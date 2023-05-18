# activate env
using RunTestsEnv
@activate_testenv

using MetXEP
using MetXNetHub
using MetXGEMs
using MetXOptim
using MetXOptim.GLPK
using MetXBase
using MetXBase.UnicodePlots
using LinearAlgebra
using Test

const LIN_SOLVER = GLPK.Optimizer

@testset "MetXEP.jl" begin
    
    ## ------------------------------------------------------------------
    # SETTER/GETTER interface (indirectly test idxmap and ixdmap_inv)
    let
        println()
        println("="^60)
        println("SETTER/GETTER INTERFACE")
        println("."^60)
        println()
        
        net0 = pull_net("toy_net")
        lep = box(net0, LIN_SOLVER)
        epm = FluxEPModelT0(lep)
        # beta

        for rxn in colids(lep)
            beta!(epm, rxn, 1)
            @test findfirst(isone, beta(epm)) == colindex(epm, rxn)
            beta!(epm, rxn, 0)
            @test count(isone, beta(epm)) |> iszero
        end

        @test isapprox(lb(epm), lb(lep))
        @test isapprox(ub(epm), ub(lep))

        println()
    end

    ## ------------------------------------------------------------------
    # EP -> FBA
    let
        println()
        println("="^60)
        println("EP -> FBA")
        println("."^60)
        println()

        model_id = "ecoli_core"
        net0 = pull_net(model_id)
        biom_id = extras(net0, "BIOM")
        glc_id = extras(net0, "EX_GLC")
        
        lep = box(net0, LIN_SOLVER; eps = 1e-4)
        M, N = size(lep)
        Srank = rank(lep.S)

        opm = fba(lep, LIN_SOLVER)

        epm = FluxEPModelT0(lep)
        config!(epm; 
            verbose = false, 
            maxiter = 3000, 
            damp = 0.9, 
            epsconv = 1e-6,
        )

        betas = 10.0.^range(-5, 9.0; length = 10)
        fba_dists = Float64[]
        biom_avs = Float64[]
        iters = Int[]
        for b in betas
            beta!(epm, biom_id, b)
            converge!(epm)
            push!(iters, state(epm, :iter))

            diffs = (solution(opm) .- mean(epm)).^2
            sort!(diffs)
            push!(fba_dists, sum(diffs[1:Srank]))

            push!(biom_avs, mean(epm, biom_id))
        end

        # Approx tests
        @test last(fba_dists) <= 1e-4
        @test mean(epm, biom_id) - solution(opm, biom_id) <= 1e-4
        
        # Plots
        # FBA tot dist
        p = lineplot(log10.(betas), log10.(fba_dists);
            title = "EP -> FBA", 
            name = "EP", 
            xlabel = "log10 (beta)", 
            ylabel = "log10 dist(ep, fba)",
            canvas = DotCanvas
        )
        println(p)
        println()

        # BIOMS tot dist
        p = lineplot(log10.(betas), fill(solution(opm, biom_id), length(betas));
            # name = "EP", 
            name = "FBA", 
            title = model_id, 
            xlabel = "log10 (beta)", 
            ylabel = "log10 biom av",
            canvas = DotCanvas
        )
        lineplot!(p, log10.(betas), biom_avs; name="EP")
        println(p)

        # Iters
        p = lineplot(log10.(betas), iters;
            title = model_id, 
            xlabel = "log10 (beta)", 
            ylabel = "iters",
            canvas = DotCanvas
        )
        println(p)
    end

    ## ------------------------------------------------------------------
    # FluxEPModelT0 Normal Join Parameters
    let
        model_id = "ecoli_core"
        net0 = pull_net(model_id)
        lep = box(net0, LIN_SOLVER; eps = 1e-4)

        epm = FluxEPModelT0(lep)
        config!(epm; 
            damp = 0.90,
            verbose = true,
            epsconv = 1e-6, 
            maxiter = Int(1e4)
        )
        converge!(epm)
        @show convergence_status(epm)

        # Join
        G = epm.G
        βi, βd = epm.betai, epm.betad
        di, dd = epm.di, epm.dd
        Di, Dd = Diagonal(inv.(di)), Diagonal(inv.(dd))
        ai, ad = epm.ai, epm.ad
        b = epm.be
        vi, Σi = epm.vi, epm.Σi

        Σi0inv = Di + G' * Dd * G
        Σi0 = similar(Σi0inv)
        MetXEP.inplaceinverse!(Σi0, Σi0inv)
        @show isposdef(Σi)
        @show isposdef(Σi0)
        MetXBase.nearPD!(Σi0)

        @show maximum(abs, vec(Σi .- Σi0))
        @test isapprox(Σi, Σi0; rtol = 1e-2)
        
        vi0 = Σi0 * (G' * Dd * (b - ad) + Di * ai - G' * βd + βi)
        @show maximum(abs, vec(vi .- vi0))
        @test isapprox(vi, vi0; rtol = 1e-2)

    end

    ## ------------------------------------------------------------------
    # Average grad desc
    let
        ## ---------------------------------------------
        model_id = "ecoli_core"
        net0 = pull_net(model_id)
        lep = box(net0, LIN_SOLVER; verbose = false)
        biom_id = extras(lep, "BIOM") 
        glc_id = extras(lep, "EX_GLC")
        lac_id = "EX_lac__D_e"

        ## ---------------------------------------------
        # EP MODEL
        epm = FluxEPModelT0(lep)
        config!(epm; epsconv = 1e-4, verbose = true)
        converge!(epm)

        ## ---------------------------------------------
        target_ids = [biom_id, glc_id, lac_id]
        target_avs = [0.2, -9.0, 0.1]
        gdmodel = average_gd!(epm, target_ids, target_avs;
            maxΔx = [1e3, 1e3, 1e4],
            gdth = 1e-2,
            maxiter = 500,
            verbose = true
        )

        ## ---------------------------------------------
        @test isapprox(mean(epm, target_ids), target_avs; rtol = 1e-2)
    end

end
