using MetXNetHub
using MetXOptim
using MetXOptim: GLPK
using MetXBase
using MetXBase.UnicodePlots
using LinearAlgebra
using MetXEP
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
        
        net = pull_net("toy_net")
        net = box(net, LIN_SOLVER)
        epm = FluxEPModelT0(net)
        # beta

        for rxn in reactions(net)
            beta!(epm, rxn, 1)
            @test findfirst(isone, beta(epm)) == rxnindex(epm, rxn)
            beta!(epm, rxn, 0)
            @test count(isone, beta(epm)) |> iszero
        end

        @test isapprox(lb(epm), lb(net))
        @test isapprox(ub(epm), ub(net))

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
        net = pull_net(model_id)
        biom_id = extras(net, "BIOM")
        glc_id = extras(net, "EX_GLC")
        
        net = box(net, GLPK.Optimizer; eps = 1e-4)
        M, N = size(net)
        Srank = rank(net.S)

        opm = fba(net, GLPK.Optimizer)

        epm = FluxEPModelT0(net)
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


end
