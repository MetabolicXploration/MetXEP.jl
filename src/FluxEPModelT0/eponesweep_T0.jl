# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

# TODO: recheck this code
# This now it has being tested externally (output testing)
# but it is lacking iternal tests

function eponesweep!(epm::FluxEPModelT0)

    @extract epm: G 
    @extract epm: ad ai dd di 
    @extract epm: μd μi sd si 
    @extract epm: avd avi vai vad
    @extract epm: betad betai
    @extract epm: siteflagave_i siteflagave_d
    @extract epm: siteflagvar_i siteflagvar_d
    
    @extract epm: vd vi Σd Σi
    @extract epm: lbd lbi ubd ubi be

    minvar = config(epm, :minvar, 1e-50)
    maxvar = config(epm, :maxvar, 1e50)
    damp = config(epm, :damp, 0.9)
    
    errav, errva, errμ, errs = fill(typemin(eltype(G)), 4)

    # independent loop
    @inbounds for i in eachindex(μi)
        # Update the marginal tilted Q^(n) parameters (μ, s) from the parameters of the join Q
        newμi, newsi = newμs(Σi[i,i], ai[i], di[i], vi[i], lbi[i], ubi[i], minvar, maxvar)
        errμ = max(errμ, abs(μi[i] - newμi))
        errs = max(errs, abs(si[i] - newsi))
        μi[i] = newμi
        si[i] = newsi
        
        # Update the marginal tilted Q^(n) moments (av, va)
        newavi, newvai = newav(si[i], μi[i], avi[i], vai[i], 
            siteflagave_i[i], siteflagvar_i[i],
            lbi[i], ubi[i], minvar, maxvar
        )
        errav = max(errav, abs(avi[i] - newavi))
        errva = max(errva, abs(vai[i] - newvai))
        avi[i] = newavi
        vai[i] = newvai

        # Force moment matching between Q^(n) and Q by recomputing (a, d)
        newai, newdi = matchmom(μi[i], si[i], avi[i], vai[i], minvar, maxvar)
        ai[i] = damp * ai[i] + (1.0 - damp) * newai # modify a in epfields
        di[i] = damp * di[i] + (1.0 - damp) * newdi # modify d in epfields
    end

    # dependent loop
    @inbounds for i in eachindex(μd)   # loop  1:M

        # Update the marginal tilted Q^(n) parameters (μi, si) from the parameters of the join Q
        newμd, newsd = newμs(Σd[i,i], ad[i], dd[i], vd[i], lbd[i], ubd[i], minvar, maxvar)
        errμ = max(errμ, abs(μd[i] - newμd))
        errs = max(errs, abs(sd[i] - newsd))
        μd[i] = newμd
        sd[i] = newsd

        # Update the marginal tilted Q^(n) moments (av, va)
        newavd, newvad = newav(sd[i], μd[i], avd[i], vad[i],
            siteflagave_d[i], siteflagvar_d[i], 
            lbd[i], ubd[i], minvar, maxvar
        )
        errav = max(errav, abs(avd[i] - newavd))
        errva = max(errva, abs(vad[i] - newvad))
        avd[i] = newavd
        vad[i] = newvad

        # Force moment matching between Q^(n) and Q by recomputing (a, d)
        newad,newbd = matchmom(μd[i], sd[i], avd[i], vad[i], minvar, maxvar)
        ad[i] = damp * ad[i] + (1.0 - damp) * newad # modify a in epfields
        dd[i] = damp * dd[i] + (1.0 - damp) * newbd # modify d in epfields
    end

    # Recompute join parameters
    Gt = G'

    # Compute Σi and vi, the parameters (of the independent variables) of the full gaussian Q
    # TODO: check alloc performance
    println("-"^50)
    Dd = Diagonal(inv.(dd))
    Di = Diagonal(inv.(di))
    Σi_inv = Di + Gt * Dd * G
    state!(epm, :elapsed_eponesweep_inv) do
        return @elapsed inplaceinverse!(Σi, Σi_inv)
    end
    vi .= Σi * (Gt * Dd * (be - ad) + Di * ai - Gt * betad + betai)
    
    # Compute Σd and vd, the parameters (of the dependent variables) of the full gaussian Q
    # mul!(Σd, G * Σi, Gt) 
    mul!(Σd, G * Σi, Gt)
    mul!(vd, G, vi)
    vd .= be - vd

    return errav, errva, errμ, errs
end