# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function eponesweep!(epm::FluxEPModelT0)

    # dep an ind prior mean (epfields)
    # dep an ind prior variance (epfields)
    # dep and ind truncated gaussian-part marginals variance 
    # dep and ind truncated gaussian-part marginals mean 
    # dep and ind truncated marginals mean (epfields)
    # dep and ind truncated marginals mean (epfields)
    # dep and ind truncated marginals variance (epfields)

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
    
    M, N = size(G)
    errav, errva, errμ, errs = fill(typemin(eltype(G)), 4)
    Gt = G'
    
    # All fields in epmat are updated from the epfields of last sweep
    # (?) covariance matrix of independent variables (epmat)
    elapsed_eponesweep_inv = @elapsed begin
        Di = Diagonal(inv.(di))
        Dd = Diagonal(inv.(dd))
        # Σi = inv(Σᵢ⁻¹)
        # Σᵢ⁻¹ = Gt * Dd * G + Di
        Σᵢ⁻¹ = Di + Gt * Dd * G
        inplaceinverse!(Σi, Σᵢ⁻¹)
    end
    state!(epm; elapsed_eponesweep_inv)
    # fast_similarity_inv!(Σi, di, dd, G)
    mul!(Σd, G * Σi, Gt) # (?) covariance matrix of dependent variables (epmat)
    # Original ep
    # mul!(vi,Σi, ai ./ di - G'*(ad ./ dd)) # (?) mean vector of independent variables (epmat)
    # ep-maxent
    vi .= Σi * (Gt * Dd * (be - ad) + Di * ai - Gt * betad + betai)

    # vd = -Gvi + b'
    mul!(vd, G, vi) # (?) mean vector of dependent variables (epmat)
    vd .= be - vd

    @inbounds for i in eachindex(μi)
        newμi, newsi = newμs(Σi[i,i], ai[i], di[i], vi[i], lbi[i], ubi[i], minvar, maxvar)
        errμ = max(errμ, abs(μi[i] - newμi))
        errs = max(errs, abs(si[i] - newsi))
        μi[i] = newμi
        si[i] = newsi
        
        newavw, newvaw = newav(si[i], μi[i], avi[i], vai[i], 
            siteflagave_i[i], siteflagvar_i[i],
            lbi[i], ubi[i], minvar, maxvar
        )
        errav = max(errav, abs(avi[i] - newavw))
        errva = max(errva, abs(vai[i] - newvaw))
        avi[i] = newavw
        vai[i] = newvaw

        newai,newdi = matchmom(μi[i], si[i], avi[i], vai[i], minvar, maxvar)
        ai[i] = damp * ai[i] + (1.0 - damp) * newai # modify a in epfields
        di[i] = damp * di[i] + (1.0 - damp) * newdi # modify d in epfields
    end

    @inbounds for i in eachindex(μd)   # loop  1:M

        newμd, newsd = newμs(Σd[i,i], ad[i], dd[i], vd[i], lbd[i], ubd[i], minvar, maxvar)
        errμ = max(errμ, abs(μd[i] - newμd))
        errs = max(errs, abs(sd[i] - newsd))
        μd[i] = newμd # modify μ in epfields
        sd[i] = newsd # modify s in epfields

        newavy, newvay = newav(sd[i],μd[i],avd[i],vad[i],
            siteflagave_d[i], siteflagvar_d[i], 
            lbd[i], ubd[i], minvar, maxvar
        )
        errav = max(errav, abs(avd[i] - newavy))
        errva = max(errva, abs(vad[i] - newvay))
        avd[i] = newavy # modify av in epfields
        vad[i] = newvay # modify va in epfields

        newad,newbd = matchmom(μd[i],sd[i],avd[i],vad[i],minvar,maxvar)
        ad[i] = damp * ad[i] + (1.0 - damp) * newad # modify a in epfields
        dd[i] = damp * dd[i] + (1.0 - damp) * newbd # modify d in epfields
    end

    return errav, errva, errμ, errs
end