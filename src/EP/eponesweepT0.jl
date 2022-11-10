# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function eponesweepT0!(epfields::EPFields, epalg::EPAlg, epmatT0::EPMatT0, stat = Dict())

    @extract epfields : av va a d μ s siteflagave siteflagvar
    @extract epalg : alpha beta_vec minvar maxvar epsconv damp
    @extract epmatT0 : Σd Σi G lb ub vd vi Y


    M = size(G,1)
    N = length(av)

    idxy = 1:M # dependent variables
    idxw = M+1:N # independent variables

    ad,ai = view(a,idxy), view(a,idxw)     # dep an ind prior mean (epfields)
    dd,di = view(d,idxy), view(d,idxw)     # dep an ind prior variance (epfields)
    sd,si = view(s,idxy), view(s,idxw)     # dep and ind untruncated marginals variance (epfields)
    μd,μi = view(μ,idxy), view(μ,idxw)     # dep and ind untruncated marginals mean (epfields)
    avd,avi = view(av,idxy), view(av,idxw) # dep and ind truncated marginals mean (epfields)
    vad,vai = view(va,idxy), view(va,idxw) # dep and ind truncated marginals variance (epfields)
    βd, βi = view(beta_vec, idxy), view(beta_vec, idxw)
    errav, errva, errμ, errs = fill(typemin(eltype(av)), 4)
    Gt = G'

    # Test
    # @show sum(beta_vec)  
    
    # All fields in epmat are updated from the epfields of last sweep
    # (?) covariance matrix of independent variables (epmat)
    stat[:elapsed_eponesweep_inv] = @elapsed begin
        Di = Diagonal(inv.(di))
        Dd = Diagonal(inv.(dd))
        # Σi = inv(Σᵢ⁻¹)
        Σᵢ⁻¹ = Gt * Dd * G + Di
        inplaceinverse!(Σi, Σᵢ⁻¹)
    end
    # fast_similarity_inv!(Σi, di, dd, G)
    mul!(Σd, G*Σi, Gt) # (?) covariance matrix of dependent variables (epmat)
    # Original ep
    # mul!(vi,Σi, ai ./ di - G'*(ad ./ dd)) # (?) mean vector of independent variables (epmat)
    # ep-maxent
    vi .= Σi * (Gt * Dd * (Y - ad) + Di * ai - Gt * βd + βi)

    # vd = -Gvi + b'
    mul!(vd, G, vi) # (?) mean vector of dependent variables (epmat)
    vd .= Y - vd

    for i in eachindex(μi)  # loop M+1:N
        newμi,newsi = newμs(Σi[i,i], ai[i], di[i], vi[i], lb[i+M], ub[i+M], minvar, maxvar)
        errμ = max(errμ, abs(μi[i]-newμi))
        errs = max(errs, abs(si[i]-newsi))
        μi[i] = newμi
        si[i] = newsi
        # println("μi[$(i+M)] = ", μi[i]," si[$(i+M)] = ", si[i], " Σi = ",Σi[i,i] )

        newavw,newvaw = newav(si[i], μi[i], avi[i], vai[i], siteflagave[i+M], siteflagvar[i+M],
                              lb[i+M], ub[i+M], minvar, maxvar)
        errav = max(errav,abs(avi[i]-newavw))
        errva = max(errva,abs(vai[i]-newvaw))
        avi[i] = newavw
        vai[i] = newvaw

        newai,newdi = matchmom(μi[i], si[i], avi[i], vai[i], minvar, maxvar)
        ai[i] = damp * ai[i] + (1.0-damp) * newai # modify a in epfields
        di[i] = damp * di[i] + (1.0-damp) * newdi # modify d in epfields
    end

    for i in eachindex(μd)   # loop  1:M

        newμd, newsd = newμs(Σd[i,i], ad[i], dd[i], vd[i], lb[i], ub[i], minvar, maxvar)
        errμ = max(errμ, abs(μd[i]-newμd))
        errs = max(errs, abs(sd[i]-newsd))
        μd[i] = newμd # modify μ in epfields
        sd[i] = newsd # modify s in epfields
#        println("μd[$i] = ", μd[i]," sd[$i] = ", sd[i], " Σ = ", Σd[i,i], " (",lb[i],":",ub[i],")"," ad[$i] = ",ad[i], " dd[$i] = ", dd[i])

        newavy,newvay = newav(sd[i],μd[i],avd[i],vad[i],
            siteflagave[i], siteflagvar[i], 
            lb[i], ub[i], minvar, maxvar)
        errav = max(errav, abs(avd[i]-newavy))
        errva = max(errva, abs(vad[i]-newvay))
        avd[i] = newavy # modify av in epfields
        vad[i] = newvay # modify va in epfields

        newad,newbd = matchmom(μd[i],sd[i],avd[i],vad[i],minvar,maxvar)
        ad[i] = damp * ad[i] + (1.0-damp) * newad # modify a in epfields
        dd[i] = damp * dd[i] + (1.0-damp) * newbd # modify d in epfields
    end
    return errav, errva, errμ, errs
end