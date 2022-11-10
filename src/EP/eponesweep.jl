# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function eponesweep!(epfields::EPFields{T}, epalg::EPAlg, epmat::EPMat, stat) where T
    @extract epfields : av va a d μ s siteflagave siteflagvar
    @extract epalg : alpha beta_vec minvar maxvar epsconv damp
    @extract epmat : αStS αStS_diag Σ lb ub αStb v


    # Parameters of the approximated join
    # We rebuild αStS and add D
    # αStS = αStS if diag(αStS) = αStS_diag
    # Σ^-1 = (αStS + D)
    Σ⁻¹ = αStS[diagind(αStS)] = αStS_diag + inv.(d)

    # Σ = inv(Σ^-1)
    stat[:elapsed_eponesweep_inv] = @elapsed inplaceinverse!(Σ, Σ⁻¹)

    # αStb = βF'y
    # v¯ = Σ(αStb + Da) (original ep)
    # mul!(v, Σ, (αStb + a./d))
    # v¯ = Σ(αStb + Da) + Σ * beta_vec (maxent-ep)
    v .= Σ * (αStb + a./d) + Σ * beta_vec

    errav = errva = errμ = errs = typemin(T)
    for i in eachindex(av)

        # Parameters of the normal part of the nth tilted
        newμ, news = newμs(Σ[i,i],a[i],d[i],v[i],lb[i],ub[i],minvar, maxvar)
        errμ = max(errμ, abs(μ[i]-newμ))
        errs = max(errs, abs(s[i]-news))
        μ[i] = newμ
        s[i] = news

        # Parameters of the whole tilted (truncated marginals)
        newave, newva = newav(s[i],μ[i],av[i],va[i],siteflagave[i],siteflagvar[i],lb[i],ub[i],minvar,maxvar)
        errav = max(errav,abs(av[i]-newave))
        errva = max(errva,abs(va[i]-newva))
        av[i] = newave
        va[i] = newva

        # Parameter of the prior
        newa,newb = matchmom(μ[i],s[i],av[i],va[i],minvar,maxvar)
        a[i] = damp * a[i] + (1.0-damp)*newa
        d[i] = damp * d[i] + (1.0-damp)*newb
    end

    return errav,errva,errμ,errs
end