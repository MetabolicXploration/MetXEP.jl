# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function newμs(Σ, a, d, v, lb, ub, minvar, maxvar)

    Σ == 0 && (Σ = minvar)
    #lΣ = clamp(Σ,minvar,maxvar)
    #s = Σ > 0 ? clamp(inv(1.0/Σ - 1.0/d), minvar, maxvar) : minvar
    s = clamp(inv(inv(Σ) - inv(d)), minvar, maxvar)
    if Σ != d
        μ = s * (v/Σ - a/d)
    else
        #@warn("I'm here: ub = ",ub," lb = ",lb, " Σ = ", Σ)
        μ = 0.5 * (ub + lb)
    end
    return μ, s
end

#
