# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function matchmom(μ,s,av,va, minvar,maxvar)
    newb = clamp(inv(1.0/va - 1.0/s), minvar, maxvar)
    newa = av + newb * (av - μ) / s
    if isnan(newa + newb)
        # TODO: handle this with good taste
        error("Moment matching failed! a = $newa b = $newb")
    end
    return newa, newb
end