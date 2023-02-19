# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function matchmom(μ,s,av,va, minvar,maxvar)
    newb = clamp(inv(1.0/va - 1.0/s),minvar,maxvar)
    newa = av + newb*(av-μ)/s
    if isnan(newa) || isnan(newb)
        @warn("a = $newa b = $newb")
    end
    return newa, newb
end