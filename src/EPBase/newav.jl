# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

# new Q_n tilted first moment
function newav(s, μ, av, va, siteflagave, siteflagvar, lb, ub, minvar, maxvar)
    sqrts = sqrt(s)
    xinf = (lb - μ) / sqrts
    xsup = (ub - μ) / sqrts
    scra1 , scra12 = compute_mom5d(xinf, xsup)
    avnew  = siteflagave ? μ + scra1 * sqrts : av # if is fixed use av
    varnew = siteflagvar ? max(minvar, s * (1.0 + scra12)) : va
    if isnan(avnew) || isnan(varnew) 
        error("Tilted moments calculation failed, avnew = $avnew varnew = $varnew")
    end
    return avnew, varnew
end