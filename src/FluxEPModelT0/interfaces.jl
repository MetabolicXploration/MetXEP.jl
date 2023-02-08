## ------------------------------------------------------------------
# ep
export beta
beta(epm::FluxEPModelT0) = [epm.betad; epm.betai][epm.idxmap_inv]
function beta(epm::FluxEPModelT0, idxs) 
    idxs = rxnindex(epm, idxs)
    if length(idxs) == 1
        Nd = length(epm.betad)
        idx = epm.idxmap_inv[idxs]
        return idx <= Nd ? epm.betad[idx] : epm.betai[idx - Nd]
    end
    return beta(epm)[idxs]
end

export beta!
function beta!(epm::FluxEPModelT0, v::Vector) 
    v1 = view(v, epm.idxmap)
    Nd = length(epm.betad)
    epm.betad .= v1[1:Nd]
    epm.betai .= v1[Nd+1:end]
    return epm
end
function beta!(epm::FluxEPModelT0, idxs, val) 
    idxs = rxnindex(epm, idxs)
    if length(idxs) == 1
        Nd = length(epm.betad)
        idx = epm.idxmap_inv[idxs]
        if idx <= Nd; epm.betad[idx] = val
            else; epm.betai[idx - Nd] = val
        end
    else
        _beta = beta(epm)
        MetXBase._setindex!(_beta, idxs, val)
        beta!(epm, _beta)
    end
    return epm
end

## ------------------------------------------------------------------
# TODO: optimize this, do not create vector unnecessarily 

# Truncated Distribution
import Statistics.mean
export mean
function mean(epm::FluxEPModelT0)
    av = [epm.avd; epm.avi]
    rmul!(av, epm.scalefact)
    return av[epm.idxmap_inv] 
end
_mean(epm::FluxEPModelT0, idxs) = mean(epm)[rxnindex(epm, idxs)]
mean(epm::FluxEPModelT0, idxs::AbstractVector) = _mean(epm, idxs)
mean(epm::FluxEPModelT0, idxs) = _mean(epm, idxs)

import Statistics.var
export var

function var(epm::FluxEPModelT0)
    va = [epm.vad; epm.vai]
    rmul!(va, epm.scalefact^2); 
    return va[epm.idxmap_inv]
end
var(epm::FluxEPModelT0, idxs) = var(epm)[rxnindex(epm, idxs)]

## ------------------------------------------------------------------
# untruncated

export untruncated_mean
function untruncated_mean(epm::FluxEPModelT0)
    v = [epm.vd; epm.vi]
    rmul!(v, epm.scalefact)
    return v[epm.idxmap_inv] 
end
untruncated_mean(epm::FluxEPModelT0, idxs) = untruncated_mean(epm)[rxnindex(epm, idxs)]

export untruncated_var
function untruncated_var(epm::FluxEPModelT0)
    s = [diag(epm.Σd); diag(epm.Σi)]
    rmul!(s, epm.scalefact^2); 
    return s[epm.idxmap_inv]
end
untruncated_var(epm::FluxEPModelT0, idxs) = untruncated_var(epm)[rxnindex(epm, idxs)]

## -----
# MetNet
import MetXBase.lb
export lb
function lb(epm::FluxEPModelT0)
    _lb = [epm.lbd; epm.lbi]
    rmul!(_lb, epm.scalefact); 
    return _lb[epm.idxmap_inv]
end
lb(epm::FluxEPModelT0, idxs) = lb(epm)[rxnindex(epm, idxs)]

import MetXBase.ub
export ub
function ub(epm::FluxEPModelT0)
    _ub = [epm.ubd; epm.ubi]
    rmul!(_ub, epm.scalefact); 
    return _ub[epm.idxmap_inv]
end
ub(epm::FluxEPModelT0, idxs) = ub(epm)[rxnindex(epm, idxs)]