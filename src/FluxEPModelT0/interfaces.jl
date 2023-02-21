## ------------------------------------------------------------------
_depind_getindex(epm::FluxEPModelT0, cold, coli) = [cold; coli][epm.idxmap_inv]
function _depind_getindex(epm::FluxEPModelT0, cold, coli, idx::Int)
    Nd = length(epm.idxd)
    inv_idx = epm.idxmap_inv[idx]
    inv_idx <= Nd ? cold[inv_idx] : coli[inv_idx - Nd]
end

_depind_getindex(epm::FluxEPModelT0, cold, coli, idxs) = 
    [_depind_getindex(epm, cold, coli, idx) for idx in rxnindex(epm, idxs)]

function _depind_setindex!(epm::FluxEPModelT0, cold, coli, v)
    v1 = view(v, epm.idxmap)
    Nd = length(epm.idxd)
    cold .= v1[1:Nd]
    coli .= v1[Nd+1:end]
    return epm
end

function _depind_setindex!(epm::FluxEPModelT0, cold, coli, idxs, val)
    _col = _depind_getindex(epm, cold, coli)
    idxs = rxnindex(epm, idxs)
    MetXBase._setindex!(_col, idxs, val)
    _depind_setindex!(epm, cold, coli, _col)
    return epm
end

# ep
export beta
beta(epm::FluxEPModelT0) = _depind_getindex(epm, epm.betad, epm.betai)
beta(epm::FluxEPModelT0, idxs) = _depind_getindex(epm, epm.betad, epm.betai, idxs)

export beta!
beta!(epm::FluxEPModelT0, v::Vector) = _depind_setindex!(epm, epm.betad, epm.betai, v)
beta!(epm::FluxEPModelT0, idxs, val) = _depind_setindex!(epm, epm.betad, epm.betai, idxs, val)

## ------------------------------------------------------------------
# internals
_a(epm::FluxEPModelT0) = _depind_getindex(epm, epm.ad, epm.ai) .* epm.scalefact
_a(epm::FluxEPModelT0, idxs) = _depind_getindex(epm, epm.ad, epm.ai, idxs) .* epm.scalefact
_d(epm::FluxEPModelT0) = _depind_getindex(epm, epm.dd, epm.di) .* epm.scalefact^2
_d(epm::FluxEPModelT0, idxs) = _depind_getindex(epm, epm.dd, epm.di, idxs) .* epm.scalefact^2
_μ(epm::FluxEPModelT0) = _depind_getindex(epm, epm.μd, epm.μi) .* epm.scalefact
_μ(epm::FluxEPModelT0, idxs) = _depind_getindex(epm, epm.μd, epm.μi, idxs) .* epm.scalefact
_s(epm::FluxEPModelT0) = _depind_getindex(epm, epm.sd, epm.si) .* epm.scalefact^2
_s(epm::FluxEPModelT0, idxs) = _depind_getindex(epm, epm.sd, epm.si, idxs) .* epm.scalefact^2
_Σi(epm::FluxEPModelT0) = epm.Σi .* epm.scalefact^2
_vi(epm::FluxEPModelT0) = epm.vi .* epm.scalefact

## ------------------------------------------------------------------
# Truncated Distribution
import Statistics.mean
export mean
mean(epm::FluxEPModelT0) = _depind_getindex(epm, epm.avd, epm.avi) .* epm.scalefact    
_mean(epm::FluxEPModelT0, idxs) = _depind_getindex(epm, epm.avd, epm.avi, idxs) .* epm.scalefact
mean(epm::FluxEPModelT0, idxs::AbstractArray) = _mean(epm, idxs)
mean(epm::FluxEPModelT0, idxs::Integer) = _mean(epm, idxs)
mean(epm::FluxEPModelT0, idxs) = _mean(epm, idxs)

import Statistics.var
export var

var(epm::FluxEPModelT0) = _depind_getindex(epm, epm.vad, epm.vai) .* epm.scalefact^2
var(epm::FluxEPModelT0, idxs) = _depind_getindex(epm, epm.vad, epm.vai, idxs) .* epm.scalefact^2

## ------------------------------------------------------------------
# untruncated

export untruncated_mean
untruncated_mean(epm::FluxEPModelT0) = _depind_getindex(epm, epm.vd, epm.vi) .* epm.scalefact
untruncated_mean(epm::FluxEPModelT0, idxs) = _depind_getindex(epm, epm.vd, epm.vi, idxs) .* epm.scalefact

export untruncated_var
untruncated_var(epm::FluxEPModelT0) = _depind_getindex(epm, diag(epm.Σd), diag(epm.Σi)) .* epm.scalefact^2
untruncated_var(epm::FluxEPModelT0, idxs) = _depind_getindex(epm, diag(epm.Σd), diag(epm.Σi), idxs) .* epm.scalefact^2

## -----
# MetNet
import MetXBase.lb
export lb
lb(epm::FluxEPModelT0) = _depind_getindex(epm, epm.lbd, epm.lbi) .* epm.scalefact
lb(epm::FluxEPModelT0, idxs) = _depind_getindex(epm, epm.lbd, epm.lbi, idxs) .* epm.scalefact

import MetXBase.ub
export ub
ub(epm::FluxEPModelT0) = _depind_getindex(epm, epm.ubd, epm.ubi) .* epm.scalefact
ub(epm::FluxEPModelT0, idxs) = _depind_getindex(epm, epm.ubd, epm.ubi, idxs) .* epm.scalefact