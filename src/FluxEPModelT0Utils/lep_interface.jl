# -------------------------------------------------------------------
# Utils
# NOTE: Do not interface with the lep the data that will be in the constraints

import MetXBase.lb
lb(epm::FluxEPModelT0) = _depind_getindex(epm, epm.lbd, epm.lbi) .* epm.scalefact
lb(epm::FluxEPModelT0, idxs) = _depind_getindex(epm, epm.lbd, epm.lbi, idxs) .* epm.scalefact

import MetXBase.ub
ub(epm::FluxEPModelT0) = _depind_getindex(epm, epm.ubd, epm.ubi) .* epm.scalefact
ub(epm::FluxEPModelT0, idxs) = _depind_getindex(epm, epm.ubd, epm.ubi, idxs) .* epm.scalefact