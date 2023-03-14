# -------------------------------------------------------------------
# Fundamentals

import MetXBase.colids
colids(opm::AbstractFluxEPModel) = extras(opm, :colids)::Union{Nothing, Vector{String}}
import MetXBase.rowids
rowids(opm::AbstractFluxEPModel) = extras(opm, :rowids)::Union{Nothing, Vector{String}}