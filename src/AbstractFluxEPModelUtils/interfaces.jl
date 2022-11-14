# -------------------------------------------------------------------
# extras interface
import MetXBase.extras
extras(m::AbstractFluxEPModel) = m.extras

# -------------------------------------------------------------------
# state interface
@extras_dict_interface AbstractFluxEPModel state

# -------------------------------------------------------------------
# config interface
@extras_dict_interface AbstractFluxEPModel config

# -------------------------------------------------------------------
# ep interface

# -------------------------------------------------------------------
# net interface
# NOTE: Do not interface with the net the data that will be in the model

# net data
@extras_val_interface AbstractFluxEPModel metnet MetNet

import MetXBase.metabolites
metabolites(m::AbstractFluxEPModel, ider...) = metabolites(metnet(m), ider...)

import MetXBase.reactions
reactions(m::AbstractFluxEPModel, ider...) = reactions(metnet(m), ider...)

import MetXBase.genes
genes(m::AbstractFluxEPModel, ider...) = genes(metnet(m), ider...)

# contraints data
import MetXBase.lb
# import MetXBase.lb!
# lb(m::AbstractFluxEPModel) = JuMP.normalized_rhs.(get_lower_bound_cons(m))
# lb(m::AbstractFluxEPModel, ridx) = lb(m)[rxnindex(m, ridx)]
# function lb!(m::AbstractFluxEPModel, idxs, lb)
#     upcons = get_lower_bound_cons(m)
#     cidxs = rxnindex(m, idxs)
#     up_con_rhs!(upcons, lb, cidxs)
#     return m
# end
# lb!(m::AbstractFluxEPModel, lb) = lb!(m, 1:_length(m), lb)

# ub
import MetXBase.ub
# import MetXBase.ub!
# ub(m::AbstractFluxEPModel) = JuMP.normalized_rhs.(get_upper_bound_cons(m))
# ub(m::AbstractFluxEPModel, ider) = ub(m)[rxnindex(m, ider)]
# function ub!(m::AbstractFluxEPModel, idxs, ub)
#     upcons = get_upper_bound_cons(m)
#     cidxs = rxnindex(m, idxs)
#     up_con_rhs!(upcons, ub, cidxs)
#     return m
# end
# ub!(m::AbstractFluxEPModel, ub) = ub!(m, 1:_length(m), ub)

import MetXBase.bounds
# import MetXBase.bounds!
# bounds(m::AbstractFluxEPModel) = (lb(m), ub(m))
# bounds(m::AbstractFluxEPModel, idx) = (lb(m, idx), ub(m, idx))
# bounds!(m::AbstractFluxEPModel, idx, lb, ub) = (lb!(m, idx, lb); ub!(m, idx, ub); nothing)
# function bounds!(m::AbstractFluxEPModel, idx; lb = nothing, ub = nothing) 
#     isnothing(lb) || lb!(m, idx, lb)
#     isnothing(ub) || ub!(m, idx, ub)
#     return nothing
# end

# import MetXBase.balance
# balance(m::AbstractFluxEPModel) = JuMP.normalized_rhs.(get_balance_cons(m))
# balance(m::AbstractFluxEPModel, ider) = balance(m)[metindex(m, ider)]

