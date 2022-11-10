# -------------------------------------------------------------------
# extras interface
import MetXBase.extras
extras(m::FluxEPModel) = m.extras

# -------------------------------------------------------------------
# config interface

# -------------------------------------------------------------------
# ep interface

# -------------------------------------------------------------------
# net interface
# NOTE: Do not interface with the net the data that will be in the constraints

# net data
metnet(opm::FluxEPModel)::Union{Nothing, MetNet} = extras(opm, :_net, nothing)
metnet!(opm::FluxEPModel, net::MetNet) = extras!(opm, :_net, net) 

import MetXBase.metabolites
metabolites(m::FluxEPModel, ider...) = metabolites(metnet(m), ider...)

import MetXBase.reactions
reactions(m::FluxEPModel, ider...) = reactions(metnet(m), ider...)

import MetXBase.genes
genes(m::FluxEPModel, ider...) = genes(metnet(m), ider...)

# contraints data
import MetXBase.lb
# import MetXBase.lb!
lb(m::FluxEPModel) = JuMP.normalized_rhs.(get_lower_bound_cons(m))
lb(m::FluxEPModel, ridx) = lb(m)[rxnindex(m, ridx)]
# function lb!(m::FluxEPModel, idxs, lb)
#     upcons = get_lower_bound_cons(m)
#     cidxs = rxnindex(m, idxs)
#     up_con_rhs!(upcons, lb, cidxs)
#     return m
# end
# lb!(m::FluxEPModel, lb) = lb!(m, 1:_length(m), lb)

# ub
import MetXBase.ub
# import MetXBase.ub!
ub(m::FluxEPModel) = JuMP.normalized_rhs.(get_upper_bound_cons(m))
ub(m::FluxEPModel, ider) = ub(m)[rxnindex(m, ider)]
# function ub!(m::FluxEPModel, idxs, ub)
#     upcons = get_upper_bound_cons(m)
#     cidxs = rxnindex(m, idxs)
#     up_con_rhs!(upcons, ub, cidxs)
#     return m
# end
# ub!(m::FluxEPModel, ub) = ub!(m, 1:_length(m), ub)

import MetXBase.bounds
# import MetXBase.bounds!
bounds(m::FluxEPModel) = (lb(m), ub(m))
bounds(m::FluxEPModel, idx) = (lb(m, idx), ub(m, idx))
# bounds!(m::FluxEPModel, idx, lb, ub) = (lb!(m, idx, lb); ub!(m, idx, ub); nothing)
# function bounds!(m::FluxEPModel, idx; lb = nothing, ub = nothing) 
#     isnothing(lb) || lb!(m, idx, lb)
#     isnothing(ub) || ub!(m, idx, ub)
#     return nothing
# end

# import MetXBase.balance
# balance(m::FluxEPModel) = JuMP.normalized_rhs.(get_balance_cons(m))
# balance(m::FluxEPModel, ider) = balance(m)[metindex(m, ider)]

