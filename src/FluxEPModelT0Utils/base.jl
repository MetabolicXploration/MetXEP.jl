# # LEPModel
# import Base.size
# size(lep::LEPModel) = size(lep.S)
# size(lep::LEPModel, dim) = size(lep.S, dim)

# import Base.==
# function ==(lep1::AbstractLEPModel, lep2::AbstractLEPModel)
#     rowids(lep1) == rowids(lep2) || return false
#     colids(lep1) == colids(lep2) || return false
#     lb(lep1) == lb(lep2) || return false
#     ub(lep1) == ub(lep2) || return false
#     balance(lep1) == balance(lep2) || return false
#     cost_matrix(lep1) == cost_matrix(lep2) || return false
#     return true
# end

# import Base.isequal
# isequal(lep1::LEPModel, lep2::LEPModel) = (lep1 == lep2)

import Base.hash
function hash(epm::FluxEPModelT0, h::UInt64 = UInt64(0)) 
    h += hash(:FluxEPModelT0)
    # structural 
    h = hash(epm.G, h)
    h = hash(epm.lbi, h)
    h = hash(epm.lbd, h)
    h = hash(epm.ubi, h)
    h = hash(epm.ubd, h)
    # convergence
    h = hash(epm.ai, h)
    h = hash(epm.ad, h)
    h = hash(epm.di, h)
    h = hash(epm.dd, h)
    return h
end
hash(epm::FluxEPModelT0, h::Integer) = hash(epm, UInt64(h)) 

# import Base.show
# show(io::IO, lep::LEPModel) = (println(io, "LEPModel ", size(lep)))

# import Base.big
# function Base.big(lep::LEPModel)
#     return LEPModel(lep; 
#         S = _collect_or_nothing(BigFloat, lep.S),
#         b = _collect_or_nothing(BigFloat, lep.b), 
#         lb = _collect_or_nothing(BigFloat, lep.lb), 
#         ub = _collect_or_nothing(BigFloat, lep.ub), 
#         c = _collect_or_nothing(BigFloat, lep.c), 
#         C = _collect_or_nothing(BigFloat, lep.C) 
#     )
# end