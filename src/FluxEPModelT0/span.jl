# DEPRECATED: The EP scaling makes this unsafe

# import MetXBase.span!
# export span!
# function span!(v::Vector, epm::FluxEPModelT0, vf::Vector)
#     v[epm.idxd] .= epm.be - epm.G * vf
#     v[epm.idxi] .= vf
#     return v
# end

# import MetXBase.span
# export span
# span(epm::FluxEPModelT0, vf::Vector) = span!(zeros(sum(size(epm.G))), epm, vf)
