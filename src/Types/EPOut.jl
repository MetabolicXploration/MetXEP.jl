# # Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

# struct EPOut{T<:Real}
#     μ::Vector{T}
#     σ::Vector{T}
#     av::Vector{T}
#     va::Vector{T}
#     # stoierr::T
#     sol
#     status::Symbol
#     iter::Integer
# end

# status(out::EPOut) = out.status
# curriter(out::EPOut) = out.iter