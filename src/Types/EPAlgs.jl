# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)
# The original type did not include beta_maxent vector
struct EPAlg{T<:Real}
    alpha::T
    beta_vec::AbstractVector{T}
    minvar::T
    maxvar::T
    epsconv::T
    damp::T
    maxiter::Int
    verbose::Bool
end