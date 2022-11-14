# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

export maxent_ep

"""
res=metabolicEP(S,b,lb,ub,...)


The output in res is of type `EPOut`: there are several fields:
-   `μ::Vector`: A parameter linked to the mean of the posterior probability
-   `σ::Vector`: A parameter linked to the std  of the posterior probability
-   `av::Vector`: The mean posterior probability
-   `va::Vector`: The variance of the posterior probability
-   `sol::EPFields`: The internal field status. From this value we can restart the sampling from a specific state.
-   `status::Symbol`: either ``:converged`` or ``:unconverged``.

Input (required)
----
- `S`: MxN matrix (either sparse or dense) please note that if you input a dense version, the algorithm is slighlty more efficient. Dense matrices can be create from sparse ones with `Matrix(S)`.
- `b`: a vector of M intakes/uptakes
- `lb`: a vector of length N of lower bounds.
- `ub`: a vector of length N of upper bounds.

Input (optional arguments).
----
- `alpha` (inverse temperature::``Real``): default 10^7;  the zero temperature algorithm is run setting ``beta=Inf``.
- `beta_vec` (MaxEnt distribution inverse temperature vector:: ``Vector{T}``): default ``T[]``
- `verbose` (``true`` or ``false``): default ``true``
- `damp` (∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield): default 0.9
- `epsconv` (convergence criterion): default 1e-6
- `maxiter` (maximum number of iterations): default 2000
- `maxvar`  (threshold on maximum variance): default 1e50
- `minvar`  (threshold on minimum variance): default 1e-50
- `solution` (start from solution. Is of type ``EPOut``): default: ``nothing``
- `expval` (fix to posterior probability of mean and/or variance to values): default ``nothing``. expval can be either at ``Tuple{Float64,Float64,Int}`` or a ``Vector{Tuple{Float64,Float64,Int}}``. Values can be fixed as``expval=(0.2,0.4,4)`` meaning that for flux index 4 the mean is set to 0.2 and the variance to 0.4. Fixing more values ``expval=[(0.2, 0.3, 4), (0.4, nothing, 5)]``: in this case, we fix the posterior of flux 4 to 0.2 (mean) and 0.3 (variance), while for flux 5 we fix the mean to 0.4 and we keep the variance free.
"""
function maxent_ep(
        S::AbstractArray{T,2}, b::AbstractArray{T,1}, lb::AbstractArray{T,1}, ub::AbstractArray{T,1};
        alpha::Real=Inf,                             # inverse temperature
        beta_vec::AbstractVector{T} = T[],           # maxent inverse temperature vector
        verbose::Bool=true,                          # output verbosity
        damp::Real=0.9,                              # damp ∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield
        epsconv::Real=1e-6,                          # convergence criterion
        maxiter::Int=2000,                           # maximum iteration count
        maxvar::Real=1e50,                           # maximum numerical variance
        minvar::Real=1e-50,                          # minimum numerical variance
        solution = nothing,                          # start from a solution
        expval=nothing,                              # fix posterior probability experimental values for std and mean
        iter0 = 0,                                   # the started iteration count
        drop_epfields = false,                       # if the final EPOut object will export the epfields
        # callbacks
        oniter::Function = (it, epmodel) -> (false, nothing)
    ) where T<:Real

    epmodel = FluxEPModel(S, b, lb, ub; alpha, beta_vec, solution, expval)
    return converge_ep!(epmodel; verbose, damp, epsconv, 
        maxiter, maxvar, minvar, iter0, drop_epfields, oniter
    )
end

maxent_ep(net::MetNet; kwargs...) = 
    maxent_ep(net.S, net.b, net.lb, net.ub; kwargs...)
