# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

# TODO remove Y unused parameter
# TODO remove alpha Inf code
function prepareinput(S, lb, ub, alpha, verbose, solution, expval)

    M,N = size(S)
    verbose && M >= N && @warn("M = $M ≥ N = $N")
    all(lb .<= ub) || error("lower bound fluxes > upper bound fluxes. Consider swapping lower and upper bounds")

    verbose && println(stderr, "Analyzing a $M × $N stoichiometric matrix.")

    updatefunction = alpha == Inf ? eponesweepT0! : eponesweep!

    scalefact = max(maximum(abs.(lb)), maximum(abs.(ub)))
    epfields = isnothing(solution) ? epfields = EPFields(N, expval, eltype(S)) :
        epfields = deepcopy(solution.sol) # preserve the original solution!

    return updatefunction, scalefact, epfields
end