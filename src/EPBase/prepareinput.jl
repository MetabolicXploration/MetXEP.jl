# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

# TODO remove Y unused parameter
# TODO remove alpha Inf code
function prepareinput(S, lb, ub, alpha, verbose, solution, expval)

    M,N = size(S)
    verbose && M >= N && @warn("M = $M ≥ N = $N")
    for (i, (l, u)) in enumerate(zip(lb, ub))
        l <= u && continue
        error("lower bound [$l] >= upper bound [$u] at index $i")
    end

    verbose && println(stderr, "Analyzing a $M x $N stoichiometric matrix.")

    updatefunction = alpha == Inf ? eponesweepT0! : eponesweep!

    scalefact = max(maximum(abs.(lb)), maximum(abs.(ub)))
    epfields = isnothing(solution) ? epfields = EPFields(N, expval, eltype(S)) :
        epfields = deepcopy(solution.sol) # preserve the original solution!

    return updatefunction, scalefact, epfields
end