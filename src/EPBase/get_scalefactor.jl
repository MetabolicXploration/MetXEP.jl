get_scalefactor(lb, ub) = max(maximum(abs.(lb)), maximum(abs.(ub)))
get_scalefactor(model::LEPModel) = get_scalefactor(lb(model), ub(model))