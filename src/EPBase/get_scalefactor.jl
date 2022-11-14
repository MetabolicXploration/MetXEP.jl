get_scalefactor(lb, ub) = max(maximum(abs.(lb)), maximum(abs.(ub)))
get_scalefactor(model::MetNet) = get_scalefactor(model.lb, model.ub)