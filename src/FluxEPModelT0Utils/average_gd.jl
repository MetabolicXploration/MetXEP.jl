# TODO: Use Optim
## ---------------------------------------------
function beta_gd!(obj_fun::Function, epm::FluxEPModelT0, target_ids::Vector, target_value::Vector;
        x0::Vector = zeros(length(target_ids)),  # initial betas
        maxΔx::Vector = ones(length(target_ids)) .* 1e3,
        minΔx::Vector = ones(length(target_ids)) .* 1e-5,
        x1::Vector = ones(length(target_ids)), # The first step
        gdth = 1e-5,
        maxiter = 5000,
        verbose = false,
        gd_kwargs...
    )

    # setup
    config!(epm; verbose = false) # silence ep
    target_ids = colindex(epm, target_ids)

    # gd
    function up_fun(gdmodel) 
        # Up betas
        _curr_betas = MetXBase.gd_value(gdmodel)
        beta!(epm, target_ids, _curr_betas)
        # converge
        converge!(epm)
        # return obj_fun
        return obj_fun(gdmodel)
    end

    gdmodel = MetXBase.grad_desc_vec(up_fun; 
        target = target_value, x0, x1, minΔx, maxΔx, gdth, maxiter, verbose, gd_kwargs...
    )

    return gdmodel

end

_beta_gd_deflt_obj_fun(epm, target_ids) = (x...) -> mean(epm, target_ids)

function beta_gd!(epm::FluxEPModelT0, target_ids::Vector, args...; kwargs...) 
    target_ids = colindex(epm, target_ids)
    beta_gd!(_beta_gd_deflt_obj_fun(epm, target_ids), epm, target_ids, args...; kwargs...)
end

# TODO: Use Optim
## ---------------------------------------------
function gamma_gd!(obj_fun::Function, epm::FluxEPModelT0, target_ids::Vector, target_value::Vector;
        x0::Vector = zeros(length(target_ids)),  # initial gammas
        maxΔx::Vector = ones(length(target_ids)) .* 1e3,
        minΔx::Vector = ones(length(target_ids)) .* 1e-5,
        x1::Vector = ones(length(target_ids)), # The first step
        gdth = 1e-5,
        maxiter = 5000,
        verbose = false,
        gd_kwargs...
    )

    # setup
    config!(epm; verbose = false) # silence ep
    target_ids = colindex(epm, target_ids)

    # gd
    function up_fun(gdmodel) 
        # Up gammas
        _curr_gammas = MetXBase.gd_value(gdmodel)
        gamma!(epm, target_ids, _curr_gammas)
        # converge
        converge!(epm)
        # return means
        # return mean(epm, target_ids).^2
        return obj_fun(gdmodel)
    end

    gdmodel = MetXBase.grad_desc_vec(up_fun; 
        target = target_value, x0, x1, minΔx, maxΔx, gdth, maxiter, verbose, gd_kwargs...
    )

    return gdmodel

end

_gamma_gd_deflt_obj_fun(epm, target_ids) = (x...) -> var(epm, target_ids)

function gamma_gd!(epm::FluxEPModelT0, target_ids::Vector, args...; kwargs...)
    target_ids = colindex(epm, target_ids)
    gamma_gd!(_gamma_gd_deflt_obj_fun(epm, target_ids), epm, target_ids, args...; kwargs...)
end