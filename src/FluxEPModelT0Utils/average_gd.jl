## ---------------------------------------------
function average_gd!(epm::FluxEPModelT0, target_ids::Vector, target_value::Vector;
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
        # return means
        sleep(0.1)
        return mean(epm, target_ids)
    end

    gdmodel = MetXBase.grad_desc_vec(up_fun; 
        target = target_value, x0, x1, minΔx, maxΔx, gdth, maxiter, verbose, gd_kwargs...
    )

    return gdmodel

end
