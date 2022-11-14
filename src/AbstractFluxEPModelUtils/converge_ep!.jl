export converge!
function converge!(epm::AbstractFluxEPModel;
        # config
        verbose::Bool=true,  # output verbosity
        # state
        damp::Real=0.9,      # damp ∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield
        epsconv::Real=1e-6,  # convergence criterion
        maxiter::Int=2000,   # maximum iteration count
        maxvar::Real=1e50,   # maximum numerical variance
        minvar::Real=1e-50,  # minimum numerical variance
        iter0 = 1,           # the starting iteration count
        # callbacks
        oniter = nothing
    )

    # TODO: move to defaults
    # config
    verbose = config(epm, :verbose, verbose)
    damp = config(epm, :damp, damp)
    epsconv = config(epm, :epsconv, epsconv)
    maxiter = config(epm, :maxiter, maxiter)
    maxvar = config(epm, :maxvar, maxvar)
    minvar = config(epm, :minvar, minvar)
    
    iter0 = config(epm, :iter0, iter0)
    iter = iter0
    state!(epm, :iter, iter)
    state!(epm, :status, UNCONVERGED_STATUS)
    state!(epm, :converge_init_time, time())
    
    # sweep ep till maxiter is reached or max(errav, errvar) < epsconv
    prog = ProgressThresh{typeof(epsconv)}(epsconv; desc =  "EP  ", dt = 0.5)
    # max_beta = findmax(epm.beta)
    for _ in iter0:maxiter

        elapsed_eponesweep = @elapsed begin
            errs = eponesweep!(epm)
        end
        
        max_err = maximum(errs)
        state!(epm; elapsed_eponesweep, max_err)

        # call back
        retflag = MetXBase.run_callbacks(oniter, epm)
        retflag === true && return epm

        # Converged
        max_err < epsconv && (state!(epm, :status, CONVERGED_STATUS); break)

        if verbose 
            sweep_time = state(epm, :elapsed_eponesweep, 0)
            inv_time = state(epm, :elapsed_eponesweep_inv, 0)
            inv_frac = round(inv_time * 100/ sweep_time; digits = 3)

            update!(prog, max_err; showvalues = [
                (:iter, state(epm, :iter, "nd")),
                (:maxiter, config(epm, :maxiter, "nd")),
                # (:alpha, alpha),
                # (:max_beta, max_beta),
                (:sweep_time, sweep_time),
                (:inv_time, string(inv_time, " [", inv_frac, " %]")),
            ])
        end

        
        state!(epm, :iter, iter += 1)

    end # for iter

    verbose && finish!(prog)

    # return
    # return produce_epout!(epm, status, iter; drop_epfields)
    return epm
end

## ------------------------------------------------------------------
function _parse_flag(ret::Tuple)
    length(ret) != 2 && return (false, nothing)
    e1, e2 = ret
    return (e1 isa Bool) ? (e1, e2) : (false, nothing)
end
_parse_flag(ret) = (false, nothing)
