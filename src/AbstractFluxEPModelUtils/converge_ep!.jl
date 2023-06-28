function converge!(epm::AbstractFluxEPModel; 
        # callbacks
        oniter = nothing
    )

    # TODO: move to defaults
    # config
    verbose = config!!(epm, :verbose, false)    # output verbosity
    damp = config!!(epm, :damp, 0.9)            # damp ∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield
    epsconv = config!!(epm, :epsconv, 1e-6)     # convergence criterion
    maxiter = config!!(epm, :maxiter, 2000)     # maximum iteration count
    maxvar = config!!(epm, :maxvar, 1e50)       # maximum numerical variance
    minvar = config!!(epm, :minvar, 1e-50)      # minimum numerical variance
    iter0 = config!!(epm, :iter0, 1)            # the starting iteration count
    
    iter = iter0
    state!(epm, :iter, iter)
    state!(epm, :status, UNCONVERGED_STATUS)
    state!(epm, :converge_init_time, time())
    
    # sweep ep till maxiter is reached or max(errav, errvar) < epsconv
    prog = ProgressThresh{typeof(epsconv)}(epsconv; desc =  "EP  ", dt = 0.5)
    # TODO: Interface this
    # minimum iters before checking convergence
    _c = 0 
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
        _c > 5 && max_err < epsconv && (state!(epm, :status, CONVERGED_STATUS); break)

        if verbose 
            sweep_time = state(epm, :elapsed_eponesweep, 0)
            inv_time = state(epm, :elapsed_eponesweep_inv, 0)
            inv_frac = round(inv_time * 100/ sweep_time; digits = 3)

            update!(prog, max_err; showvalues = () -> [
                (:iter, state(epm, :iter, "nd")),
                (:maxiter, config(epm, :maxiter, "nd")),
                (:sweep_time, sweep_time),
                (:inv_time, string(inv_time, " [", inv_frac, " %]")),
            ])
        end

        
        state!(epm, :iter, iter += 1)
        _c += 1

    end # for iter

    verbose && finish!(prog)
    verbose && println("convergence_status: ", convergence_status(epm))

    # return
    return epm
end

## ------------------------------------------------------------------
function _parse_flag(ret::Tuple)
    length(ret) != 2 && return (false, nothing)
    e1, e2 = ret
    return (e1 isa Bool) ? (e1, e2) : (false, nothing)
end
_parse_flag(ret) = (false, nothing)
