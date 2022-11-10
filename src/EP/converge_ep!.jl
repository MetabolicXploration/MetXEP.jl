
function converge_ep!(epmodel::FluxEPModel{T};
        verbose::Bool=true,  # output verbosity
        damp::Real=0.9,      # damp ∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield
        epsconv::Real=1e-6,  # convergence criterion
        maxiter::Int=2000,   # maximum iteration count
        maxvar::Real=1e50,   # maximum numerical variance
        minvar::Real=1e-50,  # minimum numerical variance
        iter0 = 1,           # the started iteration count
        drop_epfields = false,  # if the final EPOut object will export the epfields
        # callbacks
        oniter::Function = (it, epmodel) -> nothing
    ) where {T<:Real}

    @extract epmodel : epfields epmat beta_vec alpha

    stat = extras!(() -> Dict(), epmodel, :_stat)
    stat[:converge_init_time] = time()
    epalg = EPAlg(alpha, beta_vec, 
        minvar, maxvar, epsconv, damp, maxiter, verbose)

    returnstatus = UNCONVERGED_STATUS
    iter = iter0
    
    # sweep ep till maxiter is reached or max(errav, errvar) < epsconv
    prog = ProgressThresh{typeof(epsconv)}(epsconv; desc =  "EP  ", dt = 0.5)
    max_beta = findmax(beta_vec)
    for _ in iter0:maxiter

        # eponesweep! will be eponesweepT0! or eponesweep depending on alpha
        stat[:elapsed_eponesweep] = @elapsed begin
            errs = epmodel.updatealg!(epfields, epalg, epmat, stat)
        end

        max_err = stat[:max_err] = maximum(errs)

        # call back
        retflag, val = _parse_flag(oniter(iter, epmodel))
        retflag && return val

        # Converged
        max_err < epsconv && (returnstatus = CONVERGED_STATUS; break)

        if verbose 
            sweep_time = stat[:elapsed_eponesweep]
            inv_time = stat[:elapsed_eponesweep_inv]
            inv_frac = round(inv_time * 100/ sweep_time; digits = 3)

            update!(prog, max_err; showvalues = [
                (:iter, iter),
                (:maxiter, maxiter),
                (:alpha, alpha),
                (:max_beta, max_beta),
                (:sweep_time, sweep_time),
                (:inv_time, string(inv_time, " [", inv_frac, " %]")),
            ])
        end

        iter += 1
    end # for iter

    verbose && finish!(prog)

    # return
    # return produce_epout!(epmodel, returnstatus, iter; drop_epfields)
    return epmodel
end

## ------------------------------------------------------------------
function _parse_flag(ret::Tuple)
    length(ret) != 2 && return (false, nothing)
    e1, e2 = ret
    return (e1 isa Bool) ? (e1, e2) : (false, nothing)
end
_parse_flag(ret) = (false, nothing)
