# # Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

# function epconverge!(epfields::EPFields, epmat::M, epalg::EPAlg, eponesweep!::T) where {T<:Function,M<:AbstractEPMat}
    
#     @extract epalg : maxiter verbose epsconv alpha beta_vec

#     returnstatus = :unconverged
#     iter = 0
    
#     # sweep ep till maxiter is reached or max(errav, errvar) < epsconv
#     # alpha = epalg.alpha epsconv maxiter
#     prog = ProgressThresh(epsconv)# ; desc = "MaxEnt EP... ")
#     max_beta = findmax(beta_vec)
#     while iter < maxiter
#         iter += 1
#         # eponesweep! will be eponesweepT0! or eponesweep depending on alpha
#         (errav, errvar, errμ, errs) = eponesweep!(epfields, epalg, epmat)

#         max_err = max(errav, errvar)
#         if max_err < epsconv
#             returnstatus = :converged
#             break
#         end

#         verbose && update!(prog, max_err; showvalues = [
#             (:iter, iter),
#             (:maxiter, maxiter),
#             (:alpha, alpha),
#             (:max_beta, max_beta)
#         ])
#     end

#     verbose && (finish!(prog); flush(stderr))

#     return returnstatus, iter
# end