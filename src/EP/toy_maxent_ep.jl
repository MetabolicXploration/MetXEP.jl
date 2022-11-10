# ## ----------------------------------------------------------------------
# function toy_maxent_ep(model;
#         alpha::Real=1e7,                   # inverse temperature
#         beta_vec = zeros(size(model, 2)),  # maxent inverse temperature vector
#         verbose::Bool=true,                # output verbosity
#         damp::Real=0.9,                    # damp ∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield
#         epsconv::Real=1e-6,                # convergence criterion
#         maxiter::Int=2000,                 # maximum iteration count
#         maxvar::Real=1e50,                 # maximum numerical variance
#         minvar::Real=1e-50,                # minimum numerical variance
#         solution = nothing,                # start from a solution
#         iter0 = 1
#     )

#     ## ----------------------------------------------------------------------
#     # Extracting
#     S, b, lb, ub = getfield.([model], [:S, :b, :lb, :ub])
#     α = alpha
#     β = beta_vec
#     M, N = size(S);

#     if isinf(α)
#         ## ----------------------------------------------------------------------
#         # EP T0 
#         error("ToyEPT0 Not implemented yet!!!")
#     else
#         ## ----------------------------------------------------------------------
#         # EP T > 0 

#         epfileds = isnothing(solution) ? 
#             EPFields(N, nothing, eltype(S)) : # ensure compat 
#             epfileds = deepcopy(solution.sol)

#         a = epfileds.a # approx priors mean
#         d = epfileds.d # approx priors variance
#         prog = ProgressThresh(epsconv; desc = "Toy_EP... ", dt = 0.5)
        
#         status = UNCONVERGED_STATUS
#         ## ----------------------------------------------------------------------
#         # Converge
#         iter = iter0
#         for outer iter in iter:maxiter
#             err = 0.0
            
#             for (n, rxn) in enumerate(model.rxns)

#                 # Compute Qn normal part parameters
#                 D = Diagonal(inv.(d))
#                 D[n,n] = 0.0
#                 Σ = inv(α*S'*S + D)
#                 v̄ = Σ*(α*S'*b + D*a + β)

#                 # tilted 
#                 Qn = Truncated(Normal(v̄[n], sqrt(Σ[n,n])), lb[n], ub[n])
#                 avQn = mean(Qn)
#                 vaQn = clamp(var(Qn), minvar, maxvar)

#                 # update a and b
#                 old_a, old_d = a[n], d[n] # prev
#                 new_d = clamp(inv(inv(vaQn) - inv(Σ[n,n])), minvar, maxvar)
#                 new_a = new_d*(avQn*(inv(new_d) + inv(Σ[n,n])) - (v̄[n]/Σ[n,n]))
#                 a[n] = damp * old_a + (1.0 - damp) * new_a
#                 d[n] = damp * old_d + (1.0 - damp) * new_d

#                 errn = max(abs(old_a - new_a), abs(old_d - new_d))
#                 err = max(err, errn)
                
#                 # update!(prog, err; showvalues = [
#                 #         ("iter", string(iter, "/", maxiter)),
#                 #         ("n", string(n, "/", N)),
#                 #         ("epsconv", epsconv),
#                 #         ("rxn", rxn),
#                 #         ("Σ[n,n]", Σ[n,n]),
#                 #         ("v̄[n]", v̄[n]),
#                 #         ("a[n]", a[n]),
#                 #         ("d[n]", d[n]),
#                 #         ("avQn", avQn),
#                 #         ("vaQn", vaQn),
#                 #     ]
#                 # )
#                 # sleep(0.1)
#             end

#             update!(prog, err; showvalues = [
#                         ("iter", string(iter, "/", maxiter)),
#                     ]
#                 )

#             err < epsconv && (status = CONVERGED_STATUS; break)
#         end # for it in 1:maxiter

#         finish!(prog)

#         # Marginals
#         μ = zeros(N)
#         σ = zeros(N)
#         av = zeros(N)
#         va = zeros(N)
#         for (n, rxn) in enumerate(model.rxns)

#             # Compute Qn normal part parameters
#             D = Diagonal(inv.(d))
#             D[n,n] = 0.0
#             Σ⁻¹ = α*S'*S + D
#             Σ = inv(Σ⁻¹)
#             v̄ = Σ*(α*S'*b + D*a)
#             μ[n] = v̄[n]
#             σ[n] = Σ[n,n]

#             # tilted 
#             Qn = Truncated(Normal(v̄[n], sqrt(Σ[n,n])), lb[n], ub[n])
#             av[n] = mean(Qn)
#             va[n] = var(Qn)

#         end

#         EPOut(μ, σ, av, va, epfileds, status, iter)
#     end
# end

