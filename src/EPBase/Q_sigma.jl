# """
#     (This method is not efficient for intence computations)
#     In Braunstein et al paper (see README)
#     KK = αF'F
#     Σ^-1 = (KK + D) 
#     return Σ, covariance matrix of Q, the full multivariate Gaussian
# """
# function Q_sigma(model::MetNet, epout::EPOut, α)

#     # in epout sol d is scale using maxflux, we'll scale it back
#     maxflux = max(maximum(abs.(lb)), maximum(abs.(ub)))
#     a = epout.sol.a * maxflux;
#     d = epout.sol.b * maxflux^2;
#     #= a, d are the mean and variance of the univariate Gaussians
#     used as site approximations in the EP =#
#     S = model.S
#     D = Matrix(Diagonal(1.0 ./ d))
#     invΣ = α*S'*S + D
#     Σ = inv(invΣ)
# end