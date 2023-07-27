# The non truncated join normal distribution
import Distributions.MultivariateNormal
function MultivariateNormal(epm::FluxEPModelT0)
    return MultivariateNormal(_vi(epm), _Σi(epm))
end

function _normal_entropy(Σ)
    L = _cholesky(Σ).L
    N = size(Σ,1)
    S = sum(log.(diag(L))) .+ 0.5 * N * log(2 * pi * exp(1))
    return S
end

import Distributions.entropy
function entropy(epm::FluxEPModelT0)
    # TODO: check why S(Cov = basis * epm.Σi * basis') != S(Cov = epm.Σi)
    # basis * vi = v
    # basis = epm.basis
    # Cov = basis * epm.Σi * basis'
    Cov = _Σi(epm)
    return _normal_entropy(Cov)
end

function orthnorm_entropy(epm::FluxEPModelT0)
    # basis * vi = v -> basis
    basis = epm.basis
    Cov = basis * epm.Σi * basis'
    ort_basis = mgrscho(basis)
    normΣ = ort_basis' * Cov * ort_basis
    
    return _normal_entropy(normΣ)
end

# returns the truncated marginal of rxn
function trunc_marginal(epm::FluxEPModelT0, rxn::Int)
    μ = _v(epm, rxn)
    s = _s(epm, rxn)
    l, u = lb(epm, rxn), ub(epm, rxn)
    return truncated(Normal(μ, sqrt(s)), l, u)
end

trunc_marginal(epm::FluxEPModelT0, rxn) = 
    trunc_marginal(epm, colindex(epm, rxn))

import MetXBase.sample_tnorm
function sample_tnorm(epm::FluxEPModelT0, rxn; 
        xbins = 1000, digits = 15
    )
    μ = _v(epm, rxn)
    s = _s(epm, rxn)
    l, u = lb(epm, rxn), ub(epm, rxn)
    xs, ws = MetXBase.sample_tnorm(
        μ, sqrt(s), l, u; 
        xbins, digits
    )
    return xs, ws
end