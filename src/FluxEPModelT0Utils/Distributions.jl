import Distributions.MultivariateNormal
function MultivariateNormal(epm::FluxEPModelT0)
    vi = _dense(epm.vi)
    Σi = _dense(epm.Σi)
    return MultivariateNormal(
        vi * epm.scalefact,
        Σi * epm.scalefact^2
    )
end

import Distributions.entropy
function _normal_entropy(Σ)
    L = _cholesky(Σ).L
    N = size(Σ,1)
    S = sum(log.(diag(L))) .+ 0.5 * N * log(2 * pi * exp(1))
    return S
end

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