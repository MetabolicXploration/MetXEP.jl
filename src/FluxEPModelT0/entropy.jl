import Distributions.entropy
export entropy
function entropy(epm::FluxEPModelT0)
    
    # basis * vi = v
    basis = epm.basis
    Cov = basis * epm.Σi * basis'
    
    N = size(Cov,1);
    # this is just the entropy of a multivariate gaussian
    L = cholesky(Cov).L;
    S = sum(log.(diag(L))) + 0.5*N*log.(2*pi*exp(1));
    return S 
end

export orthnorm_entropy
function orthnorm_entropy(epm::FluxEPModelT0)
    basis = epm.basis
    # basis * vi = v -> basis
    ort_basis = mgrscho(basis)
    Cov = basis * epm.Σi * basis'
    normΣ = ort_basis' * Cov * ort_basis
    normΣ .= 0.5 * (normΣ + normΣ') # Correct Symetry
    
    N = size(normΣ,1);
    # this is just the entropy of a multivariate gaussian
    L = cholesky(normΣ).L;
    S = sum(log.(diag(L))) + 0.5*N*log.(2*pi*exp(1));
    return S 
end