# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function inplaceinverse!(dest::AbstractArray, source::AbstractArray)
    copyto!(dest, source)
    inv!(cholesky!(Hermitian(dest)))
end

Φ(x) = 0.5*(1.0+erf(x/sqrt(2.0)))
ϕ(x) = exp(-x.^2/2.0)/sqrt(2π)
ϕ(x, μ, σ) = inv(σ*sqrt(2π))*exp(-((x - μ)^2)/(2σ^2))