# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

# TODO: Move to MetXBase
function inplaceinverse!(dest::AbstractArray, source::AbstractArray; δ = 1e-10)
    copyto!(dest, source)
    try
        inv!(cholesky!(Hermitian(dest)))
    catch err
        if err isa PosDefException
            nearPD!(dest, δ)
            # @show isposdef(dest)
            inv!(cholesky!(Hermitian(dest)))
            return 
        end
        rethrow(err)
    end
end

Φ(x) = 0.5*(1.0+erf(x/sqrt(2.0)))
ϕ(x) = exp(-x.^2/2.0)/sqrt(2π)
ϕ(x, μ, σ) = inv(σ*sqrt(2π))*exp(-((x - μ)^2)/(2σ^2))