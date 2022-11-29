## ----------------------------------------------------------------------------------
_unbig(x; digits = 15) = round(Float64(x); digits)
function sample_tnorm(μ::Number, σ::Number, lb::Number, ub::Number; 
        xbins = 1000, digits = 15
    )
    
    local_max = clamp(μ, lb, ub)
    max_range = max(abs(local_max - lb), abs(local_max - ub))
    log_range = range(log10(max_range), log10(σ/xbins), length = div(xbins, 2))
    xs = [
        local_max .- 10.0.^(log_range);
        local_max;
        local_max .+ 10.0.^(reverse(log_range))
    ] 
    Txs = [big(x) for x in xs if lb <= x <= ub]
    Tpdf = ϕ.(Txs, big(μ), big(σ))
    Z = sum(@view(Tpdf[1:end - 1]) .* diff(Txs))
    Tpdf .= Tpdf ./ Z
    _unbig.(Txs; digits), _unbig.(Tpdf; digits)
end