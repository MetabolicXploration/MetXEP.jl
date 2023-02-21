_Ncdf(x) = 0.5*(1.0+erf(big(x)/sqrt(2.0)))
_Ncdf(x, av, sd) = _Ncdf((x - av) / sd)

export free_energy
function free_energy(epm::FluxEPModelT0{T}) where {T}

    N = length(epm.idxi)
    μs = MetXEP._μ(epm, epm.idxi)
    ss = MetXEP._s(epm, epm.idxi)
    ls = lb(epm, epm.idxi)
    us = ub(epm, epm.idxi)
    Σ = MetXEP._Σi(epm)
    
    # ZQ
    _ZQ = big(2*π)^(size(Σ, 1)/2)*sqrt(det(big.(Σ)))
    
    # ∑logZ_Qn
    _∑logZ_Qn = big(0.0)
    for i in 1:N
        # ZQn
        F_ul = _Ncdf(us[i], μs[i], sqrt(ss[i])) - _Ncdf(ls[i], μs[i], sqrt(ss[i]))
        # TODO: make this more stable
        _ZQn = _ZQ * F_ul
        _∑logZ_Qn += log(_ZQn)
    end
    
    _log_ZQ = log(_ZQ)
    
    F_EP = (N - 1)*_log_ZQ - _∑logZ_Qn

    return convert(T, F_EP), convert(T,  _log_ZQ), convert(T,  _∑logZ_Qn)
    
end
