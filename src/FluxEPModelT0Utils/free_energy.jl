# TODO: add numerical stability stuff

_Ncdf(x) = 0.5*(1.0+erf(big(x)/sqrt(2.0)))
_Ncdf(x, av, sd) = _Ncdf((x - av) / sd)

function free_energy(epm::FluxEPModelT0{T}; 
        full = false
    ) where {T}

    Ni = length(epm.idxi)
    Nd = length(epm.idxd)
    μsi = _μ(epm, epm.idxi)
    μsd = _μ(epm, epm.idxd)
    ssi = _s(epm, epm.idxi)
    ssd = _s(epm, epm.idxd)
    lsi = lb(epm, epm.idxi)
    lsd = lb(epm, epm.idxd)
    usi = ub(epm, epm.idxi)
    usd = ub(epm, epm.idxd)
    Σi = _Σi(epm)

    # ZQ
    _ZQ = big(2*π)^(size(Σi, 1)/2)*sqrt(det(big.(Σi)))
    
    # ∑logZ_Qn
    _∑logZ_Qn = big(0.0)

    for i in 1:Ni
        # ZQn
        F_ul = _Ncdf(usi[i], μsi[i], sqrt(ssi[i])) - _Ncdf(lsi[i], μsi[i], sqrt(ssi[i]))
        _ZQn = _ZQ * F_ul
        # _ZQn = F_ul
        _∑logZ_Qn += log(_ZQn)
    end

    # ∑logZ_Qn
    full && for i in 1:Nd
        # ZQn
        F_ul = _Ncdf(
            usd[i], μsd[i], sqrt(ssd[i])) - _Ncdf(lsd[i], μsd[i], sqrt(ssd[i])
        )
        _ZQn = _ZQ * F_ul
        # _ZQn = F_ul
        _∑logZ_Qn += log(_ZQn)
    end
    
    _log_ZQ = log(_ZQ)
    
    D = full ? (Ni + Nd - 1) : (Ni - 1)
    F_EP = D * _log_ZQ - _∑logZ_Qn

    return convert(T, -F_EP), convert(T,  _log_ZQ), convert(T,  _∑logZ_Qn)
end
