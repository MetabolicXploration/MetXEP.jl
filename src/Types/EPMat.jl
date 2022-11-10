struct EPMat{T<:AbstractFloat} <: AbstractEPMat

    αStS::SparseMatrixCSC{T,Int64}  # A container of KK = βS'S (Its diag will change)
    αStS_diag::Vector{T}     # Store the original diag of KK = βS'S
    αStb::Vector{T}
    Σ::Matrix{T}                    # ΣQ Enforced to be dense 
    v::Vector{T}
    lb::Vector{T}
    ub::Vector{T}

end

function EPMat(
        S::AbstractArray{T, 2}, b::Vector{T}, 
        lb::Vector{T}, ub::Vector{T}, 
        alpha
    ) where T <: AbstractFloat

    alpha = convert(T, alpha)
    isinf(alpha) && error("For Inf alpha use 'EPMatT0'")
    M, N = size(S)
    αStS = sparse(alpha * S' * S)
    return EPMat(
        #= αStS      =# αStS, 
        #= αStS_diag =# diag(αStS), 
        #= αStb      =# alpha * S' * b, 
        #= Σ         =# zeros(T,N,N), 
        #= v         =# zeros(T,N), 
        #= lb        =# lb, 
        #= ub        =# ub
    )
    
end

EPMat(S, b, lb, ub, alpha) = EPMat(S, collect(b), collect(lb), collect(ub), float(alpha))s