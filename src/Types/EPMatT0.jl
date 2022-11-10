# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

struct EPMatT0{T<:AbstractFloat} <: AbstractEPMat
    Σd::Matrix{T}
    Σi::Matrix{T}
    G::Matrix{T}
    lb::Vector{T}
    ub::Vector{T}
    vd::Vector{T}
    vi::Vector{T}
    Y::Vector{T}
    idxmap::Vector{Int} # unsorted permutation (see echelonize)
end

function EPMatT0(K::AbstractArray{T,2}, Y::Vector{T}, 
        lb::Vector{T}, ub::Vector{T}) where T <: AbstractFloat
    M, N = size(K)    
    M > N && @warn("numeber of rows M=$M larger than number of cols N=$N")
    K = K isa DenseMatrix ? K : Matrix(K)
    # idxmap ci is the inverse permutation that sends me back to the original model rxn order
    _, _, idxmap, EK, EY = echelonize(K,Y)
    Mech, Nech = size(EK)
    length(EY) == Mech || error("vector size incompatible with matrix") 
    return EPMatT0( 
        #=Σd=#  zeros(T, Mech, Mech), 
        #=Σi=#  zeros(T, Nech-Mech, Nech-Mech), 
        #=G=#   copy(EK[1:Mech, Mech+1:Nech]), 
        #=lb=#  lb[idxmap],
        #=ub=#  ub[idxmap],
        #=vd=#  zeros(T, Mech),
        #=vi=#  zeros(T, Nech-Mech),
        #=Y=#   EY,
        #=idx=# idxmap
    ) 
end

EPMatT0(K, Y, lb, ub) = EPMatT0(K, collect(Y), collect(lb), collect(ub))