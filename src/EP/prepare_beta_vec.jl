# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function prepare_beta_vec(epmat::EPMat{T}, beta_vec::AbstractVector{T}) where T
    N = length(epmat.v)
    isempty(beta_vec) && return spzeros(T, N)
    beta_len = length(beta_vec)
    beta_len != N && error("Dim missmatch, beta_vec lenght != N, ($(beta_len)) != ($N)")
    return sparse(beta_vec)
end

function prepare_beta_vec(epmat::EPMatT0{T}, beta_vec::AbstractVector{T}) where T
    M = length(epmat.vd) # Number of dependent variables
    N = M + length(epmat.vi)
    isempty(beta_vec) && return spzeros(T, N)
    beta_len = length(beta_vec)
    beta_len != N && error("Dim missmatch, beta_vec lenght != N, ($(beta_len)) != ($N)")
    return sparse(beta_vec[epmat.idxmap]) # see echelonize
end
