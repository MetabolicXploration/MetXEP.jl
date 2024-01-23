## ------------------------------------------------------------------------------
# A struct that contain all the data required for converging ep over a lep model

#TODO: rename FluxEPModelT0 -> EPT0Model
Base.@kwdef struct FluxEPModelT0{T} <: AbstractFluxEPModel where {T<:AbstractFloat}

    # T0
    Σd::Union{Nothing, Matrix{T}} = nothing
    Σi::Union{Nothing, Matrix{T}} = nothing
    G::Union{Nothing, Matrix{T}} = nothing
    vd::Union{Nothing, Vector{T}} = nothing
    vi::Union{Nothing, Vector{T}} = nothing
    be::Union{Nothing, Vector{T}} = nothing
    basis::Union{Nothing, Matrix{T}} = nothing
    # permutations (see echelonize)
    idxi::Union{Nothing, Vector{Int}}  = nothing # Map from epmi -> lepi
    idxd::Union{Nothing, Vector{Int}}  = nothing # Map from epmd -> lepd
    idxmap::Union{Nothing, Vector{Int}}  = nothing # Map from epm -> lep
    idxmap_inv::Union{Nothing, Vector{Int}}  = nothing # Map from lep -> epm

    # ep fields
    betai::Union{Nothing, Vector{T}} = nothing
    betad::Union{Nothing, Vector{T}} = nothing

    gammai::Union{Nothing, Vector{T}} = nothing
    gammad::Union{Nothing, Vector{T}} = nothing
    
    avi::Union{Nothing, Vector{T}} = nothing
    avd::Union{Nothing, Vector{T}} = nothing
    vai::Union{Nothing, Vector{T}} = nothing
    vad::Union{Nothing, Vector{T}} = nothing
    μi::Union{Nothing, Vector{T}} = nothing
    μd::Union{Nothing, Vector{T}} = nothing
    si::Union{Nothing, Vector{T}} = nothing
    sd::Union{Nothing, Vector{T}} = nothing
    ai::Union{Nothing, Vector{T}} = nothing
    ad::Union{Nothing, Vector{T}} = nothing
    di::Union{Nothing, Vector{T}} = nothing
    dd::Union{Nothing, Vector{T}} = nothing

    lbi::Union{Nothing, Vector{T}} = nothing
    lbd::Union{Nothing, Vector{T}} = nothing
    ubi::Union{Nothing, Vector{T}} = nothing
    ubd::Union{Nothing, Vector{T}} = nothing
    scalefact::Union{Nothing, T} = nothing

    siteflagave_i::Union{Nothing, BitArray{1}} = nothing
    siteflagave_d::Union{Nothing, BitArray{1}} = nothing
    siteflagvar_i::Union{Nothing, BitArray{1}} = nothing
    siteflagvar_d::Union{Nothing, BitArray{1}} = nothing

    # extras
    extras::Dict = Dict()

end

function FluxEPModelT0(
        S::AbstractArray{T,2}, b::AbstractArray{T}, 
        lb::AbstractArray{T}, ub::AbstractArray{T};
        beta::AbstractVector{T} = spzeros(T, size(S, 2)),               # maxent inverse temperature vector
        gamma::AbstractVector{T} = spzeros(T, size(S, 2)),              # maxent quadratic terms bias
        expval = nothing,                                               # fix posterior probability experimental values 
    ) where {T<:AbstractFloat}

    # Some checks
    M, N = size(S)
    M > N && @warn("Number of rows M=$M larger than number of cols N=$N")
    for (i, (l, u)) in enumerate(zip(lb, ub))
        l <= u && continue
        error("lower bound [$l] >= upper bound [$u] at index $i")
    end

    # echelonize
    S = S isa DenseMatrix ? S : Matrix(S)
    # idxmap ci is the inverse permutation that sends me back to the original model rxn order
    idxi, idxd, idxmap, G, be = echelonize(S, b)
    Mech = size(G, 1)
    Nech = sum(size(G))
    Nd, Ni = Mech, Nech - Mech
    length(be) == Nd || error("vector size incompatible with matrix") 
    Σd, Σi = zeros(T, Nd, Nd), zeros(T, Ni, Ni)
    vd, vi = zeros(T, Nd), zeros(T, Ni)
    lbd, lbi = lb[idxd], lb[idxi]
    ubd, ubi = ub[idxd], ub[idxi]

    # ep fields
    avd, avi = zeros(T, Nd), zeros(T, Ni)
    vad, vai = zeros(T, Nd), zeros(T, Ni)
    μd, μi = zeros(T, Nd), zeros(T, Ni)
    sd, si = ones(T, Nd), ones(T, Ni)
    ad, ai = zeros(T, Nd), zeros(T, Ni)
    # max_λ = (sqrt(eigmax(G * G')) + 0.5)
    max_λ = 1.0
    dd, di = ones(T, Nd) .* max_λ, ones(T, Ni) .* max_λ
    
    # expval
    # TODO: implement fixing var and ave (play with idxmap)
    siteflagvar_d, siteflagvar_i = trues(Nd), trues(Ni)
    siteflagave_d, siteflagave_i = trues(Nd), trues(Ni)
    
    # expave, expvar = parseexpval!(expval, siteflagave, siteflagvar)
    # for (k,v) in expave
    #     av[k] = v
    # end    
    # for (k,v) in expvar
    #     var[k] = v
    # end
    
    # Scale down μ, s, av, va, ub, lb and b using 
    # the maximum absolute bound (lb or ub).
    scalefact = max(maximum(abs, lb), maximum(abs, ub))
    scalefact_inv = inv(scalefact)
    rmul!(μd, scalefact_inv); rmul!(μi, scalefact_inv)
    rmul!(sd, scalefact_inv^2); rmul!(si, scalefact_inv^2)
    rmul!(avd, scalefact_inv); rmul!(avi, scalefact_inv)
    rmul!(vad, scalefact_inv^2); rmul!(vai, scalefact_inv^2)
    rmul!(lbd, scalefact_inv); rmul!(lbi, scalefact_inv);
    rmul!(ubd, scalefact_inv); rmul!(ubi, scalefact_inv);
    rmul!(be, scalefact_inv)

    # prepare beta
    betad, betai = spzeros(T, Nd), spzeros(T, Ni)
    if !isempty(beta) 
        betad .= beta[idxd]
        betai .= beta[idxi]
    end
    gammad, gammai = spzeros(T, Nd), spzeros(T, Ni)
    if !isempty(gamma) 
        gammad .= gamma[idxd]
        gammai .= gamma[idxi]
    end

    idxmap_inv = sortperm(idxmap)
    basis = basis_mat(G, idxi, idxd)

    # compute Σi and vi, the parameters (of the independent variables) of the full gaussian Q
    # TODO: check alloc performance
    Gt = G'
    Dd = Diagonal(inv.(dd))
    Di = Diagonal(inv.(di))
    gammaiD = Diagonal(gammai)
    gammadD = Diagonal(gammad)
    Σi_inv = Di - gammaiD + Gt * Dd * G - Gt * gammadD * G
    inplaceinverse!(Σi, Σi_inv)
    vi .= Σi * (Gt * Dd * (be - ad) - Gt * betad - Gt * gammadD * be + Di * ai + betai)
    
    # Compute Σd and vd, the parameters (of the dependent variables) of the full gaussian Q
    mul!(Σd, G * Σi, Gt) 
    mul!(vd, G, vi)
    vd .= be - vd
    
    return FluxEPModelT0{T}(;
        Σd, Σi, G, vd, vi, be, basis, 
        idxi, idxd, idxmap, idxmap_inv,
        betai, betad, 
        gammai, gammad, 
        avi, avd, vai, vad, 
        μi, μd, si, sd, 
        ai, ad, di, dd, 
        lbi, lbd, ubi, ubd, 
        scalefact, 
        siteflagave_i, siteflagave_d, 
        siteflagvar_i, siteflagvar_d
    )
end

# submodel
function FluxEPModelT0(template::FluxEPModelT0; to_overwrite...)
    dict = Dict{Symbol, Any}(to_overwrite)

    for field in fieldnames(typeof(template))
        haskey(dict, field) && continue # avoid use the template version
        dict[field] = getfield(template, field)
    end
    
    return FluxEPModelT0(;dict...)
end

## ------------------------------------------------------------------
function FluxEPModelT0(lep::LEPModel; 
        add_extras!::Function = (epm) -> nothing,
        ep_model_kwargs...
    )
    
    epm = FluxEPModelT0(lep.S, lep.b, lep.lb, lep.ub; ep_model_kwargs...)

    # cache lep specific stuff
    extras!(epm, :colids, colids(lep))
    extras!(epm, :rowids, rowids(lep))

    # custom estras
    add_extras!(epm)

    # default config
    config!(epm, :verbose, false)
    config!(epm, :damp, 0.9)
    config!(epm, :epsconv, 1e-6)
    config!(epm, :maxiter, 2000)
    config!(epm, :maxvar, 1e50)
    config!(epm, :minvar, 1e-50)
    config!(epm, :iter0, 1)
    
    # init state
    state!(epm, :status, UNSET_STATUS) 

    return epm
end

# Default for lep interfaced objects
FluxEPModelT0(model; kwargs...) = FluxEPModelT0(lepmodel(model); kwargs...)
