## ------------------------------------------------------------------------------
# A model for optimizing over the space of flux configurations
export FluxEPModelT0
Base.@kwdef struct FluxEPModelT0 <: AbstractFluxEPModel

    # T0
    Σd::Union{Nothing, Matrix{Float64}} = nothing
    Σi::Union{Nothing, Matrix{Float64}} = nothing
    G::Union{Nothing, Matrix{Float64}} = nothing
    vd::Union{Nothing, Vector{Float64}} = nothing
    vi::Union{Nothing, Vector{Float64}} = nothing
    be::Union{Nothing, Vector{Float64}} = nothing
    # unsorted permutation (see echelonize)
    idxmap::Union{Nothing, Vector{Int}}  = nothing

    # ep fields
    betai::Union{Nothing, Vector{Float64}} = nothing
    betad::Union{Nothing, Vector{Float64}} = nothing
    
    avi::Union{Nothing, Vector{Float64}} = nothing
    avd::Union{Nothing, Vector{Float64}} = nothing
    vai::Union{Nothing, Vector{Float64}} = nothing
    vad::Union{Nothing, Vector{Float64}} = nothing
    μi::Union{Nothing, Vector{Float64}} = nothing
    μd::Union{Nothing, Vector{Float64}} = nothing
    si::Union{Nothing, Vector{Float64}} = nothing
    sd::Union{Nothing, Vector{Float64}} = nothing
    ai::Union{Nothing, Vector{Float64}} = nothing
    ad::Union{Nothing, Vector{Float64}} = nothing
    di::Union{Nothing, Vector{Float64}} = nothing
    dd::Union{Nothing, Vector{Float64}} = nothing

    lbi::Union{Nothing, Vector{Float64}} = nothing
    lbd::Union{Nothing, Vector{Float64}} = nothing
    ubi::Union{Nothing, Vector{Float64}} = nothing
    ubd::Union{Nothing, Vector{Float64}} = nothing
    scalefact::Union{Nothing, Float64} = nothing

    siteflagave_i::Union{Nothing, BitArray{1}} = nothing
    siteflagave_d::Union{Nothing, BitArray{1}} = nothing
    siteflagvar_i::Union{Nothing, BitArray{1}} = nothing
    siteflagvar_d::Union{Nothing, BitArray{1}} = nothing

    # extras
    extras::Dict = Dict()

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


# A struct that contain all the data required for converging ep
function FluxEPModelT0(
        S::AbstractArray{T,2}, b::AbstractArray{T}, 
        lb::AbstractArray{T}, ub::AbstractArray{T};
        beta::AbstractVector{T} = spzeros(T, size(S, 2)),               # maxent inverse temperature vector
        solution = nothing,                                             # seed a solution
        expval = nothing,                                               # fix posterior probability experimental values 
    ) where {T<:AbstractFloat}

    # Some checks
    M, N = size(S)
    M > N && @warn("Number of rows M=$M larger than number of cols N=$N")
    any(lb .> ub) && error("lower bound fluxes > upper bound fluxes. Consider swapping lower and upper bounds")

    # echelonize
    S = S isa DenseMatrix ? S : Matrix(S)
    # idxmap ci is the inverse permutation that sends me back to the original model rxn order
    _, _, idxmap, IG, be = echelonize(S, b)
    Mech, Nech = size(IG)
    Ni, Nd = Mech, Nech - Mech
    Ii, Id = 1:Ni, (Ni+1):Nech
    length(be) == Ni || error("vector size incompatible with matrix") 
    Σi, Σd = zeros(T, Ni, Ni), zeros(T, Nd, Nd)
    G = IG[:, Id]
    vi, vd = zeros(T, Ni), zeros(T, Nd)
    lbi, lbd = lb[idxmap[Ii]], lb[idxmap[Id]]
    ubi, ubd = ub[idxmap[Ii]], ub[idxmap[Id]]

    # ep fields
    avi, avd = zeros(T, Ni), zeros(T, Nd)
    vai, vad = zeros(T, Ni), zeros(T, Nd)
    μi, μd = zeros(T, Ni), zeros(T, Nd)
    si, sd = ones(T, Ni), ones(T, Nd)
    ai, ad = zeros(T, Ni), zeros(T, Nd)
    di, dd = ones(T, Ni), ones(T, Nd)

    # expval
    # TODO: implement fixing var and ave (play with idxmap)
    siteflagvar_i, siteflagvar_d = trues(Ni), trues(Nd)
    siteflagave_i, siteflagave_d = trues(Ni), trues(Nd)
    
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
    rmul!(μi, scalefact_inv); rmul!(μd, scalefact_inv)
    rmul!(si, scalefact_inv^2); rmul!(sd, scalefact_inv^2)
    rmul!(avi, scalefact_inv); rmul!(avd, scalefact_inv)
    rmul!(vai, scalefact_inv^2); rmul!(vad, scalefact_inv^2)
    rmul!(lb, scalefact_inv); rmul!(ub, scalefact_inv);
    rmul!(be, scalefact_inv)

    # prepare beta
    betai = spzeros(T, Ni)
    betad = spzeros(T, Nd)
    if !isempty(beta) 
        betai .= beta[idxmap[Ii]]
        betad .= beta[idxmap[Id]]
    end

    return FluxEPModelT0(;Σd, Σi, G, vd, vi, be, idxmap, betai, betad, avi, avd, vai, vad, μi, μd, si, sd, ai, ad, di, dd, lbi, lbd, ubi, ubd, scalefact, siteflagave_i, siteflagave_d, siteflagvar_i, siteflagvar_d)
end


function FluxEPModelT0(net::MetNet; 
        netfields = [:rxns],                # fields to chache
        netcopy = false,                    # flag to make an internal copy of the net fields
        ep_model_kwargs...
    ) 
    
    epm = FluxEPModelT0(net.S, net.b, net.lb, net.ub; ep_model_kwargs...)

    # cache net
    net0 = extract_fields(net, netfields)
    net0 = netcopy ? deepcopy(net0) : net0
    net1 = MetNet(; net0...)
    metnet!(epm, reindex(net1; rxn_idxs = epm.idxmap))

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
