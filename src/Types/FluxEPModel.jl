## ------------------------------------------------------------------------------

# A model for optimizing over the space of flux configurations
export FluxEPModel
struct FluxEPModel{T} <: AbstractEPMat

    # EP data
    scalefact::T
    updatealg!::Function
    epfields::EPFields{T}
    epmat::AbstractEPMat
    alpha::T
    beta_vec::SparseVector{T, Int}
    
    # extras
    extras::Dict{Any, Any}
end

# A struct that contain all the data required for converging ep
function FluxEPModel(
        S::AbstractArray{T,2}, b::AbstractArray{T}, 
        lb::AbstractArray{T}, ub::AbstractArray{T};
        alpha::Real=Inf,                                                # inverse temperature
        beta_vec::AbstractVector{T} = spzeros(T, size(S, 2)),           # maxent inverse temperature vector
        solution::Union{FluxEPModel{T}, EPFields, Nothing} = nothing,   # seed a solution
        expval = nothing,                                               # fix posterior probability experimental values 
    ) where {T<:Real}

    # Some checks
    M, N = size(S)
    M > N && @warn("M = $M ≥ N = $N")
    any(lb .> ub) && error("lower bound fluxes > upper bound fluxes. Consider swapping lower and upper bounds")

    # The scalefactor is just the maximum absolute bound (lb or ub).
    scalefact = get_scalefactor(lb, ub)

    # Create EPFields. If a solution is not given, the EPfields will be fresh
    epfields::EPFields = isnothing(solution) ? EPFields(N, expval, eltype(S)) : 
        deepcopy(solution.sol) # preserve the original solution!

    # making a local copy to rescale
    lb, ub, b = copy.([lb, ub, b]) 

    #=
    Scale down μ, s, av, va of epfields and ub, lb and Y using the 
    previous computed scalefactor.
    If epfields is fresh, it only will have effect on ub, lb and Y
    =#
    scaleepfield!(inv(scalefact), epfields, ub, lb, b) # scaling fields to [0,1]

    epmat = (alpha < Inf) ? EPMat(S, b, lb, ub, alpha) : EPMatT0(S, b, lb, ub)

    # One iteration of EP
    updatealg! = isinf(alpha) ? eponesweepT0! : eponesweep!

    beta_vec = prepare_beta_vec(epmat, beta_vec)

    epm = FluxEPModel{T}(scalefact, updatealg!, epfields, epmat, alpha, beta_vec, Dict())

    return epm

end


function FluxEPModel(net::MetNet; 
        netfields = [:rxns],                # fields to chache
        netcopy = false,                    # flag to make an internal copy of the net fields
        ep_model_kwargs...
    ) 
    
    epm = FluxEPModel(net.S, net.b, net.lb, net.ub; ep_model_kwargs...)

    # cache net
    net0 = extract_fields(net, netfields)
    net0 = netcopy ? deepcopy(net0) : net0
    net1 = MetNet(; net0...)
    metnet!(epm, net1)

    return epm
end
