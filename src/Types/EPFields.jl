# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)
struct EPFields{T<:AbstractFloat}
    av::Vector{T}
    va::Vector{T}
    a::Vector{T}
    d::Vector{T}
    μ::Vector{T}
    s::Vector{T}

    siteflagave::BitArray{1}
    siteflagvar::BitArray{1}
end

function EPFields(N::Int, expval, T)
    
    siteflagvar = trues(N)
    siteflagave = trues(N)
    
    expave, expvar = parseexpval!(expval, siteflagave, siteflagvar)
    av = zeros(T,N)
    var = zeros(T,N)

    for (k,v) in expave
        av[k] = v
    end    
    for (k,v) in expvar
        var[k] = v
    end

    EPFields(
        #= av          =# av, 
        #= va          =# var, 
        #= a           =# zeros(T,N), 
        #= d           =# ones(T,N), 
        #= μ           =# zeros(T,N), 
        #= s           =# ones(T,N), 
        #= siteflagave =# siteflagave, 
        #= siteflagvar =# siteflagvar
    )
end

EPFields(N::Int, ::Nothing, T) = EPFields(
    #= av          =#  zeros(T,N), 
    #= va          =#  zeros(T,N), 
    #= a           =#  zeros(T,N), 
    #= d           =#  ones(T,N), 
    #= μ           =#  zeros(T,N), 
    #= s           =#  ones(T,N), 
    #= siteflagave =#  trues(N), 
    #= siteflagvar =#  trues(N)
)
