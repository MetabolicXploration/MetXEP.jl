import Base.show
show(io::IO, m::AbstractFluxEPModel) = print(io, string(typeof(m)))
