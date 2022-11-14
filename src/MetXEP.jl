module MetXEP

    
    using Reexport
    
    @reexport using MetXBase
    @reexport using MetXOptim
    @reexport using MetXNetHub

    using SparseArrays
    using ProgressMeter
    using LinearAlgebra
    using LinearAlgebra: inv!
    using Statistics
    
    import Printf: @printf
    import ExtractMacro: @extract
    import SpecialFunctions: erf
    import Distributions: Truncated, Normal, mean, var, MvNormal

    #! include Types
    include("Types/AbstractFluxEPModel.jl")
    include("Types/FluxEPModelT0.jl")
        
    #! include Utils
    include("Utils/summary.jl")
    include("Utils/trunc_sample.jl")
    include("Utils/utils.jl")
    
    # #! include EPBase
    include("EPBase/compute_mom5d.jl")
    include("EPBase/matchmom.jl")
    include("EPBase/newav.jl")
    include("EPBase/newμs.jl")
    
    #! include AbstractFluxEPModelUtils
    include("AbstractFluxEPModelUtils/base.jl")
    include("AbstractFluxEPModelUtils/converge_ep!.jl")
    include("AbstractFluxEPModelUtils/interfaces.jl")
    include("AbstractFluxEPModelUtils/status.jl")

    #! include FluxEPModelT0
    include("FluxEPModelT0/eponesweep_T0.jl")
    include("FluxEPModelT0/interfaces.jl")

    #! include .

end