module MetXEP

    
    using Reexport
    
    @reexport using MetXBase
    @reexport using MetXOptim
    @reexport using MetXNetHub

    using SparseArrays
    using ProgressMeter
    using LinearAlgebra
    using LinearAlgebra: inv!
    
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
    # include("EPBase/Q_sigma.jl")
    # include("EPBase/compute_mom5d.jl")
    # include("EPBase/converge_ep!.jl")
    # include("EPBase/epconverge.jl")
    # include("EPBase/eponesweep_T0.jl")
    # include("EPBase/eponesweep_relaxed.jl")
    # include("EPBase/fast_maxent_ep.jl")
    # include("EPBase/get_join.jl")
    # include("EPBase/get_scalefactor.jl")
    # include("EPBase/matchmom.jl")
    # include("EPBase/maxent_ep.jl")
    # include("EPBase/newav.jl")
    # include("EPBase/newμs.jl")
    # include("EPBase/prepare_beta_vec.jl")
    # include("EPBase/prepareinput.jl")
    # include("EPBase/produce_epout.jl")
    # include("EPBase/scaleepfield.jl")
    # include("EPBase/toy_maxent_ep.jl")
    # include("EPBase/utils.jl")
    
    #! include AbstractFluxEPModelUtils
    include("AbstractFluxEPModelUtils/base.jl")
    include("AbstractFluxEPModelUtils/converge_ep!.jl")
    include("AbstractFluxEPModelUtils/interfaces.jl")
    include("AbstractFluxEPModelUtils/status.jl")

    #! include FluxEPModelT0
    include("FluxEPModelT0/eponesweep_T0.jl")

    #! include .

end