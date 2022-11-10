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
    include("Types/AbstractEPMat.jl")
    include("Types/AbstractEPModel.jl")
    include("Types/EPAlgs.jl")
    include("Types/EPFields.jl")
    include("Types/EPMat.jl")
    include("Types/EPMatT0.jl")
    include("Types/EPOut.jl")
    include("Types/FluxEPModel.jl")
        
    #! include Utils
    include("Utils/summary.jl")
    include("Utils/trunc_sample.jl")
    
    #! include EP
    include("EP/Q_sigma.jl")
    include("EP/compute_mom5d.jl")
    include("EP/converge_ep!.jl")
    include("EP/epconverge.jl")
    include("EP/eponesweep.jl")
    include("EP/eponesweepT0.jl")
    include("EP/fast_maxent_ep.jl")
    include("EP/get_join.jl")
    include("EP/get_scalefactor.jl")
    include("EP/matchmom.jl")
    include("EP/maxent_ep.jl")
    include("EP/newav.jl")
    include("EP/newμs.jl")
    include("EP/parseexpval.jl")
    include("EP/prepare_beta_vec.jl")
    include("EP/prepareinput.jl")
    include("EP/produce_epout.jl")
    include("EP/scaleepfield.jl")
    include("EP/status.jl")
    include("EP/toy_maxent_ep.jl")
    include("EP/utils.jl")
    
    #! include FluxEPModelUtils
    include("FluxEPModelUtils/base.jl")
    include("FluxEPModelUtils/interfaces.jl")

    #! include .

end