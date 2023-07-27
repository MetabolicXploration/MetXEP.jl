# TODO
# Make/documment common setting interface setter(model, ider, val)
# make a fixxing-like setting interface

module MetXEP
    
    using Reexport
    
    using MetXBase
    using MetXOptim
    
    using SparseArrays
    using ProgressMeter
    using LinearAlgebra
    using LinearAlgebra: inv!
    using Statistics
    
    import Printf: @printf
    import ExtractMacro: @extract
    import Distributions: truncated, Truncated, Normal, mean, var, MvNormal

    #! include Types
    include("Types/0_AbstractFluxEPModel.jl")
    include("Types/1_FluxEPModelT0.jl")
    
    #! include EPBase
    include("EPBase/compute_mom5d.jl")
    include("EPBase/epconverge.jl")
    include("EPBase/eponesweep_relaxed.jl")
    include("EPBase/get_join.jl")
    include("EPBase/get_scalefactor.jl")
    include("EPBase/matchmom.jl")
    include("EPBase/newav.jl")
    include("EPBase/newμs.jl")
    include("EPBase/parseexpval.jl")
    include("EPBase/prepareinput.jl")
    include("EPBase/scaleepfield.jl")
    
    #! include AbstractFluxEPModelUtils
    include("AbstractFluxEPModelUtils/base.jl")
    include("AbstractFluxEPModelUtils/converge_ep!.jl")
    include("AbstractFluxEPModelUtils/interfaces.jl")
    include("AbstractFluxEPModelUtils/lep_interface.jl")
    include("AbstractFluxEPModelUtils/status.jl")

    #! include FluxEPModelT0Utils
    include("FluxEPModelT0Utils/Distributions.jl")
    include("FluxEPModelT0Utils/average_gd.jl")
    include("FluxEPModelT0Utils/base.jl")
    include("FluxEPModelT0Utils/eponesweep_T0.jl")
    include("FluxEPModelT0Utils/free_energy.jl")
    include("FluxEPModelT0Utils/interfaces.jl")
    include("FluxEPModelT0Utils/lep_interface.jl")

    #! include .

    # exports
    @_exportall_non_underscore()

end