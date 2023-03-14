# -------------------------------------------------------------------
# extras interface
import MetXBase.extras
extras(m::AbstractFluxEPModel) = m.extras

# -------------------------------------------------------------------
# state interface
@extras_dict_interface AbstractFluxEPModel state

# -------------------------------------------------------------------
# config interface
@extras_dict_interface AbstractFluxEPModel config
