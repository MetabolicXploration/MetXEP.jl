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

# -------------------------------------------------------------------
# ep interface

# -------------------------------------------------------------------
# net interface
# NOTE: Do not interface with the net the data that will be in the model

# net data
@extras_val_interface AbstractFluxEPModel metnet MetNet

import MetXBase.metabolites
metabolites(m::AbstractFluxEPModel, ider...) = metabolites(metnet(m), ider...)

import MetXBase.reactions
reactions(m::AbstractFluxEPModel, ider...) = reactions(metnet(m), ider...)

import MetXBase.genes
genes(m::AbstractFluxEPModel, ider...) = genes(metnet(m), ider...)

# contraints data
import MetXBase.lb
import MetXBase.ub
import MetXBase.bounds