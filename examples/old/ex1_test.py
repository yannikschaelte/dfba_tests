
from os.path import dirname, join, pardir

from cobra.io import read_sbml_model

from dfba import DfbaModel, ExchangeFlux, KineticVariable, Parameter

import matplotlib
matplotlib.use('TkAgg')

import pypesto
import pypesto.visualize as visualize
import numpy as np
# import scipy as sp
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import pypesto.optimize as optimize
import pandas as pd
import pickle
import pypesto.sample as sample


# define fixed parameters and instances of Parameters (to be used for estimation)
# the latter are added to model below
K_g = Parameter("K_g",0.0027)
v_gmax = Parameter("v_gmax",10.5)
K_z = Parameter("K_z",0.0165)
v_zmax = Parameter("v_zmax",6.0)
# {'K_g': 0.00027000000000000006, 'v_gmax': 105.00000000000004, 'K_z': 0.0016500000201385228, 'v_zmax': 60.0}   #t=1.8

# K_g = Parameter("K_g",0.0007499604954039659)
# v_gmax = Parameter("v_gmax",16.469289388385086)
# K_z = Parameter("K_z",0.0034085177897754857)
# v_zmax = Parameter("v_zmax",3.1255393118546695)
#
# K_g = Parameter("K_g",0.01170126028125197)
# v_gmax = Parameter("v_gmax",17.284851050108053)
# K_z = Parameter("K_z",0.0021365123322649674)
# v_zmax = Parameter("v_zmax",19.830764142515115)
#
# K_g = Parameter("K_g",0.006229950332686774)
# v_gmax = Parameter("v_gmax",25.06664906281686)
# K_z = Parameter("K_z",0.017546481661625115)
# v_zmax = Parameter("v_zmax",4.712842647909046)
#

# DfbaModel instance initialized with cobra model
fba_model = read_sbml_model(
    join(dirname(__file__), pardir, "sbml-models", "iJR904.xml.gz")
)


fba_model.solver = "glpk"
dfba_model = DfbaModel(fba_model)

# instances of KineticVariable (default initial conditions are 0.0, but can be
# set here if wanted e.g. Oxygen)
X = KineticVariable("Biomass")
Gluc = KineticVariable("Glucose")
Xyl = KineticVariable("Xylose")
Oxy = KineticVariable("Oxygen", initial_condition=0.24)
Eth = KineticVariable("Ethanol")
#
# add kinetic variables to dfba_model
dfba_model.add_kinetic_variables([X, Gluc, Xyl, Oxy, Eth])

# instances of ExchangeFlux
mu = ExchangeFlux("BiomassEcoli")
v_G = ExchangeFlux("EX_glc(e)")
v_Z = ExchangeFlux("EX_xyl_D(e)")
v_O = ExchangeFlux("EX_o2(e)")
v_E = ExchangeFlux("EX_etoh(e)")

# add exchange fluxes to dfba_model
dfba_model.add_exchange_fluxes([mu, v_G, v_Z, v_O, v_E])

# Here add parameters to dfba_model
dfba_model.add_parameters([K_g,v_gmax,K_z,v_zmax])

# add rhs expressions for kinetic variables in dfba_model
dfba_model.add_rhs_expression("Biomass", mu * X)
dfba_model.add_rhs_expression("Glucose", v_G * 180.1559 * X / 1000.0)
dfba_model.add_rhs_expression("Xylose", v_Z * 150.13 * X / 1000.0)
dfba_model.add_rhs_expression("Oxygen", v_O * 16.0 * X / 1000.0)
dfba_model.add_rhs_expression("Ethanol", v_E * 46.06844 * X / 1000.0)

# add lower/upper bound expressions for exchange fluxes in dfba_model together
# with expression that must be non-negative for correct evaluation of bounds
# v_gmax = 10.5
# K_g = 0.0027
dfba_model.add_exchange_flux_lb(
    "EX_glc(e)", v_gmax * (Gluc / (K_g + Gluc)) * (1 / (1 + Eth / 20.0)), Gluc
) #v_g glucose
dfba_model.add_exchange_flux_lb("EX_o2(e)", 15.0 * (Oxy / (0.024 + Oxy)), Oxy) # v_o, oxygen
# v_zmax = 6.0
# K_z = 0.0165
dfba_model.add_exchange_flux_lb(
    "EX_xyl_D(e)",
    v_zmax * (Xyl / (K_z + Xyl)) * (1 / (1 + Eth / 20.0)) * (1 / (1 + Gluc / 0.005)),
    Xyl,
) # v_z, xylose

# add initial conditions for kinetic variables in dfba_model biomass (gDW/L),
# metabolites (g/L)
dfba_model.add_initial_conditions(
    {
        "Biomass": 0.03,
        "Glucose": 15.5,
        "Xylose": 8.0,
        "Oxygen": 0.0,
        "Ethanol": 0.0,
    }
)

#
# dfba_model.solver_data.set_display("none")
t_out = 0.1
t_end = 25
concentrations, trajectories = dfba_model.simulate(
    0.0, t_end, t_out, ["EX_glc(e)", "EX_xyl_D(e)", "EX_etoh(e)"]
)