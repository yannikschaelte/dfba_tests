# Copyright (C) 2018, 2019 Columbia University Irving Medical Center,
#     New York, USA
# Copyright (C) 2019 Novo Nordisk Foundation Center for Biosustainability,
#     Technical University of Denmark

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


"""Anaerobic growth of *E. coli* on glucose and xylose.

Organism -> _Escherichia coli str. K-12 substr. MG1655_
Model in http://bigg.ucsd.edu/models/iJR904
"""

from os.path import dirname, join, pardir

from cobra.io import read_sbml_model

from dfba import DfbaModel, ExchangeFlux, KineticVariable, Parameter
from pypesto_dfba.optimize_dfba.objective_dfba import ObjFunction



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
from pypesto.store import (save_to_hdf5, read_from_hdf5)

# DfbaModel instance initialized with cobra model
# fba_model = read_sbml_model(
#     join(dirname(__file__), pardir, "sbml-models", "iJR904.xml.gz")
# )
# define fixed parameters and instances of Parameters (to be used for estimation)
# the latter are added to model below
K_g = Parameter("K_g",0.0027)
K_z = Parameter("K_z",0.0165)
v_gmax = Parameter("v_gmax",10.5)
v_zmax = Parameter("v_zmax",6.0)
# v_gmax = 10.5
# v_zmax =6.0

# K_g = Parameter("K_g",0.005)
# v_gmax = Parameter("v_gmax",9.5)
# K_z = Parameter("K_z",0.001)
# v_zmax = Parameter("v_zmax",7.0)

# init_biomass = Parameter("init_biomass",0.1)
# init_glucose = Parameter("init_glucose",16)
# init_xylose = Parameter("init_xylose",8.0)
# init_oxygen = Parameter("init_oxygen",0.0)
# init_ethanol = Parameter("init_ethanol",0.0)
# "Biomass": 0.1,
# "Glucose": 16,
# "Xylose": 8.0,
# "Oxygen": 0.0,
# "Ethanol": 0.0,
# K_g = 0.0027
# v_gmax = 10.5
# K_z = 0.0165
# v_zmax = 6.0

fba_model = read_sbml_model("/home/erika/Documents/Projects/DFBA/dynamic-fba/sbml-models/iJR904.xml.gz")
##
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
# dfba_model.add_parameters([K_g,K_z])#,init_biomass,init_glucose,
                          # init_xylose,init_oxygen,init_ethanol])


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
# "Biomass": 7.09720,
# "Glucose": 16.828,
# "Xylose": 0,
# "Oxygen": 0.0,
# "Ethanol": 0.0,
# "Biomass": 0.03,
# "Glucose": 15.5,
# "Xylose": 8.0,
# "Oxygen": 0.0,
# "Ethanol": 0.0,
# "Biomass": 0.1,
# "Glucose": 16,
# "Xylose": 8.0,
# "Oxygen": 0.0,
# "Ethanol": 0.0,


# simulate model across interval t = [0.0,25.0](hours) with outputs for plotting
# every 0.1h and optional list of fluxes
# concentrations, trajectories = dfba_model.simulate(
#     0.0, 25.0, 0.1, ["EX_glc(e)", "EX_xyl_D(e)", "EX_etoh(e)"]
# )
##
# dfba_model.solver_data.set_display("none")
t_out = 0.1
t_end = 25
dfba_model.solver_data.set_display("none")
concentrations, trajectories = dfba_model.simulate(
    0.0, t_end, t_out, ["EX_glc(e)", "EX_xyl_D(e)", "EX_etoh(e)"]
)



# simulate on a grid
# K_g = Parameter("K_g",0.0027)
# v_gmax = Parameter("v_gmax",10.5)
# K_z = Parameter("K_z",0.0165)
# v_zmax = Parameter("v_zmax",6.0)
# range_v_gmax = range(0,20)
# range_v_zmax = range(0,20)
# range_Kg = np.arange(1,10,1)
# range_Kz = np.arange(1,10,1)
# pd_v_var = pd.DataFrame(np.zeros([len(range_Kg),len(range_Kz)]),
#                         columns=range_Kg,index=range_Kz)
#
#
# for ikg,kg in enumerate(range_Kg):
#     for ikz,kz in enumerate(range_Kz):
#
#         par_dict = {"K_g": float(kg),
#                     "v_gmax": 10.5,
#                     "K_z": float(kz),
#                     "v_zmax": 6.0}
#         dfba_model.update_parameters(par_dict)
#         concentrations, trajectories = dfba_model.simulate(
#             0.0, t_end, t_out, ["EX_glc(e)", "EX_xyl_D(e)", "EX_etoh(e)"]
#         )
#         print(max(concentrations['time']))
#         print(kg)
#         print(kz)
#         pd_v_var.iloc[int(ikg),int(ikz)] = max(concentrations['time'])
#         print(pd_v_var)
#
# pd_v_var.to_csv("/home/erika/Documents/Projects/DFBA/results_example1/grid_kg_rows_kz_cols_1_10"
#                 ".csv")
#
# ##
# vg = 3
# vz = 1
# par_dict = {"K_g":0.0027,
#             "v_gmax":vg,
#             "K_z":0.0165,
#             "v_zmax":vz}
# print('vg ' +str(vg))
# print('vz: '+str(vz))
# dfba_model.update_parameters(par_dict)
# concentrations, trajectories = dfba_model.simulate(
#     0.0, t_end, t_out, ["EX_glc(e)", "EX_xyl_D(e)", "EX_etoh(e)"])
# print(max(concentrations['time']))
