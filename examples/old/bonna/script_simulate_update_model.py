from os.path import dirname, join, pardir

from cobra.io import read_sbml_model

from dfba import DfbaModel, ExchangeFlux, KineticVariable, Parameter
from pypesto_dfba.optimize_dfba.objective_dfba import (ObjFunction,get_t_simu)
from examples.get_dfba_model_ex1 import get_dfba_model

import matplotlib
matplotlib.use('TkAgg')
import os
import pypesto
import pypesto.visualize as visualize
import numpy as np
import matplotlib.pyplot as plt
import pypesto.optimize as optimize
import pandas as pd
import pickle
import pypesto.sample as sample
from pypesto.store import (save_to_hdf5, read_from_hdf5)
from datetime import datetime

# dir_to = "/home/erika/Documents/Projects/DFBA/results_example1/test2/"

# dfba_model, params_dict = get_dfba_model()

fba_model = read_sbml_model("/home/erika/Documents/Projects/DFBA/dynamic-fba/"
                                "sbml-models/iJR904.xml.gz")
# fba_model = read_sbml_model(join(dirname(__file__), pardir,
#                                  "sbml-models", "iJR904.xml.gz"))

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
params = {"K_g": 0.0027,
          "v_gmax": 10.5,
          "K_z": 0.0165,
          "v_zmax": 6.0}
K_g = Parameter(list(params.keys())[0],list(params.values())[0])
v_gmax = Parameter(list(params.keys())[1],list(params.values())[1])
K_z = Parameter(list(params.keys())[2], list(params.values())[2])
v_zmax = Parameter(list(params.keys())[3],list(params.values())[3])

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

##
# Simulate model
# t_out = 0.1
# t_end = 25
dfba_model.solver_data.set_display("none")
# concentrations, trajectories = dfba_model.simulate(
#     0.0, t_end, t_out, ["EX_glc(e)", "EX_xyl_D(e)", "EX_etoh(e)"]
# )

data = pd.read_csv("/home/erika/Documents/Projects/DFBA/results_example1/"
                   "simulated_data/test_8_starts_8_cores_no_maxls/"
                   "simulated_data_sigma_0_01_8starts_L-BFGS-B.csv", index_col=0)


def test_simulate(par_dict, title,save_name):
    dfba_model.update_parameters(par_dict)

    t_start = 0.0
    t_end = 20.0
    t_out = 0.1
    concentrations, trajectories = dfba_model.simulate(t_start, t_end+t_out,
                                                       t_out)
    print(concentrations)

    # print(1.9443830711793235 - concentrations['Biomass'][20])
    observables = data.columns
    colors = ['#fff7fb','#ece2f0','#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016c59','#014636']
    colors = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f']

    plt.figure()
    for i_o in range(1,len(observables)):
        plt.plot(concentrations['time'],concentrations[observables[i_o]],
                 label='simulation '+ observables[i_o], color = colors[i_o])
        plt.plot(data['time'],data[observables[i_o]], 'x',color = colors[i_o], label='data' + observables[i_o])

    plt.title(title)

    plt.legend()

    # plt.savefig(dir_to + 'simu_trajectories_' +save_name+ '.png')
p1 = [10**(-2.7855967),   10**(1.02277356), 10**(-3.20007978),  10**(0.79189569)]
p2 = [10**(-2.8940697),   10**(1.01469381), 10**(-2.13818811),  10**(0.97856977)]

# par_dict = {"K_g": 0.0027,"v_gmax": 10.5,"K_z": 0.0165, "v_zmax": 6.0 } #original
par_dict = {"K_g": p1[0],"v_gmax": p1[1],"K_z": p1[2], "v_zmax": p1[3] }
# title = 'original parameters before \n' + str(par_dict)
title = 'fval = ' + str(0.015692270)
save_name = '01_orig_params_before'
test_simulate(par_dict,title,save_name)

# par_dict = {'K_g': 0.0001, 'v_gmax': 100.0, 'K_z': 0.0004953259189234883, 'v_zmax': 99.7}
par_dict = {"K_g": p2[0],"v_gmax": p2[1],"K_z": p2[2], "v_zmax": p2[3] }
title = 'update parameters with \n' + str(par_dict)
title= 'fval = ' + str(0.4327313798)
save_name = '02_updated_params'
test_simulate(par_dict,title,save_name)

par_dict = {"K_g": 0.0027,"v_gmax": 10.5,"K_z": 0.0165, "v_zmax": 6.0 }
# title = 'original parameters after \n' + str(par_dict)
title = 'original parameters'
save_name = '03_original_after_updated_params'
test_simulate(par_dict,title,save_name)