# get cost function for different parameter values

from os.path import dirname, join, pardir
from cobra.io import read_sbml_model
from dfba import DfbaModel, ExchangeFlux, KineticVariable, Parameter
from pypesto_dfba.optimize_dfba.objective_dfba import (ObjFunction, get_t_simu)
from examples.get_dfba_model_ex1_ex6 import get_dfba_model, PicklableDFBAModel, modifun

import matplotlib
matplotlib.use('TkAgg')
import os
import pypesto
import pypesto.visualize as visualize
import numpy as np
import matplotlib.pyplot as plt
import pypesto.optimize as optimize
import pandas as pd
from pypesto_dfba.optimize_dfba.objective_dfba import (ObjFunction,
                                                       get_obs_names)
import pickle
import pypesto.sample as sample
from pypesto.store import (save_to_hdf5, read_from_hdf5)
from datetime import datetime
import time
##
# Example 1 - Real data:
model_dir = "/home/erika/Documents/Projects/DFBA/dynamic-fba/" \
                  "sbml-models/iJR904.xml.gz"
examplename = "example1_aerobic"
data_dir = "/home/erika/Documents/Projects/DFBA/results_example1/" \
                 "real_data/data_Fig1.csv"
cost_funct = "NLLH_normal" #|  "NLLH_normal" | "NLLH_laplace"
scaling_param_biomass = True

# dfba_model, params_dict = get_dfba_model(model_dir)
_, params_dict = get_dfba_model(model_dir, examplename)
# get dfba-model in PickableDFBAModel-class
dfba_model = PicklableDFBAModel(model_dir, modifun, examplename)
data = pd.read_csv(data_dir, index_col=0)

# initialize Objective Function Class
param_scale = 'lin'
# observable names
obs_names = get_obs_names(data)
# add sigmas
if "NLLH" in cost_funct:
    # define new sigma parameters for each observable
    for i_o in range(len(obs_names)):
        params_dict["sigma_" + obs_names[i_o]] = 1
# add scaling parameter for Biomass (scaling * simulated_biomass)
if scaling_param_biomass:
    params_dict["sc_biomass"] = 1

# parameter names
par_names = list(params_dict.keys())

obj_function = ObjFunction(dfba_model, data, par_names, param_scale,
                               cost_funct)
#{'K_g': 1.5473177020976163, 'v_gmax': 11.874787719859071,
# 'K_z': 0.0001, 'v_zmax': 6.718035754774346,
# 'sigma_Xylose': 0.21918679684314626, 'sigma_Glucose': 0.21877746850634938, 'sigma_Biomass': 0.18952320292283625,
# sc_biomass': 1.373112029096991}#
params = np.array([1.5473177020976163, 11.874787719859071,
          0.0001, 6.718035754774346,
          0.21918679684314626, 0.21877746850634938, 0.18952320292283625,
          1.373112029096991])

cost = obj_function(parameters=params)
