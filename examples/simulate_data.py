# produce synthetic data for example 1

from os.path import dirname, join, pardir

from cobra.io import read_sbml_model

from dfba import DfbaModel, ExchangeFlux, KineticVariable, Parameter
from pypesto_dfba.optimize_dfba.objective_dfba import (ObjFunction,get_t_simu)
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
import pickle
import pypesto.sample as sample
from pypesto.store import (save_to_hdf5, read_from_hdf5)
from datetime import datetime

dir_to = "/home/erika/Documents/Projects/DFBA/results_example1/" \
         "simulated_data/100_starts_8_cores_L-BFGS-B_broaderBounds/"
model_dir = "/home/erika/Documents/Projects/DFBA/"\
            "dynamic-fba/sbml-models/iJR904.xml.gz"
# dfba_model, params_dict = get_dfba_model()
_, params_dict = get_dfba_model(model_dir)

# params_dict = {"K_g": 0.0027,
#           "v_gmax": 10.5,
#           "K_z": 0.0165,
#           "v_zmax": 6.0}
# get dfba-model in PickableDFBAModel-class
dfba_model = PicklableDFBAModel(model_dir, modifun)

# Simulate model
t_out = 0.1
t_end = 25
# dfba_model.solver_data.set_display("none")
concentrations, trajectories = dfba_model.simulate(
    0.0, t_end, t_out, ["EX_glc(e)", "EX_xyl_D(e)", "EX_etoh(e)"]
)
##
# simulate data
mu, sigma = 0, 0.25
n_obs = 3
# create dataframe with first column 'time' and subsequent columns observable
index = [0,10,30,50,80,150,200]
noise = np.zeros((n_obs, len(index)))
for i in range(n_obs):
    noise[i, :] = np.random.normal(mu, sigma, len(index))

observables = ["Biomass","Xylose","Glucose"]
data_time = concentrations["time"].iloc[index]
data_biomass = concentrations["Biomass"].iloc[index] + noise[0, :]
data_xylose = concentrations["Xylose"].iloc[index] + noise[1, :]
data_glucose = concentrations["Glucose"].iloc[index] + noise[2, :]

data = pd.concat([data_time,data_biomass,data_xylose,data_glucose],axis = 1)

now = datetime.now()
str_now = now.strftime(format="%Y-%m-%d %H:%M:%S").replace(' ','_')

data.to_csv(os.path.join(dir_to, "simulated_data_sigma_" + str(sigma) + ".csv"))

##
colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']

plt.figure()
for i_o in range(0, len(observables)):
    plt.plot(concentrations['time'], concentrations[observables[i_o]],
             label='simulation ' + observables[i_o], color=colors[i_o])
    plt.plot(data['time'], data[observables[i_o]], 'x', color=colors[i_o],
             label='data' + observables[i_o])
