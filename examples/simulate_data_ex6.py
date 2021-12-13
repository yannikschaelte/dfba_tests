# produce synthetic data for example 2
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

dir_to = "/home/erika/Documents/Projects/DFBA/results_example6/" \
         "simulated_data/"
model_dir = "/home/erika/Documents/Projects/DFBA/dynamic-fba/" \
                  "sbml-models/iND750.xml.gz"


_, params_dict = get_dfba_model(model_dir,example_name="example6")

dfba_model = PicklableDFBAModel(model_dir, modifun, example_name="example6")

# Simulate model
t_out = 1.0
t_end = 16
# dfba_model.solver_data.set_display("none")

concentrations, trajectories = dfba_model.simulate(0.0, t_end, t_out)
##
# simulate data
mu, sigma = 0, 0.01
n_obs = 1
# create dataframe with first column 'time' and subsequent columns observable
index = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
noise = np.zeros((n_obs, len(index)))
for i in range(n_obs):
    noise[i, :] = np.random.normal(mu, sigma, len(index))

observables = ["Biomass"]
data_time = concentrations["time"].iloc[index]
data_biomass = concentrations["Biomass"].iloc[index] + noise[0, :]

data = pd.concat([data_time, data_biomass], axis=1)

now = datetime.now()
str_now = now.strftime(format="%Y-%m-%d %H:%M:%S").replace(' ', '_')

data.to_csv(os.path.join(dir_to, "simulated_data_sigma_" + str(sigma) +
                         "_each_minute_ex6.csv"))

##
colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']

plt.figure()
for i_o in range(0, len(observables)):
    plt.plot(concentrations['time'], concentrations[observables[i_o]],
             label='simulation ' + observables[i_o], color=colors[i_o])
    plt.plot(data['time'], data[observables[i_o]], 'x', color=colors[i_o],
             label='data' + observables[i_o])

plt.savefig(os.path.join(dir_to, "simulated_data_sigma_" + str(sigma) +
                         "_each_minute_ex6.png"))
