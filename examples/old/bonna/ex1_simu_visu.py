
from os.path import dirname, join, pardir

from cobra.io import read_sbml_model

from dfba import DfbaModel, ExchangeFlux, KineticVariable, Parameter
from pypesto_dfba.optimize_dfba.objective_dfba import (ObjFunction,get_t_simu)
from examples.get_dfba_model_ex1 import get_dfba_model, PicklableDFBAModel, modifun

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

dir_to = "/home/erika/Documents/Projects/DFBA/results_example1/simulated_data/" \
         "middle_40starts_TNC/"

resultspath = dir_to + "results_40starts_TNC_2.hdf5"
name_to_save = str.split(str.split(resultspath,'/')[-1],'.')[0]# 'results_40starts_TNC_'

data = pd.read_csv("/home/erika/Documents/Projects/DFBA/results_example1/"
                   "simulated_data/"
                   "simulated_data_sigma_0_01_25starts_L-BFGS-B.csv", index_col=0)
dfba_model, params_dict = get_dfba_model()

# initialize Objective Function Class
# param_scale = 'lin'
param_scale = 'log10'

par_names = list(params_dict.keys())
#par_names = ["K_g","v_gmax","K_z","v_zmax"]
# par_names = ["K_g","K_z"]
obj_function = ObjFunction(dfba_model, data, par_names, param_scale)

# create objective object for pyPESTO
objective2 = pypesto.Objective(fun=obj_function, grad=False, hess=False)

save_add = '_log'




# resultspath = "/home/erika/Documents/Projects/DFBA/results_example1/simulated_data/tightBounds_25starts_L-BGFS-B/results_25starts_L-BFGS-B.hdf5"
hdf5_reader = read_from_hdf5.OptimizationResultHDF5Reader(resultspath)
result2 = hdf5_reader.read()
# problem = read_from_hdf5.OptimizationResultHDF5Reader(resultspath[:-5])
# problem = hdf5_reader.read()

# result = pickle.load( open("/home/erika/Documents/Projects/DFBA/results_example1/"
#                            "simulated_data/middle_40starts_TNC/"
#                            "results_40starts_TNC.pickle", "rb" ) )
# pickle.dump(result, open( dir_to + 'results_' +name_to_save+ '.pickle', "wb" ) )
# results_TCN = result
##
result2.problem = problem_test
pypesto_result_writer = save_to_hdf5.OptimizationResultHDF5Writer(
            dir_to + "results_40starts_TNC_2.hdf5")
pypesto_result_writer.write(result2, overwrite=True)
pypesto_problem_writer = save_to_hdf5.ProblemHDF5Writer(
    dir_to + "results_40starts_TNC_2.hdf5")
pypesto_problem_writer.write(problem_test, overwrite=True)
##
# result2.problem.lb = np.array([-4, -1, -4, -1])
# result2.problem.lb_full = [-4, -1, -4, -1]
# result2.problem.ub = [-0.5, 2, -0.5, 2]
# result2.problem.ub_full = [-0.5, 2, -0.5, 2]
##


visualize.waterfall(result2, size=(15,6))
plt.savefig(dir_to + 'waterfall_'+name_to_save+'_'+ save_add +'.png')
##
# SIMULATE trajectories
# x_hat = result.optimize_result.list[0].x

result_nr = 0  # which opt. start should be simulated and plotted?
# x_hat = df['x'][result_nr]
x_hat = result2.optimize_result.list[result_nr]['x']

# x_hat = result_previous.optimize_result.list[result_nr]['x']
par_dict = {}
if param_scale == 'log10':
    for i_p in range(len(par_names)):
        par_dict[par_names[i_p]] = 10**x_hat[i_p]
else:
    for i_p in range(len(par_names)):
        par_dict[par_names[i_p]] = x_hat[i_p]

simu_original = False
if simu_original:
    dfba_model.update_parameters(params_dict)
else:
    dfba_model.update_parameters(par_dict)
# dfba_model.update_parameters(par_dict_original)

t_start, t_end, t_out = get_t_simu(data)
t_out = 0.1
concentrations_best, trajectories_best = dfba_model.simulate(t_start, t_end, t_out)

##
# TRAJECTORIES PLOT
observables = data.columns
colors = ['#fff7fb','#ece2f0','#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016c59','#014636']
colors = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f']

plt.figure()
for i_o in range(1,len(observables)):
    plt.plot(concentrations_best['time'],concentrations_best[observables[i_o]],
             label='simulation '+ observables[i_o], color = colors[i_o])
    plt.plot(data['time'],data[observables[i_o]], 'x',color = colors[i_o], label='data' + observables[i_o])
plt.title(str(x_hat))
if simu_original:
    plt.title('original parameters \n' + str(params_dict))

plt.legend()
#
if simu_original:
    save_add = save_add + 'original'
#
plt.savefig(dir_to + 'simu_trajectories_' +
            name_to_save +'_' + save_add +'_x' + str(result_nr)+ '.png')
# plt.close()

##
# PARAMETER PLOT
# plt.figure()
# create a reference point from it

# ref = {'x': [np.log10(par_dict_original['K_g']),
#              np.log10(par_dict_original['v_gmax']),
#              np.log10(par_dict_original['K_z']),
#              np.log10(par_dict_original['v_zmax'])],
#        'fval': df['fval'][0],
#        'color': [
#            0.2, 0.4, 1., 1.], 'legend': 'reference'}

ref = {'x': [np.log10(params_dict['K_g']),np.log10(params_dict['v_gmax']),
             np.log10(params_dict['K_z']),np.log10(params_dict['v_zmax'])],
       'fval': result2.optimize_result.as_dataframe().iloc[0,:]['fval'], 'color': [
    0.2, 0.4, 1., 1.], 'legend': 'reference'}

# ref = {'x': [np.log10(par_dict_original['K_g']),
#              np.log10(par_dict_original['K_z'])],
#        'fval': result.optimize_result.as_dataframe().iloc[0,:]['fval'], 'color': [
#     0.2, 0.4, 1., 1.], 'legend': 'reference'}
ax = visualize.parameters(result2,
                     balance_alpha=False, size=[10,6],
                     reference=ref)
ax.set_yticklabels(par_names)
#
plt.savefig(dir_to+ 'parameters_' +
            name_to_save +'_' + save_add +'_x' + str(result_nr)+ '.png')
# visualize.ReferencePoint( x=[par_dict['K_g'],par_dict['v_gmax'],
#                              par_dict['K_z'],par_dict['v_zmax']],
#                           fval=None, color='g', legend=None)
# ax.scatter(x=np.log10(par_dict['K_g']),y=4,marker='o',color='g')
# plt.scatter(x=np.log10(par_dict['v_gmax']),y=3,marker='o',color='k')
# plt.scatter(x=np.log10(par_dict['K_z']),y=2,marker='o',color='k')
# plt.scatter(x=np.log10(par_dict['v_zmax']),y=1,marker='o',color='k')
##
# plot one optimizer history
# create a reference point from it
# ref = {'x': x_hat,
#        'fval': result2.optimize_result.get_for_key('fval')[0],
#        'color': [0.2, 0.4, 1., 1.], 'legend': 'optimum'}
# visualize.optimizer_history(result2,
#                             reference=ref)
# #
# plt.savefig(dir_to+ 'optimizer_history_' +
#             name_to_save +'_' + save_add +'_x' + str(result_nr)+ '.png')
