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

dir_to = "/home/erika/Documents/Projects/DFBA/results_example1/test_parallel/"

dfba_model, params_dict = get_dfba_model()

# params_dict = {"K_g": 0.0027,
#           "v_gmax": 10.5,
#           "K_z": 0.0165,
#           "v_zmax": 6.0}
#
# erikas_dfba_model = PicklableDFBAModel("/home/erika/Documents/Projects/DFBA/dynamic-fba/sbml-models/iJR904.xml.gz", modifun)
# dfba_model = erikas_dfba_model.dfba_model
##
# Simulate model
t_out = 0.1
t_end = 25
dfba_model.solver_data.set_display("none")
concentrations, trajectories = dfba_model.simulate(
    0.0, t_end, t_out, ["EX_glc(e)", "EX_xyl_D(e)", "EX_etoh(e)"]
)

##
# generate plots of results (in this case using plotlly)
# from dfba.plot.plotly import *
#
# import plotly.io as pio
#
# pio.renderers.default = "browser"
# # in plotly version 4, default version is different than in 3
# pio.templates.default = "plotly_white"
# fig = plot_concentrations(concentrations)
# fig.show()
# fig = plot_trajectories(trajectories)
# fig.show()
# write results to file
#concentrations.to_csv("concentrations.csv")
#trajectories.to_csv("trajectories.csv")
##
# DATA
# load measurement data (from Eiteman et al., (2008), Fig 2a)
# path_data = '/home/erika/Documents/Projects/DFBA/data_example1/data_Fig2a_extracted_apprtime.csv'
# data = pd.read_csv(path_data)
# data: first column has to be time!
##
# simulate data
# mu, sigma = 0, 0.02
# n_obs = 3
# # create dataframe with first column 'time' and subsequent columns observable
# index = [0,10,30,50,80,150,200]
# noise = np.zeros((n_obs,len(index)))
# for i in range(n_obs):
#     noise[i,:] = np.random.normal(mu, sigma, len(index))
#
# data_time = concentrations["time"].iloc[index]
# data_biomass = concentrations["Biomass"].iloc[index] + noise[0,:]
# data_xylose =  concentrations["Xylose"].iloc[index] + noise[1,:]
# data_glucose =  concentrations["Glucose"].iloc[index] + noise[2,:]
#
# data = pd.concat([data_time,data_biomass,data_xylose,data_glucose],axis = 1)
#
# now = datetime.now()
# str_now = now.strftime(format = "%Y-%m-%d %H:%M:%S").replace(' ','_')
#
# data.to_csv(os.path.join(dir_to, "simulated_data_sigma_0_02_" + str_now + ".csv"))
#
##
# data = pd.read_csv("/home/erika/Documents/Projects/DFBA/results_example1/"
#                    "simulated_data/tightBounds_50starts_02_L-BFGS-B_and_sampling/"
#                    "simulated_data_sigma_0_01_50starts_L-BFGS-B.csv", index_col=0)
data = pd.read_csv("/home/erika/Documents/Projects/DFBA/results_example1/"
                   "simulated_data/tightBounds_50starts_02_L-BFGS-B_and_sampling/"
                   "simulated_data_sigma_0_01_50starts_L-BFGS-B.csv", index_col=0)

##
# import matplotlib.pyplot as plt
# # plot trajectories
# plt.figure()
# plt.plot(concentrations['time'],concentrations['Biomass'], label='Biomass')
# plt.plot(concentrations['time'],concentrations['Glucose'], label='Glucose')
# plt.plot(concentrations['time'],concentrations['Volume'], label='Volume')
# #plt.plot(concentrations['time'],measured, 'kx', label='measured')
# plt.legend()

# par_dict = {'K_g': 0.1, 'v_gmax': 13.41647121146278, 'K_z': 0.0013965913120499189, 'v_zmax': 10.507141443062393}
# # par_dict = {'K_g': 0.06490681437002449, 'v_gmax': 10.679789201848868, 'K_z': 0.011207397140323462, 'v_zmax': 2.391019700112453}
# # par_dict = {'K_g': 0.0001, 'v_gmax': 100.0, 'K_z': 0.011213266277786694, 'v_zmax': 2.7433609437779714} #IDA solve error [IDA ERROR]  IDASolve   At t = 1.55147, , mxstep steps taken before reaching tout.

#
# def test_simulate(par_dict, title,save_name):
#     dfba_model.update_parameters(par_dict)
#
#     t_start = 0.0
#     t_end = 20.0
#     t_out = 0.1
#     concentrations, trajectories = dfba_model.simulate(t_start, t_end+t_out,
#                                                        t_out)
#     print(concentrations)
#
#     # print(1.9443830711793235 - concentrations['Biomass'][20])
#     observables = data.columns
#     colors = ['#fff7fb','#ece2f0','#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016c59','#014636']
#     colors = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f']
#
#     plt.figure()
#     for i_o in range(1,len(observables)):
#         plt.plot(concentrations['time'],concentrations[observables[i_o]],
#                  label='simulation '+ observables[i_o], color = colors[i_o])
#         plt.plot(data['time'],data[observables[i_o]], 'x',color = colors[i_o], label='data' + observables[i_o])
#
#     plt.title(title)
#
#     plt.legend()
#
#     plt.savefig(dir_to + 'simu_trajectories_' +save_name+ '.png')
#
# par_dict = {"K_g": 0.0027,"v_gmax": 10.5,"K_z": 0.0165, "v_zmax": 6.0 } #original
# # par_dict = {'K_g': 0.0001, 'v_gmax': 100.0, 'K_z': 0.0004953259189234883, 'v_zmax': 100.0}
# # par_dict = {'K_g': 0.09535644985619214, 'v_gmax': 99.9999, 'K_z': 0.0008813778843277387, 'v_zmax': 99.9999}
# title = 'original parameters before \n' + str(par_dict)
# save_name = '01_orig_params_before'
# test_simulate(par_dict,title,save_name)
#
# par_dict = {'K_g': 0.0001, 'v_gmax': 100.0, 'K_z': 0.0004953259189234883, 'v_zmax': 100.0}
# title = 'update parameters with \n' + str(par_dict)
# save_name = '02_updated_params'
# test_simulate(par_dict,title,save_name)
#
# par_dict = {"K_g": 0.0027,"v_gmax": 10.5,"K_z": 0.0165, "v_zmax": 6.0 }
# title = 'original parameters after \n' + str(par_dict)
# save_name = '03_original_after_updated_params'
# test_simulate(par_dict,title,save_name)

##
# initialize Objective Function Class
# param_scale = 'lin'
param_scale = 'log10'

par_names = list(params_dict.keys())
#par_names = ["K_g","v_gmax","K_z","v_zmax"]
# par_names = ["K_g","K_z"]
obj_function = ObjFunction(dfba_model, data, par_names, param_scale)

# create objective object for pyPESTO
objective2 = pypesto.Objective(fun=obj_function, grad=False, hess=False)


import cloudpickle as pickle
dump = pickle.dumps(objective2)
pickle.loads(dump) #drauf aufrufen, mit

##
# OPTIMIZATION
# define lower and upper bound for parameters optimization
# lb = 0 * np.ones((dim_full, 1)) #log
# ub = 10 * np.ones((dim_full, 1))
# par_dict_original = {}
# par_dict_original['K_g'] = 0.0027
# par_dict_original['v_gmax'] = 10.5
# par_dict_original['K_z'] = 0.0165
# par_dict_original['v_zmax'] = 6.0

do_optimize = False

if param_scale == 'lin':
    lb = [0,5]
    ub = [1,13]
    x_sc = ['lin','lin']
    x_g = np.array([[0.05, 8.5]])
elif param_scale == 'log10':
    lb = [-4,-0.5,-4,-0.5]
    ub =[-1,2,-1,2]
    lb=[np.log10(params_dict['K_g'])-1,
        np.log10(params_dict['v_gmax']) - 1,
        np.log10(params_dict['K_z']) - 1,
        np.log10(params_dict['v_zmax'])-1]
    #
    ub = [np.log10(params_dict['K_g'])+1,
          np.log10(params_dict['v_gmax']) + 1,
          np.log10(params_dict['K_z']) + 1,
          np.log10(params_dict['v_zmax'])+1]
    # K_g = Parameter("K_g", 0.0027)
    # v_gmax = Parameter("v_gmax", 10.5)
    # K_z = Parameter("K_z", 0.0165)
    # v_zmax = Parameter("v_zmax", 6.0)
    # ub = [np.log10(10*0.5),np.log10(10*8.5)]
    x_sc = ['log10','log10','log10','log10']
    # hdf5_reader = read_from_hdf5.OptimizationResultHDF5Reader(
    #     "/home/erika/Documents/"
    #     "Projects/DFBA/results_example1/"
    #     "tightBounds_25starts_L-BGFS-B/results_25starts_L-BFGS-B.hdf5")
    # result_previous = hdf5_reader.read()
    #x000 = result_previous.optimize_result.list[0]['x0'] #-2.91552379,  1.19795398, -1.99072438,  0.76946407
    #x000 =   [[-2.49875968, 1.10079146 , -2.32763403 , 1.46871346]] #fval:0.232185780918403	fval0:206.458530563814	nfval:820

problem1 = pypesto.Problem(objective=objective2, lb=lb, ub=ub,
                           copy_objective=False, x_scales=x_sc)#,
                          # x_guesses=[[x000]])


if do_optimize:
    opt_method = ['L-BFGS-B']#'Pyswarm']#,'TNC']#],'L-BFGS-B','Fides']
        # Fides: objective function, if no hessian_update is provided, this function must return a tuple (fval, grad), otherwise this function must return a tuple (fval, grad, Hessian)
    for i_o in range(len(opt_method)):
        # opt_method = 'L-BFGS-B'
        if opt_method[i_o]=='Pyswarm':
            my_optimizer = optimize.PyswarmOptimizer(options={'swarmsize':30, 'minstep': 1e-9,'minfunc' : 1e-9})
            # nstarts = 1
        elif opt_method[i_o]=='Fides':
            my_optimizer = optimize.FidesOptimizer()
            # my_optimizer.minimize(problem1,x0 = np.array([0,0,0,0]), id=1, allow_failed_starts=10)
        else:
            my_optimizer = optimize.ScipyOptimizer(method=opt_method[i_o],
                                                   options={'eps': 1e-09,
                                                            'maxls':40})  #  'eta':0.5})  # basic optimizer'L-BFGS-B'

        nstarts = 5

        # my_optimizer = pypesto.optimize.DlibOptimizer()

        print('----- starting optimization...............')
        # record the history
        name_to_save = str(nstarts) + 'starts_' + opt_method[i_o]
        now = datetime.now()
        str_now = now.strftime(format="%Y-%m-%d %H:%M:%S").replace(' ', '_')
        file1 = open(dir_to + 'OptStart_' +str_now+'_' +name_to_save+ '.txt', "a")
        file1.close()

        history_options = pypesto.HistoryOptions(trace_record=True)
        result = optimize.minimize(problem1, optimizer=my_optimizer, n_starts=nstarts,
                                   history_options=history_options)

        # store result
        now = datetime.now()
        str_now = now.strftime(format="%Y-%m-%d %H:%M:%S").replace(' ', '_')
        file2 = open(
            dir_to + 'OptEnd_' + str_now + '_' + name_to_save + '.txt', "a")
        file2.close()
        with open(dir_to + 'results_' +name_to_save+ '.pickle',
                'wb') as result_file:
                pickle.dump(result.optimize_result, result_file)  # or the full result object
                result_file.close()
        # save optimization result as hdf5 file
        store_filename = dir_to + 'results_' +name_to_save+ '.hdf5'
        pypesto_result_writer = save_to_hdf5.OptimizationResultHDF5Writer(
            store_filename)
        pypesto_result_writer.write(result, overwrite=True)
        pypesto_problem_writer = save_to_hdf5.ProblemHDF5Writer(
            store_filename)
        pypesto_problem_writer.write(problem1, overwrite=True)

        df = result.optimize_result.as_dataframe(
            ['fval','fval0','n_fval', 'x','x0','grad','n_grad', 'n_hess', 'n_res', 'n_sres', 'time'])
        df['lb'] = str(lb)
        df['ub'] = str(ub)

        df.to_csv(dir_to + 'df_results_' + name_to_save+ ".csv")



data.to_csv(os.path.join(dir_to, "simulated_data_sigma_0_01_" + name_to_save + ".csv"))
##
if param_scale == 'lin':
    save_add = ''
elif param_scale == 'log10':
    save_add = '_log'

## plot waterfall
visualize.waterfall(result, size=(15,6))
plt.savefig(dir_to + 'waterfall_'+name_to_save+'_'+ save_add +'.png')
# plt.close()
##
# READ optimization results
# result = pickle.load( open("/home/erika/Documents/Projects/DFBA/"
#                            "results_example1/5starts_opt_Kg_Kz/"
#                            "testresults_L-BFGS-B_5starts_opt_kg_kv.pickle", "rb" ) )
# pickle.dump(result, open( dir_to + 'results_' +name_to_save+ '.pickle', "wb" ) )
# results_TCN = result
name_to_save = str(nstarts)+ 'starts_' + opt_method[0]
store_filename = dir_to + 'results_' +name_to_save+ '.hdf5'
hdf5_reader = read_from_hdf5.OptimizationResultHDF5Reader(dir_to
                                                          + 'results_' +name_to_save+ '.hdf5')

df = pd.read_csv(dir_to + 'df_results_' + name_to_save+ ".csv", index_col=0)
##
resultspath = "/home/erika/Documents/Projects/DFBA/results_example1/results_50starts_L-BFGS-B.hdf5"

# resultspath = "/home/erika/Documents/Projects/DFBA/results_example1/simulated_data/tightBounds_25starts_L-BGFS-B/results_25starts_L-BFGS-B.hdf5"
hdf5_reader = read_from_hdf5.OptimizationResultHDF5Reader(resultspath)
result = hdf5_reader.read()
##
# SIMULATE trajectories
# x_hat = result.optimize_result.list[0].x

result_nr = 0  # which opt. start should be simulated and plotted?
# x_hat = df['x'][result_nr]
x_hat = result.optimize_result.list[result_nr]['x']

# x_hat = result_previous.optimize_result.list[result_nr]['x']
par_dict = {}
if param_scale == 'log10':
    for i_p in range(len(par_names)):
        par_dict[par_names[i_p]] = 10**x_hat[i_p]
else:
    for i_p in range(len(par_names)):
        par_dict[par_names[i_p]] = x_hat[i_p]

simu_original = True
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
##
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
       'fval': result.optimize_result.as_dataframe().iloc[0,:]['fval'], 'color': [
    0.2, 0.4, 1., 1.], 'legend': 'reference'}

# ref = {'x': [np.log10(par_dict_original['K_g']),
#              np.log10(par_dict_original['K_z'])],
#        'fval': result.optimize_result.as_dataframe().iloc[0,:]['fval'], 'color': [
#     0.2, 0.4, 1., 1.], 'legend': 'reference'}
ax = visualize.parameters(result,
                     balance_alpha=False, size=[10,6],
                     reference=ref)
ax.set_yticklabels(par_names)
##
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
ref = {'x': x_hat,
       'fval': result.optimize_result.get_for_key('fval')[0],
       'color': [0.2, 0.4, 1., 1.], 'legend': 'optimum'}
visualize.optimizer_history(result,
                            reference=ref)
#
plt.savefig(dir_to+ 'optimizer_history_' +
            name_to_save +'_' + save_add +'_x' + str(result_nr)+ '.png')

##
param_numpy = np.array(list(par_dict.values()))
param_numpy = result.optimize_result.list[0]['x']

obj_function_test = ObjFunction(dfba_model, data, par_names, 'log10')
cost = obj_function_test(param_numpy)
print(cost)

# par_dict

# #########################################################################
## SAMPLING --------------------------------------------

hdf5_reader = read_from_hdf5.OptimizationResultHDF5Reader("/home/erika/Documents/"
                                                          "Projects/DFBA/"
                                                          "results_example1/"
                                                          "results_50starts_L-BFGS-B.hdf5")
result = hdf5_reader.read()
x_hat = result.optimize_result.list[0]['x']
name_to_save = '50starts_L-BFGS-B'
## start 10:26-11:06 - 1000 samples
# ---------------------------------------------------------------------------
# sampler = sample.AdaptiveParallelTemperingSampler(
#     internal_sampler=sample.AdaptiveMetropolisSampler(),
#     n_chains=3)   #-> deepcopy problem
sampler = sample.AdaptiveMetropolisSampler()
# sampler = sample.ParallelTemperingSampler(internal_sampler=sample.MetropolisSampler(),
#                                           betas=[1, 1e-1, 1e-2])  # -> deepcopy problem
#
# result = sample.sample(problem1, n_samples=100, sampler=sampler, result=result)
n_samples=10000
result = sample.sample(problem1, n_samples=n_samples, sampler=sampler,
                       x0 = x_hat)
                       #x0 = np.array([np.log10(0.5),  np.log10(8.5)]))
##
with open(dir_to+ 'sampling_' + str(n_samples) + '_'+ name_to_save+ '.pickle' ,
          'wb') as result_file:
    pickle.dump(result.sample_result,
                result_file)  # or the full result object
    result_file.close()
store_filename2 = dir_to + 'results_and_sampling_' +str(n_samples) +'_'+\
                  name_to_save+ '.hdf5'
pypesto_result_writer = save_to_hdf5.OptimizationResultHDF5Writer(
    store_filename2)
pypesto_result_writer.write(result, overwrite=True)
pypesto_problem_writer = save_to_hdf5.ProblemHDF5Writer(
    store_filename2)
pypesto_problem_writer.write(problem1, overwrite=True)
##
# plt.figure()
visualize.sampling_parameters_trace(result, use_problem_bounds=True, size=(12,5))
plt.savefig(dir_to + 'sampling_AM_' +str(n_samples) +'.png', bbox_inches = 'tight')
plt.show()
# plt.close()
##

visualize.sampling_fval_trace(result, size=(12,5))
plt.savefig(dir_to+'sampling_fval_AM_' +str(n_samples) +'.png', bbox_inches = 'tight')
# plt.close()

visualize.sampling_scatter(result, size=[13,6])
plt.savefig(dir_to+ 'sampling_scatter_AM_' +str(n_samples) +'.png', bbox_inches = 'tight')
# plt.close()
##
# read existing files
n_samples=10000
with open('/home/erika/Documents/Projects/DFBA/results_sampling_AM_' +str(n_samples) + '.pickle',
          'rb') as result_file:
    result = pickle.load(result_file)  # or the full result object
    result_file.close()
##
visualize.sampling_1d_marginals(result)
plt.savefig(dir_to+ 'sampling_1d_marginals_AM_' +str(n_samples) +'.png', bbox_inches = 'tight')

