# PLOT OPTIMIZATION RESULTS
# plots results from the setup defined in "pd_info_XXX.csv" file,
# which is located in dir_results
#
# from saved results.hdf5 files
# plot waterfall plot, trajectory plot, parameter plot
# sampling plots (for sampling plots you need hdf5-file and pickle results
# (pickle: containing the result.sampling)
# ####################################################################
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from pypesto.store import (save_to_hdf5, read_from_hdf5)
import pypesto.visualize as visualize
from examples.get_dfba_model_ex1_ex6 import get_dfba_model, PicklableDFBAModel, modifun
from pypesto_dfba.optimize_dfba.objective_dfba import (ObjFunction,get_t_simu)
import pickle
from dfba.plot.matplotlib import *
import matplotlib.pyplot as plt
import os, glob

# plots results from a "pd_info_XXX.csv" file, which is located in dir_results

dir_results = "/home/erika/Documents/Projects/DFBA/results_example6/" \
              "tests/test/"

dir_pd_info = glob.glob(os.path.join(dir_results, "pd_info*"))
if len(dir_pd_info)>1:
    raise Warning("Found more then 1 pd_info-files: " + dir_pd_info +
                  "Plot results for: " + dir_pd_info[0])
dir_pd_info_0 = dir_pd_info[0]

pd_info = pd.read_csv(dir_pd_info_0, index_col=0)
dir_to = pd_info.loc['dir_to'][0]
example_name = pd_info.loc['example_name'][0]
model_dir = pd_info.loc['model_dir'][0]
data_dir = pd_info.loc['data_dir'][0]
starts = pd_info.loc['nstarts'][0]
opt_method = pd_info.loc['opt_method'][0]
cost_func = pd_info.loc['cost_function'][0]


# folder = "SLSQP"
# starts = "100"
# nllh = True
# Example 1
# dir_to = "/home/erika/Documents/Projects/DFBA/results_example1/real_data/" \
#          + folder + "_" + starts + "/"
# dir_to = "/home/erika/Documents/Projects/DFBA/results_example1/" \
#          "real_data_laplace_noise/SLSQP_NLLH_laplace_100/"
# example_name = "example1_aerobic" # example1, example6
# model_dir = "/home/erika/Documents/Projects/DFBA/"\
#             "dynamic-fba/sbml-models/iJR904.xml.gz"     # example1
# data_dir = "/home/erika/Documents/Projects/DFBA/results_example1/" \
#                "simulated_data_sigma_0_01_25starts_L-BFGS-B.csv"
# data_dir = "/home/erika/Documents/Projects/DFBA/results_example1/" \
#            "simulated_data/simulated_data_sigma_0.25.csv"
# data_dir = "/home/erika/Documents/Projects/DFBA/results_example1/" \
#               "real_data/data_Fig1.csv"
# Example 6
# dir_to = "/home/erika/Documents/Projects/DFBA/results_example6/ex6_synthetic/" + \
#          folder + "_"+starts+"/"
# dir_to = "/home/erika/Documents/Projects/DFBA/results_example6/tests/Ausreiser/"
# example_name = "example6"
# model_dir = "/home/erika/Documents/Projects/DFBA/"\
#             "dynamic-fba/sbml-models/iND750.xml.gz"   # example6
# data_dir = "/home/erika/Documents/Projects/DFBA/results_example6/" \
#            "simulated_data/simulated_data_sigma_0.25_ex6_Ausreiser.csv"
# data_dir = "/home/erika/Documents/Projects/DFBA/results_example6/" \
#             "simulated_data/simulated_data_sigma_0.01_each_minute_ex6.csv"


res_path = "results_" + starts + "starts_" + opt_method + "_" + cost_func + "_.hdf5"
# res_path = "results_1starts_SLSQP_NLLH_laplace_.hdf5"

data = pd.read_csv(data_dir, index_col=0)
if example_name == "example6":
    plot_extra = ["Glucose", "Volume"]  # ex 6
elif example_name == "example1_aerobic":
    plot_extra = []  # "Ethanol", "Oxygen"]
    # data['error_Biomass'] = 0.014
else:
    plot_extra = ["Ethanol", "Oxygen"]


param_scale = 'log10'

name_to_save = res_path[:-5]
dfba_model, params_dict = get_dfba_model(model_dir, example_name)
par_names = list(params_dict.keys())
# dfba_model = PicklableDFBAModel(model_dir, modifun, example_name)

hdf5_reader = \
    read_from_hdf5.OptimizationResultHDF5Reader(dir_to + res_path)
result = hdf5_reader.read()
##
# plot waterfall
visualize.waterfall(result, size=(15, 6))
plt.savefig(dir_to + 'waterfall_' + name_to_save + '.png')
# plt.close()

# READ pickle optimization results
# result = pickle.load( open("/home/erika/Documents/Projects/DFBA/"
#                            "results_example1/5starts_opt_Kg_Kz/"
#                            "testresults_L-BFGS-B_5starts_opt_kg_kv.pickle", "rb" ) )
# pickle.dump(result, open( dir_to + 'results_' +name_to_save+ '.pickle', "wb" ) )
# df = pd.read_csv(dir_to + 'df_results_' + name_to_save+ ".csv", index_col=0)
#
# SIMULATE trajectories

result_nr = 0  # which opt. start should be simulated and plotted?
x_hat = result.optimize_result.list[result_nr]['x']

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

concentrations_best, trajectories_best = dfba_model.simulate(t_start, t_end,
                                                             t_out)#,
                                                           #["EX_glc__D_e",
                                                           # "EX_etoh_e",
                                                           # "EX_glyc_e"])
##
# i_t = 32
# delta_Gluc = trajectories_best.iloc[i_t, 1] * concentrations_best.iloc[i_t,1] + \
#     par_dict['D'] * (100- concentrations_best.iloc[i_t,3]) / \
#              concentrations_best.iloc[i_t,5]
# print(delta_Gluc)
# EX_glc__D_e * Biomass + d *(100-Glucose)/Volume
##
# plt.rcParams["figure.figsize"] = 9, 5
# plot_trajectories(trajectories_best)
# # plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
# plt.show()
# plt.savefig(dir_to + 'exchangeRates_' +
#             name_to_save + '_' + '_x' + str(result_nr) + '.png', bbox_inches='tight')
##
fs=12
# TRAJECTORIES PLOT
observables = data.columns
colors = ['#fff7fb', '#ece2f0', '#d0d1e6', '#a6bddb', '#67a9cf', '#3690c0',
          '#02818a', '#016c59', '#014636']
colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']

plt.figure(figsize=[7,5])
for i_o in range(1, len(observables)):
    plt.plot(concentrations_best['time'], concentrations_best[observables[i_o]],
             label='simulation ' + observables[i_o], color=colors[i_o])
    # if example_name == "example1_aerobic":
    #     plt.errorbar(data['time'], data[observables[i_o]],
    #                  yerr=data['error_Biomass'], fmt='x', color=colors[i_o],
    #                  label='data ' + observables[i_o])
    # else:
    plt.plot(data['time'], data[observables[i_o]], 'x', color=colors[i_o],
             label='data ' + observables[i_o])
if plot_extra:
    for i_e in range(len(plot_extra)):
        plt.plot(concentrations_best['time'],
                 concentrations_best[plot_extra[i_e]],
                 label='simulation ' + plot_extra[i_e], color=colors[i_o+1+i_e])
plt.title(str(x_hat))
plt.xlabel('Time [h]', fontsize=fs)
plt.ylabel('Concentration [g/L]', fontsize=fs)
if simu_original:
    plt.title('original parameters \n' + str(params_dict))

plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')

if simu_original:
    save_add = 'original'

plt.savefig(dir_to + 'simu_trajectories_' +
            name_to_save + '_' + '_x' + str(result_nr) + '.png',
            bbox_inches='tight')
# plt.close()


# PARAMETER PLOT
# create a reference point from it
if not 'NLLH' in cost_func and not example_name=="example1_aerobic":
# if not nllh and not example_name=="example1_aerobic":
    if example_name == "example1":
        ref = {'x': [np.log10(params_dict['K_g']), np.log10(params_dict['v_gmax']),
                     np.log10(params_dict['K_z']), np.log10(params_dict['v_zmax'])],
               'fval': result.optimize_result.as_dataframe().iloc[0, :]['fval'],
               'color': [0.2, 0.4, 1., 1.], 'legend': 'reference'}
    elif example_name == "example6":
        ref = {'x': [np.log10(params_dict['D']), np.log10(params_dict['Vgmax'])],
               'fval': result.optimize_result.as_dataframe().iloc[0, :]['fval'],
               'color': [0.2, 0.4, 1., 1.], 'legend': 'reference'}
    else:
        raise ValueError("Reference is not implemented yet for other then ex1 or "
                         "ex1_aerobic or ex6.")
    ax = visualize.parameters(result,
                              balance_alpha=False, size=[10, 6],
                              reference=ref)
elif example_name == "example1_aerobic":
    ax = visualize.parameters(result,
                              balance_alpha=False, size=[10, 6])
else:
    for i_o in range(1, len(observables)):
        par_names.append("sigma_" + observables[i_o])
    ax = visualize.parameters(result,
                              balance_alpha=False, size=[10, 6])


ax.set_yticklabels(par_names)

plt.savefig(dir_to + 'parameters_' +
            name_to_save + '_' + '_x' + str(result_nr) + '.png')

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

plt.savefig(dir_to + 'optimizer_history_' +
            name_to_save + '_' + '_x' + str(result_nr) + '.png')

param_numpy = np.array(list(par_dict.values()))
param_numpy = result.optimize_result.list[0]['x']

obj_function_test = ObjFunction(dfba_model, data, par_names, 'log10')
cost = obj_function_test(param_numpy)
print(cost)

##
# -------------------------------------------------------------------------
# ----------------SAMPLING ------------------------------------------------

sam_method = "AM"
n_samples = 10000
opt_method = "SLSQP"
n_multistart = "200"

example_name = "example1"

dir_to = "/home/erika/Documents/Projects/DFBA/results_example1/bonna/SLSQP_200/"

dir_result_object = dir_to +"results_after_sampling_" + str(n_samples) + \
                     "_" + sam_method + "_1ch_results_" + n_multistart + \
                     "starts_" + opt_method + "_.hdf5"

dir_sampling = dir_to + "sampling_" + str(n_samples) + "_" + sam_method + \
               "_1ch_results_" + n_multistart + "starts_" + opt_method + \
               "_.pickle"

hdf5_reader = \
    read_from_hdf5.OptimizationResultHDF5Reader(dir_result_object)
result1 = hdf5_reader.read()
result_sampling = pickle.load(open(dir_sampling, "rb"))

result1.sample_result = result_sampling
##
if example_name == "example1":
    params_dict = {"K_g": 0.0027,
               "v_gmax": 10.5,
               "K_z": 0.0165,
               "v_zmax": 6.0}
elif example_name == "example6":
    params_dict = {"D": 0.044,
               "Vgmax": 8.5}

import pypesto.sample as sample
sample.geweke_test(result=result1)


# ax = visualize.sampling_parameters_trace(result1, use_problem_bounds=True,
#                                          size=(12, 5))
# for i in range(len(params_dict)):
    # ax.figure.axes[i].set_ylabel(list(params_dict.keys())[i])

ax = visualize.sampling_parameters_trace(result1, use_problem_bounds=True,
                                         full_trace=True, size=(12,5))
for i in range(len(params_dict)):
    ax.figure.axes[i].set_ylabel(list(params_dict.keys())[i])

# ax.set_ylabel(list(params_dict.keys())[3])
plt.savefig(dir_to + 'sampling_'+sam_method+'_' + str(n_samples) + '.png',
            bbox_inches='tight')
plt.show()
# plt.close()
##

visualize.sampling_fval_trace(result1, size=(12, 5))
plt.savefig(dir_to+'sampling_fval_' + sam_method+'_' + str(n_samples) + '.png',
            bbox_inches='tight')
# plt.close()
##
ax = visualize.sampling_scatter(result1, size=[13, 6])
# ax = visualize.sampling_parameter_traces(res, use_problem_bounds=False, size=(12,5))
if example_name == "example6":
    ax.fig.axes[0].set_ylabel(list(params_dict.keys())[0])
    ax.fig.axes[2].set_ylabel(list(params_dict.keys())[1])
    ax.fig.axes[2].set_xlabel(list(params_dict.keys())[0])
    ax.fig.axes[3].set_xlabel(list(params_dict.keys())[1])
    ax.fig.axes[0].set_xlim(result1.problem.lb[0],result1.problem.ub[0])
    ax.fig.axes[2].set_xlim(result1.problem.lb[0], result1.problem.ub[0])
    ax.fig.axes[1].set_xlim(result1.problem.lb[1], result1.problem.ub[1])
    ax.fig.axes[3].set_xlim(result1.problem.lb[1], result1.problem.ub[1])
    ax.fig.axes[0].set_ylim(result1.problem.lb[0], result1.problem.ub[0])
    ax.fig.axes[1].set_ylim(result1.problem.lb[0], result1.problem.ub[0])
    ax.fig.axes[2].set_ylim(result1.problem.lb[1], result1.problem.ub[1])
    ax.fig.axes[3].set_ylim(result1.problem.lb[1], result1.problem.ub[1])
elif example_name == "example1":
    ax.fig.axes[0].set_ylabel(list(params_dict.keys())[0])
    ax.fig.axes[4].set_ylabel(list(params_dict.keys())[1])
    ax.fig.axes[8].set_ylabel(list(params_dict.keys())[2])
    ax.fig.axes[12].set_ylabel(list(params_dict.keys())[3])
    ax.fig.axes[12].set_xlabel(list(params_dict.keys())[0])
    ax.fig.axes[13].set_xlabel(list(params_dict.keys())[1])
    ax.fig.axes[14].set_xlabel(list(params_dict.keys())[2])
    ax.fig.axes[15].set_xlabel(list(params_dict.keys())[3])
    ax.fig.axes[12].set_xlim(result1.problem.lb[0]-0.5, result1.problem.ub[0]+0.5)
    ax.fig.axes[13].set_xlim(result1.problem.lb[1]-0.5, result1.problem.ub[1]+0.5)
    ax.fig.axes[14].set_xlim(result1.problem.lb[2]-0.5, result1.problem.ub[2]+0.5)
    ax.fig.axes[15].set_xlim(result1.problem.lb[3]-0.5, result1.problem.ub[3]+0.5)
    ax.fig.axes[0].set_ylim(result1.problem.lb[0], result1.problem.ub[0])
    ax.fig.axes[4].set_ylim(result1.problem.lb[1], result1.problem.ub[1])
    ax.fig.axes[8].set_ylim(result1.problem.lb[2], result1.problem.ub[2])
    ax.fig.axes[12].set_ylim(result1.problem.lb[3], result1.problem.ub[3])


plt.savefig(dir_to + 'sampling_scatter_'+sam_method+'_' + str(n_samples) +
            '.png', bbox_inches='tight')
# plt.close()
##
# read existing files
# n_samples=10000
# with open('/home/erika/Documents/Projects/DFBA/results_sampling_AM_' +str(n_samples) + '.pickle',
#           'rb') as result_file:
#     result = pickle.load(result_file)  # or the full result object
#     result_file.close()
##
ax = visualize.sampling_1d_marginals(result1)
ax[0][0].set_xlabel(list(params_dict.keys())[0])
ax[0][1].set_xlabel(list(params_dict.keys())[1])
ax[1][0].set_xlabel(list(params_dict.keys())[2])
ax[1][1].set_xlabel(list(params_dict.keys())[3])
# plt.savefig(dir_to+ 'sampling_1d_marginals_AM_' +str(n_samples) +'.png', bbox_inches = 'tight')

