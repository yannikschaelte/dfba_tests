# PLOT OPTIMIZATION RESULTS
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
##
folder = "SLSQP"
dir_to = "/home/erika/Documents/Projects/DFBA/results_example1/real_data/" \
         + folder + "_200/"
res_path = "results_200starts_" + folder + "_.hdf5"
model_dir = "/home/erika/Documents/Projects/DFBA/"\
            "dynamic-fba/sbml-models/iJR904.xml.gz"
# data_dir = "/home/erika/Documents/Projects/DFBA/results_example1/" \
#                "simulated_data_sigma_0_01_25starts_L-BFGS-B.csv"
# data_dir = "/home/erika/Documents/Projects/DFBA/results_example1/" \
#            "simulated_data/simulated_data_sigma_0.25.csv"
data_dir = "/home/erika/Documents/Projects/DFBA/results_example1/real_data/" \
           "data_Fig1.csv"

param_scale = 'log10'

name_to_save = res_path[:-5]
_, params_dict = get_dfba_model(model_dir)
par_names = list(params_dict.keys())
dfba_model = PicklableDFBAModel(model_dir, modifun)
data = pd.read_csv(data_dir, index_col=0)

hdf5_reader = \
    read_from_hdf5.OptimizationResultHDF5Reader(dir_to + res_path)
result = hdf5_reader.read()

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
                                                             t_out)

# TRAJECTORIES PLOT
observables = data.columns
colors = ['#fff7fb', '#ece2f0', '#d0d1e6', '#a6bddb', '#67a9cf', '#3690c0',
          '#02818a', '#016c59', '#014636']
colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']

plt.figure()
for i_o in range(1, len(observables)):
    plt.plot(concentrations_best['time'], concentrations_best[observables[i_o]],
             label='simulation ' + observables[i_o], color=colors[i_o])
    plt.plot(data['time'], data[observables[i_o]], 'x', color=colors[i_o],
             label='data' + observables[i_o])
plt.title(str(x_hat))
plt.xlabel('Time [h]')
plt.ylabel('Concentration [g/L]')
if simu_original:
    plt.title('original parameters \n' + str(params_dict))

plt.legend()

if simu_original:
    save_add = 'original'

plt.savefig(dir_to + 'simu_trajectories_' +
            name_to_save + '_' + '_x' + str(result_nr) + '.png')
# plt.close()


# PARAMETER PLOT
# create a reference point from it
ref = {'x': [np.log10(params_dict['K_g']), np.log10(params_dict['v_gmax']),
             np.log10(params_dict['K_z']), np.log10(params_dict['v_zmax'])],
       'fval': result.optimize_result.as_dataframe().iloc[0, :]['fval'],
       'color': [0.2, 0.4, 1., 1.], 'legend': 'reference'}

ax = visualize.parameters(result,
                          balance_alpha=False, size=[10, 6],
                          reference=ref)
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

dir_to = "/home/erika/Documents/Projects/DFBA/results_example1/" \
                    "bonna/SLSQP_200/"

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
params_dict = {"K_g": 0.0027,
               "v_gmax": 10.5,
               "K_z": 0.0165,
               "v_zmax": 6.0}

ax = visualize.sampling_parameters_trace(result1, use_problem_bounds=True,
                                         size=(12,5))
ax.set_ylabel(list(params_dict.keys())[3])
plt.savefig(dir_to + 'sampling_'+sam_method+'_' + str(n_samples) + '.png',
            bbox_inches='tight')
plt.show()
# plt.close()
##

visualize.sampling_fval_trace(result1, size=(12, 5))
plt.savefig(dir_to+'sampling_fval_' + sam_method+'_' + str(n_samples) + '.png',
            bbox_inches='tight')
# plt.close()

visualize.sampling_scatter(result1, size=[13, 6])
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

