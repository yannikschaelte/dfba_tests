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

grid = False

from os.path import dirname, join, pardir
from cobra.io import read_sbml_model
from dfba import DfbaModel, ExchangeFlux, KineticVariable, Parameter
from pypesto_dfba.optimize_dfba.objective_dfba import (ObjFunction,
                                                       get_obs_names)
from examples.get_dfba_model_ex1_ex6 import get_dfba_model, \
    PicklableDFBAModel, modifun
import matplotlib
if not grid:
    matplotlib.use('TkAgg')
import os
import pypesto
import pypesto.optimize as optimize
import pandas as pd
import pickle
import pypesto.sample as sample
from pypesto.store import (save_to_hdf5, read_from_hdf5)
# from pypesto import Objective, FD
from datetime import datetime
import tempfile
import fides


examplename = "example1_aerobic"
model_dir = "/home/erika/Documents/Projects/DFBA/dynamic-fba/" \
              "sbml-models/iJR904.xml.gz"
data_dir = "/home/erika/Documents/Projects/DFBA/results_example1/" \
             "real_data/data_Fig1.csv"
dir_to = "/home/erika/Documents/Projects/DFBA/results_example1/" \
         "test_scaling_biomass/"
lb = [-3, -1, -4, -1, -3, -3, -3, -3]
ub = [1, 2, 0, 2, 1, 1, 1, 1]
x000 = [[0.18805829, 1.07452258, -4., 0.82564202, -0.65918561,
         -0.72277778, 0.43567432, 0]]   # SLSQP

nstarts = 1
opt_method = 'Fides'
parallel = False
cost_funct='NLLH_normal'
scaling_param_biomass=True

# def run_optimization(model_dir,
#                      examplename,
#                      data_dir,
#                      dir_to,
#                      lb, ub,
#                      nstarts,
#                      opt_method,
#                      parallel,
#                      cost_funct,
#                      x000=None,
#                      scaling_param_biomass=False):

print('optimization method: ' + opt_method)
print("nstarts: " + str(nstarts))
print("Parallel: " + str(parallel))
name_to_save = str(nstarts) + 'starts_' + opt_method + '_' + cost_funct

# save example settings/inputs
pd_info = pd.DataFrame(index= ['dir_to', 'example_name',
                               'model_dir', 'data_dir', 'lb',  'ub',
                               'nstarts', 'opt_method', 'parallel',
                               'cost_function', 'x000_initialGuess'],
                       columns=['content'])
pd_info.loc['dir_to'] = dir_to
pd_info.loc['example_name'] = examplename
pd_info.loc['model_dir'] = model_dir
pd_info.loc['data_dir'] = data_dir
pd_info.loc['lb'] = str(lb)
pd_info.loc['ub'] = str(ub)
pd_info.loc['nstarts'] = nstarts
pd_info.loc['opt_method'] = opt_method
pd_info.loc['parallel'] = parallel
pd_info.loc['cost_function'] = cost_funct
pd_info.loc['x000_initialGuess'] = str(x000)
pd_info.to_csv(os.path.join(dir_to,"pd_info_" + name_to_save + ".csv"))

# dfba_model, params_dict = get_dfba_model(model_dir)
_, params_dict = get_dfba_model(model_dir, examplename)

# params_dict = {"K_g": 0.0027,
#           "v_gmax": 10.5,
#           "K_z": 0.0165,
#           "v_zmax": 6.0}
# get dfba-model in PickableDFBAModel-class
dfba_model = PicklableDFBAModel(model_dir, modifun, examplename)
#
# Simulate model
# t_out = 0.1
# t_end = 25
# print('Start Simulating....')
# dfba_model.solver_data.set_display("none")
# concentrations, trajectories = dfba_model.simulate(
#     0.0, t_end, t_out, ["EX_glc(e)", "EX_xyl_D(e)", "EX_etoh(e)"]
# )
# print(concentrations)

data = pd.read_csv(data_dir, index_col=0)

# initialize Objective Function Class
param_scale = 'log10'
# cost_funct = "LS" |  "NLLH_normal" | "NLLH_laplace"

# observable names
obs_names = get_obs_names(data)
if "NLLH" in cost_funct:
    # define new sigma parameters for each observable
    for i_o in range(len(obs_names)):
        params_dict["sigma_" + obs_names[i_o]] = 1
# add scaling parameter for Biomass (scaling * simulated_biomass)
if scaling_param_biomass:
    params_dict["sc_biomass"] = 1

# parameter names
par_names = list(params_dict.keys())  # ["K_g","v_gmax","K_z","v_zmax"]

obj_function = ObjFunction(dfba_model, data, par_names, param_scale,
                           cost_funct)

# create objective object for pyPESTO
objective_pypesto = pypesto.Objective(fun=obj_function, grad=False, hess=False)
##
if opt_method == 'Fides':
# Fides: objective function, if no hessian_update is provided, this
# function must return a tuple (fval, grad), otherwise this function must
# return a tuple (fval, grad, Hessian)
#     grad2 = FD(Objective(fun=fun, res=res))
    obj_grad = pypesto.FD(objective_pypesto, delta_fun=1e-1, delta_grad=0.1, delta_res=0.1)

# param_opt=[-0.42959639, -0.67591189, -0.32051595,  1.49677907,  0.49781245 , 1.,
#   0.96267552]
param_opt = [0.18999292,  1.07469067, -3.59855968,  0.82761138,
             -0.65918561, -0.65993364, -0.72224595,  0.1376983]

result = obj_grad(param_opt, (0,1), delta_fun=1e-1, delta_grad=0.1, delta_res=0.1)  # (0,1)-> sensi orders: (fun,grad)

##
# test obj-function
param_numpy = [0.0028, 10.5, 0.0165, 6.0]
obj_function_test = ObjFunction(dfba_model, data, par_names, 'lin', "NLLH")
cost = obj_function_test(param_numpy)
print(cost)
# # Test, ob das Pickling, was pyPESTO intern verwendet,
# # um die Objective für die MultiProcessEngine zu kopieren
# import cloudpickle as pickle
# dump = pickle.dumps(objective2)
# pickle.loads(dump)

#
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
    x_sc = ['lin']*len(params_dict)
elif param_scale == 'log10':
    x_sc = ['log10']*len(params_dict)

if parallel:
    # engine = pypesto.engine.MultiThreadEngine(n_threads=2)
    engine = pypesto.engine.MultiProcessEngine() #n_procs=2
else:
    engine = pypesto.engine.SingleCoreEngine()
print(engine)

if len(par_names) != len(lb):
    raise ValueError("Defined lower/upper bound array does not coincide "
                     "with number of parameters! \n Parameters: (" +
                     str(len(par_names)) + ") " + str(par_names)
                     + "\n dimension(lb) = " + str(len(ub)))

# problem1 = pypesto.Problem(objective=objective2, lb=lb, ub=ub,
#                            copy_objective=False, x_scales=x_sc,
#                            x_guesses=x000)
if opt_method == 'Fides':
    problem1 = pypesto.Problem(objective=obj_grad, lb=lb, ub=ub,
                               copy_objective=False, x_scales=x_sc,
                               x_guesses=x000)
else:
    problem1 = pypesto.Problem(objective=objective_pypesto, lb=lb, ub=ub,
                               copy_objective=False, x_scales=x_sc,
                               x_guesses=x000)

# maxls = 40  # maxiter = 10 # maxfun = 10

if do_optimize:
    # Fides: objective function, if no hessian_update is provided, this
    # function must return a tuple (fval, grad), otherwise this function must
    # return a tuple (fval, grad, Hessian)
    #for i_o in range(len(opt_method)):
    if opt_method == 'Pyswarm':
        my_optimizer = optimize.PyswarmOptimizer(options={'swarmsize': 30,
                                                          'minstep': 1e-9,
                                                          'minfunc': 1e-9})
    elif opt_method == 'Fides':
        # trust region optimizer fides
        # my_hess_update = fides.hessian_approximation.HessianApproximation()
        # my_optimizer = optimize.FidesOptimizer(hessian_update=my_hess_update)
        my_optimizer = optimize.FidesOptimizer()
    else:
        my_optimizer = optimize.ScipyOptimizer(method=opt_method,
                                               options={'eps': 1e-09})#,
                                                        #'maxiter': maxiter})
                                                        # 'maxls':maxls,
                                                        # 'maxfun':maxfun}
                                                        # )  #  'eta':0.5})  # basic optimizer'L-BFGS-B'

    # my_optimizer = pypesto.optimize.DlibOptimizer()

    print('----- starting optimization...............')
    # record the history
    now = datetime.now()
    str_now = now.strftime(format="%Y-%m-%d %H:%M:%S").replace(' ', '_')
    file1 = open(os.path.join(dir_to, 'OptStart_' + str_now + '_' +
                              name_to_save + '.txt'), "a")
    file1.close()


    # save optimizer trace (to temporary file fn)
    store_filename = os.path.join(dir_to, 'results_' + name_to_save + '_' +
                                  '.hdf5')
    history_options = pypesto.HistoryOptions(trace_record=True)#,
                                             #storage_file=store_filename) #funktioniert parallel nicht#

    result = optimize.minimize(problem1, optimizer=my_optimizer,
                               n_starts=nstarts,
                               history_options=history_options,
                               engine=engine,
                               options={'allow_failed_starts': False})

    # save optimization result as hdf5 file
    pypesto_result_writer = save_to_hdf5.OptimizationResultHDF5Writer(
        store_filename)
    pypesto_result_writer.write(result, overwrite=True)
    pypesto_problem_writer = save_to_hdf5.ProblemHDF5Writer(
        store_filename)
    pypesto_problem_writer.write(problem1, overwrite=True)

    # store result in pickle
    now = datetime.now()
    str_now = now.strftime(format="%Y-%m-%d %H:%M:%S").replace(' ', '_')
    file2 = open(os.path.join(dir_to, 'OptEnd_' + str_now + '_' +
                              name_to_save + '.txt'), "a")
    file2.close()
    with open(os.path.join(dir_to, 'result_optimize_result_' +
                                   name_to_save + '.pickle'), 'wb') as result_file:
        pickle.dump(result.optimize_result, result_file)  # or the full result object
        result_file.close()

    # save main results in csv
    df = result.optimize_result.as_dataframe(
        ['fval', 'fval0', 'n_fval', 'x', 'x0', 'grad', 'n_grad', 'n_hess',
         'n_res', 'n_sres', 'time'])
    df['lb'] = str(lb)
    df['ub'] = str(ub)

    df.to_csv(os.path.join(dir_to, 'df_results_' + name_to_save+ ".csv"))

    # data.to_csv(os.path.join(dir_to, "simulated_data_sigma_0_01_" +
    #                          name_to_save + ".csv"))
    print("Successfully finished Optimization! ")


if not grid:
    # Example 1 - Synthetic Data:
    # name_ex = "example1"
    # model_direc = "/home/erika/Documents/Projects/DFBA/dynamic-fba/" \
    #             "sbml-models/iJR904.xml.gz"
    # lo_b = [-4, -1, -4, -1, -2,-2,-2]
    # up_b = [-0.5, 2, -0.5, 2,1,1,1]
    # data_direc = "/home/erika/Documents/Projects/DFBA/results_example1/" \
    #              "simulated_data_sigma_0_01_25starts_L-BFGS-B.csv"
    # direc_to = "/home/erika/Documents/Projects/DFBA/results_example1/tests/"#

    # Example 1 - Real data:
    name_ex = "example1_aerobic"
    model_direc = "/home/erika/Documents/Projects/DFBA/dynamic-fba/" \
                  "sbml-models/iJR904.xml.gz"
    data_direc = "/home/erika/Documents/Projects/DFBA/results_example1/" \
                 "real_data/data_Fig1.csv"
    direc_to = "/home/erika/Documents/Projects/DFBA/results_example1/" \
             "test_scaling_biomass/"
    lo_b = [-3, -1, -4, -1, -3, -3, -3, -3]
    up_b = [1, 2, 0, 2, 1, 1, 1, 1]
    # x000 = [[-2.78442879,  1.12691846, -3.20043592,  0.87725271,-2,-2,-2]]
    # x000 = [[-1.95850078, -0.42881366, -1.72887303, -0.58208385],
    #         [-1.95850078, -0.42881366, -1.72887303, -0.58208385]] # ends fast
    # x000 = [[-2.78442879,  1.12691846, -3.20043592,  0.87725271]] # good one
    # x000 = [[-2.78442879 , 1.12691846, -3.20043592,  0.87725271],
    #         [-2.89694802  ,0.56741186 ,-2.13841139 , 1.06010259],
    #         [-3.58966931  ,0.72260446 ,-2.59264231 , 0.96568738],
    #         [-3.9995114 ,1.67331706 ,-3.87885693 , 1.29158789],
    #         [-2.535844    ,1.14214638 ,-2.86931482 , 0.72088211],
    #         [-2.4933603   ,1.19034203 ,-3.56790169 , 1.95862843],
    #         [-3.98467895  ,1.45992462 ,-1.2082032  , 1.5633985 ],
    #         [-1.33089788  ,1.70595286 ,-1.8141665  , 0.74148177]]
    # x000 = [[-2.78442879,  1.12691846, -3.20043592,  0.87725271,
    #          -2,-2,-2,-1]]#good one, with sigma = 0.01
    x000 = [[0.18805829, 1.07452258, -4., 0.82564202, -0.65918561,
             -0.72277778, 0.43567432, 0]]   # SLSQP
#
    # Example 6:
    # name_ex = "example6"
    # model_direc = "/home/erika/Documents/Projects/DFBA/dynamic-fba/" \
    #               "sbml-models/iND750.xml.gz"
    # lo_b = [-3, -1, -3]
    # up_b = [0, 3, 1]
    # data_direc = "/home/erika/Documents/Projects/DFBA/results_example6/" \
    #              "simulated_data/simulated_data_sigma_0.25_ex6_Ausreiser.csv"
    # direc_to = "/home/erika/Documents/Projects/DFBA/results_example6/tests/test/"

    n_starts = 1
    # optimization_method = 'TNC'  # Pyswarm']#'Pyswarm']#,'TNC']#],'L-BFGS-B','SLSQP']
    optimization_method = 'SLSQP'
    # optimization_method = 'Fides'
    run_parallel = False
    cost_funct='NLLH_normal'
    # x000 = [[[0.08348173, 1.23262928, 0.85721095]]] #example 6:

    run_optimization(model_direc, name_ex, data_direc, direc_to, lo_b, up_b,
                     n_starts,
                     optimization_method, run_parallel,
                     cost_funct, x000=x000, scaling_param_biomass=True)#, x000=x000)
