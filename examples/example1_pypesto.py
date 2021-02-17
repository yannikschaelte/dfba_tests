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
##
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

##
# generate plots of results (in this case using plotlly)
from dfba.plot.plotly import *

import plotly.io as pio

pio.renderers.default = "browser"
# in plotly version 4, default version is different than in 3
pio.templates.default = "plotly_white"
fig = plot_concentrations(concentrations)
fig.show()
# fig = plot_trajectories(trajectories)
# fig.show()
# write results to file
#concentrations.to_csv("concentrations.csv")
#trajectories.to_csv("trajectories.csv")
##
# DATA
# load measurement data (from Eiteman et al., (2008), Fig 2a)
# path_data = '/home/erika/Documents/Projects/DFBA/data_Fig2a_extracted_apprtime.csv'
# data = pd.read_csv(path_data)
# data: first column has to be time!

# simulate data
mu, sigma = 0, 0.02
n_obs = 3
# create dataframe with first column 'time' and subsequent columns observable
index = [0,10,30,50,80,150,200]
noise = np.zeros((n_obs,len(index)))
for i in range(n_obs):
    noise[i,:] = np.random.normal(mu, sigma, len(index))

data_time = concentrations["time"].iloc[index] 
data_biomass = concentrations["Biomass"].iloc[index] + noise[0,:]
data_xylose =  concentrations["Xylose"].iloc[index] + noise[1,:]
data_glucose =  concentrations["Glucose"].iloc[index] + noise[2,:]

data = pd.concat([data_time,data_biomass,data_xylose,data_glucose],axis = 1)

data.to_csv("/home/erika/Documents/Projects/DFBA/results_example1/"
            "simulated_data_sigma_0_01.csv")





##
# import matplotlib.pyplot as plt
# # plot trajectories
# plt.figure()
# plt.plot(concentrations['time'],concentrations['Biomass'], label='Biomass')
# plt.plot(concentrations['time'],concentrations['Glucose'], label='Glucose')
# plt.plot(concentrations['time'],concentrations['Volume'], label='Volume')
# #plt.plot(concentrations['time'],measured, 'kx', label='measured')
# plt.legend()

##
# get t_start, t_end, t_out from measured time in data


def get_t_simu(data):
    measured_time = data.iloc[:,0]
    t_start = measured_time.iloc[0]
    t_end = np.round(measured_time.iloc[-1],decimals=5)
    # minimal time steps
    t_steps = np.zeros(len(measured_time)-1)
    for i_st in range(len(measured_time)-1):
        t_steps[i_st] = measured_time.iloc[i_st+1]-measured_time.iloc[i_st]
    t_out = np.round(min(t_steps),decimals=5)

    return t_start, t_end, t_out


class ObjFunction:
    # define objective function as class: initialization with: model,
    # measured_observables, parameter_names,
    # call the class with the parameters (np.array)
    # return cost of the objective function

    def __init__(self, model, data, par_names, param_scale):
        self.model = model
        # self.measured_observables = measured_observables
        self.data = data
        self.par_names = par_names
        self.param_scale = param_scale

    def __call__(self,parameters):
        # transform log parameters to lin
        if self.param_scale == 'log10':
            self.parameters = 10**parameters
        elif self.param_scale == 'lin':
            self.parameters = parameters

        # numpy parameters to dict
        par_dict = {}
        for i_p in range(len(self.par_names)):
            par_dict[self.par_names[i_p]] = self.parameters[i_p]

        # get t_start, t_end, t_out from measured time
        t_start, t_end, t_out = get_t_simu(self.data)

        self.model.update_parameters(par_dict)
        concentrations, trajectories = self.model.simulate(t_start, t_end,
                                                           t_out)

        # get subset of times, that are present in measurement times.
        indx = []
        for i in range(len(self.data['time'])):
            for j in range(len(concentrations['time'])):
                if np.isclose(concentrations['time'][j],
                              self.data['time'].values[i]):
                    indx.append(j)
        # check if length of found matched indices equals length of measurement
        # time
        if len(indx) != len(self.data['time']):
            ValueError('Some time points are lost. Please check if inferred '
                       't_out is correct, or why the found time subset does '
                       'not match all measurement times.')
        conc_subset = concentrations.iloc[indx]

        obs_names = self.data.columns[1:]  # first column needs to be 'Time'!! obs_names = ['Biomass','Xylose','Glucose']
        # to check!
        # check for exemplary observable (obs_names[0]) if simulation values
        # exist in last time points, if not copy last entries
        if len(conc_subset[obs_names[0]]) < len(self.data):
            # copy last results into not simulated time points
            while len(conc_subset[obs_names[0]]) < len(self.data):
                row_next = conc_subset[-1:].copy()
                row_next['time'] = conc_subset['time'][-1:] + 1
                conc_subset = conc_subset.append(row_next,ignore_index=True)

            # TODO: add nan?



        # Least squares
        cost = 0
        for i_obs in range(len(obs_names)):
            difference = np.asarray(self.data[obs_names[i_obs]]) - np.asarray(
                conc_subset[obs_names[i_obs]])
            cost = cost + np.sum(np.power(difference, 2))
           # print(cost)


        # negloglikelihood


        return cost

##
# initialize class
# param_scale = 'lin'
param_scale = 'log10'
do_optimize = True

par_names = ["K_g","v_gmax","K_z","v_zmax"]
# par_names = ["K_g","K_z"]
obj_function = ObjFunction(dfba_model, data, par_names, param_scale)


# create objective object for pyPESTO
objective2 = pypesto.Objective(fun=obj_function, grad=False, hess=False)

##

where_to_save = "/home/erika/Documents/Projects/DFBA/results_example1/"

# OPTIMIZATION
# define lower and upper bound for parameters optimization
# dim_full = 2
# lb = 0 * np.ones((dim_full, 1)) #log
# ub = 10 * np.ones((dim_full, 1))
par_dict_original = {}
par_dict_original['K_g'] = 0.0027
par_dict_original['v_gmax'] = 10.5
par_dict_original['K_z'] = 0.0165
par_dict_original['v_zmax'] = 6.0


if param_scale == 'lin':
    lb = [0,5]
    ub = [1,13]
    x_sc = ['lin','lin']
    x_g = np.array([[0.05, 8.5]])
elif param_scale == 'log10':
    # lb = [np.log10(1/10*0.5),np.log10(1/10*8.5)]
    # lb = [-3,-3,-3,-3]
    # ub =[3,3,3,3]
    lb=[np.log10(par_dict_original['K_g'])-1,
        np.log10(par_dict_original['v_gmax'])-1,
        np.log10(par_dict_original['K_z'])-1,
        np.log10(par_dict_original['v_zmax'])-1]
    # lb = [np.log10(par_dict_original['K_g']) - 1,
    #       np.log10(par_dict_original['K_z']) - 1]
    ub = [np.log10(par_dict_original['K_g'])+1,
          np.log10(par_dict_original['v_gmax'])+1,
          np.log10(par_dict_original['K_z'])+1,
          np.log10(par_dict_original['v_zmax'])+1]
    # ub = [np.log10(par_dict_original['K_g']) + 1,
    #       np.log10(par_dict_original['K_z']) + 1]
    # K_g = Parameter("K_g", 0.0027)
    # v_gmax = Parameter("v_gmax", 10.5)
    # K_z = Parameter("K_z", 0.0165)
    # v_zmax = Parameter("v_zmax", 6.0)
    # ub = [np.log10(10*0.5),np.log10(10*8.5)]
    x_sc = ['log10','log10','log10','log10']
    # x_g = np.log10(np.array([[0.05, 8.5]]))


problem1 = pypesto.Problem(objective=objective2, lb=lb, ub=ub,
                           copy_objective=False, x_scales=x_sc)#,
                          # x_guesses=x_g)
if do_optimize:
    opt_method = ['L-BFGS-B']#'Pyswarm']#,'TNC']#],'L-BFGS-B','Fides']
    for i_o in range(len(opt_method)):
        # opt_method = 'L-BFGS-B'
        if opt_method[i_o]=='Pyswarm':
            my_optimizer = optimize.PyswarmOptimizer(options={'swarmsize':30, 'minstep': 1e-9,'minfunc' : 1e-9})
            # nstarts = 1
        elif opt_method[i_o]=='Fides':
            my_optimizer = optimize.FidesOptimizer()
        else:
            my_optimizer = optimize.ScipyOptimizer(method=opt_method[i_o],
                                                   options={'eps': 1e-06})  #  'eta':0.5})  # basic optimizer'L-BFGS-B'

        nstarts = 25

        # my_optimizer = pypesto.optimize.DlibOptimizer()

        print('................starting optimization...............')
        # record the history
        history_options = pypesto.HistoryOptions(trace_record=True)
        result = optimize.minimize(problem1, optimizer=my_optimizer, n_starts=nstarts,
                                   history_options=history_options)
        #
        # plot one optimizer history
        # create a reference point from it
        # ref = {'x': result.optimize_result.get_for_key('x')[0],
        #        'fval': result.optimize_result.get_for_key('fval')[0],
        #        'color': [0.2, 0.4, 1., 1.], 'legend': 'optimum'}
        # plt.figure()
        # visualize.optimizer_history(result,
        #                             reference=ref)
        #
        # plt.savefig('/home/erika/Documents/Projects/DFBA/optimizer_history_'+opt_method+'.png')
        #
        # store result
        name_to_save = str(nstarts)+ 'starts_' + opt_method[i_o]
        with open(where_to_save + 'results_' +name_to_save+ '.pickle',
                'wb') as result_file:
                pickle.dump(result.optimize_result, result_file)  # or the full result object
                result_file.close()
        # save optimization result as hdf5 file
        store_filename = where_to_save + 'results_' +name_to_save+ '.hdf5'
        pypesto_result_writer = save_to_hdf5.OptimizationResultHDF5Writer(
            store_filename)
        pypesto_result_writer.write(result)#, overwrite=True)
        pypesto_problem_writer = save_to_hdf5.ProblemHDF5Writer(
            store_filename)
        pypesto_problem_writer.write(problem1)#, overwrite=True)


        # with open('/home/erika/Documents/Projects/DFBA/results_problem.pickle',
        #           'wb') as result_file:
        #     pickle.dump(result.problem,
        #                 result_file)  # or the full result object
        #     result_file.close()
        #
        #
        df = result.optimize_result.as_dataframe(
            ['fval','fval0','n_fval', 'x','x0','n_grad', 'n_hess', 'n_res', 'n_sres', 'time'])
        df['lb'] = str(lb)
        df['ub'] = str(ub)
        #
        df.to_csv(where_to_save + 'df_results_' + name_to_save+ ".csv")
##

##
if param_scale == 'lin':
    save_add = ''
elif param_scale == 'log10':
    save_add = '_log'

## plot waterfall
visualize.waterfall(result, size=(15,6))
plt.savefig(where_to_save + 'waterfall_'+name_to_save+'_'+ save_add +'.png')
# plt.close()
##
# READ optimization results
# result = pickle.load( open("/home/erika/Documents/Projects/DFBA/"
#                            "results_example1/5starts_opt_Kg_Kz/"
#                            "testresults_L-BFGS-B_5starts_opt_kg_kv.pickle", "rb" ) )
# pickle.dump(result, open( where_to_save + 'results_' +name_to_save+ '.pickle', "wb" ) )
##

hdf5_reader = read_from_hdf5.OptimizationResultHDF5Reader(where_to_save
                                                          + 'results_' +name_to_save+ '.hdf5')
result = hdf5_reader.read()
##
load_opt_results=True
name_to_save = str(nstarts)+ 'starts_' + opt_method[0]
store_filename = where_to_save + 'results_' +name_to_save+ '.hdf5'
result = pypesto.store.OptimizationResultHDF5Reader(store_filename)
df = pd.read_csv(where_to_save + 'df_results_' + name_to_save+ ".csv", index_col=0)
##
# simulate trajectories
# x_hat = result.optimize_result.list[0].x
result_nr = 0  # which opt. start should be simulated and plotted?
# x_hat = df['x'][result_nr]
x_hat = result.optimize_result.list[result_nr]['x']
# if load_opt_results:
#     x_hat = x_hat.replace('  ', ',')
#     x_hat = x_hat.replace(' ', ',')
#     x_hat = x_hat.replace(',,',',')
#     x_hat = eval(x_hat)
##

# x_hat = np.array([0.044,8.5])

par_dict = {}
if param_scale == 'log10':
    for i_p in range(len(par_names)):
        par_dict[par_names[i_p]] = 10**x_hat[i_p]
else:
    for i_p in range(len(par_names)):
        par_dict[par_names[i_p]] = x_hat[i_p]



#
dfba_model.update_parameters(par_dict)
# dfba_model.update_parameters(par_dict_original)

t_start, t_end, t_out = get_t_simu(data)
t_out = 0.1
concentrations_best, trajectories_best = dfba_model.simulate(t_start, t_end, t_out)

##
# plot trajectories
observables = data.columns
colors = ['#fff7fb','#ece2f0','#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016c59','#014636']
colors = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f']

plt.figure()
for i_o in range(1,len(observables)):
    plt.plot(concentrations_best['time'],concentrations_best[observables[i_o]],
             label='simulation '+ observables[i_o], color = colors[i_o])
    plt.plot(data['time'],data[observables[i_o]], 'x',color = colors[i_o], label='data' + observables[i_o])

plt.legend()
plt.savefig(where_to_save + 'simu_trajectories_' +
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

ref = {'x': [np.log10(par_dict_original['K_g']),np.log10(par_dict_original['v_gmax']),
             np.log10(par_dict_original['K_z']),np.log10(par_dict_original['v_zmax'])],
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
plt.savefig(where_to_save+ 'parameters_' +
            name_to_save +'_' + save_add +'_x' + str(result_nr)+ '.png')
# visualize.ReferencePoint( x=[par_dict['K_g'],par_dict['v_gmax'],
#                              par_dict['K_z'],par_dict['v_zmax']],
#                           fval=None, color='g', legend=None)
# ax.scatter(x=np.log10(par_dict['K_g']),y=4,marker='o',color='g')
# plt.scatter(x=np.log10(par_dict['v_gmax']),y=3,marker='o',color='k')
# plt.scatter(x=np.log10(par_dict['K_z']),y=2,marker='o',color='k')
# plt.scatter(x=np.log10(par_dict['v_zmax']),y=1,marker='o',color='k')