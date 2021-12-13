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

"""Aerobic growth of *S. cerevisiae* on glucose.

Organism -> Saccharomyces cerevisiae S288C
Model stored in http://bigg.ucsd.edu/models/iND750
"""

from os.path import dirname, join, pardir

import numpy as np
from cobra.io import read_sbml_model

from dfba import DfbaModel, ExchangeFlux, KineticVariable, Parameter

import matplotlib
matplotlib.use('TkAgg')
#
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
##

do_optimize = True # run optimization?

# param_scale = 'lin'
param_scale = 'log10'



# define fixed parameters and instances of Parameters (to be used for estimation)
# the latter are added to model below (see line 69)
Vgmax = Parameter("Vgmax",8.5)      # maximum uptake rate of glucose
# Vgmax = Parameter("Vgmax",7.3)#changed
# Kg = Parameter("Kg",1.026)                 # saturation constant of glucose
Kg = Parameter("Kg",0.5)
# Kg = 1.026#changed
# D = Parameter("D",0.044)
D = 0.044
Gin = 100.0

# create DfbaModel instance initialized with cobra model
# this model is related to example 3 without control parameter
# fba_model = read_sbml_model(join(dirname(__file__), pardir, "sbml-models", "iND750.xml.gz"))
fba_model = read_sbml_model("/home/erika/Documents/Projects/DFBA/dynamic-fba/sbml-models/iND750.xml.gz")
##
fba_model.solver = "glpk"
dfba_model = DfbaModel(fba_model)
V = KineticVariable("Volume")
X = KineticVariable("Biomass")
Gluc = KineticVariable("Glucose")
Eth = KineticVariable("Ethanol")
Glyc = KineticVariable("Glycerol")
dfba_model.add_kinetic_variables([V, X, Gluc, Eth, Glyc])
mu = ExchangeFlux("BIOMASS_SC4_bal")
v_G = ExchangeFlux("EX_glc__D_e")
v_E = ExchangeFlux("EX_etoh_e")
v_H = ExchangeFlux("EX_glyc_e")
dfba_model.add_exchange_fluxes([mu, v_G, v_E, v_H])

dfba_model.add_parameters([Kg,Vgmax]) # Here add parameters to dfba_model
dfba_model.add_rhs_expression("Volume", D)
dfba_model.add_rhs_expression("Biomass", mu * X - D * X / V)
dfba_model.add_rhs_expression("Glucose", v_G * X + D * (Gin - Gluc) / V)
dfba_model.add_rhs_expression("Ethanol", v_E * X - D * Eth / V)
dfba_model.add_rhs_expression("Glycerol", v_H * X - D * Glyc / V)
dfba_model.add_exchange_flux_lb("EX_glc__D_e", Vgmax * (Gluc / (Kg + Gluc)), Gluc)  #lower bound
dfba_model.add_initial_conditions(
    {
        "Volume": 0.5,
        "Biomass": 0.05,
        "Glucose": 16.0,
        "Ethanol": 0.0,
        "Glycerol": 0.0,
    }
)   #"Glucose": 10.0,
##
# perform simulation with initial parameter values
# add noise to Biomass to simulate measurement
t_end = 16.0
dfba_model.solver_data.set_display("none")
concentrations, trajectories = dfba_model.simulate(0.0, t_end, 1.0,)

mu, sigma = 0, 0.01
noise = np.random.normal(mu, sigma, [17])
# measured = (concentrations["Biomass"] + noise).tolist()
measured = (concentrations["Biomass"] + noise)
# create dataframe with first column 'time' and subsequent columns observable
data = pd.DataFrame(np.arange(0,t_end),columns=['time'])
data['Biomass'] = measured

#
# create parameters dict with original parameter values
# use to calculate objective function value in grid scan
# parameters = {"D": 0.044, "Vgmax": 8.5} # -> write in numpy array
parameters = {"Kg": 0.5, "Vgmax": 8.5} # -> write in numpy array

## plot trajectories
from dfba.plot.plotly import *
import plotly.io as pio

pio.renderers.default = "browser"
# in plotly version 4, default version is different than in 3
pio.templates.default = "plotly_white"
fig = plot_concentrations(concentrations)
fig.show()

# fig = plot_trajectories(trajectories)
# fig.show()


##
# plot trajectories
plt.figure()
plt.plot(concentrations['time'],concentrations['Biomass'], label='Biomass')
plt.plot(concentrations['time'],concentrations['Glucose'], label='Glucose')
plt.plot(concentrations['time'],concentrations['Volume'], label='Volume')
plt.plot(concentrations['time'],measured, 'kx', label='measured')
plt.legend()
# plt.savefig('/home/erika/Documents/Projects/DFBA/simu_trajectories_basis_.png')

##
# # define objective function as class: initialization with: model,
# # measured_observables, parameter_names,
# # call the class with the parameters (np.array)
# # return cost of the objective function
#
#
# class ObjFunction:
#
#     def __init__(self, model, measured_observables, par_names, param_scale, t_end):
#         self.model = model
#         self.measured_observables = measured_observables
#         self.par_names = par_names
#         self.param_scale = param_scale
#
#     def __call__(self,parameters):
#         # transform log parameters to lin
#         if self.param_scale == 'log10':
#             self.parameters = 10**parameters
#         elif self.param_scale == 'lin':
#             self.parameters = parameters
#
#         # numpy to dict
#         par_dict = {self.par_names[0]: self.parameters[0],
#                     self.par_names[1]: self.parameters[1]}
#
#         self.model.update_parameters(par_dict)
#         concentrations, trajectories = self.model.simulate(0.0, t_end, 1.0)
#
#         if len(concentrations["Biomass"]) < len(measured):
#             # copy last results into not simulated time points
#             while len(concentrations["Biomass"]) < len(measured):
#                 row_next = concentrations[-1:].copy()
#                 row_next['time'] = concentrations['time'][-1:] + 1
#                 concentrations = concentrations.append(row_next,ignore_index=True)
#
#         # Least squares
#         difference = np.asarray(measured) - np.asarray(
#             concentrations["Biomass"])
#         cost = np.sum(np.power(difference, 2))
#         # negloglikelihood
#
#
#         return cost

##


# get t_start, t_end, t_out from measured time in data
def get_t_simu(data):
    measured_time = data.iloc[:,0]
    t_start = measured_time.iloc[0]
    t_end = measured_time.iloc[-1]
    # minimal time steps
    t_steps = np.zeros(len(measured_time)-1)
    for i_st in range(len(measured_time)-1):
        t_steps[i_st] = measured_time.iloc[i_st+1]-measured_time.iloc[i_st]
    t_out = min(t_steps)

    return t_start, t_end, t_out

# define objective function as class: initialization with: model,
# measured_observables, parameter_names,
# call the class with the parameters (np.array)
# return cost of the objective function


class ObjFunction:

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
        concentrations, trajectories = self.model.simulate(t_start, t_end, t_out)

        conc_subset = concentrations[concentrations['time'].isin(self.data['time'])]

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

        # Least squares
        cost = 0
        for i_obs in range(len(obs_names)):
            difference = np.asarray(self.data[obs_names[i_obs]]) - np.asarray(
                conc_subset[obs_names[i_obs]])
            cost = cost + np.sum(np.power(difference, 2))

        # negloglikelihood

        return cost
##
# initialize class
# par_names = ['D','Vgmax']
par_names = ['Kg','Vgmax']
obj_function = ObjFunction(dfba_model, data, par_names, param_scale)
# obj_function = ObjFunction(dfba_model, par_names, param_scale)
#test
# obj_function22 = ObjFunction(dfba_model, measured, par_names, 'lin', t_end)
# obj_function22([0.044,8.5])
# obj_function22([1.1,10.5])    #-> leads to inf

# create objective object for pyPESTO
objective2 = pypesto.Objective(fun=obj_function, grad=False, hess=False)

# call class
# params_mat = 5 * np.random.random_sample((20, 2))
# params = np.array([0.044,8.5])
# for i in range(20):
#     print('-----------------------' + str(i) + '---------------')
#     params = params_mat[i]
#     obj_function(params)
## TEST ON PARAMTER GRID
#test on grid how often fails
# grid_x = np.arange(0.45, 0.55, 0.01)
# grid_y = np.arange(7.5, 9.6, 0.2)
# #
# df_cost=pd.DataFrame(np.zeros([len(grid_y),len(grid_x)]), index=grid_y,
#                      columns=grid_x)
# for ix in range(0,len(grid_x)):
#     for iy in range(0,len(grid_y)):
#         obj_function22 = ObjFunction(dfba_model, measured, par_names, 'lin', t_end)
#         df_cost.iloc[iy,ix] = obj_function22([grid_x[ix],grid_y[iy]])
#         print(df_cost)

## #
# df_cost.to_csv("/home/erika/Documents/Projects/DFBA/df_cost_grid_Kg_Vgmax_05_85_10percent.csv")
##
# OPTIMIZATION
# define lower and upper bound for parameters optimization
# dim_full = 2
# lb = 0 * np.ones((dim_full, 1)) #log
# ub = 10 * np.ones((dim_full, 1))
if param_scale == 'lin':
    lb = [0,5]
    ub = [1,13]
    x_sc = ['lin','lin']
    x_g = np.array([[0.05, 8.5]])
elif param_scale == 'log10':
    lb = [np.log10(1/10*0.5),np.log10(1/10*8.5)]
    lb = [-3,-1]
    ub =[1,3]
    # ub = [np.log10(10*0.5),np.log10(10*8.5)]
    x_sc = ['log10','log10']
    # x_g = np.log10(np.array([[0.05, 8.5]]))


problem1 = pypesto.Problem(objective=objective2, lb=lb, ub=ub,
                           copy_objective=False, x_scales=x_sc)#,
                          # x_guesses=x_g)

# # def standard_sampling():
# #     objective = pypesto.Objective(fun=nllh)
# #     problem = pypesto.Problem(objective=objective,
# #                               lb=[-5, -5, -np.inf, -np.inf], ub=[5, 5, np.inf, np.inf],
# #                               x_names=['k1', 'k2', 'scaling', 'sigma'],
# #                               x_scales=['log10', 'log10', 'log10', 'log10'])
# #     return problem
# #
#
if do_optimize:
    opt_method = ['TNC']#'Pyswarm']#,'TNC']#],'L-BFGS-B']
    for i_o in range(len(opt_method)):
        # opt_method = 'L-BFGS-B'

        if opt_method[i_o]=='Pyswarm':
            my_optimizer = optimize.PyswarmOptimizer(options={'swarmsize':30, 'minstep': 1e-9,'minfunc' : 1e-9})
            nstarts = 1
        else:
            my_optimizer = optimize.ScipyOptimizer(method=opt_method[i_o],
                                                   options={'eps': 1e-06})  #  'eta':0.5})  # basic optimizer'L-BFGS-B'
            nstarts = 3

        # my_optimizer = pypesto.optimize.DlibOptimizer()

        print('................starting optimization...............')
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
        with open('/home/erika/Documents/Projects/DFBA/testresults_'+opt_method[i_o]+'.pickle',
                'wb') as result_file:
                pickle.dump(result.optimize_result, result_file)  # or the full result object
                result_file.close()

        # with open('/home/erika/Documents/Projects/DFBA/results_problem.pickle',
        #           'wb') as result_file:
        #     pickle.dump(result.problem,
        #                 result_file)  # or the full result object
        #     result_file.close()
        #
        #
        df = result.optimize_result.as_dataframe(
            ['fval','n_fval', 'x','x0','n_grad', 'n_hess', 'n_res', 'n_sres', 'time'])
        df['lb'] = str(lb)
        df['ub'] = str(ub)
        #
        df.to_csv("/home/erika/Documents/Projects/DFBA/df_results_"+opt_method[i_o]+"_"+str(nstarts) + "runs.csv")
##
##
# if param_scale == 'lin':
#     save_add = ''
# elif param_scale == 'log10':
#     save_add = '_log'
#
# ## plot waterfall
# visualize.waterfall(result, size=(15,6))
# plt.savefig('/home/erika/Documents/Projects/DFBA/waterfall_'+str(nstarts) +'runs_'+ save_add +'.png')
# plt.close()
# ##
# # simulate trajectories
# # x_hat = result.optimize_result.list[0].x
# x_hat = df['x'][0]
# # x_hat = np.array([0.044,8.5])
# # par_dict = {"D": 0.044, "Vgmax": 8.5}
#
# if param_scale == 'log10':
#     par_dict = {"Kg": 10**x_hat[0],
#                 "Vgmax": 10**x_hat[1]}
# else:
#     par_dict = {"Kg": x_hat[0],
#                 "Vgmax": x_hat[1]}
# # par_dict = {"D": 0.044, "Vgmax": 8.5}
# dfba_model.update_parameters(par_dict)
# concentrations_best, trajectories_best = dfba_model.simulate(0.0, 16.0, 1.0)

##
# obj_function = ObjFunction(dfba_model, measured, par_names, 'lin')
# obj_function([0.044,8.5])
##
# plot trajectories
# plt.figure()
# plt.plot(concentrations_best['time'],concentrations_best['Biomass'], label='simulation Biomass')
# plt.plot(concentrations_best['time'],measured, 'kx', label='measured')
# plt.plot(concentrations_best['time'],concentrations_best['Glucose'], label='Glucose')
# plt.plot(concentrations_best['time'],concentrations_best['Glycerol'], label='Glycerol')
# plt.plot(concentrations_best['time'],concentrations_best['Volume'], label='Volume')
# plt.legend()
# plt.savefig('/home/erika/Documents/Projects/DFBA/simu_trajectories_' +str(nstarts) +'_' + save_add +'.png')
# plt.close()
# # ##
# # #
# ## #
# open result
# folder = 'results_optimization_TNC_10starts_05_2_lbub'
# opt_method = ['TNC']
# file = open('/home/erika/Documents/Projects/DFBA/'+folder+'/testresults_'+opt_method[i_o]+'.pickle', "rb")
# result.optimize_result = pickle.load(file)



# plt.figure()
# visualize.parameters(result,
#                      balance_alpha=False, size=[10,6])
# plt.scatter(x=np.log10(parameters['Kg']),y=2,marker='o',color='k')
# plt.scatter(x=np.log10(parameters['Vgmax']),y=1,marker='o',color='k')
#
# plt.savefig('/home/erika/Documents/Projects/DFBA/parameters' + save_add +'.svg',bbox_inches='tight')
# plt.savefig('/home/erika/Documents/Projects/DFBA/parameters' + save_add +'.png',bbox_inches='tight')
# plt.close()
# ##
# # plot x0 and x
# i_run = 0
# # param_scale='log10'
# colors = ['#fff7fb','#ece2f0','#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016c59','#014636']
# colors.reverse()
#
# plt.figure()
#
# for i_run in range(nstarts):
#
#     if param_scale=='log10':
#         lb_list = result.problem.lb.astype(float)
#         ub_list = result.problem.ub.astype(float)
#         x_est = result.optimize_result.list[i_run].x
#         x_0 = result.optimize_result.list[i_run].x0
#         param1 = np.log10(parameters['Kg'])
#         param2 = np.log10(parameters['Vgmax'])
#         xlab = 'log10(Kg)'
#         ylab='log10(Vgmax)'
#     else:
#         lb_list = result.problem.lb.astype(float)
#         ub_list = result.problem.ub.astype(float)
#         x_est = result.optimize_result.list[i_run].x
#         x_0 = result.optimize_result.list[i_run].x0
#         param1 = parameters['Kg']
#         param2 = parameters['Vgmax']
#         xlab = 'Kg'
#         ylab = 'Vgmax'
#
#     plt.scatter(x=parameters['Kg'],y=parameters['Vgmax'],marker='x',s=20,c='r',label='best')
#     plt.scatter(x=x_est[0],y=x_est[1],label='x_hat',c=colors[i_run])
#     plt.scatter(x=x_0[0],y=x_0[1],label='x_0',c='w',edgecolors=colors[i_run])
#     plt.plot([x_est[0], x_0[0]], [x_est[1],x_0[1]], color=colors[i_run], linestyle='-', linewidth=2)
# plt.scatter(x=param1,y=param2,marker='x',s=40,c='r')
#
# plt.xlim([lb_list[0],ub_list[0]])
# plt.ylim([lb_list[1],ub_list[1]])
# plt.xlabel(xlab)
# plt.ylabel(ylab)
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.savefig('/home/erika/Documents/Projects/DFBA/x0_x' +str(nstarts) +'_' + save_add +'.svg', bbox_inches = 'tight')
# plt.savefig('/home/erika/Documents/Projects/DFBA/x0_x' +str(nstarts) +'_' + save_add +'.png', bbox_inches = 'tight')
# plt.close()














## ##
# # ##
# sampler = sample.AdaptiveParallelTemperingSampler(
#     internal_sampler=sample.AdaptiveMetropolisSampler(),
#     n_chains=3)   #-> deepcopy problem
sampler = sample.AdaptiveMetropolisSampler()
# sampler = sample.ParallelTemperingSampler(internal_sampler=sample.MetropolisSampler(),
#                                           betas=[1, 1e-1, 1e-2])  # -> deepcopy problem
#
# result = sample.sample(problem1, n_samples=100, sampler=sampler, result=result)
n_samples=1000
result = sample.sample(problem1, n_samples=n_samples, sampler=sampler,
                       x0 = np.array([np.log10(0.5),  np.log10(8.5)]))
##
with open('/home/erika/Documents/Projects/DFBA/results_sampling_AM_' +str(n_samples) + '.pickle',
          'wb') as result_file:
    pickle.dump(result.sample_result,
                result_file)  # or the full result object
    result_file.close()
##
# plt.figure()
visualize.sampling_parameters_trace(result, use_problem_bounds=True, size=(12,5))
# plt.savefig('/home/erika/Documents/Projects/DFBA/sampling_AM_' +str(n_samples) +'.png', bbox_inches = 'tight')
plt.show()
# plt.close()
##
plt.figure()
visualize.sampling_fval_trace(result, size=(12,5))
# plt.savefig('/home/erika/Documents/Projects/DFBA/sampling_fval_AM_' +str(n_samples) +'.png', bbox_inches = 'tight')
# plt.close()
plt.figure()
visualize.sampling_scatter(result, size=[13,6])
# plt.savefig('/home/erika/Documents/Projects/DFBA/sampling_scatter_AM_' +str(n_samples) +'.png', bbox_inches = 'tight')
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
