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

import pypesto
import pypesto.visualize as visualize
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# define example objective function.
# for now takes model instance, list of measured observables, and dict of parameter values
# calculates simple point-wise least squares cost function for Biomass
def objective_function(model, measured_observables, parameters):
    dfba_model.update_parameters(parameters)
    concentrations, trajectories = dfba_model.simulate(0.0, 16.0, 1.0)

    cost = 0.0
    difference = measured - concentrations["Biomass"]
    for diff in difference:
        cost += diff**2
    return cost



# define fixed parameters and instances of Parameters (to be used for estimation)
# the latter are added to model below (see line 69)
Vgmax = Parameter("Vgmax",8.5)
Kg = 0.5
D = Parameter("D",0.044)
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
dfba_model.add_parameters([D,Vgmax]) # Here add parameters to dfba_model
dfba_model.add_rhs_expression("Volume", D)
dfba_model.add_rhs_expression("Biomass", mu * X - D * X / V)
dfba_model.add_rhs_expression("Glucose", v_G * X + D * (Gin - Gluc) / V)
dfba_model.add_rhs_expression("Ethanol", v_E * X - D * Eth / V)
dfba_model.add_rhs_expression("Glycerol", v_H * X - D * Glyc / V)
dfba_model.add_exchange_flux_lb("EX_glc__D_e", Vgmax * (Gluc / (Kg + Gluc)), Gluc)
dfba_model.add_initial_conditions(
    {
        "Volume": 0.5,
        "Biomass": 0.05,
        "Glucose": 10.0,
        "Ethanol": 0.0,
        "Glycerol": 0.0,
    }
)
##
# perform simulation with initial parameter values
# add noise to Biomass to simulate measurement
concentrations, trajectories = dfba_model.simulate(0.0, 16.0, 1.0)

mu, sigma = 0, 0.01
noise = np.random.normal(mu, sigma, [17])
measured = (concentrations["Biomass"] + noise).tolist()
##
# create parameters dict with original parameter values
# use to calculate objective function value in grid scan
parameters = {"D": 0.044, "Vgmax": 8.5} # -> write in numpy array

##
for i in range(1,5):
    D_param = 0.04 + i*0.002
    parameters["D"] = D_param
    for j in range(1,5):
        Vgmax_param = 8 + j*0.25
        parameters["Vgmax"] = Vgmax_param
        obj_val = objective_function(dfba_model,measured,parameters)
        print("Objective val with D = " + str(D_param) + ", Vgmax = " + str(Vgmax_param) + " is " + str(obj_val))


##
# obj function ggf. als klasse mit __init__(model, observables) und __call__(parameters)
class test():

    def __init__(self, model, measured_observables):
        self.model = model
        self.measured_observables = measured_observables

    def __call__(self,parameters):
        self.parameters = parameters

        # numpy to dict
        print(self.parameters[0])
        print(self.parameters[1])
        par_dict = {"D": self.parameters[0], "Vgmax": self.parameters[1]}

        self.model.update_parameters(par_dict)
        concentrations, trajectories = self.model.simulate(0.0, 16.0, 1.0)

        cost = 0.0
        difference = self.measured_observables - concentrations["Biomass"]
        for diff in difference:
            cost += diff ** 2

        print(cost)
        return cost




##


o2 = test(dfba_model,measured)  #initialize class


params = np.array([0.044,8.5])
# o2.objective_function(params)
o2(params)


##
# transfor parameters-dict into numpy-array
# for key, value in zip(ordered_parameter_keys, parameter_array):
#     parameter[key] = value

# objective2 = pypesto.Objective(fun=o2.objective_function, grad=False, hess=False)

objective2 = pypesto.Objective(fun=o2, grad=False, hess=False)
##
dim_full = 2
lb = 1 * np.ones((dim_full, 1))
ub = 5 * np.ones((dim_full, 1))
##
problem1 = pypesto.Problem(objective=objective2, lb=lb, ub=ub)
##
# import pypesto.optimize as optimize
# result = optimize.minimize(problem1, n_starts=10)
# ##
# import pypesto.sample as sample
# sampler = sample.AdaptiveParallelTemperingSampler(
#     internal_sampler=sample.AdaptiveMetropolisSampler(),
#     n_chains=3)
#
# result = sample.sample(problem1, n_samples=10000, sampler=sampler, result=result)
