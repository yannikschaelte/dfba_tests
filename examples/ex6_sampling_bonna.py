# STAND: 06.Juni.2021 (not tested on newer functions of ex1_pypesto_bonna.py)
#
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

grid = True

from os.path import dirname, join, pardir

from cobra.io import read_sbml_model

from dfba import DfbaModel, ExchangeFlux, KineticVariable, Parameter
from pypesto_dfba.optimize_dfba.objective_dfba import (ObjFunction, get_t_simu)
from examples.get_dfba_model_ex1_ex6 import get_dfba_model, PicklableDFBAModel, modifun

import matplotlib
if not grid:
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
import tempfile

# activate debugging
# import logging
# import sys
# root = logging.getLogger()
# root.setLevel(logging.DEBUG)
# handler = logging.StreamHandler(sys.stdout)
# handler.setLevel(logging.DEBUG)
# # create formatter
# formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# handler.setFormatter(formatter)
# root.addHandler(handler)


def init_pypesto_problem(model_dir: str, data_dir: str,
                         lb: np.array, ub: np.array, parallel: bool,
                         examplename:str):
    # dfba_model, params_dict = get_dfba_model(model_dir)
    _, params_dict = get_dfba_model(model_dir, examplename)

    # params_dict = {"K_g": 0.0027,
    #           "v_gmax": 10.5,
    #           "K_z": 0.0165,
    #           "v_zmax": 6.0}
    # get dfba-model in PickableDFBAModel-class
    dfba_model = PicklableDFBAModel(model_dir, modifun, examplename)

    data = pd.read_csv(data_dir, index_col=0)

    #
    # initialize Objective Function Class
    param_scale = 'log10'
    par_names = list(params_dict.keys())  # ["K_g","v_gmax","K_z","v_zmax"]

    cost_funct = "LS"   # TODO: initialize funtion with cost-funct
    obj_function = ObjFunction(dfba_model, data, par_names, param_scale,
                               cost_funct)

    # create objective object for pyPESTO
    objective2 = pypesto.Objective(fun=obj_function, grad=False, hess=False)

    # pyPESTO-PROBLEM
    # define lower and upper bound for parameters optimization
    if param_scale == 'lin':
        x_sc = ['lin']*len(params_dict)
    elif param_scale == 'log10':
        x_sc = ['log10']*len(params_dict)

    if parallel:
        # engine = pypesto.engine.MultiThreadEngine(n_threads=2)
        engine = pypesto.engine.MultiProcessEngine()  # n_procs=2
    else:
        engine = pypesto.engine.SingleCoreEngine()

    # x000 = [[-1.95850078, -0.42881366, -1.72887303, -0.58208385],
    #         [-1.95850078, -0.42881366, -1.72887303, -0.58208385]] # ends fast
    # x000 = [[-2.78442879,  1.12691846, -3.20043592,  0.87725271]] # good one

    problem1 = pypesto.Problem(objective=objective2, lb=lb, ub=ub,
                               copy_objective=False, x_scales=x_sc)  # ,
    # x_guesses=x000)
    return problem1, engine, objective2


# ###############################################################
# -------------------SAMPLING ----------------------------------#


def run_sampling(dir_hdf5_file, model_dir, data_dir, dir_to, n_samples,
                 which_sampler, examplename,
                 n_chains=None
                 ):

    # read optimization result
    hdf5_reader = read_from_hdf5.OptimizationResultHDF5Reader(dir_hdf5_file)
    result = hdf5_reader.read()
    problem_hdf5 = result.problem

    problem1, _, objective2 = init_pypesto_problem(model_dir, data_dir,
                                                   problem_hdf5.lb,
                                                   problem_hdf5.ub,
                                                   parallel=False,
                                                   examplename=examplename)
    x_hat = result.optimize_result.list[0]['x']
    name_to_save = os.path.split(dir_hdf5_file)[1][:-5]

    # ---------------------------------------------------------------------------
    if which_sampler == "AM":
        sampler = sample.AdaptiveMetropolisSampler()
    elif which_sampler == "PT_AM":
        sampler = sample.AdaptiveParallelTemperingSampler(
            internal_sampler=sample.AdaptiveMetropolisSampler(),
            n_chains=n_chains)
    elif which_sampler == "PT_M":
        sampler = sample.ParallelTemperingSampler(internal_sampler=sample.MetropolisSampler(),
                                                  n_chains=n_chains)#betas=[1, 1e-1, 1e-2],
    else:
        raise ValueError("Choose one of the following samplers: which_sampler="
                         "'AM' (Adaptive Metropolis),"
                         " 'PT_AM' (Parallel Tempering w internal: Adaptive "
                         "Metropolis), "
                         "'PT_M' (Parallel Tempering w internal: Metropolis).")

    now = datetime.now()
    str_now = now.strftime(format="%Y-%m-%d %H:%M:%S").replace(' ', '_')
    file1 = open(os.path.join(dir_to, 'SamplStart_' + str_now + '_' +
                              name_to_save + '.txt'), "a")
    file1.close()

    # start SAMPLING
    result = sample.sample(problem1, n_samples=n_samples, sampler=sampler,
                           x0=x_hat)
    # end SAMPLING

    now = datetime.now()
    str_now = now.strftime(format="%Y-%m-%d %H:%M:%S").replace(' ', '_')
    file2 = open(os.path.join(dir_to, 'SamplEnd_' + str_now + '_' +
                              name_to_save + '.txt'), "a")
    file2.close()

    # save results
    with open(os.path.join(dir_to, 'sampling_' + str(n_samples) + '_' +
                                   which_sampler + '_' + str(n_chains) +
                                   'ch_' + name_to_save + '.pickle'),
              'wb') as result_file:
        pickle.dump(result.sample_result, result_file)  # or the full result object
        result_file.close()

    store_filename2 = os.path.join(dir_to, 'results_after_sampling_' +
                                   str(n_samples) + '_' + which_sampler +
                                   '_' + str(n_chains) + 'ch_' + name_to_save +
                                   '.hdf5')
    pypesto_result_writer = save_to_hdf5.OptimizationResultHDF5Writer(
        store_filename2)
    pypesto_result_writer.write(result, overwrite=True)
    pypesto_problem_writer = save_to_hdf5.ProblemHDF5Writer(
        store_filename2)
    pypesto_problem_writer.write(problem1, overwrite=True)


if not grid:
    model_dir = "/home/erika/Documents/Projects/DFBA/dynamic-fba/" \
                    "sbml-models/iJR904.xml.gz"
    data_dir = "/home/erika/Documents/Projects/DFBA/results_example1/" \
                   "simulated_data_sigma_0_01_25starts_L-BFGS-B.csv"
    dir_hdf5_file = "/home/erika/Documents/Projects/DFBA/results_example1/tests/" \
                    "results_2starts_SLSQP_.hdf5"
    n_samples = 5
    dir_to = "/home/erika/Documents/Projects/DFBA/results_example1/" \
                 "tests/"
    which_sampler = "PT_AM"    #"PT_AM", PT_M
    n_chains = 3
    run_sampling(dir_hdf5_file, model_dir, data_dir, dir_to, n_samples,
                 which_sampler, n_chains)
