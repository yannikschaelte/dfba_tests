from os.path import dirname, join, pardir
from cobra.io import read_sbml_model
from dfba import DfbaModel, ExchangeFlux, KineticVariable, Parameter
from pypesto_dfba.optimize_dfba.objective_dfba import (ObjFunction,
                                                       get_obs_names)
from examples.get_dfba_model_ex1_ex6 import get_dfba_model, \
    PicklableDFBAModel, modifun
grid = False
import matplotlib
if not grid:
    matplotlib.use('TkAgg')
import os
import pypesto
import pypesto.optimize as optimize
import pandas as pd
import pickle
import numpy as np
import pypesto.sample as sample
from pypesto.store import (save_to_hdf5, read_from_hdf5)
# from pypesto import Objective, FD
from datetime import datetime
import tempfile
import fides
import csv

history_path = "/home/erika/Documents/Projects/DFBA/results_example1/" \
               "real_data/Fides2/Fides_NLLH_normal_32/" \
               "history_32starts_Fides_NLLH_normal/"

filename = "history_32starts_Fides_NLLH_normal_"

# write hdf5-result to pandasDataframe
pd_hist = pd.DataFrame(columns=['id','x','fval','n_fval','fval0','x0','grad',
                                    'hess','res','sres','n_grad','n_hess','n_res',
                                    'n_sres','exitflag','time','message'], index=[0])
# Read result
for ii in range(32):
    hist = pypesto.store.optimization_result_from_history(history_path + filename +
                                                          str(ii) +".hdf5")

    pd_hist.loc[ii,'id'] = hist.optimize_result.list[0]['id']
    pd_hist.loc[ii,'x'] = hist.optimize_result.list[0]['x']
    pd_hist.loc[ii,'fval'] = hist.optimize_result.list[0]['fval']
    pd_hist.loc[ii,'n_fval'] = hist.optimize_result.list[0]['n_fval']
    pd_hist.loc[ii,'fval0'] = hist.optimize_result.list[0]['fval0']
    pd_hist.loc[ii,'x0'] = hist.optimize_result.list[0]['x0']
    pd_hist.loc[ii,'grad'] = hist.optimize_result.list[0]['grad']
    pd_hist.loc[ii,'hess'] = hist.optimize_result.list[0]['hess']
    pd_hist.loc[ii,'res'] = hist.optimize_result.list[0]['res']
    pd_hist.loc[ii,'sres'] = hist.optimize_result.list[0]['sres']
    pd_hist.loc[ii,'n_grad'] = hist.optimize_result.list[0]['n_grad']
    pd_hist.loc[ii,'n_hess'] = hist.optimize_result.list[0]['n_hess']
    pd_hist.loc[ii,'n_res'] = hist.optimize_result.list[0]['n_res']
    pd_hist.loc[ii,'n_sres'] = hist.optimize_result.list[0]['n_sres']
    pd_hist.loc[ii,'exitflag'] = hist.optimize_result.list[0]['exitflag']
    pd_hist.loc[ii,'time'] = hist.optimize_result.list[0]['time']
    pd_hist.loc[ii,'message'] = hist.optimize_result.list[0]['message']

pd_hist.to_csv(history_path + filename + '_SUMMARY.csv')