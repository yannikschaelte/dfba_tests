# plot initial x0 values, compared to optimized x-values
# in order to check, if optimization changed initial x0-values
import pandas
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from pypesto.store import (save_to_hdf5, read_from_hdf5)

##
hdf5_reader = \
    read_from_hdf5.OptimizationResultHDF5Reader("/home/erika/Documents/"
                                                "Projects/DFBA/results_example1/"
                                                "results_50starts_L-BFGS-B.hdf5")
result = hdf5_reader.read()

plt.figure()
colors = ['red','darkblue','darkblue','darkblue','grey','lightblue','lightblue'
          ,'lightblue','lightblue']
#Kg, Kz
for i_optrun in range(len(colors)):
    x0 = [result.optimize_result.list[i_optrun]['x0'][0],
          result.optimize_result.list[i_optrun]['x0'][2]]
    x = [result.optimize_result.list[i_optrun]['x'][0],
         result.optimize_result.list[i_optrun]['x'][2]]
    y = range(len(x))



    plt.scatter(x=x0, y=y, marker='o',c='white',edgecolors=colors[i_optrun])
    plt.scatter(x=x, y=y, marker='o',c=colors[i_optrun],edgecolors=colors[i_optrun])

    for i_p in range(len(x)):
        ax = plt.arrow(x0[i_p], y[i_p], x[i_p]-x0[i_p],0,length_includes_head=True,
                      head_width=0.06, head_length=0.08, color='k')

plt.yticks(np.array([0,1]),['Kg','Kz'])

#Kg, Kz
plt.figure()
for i_optrun in range(len(colors)):
    x0 = [result.optimize_result.list[i_optrun]['x0'][1],
          result.optimize_result.list[i_optrun]['x0'][3]]
    x = [result.optimize_result.list[i_optrun]['x'][1],
         result.optimize_result.list[i_optrun]['x'][3]]
    y = range(len(x))



    plt.scatter(x=x0, y=y, marker='o',c='white',edgecolors=colors[i_optrun])
    plt.scatter(x=x, y=y, marker='o',c=colors[i_optrun],edgecolors=colors[i_optrun])

    for i_p in range(len(x)):
        ax = plt.arrow(x0[i_p], y[i_p], x[i_p]-x0[i_p],0,length_includes_head=True,
                      head_width=0.06, head_length=0.08, color='k')

plt.yticks(np.array([0,1]),['vgmax','vzmax'])

