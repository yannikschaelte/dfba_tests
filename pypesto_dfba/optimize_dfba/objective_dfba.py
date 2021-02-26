import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

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

        #original
        # dir_to ='/home/erika/Documents/Projects/DFBA/results_example1/test2/'
        # par_dict_orig =  {"K_g": 0.0027,"v_gmax": 10.5,"K_z": 0.0165, "v_zmax": 6.0 }
        # self.model.update_parameters(par_dict_orig)
        # concentrations_orig, trajectories_orig = self.model.simulate(t_start,
        #                                                    t_end + t_out,
        #                                                    t_out)
        # diff = 1.9443830711793235 - concentrations_orig['Biomass'].iloc[20]
        # if diff == 0:
        #     writet = 'True'
        # else: writet = 'False'
        # Append-adds at last
        # file1 = open(dir_to + "output.txt", "a")  # append mode
        # file1.write(str(diff = 1.9443830711793235 - concentrations_orig['Biomass'].iloc[20]) +
        #             " \n" + writet + "\n" + )
        # file1.close()

        # now = datetime.now()
        # str_now = now.strftime(format="%Y-%m-%d %H:%M:%S").replace(' ', '_')
        # observables = self.data.columns
        # colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854',
        #           '#ffd92f']
        # plt.figure()
        # for i_o in range(1, len(observables)):
        #     plt.plot(concentrations_orig['time'], concentrations_orig[observables[i_o]],
        #              label='simulation ' + observables[i_o], color=colors[i_o])
        #     plt.plot(self.data['time'], self.data[observables[i_o]], 'x',
        #              color=colors[i_o], label='data' + observables[i_o])
        # plt.title('original parameters before \n' + str(par_dict_orig) + '\n updated: '+
        #           str(par_dict) + '\n difference: ' + str(diff))
        # plt.legend()
        # plt.savefig(dir_to + 'simu_trajectories_before_' + str_now +'.png')
        # plt.close()
        # print(diff)

        self.model.update_parameters(par_dict)

        concentrations, trajectories = self.model.simulate(t_start, t_end+t_out,
                                                           t_out)
        print(par_dict)
        print(concentrations)
        # now = datetime.now()
        # str_now = now.strftime(format="%Y-%m-%d %H:%M:%S").replace(' ', '_')
        # observables = self.data.columns
        # colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854',
        #           '#ffd92f']
        # plt.figure()
        # for i_o in range(1, len(observables)):
        #     plt.plot(concentrations['time'], concentrations[observables[i_o]],
        #              label='simulation ' + observables[i_o], color=colors[i_o])
        #     plt.plot(self.data['time'], self.data[observables[i_o]], 'x',
        #              color=colors[i_o], label='data' + observables[i_o])
        # plt.title('parameters update \n' + str(par_dict))
        # plt.legend()
        # plt.savefig(dir_to + 'simu_trajectories_update_' +str_now+ '.png')
        # plt.close()



        # get subset of times, that are present in measurement times.
        indx = []
        for i in range(len(self.data['time'])):
            for j in range(len(concentrations['time'])):
                if np.isclose(concentrations['time'][j],
                              self.data['time'].values[i]):
                    indx.append(j)

        conc_subset = concentrations.iloc[indx]

        # check that first column is time:
        if self.data.columns[0] != 'time':
            raise ValueError('First column of the data needs to be "time"-column!')

        obs_names = self.data.columns[1:]  # first column needs to be 'Time'!! obs_names = ['Biomass','Xylose','Glucose']
        # to check!
        # check for exemplary observable (obs_names[0]) if simulation values
        # exist in last time points, if not copy last entries
        if len(conc_subset) < len(self.data):
            # copy last results into not simulated time points
            while len(conc_subset) < len(self.data):
                row_next = conc_subset[-1:].copy()
                row_next['time'] = self.data['time'].iloc[len(conc_subset[obs_names[0]])]
                    # conc_subset['time'][-1:] + 1
                conc_subset = conc_subset.append(row_next,ignore_index=True)
            # TODO: add nan?

        # check if length of found matched indices equals length of measurement
        # time
        if len(conc_subset) != len(self.data['time']):
            raise ValueError('Some time points are lost. Please check if inferred '
                       't_out is correct, or why the found time subset does '
                       'not match all measurement times.')


        # Least squares
        cost = 0
        for i_obs in range(len(obs_names)):
            difference = np.asarray(self.data[obs_names[i_obs]]) - np.asarray(
                conc_subset[obs_names[i_obs]])
            cost = cost + np.sum(np.power(difference, 2))
           # print(cost)


        # negloglikelihood


        return cost