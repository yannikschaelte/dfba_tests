import sys
from ex1_pypesto_bonna import run_optimization
import os
import numpy as np

# ###########################################
# Script to run DFBA-OPTIMIZATIONS on grid
# function parameters:
# nstarts_str = number of starts for multistart-optimizations  
# opt_method_str = optimization method (TNC; L-BFGS-B,...)
# ###################################################

# number of string arguments (in order to calculate number of
# parameters, to infer size of lower/upper bound array)
num_str_arg = 10    # !!!!! >>> CHANGE HERE IF YOU ADD NEW ARGVS <<< !!!!!!!

nstarts_str = sys.argv[1]
print("Number of arguments: " + str(len(sys.argv)))
print("Number of parameters are: " + str(int((len(sys.argv) - num_str_arg)/2)))
print("The arguments are: " + str(sys.argv))
opt_method_str = sys.argv[2]
parallel_str = sys.argv[3]
model_dir_str = sys.argv[4]
data_dir_str = sys.argv[5]
examplename_str = sys.argv[6]
dir_to_str = sys.argv[7]
cost_func = sys.argv[8]
scaling_param_biomass_str = sys.argv[9]

# get lower bound, upper bound arrays
num_params = int((len(sys.argv) - num_str_arg)/2)
lb_l = []
ub_l = []
for i_b in range(num_params):
    lb_l.append(sys.argv[num_str_arg+i_b])
    ub_l.append(sys.argv[num_str_arg+num_params+i_b])
lb_a = np.asarray(lb_l).astype(float)
ub_a = np.asarray(ub_l).astype(float)
print(lb_a)
print(ub_a)


def run_dfba(nstarts_str, opt_method, parallel_str, model_dir, data_dir,
             examplename, dir_to_folder, cost_funct, lb, ub,
             scaling_param_biomass_str):
    dir_to = os.path.join("./results/", dir_to_folder+ opt_method_str + '_' + cost_funct + '_' +
                          nstarts_str)
    if not os.path.exists(dir_to):
        os.makedirs(dir_to)

    # convert input strings to integer
    nstarts = int(nstarts_str)
    print("lb:")
    print(lb)
    print("ub:")
    print(ub)
    print(type(parallel_str))
    if parallel_str == "True":
        parallel = True
        print('parallel is true, in run_dfba_on_bonna.py')
    else:
        parallel = False
        print("parallel is false, in run_dfba_on_bonna.py")
    parallel = True
    if scaling_param_biomass_str == "True":
        scaling_param_biomass = True
    else:
        scaling_param_biomass = False

    # if opt_method_str=="Fides":
    #     opt_options = {"fatol": 1E-05, "frtol": 1E-05}  # Fides
    # else:
    #     opt_options = None

    # RUN OPTIMIZATION FOR DFBA MODEL
    run_optimization(model_dir, examplename, data_dir, dir_to, lb, ub,
                     nstarts, opt_method, parallel, cost_funct,
                     scaling_param_biomass=scaling_param_biomass,
                     opt_options={"fatol": 1E-05, "frtol": 1E-05})

    return


if __name__ == '__main__':
    # Map command line arguments to function arguments.
    run_dfba(nstarts_str, opt_method_str, parallel_str, model_dir_str,
             data_dir_str, examplename_str, dir_to_str, cost_func,
             lb_a, ub_a, scaling_param_biomass_str)
