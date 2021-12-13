import sys
from ex6_sampling_bonna import run_sampling
import os

# ###########################################
# Script to run DFBA-SAMPLING on grid
# ###################################################
print("Number of arguments: " + str(len(sys.argv)))
print("The arguments are: " + str(sys.argv))
dir_hdf5_folder_arg = sys.argv[1] #dir_results_folder="SLSQP_200"
n_samples_arg = sys.argv[2]
which_sampler_arg = sys.argv[3]
n_chains_arg = sys.argv[4]


def run_dfba_sampling(dir_hdf5_folder_str, n_samples_str, which_sampler_str,
                      n_chains_str, examplename="example6"):
    # dir_hdf5_folder_str = "TNC_200"
    opt_method = str.split(dir_hdf5_folder_str, '_')[0]
    nstarts = str.split(dir_hdf5_folder_str, '_')[1]
    dir_hdf5_file = os.path.join(".", "results/ex6_synthetic", dir_hdf5_folder_str,
                                 "results_" + str(nstarts) + "starts_" +
                                 opt_method + '_.hdf5')
    model_dir = "../dynamic-fba/sbml-models/iND750.xml.gz"
    data_dir = "./data/simulated_data_sigma_0.25_ex6.csv"

    dir_to = os.path.join(".", "results", dir_hdf5_folder_str)
    # if not os.path.exists(dir_to):
    #     os.makedirs(dir_to)

    # convert input strings to integer
    n_samples = int(n_samples_str)
    n_chains = int(n_chains_str)
    which_sampler = which_sampler_str   # 'AM', 'PT_AM', 'PT_M'

    # RUN OPTIMIZATION FOR DFBA MODEL
    run_sampling(dir_hdf5_file, model_dir, data_dir, dir_to, n_samples,
                 which_sampler, examplename, n_chains)

    return


if __name__ == '__main__':
    # Map command line arguments to function arguments.
    run_dfba_sampling(dir_hdf5_folder_arg, n_samples_arg, which_sampler_arg,
                      n_chains_arg)
