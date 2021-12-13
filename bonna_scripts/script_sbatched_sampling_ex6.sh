#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5800
#SBATCH --output=%j.slurm.log
#SBATCH --time=47:00:00

echo Begin submission

date
hostname
uname -r

module load GCCcore/7.3.0
source ~/miniconda/etc/profile.d/conda.sh
conda activate dfba-env39
#source ~/env/bin/activate

# start dfba analysis
export PATH=/home/edudkin/dfba:$PATH
export PATH=/home/edudkin/dfba/dfba_tests:$PATH
cd /home/edudkin/dfba/dfba_tests

echo $dir_results_folder
echo $nsamples
echo $which_samplers
echo $nchains_l

srun --ntasks=1 python run_sampling_ex6_on_bonna.py $dir_results_folder $nsamples $which_samplers $nchains_l

echo Done submission




