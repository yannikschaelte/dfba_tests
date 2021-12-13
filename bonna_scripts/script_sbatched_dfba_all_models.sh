#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00

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

allargs=("$@")
echo ${allargs[@]}
echo -----------------------------------

srun --ntasks=1 python run_dfba_bonna_all_models.py ${allargs[@]} 
#$nstarts $opt_method $parallel $model_dir $data_dir $examplename $dir_to $costfunct ${lb[@]} ${ub[@]} $scaling_param_biomass
echo Done submission




