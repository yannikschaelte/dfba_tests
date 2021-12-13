echo qsub start

# Declare an array of string with cell types
#declare -a opt_method=("TNC")

declare -x dir_results_folder="SLSQP_200"	#in folder: ex6_synthetic/ defined in run_sampling_ex6_on_bonna.py
declare -x nsamples=10000
#declare -a which_samplers=("AM") #,"PT_AM", "PT_AM")
#declare -a nchains_l=("1") #,"3","10")
declare -x which_samplers="AM" #,"PT_AM", "PT_AM")
declare -x nchains_l=1  #,"3","10")


#for cell_type in ${cellNames[@]}; do
#  start_index=$((gene_index*gene_range))
#done


#for which_sampler in ${which_samplers[@]}; do
echo $dir_results_folder
echo $nsamples
echo $which_samplers
echo $nchains_l

sbatch /home/edudkin/dfba/dfba_tests/script_sbatched_sampling_ex6.sh $dir_results_folder $nsamples $which_samplers $nchains_l

#done


echo qsub done

