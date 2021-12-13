echo qsub start

# Declare an array of string with cell types
#declare -a opt_method=("TNC")

declare -x opt_method="SLSQP"
declare -x nstarts=1
declare -x parallel="True"
declare -x model_dir="../dynamic-fba/sbml-models/iJR904.xml.gz"
declare -x data_dir="./data/simulated_data_sigma_0_01_25starts_L-BFGS-B.csv"
declare -x examplename="example1"
declare -x dir_to="ex1_synthetic" 	#foldername
declare -x costfunct="LS"	# "LS" or "NLLH"
declare -a lb=(-4 -1 -4 -1)
declare -a ub=(0 2 0 2)


#for cell_type in ${cellNames[@]}; do
#for gene_index in {0..1}; do
#  start_index=$((gene_index*gene_range))

echo $nstarts
echo $opt_method
echo $parallel
echo $model_dir
echo $data_dir
echo $examplename
echo $dir_to
echo ${lb[@]}
echo ${ub[@]}
 
sbatch /home/edudkin/dfba/dfba_tests/script_sbatched_dfba_all_models.sh $nstarts $opt_method $parallel $model_dir $data_dir $examplename $dir_to $costfunct ${lb[@]} ${ub[@]}


echo qsub done

