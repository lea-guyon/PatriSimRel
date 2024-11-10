#!/bin/bash

STARTTIME=$(date +%s)

#############################
######## PARAMETERS #########
#############################

dir=~ # path

chr_size="c(1e6, 1e6, 1e4)"
random_fission=false
transmission="full"
fission_threshold=150
pM=0
nb_villages=5
nb_group=1 # nb of lineages -> doesn't make sense for bilateral descent but usefull to normalize villages' sizes
K=300 # carrying capacity per lineage
declare -i K_total=$nb_villages*$nb_group*$K # total carrying capacity
mf=0.5
mm=0
sigma=0.1 # parameter of the gamma law used to draw growth rates
growth_rate=0 #1.004 or 1.01 # growth rate of villages and outgroup, if 1 : population has a constant size
sample_size=20
pseudohap=False
nbsimu=100 # nb of simulations
cores=30

############################
####### SIMULATIONS ########
############################

echo "Starting simulations"

cd $dir

vl="F"
rf="F"
descent="unilineal"
descent_rule="patrilineal"
nameDir="patrilineal_2_bilateral"
path=$descent/extended/r=$growth_rate/sigma=$sigma/FT=$fission_threshold
rm -rf $dir/Tables/Simulations_metrics/$path/$nameDir
mkdir -p $dir/Tables/Simulations_metrics/$path/$nameDir
echo "'Replicat'    'Generation'    'N_ind'    'Nb_of_fissions'    'Nb_of_extinctions'    'Nb_of_groups'    'fathers'    'Nb_indiv_per_group'    'var_nb_ind_per_group'    'Nb_women_per_group'    'mothers'    'failed_couples'    'singleInd'	'Nb_children_per_couple'    'var_nb_children_per_couple'    'mean_group_depth'    'var_group_depth'    'mean_migrant_ratio'    'var_migrant_ratio' 'meanFissionTime' 'varFissionTime'" > $dir/Tables/Simulations_metrics/$path/$nameDir/metrics.txt
echo "'Replicat'	'Generation'	'Step'	'nMales'	'nInds'	'sexRatio'" > $dir/Tables/Simulations_metrics/$path/$nameDir/sexRatio.txt

rm -rf Simulations_folders/$path/$nameDir
mkdir -p Simulations_folders/$path/$nameDir

cd Simulations_folders/$path/$nameDir/

## Replace parameters in the slim file ##

cat $dir/SLiM_scripts/patrilineal2bilateral_relatedness.slim | sed "s/bash_wd/${dir}/g;s/bash_Num_villages/${nb_villages}/g;s/bash_carrying_capacity/${K}/g;s/bash_chr_size/${chr_size}/g;s/bash_random_fission/${rf}/g;s/bash_fission_threshold/${fission_threshold}/g;s/bash_pM/${pM}/g;s/bash_descent_rule/${descent_rule}/g;s/bash_residence_rule/${residence_rule}/g;s/bash_sigma/${sigma}/g;s/bash_growth_rate/${growth_rate}/g;s/bash_transmission/${transmission}/g;s/bash_dir_table/${nameDir}/g" > "islandmodel.slim"

## Create a new file for each simulation ##
for i in $(seq 1 1 $nbsimu)
do	
	mkdir $dir/simulations/$path/$nameDir/$i
done

cd $dir/simulations/$path/$nameDir/
if $burnin; then
	echo "burnin"
	for i in $(seq 1 1 $nbsimu)
	do	
		cd $dir/simulations/$path/$nameDir/$i
		cat $dir/SLiM_scripts/burnin.slim | sed "s/bash_wd/${dir}/g;s/bash_Num_villages/${nb_villages}/g;s/bash_nGroupsPerVillage/${nb_groups}/g;s/bash_total_carrying_capacity/${K_total}/g;s/bash_carrying_capacity/${K}/g;s/bash_chr_size/${chr_size}/g;s/bash_descent/${descent}/g;s/bash_Num_replicat/${i}/g" > "burnin_${i}.slim"
		echo "slim $i/burnin_${i}.slim"
		cd ..
	done > launcher.txt
	parallel -a launcher.txt -j $cores
fi

for i in $(seq 1 1 $nbsimu)
do
	cd $dir/Simulations_folders/$path/$nameDir/$i
	cat ../islandmodel.slim | sed "s/bash_Num_replicat/${i}/g" > "islandmodel_${i}.slim"

	echo "echo $i"
    echo "slim $i/islandmodel_${i}.slim > $i/outputSlim${i}.slim"
		
	cd ..
done > launcher.txt
parallel -a launcher.txt -j $cores > outputErrors.slim

ENDTIME=$(date +%s)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"

##################################
######### SUBSET TREES ###########
##################################
cd $dir/Simulations_folders/$path/$nameDir/
echo "Subset trees"
STARTTIME2=$(date +%s)

cd $dir/Simulations_folders/$path/$nameDir/

gen_patrilineal=$(seq 0 20 100)
gen_patrilocal=$(seq 120 20 200)

for i in $(seq 1 1 $nbsimu)
do
	echo "python $dir/Subset_trees/subset_trees_bilateral_descent_relatedness.py -s $dir/Simulations_folders/$path/$nameDir/ -rep $i -g $gen_patrilocal --sample-size $sample_size -K $K_total -o $dir/Simulations_folders/$path/$nameDir/$i/ > $i/outputPy${i}.txt"
	echo "python $dir/Subset_trees/subset_trees_unilineal_descent_relatedness.py -s $dir/Simulations_folders/$path/$nameDir/ -rep $i -g $gen_patrilineal --sample-size $sample_size -K $K_total -og False -d $descent_rule -o $dir/Simulations_folders/$path/$nameDir/$i/ -t $dir/Tables/Simulations_metrics/$path/$nameDir > $i/outputPy${i}.txt"
done > launcher.txt
parallel -a launcher.txt -j $cores

ENDTIME=$(date +%s)
echo "It takes $(($ENDTIME - $STARTTIME2)) seconds to complete this task"

##################################
######## COMPUTE METRICS #########
##################################
echo "Compute metrics"
STARTTIME3=$(date +%s)

cd $dir/

rm -rf $dir/Tables/Pi/$path/$DirTable
mkdir -p $dir/Tables/Pi/$path/$DirTable/village
python ~/Documents/SLiM_model/Python_scripts/Pi_villages.py -s $dir/Simulations_folders/$path/$nameDir/village/ -c "1e6, 1e6, 1e6, 1e4" -ph $pseudohap -o $dir/Tables/Pi/$path/$DirTable/village/