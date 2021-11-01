#!/bin/bash
#SBATCH -p scavenge                # partition (queue)
#SBATCH -n 1                      # number of cores
#SBATCH --output=scheduler_stdout.txt
#SBATCH --error=setup_stderr.txt
#SBATCH --job-name=gamma_bouts

export FOLDER=/gpfs/ysm/home/ahf38/Documents/gamma_bouts/Figures_and_Sandbox
FOLDER_LOGS="$FOLDER/Logs"
# FOLDER_LOGS=/gpfs/ysm/home/ahf38/Documents/gamma_bouts/Demo/Logs
# export FOLDER=${1}
folder_name=$(basename $FOLDER)

if [ ! -d "$FOLDER_LOGS" ]; then
    echo "creating $FOLDER_LOGS."
    mkdir -p $FOLDER_LOGS
fi

# cd ${FOLDER_LOGS}
# pwd


# grab the files, and export it so the 'child' sbatch jobs can access it
# ls *.WAV *.wav> filenames.txt

gammaBouts="gamma_bouts_$folder_name"
echo "running job $gammaBouts"
rm -f joblist.txt

# Declare an array of string with type
declare -a StringBandLabel=("Beta" "Gamma")
declare -a StringBands=("[15,30]" "[30,80]")
declare -a StringStateLabel=("Stim" "Running");  
declare -a StringcFormat=("complex" "polar");  
declare -a Stringcbl1ZScore=("true" "false");  
 
# Iterate the string array using for loop
INDEX=0
for bandLabel in ${StringBandLabel[@]}; do
    for stateLabel in ${StringStateLabel[@]}; do
        for cFormat in ${StringcFormat[@]}; do
        	for bl1ZScore in ${Stringcbl1ZScore[@]}; do
#         echo $bandLabel
	            OUT0=$(sbatch -p scavenge --requeue --chdir=${FOLDER_LOGS} --parsable -J $bandLabel$stateLabel$cFormat$bl1ZScore -o $bandLabel$stateLabel$cFormat$bl1ZScore.stdout.txt -e $bandLabel$stateLabel$cFormat$bl1ZScore.stderr.txt cluster_test.sh ${bandLabel} ${StringBands[$INDEX]} ${stateLabel} $cFormat $bl1ZScore ${folder_name})
	        done
        done
    done
    let INDEX=${INDEX}+1
done


# while read p; do
  

  sleep 1 # pause to be kind to the scheduler
# done <filenames.txt

ID_LIST=$(squeue -u ahf38|  awk '{ printf (":%i", $1 )}')
ID_LIST=${ID_LIST:3}
#C=("${A[@]:1}")
echo $ID_LIST



