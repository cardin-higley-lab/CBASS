#!/bin/bash
#SBATCH --partition=scavenge
#SBATCH --requeue
#SBATCH --ntasks=1 --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6000 
#SBATCH --time=2:00:00

<<COMMENT
This script creates and assigns a job to a node in the cluster with the specifications given above (currently using one node with 4 cpus and 24GB of RAM). Also, if intending 
to use the embedding and plotting routines, the python enviroment needs to be activated in addition to the Matlab.
All the parameters in this script are specified by the scheduler (scheduler_CBASS.sh)
COMMENT

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Loading Matlab, conda and our customized environment to the node.
module load MATLAB/2020a
module load miniconda
source activate gammaBouts_env

# matlab -nodisplay -nosplash -nodesktop -r "root_path='"$FOLDER"', cd(root_path); cBAND_LABEL = '"${1}"', cBAND = eval('"${2}"'), cSTATE_LBL = '"${3}"', cFORMAT = '"${4}"', bl1ZScore = "${5}", run('CBASS_Call_Main_AF.m') "
matlab -nodisplay -nosplash -nodesktop -r "root_path='"$FOLDER"', cd(root_path); cBAND_LABEL = '"${1}"', cBAND = eval('"${2}"'), cSTATE_LBL = '"${3}"', cFORMAT = '"${4}"', bl1ZScore = "${5}", run('"$ScriptName"') "
