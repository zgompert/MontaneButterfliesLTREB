#!/bin/sh 
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=stan22
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load R/4.1.1

cd /uufs/chpc.utah.edu/common/home/gompert-group3/projects/butterfly_predict

perl forkRunModel.pl

#R CMD BATCH --no-save --no-restore run_predict_2022.R
