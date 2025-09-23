#!/bin/bash
#SBATCH --job-name=rgcca # Set a partition or queue to submit to
#SBATCH --partition=long # Set a partition or queue to submit to
#SBATCH --account=generic # Set a user account to submit as
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adria.seto@isglobal.org # Where to send mail events
#SBATCH --ntasks=1 # Run a single tasks
#SBATCH --cpus-per-task=24 #Use 4 cpus for each task
#SBATCH --mem=64gb # Job memory request
#SBATCH --output=/home/isglobal.lan/aseto/Denoising_Slurm_%j.out
#SBATCH --error=/home/isglobal.lan/aseto/Denoising_Slurm_%j.err

# Clear the environment from any previously loaded modules
module purge > /dev/null 2>&1

# Load the bashrc and conda environment needed
source ~/.bashrc
source activate TFM_env

# And finally run the job (in this case it’s a R script called from bash)
Rscript /home/isglobal.lan/aseto/RGCCA_methyl_v1.R