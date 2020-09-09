#!/bin/bash
#####  Constructed by HPC everywhere #####
#SBATCH --mail-user=lauraliu@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=3:00:00
#SBATCH --partition=general
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=hh-lik
#SBATCH --array=2-5

######  Module commands #####



######  Job commands go below this line #####
cd /N/u/lauraliu/BigRed3/het_agents_bayes/program/
matlab -r check_likelihood_hh_array
