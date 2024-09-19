#!/bin/bash

#SBATCH --partition=fn_medium,cpu_medium,fn_long,cpu_long
#SBATCH --tasks=10
#SBATCH --nodes=1
#SBATCH --tasks-per-node=10
#SBATCH --cpus-per-task=1
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=axel.martin@nyulangone.org
#SBATCH --out=Mult_n5000.txt

module load r/4.1.2
Rscript R/Mult_n5000.R
