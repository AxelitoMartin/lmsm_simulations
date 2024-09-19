#!/bin/bash

#SBATCH --partition=fn_short,cpu_short,fn_medium,cpu_medium
#SBATCH --tasks=5
#SBATCH --nodes=1
#SBATCH --tasks-per-node=5
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=axel.martin@nyulangone.org
#SBATCH --out=Mult_n1000.txt

module load r/4.1.2
Rscript R/Mult_n1000.R
