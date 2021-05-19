#!/bin/bash

#SBATCH --job-name="A_snk"
#SBATCH --partition=long
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH --ntasks=1

snakemake --keep-going --cluster "sbatch -p medium --mem 60000 -N 1 --cpus-per-task 16 -J Assemblies" --jobs 10