#!/bin/bash
#SBATCH --name="rpsh"
#SBATCH --nodes=%%nodes
#SBATCH --ntasks-per-node=%%tasks
#SBATCH --time=%%time
#SBATCH --mem-per-cpu=%%mem
#SBATCH --cpus-per-task=%%cpus
#SBATCH --output=%%output
#SBATCH --chdir=%%fullpath

# Path to executable
x=/home/chem/msrkhg/bin/julia

# Path to input script
f=full_simulation.jl

# Path to output file
o=%%file_output


# Any further definitions
export JULIA_PROJECT=%%project


# Run everything
$x $f --julia-project $JULIA_PROJECT --slurm
