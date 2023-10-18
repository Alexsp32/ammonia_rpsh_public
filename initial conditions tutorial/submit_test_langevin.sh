#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=3700
#SBATCH --cpus-per-task=1
#SBATCH --output=/storage/mssgwp_grp/msrkhg/ammonia_rpsh/new_initial_conditions/langevin_test.log
#SBATCH --chdir=/storage/mssgwp_grp/msrkhg/ammonia_rpsh/new_initial_conditions

# Executable to call
x="julia +release"

# Path to input script
f=langevin_initial_conditions.jl

# Path to output file
# Parallelisation options

PARALLEL_OPTIONS="-N 1 --delay=2 -j ${SLURM_NTASKS} --joblog parallel_MC.log"
PARALLEL_EXEC="$x $f --params="langevin_settings.jld2" --worker"

# Any further definitions
THREADING=1

export JULIA_PROJECT=/storage/mssgwp_grp/msrkhg/ammonia_rpsh/
export MKL_DYNAMIC=FALSE # ToDo: Test if this makes a difference here
export OMP_NUM_THREADS=$THREADING
export MKL_NUM_THREADS=$THREADING

echo "(1/2) Creating job queue from $f:";
$x --project=$JULIA_PROJECT $f --params="langevin_settings.jld2"

echo "(2/2) Running parallel jobs from $f:"
$x --project=$JULIA_PROJECT $f --params="langevin_settings.jld2" --worker 1