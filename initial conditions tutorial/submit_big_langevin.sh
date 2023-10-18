#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=3700
#SBATCH --cpus-per-task=1

module load GCCcore/11.2.0 parallel/20210722;

# The script takes three arguments: 1:2:3 define the start, step and end of the range of things to run.
# Build the range of simulation queue items to run 
job_array=$(seq $1 $2 $3)

# Executable to call
x="julia +release"

# Path to input script
f=langevin_initial_conditions.jl

# Path to output file
# Parallelisation options

test ! $SLURM_NTASKS && export SLURM_NTASKS=$SLURM_NTASKS_PER_NODE
echo $SLURM_NTASKS
echo $SLURM_NTASKS_PER_NODE
echo $SLURM_JOBID

SRUN_OPTIONS="-N 1 -n 1 --exclusive"
PARALLEL_OPTIONS="-N 1 --delay=2 -j ${SLURM_NTASKS} --joblog parallel.log"
PARALLEL_EXEC="$x $f --params="langevin_settings.jld2" --worker {1}"

# Any further definitions
THREADING=1

export JULIA_PROJECT=.
export MKL_DYNAMIC=FALSE # ToDo: Test if this makes a difference here
export OMP_NUM_THREADS=$THREADING
export MKL_NUM_THREADS=$THREADING

echo "(1/2) Creating job queue from $f:";
$x --project=$JULIA_PROJECT $f --params="langevin_settings.jld2"



echo "(3/3) Running parallel jobs from $f:"
parallel $PARALLEL_OPTIONS srun $SRUN_OPTIONS $PARALLEL_EXEC ::: ${job_array[@]}