module load GCCcore/11.2.0 parallel/20210722;

# The script takes three arguments: 1:2:3 define the start, step and end of the range of things to run.
# Build the range of simulation queue items to run 
job_array=$(seq 1 1 114)

# Executable to call
x="julia +release"

# Path to input script
f=mc_rpmd_initial_conditions.jl

# Path to output file
# Parallelisation options
export SLURM_JOBID=rpmontecarlo-ded2

SRUN_OPTIONS="-N 1 -n 1 --exclusive"
PARALLEL_OPTIONS="-N 1 --delay=2 -j 24 --joblog parallel_rpmc.log"
PARALLEL_EXEC="$x $f --params="rpmc_settings.jld2" --worker {1}"

# Any further definitions
THREADING=1

export JULIA_PROJECT=.
export MKL_DYNAMIC=FALSE # ToDo: Test if this makes a difference here
export OMP_NUM_THREADS=$THREADING
export MKL_NUM_THREADS=$THREADING

echo "(1/2) Creating job queue from $f:";
$x --project=$JULIA_PROJECT $f --params="rpmc_settings.jld2"



echo "(3/3) Running parallel jobs from $f:"
parallel $PARALLEL_OPTIONS $PARALLEL_EXEC ::: ${job_array[@]}