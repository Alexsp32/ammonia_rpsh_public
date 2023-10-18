# ammonia_rpsh repo
This repository contains scripts for the simulation of hydrogen roaming in photo-excited ammonia. 

## Files
**initial_distribution.jl**: Contains routines for different dynamics methods, each of which takes a Dictionary of parameters. 

**full_simulation.jl**: Defines a dictionary of simulation parameters (or parameter screenings) and can be given arguments to run certain sub-jobs in the simulation defined here. This is what you will want to launch using SLURM. 

**jobscipt_slurm.sh**: Example submission script for a larger simulation on Taskfarm. 

**analysis_functions.jl**: Old script, contains some functions useful for analysis that might be broken by now. 


## Julia Environment
**Project.toml, Manifest.toml**: These define the Julia environment required to use the Julia scripts included in the package. 

Activate the environment in a Julia Terminal using `]activate $FOLDER_WITH_PROJECT_TOML_IN_IT`. 
Then download and precompile all necessary packages with `instantiate`. 
You may need to `rm AmmoniaSingletsModel` and re-install it with `add $GIT_REPO_LINK`. 
Make sure you have the NQCRegistry installed or NQCDynamics and associated packages will not be found. 
