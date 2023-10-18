using AmmoniaSingletsModel
using ArgParse
using Distributed
using JLD2
using Libdl: DL_LOAD_PATH
using ClusterScripts

function argparse()
    s=ArgParseSettings()
    @add_arg_table s begin
        "--worker"
            help="Julia processes called by parallel should use this option to indicate they need to run something."
            arg_type=Int
            default=0
        "--params"
            help="The location of the queue file. "
            default="simulation_parameters.jld2"
    end
    return parse_args(s)
end
args=argparse()
# Execution block for a worker process. 
if args["worker"]!=0
    if !isdir("$(ENV["SLURM_JOBID"])")
        mkdir("$(ENV["SLURM_JOBID"])")
    end
    @everywhere include("/storage/mssgwp_grp/msrkhg/ammonia_rpsh/script_repo/initial_distribution.jl");
    all_params=jldopen(args["params"])
    run_params=all_params["queue"][args["worker"]][1,1]
    run_params["jobid"]=args["worker"]
    kick_params=copy(run_params) # Absolutely DO NOT REMOVE copy() under any circumstances - You will bork your simulation parameters by not copying the dict. 
    kick_params["MC_steps"]=1
    println("Running a single time step to speed things up.")
    driver(kick_params)
    println("Now running the requested simulation.")
    results=driver(run_params)
    jldsave(ENV["SLURM_JOBID"]*"/"*string(args["worker"])*".jld2"; results=results)
    exit()
end
using NQCDynamics
using Unitful
using DiffEqBase
include("/storage/mssgwp_grp/msrkhg/ammonia_rpsh/script_repo/initial_distribution.jl")
######### SIMULATION PARAMETER BLOCK

# Set fixed parameters
fixed_parameters=Dict(
    "model" => "AmmoniaSingletsModel",
    "atoms" => Atoms([:N, :H, :H, :H]), # Atomic structure to consider
    "initial_positions" => [0.0 1.91569 -0.901583 -0.901583; 0.0 -0.0 -1.62123 1.62123; 0.0 0.0 0.0 0.0], # Initial ammonia positions (transposed with atoms columns, cartesian coordinate rows)
    "trajectories" => 1, #  
    "adiabatic_state_number" => 1, # Run in the first excited state
    "ensemble_algorithm" => EnsembleSerial(), # Parallelisation method to use
    "md_type" => "MC",
    "MC_move_ratio" => 0.65,
    "MC_stepsize_N" => 0.02,
    "MC_stepsize_H" => 0.1,
    "MC_steps" => 8e5,
)


# Set variables

variables=Dict(
    "temperature" => collect(Float64, 100:50:1000).*u"K",
)
# Total trajectories to run: 18*18*4=1296

if isfile(args["params"])
    @warn "Simulation queue already exists and will not be regenerated. Delete "*args["params"]*" to regenerate the queue."
    exit()
else
    job_queue=build_job_queue(fixed_parameters, variables)
    serialise_queue!(job_queue; filename=args["params"])
    @info "Generated simulation queue at "*args["params"]
    exit()
end

