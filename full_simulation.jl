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
    if !isdir("tmp")
        mkdir("tmp")
    end
    @everywhere include(%%simulation_functions);
    all_params=jldopen(args["params"])
    run_params=all_params["queue"][args["worker"]][1,1]
    run_params["jobid"]=args["worker"]
    kick_params=copy(run_params) # Absolutely DO NOT REMOVE copy() under any circumstances - You will bork your simulation parameters by not copying the dict. 
    kick_params["runtime"]=kick_params["timestep"]
    kick_params["saveat"]=kick_params["timestep"]
    println("Running a single time step to speed things up.")
    driver(kick_params)
    println("Now running the requested simulation.")
    results=driver(run_params)
    jldsave("tmp/"*ENV["SLURM_JOBID"]*"_"*string(args["worker"])*".jld2"; results=results)
    exit()
end
using NQCDynamics
using Unitful
using DiffEqBase
include(%%simulation_scripts)
######### SIMULATION PARAMETER BLOCK

# Set fixed parameters
fixed_parameters=Dict(
    "model" => SurfgenAmmoniaModel(),
    "atoms" => Atoms([:N, :H, :H, :H]), # Atomic structure to consider
    "initial_positions" => [0.0 1.91569 -0.901583 -0.901583; 0.0 -0.0 -1.62123 1.62123; 0.0 0.0 0.0 0.0], # Initial ammonia positions (transposed with atoms columns, cartesian coordinate rows)
    "trajectories" => 10, #Trajectories for RPMD
    "runtime" => 5.0u"ps", # Runtime for RPMD
    "timestep" => 0.5u"fs", # Timestep 
    "adiabatic_state_number" => 2, # Run in the first excited state
    # "equilibration_time" => 50u"fs",
    "ensemble_algorithm" => EnsembleThreads(), # Parallelisation method to use
    "saveat" => 10u"fs", # Sampling time interval for RPSH
    "gamma" => 1.0, # Determines the initial distribution to load from. (Shouldn't matter too much, since we have enough time to equilibrate. )
    "name" => "rpsh",
)


# Set variables

variables=Dict(
    # "gamma" => collect(0.0:0.1:2.0),
    "beads" => [4,8,16,32],
    "temperature" => collect(Float64, 100:100:1000).*u"K",
)
# Total trajectories to run: 
# Estimated time per trajectory: 3


function input_preprocess(input_dict)
    input_dict["nuclear_distrn"]="../rpmd_initial_surfgen/initial/rpmd_"*string(ustrip(input_dict["temperature"]))*"K_"*string(input_dict["beads"])*"beads_1ps.jld2"
end

if isfile(args["params"])
    @warn "Simulation queue already exists and will not be regenerated. Delete "*args["params"]*" to regenerate the queue."
    exit()
else
    job_queue=build_job_queue(fixed_parameters, variables, input_preprocess)
    serialise_queue!(job_queue)
    @info "Generated simulation queue at "*args["params"]
    exit()
end

