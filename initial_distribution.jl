using AmmoniaSingletsModel
using Combinatorics
using DiffEqBase
using Distributed
using JLD2
using LinearAlgebra
using NQCDynamics
using RingPolymerArrays
using Unitful
using UnitfulAtomic

#! Utility functions

"""
    distribution_from_single_datapoint((trajectories, parameters), electronic_state, filename::String)

Creates a [DynamicalDistribution](https://nqcd.github.io/NQCDynamics.jl/stable/NQCDistributions/overview/) of initial conditions from results generated using `langevin_ammonia, langevin_rpmd`, ...

Input a (results, parameters) tuple, the electronic state in question (e.g. `PureState(2, Adiabatic())`) and the filename to save it as. 

`results`: Vector containing output dictionaries of `NQCDynamics.run_dynamics`
`parameters`: Dictionary containing simulation parameters. 
"""
function distribution_from_single_datapoint((trajectories, parameters), electronic_state, filename::String)
    if parameters["md_type"]=="langevin"
        # Package results from Langevin dynamics into an initial distribution, using both velocities and positions. 
        # This slightly stupid way of defining my "long list of ..." matrices improves compatibility with different NQCD methods. 
        all_positions=map(get_positions, trajectories[1][:OutputDynamicsVariables])
        all_velocities=map(get_velocities, trajectories[1][:OutputDynamicsVariables])
        for trajectory in trajectories[2:end]
            append!(all_positions, map(get_positions, trajectory[:OutputDynamicsVariables]))
            append!(all_velocities, map(get_velocities, trajectory[:OutputDynamicsVariables]))
        end
        distribution=NQCDistributions.DynamicalDistribution(all_velocities, all_positions, size(all_positions[1]))
    elseif occursin("rpmd", parameters["md_type"]) # catch-all for rpmd initialisation methods. 
        all_positions=map(get_positions, trajectories[1][:OutputDynamicsVariables])
        for trajectory in trajectories[2:end]
            append!(all_positions, map(get_positions, trajectory[:OutputDynamicsVariables]))
        end
        velocities=NQCDistributions.RingPolymerWrapper(VelocityBoltzmann(parameters["temperature"], parameters["atoms"].masses, size(all_positions)[1:2]), size(all_positions[1]), Int[])
        distribution=DynamicalDistribution(velocities, all_positions, size(all_positions))
    elseif parameters["md_type"]=="MC"
        # This slightly stupid way of defining my "long list of ..." matrices improves compatibility with different NQCD methods. 
        all_positions=trajectories[1][:OutputPosition]
        for chain in trajectories[2:end]
            append!(all_positions, chain[:OutputPosition])
        end
        velocities=VelocityBoltzmann(parameters["temperature"], parameters["atoms"].masses, size(all_positions[1]))
        distribution=DynamicalDistribution(velocities, all_positions, size(all_positions[1]))
    end
    distribution*=electronic_state
    jldsave(filename; distribution)
end

"""
    generate_starting_configuration(parameters)

This function is used to avoid loading a large initial distribution every time by sampling the initial distribution beforehand. 
Use in the simulation queue postprocessing step. 

"""
function generate_starting_configuration!(parameters)
    if haskey(parameters, "initial_configuration")
        # An initial configuration was already sampled, pass through. 
        return parameters
    elseif haskey(parameters, "initial_distribution")
        # We have a DynamicalDistribution that needs sampling. 
        distribution=jldopen(parameters["initial_distribution"])["distribution"]
        if hasfield(typeof(distribution), :nuclear)
            # This is a product distribution, so we need to randomly sample the nuclear part
            config=rand(distribution.nuclear)
            elec=PureState(distribution.electronic.state, distribution.electronic.statetype)
        else
            # Or we just have a nuclear distribution, which we can sample directly
            config=rand(distribution)
            elec=nothing
        end
        # Put everything into a Dict
        parameters["initial_configuration"]=Dict(:v => config.v, :r => config.r, :elec => elec)
        return parameters
    end
end

"""
    generate_dynamicsvariables(parameters, sim::Simulation)

Function to load the initial configuration as a DynamicsVariables object for NQCDynamics. 
This is in a separate function since all MD methods need it. 
"""
function generate_dynamicsvariables(parameters, sim::Simulation)
    if !haskey(parameters, "initial_configuration")
        parameters=generate_starting_configuration!(parameters)
    end
    initial_config=parameters["initial_configuration"]
    if isa(initial_config[:elec], Nothing)
        return DynamicsVariables(sim, initial_config[:v], initial_config[:r])
    else
        return DynamicsVariables(sim, initial_config[:v], initial_config[:r], initial_config[:elec])
    end
end

"""
    OutputInitialConfiguration(sol, u)

Outputs the initial DynamicsVariables of a simulation so we have a record of which configuration a trajectory started with. 
(Useful for random sampling from distributions)
"""
function OutputInitialConfiguration(sol, u)
    return first(sol.u)
end

struct DesorptionTerminator
    dissociation_distance::Float64
end

"""
    (desorption_terminator::DesorptionTerminator)(u, t, integrator)::Bool

Checks if any atoms are more than `dissociation_distance` (in atomic units!) apart to terminate a simulation. 
"""
function (desorption_terminator::DesorptionTerminator)(u, t, integrator)::Bool
    positions=DynamicsUtils.get_positions(u)
    distance(ind1, ind2)=norm(positions[:,ind2]-positions[:,ind1])
    combinations=collect(multiset_combinations(1:size(positions)[2], 2))
    distances=map(x->distance(x...), combinations)
    if length(distances[distances.>desorption_terminator.dissociation_distance])>0
        return true
    else
        return false
    end
end

"""
    driver(params)

This function is called by your simulation script to run the correct function based on the input dictionary provided. 
"""
function driver(params)
    # Execution block, run correct function for the simulation job. 
    if params["md_type"]=="langevin"
        results=langevin_ammonia(params)
    elseif params["md_type"]=="langevin_rpmd"
        results=langevin_rpmd(params)
    elseif params["md_type"]=="MC"
        results=thermal_montecarlo(params)
    elseif params["md_type"]=="rpsh"
        results=run_rpsh(params)
    elseif params["md_type"]=="fssh"
        results=run_fssh(params)
    elseif params["md_type"]=="classical"
        results=run_classical_adiabatic_md(params)
    elseif params["md_type"]=="classical_rpmd"
        results=run_adiabatic_rpmd(params)
    end
    return results
end 

#! MD functions

"""
    activate_pes_model(input_params)

AmmoniaSingletsModel.jl tries to activate the fortran wrapper for the model as soon as it's called for the first time. 
This function pushes loading the model as late in the simulation process to avoid crashes on incompatible systems. 

I really couldn't be bothered to get AmmoniaSingletsModel.jl working on macOS and this is my half-arsed attempt at avoiding issues when I don't need the model. 
"""
function activate_pes_model(input_params)
    if input_params["model"]=="AmmoniaSingletsModel"
        return SurfgenAmmoniaModel()
    elseif input_params["model"]=="AlternateAmmoniaModel"
        return AlternateAmmoniaModel()
    end
end

"""
    langevin_ammonia(input_params)

Runs [Langevin dynamics](https://nqcd.github.io/NQCDynamics.jl/stable/dynamicssimulations/dynamicsmethods/langevin/) using input dictionary. 
"""
function langevin_ammonia(input_params)
    # Fix the model in the ground state. 
    input_model=activate_pes_model(input_params)
    ground_state_model=AdiabaticStateSelector(input_model, 1)
    # Initialise simulation object
    sim=Simulation{Langevin}(
        get!(input_params, "atoms", Atoms([:N, :H, :H, :H])),
        ground_state_model,
        temperature=get!(input_params, "temperature", 300u"K"),
        γ=get!(input_params, "gamma", 1.0)
    )
    # Set initial positions and zero velocities. 
    initial=DynamicsVariables(sim, zeros(size(sim)), input_params["initial_positions"])
    # Decide when run_dynamics should save
    # Save at every timestep if not already differently defined.
    if !haskey(input_params, "saveat")
        input_params["saveat"]=get!(input_params, "timestep", 0.1u"fs")
    end
    # Allow for equilibration time if defined. 
    if get!(input_params, "equilibration_time", 0.0)==0.0
        saveat_arg=input_params["saveat"]
    else
        saveat_arg=map(austrip, collect(input_params["equilibration_time"]:input_params["timestep"]:input_params["runtime"]))
    end
    # Now run dynamics
    traj=run_dynamics(
        sim, 
        (0.0u"fs", get!(input_params, "runtime", 0.3u"ps")),
        initial,
        dt=get!(input_params, "timestep", 0.1u"fs"), # Need to double-check because previous if-statement doesn't cover every possibility
        trajectories=get!(input_params, "trajectories", 1),
        output=(OutputDynamicsVariables, OutputPotentialEnergy, OutputKineticEnergy),
        saveat=saveat_arg, 
        ensemble_algorithm=get!(input_params, "ensemble_algorithm", EnsembleSerial()), # Parallelise if possible
    )
    # Return single trajectories as size 1 vector for a consistent output format with ensemble simulations. 
    if input_params["trajectories"]==1
        traj=[traj]
    end
    return (traj, input_params) # Return modified input parameters as well, to give additional information. 
end

"""
    langevin_rpmd(parameters)

Runs ring-polymer Langevin dynamics. 
"""
function langevin_rpmd(parameters)
    # Fix the model in the ground state. 
    input_model=activate_pes_model(parameters)
    ground_state_model=AdiabaticStateSelector(input_model, 1)
    # Initialise simulation object
    sim=RingPolymerSimulation{ThermalLangevin}(
        get!(parameters, "atoms", Atoms([:N, :H, :H, :H])),
        ground_state_model,
        parameters["beads"];
        temperature=get!(parameters, "temperature", 300u"K"),
        γ=get!(parameters, "gamma", 1.0)
    )
    # Set initial positions and zero velocities. 
    initial=DynamicsVariables(sim, RingPolymerArray(zeros(size(sim))), RingPolymerArray(repeat(parameters["initial_positions"], 1,1,parameters["beads"])))
    # Decide when run_dynamics should save
    # Save at every timestep if not already differently defined.
    if !haskey(parameters, "saveat")
        parameters["saveat"]=get!(parameters, "timestep", 0.1u"fs")
    end
    # Allow for equilibration time if defined. 
    if get!(parameters, "equilibration_time", 0.0)==0.0
        saveat_arg=parameters["saveat"]
    else
        saveat_arg=map(austrip, collect(parameters["equilibration_time"]:parameters["timestep"]:parameters["runtime"]))
    end
    # Now run dynamics
    traj=run_dynamics(
        sim, 
        (0.0u"fs", get!(parameters, "runtime", 0.3u"ps")),
        initial,
        dt=get!(parameters, "timestep", 0.1u"fs"), # Need to double-check because previous if-statement doesn't cover every possibility
        trajectories=get!(parameters, "trajectories", 1),
        output=(OutputDynamicsVariables, OutputPotentialEnergy, OutputKineticEnergy, OutputCentroidPosition, OutputCentroidVelocity),
        saveat=saveat_arg, 
        ensemble_algorithm=get!(parameters, "ensemble_algorithm", EnsembleSerial()), # Parallelise if possible
    )
    # Return single trajectories as size 1 vector for a consistent output format with ensemble simulations. 
    if parameters["trajectories"]==1
        traj=[traj]
    end
    return (traj, parameters) # Return modified input parameters as well, to give additional information.
end

function thermal_montecarlo(parameters)
    # Fix the model in the ground state. 
    input_model=activate_pes_model(parameters)
    ground_state_model=AdiabaticStateSelector(input_model, 1)
    # Initialise simulation object
    if !haskey(parameters, "beads")
        # No beads parameters --> classical simulation
        sim=Simulation{Classical}(
            get!(parameters, "atoms", Atoms([:N, :H, :H, :H])),
            ground_state_model,
            temperature=get!(parameters, "temperature", 300u"K"),
        )
    else
        # Beads parameter --> Ring polymer simulation
        sim=RingPolymerSimulation{Classical}(
            get!(parameters, "atoms", Atoms([:N, :H, :H, :H])),
            ground_state_model,
            parameters["beads"];
            temperature=get!(parameters, "temperature", 300u"K"),
        )
    end
    # Set initial positions and zero velocities. 
    step_sizes=Dict(
        :H => get!(parameters, "MC_stepsize_H", 0.05),
        :N => get!(parameters, "MC_stepsize_N", 0.0), # Fix N in position to avoid translation of the entire molecule. 
    )
    rpmd_positions=isa(sim, RingPolymerSimulation) ? jldopen(parameters["initial_distribution"])["distribution"].nuclear : nothing
    mc_chain=InitialConditions.ThermalMonteCarlo.run_advancedmh_sampling(
        sim, # Simulation to run HMC sampling with
        isa(sim, RingPolymerSimulation) ? cat([rand(rpmd_positions).r for i in 1:parameters["beads"]]...;dims=3) : parameters["initial_positions"], # Initial configuration
        get!(parameters, "MC_steps", 10), # Number of MonteCarlo steps to do. 
        step_sizes; # Step sizes per atom species. 
        move_ratio=get!(parameters, "MC_move_ratio", 0.0)
    )
    potential_energy=NQCModels.potential.(ground_state_model, mc_chain)
    return([Dict(:OutputPosition => mc_chain, :OutputPotentialEnergy => potential_energy)], parameters)
end

"""
    run_rpsh(input_params)

Performs MD using [Ring-polymer Surface Hopping](https://nqcd.github.io/NQCDynamics.jl/stable/dynamicssimulations/dynamicsmethods/rpsh/) based on input parameters an initial conditions file. 

This is the ring-polymer version of Fewest-switches Surface Hopping - so using a ring polymer to further approximate quantum nuclear and electronic dynamics based on connected copies of the system in question. 
"""
function run_rpsh(input_params)
    model=activate_pes_model(input_params)
    sim=RingPolymerSimulation{FSSH}(
        get!(input_params, "atoms", Atoms([:N, :H, :H, :H])),
        model,
        get!(input_params, "beads", 4);
        temperature=get!(input_params, "temperature", 300u"K"),
    )
    # Load the initial conditions (positions, velocity, electronic state) from file. 
    distribution=jldopen(parameters["initial_conditions"])["distribution"]
    # Decide when run_dynamics should save
    # Save at every timestep if not already differently defined.
    if !haskey(input_params, "saveat")
        input_params["saveat"]=get!(input_params, "timestep", 0.1u"fs")
    end
    # Allow for equilibration time if defined. 
    if get!(input_params, "equilibration_time", 0.0)==0.0
        saveat_arg=input_params["saveat"]
    else
        saveat_arg=map(austrip, collect(input_params["equilibration_time"]:input_params["timestep"]:input_params["runtime"]))
    end
    results=run_dynamics(
        sim, 
        (0.0u"fs", runtime),
        distribution;
        dt=timestep,
        output=(OutputDynamicsVariables, OutputKineticEnergy, OutputPotentialEnergy, OutputDiscreteState, OutputDiabaticPopulation, OutputCentroidPosition, OutputCentroidVelocity),
        trajectories=get!(input_params, "trajectories", 1),
        ensemble_algorithm=get!(input_params, "ensemble_algorithm", EnsembleThreads()),
        saveat=saveat_arg,
        #reduction=get!(input_params, "reduction_type", MeanReduction())
    )
    return results, input_params
end

"""
    run_fssh(input_params)

Runs Tully's [Fewest-switches Surface Hopping](https://nqcd.github.io/NQCDynamics.jl/stable/dynamicssimulations/dynamicsmethods/fssh/). 

This mixed quantum-classical method treats electronic populations across states, whereby the potential energy surface experienced by nuclei is determined through stochastic hopping between states. 
"""
function run_fssh(input_params)
    model=activate_pes_model(input_params)
    sim=Simulation{FSSH}(
        get!(input_params, "atoms", Atoms([:N, :H, :H, :H])),
        model;
        temperature=get!(input_params, "temperature", 300u"K")
    )
    # Load initial distribution file to sample starting conditions from. 
    distribution=generate_dynamicsvariables(input_params, sim)
    # Decide when run_dynamics should save
    # Save at every timestep if not already differently defined.
    if !haskey(input_params, "saveat")
        input_params["saveat"]=get!(input_params, "timestep", 0.1u"fs")
    end
    # Allow for equilibration time if defined. 
    if get!(input_params, "equilibration_time", 0.0)==0.0
        saveat_arg=input_params["saveat"]
    else
        saveat_arg=map(austrip, collect(input_params["equilibration_time"]:input_params["timestep"]:input_params["runtime"]))
    end
    # And now run the dynamics
    results=run_dynamics(
        sim, 
        (0.0u"fs", input_params["runtime"]),
        distribution;
        output=(OutputDynamicsVariables, OutputPotentialEnergy, OutputKineticEnergy, OutputDiscreteState, OutputDiabaticPopulation),
        trajectories=get!(input_params, "trajectories", 1),
        ensemble_algorithm=get!(input_params, "ensemble_algorithm", EnsembleThreads()),
        saveat=saveat_arg,
        callback=DynamicsUtils.TerminatingCallback(DesorptionTerminator(input_params["dissociation_distance"])) # Stop simulation if something has dissociated
    )
    # Give the simulation object to enable Estimators.jl
    input_params["simulation"]=sim
    return results, input_params
end

"""
    run_classical_adiabatic_md(input_params)

Simulation function for [Classical MD](https://nqcd.github.io/NQCDynamics.jl/stable/dynamicssimulations/dynamicsmethods/classical/). 

These are entirely classical nuclear dynamics, bringing us closest to the Born-Oppenheimer approximation. 
Any PES (ground or excited state) will be explored classically, with no transition between states. 
"""
function run_classical_adiabatic_md(input_params)
    model=activate_pes_model(input_params)
    adiabatic_model=AdiabaticStateSelector(model, get!(input_params, "adiabatic_state", 1))
    sim=Simulation{Classical}(
        get!(input_params, "atoms", Atoms([:N, :H, :H, :H])),
        adiabatic_model;
        temperature=get!(input_params, "temperature", 300u"K")
    )
    # Load initial conditions from file. 
    distribution=jldopen(parameters["initial_conditions"])["distribution"]
    # Decide when run_dynamics should save
    # Save at every timestep if not already differently defined.
    if !haskey(input_params, "saveat")
        input_params["saveat"]=get!(input_params, "timestep", 0.1u"fs")
    end
    # Allow for equilibration time if defined. 
    if get!(input_params, "equilibration_time", 0.0)==0.0
        saveat_arg=input_params["saveat"]
    else
        saveat_arg=map(austrip, collect(input_params["equilibration_time"]:input_params["timestep"]:input_params["runtime"]))
    end
    run_dynamics(
        sim, 
        (0.0u"fs", runtime),
        distribution;
        dt=timestep,
        output=(OutputDynamicsVariables, OutputPotentialEnergy, OutputKineticEnergy),
        trajectories=get!(input_params, "trajectories", 1),
        ensemble_algorithm=get!(input_params, "ensemble_algorithm", EnsembleThreads()),
        saveat=saveat_arg,
    )
    return results, input_params
end

"""
    run_adiabatic_rpmd(input_params)

Simulation function for [Classical RPMD](https://nqcd.github.io/NQCDynamics.jl/stable/dynamicssimulations/dynamicsmethods/rpmd/). 

This would be the equivalent of running adiabatic classical MD, but in a ring polymer. 
"""
function run_adiabatic_rpmd(input_params)
    model=activate_pes_model(input_params)
    adiabatic_model=AdiabaticStateSelector(model, get!(input_params, "adiabatic_state_number", 1))
    sim=RingPolymerSimulation{Classical}(
        get!(input_params, "atoms", Atoms([:N, :H, :H, :H])),
        adiabatic_model,
        get!(input_params, "beads", 4);
        temperature=get!(input_params, "temperature", 300u"K"),
    )
    # Load initial distribution from file. 
    distribution=jldopen(parameters["initial_conditions"])["distribution"]
    # Decide when run_dynamics should save
    # Save at every timestep if not already differently defined.
    if !haskey(input_params, "saveat")
        input_params["saveat"]=get!(input_params, "timestep", 0.1u"fs")
    end
    # Allow for equilibration time if defined. 
    if get!(input_params, "equilibration_time", 0.0)==0.0
        saveat_arg=input_params["saveat"]
    else
        saveat_arg=map(austrip, collect(input_params["equilibration_time"]:input_params["timestep"]:input_params["runtime"]))
    end
    results=run_dynamics(
        sim, 
        (0.0u"fs", runtime),
        distribution;
        dt=timestep,
        output=(OutputDynamicsVariables, OutputCentroidPosition, OutputCentroidVelocity, OutputKineticEnergy, OutputPotentialEnergy),
        trajectories=get!(input_params, "trajectories", 1),
        ensemble_algorithm=get!(input_params, "ensemble_algorithm", EnsembleThreads()),
        saveat=saveat_arg,
    )
    return results, input_params
end
