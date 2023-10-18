using NQCDynamics
using Unitful
using JLD2
using DiffEqBase
using Measurements
using Plots
using Statistics
using ExtXYZ
using PyCall

function calculate_gyration_radius(results_file)
    #=
        radius_of_gyration essentially flattens a ring polymer snapshot with the RMSD from centroid positions for each atom and DOF. #!1
        Each element in the matrix needs to be averaged over the whole trajectory. #!2
        Each element in the trajectory average needs to be averaged #!3
        This total average could be turned into the total component over three DOFs (i.e. sqrt(x^2+y^2+z^2)) #!4
    =#
    if isa(results_file, String)
        results=jldopen(results_file, "r")["results"]
    elseif isa(results_file, Array)
        results=results_file
    end
    r_gyr_reduction(sim,x)=Estimators.radius_of_gyration(sim, x) #1 needs to be passed a matrix x of size sim
    trajectory_mean(x)=mean(x;dims=3) #2 needs to be passed a collection/generator of centroid positions
    all_trajectories_mean(x)=mean(x;dims=3) #3 needs to be passed a collection/generator of trajectory centroid averages
    dimensions_mean(x)=dropdims(sqrt.(sum(.^(x,2);dims=1));dims=1) #4 reduces into a single component amplitude
    #!5 Need to decide what I want to do. 
    r_gyr(datapoint, sim_params)=au_to_ang.(
        dimensions_mean(
            dropdims(all_trajectories_mean(
                cat([
                    dropdims(trajectory_mean(
                        cat([
                            r_gyr_reduction(
                            sim_params["simulation"], 
                            get_positions(frame)
                            ) 
                        for frame in trajectory[:OutputDynamicsVariables]
                        ]...;dims=3)
                    );dims=3) 
                    for trajectory in datapoint
                ]...;dims=3)
            );dims=3)
        ) # This looks absolutely horrible
    )
    results=map(x->r_gyr(x...), results)
    # results=au_to_ang.(results)
    return results
end

function NH_distance_centroid(datapoint)
    if isa(datapoint, String)
        datapoint=jldopen(datapoint, "r")["results"]
        println("Loading results. ")
    elseif isa(datapoint, Array)
        println("Results already loaded. ")
    end
    #=
    We are interested in the centroid positions, which are stored in ReplacementCentroidOutputPosition. 
    For each trajectory in a datapoint #!1, 
    for each frame in a trajectory #!2, 
    output the list of distances between N and the three H #!3
    =#
    nh_distance(snapshot)=au_to_ang.(sqrt.(dropdims(sum(.^(snapshot[:,2:end].-snapshot[:,1], 2);dims=1);dims=1))) #!3
    nh_distance_trajectory(traj)=cat([nh_distance(snapshot) for snapshot in traj[:ReplacementOutputCentroidPosition]]...;dims=2) #!2
    nh_distance_trajectories(data)=cat([nh_distance_trajectory(traj) for traj in data[1]]...;dims=3) #!1
    results=map(nh_distance_trajectories, datapoint)
    return results
end

function NH_distance_beads(datapoint)
    if isa(datapoint, String)
        datapoint=jldopen(datapoint, "r")["results"]
        println("Loading results. ")
    elseif isa(datapoint, Array)
        println("Results already loaded. ")
    end
    #=
    First, compress atomic coordinates into a matrix of NH distances #!1
    The atomic coordinates come from the positions of each bead, which are stored in a ndofs*natoms*nbeads matrix #!2
    Atomic coordinates are stored in :OutputDynamicsVariables of a trajectory #!3
    Each datapoint contains a vector of NQCD output Dicts #!4
    =# 
    convert_to_nh_distance(rp_snapshot)=sqrt.(mapreduce(
        x->x^2, 
        +,
        mapslices(
            x->au_to_ang.(x[2:end].-x[1]),#!1 
            rp_snapshot;
            dims=2
        ) #!2
        ;dims=1
    ))
    all_bead_trajectories(rp_trajectory)=cat(map(convert_to_nh_distance, [get_positions(frame) for frame in rp_trajectory])...;dims=1)#!3
    all_trajectories_for_datapoint(data)=map(all_bead_trajectories, [trj[:OutputDynamicsVariables] for trj in data[1]]) #!4
    return map(all_trajectories_for_datapoint, datapoint) # Returns a Matrix{Vector{Matrix{time, nh-distance, bead}}}
end

function NH_distance(datapoint)
    if isa(datapoint, String)
        datapoint=jldopen(datapoint, "r")["results"]
        println("Loading results. ")
    elseif isa(datapoint, Array)
        println("Results already loaded. ")
    end
    #=
    For each trajectory in a datapoint #!1, 
    for each frame in a trajectory #!2, 
    output the list of distances between N and the three H. I can't mash it all into a single matrix due to adaptive timestepping and dissociation #!3
    =#
    nh_distance(snapshot)=au_to_ang.(sqrt.(dropdims(sum(.^(snapshot[:,2:end].-snapshot[:,1], 2);dims=1);dims=1))) #!3
    nh_distance_trajectory(traj)=cat([nh_distance(get_positions(snapshot)) for snapshot in traj[:OutputDynamicsVariables]]...;dims=2) #!2
    nh_distance_trajectories(data)=[nh_distance_trajectory(traj) for traj in data[1]] #!1
    results=map(nh_distance_trajectories, datapoint)
    return results
end

function centroid_trajectories_to_xyz(datapoint, filename_base;center_on_atom=0)
    trj_atoms=datapoint[2]["atoms"]
    ase=pyimport("ase")
    ase_io=pyimport("ase.io")
    n_traj=length(datapoint[1])
    for i in 1:n_traj
        c_traj=datapoint[1][i][:ReplacementOutputCentroidPosition]
        if center_on_atom!=0
            modified_trj=map(x->x.-x[:,center_on_atom], c_traj)
        else
            modified_trj=c_traj
        end
        final_trj=NQCBase.convert_to_ase_atoms(
            trj_atoms, 
            modified_trj, 
            datapoint[2]["simulation"].cell
            )
        ase.io.write(
            filename_base*"_"*string(i;pad=length(digits(n_traj)))*".xyz", 
            final_trj
        )
    end
    return true
end

function trajectories_to_xyz(datapoint, filename_base;center_on_atom=0)
    trj_atoms=datapoint[2]["atoms"]
    ase=pyimport("ase")
    ase_io=pyimport("ase.io")
    n_traj=length(datapoint[1])
    for i in 1:n_traj
        c_traj=[get_positions(frame) for frame in datapoint[1][i][:OutputDynamicsVariables]]
        if center_on_atom!=0
            modified_trj=map(x->x.-x[:,center_on_atom], c_traj)
        else
            modified_trj=copy(c_traj)
        end
        final_trj=NQCBase.convert_to_ase_atoms(
            trj_atoms, 
            modified_trj, 
            datapoint[2]["simulation"].cell
            )
        ase.io.write(
            filename_base*"_"*string(i;pad=length(digits(n_traj)))*".xyz", 
            final_trj
        )
    end
    return true
end
