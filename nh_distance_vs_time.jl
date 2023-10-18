include("analysis_functions.jl")
using LaTeXStrings
using Plots
using Unitful

base_path="/storage/mssgwp_grp/msrkhg/ammonia_rpsh/rpsh_surfgen_4beads_temperature/centroid_nh_distances/"
beads=[4,8,16,32]
temperature=collect(100:100:1000).*u"K"
time=collect(0:0.01:5.01).*u"ps"
colors=["royalblue", "rebeccapurple", "forestgreen", "firebrick"]

results=NH_distance("/storage/mssgwp_grp/msrkhg/ammonia_rpsh/rpsh_surfgen_4beads_temperature/results_surfgen_rpsh_beads_temperature.jld2") # Import all data
for i in eachindex(beads)
    for j in eachindex(temperature)
        plot([
            plot(
            time, 
            permutedims(results[i,j][:,:,k], [2,1]);
            label=cat(["H1" "H2" "H3"]...;dims=2),
            color=cat(colors[1:3]...;dims=2),
            legend=k==1 ? true : false,
            xlabel="Time / ps",
            ylabel=L"R_{H}-R_{N}"*" / Å", 
            ylims=(1.6,2.4),
            minorticks=true, 
            xticks=collect(0.0:0.5:5).*u"ps",
            size=(600,400*10),
            left_margin=10Plots.mm,
            bottom_margin=10Plots.mm,
        );
        for k in 1:size(results[1,1])[end]
        ]...;layout=(:,1));
        
        savefig(base_path*"T_"*string(ustrip(temperature[j]))*"K_"*string(beads[i])*"beads.pdf");
        savefig(base_path*"T_"*string(ustrip(temperature[j]))*"K_"*string(beads[i])*"beads.svg");
    end
end


