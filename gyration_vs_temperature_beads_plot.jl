include("analysis_functions.jl")
using LaTeXStrings
using Plots

names=["N", "H"]
beads=[4,8,16,32]
temperatures=collect(100:100:1000).*u"K"
colors=["royalblue", "rebeccapurple", "forestgreen", "firebrick"]

results=jldopen("../rpsh_surfgen_4beads_temperature/results_surfgen_rpsh_beads_temperature.jld2")["results"]
results_flat=sort(reshape(results, :);by=x->x[2]["temperature"])
results_fixed=reshape(results_flat, 4,10)

r_gyr=calculate_gyration_radius(results_fixed)
r_gyr_N=map(x->x[1], r_gyr)
r_gyr_H=map(x->mean(x[2:end]), r_gyr)
outputs=[r_gyr_N, r_gyr_H]

for i in 1:2
    plot(;
        xlabel="Temperature / K", 
        ylabel=L"R_{gyr}"*" / Å", 
        minorticks=true, 
        xticks=collect(0:100:1000).*u"K",
        size=(600,400),
        left_margin=10Plots.mm,
        bottom_margin=10Plots.mm,
        );
    for j in eachindex(beads)
        plot!(temperatures, outputs[i][j,:]; label=string(beads[j])*" beads", color=colors[j], ms=5, marker=:x);
    end
    plot!(;legend=:topright);
    savefig("/storage/mssgwp_grp/msrkhg/ammonia_rpsh/rpsh_surfgen_4beads_temperature/figure_r_gyr_"*names[i]*".pdf")
    savefig("/storage/mssgwp_grp/msrkhg/ammonia_rpsh/rpsh_surfgen_4beads_temperature/figure_r_gyr_"*names[i]*".svg")
end


