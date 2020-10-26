#Because the wave interface can get kinda slow in Julia, we have to run it in atom
using ProgressLogging
using DataFrames, XLSX

using RetinalChaos
import RetinalChaos: read_JSON, Network, noise_2D, extract_dict
import RetinalChaos: SDEProblem, SOSRI, solve
import RetinalChaos: savefig, trace_plot, raster_plot

cd("C:\\Users\\mtarc\\JuliaScripts\\RetinalChaos\\Parameters")
#Defining Initial conditions and parameters
tspan = (0.0, 60e3)
p_dict = read_JSON("params.json", is_type = Dict{Symbol,Float64});
u_dict = read_JSON("conds.json", is_type = Dict{Symbol,Float64});
nx = ny = 96
SACnet = Network(nx, ny; μ = 0.60, nullout = :gACh)
u0_mat = cat(map(x -> fill(u_dict[x], (ny, nx)), model_conds)..., dims = 3)
u0 = extract_dict(u_dict, model_conds)
p0 = extract_dict(p_dict, model_pars)
#warm up the model
println("Warming up the model for 60s")
SDE_mat_prob = SDEProblem(SACnet, noise_2D, u0_mat, (0.0, 60e3), p0);
@time SDE_mat_sol = solve(
    SDE_mat_prob,
    SOSRI(),
    abstol = 0.2,
    reltol = 2e-2,
    maxiters = 1e7,
    progress = true,
    save_everystep = false,
    #saveat = 100
);
#get the last solution from the warmup
println("Running the model for $(tspan[end]/1000)s")
u0_new = SDE_mat_sol[end]
SDE_mat_prob = SDEProblem(SACnet, noise_2D, u0_new, tspan, p0);
@time SDE_mat_sol = solve(
    SDE_mat_prob,
    SOSRI(),
    abstol = 0.2,
    reltol = 2e-2,
    maxiters = 1e7,
    progress = true,
    saveat = 10.0,
);
sol_array = Array(SDE_mat_sol)
#run the model
mkdir("kd_$k_d")
cd("kd_$k_d\\")
tp = RetinalChaos.trace_plot(sol_array)
savefig(tp, "Trace_plot.png")
rp = RetinalChaos.raster_plot(sol_array)
savefig(rp, "Raster.png")
ap = RetinalChaos.frame_draw(sol_array)



#playing around
