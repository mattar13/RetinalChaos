using CuArrays
using Dates
import RetinalChaos
import RetinalChaos: Network, extract_dict, params_to_datasheet, read_JSON
import RetinalChaos: tar_conds, tar_pars
import RetinalChaos: SDEProblem, noise_2D, SOSRI, solve
cd("C:\\Users\\RennaLabSA1\\Documents\\JuliaScripts\\RetinalChaos\\Parameters")
#Defining Initial conditions and parameters


dt = 10.0; tspan = (0.0, 300e3)
nx = 96; ny = 96; μ = 1.0
gpu = true
SACnet = Network(nx, ny; μ = μ, gpu = gpu, version = :gACh)
p_dict = read_JSON("params.json");
u_dict = read_JSON("conds.json");
u0 = extract_dict(u_dict, tar_conds, (nx, ny)) |> cu;
p0 = extract_dict(p_dict, tar_pars);

#warm up the model
println("[$(now())]: Warming up the model for 60s")
SDE_mat_prob = SDEProblem(SACnet, noise_2D, u0, (0.0, 60e3), p0);
CuArrays.allowscalar(false)
@time SDE_mat_sol = solve(
   SDE_mat_prob,
   SOSRI(),
   abstol = 0.2,
   reltol = 2e-2,
   maxiters = 1e7,
   progress = true,
   save_everystep = false,
);
#get the last solution from the warmup
println("[$(now())]: Running the model for $(tspan[end]/1000)s")
u0_new = SDE_mat_sol[end]
SDE_mat_prob = SDEProblem(SACnet, noise_2D, u0_new, tspan, p0);
@time SDE_mat_sol = solve(
   SDE_mat_prob,
   SOSRI(),
   abstol = 0.2,
   reltol = 2e-2,
   maxiters = 1e7,
   progress = true,
   saveat = dt,
);
println("[$(now())]: Model completed")
#Running repeated trials 1
p_rng = LinRange(7000, 10000, 25)
for x in p_rng
    println("Simulating tau_a = $x")
    cd("C:\\Users\\RennaLabSA1\\Documents\\JuliaScripts\\RetinalChaos\\Parameters")
    p_dict = read_JSON("params.json", is_type = Dict{Symbol,Float64});
    u_dict = read_JSON("conds.json", is_type = Dict{Symbol,Float64});
    p_dict[:τa] = x
    #If not on a GPU active computer, set gpu = false
    sol_array, df_params, df_stats = run_model(p_dict, u_dict, tspan)
    append_modeldata("DataSheet.xlsx", df_stats, df_params)
    #run the model
    mkdir("T7_10_tau_a_$x")
    cd("T7_10_tau_a_$x\\")
    tp = RetinalChaos.trace_plot(sol_array)
    savefig(tp, "Trace_plot.png")
    rp = RetinalChaos.raster_plot(sol_array)
    savefig(rp, "Raster.png")
    ap = RetinalChaos.frame_draw(sol_array)
end

#Trial 3
for x in p_rng
    println("Simulating tau_A = $x")
    cd("C:\\Users\\RennaLabSA1\\Documents\\JuliaScripts\\RetinalChaos\\Parameters")
    p_dict = read_JSON("params.json", is_type = Dict{Symbol,Float64});
    u_dict = read_JSON("conds.json", is_type = Dict{Symbol,Float64});
    p_dict[:k_d] = k_d
    sol_array, df_params, df_stats = run_model(p_dict, u_dict, tspan)
    append_modeldata("DataSheet.xlsx", df_stats, df_params)
    #run the model
    mkdir("T2_tau_A_$x")
    cd("T2_tau_A_$x\\")
    tp = RetinalChaos.trace_plot(sol_array)
    savefig(tp, "Trace_plot.png")
    rp = RetinalChaos.raster_plot(sol_array)
    savefig(rp, "Raster.png")
    ap = frame_draw(sol_array)
end

for x in p_rng
    println("Simulating tau_a = $x")
    cd("C:\\Users\\RennaLabSA1\\Documents\\JuliaScripts\\RetinalChaos\\Parameters")
    p_dict = read_JSON("params.json", is_type = Dict{Symbol,Float64});
    u_dict = read_JSON("conds.json", is_type = Dict{Symbol,Float64});
    p_dict[:k_d] = k_d
    sol_array, df_params, df_stats = run_model(p_dict, u_dict, tspan)
    append_modeldata("DataSheet.xlsx", df_stats, df_params)
    #run the model
    mkdir("T3_tau_A_$x")
    cd("T3_tau_A_$x\\")
    tp = RetinalChaos.trace_plot(sol_array)
    savefig(tp, "Trace_plot.png")
    rp = RetinalChaos.raster_plot(sol_array)
    savefig(rp, "Raster.png")
    ap = frame_draw(sol_array)
end
