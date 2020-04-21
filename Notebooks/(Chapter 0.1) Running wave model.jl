#Because the wave interface can get kinda slow in Julia, we have to run it in atom
using DifferentialEquations, Sundials
using Plots, Colors, LaTeXStrings
using ProgressLogging
using RetinalChaos
using DataFrames, XLSX
cd("C:\\Users\\mtarc\\JuliaScripts\\RetinalChaos\\Parameters")
#Defining Initial conditions and parameters
p_dict = read_JSON("params.json", is_type = Dict{Symbol,Float64});
u_dict = read_JSON("conds.json", is_type = Dict{Symbol,Float64});
p_dict[:g_ACh] = 1.0
tspan = (0.0, 60e3)
p_dict = read_JSON("params.json", is_type = Dict{Symbol,Float64});
u_dict = read_JSON("conds.json", is_type = Dict{Symbol,Float64});
p_dict[:k_d] = 0.001
sol_array, df_params, df_stats = run_model(p_dict, u_dict, tspan)
#Run 3 seperate trials
#Trial 1
k_d_rng = LinRange(0.001, 1.0, 10)
tspan = (0.0, 60e3)
for k_d in k_d_rng
    println("Simulating kd = $k_d")
    cd("C:\\Users\\mtarc\\JuliaScripts\\RetinalChaos\\Parameters")
    p_dict = read_JSON("params.json", is_type = Dict{Symbol,Float64});
    u_dict = read_JSON("conds.json", is_type = Dict{Symbol,Float64});
    p_dict[:k_d] = k_d
    sol_array, df_params, df_stats = run_model(p_dict, u_dict, tspan)
    append_modeldata("DataSheet.xlsx", df_stats, df_params)
    #run the model
    mkdir("kd_$k_d")
    cd("kd_$k_d\\")
    tp = RetinalChaos.trace_plot(sol_array)
    savefig(tp, "Trace_plot.png")
    rp = RetinalChaos.raster_plot(sol_array)
    savefig(rp, "Raster.png")
    ap = RetinalChaos.frame_draw(sol_array)
end
