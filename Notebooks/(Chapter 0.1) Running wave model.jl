#Because the wave interface can get kinda slow in Julia, we have to run it in atom
using CuArrays
using RetinalChaos
cd("C:\\Users\\RennaLabSA1\\Documents\\JuliaScripts\\RetinalChaos\\Parameters")
#Defining Initial conditions and parameters
p_dict = read_JSON("params.json", is_type = Dict{Symbol,Float64});
u_dict = read_JSON("conds.json", is_type = Dict{Symbol,Float64});
p_dict[:g_ACh] = 1.0
tspan = (0.0, 60e3)

#Trial 1
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
