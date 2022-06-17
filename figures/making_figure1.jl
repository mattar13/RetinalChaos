using Revise
using RetinalChaos

include("figure_setup.jl")

#%% Run a bare ODE model
#Step 1: Import the initial conditions
conds_dict = read_JSON("params\\conds.json")
u0 = conds_dict |> extract_dict
#Step 2: Import the parameters
pars_dict = read_JSON("params\\params.json")
pars_dict[:I_app] = 0.0
pars_dict[:ρi] = 0.0
pars_dict[:ρe] = 5.0
p = pars_dict |> extract_dict
#Step 3: determine the timespan
tspan = (0.0, 5e3)
#Step 4: set up the problem
prob = ODEProblem(GABA_ODE, u0, tspan, p)
#Step 5: Solve the problem
@time sol = solve(prob, progress=true, progress_steps = 1);
Plots.plot(sol, vars = [1])
#%%
# Run the analysis 
dt = 0.01
#v_thresh = calculate_threshold(sol; dt = dt)
#timestamps, data = timeseries_analysis(sol; dt = dt)
#timestamps["Bursts"]

t_series = tspan[1]:dt:tspan[end]
Vt = sol(t_series, idxs = [1])'


nothing #Really annoying how this keeps jumping to the bottom of the page

#=
# ╔═╡ 4e11ec08-cbb3-462d-8c5f-6eb37c5c79f9
begin 	
	burst_lims = timestamps["Bursts"][1][2,:]#lets set the limits for area of interest
	t_rng = burst_lims[1]:dt:burst_lims[2] #Set up the plotting range
	#Plot A tick intervals and limits
	A_dx = 20 #the tick interval is 20ms
	A_xlims = (burst_lims[1], burst_lims[1]+200)
	A_xticks = collect(A_xlims[1]:A_dx:A_xlims[2])
	A_trng = A_xlims[1]:dt:A_xlims[2]
	#Plot B tick intervals and limits
	B_dx = 500
	B_xlims = (burst_lims[1]-1500, burst_lims[2]+1500)
	B_xticks = collect(B_xlims[1]:B_dx:B_xlims[2])
	B_trng = B_xlims[1]:1.0:B_xlims[2]
	#Plot C tick intervals and limits
	C_dx = 10e3
	C_xlims = (timestamps["Bursts"][1][2,2], timestamps["Bursts"][1][4,2])
	C_xticks = collect((C_xlims[1]|>Float64):C_dx:(C_xlims[2]|>Float64))
	C_trng = C_xlims[1]:1.0:C_xlims[2]
end

# ╔═╡ abb10995-b168-4c2a-83a8-5cc386bca927
begin
	fig1_Aa1 = plot(A_trng, sol(A_trng, idxs = 1), label = "",
		ylabel = "Vₜ (mV)", xlabel = "", 
		lw = 4.0, c = v_color, grid = false,
		xlims = A_xlims, xticks = false, 
		yticks = (LinRange(-50, -20, 3)),
		margin = 0.0mm
	)
	
	fig1_Aa2 = plot(A_trng, sol(A_trng, idxs = 2), label = "",
		ylabel = "Nₜ", xlabel = "Time (ms)",
		lw = 4.0, c = n_color, grid = false,
		xlims = A_xlims, xticks = A_xticks, 
		x_formatter = x -> round(Int64,x-A_xlims[1]),
		yticks = (LinRange(0.0, 0.5, 3)),
		margin = 0.0mm	
	)
	
	fig1_Aa = plot(fig1_Aa1, fig1_Aa2, layout = grid(2,1))
		
	plot!(fig1_Aa[1], xticks = false)
	#hline!(fig1_Aa, [v_thresh], 
	#	label = "Spike threshold", c = :red, linewidth = 2.0, legend = :bottomright
	#)
	
	fig1_Ab = plot(sol(A_trng, idxs = 2), sol(A_trng, idxs = 1), 
		yticks = LinRange(-50, -20, 3), xticks = LinRange(0.0, 0.5, 3),
		xlabel = "Nₜ", ylabel = "Vₜ (mV)", lw = 3.0, label = "", c = :black
	)
	#hline!(fig1_Ab, [v_thresh], seriestype = :scatter,
	#	label = "Spike threshold",  c = :red, legend = :bottomright
	#)

	fig1_A = plot(fig1_Aa, fig1_Ab, 
		layout = grid(1,2, widths = (0.75, 0.25)), margin = 0mm
	)
	#title!(fig1_A[1], "A", titlepos = :left)
end;

# ╔═╡ 4e2a0bb6-4999-4519-931d-52c415205bac
begin
	fig1_Ba1 = plot(B_trng, sol(B_trng, idxs = 1), label = "",
		ylabel = "Vₜ (mV)", xlabel = "", 
		lw = 4.0, c = v_color, grid = false,
		xlims = B_xlims, xticks = false, 
		yticks = LinRange(-70.0, -20.0, 3),
		margin = 0.0mm
	)
	
	fig1_Ba2 = plot(B_trng, sol(B_trng, idxs = 3), label = "",
		ylabel = "[Cₜ] (μM)", xlabel = "Time (s)",
		lw = 4.0, c = c_color, grid = false,
		xlims = B_xlims, xticks = B_xticks, 
		yticks = LinRange(0.1, 0.35, 3), y_formatter = y -> round(y, digits = 2),
		x_formatter = x -> round((x-B_xlims[1])/1000, digits = 3),
		margin = 0.0mm	
	)
	
	fig1_Ba = plot(fig1_Ba1, fig1_Ba2, layout = grid(2,1))	
	
	f(v) = p[:δ]*(-p[:g_Ca] * M_INF(v, p[:V1], p[:V2]) * (v - p[:E_Ca]))
	fig1_Bb = plot(f, -50, 10.0, 
		xlabel = "Vₜ (mV)", ylabel = "[δC] (mM)",
		xticks = LinRange(-50.0, 10.0, 5),
		lw = 3.0, c = :green, label = "", linestyle = :dash, grid = false
	)
	fig1_B = plot(fig1_Ba, fig1_Bb, layout = grid(1,2, widths = (0.75, 0.25)))
	#title!(fig1_B[1], "B", titlepos = :left)
end;

# ╔═╡ b516f41c-4e30-47d2-aab7-acbfb561a7d4
begin	
	fig1_Ca1 = plot(C_trng, sol(C_trng, idxs = 3), label = "",
		ylabel = "[Cₜ] (μM)", xlabel = "", 
		lw = 4.0, c = c_color, grid = false,
		xlims = C_xlims, xticks = false, 
		yticks = LinRange(0.05, 0.35,3)
	)
	
	fig1_Ca2 = plot(C_trng, sol(C_trng, idxs = 4), label = "",
		ylabel = "Aₜ", xlabel = "",
		lw = 4.0, c = a_color, grid = false,
		xlims = C_xlims, xticks =false,
		yticks = LinRange(0.1, 0.6, 3)
	)

	fig1_Ca3 = plot(C_trng, sol(C_trng, idxs = 5), label = "",
		ylabel = "Bₜ", xlabel = "", 
		lw = 4.0, c = b_color, grid = false,
		xlims = C_xlims, xticks = false, 
		yticks = LinRange(0.10, 0.5, 3)
	)
	
	fig1_Ca4 = plot(C_trng, sol(C_trng, idxs = 1), label = "",
		ylabel = "Vₜ (mV)", xlabel = "Time (s)",
		lw = 4.0, c = v_color, grid = false,
		xlims = C_xlims, xticks = C_xticks, 
		x_formatter = x -> round(Int64, (x-C_xlims[1])/1000),
		yticks = LinRange(-70.0, -20.0, 3),
		margin = 0.0mm	
	)
	
	fig1_Ca = plot(fig1_Ca1, fig1_Ca2, fig1_Ca3, fig1_Ca4,layout = grid(4,1))
	fig1_Cb = plot(sol(C_trng, idxs = 4), sol(C_trng, idxs = 3), label = "", 
		ylabel = "[Cₜ] (μM)", xlabel = "Aₜ", 
		lw = 4.0, c = a_color, grid = false,
		yticks = LinRange(0.05, 0.35,3), xticks = LinRange(0.1, 0.6, 3)
	)
	
	fig1_Cc = plot(sol(C_trng, idxs = 5), sol(C_trng, idxs = 4), label = "", 
		ylabel = "Aₜ", xlabel = "Bₜ", 
		lw = 4.0, c = b_color, linewidth = 3.0, grid = false, 
		yticks = LinRange(0.1, 0.6, 3), xticks = LinRange(0.10, 0.5, 3)
	)
	
	fig1_Cd = plot(sol(C_trng, idxs = 1), sol(C_trng, idxs = 5), label = "", 
		ylabel = "Bₜ", xlabel = "Vₜ (mV)", 
		lw = 4.0, c = v_color, grid = false,
		yticks = LinRange(0.10, 0.5, 3), xticks = LinRange(-70.0, -20.0, 3)
	)
	
	fig1_C = plot(
		fig1_Ca, 
		plot(fig1_Cb, fig1_Cc, fig1_Cd, layout = grid(3,1)), 
		layout = grid(1,2, widths = (0.75, 0.25))
	)
	#title!(fig1_C[1], "C", titlepos = :left)
end;

# ╔═╡ 41af4f36-94a4-4d18-be6a-ca402e98801f
begin
	fig1_labels = plot(
		layout = grid(3,1, heights = (0.2, 0.3, 0.5)), 
		xaxis = false, yaxis = false, xticks = false, yticks = false
	)
	annotate!(fig1_labels[1], [0.5], [0.99], "A", font("Sans",24))
	annotate!(fig1_labels[2], [0.5], [0.99], "B", font("Sans",24))
	annotate!(fig1_labels[3], [0.5], [0.99], "C", font("Sans",24))
	
	
	
	fig1_boxes = plot(fig1_A, fig1_B, fig1_C, 
		layout = grid(3,1, heights = [0.2, 0.3, 0.5]), 
		 grid = false
		)
	fig1 = plot(
		fig1_labels, fig1_boxes, 
		layout = grid(1,2, widths = (0.05, 0.95)), size = (1000, 800),
	)
end

# ╔═╡ 327ff679-53b5-4b48-a579-327ea49312e1
savefig(fig1, "E:\\Projects\\2021_Modelling_Paper\\Figures\\Fig1_Model_Dynamics.png")
=#