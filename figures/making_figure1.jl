### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ 1f3a5020-9a18-11eb-3b98-2fa7c5574805
using RetinalChaos

# ╔═╡ 39edbae4-3774-4de8-aa01-7597be0cc8d7
import RetinalChaos.M_INF

# ╔═╡ 80ce4218-d0b0-4fa6-9784-74fb9b68bc4e
begin
	#Set up fonts 
	font_title = font("Arial", 24) #title font and size
	font_axis = font("Arial", 12) #axis font and size
	font_legend = font("Arial", 8) #legend font and size
	gr(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)
end

# ╔═╡ e908115c-8389-49a1-8d0a-530d292f70c0
begin
	#Set up the parameter files and save files
	param_root = joinpath(splitpath(pwd())[1:end-1]..., "params\\")
	params_file = joinpath(param_root, "params.json")
	conds_file = joinpath(param_root, "conds.json")

	#save everything in the figures folder
	save_figs = "figures\\"
	if isdir(save_figs) == false
		#The directory does not exist, we have to make it 
		mkdir(save_figs)
	end	
end

# ╔═╡ 1cb6ef5d-ddb0-4c90-9bc0-0c3924399fbe
begin
	#Set up the parameters 
	p = read_JSON(params_file) 
	p[:I_app] = 10.0
	p[:g_ACh] = 0.0
	#Set up the initial conditions
	u0 = read_JSON(conds_file);
	#Set up the Time span
	tspan = (0.0, 120e3)
	#set up the problem to solve
	prob = ODEProblem(T_ode, u0 |> extract_dict, tspan, p |> extract_dict);
end

# ╔═╡ b215888a-6acf-4b9e-aa7e-4dec87b8738b
#Solve the problem
@time sol = solve(prob, progress = true); 

# ╔═╡ 578be207-80cb-4738-aa39-3db181b7c058
begin 
	#Run analysis 
	dt = 0.1 #set the time differential
	v_thresh = calculate_threshold(sol, dt = dt)
	timestamps = get_timestamps(sol, dt = dt)
	burst_idxs, dur_list, spb_list, ibi_list = max_interval_algorithim(sol, dt = dt)
end

# ╔═╡ 4e11ec08-cbb3-462d-8c5f-6eb37c5c79f9
begin 	
	burst_lims = burst_idxs[2]	#lets set the limits for area of interest
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
	C_xlims = (burst_idxs[2][1], burst_idxs[4][1])
	C_xticks = collect(C_xlims[1]:C_dx:C_xlims[2]) 
	C_trng = C_xlims[1]:1.0:C_xlims[2]
end

# ╔═╡ abb10995-b168-4c2a-83a8-5cc386bca927
begin
	fig1_Aa1 = plot(A_trng, sol(A_trng, idxs = 1), label = "",
		ylabel = "Vₜ (mV)", xlabel = "", 
		lw = 2.0, c = v_color, grid = false,
		xlims = A_xlims, xticks = false, 
		margin = 0.0mm
	)
	
	fig1_Aa2 = plot(A_trng, sol(A_trng, idxs = 2), label = "",
		ylabel = "Nₜ", xlabel = "Time (ms)",
		lw = 2.0, c = n_color, grid = false,
		xlims = A_xlims, xticks = A_xticks, 
		x_formatter = x -> round(Int64,x-A_xlims[1]),
		margin = 0.0mm	
	)
	
	fig1_Aa = plot(fig1_Aa1, fig1_Aa2, layout = grid(2,1))
		
	plot!(fig1_Aa[1], xticks = false)
	hline!(fig1_Aa, [v_thresh], 
		label = "Spike threshold", c = :red, linewidth = 2.0, legend = :bottomright
	)
	
	fig1_Ab = plot(sol(A_trng, idxs = 2), sol(A_trng, idxs = 1), 
		xlabel = "Nₜ", ylabel = "Vₜ (mV)", label = "", c = :black
	)
	hline!(fig1_Ab, [v_thresh], seriestype = :scatter,
		label = "Spike threshold",  c = :red, legend = :bottomright
	)

	fig1_A = plot(fig1_Aa, fig1_Ab, 
		layout = grid(1,2, widths = (0.75, 0.25)), margin = 0mm
	)
	title!(fig1_A[1], "A", titlepos = :left)
end

# ╔═╡ 4e2a0bb6-4999-4519-931d-52c415205bac
begin
	fig1_Ba1 = plot(B_trng, sol(B_trng, idxs = 1), label = "",
		ylabel = "Vₜ (mV)", xlabel = "", 
		lw = 2.0, c = v_color, grid = false,
		xlims = B_xlims, xticks = false, 
		margin = 0.0mm
	)
	
	fig1_Ba2 = plot(B_trng, sol(B_trng, idxs = 3), label = "",
		ylabel = "[Cₜ] (μM)", xlabel = "Time (s)",
		lw = 2.0, c = c_color, grid = false,
		xlims = B_xlims, xticks = B_xticks, 
		x_formatter = x -> round((x-B_xlims[1])/1000, digits = 3),
		margin = 0.0mm	
	)
	
	fig1_Ba = plot(fig1_Ba1, fig1_Ba2, layout = grid(2,1))	
	
	f(v) = p[:δ]*(-p[:g_Ca] * M_INF(v, p[:V1], p[:V2]) * (v - p[:E_Ca]))
	fig1_Bb = plot(f, -50, 10.0, 
		xlabel = "Vₜ (mV)", ylabel = "[δC] (mM)",
		lw = 3.0, c = :green, label = "", linestyle = :dash, grid = false
	)
	fig1_B = plot(fig1_Ba, fig1_Bb, layout = grid(1,2, widths = (0.75, 0.25)))
	title!(fig1_B[1], "B", titlepos = :left)
end

# ╔═╡ b516f41c-4e30-47d2-aab7-acbfb561a7d4
begin	
	fig1_Ca1 = plot(C_trng, sol(C_trng, idxs = 3), label = "",
		ylabel = "[Cₜ] (μM)", xlabel = "", 
		lw = 2.0, c = c_color, grid = false,
		xlims = C_xlims, xticks = false, 
	)
	
	fig1_Ca2 = plot(C_trng, sol(C_trng, idxs = 4), label = "",
		ylabel = "Aₜ", xlabel = "",
		lw = 2.0, c = a_color, grid = false,
		xlims = C_xlims, xticks =false,
	)

	fig1_Ca3 = plot(C_trng, sol(C_trng, idxs = 5), label = "",
		ylabel = "Bₜ", xlabel = "", 
		lw = 2.0, c = b_color, grid = false,
		xlims = C_xlims, xticks = false, 
	)
	
	fig1_Ca4 = plot(C_trng, sol(C_trng, idxs = 1), label = "",
		ylabel = "Vₜ (mV)", xlabel = "Time (s)",
		lw = 2.0, c = v_color, grid = false,
		xlims = C_xlims, xticks = C_xticks, 
		x_formatter = x -> round(Int64, (x-C_xlims[1])/1000),
		margin = 0.0mm	
	)
	
	fig1_Ca = plot(fig1_Ca1, fig1_Ca2, fig1_Ca3, fig1_Ca4,layout = grid(4,1))
	fig1_Cb = plot(sol(C_trng, idxs = 4), sol(C_trng, idxs = 3), label = "", 
		ylabel = "[Cₜ] (μM)", xlabel = "Aₜ", 
		lw = 3.0, c = a_color, grid = false
	)
	
	fig1_Cc = plot(sol(C_trng, idxs = 5), sol(C_trng, idxs = 4), label = "", 
		ylabel = "Aₜ", xlabel = "Bₜ", 
		lw = 3.0, c = b_color, linewidth = 3.0, grid = false
	)
	
	fig1_Cd = plot(sol(C_trng, idxs = 1), sol(C_trng, idxs = 5), label = "", 
		ylabel = "Bₜ", xlabel = "Vₜ (mV)", 
		lw = 3.0, c = v_color, grid = false
	)
	
	fig1_C = plot(
		fig1_Ca, 
		plot(fig1_Cb, fig1_Cc, fig1_Cd, layout = grid(3,1)), 
		layout = grid(1,2, widths = (0.75, 0.25))
	)
	title!(fig1_C[1], "C", titlepos = :left)
end

# ╔═╡ 41af4f36-94a4-4d18-be6a-ca402e98801f
fig1 = plot(fig1_A, fig1_B, fig1_C, 
    layout = grid(3,1, heights = [0.2, 0.3, 0.5]), 
    size = (1000, 1000), grid = false
    )

# ╔═╡ 327ff679-53b5-4b48-a579-327ea49312e1
savefig(fig1, "Fig1_Model_Dynamics.png")

# ╔═╡ Cell order:
# ╠═1f3a5020-9a18-11eb-3b98-2fa7c5574805
# ╠═39edbae4-3774-4de8-aa01-7597be0cc8d7
# ╠═80ce4218-d0b0-4fa6-9784-74fb9b68bc4e
# ╠═e908115c-8389-49a1-8d0a-530d292f70c0
# ╠═1cb6ef5d-ddb0-4c90-9bc0-0c3924399fbe
# ╠═b215888a-6acf-4b9e-aa7e-4dec87b8738b
# ╠═578be207-80cb-4738-aa39-3db181b7c058
# ╠═4e11ec08-cbb3-462d-8c5f-6eb37c5c79f9
# ╠═abb10995-b168-4c2a-83a8-5cc386bca927
# ╠═4e2a0bb6-4999-4519-931d-52c415205bac
# ╠═b516f41c-4e30-47d2-aab7-acbfb561a7d4
# ╠═41af4f36-94a4-4d18-be6a-ca402e98801f
# ╠═327ff679-53b5-4b48-a579-327ea49312e1
