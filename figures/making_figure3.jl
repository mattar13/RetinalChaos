### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ 7db66360-9c6f-11eb-0c9a-27ad7a4cc638
using RetinalChaos

# ╔═╡ bc5c19fa-a5a6-40f9-b081-10b8a7d2ccdc
using StatsPlots

# ╔═╡ b8a24735-3ad0-45c5-8584-6fe27dc67162
md"
# Figure 3: Spontaneous Activity
"
	
# ╔═╡ 85b88d7a-7e08-4b31-9438-647073d00ae2
begin
	#Set up fonts 
	font_title = font("Arial", 24) #title font and size
	font_axis = font("Arial", 12) #axis font and size
	font_legend = font("Arial", 8) #legend font and size
	gr(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)
end

# ╔═╡ 5d84dba7-a262-43f9-9364-5e13d322d832
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

# ╔═╡ 43050c9b-76c1-47df-88c4-3b1ff215497b
begin
	#Set up the parameters
	p = read_JSON(params_file)
	p[:σ] = 0.1 #Adjust noise amplitude paramaters
	p[:τw] = 1000.0 #Adjust the drift parameters
	u0 = read_JSON(conds_file);
	#Set up the Time span
	tspan = (0.0, 120e3)
	#set up the problem to solve
	prob = SDEProblem(T_sde, u0 |> extract_dict, tspan, p |> extract_dict);

end

# ╔═╡ ea8baa30-f4b4-4f56-87ba-5b03aec8c574
#Solve the problem using SOSRI algorithim
@time sol = solve(prob, 
		SOSRI(), abstol = 2e-2, reltol = 2e-2, maxiters = 1e7, 
		progress = true
	); 

# ╔═╡ f57aa1f7-41b4-435f-b6fc-24815a655b93
begin 
	#Run analysis 
	dt = 0.1 #set the time differential
	v_thresh = calculate_threshold(sol, dt = dt)
	timestamps = get_timestamps(sol, dt = dt)
	burst_idxs, dur_list, spb_list, ibi_list = max_interval_algorithim(sol, dt = dt)
end

# ╔═╡ 27f07f0a-bff3-4ab1-8997-8d56ef4f84fc
begin
	fig3_A = plot(sol, vars = [:v, :W], layout = grid(2,1),label = ["" ""],
		c = [v_color w_color],lw = 3.0, grid = false,
		xlabel = ["" "Time (s)"], 
		ylabel = ["Vₜ (mV)" "Noise (pA)"], 
		xformatter = x -> x/1000, 
	
	)
	title!(fig3_A[1], "A", title_pos = :left)
end

# ╔═╡ bdc1f75d-d992-4b03-a6ad-b61650d8db94
begin
	#Set up bifurcation limits
	c1_lims = (minimum(sol(sol.t, idxs=7)), maximum(sol(sol.t, idxs=7)))
	#Set up a ideal condition for the bifurcation analysis 
	pEQ = read_JSON(params_file) 
	pEQ[:I_app] = 0.0 #Set initial applied current to 0
	pEQ[:g_ACh] = 0.0 #Remove g_ACh influence
	pEQ[:g_TREK] = 0.0 #Remove the sAHP
	#Set up the equilibrium problem
	prob_eq = ODEProblem(T_ode, u0|>extract_dict, (0.0, 30e3), pEQ|>extract_dict)
	codim1 = (:I_app)
	@time c1_map = codim_map(prob_eq, codim1, c1_lims, equilibrium_resolution = 10)
	bif_val, bif_eq = find_bifurcation(c1_map)
	saddle_vs = map(x -> x.saddle[1][1], bif_eq)
end

# ╔═╡ 79fd10a2-b662-469e-ab24-508b78ec53b4
begin
	fig3_B = vspan([bif_val[end], c1_lims[end]], label = "Spiking Range",
			c = :green, alpha = 0.5
		) 
	vspan!(fig3_B, [c1_lims[1], bif_val[end]], label = "Non-spiking range", 
		c = :red, alpha = 0.5
	)
	plot!(fig3_B, c1_map, 
		lw = 3.0,
		xlabel = "Injected Current", ylabel = "Membrane Voltage"
	)

	plot!(fig3_B, bif_val, saddle_vs, label= "Saddle node bif",
		marker = :square, st = :scatter, grid = false, 
		xlims = c1_lims
	)
	
end

# ╔═╡ df7e26fa-a550-4594-bc5b-f03c38971381
begin
	#Draw a histogram of noise
	fig3_C = vspan([bif_val[end], c1_lims[end]], label = "Spiking Range",
			c = :green, alpha = 0.5, grid = false
		) 
	vspan!(fig3_C, [c1_lims[1], bif_val[end]], label = "Non-spiking range", 
		c = :red, alpha = 0.5
	)
	histogram!(fig3_C, sol(sol.t, idxs = 7), label = "",
		c = :blue,
		normalize = :pdf,
		xlabel = "Noisy current (pA)", ylabel = "Probability of Occurrance (PDF)", 
		xlims = c1_lims
	)
	
end

# ╔═╡ 976de364-ca10-4c91-8d9d-143926e48613
fig3 = plot(fig3_A, fig3_B, fig3_C, layout = grid(3,1), size  = (1000,1000))

# ╔═╡ b5d23852-523d-4c6f-b374-316a4052c87d
savefig(fig3, "Fig3_Spontneous_Activity.png")

# ╔═╡ Cell order:
# ╠═7db66360-9c6f-11eb-0c9a-27ad7a4cc638
# ╠═bc5c19fa-a5a6-40f9-b081-10b8a7d2ccdc
# ╟─b8a24735-3ad0-45c5-8584-6fe27dc67162
# ╠═85b88d7a-7e08-4b31-9438-647073d00ae2
# ╟─5d84dba7-a262-43f9-9364-5e13d322d832
# ╠═43050c9b-76c1-47df-88c4-3b1ff215497b
# ╠═ea8baa30-f4b4-4f56-87ba-5b03aec8c574
# ╠═f57aa1f7-41b4-435f-b6fc-24815a655b93
# ╟─27f07f0a-bff3-4ab1-8997-8d56ef4f84fc
# ╟─bdc1f75d-d992-4b03-a6ad-b61650d8db94
# ╟─79fd10a2-b662-469e-ab24-508b78ec53b4
# ╟─df7e26fa-a550-4594-bc5b-f03c38971381
# ╟─976de364-ca10-4c91-8d9d-143926e48613
# ╠═b5d23852-523d-4c6f-b374-316a4052c87d
