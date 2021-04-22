### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ a2755664-125e-4b21-9784-1b2e50e45f50
using Revise

# ╔═╡ 0bfc7db2-9c77-11eb-3262-f1f062593489
using RetinalChaos

# ╔═╡ bc8ecd2f-2d2d-40a5-99c2-bb5377689c5f
#We need to import my other package NeuroPhys
using NeuroPhys

# ╔═╡ 4886869d-19d7-43c9-90c0-2c8a1e5a18f1
using Dates,JLD2

# ╔═╡ 926d1c0c-171b-4001-991c-1c822b3e2fb8
begin
	#Set up fonts 
	font_title = font("Arial", 24) #title font and size
	font_axis = font("Arial", 12) #axis font and size
	font_legend = font("Arial", 8) #legend font and size
	gr(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)
end

# ╔═╡ bd8eddc8-d266-4a63-9151-4782a1e72be8
begin
	#Set up the parameter files and save files
	param_root = joinpath(splitpath(pwd())[1:end-1]..., "params\\")
	params_file = joinpath(param_root, "params.json")
	conds_file = joinpath(param_root, "conds.json")
end

# ╔═╡ 0fd4fa03-d96c-42b4-b78e-747e716588d6
md"
#### Run a simulation for many single traces

"

# ╔═╡ 0bc9c0b4-74bd-444e-bc1b-27734e760d7d
begin
	#Set up the parameters 
	p = read_JSON(params_file) 
	p[:σ] = 0.1
	p[:τw] = 800.0
	#Set up the initial conditions
	u0 = read_JSON(conds_file);
	#Set up the Time span
	tspan = (0.0, 100e3)
	#set up the problem to solve
	prob = SDEProblem(T_sde, u0 |> extract_dict, tspan, p |> extract_dict);
	#Load the problem into a ensemble simulation
	n_sims = 50
	test_rng = repeat([p[:σ]], n_sims)
	par_idx = findall(x -> x==:σ, Symbol.(T_sde.ps))
	prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par_idx, test_rng)
	#make the ensemble problem
	simProb = EnsembleProblem(prob, prob_func = prob_func);
end

# ╔═╡ b93d9c06-126f-4e67-adb4-27027fd6bad4
@time sim = solve(simProb, trajectories = n_sims, EnsembleThreads());

# ╔═╡ a82156f5-ff66-4e64-91d1-8d126b0acf46
begin
	iso_spike_durs = Float64[]
	iso_burst_durs = Float64[]
	iso_ibis = Float64[]
	for (sim_idx, sol) in enumerate(sim)
		println(sim_idx)
		#Run the timescale analysis now for real
		dt = 0.1 #set the time differential
		v_thresh = RetinalChaos.calculate_threshold(sol, dt = dt)
		println(v_thresh)
		timestamps = RetinalChaos.get_timestamps(sol, dt = dt)
		sd, bd, ibi = RetinalChaos.timescale_analysis(sol, dt = dt)
		push!(iso_spike_durs, sd...)
		push!(iso_burst_durs, bd...)
		push!(iso_ibis, ibi...)
	end
end

# ╔═╡ bd50617a-ad3a-40f1-9d80-3fea8facf83e
#Load the solution\
JLD2.@load "sol.jld2" NetSol

# ╔═╡ 0379e9a9-d8e7-48f8-aaf6-5e04461e08eb
JLD2.@load "thresholds.jld2" thresholds

# ╔═╡ 60197d1a-fff7-4c55-bdbd-edcf7f534710
JLD2.@load "C:\\Users\\mtarc\\OneDrive\\Documents\\GithubRepositories\\RetinalChaos\\figures\\ts_analysis.jld2" ts

# ╔═╡ 280a8f63-7746-407a-9a28-c5e83e15273c
begin
	plt_sd = plot()
	
	#jitter1 = randn(length(iso_spike_durs))/10.0 .+ 1.0
	#plot!(plt_sd, jitter1, iso_spike_durs, 
	#	st = :scatter, ms = 0.4
	#)
	plot!(plt_sd, iso_spike_durs, st = :violin, 
		xticks = ([1,2], ["Isolated", "Network"]), 
		legend = false
	)
	plot!(plt_sd, ts[1], st = :violin)
	
	
	plt_bd = plot(iso_burst_durs, st = :violin, 
		xticks = ([1,2], ["Isolated", "Network"]),
		legend = false, 
		yformatter = x -> x/1000
	)
	plot!(plt_bd, ts[2], st = :violin)
	
	plt_ibi = plot(iso_ibis, st = :violin, 
		xticks = ([1,2], ["Isolated", "Network"]),
		legend = false, 
		yformatter = x -> x/1000
	)
	plot!(plt_ibi, ts[3], st = :violin)
	
	plot(plt_sd, plt_bd, plt_ibi, layout = grid(1,3))
end

# ╔═╡ a422bb0c-f779-44f6-9141-9d8694dd3239
begin
	#Extract the physiological data file
	target_file = 				"E:\\Data\\Patching\\2019_11_03_Patch\\Animal_2\\Cell_3\\19n03042.abf"
	data = extract_abf(target_file)
	#Conduct the timescale analysis on this file
	#reduce the size of the file
	truncate_data!(data, 
		t_pre = 140.0, t_post = 240.0, 
		truncate_based_on = :time_range)
	data.data_array .-= 25
end

# ╔═╡ 4e6c5aee-c93d-4752-ba70-c97d828381d0
begin
	plot(NetSol(NetSol.t, idxs = 1), xformatter = x -> x/1000)
	plot!(data.t[1:10:length(data.t)].*1000, data.data_array[1:10:length(data.t)])
end

# ╔═╡ abc98139-8c4b-488e-ae88-961b6edbb5a6


# ╔═╡ 164b2668-f893-46eb-8bb7-025391163840
md"
# Analyze multiple files

"

# ╔═╡ a3266346-58f0-4ef5-b747-c88567da7a0d
begin
	#Open multiple files
	target_folder = "E:\\Data\\Jordans_Patch_Data\\UsuableData\\"
	paths = target_folder |> parse_abf;	
end

# ╔═╡ bcba180e-fc70-4c55-8e3f-59979216e146
begin
	#Walk through each file analyzing using timescale analysis
end

# ╔═╡ Cell order:
# ╠═a2755664-125e-4b21-9784-1b2e50e45f50
# ╠═0bfc7db2-9c77-11eb-3262-f1f062593489
# ╠═bc8ecd2f-2d2d-40a5-99c2-bb5377689c5f
# ╠═4886869d-19d7-43c9-90c0-2c8a1e5a18f1
# ╠═926d1c0c-171b-4001-991c-1c822b3e2fb8
# ╠═bd8eddc8-d266-4a63-9151-4782a1e72be8
# ╟─0fd4fa03-d96c-42b4-b78e-747e716588d6
# ╠═0bc9c0b4-74bd-444e-bc1b-27734e760d7d
# ╠═b93d9c06-126f-4e67-adb4-27027fd6bad4
# ╠═a82156f5-ff66-4e64-91d1-8d126b0acf46
# ╠═bd50617a-ad3a-40f1-9d80-3fea8facf83e
# ╠═0379e9a9-d8e7-48f8-aaf6-5e04461e08eb
# ╠═60197d1a-fff7-4c55-bdbd-edcf7f534710
# ╟─280a8f63-7746-407a-9a28-c5e83e15273c
# ╠═a422bb0c-f779-44f6-9141-9d8694dd3239
# ╠═4e6c5aee-c93d-4752-ba70-c97d828381d0
# ╠═abc98139-8c4b-488e-ae88-961b6edbb5a6
# ╟─164b2668-f893-46eb-8bb7-025391163840
# ╠═a3266346-58f0-4ef5-b747-c88567da7a0d
# ╠═bcba180e-fc70-4c55-8e3f-59979216e146
