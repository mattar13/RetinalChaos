### A Pluto.jl notebook ###
# v0.14.2

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
using Dates,JLD2, Statistics, StatsPlots

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
	prob = SDEProblem(T_sde, noise, u0 |> extract_dict, tspan, p |> extract_dict);
	#Load the problem into a ensemble simulation
	n_sims = 50
	test_rng = repeat([p[:σ]], n_sims)
	prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, :σ, test_rng)
	#make the ensemble problem
	simProb = EnsembleProblem(prob, prob_func = prob_func);
end

# ╔═╡ b93d9c06-126f-4e67-adb4-27027fd6bad4
@time sim = solve(simProb, trajectories = n_sims, EnsembleThreads());

# ╔═╡ a82156f5-ff66-4e64-91d1-8d126b0acf46
begin
	iso_baseline = Float64[]
	iso_max = Float64[]
	iso_min = Float64[]
	iso_spike_durs = Float64[]
	iso_burst_durs = Float64[]
	iso_ibis = Float64[]
	for (sim_idx, sol) in enumerate(sim)
		println(sim_idx)
		#Run the timescale analysis now for real
		dt = 0.1 #set the time differential
		v_thresh = RetinalChaos.calculate_threshold(sol, dt = dt)
		v_base = RetinalChaos.calculate_threshold(sol, Z = 0, dt = dt)
		println(v_thresh)
		timestamps = RetinalChaos.get_timestamps(sol, dt = dt)
		sd, bd, ibi = RetinalChaos.timescale_analysis(sol, dt = dt)
		push!(iso_baseline, v_base[1])
		push!(iso_max, maximum(sol(sol.t, idxs = 1)))
		push!(iso_min, minimum(sol(sol.t, idxs = 1)))
		push!(iso_spike_durs, sd...)
		push!(iso_burst_durs, bd...)
		push!(iso_ibis, ibi...)
	end
end

# ╔═╡ 5fcd8df1-9889-41f8-8f3e-c53b00edcf75
md"
## Isolated
- Baseline $(sum(iso_baseline)/length(iso_baseline)) +- $(std(iso_baseline)/sqrt(length(length(iso_baseline))))

- Max Amp $(sum(iso_max)/length(iso_max)) +- $(std(iso_max)/sqrt(length(length(iso_max))))

- Min Amp $(sum(iso_min)/length(iso_min)) +- $(std(iso_min)/sqrt(length(length(iso_min))))

- Spike Duration = $(sum(iso_spike_durs)/length(iso_spike_durs)) +- $(std(iso_spike_durs)/sqrt(length(length(iso_spike_durs))))

- Burst Duration = $((sum(iso_burst_durs)/length(iso_burst_durs))/1000) +- $((std(iso_burst_durs)/sqrt(length(length(iso_burst_durs))))/1000)

- Interburst Interval = $((sum(iso_ibis)/length(iso_ibis))/1000) +- $((std(iso_ibis)/sqrt(length(length(iso_ibis))))/1000)

"

# ╔═╡ bd50617a-ad3a-40f1-9d80-3fea8facf83e
#Load the solution\
JLD2.@load "sol.jld2" NetSol

# ╔═╡ 0379e9a9-d8e7-48f8-aaf6-5e04461e08eb
#JLD2.@load "thresholds.jld2" thresholds

# ╔═╡ 1450b118-bba0-49eb-b4ca-630ee65074fa
n_points = size(NetSol, 1) * size(NetSol,2)

# ╔═╡ b96e7e0b-945a-4b3b-8549-c2412c0be3a6
begin
	net_baselines = Float64[]
	net_mins = Float64[]
	net_maxs = Float64[]
	for sim_idx = 1:size(NetSol,1)
		println(sim_idx)
		vt_sim = NetSol(NetSol.t,idxs = sim_idx)
		push!(net_baselines, sum(vt_sim)/length(vt_sim))
		push!(net_mins, minimum(vt_sim))
		push!(net_maxs, maximum(vt_sim))
	end	
end

# ╔═╡ b69cdf52-3dc2-4cbc-989e-4ea838dd70cd
#Print out the timescale characteristics

# ╔═╡ ed3f2c93-34a8-454b-9e5a-193ce88e0fe3
#Calculate the baseline
#RetinalChaos.calculate_threshold(NetSol; Z = 0)

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
ts = NeuroPhys.timescale_analysis(data)

# ╔═╡ 60197d1a-fff7-4c55-bdbd-edcf7f534710
JLD2.@load "ts_analysis.jld2" ts

# ╔═╡ cbf79a6d-981d-49da-9a5f-1bdd9bc53e55
md"
## Network
- Baseline $(sum(net_baselines)/length(net_baselines)) +- $(std(net_baselines)/sqrt(length(length(net_baselines))))

- Max Amp $(sum(net_maxs)/length(net_maxs)) +- $(std(net_maxs)/sqrt(length(length(net_maxs))))

- Min Amp $(sum(net_mins)/length(net_mins)) +- $(std(net_mins)/sqrt(length(length(net_mins))))

- Spike Duration = $(sum(ts[1])/length(ts[1])) +- $(std(ts[1])/sqrt(length(length(ts[1]))))

- Burst Duration = $((sum(ts[2])/length(ts[2]))/1000) +- $((std(ts[2])/sqrt(length(length(ts[2]))))/1000)

- Interburst Interval = $((sum(ts[3])/length(ts[3]))/1000) +- $((std(ts[3])/sqrt(length(length(ts[3]))))/1000)

"

# ╔═╡ 164b2668-f893-46eb-8bb7-025391163840
md"
# Analyze multiple files

"

# ╔═╡ a3266346-58f0-4ef5-b747-c88567da7a0d
begin
	#Open multiple files
	target_folder = "E:\\Data\\Jordans_Patch_Data\\UsuableData"
	paths = target_folder |> parse_abf;	
	
	#Walk through each file analyzing using timescale analysis
	phys_baseline = Float64[]
	phys_max = Float64[]
	phys_min = Float64[]
	phys_spike_durs = Float64[]
	phys_burst_durs = Float64[]
	phys_ibis = Float64[]
	for path in paths
		print("Opening $path ...")
		try
			data_i = extract_abf(path, chs = ["Im_scaled"], continuous = true)
			
			push!( 
				phys_baseline, 
				sum(data_i.data_array[1, :, 1])/length(data_i.data_array[1, :, 1])
			)
			
			push!(phys_min, minimum(data_i))
			push!(phys_max, maximum(data_i))
			#println(min_phys)
			println("Success")
			
			ts_i = NeuroPhys.timescale_analysis(data_i, DURmax = 25)
			
			push!(phys_spike_durs, ts_i[1]...)
			push!(phys_burst_durs, ts_i[2]...)
			push!(phys_ibis, ts_i[3]...)
		catch error
			#println(error)
			println("failed")
		end
	end
	paths
end

# ╔═╡ 99a1108e-9978-4726-9f87-96264794ab1e
md"
## Isolated
- Baseline $(sum(phys_baseline)/length(phys_baseline)) +- $(std(phys_baseline)/sqrt(length(length(phys_baseline))))

- Max Amp $(sum(phys_max)/length(phys_max)) +- $(std(phys_max)/sqrt(length(length(phys_max))))

- Min Amp $(sum(phys_min)/length(phys_min)) +- $(std(phys_min)/sqrt(length(length(phys_min))))

- Spike Duration = $(sum(phys_spike_durs)/length(phys_spike_durs)) +- $(std(phys_spike_durs)/sqrt(length(length(iso_spike_durs))))

- Burst Duration = $((sum(phys_burst_durs)/length(phys_burst_durs))/1000) +- $((std(phys_burst_durs)/sqrt(length(length(phys_burst_durs))))/1000)

- Interburst Interval = $((sum(phys_ibis)/length(phys_ibis))/1000) +- $((std(phys_ibis)/sqrt(length(length(phys_ibis))))/1000)

"

# ╔═╡ d9139dd0-54f3-42cc-b1af-2a8cd809c37e
phys_spike_durs

# ╔═╡ 280a8f63-7746-407a-9a28-c5e83e15273c
begin
	plt_sd = plot(iso_spike_durs, st = :violin, 
		xticks = ([1,2, 3], ["Iso", "Net", "Phys"]), 
		legend = false
	)
	plot!(plt_sd, ts[1], st = :violin)
	plot!(plt_sd, phys_spike_durs, st = :violin)
	
	plt_bd = plot(iso_burst_durs, st = :violin, 
		xticks = ([1,2, 3], ["Iso", "Net", "Phys"]),
		legend = false, 
		yformatter = x -> x/1000
	)
	plot!(plt_bd, ts[2], st = :violin)
	plot!(plt_bd, phys_burst_durs, st = :violin)
	
	plt_ibi = plot(iso_ibis, st = :violin, 
		xticks = ([1,2, 3], ["Iso", "Net", "Phys"]),
		legend = false, 
		yformatter = x -> x/1000
	)
	plot!(plt_ibi, ts[3], st = :violin)
	plot!(plt_ibi, phys_ibis, st = :violin)
	
	plot(plt_sd, plt_bd, plt_ibi, layout = grid(1,3))
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
# ╟─5fcd8df1-9889-41f8-8f3e-c53b00edcf75
# ╠═bd50617a-ad3a-40f1-9d80-3fea8facf83e
# ╠═0379e9a9-d8e7-48f8-aaf6-5e04461e08eb
# ╠═1450b118-bba0-49eb-b4ca-630ee65074fa
# ╠═b96e7e0b-945a-4b3b-8549-c2412c0be3a6
# ╠═60197d1a-fff7-4c55-bdbd-edcf7f534710
# ╟─cbf79a6d-981d-49da-9a5f-1bdd9bc53e55
# ╟─b69cdf52-3dc2-4cbc-989e-4ea838dd70cd
# ╟─ed3f2c93-34a8-454b-9e5a-193ce88e0fe3
# ╟─a422bb0c-f779-44f6-9141-9d8694dd3239
# ╠═4e6c5aee-c93d-4752-ba70-c97d828381d0
# ╟─abc98139-8c4b-488e-ae88-961b6edbb5a6
# ╟─164b2668-f893-46eb-8bb7-025391163840
# ╟─a3266346-58f0-4ef5-b747-c88567da7a0d
# ╟─99a1108e-9978-4726-9f87-96264794ab1e
# ╠═d9139dd0-54f3-42cc-b1af-2a8cd809c37e
# ╟─280a8f63-7746-407a-9a28-c5e83e15273c
