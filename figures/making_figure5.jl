### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ 2756958e-e5b1-11eb-3da6-6f800a6d618c
using Revise

# ╔═╡ 09475831-f6a1-4989-8884-2369e29e16fc
using RetinalChaos

# ╔═╡ 63786ca7-4b1a-4805-be1e-eae87fedab19
using BSON, JLD2, StatsPlots, StatsBase

# ╔═╡ aa88f543-633a-476b-b6cb-441d020dd649
md"
For this one we need to open all of the simulation files I previously opened
"

# ╔═╡ a17974ca-1936-4cc9-bba9-9587012e9b3a
#Lets build a temporary function to convert the timestamps into the correct format
function convert_timestamps(timestamps::Vector{Vector{Tuple{T,T}}}) where T <: Real
	tstamps = Vector{Matrix{Float32}}(undef, size(timestamps))
	for (i, ts_i) in enumerate(timestamps)
		ts_start = map(ts -> ts[1], ts_i)
		ts_end = map(ts -> ts[2], ts_i)
		tstamps[i] = hcat(ts_start, ts_end)
	end
	#Because we dont have the version of the max_interval_algorithim saved 
	burst_tstamps, spd = max_interval_algorithim(tstamps)
	return tstamps, burst_tstamps
end

# ╔═╡ c8d0590c-b9f5-4b2d-95eb-19ecaf816cee
directory = "E:\\Data\\Modelling\\mu_experiment"

# ╔═╡ ac30fece-ca25-403d-b1a3-e4005c26b910
begin
	mus = Float64[]
	avg_spike_dur = Float64[]
	avg_burst_dur = Float64[]
	avg_IBI = Float64[]
	full_spike_dur = []
	full_ISI = []
	full_burst_dur = []
	full_IBI = []
	for (dir_idx, dir) in enumerate(readdir(directory))
		println(dir_idx)
		tstamps_file = "$(directory)\\$(dir)\\timestamps.jld2"
		data_file = "$(directory)\\$(dir)\\data.jld2"
		data_alt = "$(directory)\\$(dir)\\data.bson"
		params_file = "$(directory)\\$(dir)\\params.json"
		p_dict= read_JSON(params_file)
		push!(mus, p_dict[:μ])
		if !isfile(data_alt)
			timestamps = JLD2.load(tstamps_file) #files are accessible
			#Extract the spikes and bursts
			spike_tstamps, burst_tstamps = convert_timestamps(timestamps["Spikes"])
			spike_durs, isi = extract_interval(spike_tstamps, max_duration = 100.0)
			burst_durs, ibi = extract_interval(burst_tstamps, max_duration = 10e5)
		else
			data = BSON.load(data_alt)
			#println(data)
			spike_durs = data["SpikeDurs"]
			burst_durs = data["BurstDurs"]
			isi = data["ISIs"]
			ibi = data["IBIs"]
		end
		push!(full_spike_dur, spike_durs)
		push!(full_ISI, isi)
		push!(full_burst_dur, burst_durs)
		push!(full_IBI, ibi)
		push!(avg_spike_dur, sum(spike_durs)/length(spike_durs))
		push!(avg_burst_dur, sum(burst_durs)/length(burst_durs))
	end
	#we need to sort the mus
	n_plots = length(mus)
	idxs = sortperm(mus)
	mus = mus[idxs]
	avg_spike_dur = avg_spike_dur[idxs]
	avg_burst_dur = avg_burst_dur[idxs]
	full_spike_dur = full_spike_dur[idxs]
	full_burst_dur = full_burst_dur[idxs]
	full_ISI = full_ISI[idxs]
	full_IBI = full_IBI[idxs]
end

# ╔═╡ c9db69d3-cc2a-40e6-9f88-c28615868f9d
begin
	plt_spike_ridge = plot(
		size = (2500,1000),
		bottom_margin = 10.0mm,
		cbar = false,
		yaxis = false, grid = false,
		xlims = (0.0, Inf),
		xticks = (0:0.5:n_plots/2, round.(Int64, mus.*100)), xlabel = "Mu (%)"
	)

	for i in reverse(1:n_plots)
		hfit = fit(Histogram, full_spike_dur[i], LinRange(0, 20, 100))
		weights = hfit.weights/maximum(hfit.weights)
		edges = collect(hfit.edges[1])[1:length(weights)]
	
		if i == 1
			plot!(plt_spike_ridge[1], c = :black, 
				left_margin = 20.0mm,
				weights, edges, 
				fill_z = i,
				fillrange = [1.0, 1.0], fillcolor = :jet, 
				ylabel = "Spike Durations (ms)", label = ""
			)
		else
			plot!(plt_spike_ridge[1], c = :black, label = "",
				#right_margin = -10.0mm, 
				weights.+((i-1)/2), edges,
				fill_z = i,
				fillrange = [0.0, 0.0], fillcolor = :jet, 
			)
		end
	end
	#plot!(plt_spike_ridge, 
	#	collect(0:n_plots-1)./2, avg_spike_dur, 
	#	marker = :circle, c = :red
	#)
	plt_spike_ridge
end

# ╔═╡ a32fbb84-c681-492c-a37c-3867db567b7a
begin
	plt_ISI_ridge = plot(
		size = (2500,1000),
		bottom_margin = 10.0mm,
		cbar = false,
		yaxis = false, grid = false,
		xlims = (0.0, Inf),
		xticks = (0:0.5:n_plots/2, round.(Int64, mus.*100)), xlabel = "Mu (%)"
	)

	for i in reverse(1:n_plots)
		hfit = fit(Histogram, full_ISI[i], LinRange(0, 200, 100))
		weights = hfit.weights/maximum(hfit.weights)
		edges = collect(hfit.edges[1])[1:length(weights)]
	
		if i == 1
			plot!(plt_ISI_ridge[1], c = :black, 
				left_margin = 20.0mm,
				weights, edges, 
				fill_z = i,
				fillrange = [1.0, 1.0], fillcolor = :jet, 
				ylabel = "Interspike Interval (ISI) (ms)", label = ""
			)
		else
			plot!(plt_ISI_ridge[1], c = :black, label = "",
				#right_margin = -10.0mm, 
				weights.+((i-1)/2), edges,
				fill_z = i,
				fillrange = [0.0, 0.0], fillcolor = :jet, 
			)
		end
	end
	#plot!(plt_spike_ridge, 
	#	collect(0:n_plots-1)./2, avg_spike_dur, 
	#	marker = :circle, c = :red
	#)
	plt_ISI_ridge
end

# ╔═╡ 9d2d5214-9d52-4048-a99f-71fc5403cdff
begin
	plt_burst_ridge = plot(
		size = (2500,1000),
		bottom_margin = 10.0mm,
		yaxis = false, grid = false,
		cbar = false,
		xlims = (0.0, Inf),
		xticks = (0:0.5:n_plots/2, round.(Int64, mus.*100)), xlabel = "Mu (%)"
	)
	max_val = map(fbd -> maximum(fbd), full_burst_dur) |> maximum
	for i in reverse(1:n_plots)
		hfit = fit(Histogram, full_burst_dur[i], LinRange(0, 2000, 100))
		weights = hfit.weights/maximum(hfit.weights)
		edges = collect(hfit.edges[1])[1:length(weights)]./1000
	
		if i == 1
			plot!(plt_burst_ridge[1], c = :black, 
				#left_margin = 20.0mm,
				weights, edges, 
				fill_z = i,
				fillrange = [2.0, 2.0], fillcolor = :jet, 
				ylabel = "Burst Durations (s)", label = ""
			)
		else
			plot!(plt_burst_ridge[1], c = :black, label = "",
				#right_margin = -10.0mm, 
				weights.+((i-1)/2), edges,
				fill_z = i,
				fillrange = [2.0, 2.0], fillcolor = :jet, 
			)
		end
	end
	#plot!(plt_burst_ridge, 
	#	collect(0:n_plots-1)./2, avg_burst_dur./1000, 
	#	marker = :circle, c = :red
	#)
	plt_burst_ridge
end

# ╔═╡ 29487316-45cd-4630-8472-72c1e8cae38a
begin
	plt_IBI_ridge = plot(
		size = (2500,1000),
		bottom_margin = 10.0mm,
		yaxis = false, grid = false,
		cbar = false,
		xlims = (0.0, Inf),
		xticks = (0:0.5:n_plots/2, round.(Int64, mus.*100)), xlabel = "Mu (%)"
	)

	for i in reverse(1:n_plots)
		hfit = fit(Histogram, full_IBI[i], LinRange(0, 30e3, 100))
		weights = hfit.weights/maximum(hfit.weights)
		edges = collect(hfit.edges[1])[1:length(weights)]./1000
	
		if i == 1
			plot!(plt_IBI_ridge[1], c = :black, 
				left_margin = 20.0mm,
				weights, edges, 
				fill_z = i,
				fillrange = [2.0, 2.0], fillcolor = :jet, 
				ylabel = "Interburst Interval (s)", label = ""
			)
		else
			plot!(plt_IBI_ridge[1], c = :black, label = "",
				#right_margin = -10.0mm, 
				weights.+((i-1)/2), edges,
				fill_z = i,
				fillrange = [2.0, 2.0], fillcolor = :jet, 
			)
		end
	end
	#plot!(plt_burst_ridge, 
	#	collect(0:n_plots-1)./2, avg_burst_dur./1000, 
	#	marker = :circle, c = :red
	#)
	plt_IBI_ridge
end

# ╔═╡ c7a91138-02db-42a6-a42d-138fef07836d
begin #Load all the paths
	load_path0 = "E:\\Data\\Modelling\\mu_0\\"
	p_dict = read_JSON("$(load_path0)\\params.json", is_type = Dict{Symbol, Float32})
	u_dict = read_JSON("$(load_path0)\\iconds.json", is_type = Dict{Symbol, Float32}) 
	#Mu = 0%
	sol_mu0 = load_model(load_path0, p_dict, u_dict, gpu = false)
	#Mu = 12%
	load_path12 = "E:\\Data\\Modelling\\mu_12\\"
	sol_mu12 = load_model(load_path12, p_dict, u_dict, gpu = false)
	#Mu = 25%
	load_path25 = "E:\\Data\\Modelling\\mu_25\\"
	sol_mu25 = load_model(load_path25, p_dict, u_dict, gpu = false)
	#Mu = 50%
	load_path50 = "E:\\Data\\Modelling\\mu_50\\"
	sol_mu50 = load_model(load_path50, p_dict, u_dict, gpu = false)
end	

# ╔═╡ 6c2abaec-8214-4c20-94ed-6cda1f216ee0
begin #Heatmap range
	nx, ny = (125, 125)
	plt_grid = plot(layout = grid(4,4), 
		yticks = false, xticks = false, yaxis = false, xaxis = false,
		ratio = :equal, cbar = false,
		xlims = (0, nx), ylims = (0, ny), margin = 0.0mm, 
		size = (500,500)
	)
	plt_trace = plot(layout = grid(4,1), grid = false, xformatter = x -> x/1000)
	#Plot the 0% mu
	
	plot!(plt_trace[1], sol_mu0, vars = collect(1:7), label = false, 
		c = :rainbow, line_z = collect(1:7)', clims = (1,7), cbar = false,
		xlabel = "", xticks = false, 
		ylabel = "Vt (mV)", ylims = (-80.0, 0.0), 
		yticks = LinRange(-80.0,0.0, 4), yformatter = y -> round(Int64, y)
	)
	for (idx, frame) in enumerate(LinRange(40e3, 50e3, 4))
		grid_sim = reshape(sol_mu0(frame), nx, ny)
		plot!(plt_grid[1,idx], grid_sim, 
			c = :curl, st = :heatmap, clims = (-70.0, 0.0)
		)
		annotate!(plt_grid[1,idx], round(Int64,125/2), 20, "t", :white)
	end
	
	plot!(plt_trace[2], sol_mu12, vars = collect(1:7), label = false, 
		c = :rainbow, line_z = collect(1:7)', clims = (1,7), cbar = false,
		xlabel = "", xticks = false, 
		ylabel = "Vt (mV)", ylims = (-80.0, 0.0), 
		yticks = LinRange(-80.0,0.0, 4), yformatter = y -> round(Int64, y)
	)
	for (idx, frame) in enumerate(LinRange(10e3, 20e3, 4))
		grid_sim = reshape(sol_mu12(frame), nx, ny)
		plot!(plt_grid[2,idx], grid_sim, 
			c = :curl, st = :heatmap, clims = (-70.0, 0.0)
		)
		annotate!(plt_grid[2,idx], round(Int64,125/2), 20, "t", :white)
	end
	
	plot!(plt_trace[3], sol_mu25, vars = collect(1:7), label = false, 
		c = :rainbow, line_z = collect(1:7)', clims = (1,7), cbar = false,
		xlabel = "", xticks = false, 
		ylabel = "Vt(mV)", ylims = (-80.0, 0.0), 
		yticks = LinRange(-80.0,0.0, 4), yformatter = y -> round(Int64, y)
	)
	for (idx, frame) in enumerate(LinRange(4e4, 5e4, 4))
		grid_sim = reshape(sol_mu25(frame), nx, ny)
		plot!(plt_grid[3,idx], grid_sim, 
			c = :curl, st = :heatmap, clims = (-70.0, 0.0)
		)
		annotate!(plt_grid[3,idx], round(Int64,125/2), 20, "t", :white)
	end
	
	plot!(plt_trace[4], sol_mu50, vars = collect(1:7), label = false, 
		c = :rainbow, line_z = collect(1:7)', clims = (1,7), cbar = false,
		xlabel = "Time (s)", 
		ylabel = "Vt (mV)", ylims = (-80.0, 0.0), 
		yticks = LinRange(-80.0,0.0, 4), yformatter = y -> round(Int64, y)
	)
	
	for (idx, frame) in enumerate(LinRange(5e3, 15e3, 4))
		grid_sim = reshape(sol_mu50(frame), nx, ny)
		plot!(plt_grid[4,idx], grid_sim, 
			c = :curl, st = :heatmap, clims = (-70.0, 0.0)
		)
		annotate!(plt_grid[4,idx], round(Int64,125/2), 20, "t", :white)
	end
	
	fig5_BC = plot(plt_trace, plt_grid, size = (1000, 500))
end

# ╔═╡ 19cc489d-8df4-4d32-93f0-a0c842fc1f64
fig5 = plot(
	plt_burst_ridge, fig5_BC, 
	
	bottom_margin = 0.0mm,
	layout = grid(2,1), size = (1000, 1000), dpi = 500
	) 

# ╔═╡ e3d313db-b96f-4213-b695-16ddec786806
begin
	save_fig = "E:\\Projects\\2021_Modelling_Paper\\Figures\\"
	savefig(fig5, "$(save_fig)\\Fig5_Bursts_and_nullification.png")
end

# ╔═╡ Cell order:
# ╠═2756958e-e5b1-11eb-3da6-6f800a6d618c
# ╠═09475831-f6a1-4989-8884-2369e29e16fc
# ╠═63786ca7-4b1a-4805-be1e-eae87fedab19
# ╟─aa88f543-633a-476b-b6cb-441d020dd649
# ╠═a17974ca-1936-4cc9-bba9-9587012e9b3a
# ╟─c8d0590c-b9f5-4b2d-95eb-19ecaf816cee
# ╠═ac30fece-ca25-403d-b1a3-e4005c26b910
# ╟─c9db69d3-cc2a-40e6-9f88-c28615868f9d
# ╟─a32fbb84-c681-492c-a37c-3867db567b7a
# ╟─9d2d5214-9d52-4048-a99f-71fc5403cdff
# ╠═29487316-45cd-4630-8472-72c1e8cae38a
# ╟─c7a91138-02db-42a6-a42d-138fef07836d
# ╟─6c2abaec-8214-4c20-94ed-6cda1f216ee0
# ╠═19cc489d-8df4-4d32-93f0-a0c842fc1f64
# ╠═e3d313db-b96f-4213-b695-16ddec786806
