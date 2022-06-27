using Revise, Dates
using RetinalChaos
using ABFReader
using DataFrames, XLSX, BSON
using Statistics, StatsPlots, StatsBase



#==
font_title = font("Arial", 24) #title font and size
font_axis = font("Arial", 12) #axis font and size
font_legend = font("Arial", 8) #legend font and size
gr(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)

params_file = joinpath(param_path, "params.json")
conds_file = joinpath(param_path, "conds.json")

#Set up the initial conditions
u0 = read_JSON(conds_file)|> extract_dict;
#Set up the parameters 
p_dict = read_JSON(params_file) 
#p_dict[:I_app] = 10.0
p_dict[:g_ACh] = 0.0
p0 = p_dict |> extract_dict
#Set up the Time span
tspan = (0.0, 300e3)
#set up the problem to solve
prob = SDEProblem(T_sde, noise, u0 , tspan, p0);
#Load the problem into a ensemble simulation
n_sims = 20
test_rng = LinRange(1.0, 40.0, n_sims)
prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, :g_TREK, test_rng)
#make the ensemble problem and only save the voltage
simProb = EnsembleProblem(prob, prob_func = prob_func);

# ╔═╡ 2c12afea-8238-45d9-9d00-4aa5692a6e14
@time sim = solve(simProb, SOSRI(), trajectories = n_sims, save_idxs = 1, EnsembleThreads());

var_baseline = Float64[]
var_max = Float64[]
var_min = Float64[]
var_spike_durs = Float64[]
var_isis = Float64[]
var_burst_durs = Float64[]
var_ibis = Float64[]

p_var_trace = plot()
p_varSD_hist = plot()
p_varISD_hist = plot()
p_varBD_hist = plot()
p_varIBD_hist = plot()

for (i, sol) in enumerate(sim)
	thresh = RetinalChaos.calculate_threshold(sol)
	spike_tstamps_var = RetinalChaos.get_timestamps(sol, thresh)
	burst_tstamps_var = RetinalChaos.max_interval_algorithim(spike_tstamps_var)
	println(!isempty(spike_tstamps_var))
	println(!isempty(burst_tstamps_var))
	if !isempty(spike_tstamps_var)
		spike_durs, isi = extract_interval(spike_tstamps_var)
		
		hfit = fit(Histogram, spike_durs, LinRange(5, 50, 100))
		weights = hfit.weights/maximum(hfit.weights)
		edges = collect(hfit.edges[1])[1:length(hfit.weights)]
		plot!(p_varSD_hist, c = :black, 
				weights.+(i-1), edges,
				fill_z = i,
				fillrange = [0.0, 0.0], fillcolor = :jet, 
				ylims = (5.0, 50.0),
				xticks = (
					1:5:length(test_rng), 
					round.(test_rng[1:5:end])
					),
				ylabel = "Spikes (ms)", label = "", cbar = false
			)
		
		hfitI = fit(Histogram, isi, LinRange(1, 100, 100))
		weightsI = hfitI.weights/maximum(hfitI.weights)
		edgesI = collect(hfitI.edges[1])[1:length(hfitI.weights)]
		plot!(p_varISD_hist, c = :black, 
				weightsI.+(i-1), edgesI,
				fill_z = i,
				fillrange = [0.0, 0.0], fillcolor = :jet, 
				xticks = (
					1:5:length(test_rng), 
					round.(test_rng[1:5:end])
					),
				ylabel = "ISI (ms)", label = "", cbar = false
			)
		
	end
	if !isempty(burst_tstamps_var)
		burst_durs, ibi = extract_interval(burst_tstamps_var[1])
		
		hfit = fit(Histogram, burst_durs, LinRange(100, 3000, 50))
		weights = hfit.weights/maximum(hfit.weights)
		edges = collect(hfit.edges[1])[1:length(hfit.weights)]
		plot!(p_varBD_hist, c = :black, 
				weights.+(i-1), edges,
				yformatter = y -> y/1000,
				fill_z = i,
				fillrange = [0.0, 0.0], fillcolor = :jet, 
				xticks = (
					1:5:length(test_rng), 
					round.(test_rng[1:5:end])
					),
				ylabel = "Bursts (s)", label = "", cbar = false
			)
	
		hfitI = fit(Histogram, ibi, LinRange(10, 100e3, 100))
		weightsI = hfitI.weights/maximum(hfitI.weights)
		edgesI = collect(hfitI.edges[1])[1:length(hfitI.weights)]
		plot!(p_varIBD_hist, c = :black, 
				weightsI.+(i-1), edgesI,
				yformatter = y -> y/1000,
				fill_z = i,
				ylims = (10, 100e3),
				fillrange = [10.0, 10.0], fillcolor = :jet, 
				xticks = (
					1:5:length(test_rng), 
					round.(test_rng[1:5:end])
					),
				ylabel = "IBI (s)", label = "", cbar = false
			)
	end

	plot!(p_var_trace, sol, vars = 1, 
		line_z = test_rng[i], c = :jet, label = ""
	)
end
p_hists = plot(
	p_varSD_hist, p_varISD_hist, 
	p_varBD_hist, p_varIBD_hist, 
	layout = grid(2,2)
	)
plot(p_var_trace, p_hists, layout = grid(2,1), size = (1000,1000))

data_loc = "E:\\Data\\Modelling"

load_path_iso = "$(data_loc)\\μ_0\\"
#This takes a really long time
#IsoSol = RetinalChaos.load_solution(load_path_iso)
iso_data = BSON.load("$(load_path_iso)\\data.bson")
iso_spike_durs = iso_data["SpikeDurs"]
iso_isi = iso_data["ISIs"]
iso_burst_durs = iso_data["BurstDurs"]
iso_ibi = iso_data["IBIs"]

#iso_baseline = map(i -> sum(IsoSol[i, :])/length(IsoSol), 1:size(IsoSol,1))
#iso_max = map(i -> maximum(IsoSol[i, :]), 1:size(IsoSol,1))
#iso_min = map(i -> minimum(IsoSol[i, :]), 1:size(IsoSol,1))

#Load the simulation data (This may take really long for such large sims)
load_path = "$(data_loc)\\μ_18\\" #run one with 18%
#NetSol = RetinalChaos.load_solution(load_path)
net_data = BSON.load("$(load_path)\\data.bson")
net_spike_durs = net_data["SpikeDurs"]
net_isi = net_data["ISIs"]
net_burst_durs = net_data["BurstDurs"]
net_ibi = net_data["IBIs"]

#net_baseline = map(i -> sum(NetSol[i, :])/length(NetSol), 1:size(NetSol,1))
#net_max = map(i -> maximum(NetSol[i, :]), 1:size(NetSol,1))
#net_min = map(i -> minimum(NetSol[i, :]), 1:size(NetSol,1))

#Load the example trace I collected
target_file = "$(data_loc)\\Patching\\2019_11_03_Patch\\Animal_2\\Cell_3\\19n03042.abf"
data = readABF(target_file, channels = ["Vm_prime4"], stimulus_name = nothing)
data - 25.0
spike_timestamps = NeuroPhys.get_timestamps(data) #account for s -> ms
burst_timestamps = NeuroPhys.max_interval_algorithim(spike_timestamps)
spike_durs, isi = extract_interval(spike_timestamps)
burst_durs, ibi = extract_interval(burst_timestamps[1])

#%% Load all physiological data 
target_folder = "$(data_loc)\\Patching\\Jordans_Patch_Data\\UsuableData"
paths = target_folder |> parse_abf;	

#Walk through each file analyzing using timescale analysis
phys_baseline = Float64[]
phys_max = Float64[]
phys_min = Float64[]
phys_spike_durs = Float64[]
phys_isis = Float64[]
phys_burst_durs = Float64[]
phys_ibis = Float64[]
for path in paths
	print("Opening $path ...")
	data_i = readABF(path, 
		channels = ["Im_scaled"], 
		stimulus_name = nothing, flatten_episodic = true
	)
	spike_timestampsi = NeuroPhys.get_timestamps(data_i)
	burst_timestampsi = NeuroPhys.max_interval_algorithim(spike_timestampsi, DURmin = 100.0)
	spike_durs, isi = extract_interval(spike_timestampsi)
	burst_durs, ibi = extract_interval(burst_timestampsi[1])	
	push!( 
		phys_baseline, 
		sum(data_i.data_array[1, :, 1])/length(data_i.data_array[1, :, 1])
	)

	push!(phys_min, minimum(data_i))
	push!(phys_max, maximum(data_i))
	#println(min_phys)
	
		
	push!(phys_spike_durs, (spike_durs)...)
	push!(phys_isis, (isi)...)
	push!(phys_burst_durs, (burst_durs)...)
	push!(phys_ibis, (ibi)...)
	println("Success")
end
end

phys_color = :green
iso_color = :blue
net_color = :orange


thresh = ABFReader.calculate_threshold(data)
tidxs = round(Int64, 150e3/data.dt):10:round(Int64, 210e3/data.dt)
tseries = (data.t[tidxs].-data.t[tidxs[1]])

fig_traces1 = plot(tseries, data[1, tidxs, 1], 
	legend = :bottomright, label = "Physiological Data", 
	c = phys_color, xticks = false, 
	ylims = (-120.0, 0.0), yticks = LinRange(-100.0, 0.0, 3)
	
)
fig_traces2 = plot(IsoSol.t, IsoSol[30,:], legend = :bottomright, 
	label = "Isolated Simulation",
	xlims = (0.0, 60e3), xticks = false,
	ylims = (-120.0, 0.0), yticks = LinRange(-100.0, 0.0, 3), ylabel = "Vₜ (mV)",
	c = iso_color, 
)
fig_traces3 = plot(NetSol.t, NetSol[25,:], label = "Network Simulation",
	xlims = (0.0, 60e3),xlabel = "Time (s)",
	ylims = (-120.0, 0.0), yticks = LinRange(-100.0, 0.0, 3),
	c = net_color, 
)
fig_traces = plot(fig_traces1, fig_traces2, fig_traces3, 
	layout = grid(3,1), grid = false, lw = 2.0,
	xformatter = x -> x/1000, 
)

n_hist_bins = 10

# ╔═╡ 97a95cba-d8bd-4a39-875b-1123f83d6889
begin
	#burst phys
	hfit_sp = fit(Histogram, phys_spike_durs, LinRange(5, 20, n_hist_bins))
	weights_sp = hfit_sp.weights/maximum(hfit_sp.weights)
	edges_sp = collect(hfit_sp.edges[1])[1:length(hfit_sp.weights)]
	
	p1s = plot(edges_sp, weights_sp, label = "Physiological Data", 
		st = :bar, c = phys_color, xticks = false,
	)
	
	#burst isolated
	hfit_si = fit(Histogram, iso_spike_durs, LinRange(5, 10, 100))
	weights_si = hfit_si.weights/maximum(hfit_si.weights)
	edges_si = collect(hfit_si.edges[1])[1:length(hfit_si.weights)]
	
	p2s = plot(edges_si, weights_si, label = "Isolated Simulation",legend = :topleft,
		ylabel = "Probability of Occurrance",
		c = :black, xticks = false,
		fillrange = [0.0, 0.0], fillcolor = iso_color, 
	)
	
	#burst network
	hfit_sn = fit(Histogram, net_spike_durs, LinRange(5, 10, 100))
	weights_sn = hfit_sn.weights/maximum(hfit_sn.weights)
	edges_sn = collect(hfit_sn.edges[1])[1:length(hfit_sn.weights)]
	
	p3s = plot(edges_sn, weights_sn, label = "Network Simulation", legend = :topleft,
		xlabel = "Spike Duration (ms)",
		c = :black, 
		fillrange = [0.0, 0.0], fillcolor = net_color, 
	)
	
	pSD = plot(p1s, p2s, p3s, sharex = true, legend = false,
		dpi = 300, layout = grid(3,1), grid = false
	)
end

# ╔═╡ cf21cd4f-8fb7-4328-b781-6b965abb159d
begin
	#burst phys
	hfit_ISIp = fit(Histogram, phys_isis, LinRange(1, 150.0, n_hist_bins))
	weights_ISIp = hfit_ISIp.weights/maximum(hfit_ISIp.weights)
	edges_ISIp = collect(hfit_ISIp.edges[1])[1:length(hfit_ISIp.weights)]
	
	p1ISI = plot(edges_ISIp, weights_ISIp, label = "Physiological Data", 
		st = :bar, c = phys_color, xticks = false,
		#xformatter = x -> x/1000
	)
	
	#burst isolated
	hfit_ISIi = fit(Histogram, iso_isi, LinRange(1, 150.0, 100))
	weights_ISIi = hfit_ISIi.weights/maximum(hfit_ISIi.weights)
	edges_ISIi = collect(hfit_ISIi.edges[1])[1:length(hfit_ISIi.weights)]
	
	p2ISI = plot(edges_ISIi, weights_ISIi, label = "Isolated Simulation", 
		ylabel = "Probability of Occurrance",
		legend = :topright, xticks = false,
		c = :black, 
		fillrange = [0.0, 0.0], fillcolor = iso_color, 
		#xformatter = x -> x/1000
	)
	
	#burst network
	hfit_ISIn = fit(Histogram, net_isi, LinRange(1, 150.0, 100))
	weights_ISIn = hfit_ISIn.weights/maximum(hfit_ISIn.weights)
	edges_ISIn = collect(hfit_ISIn.edges[1])[1:length(hfit_ISIn.weights)]
	
	p3ISI = plot(edges_ISIn, weights_ISIn, label = "Network Simulation", 
		legend = :topright, 
		xlabel = "Interspike Interval (s)",
		c = :black, 
		fillrange = [0.0, 0.0], fillcolor = net_color, 
		#xformatter = x -> x/1000
	)
	
	pISI = plot(p1ISI, p2ISI, p3ISI, sharex = true, legend = false,
		dpi = 300, layout = grid(3,1), grid = false
	)
end

# ╔═╡ c76e4565-637b-42bc-9eac-cd036ca9ba1d
begin
	#burst phys
	hfit_bp = fit(Histogram, phys_burst_durs, LinRange(5, 2500, n_hist_bins))
	weights_bp = hfit_bp.weights/maximum(hfit_bp.weights)
	edges_bp = collect(hfit_bp.edges[1])[1:length(hfit_bp.weights)]
	
	p1b = plot(edges_bp, weights_bp, label = "Physiological Data", 
		st = :bar, c = phys_color, xticks = false,
		xformatter = x -> x/1000
	)
	
	#burst isolated
	hfit_bi = fit(Histogram, iso_burst_durs, LinRange(5, 2500, 100))
	weights_bi = hfit_bi.weights/maximum(hfit_bi.weights)
	edges_bi = collect(hfit_bi.edges[1])[1:length(hfit_bi.weights)]
	
	p2b = plot(edges_bi, weights_bi,label = "Isolated Simulation",legend = :topright,
		ylabel = "Probability of Occurrance",
		c = :black, xticks = false,
		fillrange = [0.0, 0.0], fillcolor = iso_color, 
		xformatter = x -> x/1000
	)
	
	#burst network
	hfit_bn = fit(Histogram, net_burst_durs, LinRange(5, 2500, 100))
	weights_bn = hfit_bn.weights/maximum(hfit_bn.weights)
	edges_bn = collect(hfit_bn.edges[1])[1:length(hfit_bn.weights)]
	
	p3b = plot(edges_bn, weights_bn, label = "Network Simulation", legend = :topright, 
		xlabel = "Burst Duration (s)",
		c = :black,
		fillrange = [0.0, 0.0], fillcolor = net_color, 
		xformatter = x -> x/1000
	)
	
	pBD = plot(p1b, p2b, p3b,
		dpi = 300, legend = false, 
		grid = false, layout = grid(3,1), margins = 0.0mm
	)
end

# ╔═╡ 7d86641c-76b8-44d5-85fb-23946c326c1d
begin
	#burst phys
	hfit_IBIp = fit(Histogram, phys_ibis, LinRange(10, 85e3, n_hist_bins))
	weights_IBIp = hfit_IBIp.weights/maximum(hfit_IBIp.weights)
	edges_IBIp = collect(hfit_IBIp.edges[1])[1:length(hfit_IBIp.weights)]
	
	p1IBI = plot(edges_IBIp, weights_IBIp, label = "Physiological Data", 
		st = :bar, c = phys_color, 
		#xlims = (0.0, 80e3),
		xticks = false,
		xformatter = x -> x/1000
		
	)
	
	#burst isolated
	hfit_IBIi = fit(Histogram, iso_ibi, LinRange(1, 80e3, 100))
	weights_IBIi = hfit_IBIi.weights/maximum(hfit_IBIi.weights)
	edges_IBIi = collect(hfit_IBIi.edges[1])[1:length(hfit_IBIi.weights)]
	
	p2IBI = plot(edges_IBIi, weights_IBIi, label = "Isolated Simulation", 
		ylabel = "Probability of Occurrance",
		legend = :topright, 
		c = :black, 
		fillrange = [0.0, 0.0], fillcolor = iso_color,
		xticks = false,
		xformatter = x -> x/1000, margins = 0.0mm,
	)
	
	#burst network
	hfit_IBIn = fit(Histogram, net_ibi, LinRange(1, 80e3, 100))
	weights_IBIn = hfit_IBIn.weights/maximum(hfit_IBIn.weights)
	edges_IBIn = collect(hfit_IBIn.edges[1])[1:length(hfit_IBIn.weights)]
	
	p3IBI = plot(edges_IBIn, weights_IBIn, label = "Network Simulation", 
		legend = :topright, 
		xlabel = "Interburst Interval (s)",
		c = :black, 
		fillrange = [0.0, 0.0], fillcolor = net_color, 
		xformatter = x -> x/1000
	)
	
	pIBI = plot(p1IBI, p2IBI, p3IBI, sharex = true, 
		dpi = 300, legend = false, grid = false, layout = grid(3,1), margins = 0.0mm
	)
end

# ╔═╡ 1afbd83f-424d-4bf0-a53f-79c1c163aff8
begin
	fig4_labels = plot(
		layout = grid(3,1, heights = (0.5, 0.25, 0.25)), 
		xaxis = false, yaxis = false, xticks = false, yticks = false
	)
	annotate!(fig4_labels[1], [0.5], [0.99], "A", font("Sans",24))
	annotate!(fig4_labels[2], [0.5], [0.99], "B", font("Sans",24))
	annotate!(fig4_labels[3], [0.5], [0.99], "C", font("Sans",24))
	
	fig_hists = plot(
			pSD, pISI,
			pBD, pIBI, 
			layout = grid(2,2)
		)
	fig4_graphs = plot(fig_traces, fig_hists, layout = grid(2,1), size = (1000,1000))
	
	fig4 = plot(fig4_labels, fig4_graphs, 
		layout = grid(1,2, widths = (0.05, 0.95), dpi = 300)
	)
end

# ╔═╡ b71f5712-7e95-4ff4-a153-3257ef844b4e
#savefig(fig4, "F:\\Projects\\2021_Modelling_Paper\\Figures\\Fig4_Phys_data.png")

# ╔═╡ 63eebfe4-75b3-49ef-a2c4-2aa832cfa25b
begin 
	#Make a dataframe for all of the data
	df = DataFrame(
		DataFrom = String[], 
		Baseline = Float64[], Baseline_SEM = Float64[], 
		Max_Amp = Float64[], Max_Amp_SEM = Float64[], 
		Min_Amp = Float64[], Min_Amp_SEM = Float64[], 
		Spike_Duration = Float64[], Spike_Duration_SEM = Float64[], 
		Interspike_Interval = Float64[], Interspike_Interval_SEM = Float64[], 
		Burst_Duration = Float64[], Burst_Duration_SEM = Float64[], 
		Interburst_Interval = Float64[], Interburst_Interval_SEM = Float64[]
	)
	push!(df, [
			"Isolated Simulation", 
			sum(iso_baseline)/length(iso_baseline),
			std(iso_baseline)/sqrt(length(iso_baseline)),

			sum(iso_max)/length(iso_max),
			std(iso_max)/sqrt(length(iso_max)),

			sum(iso_min)/length(iso_min),
			std(iso_min)/sqrt(length(iso_min)),

			sum(iso_spike_durs)/length(iso_spike_durs) , 		
			std(iso_spike_durs)/sqrt(length(iso_spike_durs)), 
			
			sum(iso_isi)/length(iso_isi),  	
			std(iso_isi)/sqrt(length(iso_isi)),
			
			(sum(iso_burst_durs)/length(iso_burst_durs))/1000, 	
			(std(iso_burst_durs)/sqrt(length(iso_burst_durs)))/1000, 

			(sum(iso_ibi)/length(iso_ibi))/1000,  	
			(std(iso_ibi)/sqrt(length(iso_ibi)))/1000
			]
		)
	push!(df, [
			"Network Simulation", 
			sum(net_baseline)/length(net_baseline),
			std(net_baseline)/sqrt(length(net_baseline)),

			sum(net_max)/length(net_max),
			std(net_max)/sqrt(length(net_max)),

			sum(net_min)/length(net_min),
			std(net_min)/sqrt(length(net_min)),

			sum(net_spike_durs)/length(net_spike_durs) , 		
			std(net_spike_durs)/sqrt(length(net_spike_durs)), 

			sum(net_isi)/length(net_isi),  	
			std(net_isi)/sqrt(length(net_isi)),

			(sum(net_burst_durs)/length(net_burst_durs))/1000, 	
			(std(net_burst_durs)/sqrt(length(net_burst_durs)))/1000, 

			(sum(net_ibi)/length(net_ibi))/1000,  	
			(std(net_ibi)/sqrt(length(net_ibi)))/1000
			]
		)
	push!(df, [
			"Physiological Data", 
			sum(phys_baseline)/length(phys_baseline),
			std(phys_baseline)/sqrt(length(phys_baseline)),

			sum(phys_max)/length(phys_max),
			std(phys_max)/sqrt(length(phys_max)),

			sum(phys_min)/length(phys_min),
			std(phys_min)/sqrt(length(phys_min)),

			sum(phys_spike_durs)/length(phys_spike_durs) , 		
			std(phys_spike_durs)/sqrt(length(phys_spike_durs)), 

			sum(phys_isis)/length(phys_isis),  	
			std(phys_isis)/sqrt(length(phys_isis)),

			(sum(phys_burst_durs)/length(phys_burst_durs))/1000, 	
			(std(phys_burst_durs)/sqrt(length(phys_burst_durs)))/1000, 

			(sum(phys_ibis)/length(phys_ibis))/1000,  	
			(std(phys_ibis)/sqrt(length(phys_ibis)))/1000
			]
		)
end

# ╔═╡ d4e7bca9-0605-4234-91cb-e4d8d7eade3d
begin
	save_path = "F:\\Projects\\2021_Modelling_Paper\\data.xlsx"
	if !ispath(save_path)
		XLSX.writetable("F:\\Projects\\2021_Modelling_Paper\\data.xlsx", 
			collect(DataFrames.eachcol(df)), 
			DataFrames.names(df)
		)
	else
		rm(save_path)
		XLSX.writetable("F:\\Projects\\2021_Modelling_Paper\\data.xlsx", 
			collect(DataFrames.eachcol(df)), 
			DataFrames.names(df)
		)
	end
end
==#