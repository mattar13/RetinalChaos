using Revise
using RetinalChaos
import RetinalChaos: calculate_threshold, get_timestamps, max_interval_algorithim

#We can do the data analysis on a wave sim
#1) Set up the model
print("[$(now())]: Loading Parameters... ")
conds_dict = read_JSON("params/conds.json")
u0 = extract_dict(conds_dict, t_conds, dims=(nx, ny))
pars_dict = read_JSON("params/params.json")
p = pars_dict |> extract_dict #Extract the parameters
tspan = (0.0, 120e3) #Set the tspan
println(" [$(now())]: Completed")
#2) Warmup the model
print("[$(now())]: Warming up model... ")
prob = SDEProblem(T_PDE, noise, u0, 60e3, p) #run the warmup for 60s
@time warmup = solve(prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, save_everystep=false, progress=true, progress_steps=1)
println("  [$(now())]: Completed")
#3) Run the model and save the solution
print("[$(now())]: Running model... ")
prob = SDEProblem(T_PDE, noise, warmup[end], tspan, p) #Use the last solution from the warmup
@time sol = solve(prob, SROCK1(), dt=1.0, abstol=2e-2, reltol=0.2, maxiters=1e7, progress=true, progress_steps=1)#, save_idxs=[1:(nx*ny)...])
println("  [$(now())]: Completed")

#%% Step 2: Conduct the analysis
thresh = calculate_threshold(sol)[1]
spike_tstamps = get_timestamps(sol; threshold=thresh)
spike_durs, isi = extract_interval(spike_tstamps)
burst_tstamps, SPB = max_interval_algorithim(spike_tstamps)
burst_durs, ibi = extract_interval(burst_tstamps)

#%% Lets plot things to see more
plt_a = plot(sol, vars=[1]) #Plot the solution

#%% Plot each parameter in a histogram
plt_b = plot(layout=(2, 2))

#Plot spike durations
hfit = fit(Histogram, spike_durs, LinRange(1, 10, 50))
weights = hfit.weights / maximum(hfit.weights)
edges = collect(hfit.edges[1])[1:length(hfit.weights)]
plot!(plt_b[1, 1], edges, weights)

#Plot interspike interval
hfit = fit(Histogram, isi, LinRange(1, 10, 50))
weights = hfit.weights / maximum(hfit.weights)
edges = collect(hfit.edges[1])[1:length(hfit.weights)]
plot!(plt_b[1, 2], edges, weights)

#Plot burst durations
hfit = fit(Histogram, burst_durs, LinRange(1, 10, 50))
weights = hfit.weights / maximum(hfit.weights)
edges = collect(hfit.edges[1])[1:length(hfit.weights)]
plot!(plt_b[1, 2], edges, weights)

#Plot burst durations
hfit = fit(Histogram, ibi, LinRange(1, 10, 50))
weights = hfit.weights / maximum(hfit.weights)
edges = collect(hfit.edges[1])[1:length(hfit.weights)]
plot!(plt_b[1, 2], edges, weights)

plot(plt_a, plt_b, layout=(1, 2))