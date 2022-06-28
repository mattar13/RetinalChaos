using Revise
using RetinalChaos
import RetinalChaos: calculate_threshold, get_timestamps, max_interval_algorithim

#%% Running the data analysis on a single cell simulation
#Step 1: Import the initial conditions, paramaters and tspan
conds_dict = read_JSON("params\\conds.json")
u0 = conds_dict |> extract_dict
pars_dict = read_JSON("params\\params.json")
p = pars_dict |> extract_dict
tspan = (0.0, 300e3)
#Step 2: set up and solve the problem
prob = SDEProblem(T_SDE, noise, u0, tspan, p)
@time solSDE = solve(prob, SOSRI()); #So far the best method is SOSRI
#plot(sol, vars = [1,5,8], layout=(3, 1))

#%% Step 2: Conduct the analysis
thresh = calculate_threshold(solSDE, idxs=-1)
spike_tstamps = get_timestamps(solSDE, thresh, (solSDE.t[1], solSDE.t[end]), idxs=-1)
spike_durs, isi = extract_interval(spike_tstamps, max_duration=100, max_interval=100)
burst_tstamps, SPB = max_interval_algorithim(spike_tstamps)
burst_durs, ibi = extract_interval(burst_tstamps)

#%% Lets plot things to see more
plt_a = plot(solSDE, vars=[1]) #Plot the solution

# Plot each parameter in a histogram
plt_b = plot(layout=(2, 2))

#Plot spike durations (Best range is )
hfit = fit(Histogram, spike_durs, LinRange(1, 25, 50))
weights = hfit.weights / maximum(hfit.weights)
edges = collect(hfit.edges[1])[1:length(hfit.weights)]
plot!(plt_b[1, 1], edges, weights)

#Plot interspike interval
hfit = fit(Histogram, isi, LinRange(1, 100, 50))
weights = hfit.weights / maximum(hfit.weights)
edges = collect(hfit.edges[1])[1:length(hfit.weights)]
plot!(plt_b[1, 2], edges, weights)

#Plot burst durations
hfit = fit(Histogram, burst_durs, LinRange(100, 2000, 100))
weights = hfit.weights / maximum(hfit.weights)
edges = collect(hfit.edges[1])[1:length(hfit.weights)]
plot!(plt_b[2, 1], edges, weights)

#Plot burst durations
hfit = fit(Histogram, ibi, LinRange(10e3, 100e3, 50))
weights = hfit.weights / maximum(hfit.weights)
edges = collect(hfit.edges[1])[1:length(hfit.weights)]
plot!(plt_b[2, 2], edges, weights)

plot(plt_a, plt_b, layout=(1, 2))

#%% Running the data analysis on a 2D wave sim
#1) Set up the model
nx = ny = 50
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
@time solPDE = solve(prob, SOSRI(), abstol=2e-2, reltol=0.2, maxiters=1e7, progress=true, progress_steps=1)#, save_idxs=[1:(nx*ny)...])
println("  [$(now())]: Completed")

#%% Step 2: Conduct the analysis
thresh = calculate_threshold(solPDE, idxs=1)
spike_tstamps = get_timestamps(solPDE, idxs=1)
spike_durs, isi = extract_interval(spike_tstamps, max_duration=100, max_interval=100)
burst_tstamps, SPB = max_interval_algorithim(spike_tstamps)
burst_durs, ibi = extract_interval(burst_tstamps)

#Pull out an example trace
t = solPDE.t
vt1 = solPDE(t, idxs=1)
vt2 = solPDE(t, idxs=10)
vt3 = solPDE(t, idxs=100)
vt4 = solPDE(t, idxs=1000)
#%% Lets plot things to see more
plt_a = plot(t, vt1) #Plot the solution
plot!(t, vt2)
plot!(t, vt3)
plot!(t, vt4)

# Plot each parameter in a histogram
plt_b = plot(layout=(2, 2))

#Plot spike durations (Best range is )
hfit = fit(Histogram, spike_durs, LinRange(1, 100, 50))
weights = hfit.weights / maximum(hfit.weights)
edges = collect(hfit.edges[1])[1:length(hfit.weights)]
plot!(plt_b[1, 1], edges, weights)

#Plot interspike interval
hfit = fit(Histogram, isi, LinRange(1, 100, 50))
weights = hfit.weights / maximum(hfit.weights)
edges = collect(hfit.edges[1])[1:length(hfit.weights)]
plot!(plt_b[1, 2], edges, weights)

#Plot burst durations
hfit = fit(Histogram, burst_durs, LinRange(100, 2000, 100))
weights = hfit.weights / maximum(hfit.weights)
edges = collect(hfit.edges[1])[1:length(hfit.weights)]
plot!(plt_b[2, 1], edges, weights)

#Plot burst durations
hfit = fit(Histogram, ibi, LinRange(10e3, 100e3, 50))
weights = hfit.weights / maximum(hfit.weights)
edges = collect(hfit.edges[1])[1:length(hfit.weights)]
plot!(plt_b[2, 2], edges, weights)

plot(plt_a, plt_b, layout=(1, 2))


