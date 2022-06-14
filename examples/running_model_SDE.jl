using Revise
using RetinalChaos

#Step 1: Import the initial conditions
conds_dict = read_JSON("params\\conds.json")
u0 = conds_dict |> extract_dict

#Step 2: Import the parameters
pars_dict = read_JSON("params\\params.json")
p = pars_dict |> extract_dict

#Step 3: determine the timespan
tspan = (0.0, 300e3)

#Step 4: set up the problem
prob = SDEProblem(T_SDE, noise, u0, tspan, p)

#Step 5: Solve the problem
@time sol = solve(prob, SOSRI()); #So far the best method is SOSRI
plot(sol, vars=[1, 6, 7], layout=(2, 1))

#%% Figure: Comparisons of different SDE algorithims
@time sol = solve(prob, SOSRI());
fig = plot(sol, vars=[1, 7], layout=(2, 1), label="SOSRI")

@time sol = solve(prob, SKenCarp());
plot!(fig, sol, vars=[1, 7], layout=(2, 1), label="SKenCarp")

@time sol = solve(prob, SRIW1());
plot!(fig, sol, vars=[1, 7], layout=(2, 1), label="SRIW1")

@time sol = solve(prob, SROCK1(), dt=0.5);
plot!(fig, sol, vars=[1, 7], layout=(2, 1), label="SROCK1 dt = 0.5ms")

fig
#%% Lets test the efficiency of several algorithims
a = plot()
for dt in 0.1:0.1:1.0
     println("Running for dt = $dt")
     @time sol = solve(prob, SROCK1(), dt=dt)
     plot!(a, sol, vars=[1, 7], layout=(2, 1), label="dt = $dt")
end
a