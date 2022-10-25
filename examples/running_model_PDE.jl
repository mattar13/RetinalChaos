using RetinalChaos
using Plots

import RetinalChaos.parameters
import RetinalChaos.u0
#Step 1 setup the dimension parameters
dx = 1.0
dy = 1.0
dt = 1.0

tmin = 0.0
tmax = 10000.0
tstops = tmin:dt:tmax
reload_parameters()
parameters[I_app] = 0.0
parameters[g_Na] = 0.0
parameters[De] = 0.1
#For ODE and SDE models these can be set manually
GRID, probSPDE = loadSPDE(u0, parameters,
     tmax=tmax,
     xmax=10.0, ymax=10.0,
     dx=dx, dy=dy
);

# Run the model
sol = solve(probSPDE, saveat=tstops, progress=true, progress_steps=1)

#%% Plot the animation
discrete_t = sol.t
discrete_x = GRID[x]
discrete_y = GRID[y]
anim = @animate for i in 1:100:length(discrete_t)
     println(i)
     v_map = map(d -> sol[d][i], GRID[ê(x, y, t)])
     h1 = heatmap(discrete_x, discrete_y, v_map, clim=(0.0, 0.6), title="$(discrete_t[i])", aspect_ratio=:equal)
end
gif(anim, "sync_voltage.gif", fps=60.0)

#%% Plot the 
v_traces = reshape(sol[GRID[v̂(x, y, t)]], length(discrete_x) * length(discrete_y))
n_traces = reshape(sol[GRID[n̂(x, y, t)]], length(discrete_x) * length(discrete_y))
m_traces = reshape(sol[GRID[m̂(x, y, t)]], length(discrete_x) * length(discrete_y))
h_traces = reshape(sol[GRID[ĥ(x, y, t)]], length(discrete_x) * length(discrete_y))
c_traces = reshape(sol[GRID[ĉ(x, y, t)]], length(discrete_x) * length(discrete_y))
a_traces = reshape(sol[GRID[â(x, y, t)]], length(discrete_x) * length(discrete_y))
b_traces = reshape(sol[GRID[b̂(x, y, t)]], length(discrete_x) * length(discrete_y))
e_traces = reshape(sol[GRID[ê(x, y, t)]], length(discrete_x) * length(discrete_y))
i_traces = reshape(sol[GRID[ê(x, y, t)]], length(discrete_x) * length(discrete_y))

p1 = Plots.plot(discrete_t, v_traces)
p2 = Plots.plot(discrete_t, e_traces)
p3 = Plots.plot(discrete_t, i_traces)
Plots.plot(p1, p2, layout = (2,1))