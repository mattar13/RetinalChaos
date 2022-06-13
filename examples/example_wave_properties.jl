using Revise
using RetinalChaos
import RetinalChaos.∇
#Extract the parameters
pars_dict = read_JSON("params/GABA_params.json")
De = pars_dict[:De]
Di = pars_dict[:Di]

#%% Construct the lattice sim
tseries = 1:100
nx, ny = (50, 50)
simE = zeros(nx, ny, tseries);
simI = zeros(nx, ny, tseries);
dE = zeros(nx, ny)
dI = zeros(nx, ny)
dXi = (1.9, 0.1)
dYi = (1.0, 1.0)
# Run the sim lattice
for t in tseries
     println(t)
     if t == 1 #Boundary condition I 
          simE[25, 25, 1] = 100.0
          simI[25, 25, 1] = 100.0
     else
          E = simE[:,:,t-1]
          I = simI[:,:,t-1]
          ∇(dE, E, De)
          ∇(dI, I, Di; dX = dXi, dY = dYi)
          simE[:,:,t] = E+dE
          simI[:,:,t] = I+dI
     end
end

#%% Animate the sim lattice
hmpE = plot(ratio=:equal, grid=false,
     xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny),
     c=:curl, clims = (0.0, 10.0))
hmpI = plot(ratio=:equal, grid=false,
     xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny),
     c=:curl, clims = (0.0, 10.0))

plt_A = plot(hmpE, hmpI, layout = 2)
anim = @animate for t = tseries
     println("[$(now())]: Animating simulation...")
     println(t)
     frameE = simE[:, :, t]
     frameI = simI[:, :, t]
     plot!(plt_A[1], frameE, st = :heatmap)
     plot!(plt_A[2], frameI, st = :heatmap)
     plt_A
end
gif(anim, "animation.gif", fps=1000.0 / animate_dt)