using RetinalChaos
include("../figures/figure_setup.jl")
import BifurcationKit.@set
#reload_parameters()
#Step 1: Import:
import RetinalChaos.SDEModel #the SDEModel
import RetinalChaos.ODEModel #the ODEModel
import RetinalChaos.u0 #the initial condutions 
import RetinalChaos.parameters #the parameters

#Step 2: determine the timespan
tmin = 0.0
tmax = 500.0

#%% Run the analysis for I_app
#Step 3: set up the problem making sure enable jacobian is true
reload_parameters()
parameters[I_app] = 0.0
#parameters[g_ACh] = 0.0
#parameters[g_GABA] = 0.0
parameters[g_TREK] = 0.0 #Remove the sAHP
parameters[ρi] = u0[i] = 0.0
parameters[ρe] = u0[e] = 0.0

@time ODEProb = ODEProblem(ODEModel, u0, (tmin, tmax), parameters, jac=true)
sol = solve(ODEProb)
initial_plot = plot(sol, idxs=[v, e, i], layout = 3)

# Part A: Determining jacobians and gradients
#Step 1: Make the function and Jacobian
ODEFunc = ODEProb.f
F = (u, p) -> ODEFunc(u, p, 0.0)
J = (u, p) -> ODEFunc.jac(u, p, 0.0)

id_par = indexof(I_app)
par_tm = ODEProb.p

prob = BK.BifurcationProblem(
    F, ODEProb.u0, par_tm, (@lens _[id_par]);
    J=J,
    recordFromSolution=(x, p) -> (V=x[2], N=x[3])
)

# continuation options (Theses need to be tuned)
opts_br = BK.ContinuationPar(
    pMin=-100.0, pMax=100.0, # parameters to have a smooth result
    ds=0.004, dsmax=0.5,# this is to detect bifurcation points precisely with bisection
    detectBifurcation=3, # Optional: bisection options for locating bifurcations
    nInversion=8, maxBisectionSteps=50,
    maxSteps=500,
    nev=3
)

continuation_method = BK.PALC(tangent=BK.Bordered())

br = BK.continuation(prob, continuation_method, opts_br;
    normC=norminf,
    verbosity=3,
    bothside=true, plot=true
)

plot(br)
# This section determines the periodic orbits 
optn_po = BK.NewtonPar(verbose=true, tol=1e-8, maxIter=100)
opts_po_cont = BK.ContinuationPar(
    dsmin=0.04, ds=0.04, dsmax=0.5,
    pMin=0.0, pMax=50.0,
    maxSteps=110,
    newtonOptions=(BK.@set optn_po.tol = 1e-7),
    nev=3, plotEveryStep=2, detectBifurcation=0,
    saveSolEveryStep=1
)

#This is for printing
args_po = (recordFromSolution=(x, p) -> begin
        xtt = BK.getPeriodicOrbit(p.prob, x, @set par_tm[id_par] = p.p)
        return (max=maximum(xtt[2, :]),
            min=minimum(xtt[2, :]),
            period=BK.getPeriod(p.prob, x, @set par_tm[id_par] = p.p))
    end,
    plotSolution=(x, p; k...) -> begin
        xtt = BK.getPeriodicOrbit(p.prob, x, @set par_tm[id_par] = p.p)
        plot!(xtt.t, xtt[2, :]; label="V", k...)
        plot!(xtt.t, xtt[3, :]; label="N", k...)
        plot!(br; subplot=1, putspecialptlegend=false)
    end,
    normC=norminf
)

Mt = 15 # number of time sections
@time br_pocoll = BK.continuation(
    br, 4, opts_po_cont, # we want to branch form the 4th bif. point
    BK.PeriodicOrbitOCollProblem(Mt, 5), # we want to use the Collocation method to locate PO, with polynomial degree 5
    plot=true, verbosity=3;
    args_po... # regular continuation options
)
scene = plot(br, br_pocoll, markersize=3)
plot!(scene, br_pocoll.param, br_pocoll.min, label="")

# fetch the saved solutions
for sol in br_pocoll.sol
    # periodic orbit
    po = sol.x
    # get the mesh and trajectory
    traj = BK.getPeriodicOrbit(br_pocoll.prob, po, @set par_tm[id_par] = sol.p)
    plot!(scene, traj[1, :], traj[2, :], xlabel="E", ylabel="x", label="")
end
scene
