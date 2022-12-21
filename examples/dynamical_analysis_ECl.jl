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
scene = plot()
for x = LinRange(1.0, 10.0, 20)
    reload_parameters()
    parameters[I_app] = -10.0
    parameters[g_ACh] = x
    parameters[g_TREK] = 0.0 #Remove the sAHP
    u0[i] = parameters[ρi]
    u0[e] = parameters[ρe]

    @time ODEProb = ODEProblem(ODEModel, u0, (tmin, tmax), parameters, jac=true)
    sol = solve(ODEProb)
    #initial_plot = plot(sol, idxs=[v, e, i], layout = 3)

    ODEFunc = ODEProb.f
    F = (u, p) -> ODEFunc(u, p, 0.0) #Make the ODE function
    J = (u, p) -> ODEFunc.jac(u, p, 0.0) #Make the jacobian

    id_par = indexof(I_app)
    par_tm = ODEProb.p

    prob = BK.BifurcationProblem(F, ODEProb.u0, par_tm, (@lens _[id_par]); #Setup the BifurcationProblem 
        J=J,#Setup the Jacobian
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
    plot!(scene, br, label = "")
    #=
    # This section determines the periodic orbits 
    optn_po = BK.NewtonPar(verbose=true, tol=1e-8, maxIter=100)
    opts_po_cont = BK.ContinuationPar(
        dsmin=0.04, ds=0.04, dsmax=0.5,
        pMin=0.0, pMax=100.0,
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
    plot!(scene, br_pocoll, markersize=3)
    plot!(scene, br_pocoll.param, br_pocoll.min, label="")

    # fetch the saved solutions
    for sol in br_pocoll.sol
        # periodic orbit
        po = sol.x
        # get the mesh and trajectory
        traj = BK.getPeriodicOrbit(br_pocoll.prob, po, @set par_tm[id_par] = sol.p)
        plot!(scene, traj[1, :], traj[2, :], xlabel="E", ylabel="x", label="")
    end
    =#
end
scene
#%% Continuation of fold point
opts_br_fold = BK.ContinuationPar(opts_br, pMin = -75.0, pMax = -55.0,  dsmin = 0.001, ds = 0.01, dsmax = 0.5)
codim2_idx = indexof(E_Cl)

sn_codim2 = BK.continuation(br, 2, (@lens _[codim2_idx]), 
    normC=norminf, # define the sup norm
    opts_br_fold, #fold branching options
    recordFromSolution=(x, p) -> (V=x[2]),    
    detectCodim2Bifurcation = 2, # detection of codim 2 bifurcations with bisection
	updateMinAugEveryStep = 1, # we update the Fold problem at every continuation step
    bothside = true, plot=true, 
    verbosity=3,
)
sn_codim2
plot!(scene, sn_codim2, vars = (:V, E_Cl), branchlabel = "Fold")