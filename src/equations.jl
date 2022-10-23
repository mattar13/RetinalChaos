print("Loading ODE models... ")
@named SAC_Algebraic = ODESystem([
     Dt(I_ext) ~ -I_ext + I_app, #This parameter is controlled by an outside      
     Dt(v) ~ (ILeak(v) + ICa(v) + IK(v) + ITREK(v) + IACh(v) + IGABA(v) + INa(v) + I_ext + W) / C_m,
     I_Ca ~ ICa(v),
     I_Na ~ INa(v),
     I_K ~ IK(v),
     Dt(n) ~ (Λ(v) * ((N∞(v) - n))) / τn,
     Dt(m) ~ α_M(v) * (1 - m) - β_M(v) * m,
     Dt(h) ~ α_H(v) * (1 - h) - β_H(v) * h,
     Dt(c) ~ (C_0 + δ * (ICa(v)) - λ * c) / τc,
     Dt(a) ~ (α * c^4 * (1 - a) - a) / τa,
     Dt(b) ~ (β * a^4 * (1 - b) - b) / τb,
     Dt(e) ~ (ρe * Φe(v) - e) / τACh,
     Dt(i) ~ (ρi * Φi(v) - i) / τGABA,
     Dt(W) ~ -W / τw
])
ODEModel = structural_simplify(SAC_Algebraic)
println("Complete")

print("Loading SDE models... ")
#noise_eqs = zeros(length(ODEModel.eqs))
noise_eqs = zeros(length(ODEModel.eqs))
noise_eqs[end] = 0.1
#Add in the observed variables
@named SDEModel = SDESystem(ODEModel, noise_eqs)
println("Complete")

println("Loading PDE models")
PDEeqs = [
     Dt(v̂(x, y, t)) ~ (ÎLeak(v̂(x, y, t)) + ÎCa(v̂(x, y, t)) + ÎK(v̂(x, y, t)) + ÎTREK(v̂(x, y, t)) + ÎACh(v̂(x, y, t)) + ÎGABA(v̂(x, y, t)) + ÎNa(v̂(x, y, t)) + I_app + Ŵ(x, y, t)) / C_m,
     Dt(n̂(x, y, t)) ~ (Λ(v̂(x, y, t)) * ((N∞(v̂(x, y, t)) - n̂(x, y, t)))) / τn,
     Dt(m̂(x, y, t)) ~ α_M(v̂(x, y, t)) * (1 - m̂(x, y, t)) - β_M(v̂(x, y, t)) * m̂(x, y, t),
     Dt(ĥ(x, y, t)) ~ α_H(v̂(x, y, t)) * (1 - ĥ(x, y, t)) - β_H(v̂(x, y, t)) * ĥ(x, y, t),
     Dt(ĉ(x, y, t)) ~ (C_0 + δ * ÎCa(v̂(x, y, t)) - λ * ĉ(x, y, t)) / τc,
     Dt(â(x, y, t)) ~ (α * ĉ(x, y, t)^4 * (1 - â(x, y, t)) - â(x, y, t)) / τa,
     Dt(b̂(x, y, t)) ~ (β * â(x, y, t)^4 * (1 - b̂(x, y, t)) - b̂(x, y, t)) / τb,
     Dt(ê(x, y, t)) ~ (De * ∇²(ê(x, y, t)) + ρe * Φe(v̂(x, y, t)) - ê(x, y, t)) / τACh,
     Dt(î(x, y, t)) ~ (Di * ∇²(î(x, y, t)) + ρi * Φi(v̂(x, y, t)) - î(x, y, t)) / τGABA,
     Dt(Ŵ(x, y, t)) ~ -Ŵ(x, y, t) / τw
]

#Load the PDE system
@named PDEModel = PDESystem(PDEeqs, bcs, domains, dimensions, states, ps) #Create the undiscretized PDE system

println("Discretizing the model")
#Create the Discretization scheme
discretization = MOLFiniteDifference([x => dx, y => dy], t)
grid = get_discrete(PDEModel, discretization) #Make a representation of the discrete map
@time probPDE = discretize(PDEModel, discretization) # This gives an ODEProblem since it's time-dependent

#Convert the PDE model into 