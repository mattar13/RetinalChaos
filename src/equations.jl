println("Loading model")
@named SAC_seperate = ODESystem([
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
ODEModel = structural_simplify(SAC_seperate)

#noise_eqs = zeros(length(ODEModel.eqs))
noise_eqs = zeros(length(ODEModel.eqs))
noise_eqs[end] = 0.1
#Add in the observed variables

@named SDEModel = SDESystem(ODEModel, noise_eqs)