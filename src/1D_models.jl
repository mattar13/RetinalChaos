######################Modelling Toolkit for Tarchick Model
@parameters t 
#Register parameters for ionic currents 
@parameters g_leak E_leak #Leaky currents
@parameters g_Ca V1 V2 E_Ca #Cav currents
@parameters g_K V3 V4 E_K #Kv currents
@parameters g_TREK g_sAHP #TREK and sAHP currents
@parameters g_ACh k_d E_ACh #Acetylcholine current parameters
@parameters g_HCN V5 V6 E_HCN
@parameters I_app C_m  #Extras
@parameters τn τc τa τb τACh τw τr τs #Time constants
@parameters C_0 V0
@parameters λ δ α β ρ k σ D #Relationships in my model
@parameters δc αs αc αr H_x VS VH g_n R E_n #Relationships in other models
@variables v(t) n(t) c(t) a(t) b(t) e(t) W(t) r(t) s(t)
@derivatives d'~t

T_model_eqs = [
          d(v)~ (
                    -g_leak*(v-E_leak)
                  + -g_Ca*M_INF(v, V1, V2)*(v-E_Ca)
                  + -g_K*n*(v-E_K)
                  + -g_TREK*b*(v-E_K)
                  + -g_ACh*ħ(e, k_d)*(v-E_ACh)
                  + -g_HCN*H_INF(v, V5, V6)*(v-E_HCN) #Add the HCN channels here
                  + I_app
                  + W + 0*σ #This 
                  )/C_m ,
          d(n) ~ (Λ(v, V3, V4) * (N_INF(v, V3, V4) - n))/τn,
          d(c) ~ (C_0 + δ*(-g_Ca*M_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc,
          d(a) ~ (α*c^4*(1-a)-a)/τa,
          d(b) ~ (β*a^4*(1-b)-b)/τb,
          d(e) ~ (-2*D*e) + (ρ*Φ(v, k, V0) - e)/τACh,
          d(W) ~ -W/τw
]
T_noise_eqs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, σ]
T_ode = ODESystem(T_model_eqs)
T_sde = SDESystem(T_model_eqs, T_noise_eqs, t, T_ode.states, T_ode.ps)
sym_ps = map(x -> x |> Symbol, T_ode.ps)
sym_cs = map(x -> x |> Symbol, T_ode.states)

#Karvourniari Model
K_model_eqs = [
          d(v) ~ (-g_leak*(v-E_leak)
                    + -g_Ca*M_INF(v, V1, V2)*(v-E_Ca)
                    + -g_K*n*(v-E_K)
                    + -g_sAHP*r^4*(v-E_K)
                    + I_app
                    + W + 0*σ
                    )/C_m,
          d(n)~ (Λ(v, V3, V4) * ((N_INF(v, V3, V4) - n)))/τn,
          d(r) ~ (αr*s*(1-r)-r)/τr,
          d(s) ~ (αs*c^4*(1-s)-s)/τs,
          d(c) ~ (-(αc/H_x)*c + C_0 + δc*(-g_Ca*M_INF(v, V1, V2)*(v-E_Ca)))/τc,
          d(W) ~ -W
]
K_noise_eqs = [0.0, 0.0, 0.0, 0.0, 0.0, σ]
K_ode = ODESystem(K_model_eqs)
K_sde = SDESystem(K_model_eqs, K_noise_eqs, t, K_ode.states, K_ode.ps)
