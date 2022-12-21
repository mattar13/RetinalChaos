println("Loading dimensions")
#These parameters are about the x, y, and time dimension
@parameters t x y #Dimensions

println("Loading parameters")
@variables v(t) n(t) m(t) h(t) c(t) a(t) b(t) e(t) i(t) W(t) #Initial conditions
@variables v̂(..) n̂(..) m̂(..) ĥ(..) ĉ(..) â(..) b̂(..) ê(..) î(..) Ŵ(..)
#Observables
@variables I_Ca(t) I_Na(t) I_K(t)
@variables I_ext(t)
@parameters g_leak E_leak g_Ca V1 V2 E_Ca g_K E_K g_TREK g_ACh k_ACh E_ACh g_GABA k_GABA E_Cl I_app C_m
@parameters V3 V4 τn
@parameters C_0 λ δ τc
@parameters α τa β τb ρe ρi τACh τGABA VSe VSi V0e V0i
@parameters De Di
@parameters τw σ
@parameters g_Na E_Na
@parameters V7 V8 V9
@parameters V10 V11 V12
@parameters V13 V14 V15
@parameters V16 V17 V18

Dt = Differential(t)
Dx = Differential(x)
Dy = Differential(y)
Dxx = Differential(x)^2
Dyy = Differential(y)^2

parameters = Dict(
     g_leak => 2.0,
     E_leak => -70.0,
     g_Ca => 8.5,
     V1 => -20.0,
     V2 => 20.0,
     E_Ca => 50.0,
     g_K => 4.0,
     E_K => -90.0,
     g_TREK => 3.0,
     I_app => 0.0,
     C_m => 13.6,
     V3 => -25.0,
     V4 => 7.0,
     τn => 5.0,
     C_0 => 0.088,
     λ => 2.702,
     δ => 0.010503,
     τc => 2000.0,
     α => 625.0,
     τa => 8300.0,
     β => 34.0,
     τb => 8300.0,
     ρe => 6.0,
     ρi => 5.0,
     τACh => 540.0,
     τGABA => 1000.0,
     VSe => 0.2,
     VSi => 0.2,
     V0e => -40.0,
     V0i => -40.0,
     g_ACh => 0.215,
     k_ACh => 0.1,
     E_ACh => 0.0,
     g_GABA => 0.9,
     k_GABA => 0.1,
     E_Cl => -65.0,
     g_Na => 2.0,
     E_Na => 55.0, #55.0
     De => 0.01,
     Di => 0.01,
     τw => 800.0,
     V7 => 10.0,
     V8 => -40.0,
     V9 => 10.0,
     V10 => 4.0,
     V11 => -65.0,
     V12 => 18.0,
     V13 => 0.07,
     V14 => -65.0,
     V15 => 20.0,
     V16 => 1.0,
     V17 => -35.0,
     V18 => 10.0,
     σ => 0.1,
)

u0 = Dict(
     I_ext => 0.0,
     v => -63.6,
     I_Na => 0.0,
     I_Ca => 0.0,
     I_K => 0.0,
     n => 0.0,
     m => 0.062,
     h => 0.550,
     c => 0.085,
     a => 0.026,
     b => 0.0,
     e => 0.066,
     i => 0.053,
     W => 0.0
)

#Load the time and space domains
#if the parameters have been reset you can use this to reload their defaults
function reload_parameters()
     parameters[g_leak] = 2.0
     parameters[E_leak] = -70.0
     parameters[g_Ca] = 8.5
     parameters[V1] = -20.0
     parameters[V2] = 20.0
     parameters[E_Ca] = 50.0
     parameters[g_K] = 4.0
     parameters[E_K] = -90.0
     parameters[g_TREK] = 3.0
     parameters[I_app] = 0.0
     parameters[C_m] = 13.6
     parameters[V3] = -25.0
     parameters[V4] = 7.0
     parameters[τn] = 5.0
     parameters[C_0] = 0.088
     parameters[λ] = 2.702
     parameters[δ] = 0.010503
     parameters[τc] = 2000.0
     parameters[α] = 625.0
     parameters[τa] = 8300.0
     parameters[β] = 34.0
     parameters[τb] = 8300.0
     parameters[ρe] = 6.0
     parameters[ρi] = 5.0
     parameters[τACh] = 540.0
     parameters[τGABA] = 1000.0
     parameters[VSe] = 0.2
     parameters[VSi] = 0.2
     parameters[V0e] = -40.0
     parameters[V0i] = -40.0
     parameters[g_ACh] = 0.215
     parameters[k_ACh] = 0.1
     parameters[E_ACh] = 0.0
     parameters[g_GABA] = 0.9
     parameters[k_GABA] = 0.1
     parameters[E_Cl] = -65.0
     parameters[g_Na] = 2.0
     parameters[E_Na] = 55.0 #55.
     parameters[De] = 0.005
     parameters[Di] = 0.005
     parameters[τw] = 800.0
     parameters[σ] = 0.1
     parameters[V7] = 10.0
     parameters[V8] = -40.0
     parameters[V9] = 10.0
     parameters[V10] = 4.0
     parameters[V11] = -65.0
     parameters[V12] = 18.0
     parameters[V13] = 0.07
     parameters[V14] = -65.0
     parameters[V15] = 20.0
     parameters[V16] = 1.0
     parameters[V17] = -35.0
     parameters[V18] = 10.0
end