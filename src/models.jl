mutable struct BurstPDE{A} <: Function
        Mx::Union{Tridiagonal, A}
        My::Union{Tridiagonal, A}
        MyA::A 
        AMx::A 
        DA::A
        null::A
        nullout::Symbol
end


#Auxillary functions
M_INF(V, V1, V2) = (1 + tanh((V - V1)/V2))/2;
N_INF(V, V3, V4) = (1 + tanh((V - V3)/V4))/2;
LAM(V, V3, V4) = cosh((V-V3)/(2*V4));
Φ(v, κ, V_0) = 1/(1 + exp(-κ * (v - V_0)))
ħ(a, K_d) = (a^2)/(a^2 + K_d)
#This is the generalized current calculation function
fI(g::Float64, r, v, e::Float64) = -g*r*(v-e)
#This file will include several models
#1) Starburst Amacrine Cell Burst models

BurstModel = @ode_def begin
    dv = (
          fI(g_leak,  1.0, v, E_leak)
        + fI(g_Ca, M_INF(v, V1, V2), v, E_Ca)
        + fI(g_K, n , v, E_K)
        + fI(g_TREK, b, v, E_K)
        + fI(g_ACh, ħ(ACh, k_d), v, E_ACh)
        + I_app
        + W
        )/C_m
    dn = (LAM(v, V3, V4) * ((N_INF(v, V3, V4) - n)))/τn
    dc = (C_0 + δ*(-g_Ca * M_INF(v, V1, V2)* (v - E_Ca)) - λ*c)/τc
    da =  (α*c^4*(1-a) - a)/τa
    db =  (β*a^4*(1-b) - b)/τb
    dACh = (-2*D*ACh) + (ρ*Φ(v, k, V0) - ACh)/τACh
    dW = -W/τw
end g_leak E_leak g_Ca V1 V2 E_Ca g_K E_K g_TREK g_ACh k_d E_ACh I_app C_m V3 V4 τn C_0 λ δ τc α τa β τb ρ τACh k V0 σ D τw;
model_pars = BurstModel.params
model_conds = BurstModel.syms

function diffuse(lattice, D, PDE::BurstPDE)
    mul!(PDE.MyA, PDE.My, lattice)
    mul!(PDE.AMx, lattice, PDE.Mx)
    DA = D*(PDE.MyA + PDE.AMx)
    return lattice + DA
end

function (PDE::BurstPDE)(dU, U, p, t)
    v = view(U, :, :, 1)
    n = view(U, :, :, 2)
    c = view(U, :, :, 3)
    a = view(U, :, :, 4)
    b = view(U, :, :, 5)
    ACh = view(U, :, :, 6)
    W = view(U, :, :, 7)

    dv = view(dU, :, :, 1)
    dn = view(dU, :, :, 2)
    dc = view(dU, :, :, 3)
    da = view(dU, :, :, 4)
    db = view(dU, :, :, 5)
    dACh = view(dU, :, :, 6)
    dW = view(dU, :,:,7)

    (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, σ, D, τw) = p
    if PDE.nullout == :g_TREK
        @. dv = (
                  fI(g_leak,  1.0, v, E_leak)
                + fI(g_Ca, M_INF(v, V1, V2), v, E_Ca)
                + fI(g_K, n , v, E_K)
                + fI(g_TREK, b, v, E_K) .* PDE.null
                + fI(g_ACh, ħ(ACh, k_d), v, E_ACh) 
                + I_app
                + W
            )/C_m
    elseif PDE.nullout == :g_ACh
        @. dv = (
                  fI(g_leak,  1.0, v, E_leak)
                + fI(g_Ca, M_INF(v, V1, V2), v, E_Ca)
                + fI(g_K, n , v, E_K)
                + fI(g_TREK, b, v, E_K) 
                + fI(g_ACh, ħ(ACh, k_d), v, E_ACh) .* PDE.null
                + I_app
                + W
            )/C_m
    else
        @. dv = (
                  fI(g_leak,  1.0, v, E_leak)
                + fI(g_Ca, M_INF(v, V1, V2), v, E_Ca)
                + fI(g_K, n , v, E_K)
                + fI(g_TREK, b, v, E_K)
                + fI(g_ACh, ħ(ACh, k_d), v, E_ACh)
                + I_app
                + W
            )/C_m
    end

    @. dn = (LAM(v, V3, V4) * ((N_INF(v, V3, V4) - n)))/τn
    @. dc = (C_0 + δ*(-g_Ca*M_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc
    @. da =  (α*c^4*(1-a) - a)/τa
    @. db =  (β*a^4*(1-b) - b)/τb
    mul!(PDE.MyA, PDE.My, ACh)
    mul!(PDE.AMx, ACh, PDE.Mx)
    @. PDE.DA = D*(PDE.MyA + PDE.AMx)
    if PDE.nullout == :ρ
        @. dACh = PDE.DA + (ρ*Φ(v, k, V0)*PDE.null - ACh)/τACh
    else
        @. dACh = PDE.DA + (ρ*Φ(v, k, V0) - ACh)/τACh
    end
    @. dW = -W/τw
    nothing
end


"""
This constructs the PDE function so that it can be called
"""
function BurstPDE(nx::Int64, ny::Int64; gpu::Bool = false, μ::Float64 = 0.75, nullout = :g_ACh,
        DX::Tuple{Float64, Float64} = (-2.0, 1.0), DY::Tuple{Float64, Float64} = (-2.0, 1.0))
    #Set up x diffusion steps
    x_dv = repeat([DX[1]], nx)
    x_uv = repeat([DX[2]], nx-2)
    x_lv = repeat([DX[2]], nx-2)
    push!(x_uv, abs(DX[1]))
    pushfirst!(x_lv, abs(DX[1]))

    #Set up y diffusion steps
    y_dv = repeat([DY[1]], ny)
    y_uv = repeat([DY[2]], ny-2)
    y_lv = repeat([DY[2]], ny-2)
    push!(y_lv, abs(DY[1]))
    pushfirst!(y_uv, abs(DY[1]))
    Mx = Tridiagonal(x_lv, x_dv, x_uv)
    My = Tridiagonal(y_lv, y_dv, y_uv)
    MyA = zeros(ny, nx);
    AMx = zeros(ny, nx);
    DA = zeros(ny, nx);
    d = Binomial(1, μ)
    null = rand(d, (ny, nx))
    if gpu
        gMx = CuArray(Float32.(Mx))
        gMy = CuArray(Float32.(Mx))
        gAMx = CuArray(Float32.(Mx))
        gMyA = CuArray(Float32.(Mx))
        gDA = CuArray(Float32.(Mx))
        gNull = CuArray(Float32.(null)
        return BurstPDE(gMx, gMy, gMyA, gAMx, gDA, null, nullout)
    else
        return BurstPDE(Mx, My, Array(MyA), Array(AMx), Array(DA), null, nullout)
    end
end

#Noise models
noise(du, u, p, t) = du[7] = p[30]
noise_2D(du, u, p, t) = du[:,:,7] .= p[30]

p_find(p) = findall(x -> x == p, BurstModel.params)
u_find(u) = findall(x -> x == u, BurstModel.syms)
