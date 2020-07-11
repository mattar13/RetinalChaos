#Version 0: Model without nullout
function (PDE::Network{T, :Default} where T)(dU, U, p, t::Float64)
    v = view(U, :, :, 1)
    n = view(U, :, :, 2)
    c = view(U, :, :, 3)
    a = view(U, :, :, 4)
    b = view(U, :, :, 5)
    e = view(U, :, :, 6)
    W = view(U, :, :, 7)

    dv = view(dU, :, :, 1)
    dn = view(dU, :, :, 2)
    dc = view(dU, :, :, 3)
    da = view(dU, :, :, 4)
    db = view(dU, :, :, 5)
    de = view(dU, :, :, 6)
    dW = view(dU, :,:,7)

    (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, σ, D, τw) = p

    @. dv = (
              fI(g_leak,  1.0, v, E_leak)
            + fI(g_Ca, M_INF(v, V1, V2), v, E_Ca)
            + fI(g_K, n , v, E_K)
            + fI(g_TREK, b, v, E_K)
            + fI(g_ACh, ħ(e, k_d), v, E_ACh)
            + I_app
            + W
        )/C_m

    @. dn = (Λ(v, V3, V4) * ((N_INF(v, V3, V4) - n)))/τn
    @. dc = (C_0 + δ*(-g_Ca*M_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc
    @. da =  (α*c^4*(1-a) - a)/τa
    @. db =  (β*a^4*(1-b) - b)/τb
    mul!(PDE.MyA, PDE.My, e)
    mul!(PDE.AMx, e, PDE.Mx)
    @. PDE.DA = D*(PDE.MyA + PDE.AMx)
    @. de = PDE.DA + (ρ*Φ(v, k, V0) - e)/τACh
    @. dW = -W/τw
    nothing
end

#Version 1: gTREK nullout
function (PDE::Network{T, :gTREK} where T)(dU, U, p, t::Float64)
    v = view(U, :, :, 1)
    n = view(U, :, :, 2)
    c = view(U, :, :, 3)
    a = view(U, :, :, 4)
    b = view(U, :, :, 5)
    e = view(U, :, :, 6)
    W = view(U, :, :, 7)

    dv = view(dU, :, :, 1)
    dn = view(dU, :, :, 2)
    dc = view(dU, :, :, 3)
    da = view(dU, :, :, 4)
    db = view(dU, :, :, 5)
    de = view(dU, :, :, 6)
    dW = view(dU, :,:,7)

    (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, σ, D, τw) = p

    @. dv = (
              fI(g_leak,  1.0, v, E_leak)
            + fI(g_Ca, M_INF(v, V1, V2), v, E_Ca)
            + fI(g_K, n , v, E_K)
            + fI(g_TREK, b, v, E_K) .* PDE.null
            + fI(g_ACh, ħ(e, k_d), v, E_ACh)
            + I_app
            + W
        )/C_m

    @. dn = (Λ(v, V3, V4) * ((N_INF(v, V3, V4) - n)))/τn
    @. dc = (C_0 + δ*(-g_Ca*M_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc
    @. da =  (α*c^4*(1-a) - a)/τa
    @. db =  (β*a^4*(1-b) - b)/τb
    mul!(PDE.MyA, PDE.My, e)
    mul!(PDE.AMx, e, PDE.Mx)
    @. PDE.DA = D*(PDE.MyA + PDE.AMx)
    @. de = PDE.DA + (ρ*Φ(v, k, V0) - e)/τACh
    @. dW = -W/τw
    nothing
end

#Version 2: gACh nullout
function (PDE::Network{T, :gACh} where T)(dU, U, p, t::Float64)
    v = view(U, :, :, 1)
    n = view(U, :, :, 2)
    c = view(U, :, :, 3)
    a = view(U, :, :, 4)
    b = view(U, :, :, 5)
    e = view(U, :, :, 6)
    W = view(U, :, :, 7)

    dv = view(dU, :, :, 1)
    dn = view(dU, :, :, 2)
    dc = view(dU, :, :, 3)
    da = view(dU, :, :, 4)
    db = view(dU, :, :, 5)
    de = view(dU, :, :, 6)
    dW = view(dU, :,:,7)

    (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, σ, D, τw) = p

    @. dv = (
              fI(g_leak,  1.0, v, E_leak)
            + fI(g_Ca, M_INF(v, V1, V2), v, E_Ca)
            + fI(g_K, n , v, E_K)
            + fI(g_TREK, b, v, E_K)
            + fI(g_ACh, ħ(e, k_d), v, E_ACh) .* PDE.null
            + I_app
            + W
        )/C_m

    @. dn = (Λ(v, V3, V4) * ((N_INF(v, V3, V4) - n)))/τn
    @. dc = (C_0 + δ*(-g_Ca*M_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc
    @. da =  (α*c^4*(1-a) - a)/τa
    @. db =  (β*a^4*(1-b) - b)/τb
    mul!(PDE.MyA, PDE.My, e)
    mul!(PDE.AMx, e, PDE.Mx)
    @. PDE.DA = D*(PDE.MyA + PDE.AMx)
    @. de = PDE.DA + (ρ*Φ(v, k, V0) - e)/τACh
    @. dW = -W/τw
    nothing
end

#Version 3 :Acetylcholine release nullout
function (PDE::Network{T, :ρ} where T)(dU, U, p, t::Float64)
    v = view(U, :, :, 1)
    n = view(U, :, :, 2)
    c = view(U, :, :, 3)
    a = view(U, :, :, 4)
    b = view(U, :, :, 5)
    e = view(U, :, :, 6)
    W = view(U, :, :, 7)

    dv = view(dU, :, :, 1)
    dn = view(dU, :, :, 2)
    dc = view(dU, :, :, 3)
    da = view(dU, :, :, 4)
    db = view(dU, :, :, 5)
    de = view(dU, :, :, 6)
    dW = view(dU, :,:,7)

    (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, σ, D, τw) = p

    @. dv = (
              fI(g_leak,  1.0, v, E_leak)
            + fI(g_Ca, M_INF(v, V1, V2), v, E_Ca)
            + fI(g_K, n , v, E_K)
            + fI(g_TREK, b, v, E_K)
            + fI(g_ACh, ħ(e, k_d), v, E_ACh)
            + I_app
            + W
        )/C_m
    @. dn = (Λ(v, V3, V4) * ((N_INF(v, V3, V4) - n)))/τn
    @. dc = (C_0 + δ*(-g_Ca*M_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc
    @. da =  (α*c^4*(1-a) - a)/τa
    @. db =  (β*a^4*(1-b) - b)/τb
    mul!(PDE.MyA, PDE.My, e)
    mul!(PDE.AMx, e, PDE.Mx)
    @. PDE.DA = D*(PDE.MyA + PDE.AMx)
    @. de = PDE.DA + (ρ*Φ(v, k, V0)*PDE.null - e)/τACh
    @. dW = -W/τw
    nothing
end

lansdell_pars = [:I_app, :E_Ca, :E_K, :E_Leak, :E_ACh, :V1, :V2, :V3, :V4, :g_Ca, :g_K, :g_Leak, :λ, :g_ACh, :δ, :C_m, :τr, :τs, :τACh, :γ, :α, :β, :k, :V0, :D, :μ]
lansdell_conds = [:v, :r, :s, :a, :W]

#This is the Lansdell version of the SAC model
function (PDE::Network{T, :Lansdell} where T)(dU, U, p, t::Float64)
    v = view(U, :, :, 1)
    r = view(U, :, :, 2)
    s = view(U, :, :, 3)
    a = view(U, :, :, 4)
    W = view(U, :, :, 5)

    dv = view(dU, :, :, 1)
    dr = view(dU, :, :, 2)
    ds = view(dU, :, :, 3)
    da = view(dU, :, :, 4)
    dW = view(dU, :, :, 5)

    (I_app, E_Ca, E_K, E_leak, E_ACh, V1, V2, V3, V4, g_Ca, g_K, g_leak, λ, g_ACh, δ, C_m, τr, τs, τACh, γ, α, β, k, V0, D, μ) = p

    @. dv = (
              fI(g_leak,  1.0, v, E_leak)
            + fI(g_Ca, M_INF(v, V1, V2), v, E_Ca)
            + fI(g_K, r , v, E_K)
            + fI(g_ACh, h(a, δ), v, E_ACh)
            + fIn(abs(W), λ, v, E_Ca)
            + I_app
        )/C_m

    @. dr = (Λ(v, V3, V4)*(N_INF(v, V3, V4) - r) + α*s*(1-r))/τr
    @. ds = γ*Φ(v, k, V0)-s/τs
    mul!(PDE.MyA, PDE.My, a)
    mul!(PDE.AMx, a, PDE.Mx)
    @. PDE.DA = D*(PDE.MyA + PDE.AMx)
    @. da = PDE.DA + β*Φ(v, k, V0)-a/τACh
    @. dW = -W
    nothing
end
