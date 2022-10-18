PDEeqs = [
     Dt(v(t, x, y)) ~ (ILeak(v(t, x, y)) + ICa(v(t, x, y)) + IK(v(t, x, y)) + ITREK(v(t, x, y)) + IACh(v(t, x, y)) + IGABA(v(t, x, y)) + INa(v(t, x, y)) + I_app + W(t, x, y)) / C_m,
     Dt(n(t, x, y)) ~ (Λ(v(t, x, y)) * ((N∞(v(t, x, y)) - n(t, x, y)))) / τn,
     Dt(m(t, x, y)) ~ α_M(v(t, x, y)) * (1 - m(t, x, y)) - β_M(v(t, x, y)) * m(t, x, y),
     Dt(h(t, x, y)) ~ α_H(v(t, x, y)) * (1 - h(t, x, y)) - β_H(v(t, x, y)) * h(t, x, y),
     Dt(c(t, x, y)) ~ (C_0 + δ * ICa(v(t, x, y)) - λ * c(t, x, y)) / τc,
     Dt(a(t, x, y)) ~ (α * c(t, x, y)^4 * (1 - a(t, x, y)) - a(t, x, y)) / τa,
     Dt(b(t, x, y)) ~ (β * a(t, x, y)^4 * (1 - b(t, x, y)) - b(t, x, y)) / τb,
     Dt(e(t, x, y)) ~ (De * ∇²(e(t, x, y)) + ρe * Φe(v(t, x, y)) - e(t, x, y)) / τACh,
     Dt(i(t, x, y)) ~ (Di * ∇²(i(t, x, y)) + ρi * Φi(v(t, x, y)) - i(t, x, y)) / τGABA,
     Dt(W(t, x, y)) ~ -W(t, x, y) / τw
]