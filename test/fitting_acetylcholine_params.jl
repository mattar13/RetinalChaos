using Revise
using RetinalChaos
import RetinalChaos.ħ

#%% Take info from Zhous 2004 paper
#extract the function for ACh activation
p = read_JSON(params_file)
#x -> voltage
#p -> [gACh, k_d]
model1(x, par) = map(v -> par[1] * ħ(1.0, par[2]) * (v - p[:E_ACh]), x)
v_known = [35.0, 15.0, -5.0, -25.0, -45.0, -75.0] #Zhou known v's
i_known = [0.0, 0.0, -25.0, -50.0, -100.0, -200.0] #Zhou known I's
p0 = [0.75, 0.1]
fit = curve_fit(model, v_known, i_known, p0, lower=[0.0, 0.0], upper=[Inf, Inf])
p1 = plot(v -> model1(v, fit.param), -75, 35)
p2 = plot(x -> ħ(x, fit.param[2]), 0.001, 0.1)
plot(p1, p2, layout=grid(2, 1))
