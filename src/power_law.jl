struct PowerLaw <: ContinuousUnivariateDistribution
    α::Float64
    θ::Float64
end

PowerLaw(α::Real) = PowerLaw(α, 1.0)
PowerLaw() = PowerLaw(1.0, 1.0)
#@distr_support PowerLaw d.θ Inf
#### Parameters

shape(d::PowerLaw) = d.α
scale(d::PowerLaw) = d.θ
params(d::PowerLaw) = (d.α, d.θ)
#### Statistics

mean(d::PowerLaw) = ((α, θ) = params(d); α > 2.0 ? θ *((α - 1.0)/(α - 2.0)) : Inf)

median(d::PowerLaw) = ((α, θ) = params(d);α > 1.0 ? 2.0^(1.0/(α - 1.0)) * θ : NaN)
mode(d::PowerLaw) = d.θ

function var(d::PowerLaw)
    (α, θ) = params(d)
    α > 3.0 ? (θ^2 * (α-1)) / ((α - 2.0)^2 * (α - 3.0)) : Inf
end

function skewness(d::PowerLaw)
    α = shape(d)
    α > 4.0 ? ((2.0 * (α)) / (α - 4.0)) * sqrt((α - 3.0) / (α-1)) : NaN
end

function kurtosis(d::PowerLaw)
    α = shape(d)
    α > 5.0 ? (6.0 * ((α-1)^3 + (α-1)^2 - 6.0 * (α-1) - 2.0)) / ((α-1) * (α - 4.0) * (α - 5.0)) : NaN
end

entropy(d::PowerLaw) = ((α, θ) = params(d); log(θ / (α-1)) + 1.0 / (α-1) + 1.0)

import Plots.plot
plot(d::PowerLaw, rng) = plot(rng, map(x -> pdf(d, x), rng))

#### Evaluation
export pdf, logpdf, cdf, ccdf, logcdf, logccdf

function pdf(d::PowerLaw, x::Float64)
    (α, θ) = params(d)
    x >= θ ? ((α-1.0)/θ) * ((x/θ)^(-α)) : 0.0
end

function pdf(d::PowerLaw, x::AbstractArray)
    (α, θ) = params(d)
    cons = ((α-1.0)/θ)
    pdfs = Array(Float64,0)
    for num in x
        push!(pdfs,(num >= θ ? cons * ((x/θ)^(-α)) : 0.0))
    end
    return pdfs
end

function logpdf(d::PowerLaw, x::Float64)
    (α, θ) = params(d)
    x >= θ ? log(α-1.0) - log(θ) - α * log(x/θ) : -Inf
end

function logpdf(d::PowerLaw, x::AbstractArray)
    (α, θ) = params(d)
    l_const = log(α-1.0) - log(θ)
    lpdfs = Array(Float64,0)
    for num in x
        push!(lpdfs,(num >= θ ? l_const - α * log(num/θ) : -Inf))
    end
    return lpdfs
end


function ccdf(d::PowerLaw, x::Float64)
    (α, θ) = params(d)
    x >= θ ? (x/θ)^(-α +1.0) : 1.0
end

cdf(d::PowerLaw, x::Float64) = 1.0 - ccdf(d, x)

logccdf(d::PowerLaw, x::Float64) = log(ccfd(d,x))

logcdf(d::PowerLaw, x::Float64) = log(cdf(d, x))

cquantile(d::PowerLaw, p::Float64) = d.θ *((p) ^ (-1.0 / (d.α - 1.0)))
quantile(d::PowerLaw, p::Float64) = cquantile(d, 1.0 - p)


#### Sampling
import Base.rand
rand(d::PowerLaw) = quantile(d,rand())
rand(d::PowerLaw, n::Int64) = quantile.(d, rand(n::Int64))
## Fitting

import StatsBase.fit
function fit(::Type{PowerLaw}, x::AbstractArray{T} where T<:Real)
    θ = minimum(x)

    n = length(x)
    lθ = log(θ)
    temp1 = 0.0
    for i=1:n
        temp1 += log(x[i]) - lθ
    end
    α = 1.0+n*(temp1)^(-1.0)

    return PowerLaw(α, θ)
end
