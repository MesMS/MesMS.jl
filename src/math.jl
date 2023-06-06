using Statistics

error_rel(a, b) = (b - a) / max(abs(a), abs(b))
error_ppm(a, b) = error_rel(a, b) * 1e6
error_ppb(a, b) = error_rel(a, b) * 1e9

# moe: margin of error
in_moe(a, b, ε) = abs(a - b) <= max(abs(a), abs(b)) * ε

m_to_mz(m, z) = (m + mₚ * z) / abs(z)
mz_to_m(mz, z) = mz * abs(z) - mₚ * z

mh_to_mz(mh, z) = begin
    if z > 0 return (mh + mₚ * (z - 1)) / z end
    if z < 0 return (mh - mₚ * (-z - 1)) / -z end
    return 0.0
end

mz_to_mh(mz, z) = begin
    if z > 0 return (mz - mₚ) * z + mₚ end
    if z < 0 return (mz + mₚ) * -z - mₚ end
    return 0.0
end

calc_centroid(xs, ws) = begin
    # TODO: the `@assert` may be skipped at runtime
    @assert length(xs) == length(ws) "`xs` and `ws` are not of the same length"
    @assert !isempty(xs) "`xs` and `ws` are both empty"
    if xs[begin] == xs[end]
        return xs[begin], mean(ws)
    else
        ss2 = (ws[begin:end-1] .+ ws[begin+1:end]) .* (xs[begin+1:end] .- xs[begin:end-1])
        s2 = sum(ss2)
        x2 = sum(ss2 .* (xs[begin:end-1] .+ xs[begin+1:end])) / 2
        return x2 / s2, s2 / 2 / (xs[end] - xs[begin])
    end
end

log_softer(s=1) = x -> copysign(s * (log(abs(x) + s) - log(s)), x)
exp_softer(s=1) = x -> x * exp(-abs(x)/s)

"linear interpolation between `a` and `b` with distance `δa` (and `δb`)"
lerp(a, b, δa) = a + (b - a) * δa
lerp(a, b, δa, δb) = lerp(a, b, δa / (δa + δb)) # == (δb * a + δa * b) / (δa + δb)

calc_tda_fdr_xl(xs, r=1.0) = begin
    tt, td, dd = 0, 0, 0
    fdr = zeros(length(xs))
    for (i, x) in enumerate(xs)
        if x == :TT tt += 1
        elseif x == :TD td += 1
        elseif x == :DD dd += 1
        else error("unexcepted T/D type: $(x)")
        end
        fdr[i] = max(dd / r^2, (td / r - dd / r^2)) / tt
    end
    v = 1.0
    for i in reverse(eachindex(fdr))
        v = min(v, fdr[i])
        fdr[i] = v
    end
    return fdr
end

include("esimate_td_xl.jl")
