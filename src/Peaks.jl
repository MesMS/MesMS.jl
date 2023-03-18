abstract type AbstractPeak end

Base.isless(a::AbstractPeak, b::AbstractPeak) = a.mz < b.mz
Base.isless(peak::AbstractPeak, mz::Real) = peak.mz < mz
Base.isless(mz::Real, peak::AbstractPeak) = mz < peak.mz

struct Peak <: AbstractPeak
    mz::Float64
    inten::Float64
end

merge_peaks(peaks, ε) = begin
    ans = empty(peaks)
    if isempty(peaks)
        return ans
    end
    mz, inten, weight = peaks[begin].mz, peaks[begin].inten, peaks[begin].inten
    for p in peaks[begin+1:end]
        if in_moe(p.mz, mz, ε)
            mz = (mz * weight + p.mz * p.inten) / (weight + p.inten)
            inten = max(inten, p.inten)
            weight += p.inten
        else
            push!(ans, Peak(mz, inten))
            mz, inten, weight = p.mz, p.inten, p.inten
        end
    end
    push!(ans, Peak(mz, inten))
    return ans
end

max_inten(ps, lower, upper) = maximum(p -> p.inten, query(ps, lower, upper); init=0.0)
max_inten_δ(ps, x, δ) = max_inten(ps, x - δ, x + δ)
max_inten_ε(ps, x, ε) = max_inten_δ(ps, x, ε * x)

pick_by_inten(ps, n) = begin
    if length(ps) <= n return ps end
    τ = partialsort!(map(p -> p.inten, ps), n; rev=true)
    return filter(p -> p.inten ≥ τ, ps)
end
