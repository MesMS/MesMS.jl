import ProgressMeter: @showprogress

aa_elements = (  # H,  C, N, O, S
    (sym=:Alg, n=[ 5,  3, 1, 1, 0], w=0.0876),
    (sym=:Arg, n=[12,  6, 4, 1, 0], w=0.0578),
    (sym=:Asn, n=[ 6,  4, 2, 2, 0], w=0.0393),
    (sym=:Asp, n=[ 5,  4, 1, 3, 0], w=0.0549),
    (sym=:Cys, n=[ 5,  3, 1, 1, 1], w=0.0138),
    (sym=:Gln, n=[ 8,  5, 2, 2, 0], w=0.039),
    (sym=:Glu, n=[ 7,  5, 1, 3, 0], w=0.0632),
    (sym=:Gly, n=[ 3,  2, 1, 1, 0], w=0.0703),
    (sym=:His, n=[ 7,  6, 3, 1, 0], w=0.0226),
    (sym=:Ile, n=[11,  6, 1, 1, 0], w=0.0549),
    (sym=:Leu, n=[11,  6, 1, 1, 0], w=0.0968),
    (sym=:Lys, n=[12,  6, 2, 1, 0], w=0.0519),
    (sym=:Met, n=[ 9,  5, 1, 1, 1], w=0.0232),
    (sym=:Phe, n=[ 9,  9, 1, 1, 0], w=0.0387),
    (sym=:Pro, n=[ 7,  5, 1, 1, 0], w=0.0502),
    (sym=:Ser, n=[ 5,  3, 1, 2, 0], w=0.0714),
    (sym=:Thr, n=[ 7,  4, 1, 2, 0], w=0.0553),
    (sym=:Trp, n=[10, 11, 2, 1, 0], w=0.0125),
    (sym=:Tyr, n=[ 9,  9, 1, 2, 0], w=0.0291),
    (sym=:Val, n=[ 9,  5, 1, 1, 0], w=0.0673),
)

Iso_H = (w=[0.99985, 0.00015], m=[1.007825, 2.0140])
Iso_C = (w=[0.9889, 0.0111], m=[12., 13.00335])
Iso_N = (w=[0.9964, 0.0036], m=[14.00307, 15.00011])
Iso_O = (w=[0.9976, 0.0004, 0.002], m=[15.99491, 16.99913, 17.99916])
Iso_S = (w=[0.950, 0.0076, 0.0422], m=[31.97207, 32.97146, 33.96786])

"""
XW, YW: vectors of weights
XV, YV: vectors of values
"""
conv_ipv((XW, XV), (YW, YV)) = begin
    WX = map(1:(length(XW) + length(YW) - 1)) do i
        l, r = max(1, i - length(YW) + 1), min(i, length(XW))
        ws = XW[l:r] .* YW[i+1-l:-1:i+1-r]
        vs = XV[l:r] .+ YV[i+1-l:-1:i+1-r]
        v, w = sum(vs .* ws), sum(ws)
        return w, (w != 0.0 ? v / w : 0.0)
    end
    return first.(WX), last.(WX)
end

read_ipv(io) = begin
    map(1:2) do _
        n = read(io, Int64)
        return map(1:n) do _
            m = read(io, Int64)
            return read!(io, Array{Float64}(undef, m))
        end
    end |> Tuple
end

write_ipv(io, (W, V)) = begin
    for X in (W, V)
        write(io, Int64(length(X)))
        for xs in X
            write(io, Int64(length(xs)), Float64.(xs))
        end
    end
end

calc_ipv(ms, E, n_mean, trunc=0.99) = begin
    m_mono = map(e -> first(e.m), E)
    n_mean = n_mean ./ sum(n_mean .* m_mono)
    E = map(e -> [(e.w ./ sum(e.w), e.m .- first(e.m))], E)
    V = @showprogress map(ms) do m
        ns = round.(Int, n_mean .* m)
        ns[begin] += round(Int, (m - sum(ns .* m_mono)) / m_mono[begin]) # fill the rest part with `H`
        for (n, e) in zip(ns, E)
            while length(e) < n push!(e, conv_ipv(e[begin], e[end])) end
        end
        elements = [e[n] for (n, e) in zip(ns, E) if n > 0]
        v = foldl(conv_ipv, sort(elements, by=length∘first))
        v = (v[1] ./ sum(v[1]), v[2])
        s = 0.0
        i = 0
        while s < trunc
            i += 1
            s += v[1][i]
        end
        return v[1][1:i], v[2][1:i]
    end
    return first.(V), last.(V)
end

build_ipv(path=joinpath(homedir(), ".MesMS", "default.ipv"), r=1:20000, trunc=0.99; verbose=true) = begin
    if isfile(path)
        verbose && @info "IPV loading from " * path
        V = open(read_ipv, path)
    else
        verbose && @info "IPV building"
        E = [Iso_H, Iso_C, Iso_N, Iso_O, Iso_S]
        n_mean = sum(a -> a.n .* a.w, aa_elements)
        V = calc_ipv(r, E, n_mean, trunc)
        safe_save(p -> open(io -> write_ipv(io, V), p; write=true), path, "IPV")
    end
    return V
end

ipv_w(m::Number, V) = V[1][round(Int, m)]
ipv_w(ion, V) = ipv_w(ion.mz * ion.z, V)
ipv_m(m::Number, V) = V[2][round(Int, m)]
ipv_m(ion, V) = ipv_m(ion.mz * ion.z, V)
ipv_mz(mz::Number, z, n, V) = mz + V[2][round(Int, mz * z)][n] / z
ipv_mz(mz::Number, z, V) = mz .+ V[2][round(Int, mz * z)] ./ z
ipv_mz(ion, n, V) = ipv_mz(ion.mz, ion.z, n, V)
ipv_mz(ion, V) = ipv_mz(ion.mz, ion.z, V)
