import BSON
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

Iso_H = [(m=1.007825, w=0.99985), (m=2.0140, w=0.00015)]
Iso_C = [(m=12., w=0.9889), (m=13.00335, w=0.0111)]
Iso_N = [(m=14.00307, w=0.9964), (m=15.00011, w=0.0036)]
Iso_O = [(m=15.99491, w=0.9976), (m=16.99913, w=0.0004), (m=17.99916, w=0.002)]
Iso_S = [(m=31.97207, w=0.950), (m=32.97146, w=0.0076), (m=33.96786, w=0.0422)]

conv_ipv(a, b) = begin
    map(1:(length(a) + length(b) - 1)) do i
        l, r = max(1, i - length(b) + 1), min(i, length(a))
        ms = [x.m + y.m for (x, y) in zip(a[l:r], b[i+1-l:-1:i+1-r])]
        ws = [x.w * y.w for (x, y) in zip(a[l:r], b[i+1-l:-1:i+1-r])]
        m, w = sum(ms .* ws), sum(ws)
        return (m=(w != 0.0 ? m / w : 0.0), w=w)
    end
end

calc_ipv(ms, E, n_mean, trunc=0.99) = begin
    m_mono = map(e -> e[begin].m, E)
    n_mean = n_mean ./ sum(n_mean .* m_mono)
    E = map(e -> map(i -> (m=i.m, w=i.w/sum(x -> x.w, e)), e), E)
    E = map(e -> map(i -> (m=i.m - e[begin].m, w=i.w), e), E)
    E = map(e -> [e], E)
    V = @showprogress map(ms) do m
        ns = round.(Int, n_mean .* m)
        ns[begin] += round(Int, (m - sum(ns .* m_mono)) / m_mono[begin])
        elements = []
        for (n, e) in zip(ns, E)
            while length(e) < n push!(e, conv_ipv(e[begin], e[end])) end
            n > 0 && push!(elements, e[n])
        end
        v = foldl(conv_ipv, sort(elements, by=length))
        s = sum(x -> x.w, v)
        v = map(x -> (m=x.m, w=x.w/s), v)
        s = 0.0
        i = 0
        while s < trunc
            i += 1
            s += v[i].w
        end
        return map(x -> (m=x.m, w=x.w/s), v[begin:i])
    end
    return V
end

build_ipv(fname=joinpath(homedir(), ".PepIso", "IPV.bson"), r=1:20000, trunc=0.99) = begin
    if isfile(fname)
        @info "IPV loading from " * fname
        BSON.@load fname V
    else
        @info "IPV building"
        E = [Iso_H, Iso_C, Iso_N, Iso_O, Iso_S]
        n_mean = sum(a -> a.n .* a.w, aa_elements)
        V = calc_ipv(r, E, n_mean, trunc)
        @info "IPV caching as " * fname
        mkpath(dirname(fname))
        BSON.@save fname V
    end
    return V
end

ipv_m(m::Number, V) = map(x -> x.m, V[round(Int, m)])
ipv_m(ion, V) = ipv_m(ion.mz * ion.z, V)
ipv_w(m::Number, V) = map(x -> x.w, V[round(Int, m)])
ipv_w(ion, V) = ipv_w(ion.mz * ion.z, V)
ipv_mz(mz, z, n, V) = mz + V[round(Int, mz * z)][n].m / z
ipv_mz(ion, n, V) = ipv_mz(ion.mz, ion.z, n, V)
ipv_mz(ion, V) = [ion.mz + x.m / ion.z for x in V[round(Int, ion.mz * ion.z)]]
