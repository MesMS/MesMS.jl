module Plot

import ..MesMS

include("api.jl")

const ion_a = (type=:a, sym=:a, part=:l, Δ=MesMS.Formula(C=-1, O=-1), color=:lightgreen)
const ion_a_NH3 = (type=:a_NH3, sym=Symbol("a*"), part=:l, Δ=MesMS.Formula(H=-3, C=-1, N=-1, O=-1), color=:lightgreen)
const ion_a_H2O = (type=:a_H2O, sym=Symbol("a∘"), part=:l, Δ=MesMS.Formula(H=-2, C=-1, O=-2), color=:lightgreen)
const ion_b = (type=:b, sym=:b, part=:l, Δ=MesMS.Formula(), color=:green)
const ion_b_NH3 = (type=:b_NH3, sym=Symbol("b*"), part=:l, Δ=MesMS.Formula(H=-3, N=-1), color=:green)
const ion_b_H2O = (type=:b_H2O, sym=Symbol("b∘"), part=:l, Δ=MesMS.Formula(H=-2, O=-1), color=:green)
const ion_c = (type=:c, sym=:c, part=:l, Δ=MesMS.Formula(H=3, N=1), color=:darkgreen)
const ion_x = (type=:x, sym=:x, part=:r, Δ=MesMS.Formula(C=1, O=2), color=:lightred)
const ion_y = (type=:y, sym=:y, part=:r, Δ=MesMS.Formula(H=2, O=1), color=:red)
const ion_y_NH3 = (type=:y_NH3, sym=Symbol("y*"), part=:r, Δ=MesMS.Formula(H=-1, N=-1, O=1), color=:red)
const ion_y_H2O = (type=:y_H2O, sym=Symbol("y∘"), part=:r, Δ=MesMS.Formula(), color=:red)
const ion_z = (type=:z, sym=:z, part=:r, Δ=MesMS.Formula(H=-1, N=-1, O=1), color=:darkred)

lsum_seq(seq, mods, tab_aa, tab_mod) = begin
    seq = map(x -> tab_aa[Symbol(x)], collect(seq))
    foreach(m -> seq[m[2]] += tab_mod[Symbol(m[1])], mods)
    return cumsum(seq)
end

calc_base(lsum, var=[]) = begin
    n = length(lsum)
    site_min, site_max = isempty(var) ? (0, n + 1) : extrema(v -> v[2], var)
    rsum = push!(reverse(lsum[end] .- vcat(lsum[begin:end-1])), lsum[end])
    l = [(; part=:l, idx=i, var="", mass=m) for (i, m) in enumerate(lsum[begin:site_max-1])]
    r = [(; part=:r, idx=i, var="", mass=m) for (i, m) in enumerate(rsum[begin:n-site_min])]
    for (mass, site, sym) in var
        append!(l, [(; part=:l, idx=i+site-1, var=sym, mass=m) for (i, m) in enumerate(lsum[site:end] .+ mass)])
        append!(r, [(; part=:r, idx=i+n-site, var=sym, mass=m) for (i, m) in enumerate(rsum[n-site+1:end] .+ mass)])
    end
    return (; l=map(i -> (; i..., loc=i.idx), l), r=map(i -> (; i..., loc=n - i.idx), r))
end

calc_ion(base, δ, z, type, sym, color=:green; sym_abbr=sym) = begin
    (;
        base..., type, sym, mz=(base.mass + δ + z * MesMS.mₚ) / z, z, color,
        text="\$$(sym)_{$(base.idx)}^{$(z == 1 ? "" : z)+}\$",
        text_abbr="\$$(sym_abbr)_{$(base.idx)}^{$(z == 1 ? "" : z)+}\$",
    )
end

match_peak(spec, ion, ε; abu=false) = begin
    r = MesMS.argquery_ε(spec, ion.mz, ε)
    if isempty(r)
        return (; ion..., peak=0, mz_exp=0.0, error=0.0, inten=0.0)
    else
        p = abu ? argmax(i -> spec[i].inten, r) : argmin(i -> abs(spec[i].mz - ion.mz), r)
        return (; ion..., peak=p, mz_exp=spec[p].mz, error=MesMS.error_ppm(ion.mz, spec[p].mz), inten=spec[p].inten)
    end
end

check_dup(ionss...; color=:grey) = begin
    c = Dict{Int, Int}()
    foreach(x -> foreach(i -> c[i.peak] = get(c, i.peak, 0) + 1, x), ionss)
    ionss = map(x -> [(; i..., color=(c[i.peak] == 1 ? i.color : color)) for i in x], ionss)
    return length(ionss) > 1 ? ionss : ionss[begin]
end

build_ions(spec, seq, mods, ε, tab_ele, tab_aa, tab_mod; types=[(ion_b, 1:2), (ion_y, 1:2)], abu=false) = begin
    base = calc_base(lsum_seq(seq, mods, tab_aa, tab_mod))
    ions = [calc_ion(i, MesMS.calc_mass(t.Δ, tab_ele), z, t.type, t.sym, t.color) for (t, zs) in types for z in zs for i in base[t.part]]
    ions = [match_peak(spec, ion, ε; abu) for ion in ions]
    return ions
end

build_ions_xl(spec, seqs, modss, linker, sites, ε, tab_ele, tab_aa, tab_mod; types=[(ion_b, 1:3), (ion_y, 1:3)], abu=false) = begin
    lsums = [lsum_seq(seq, mod, tab_aa, tab_mod) for (seq, mod) in zip(seqs, modss)]
    anothers = reverse([s[end] for s in lsums])
    ions = map(zip((:α, :β), lsums, sites, anothers, (0.25, -0.25), (:circle, :diamond))) do (pep_sym, seq, site, another, offset, shape)
        if linker.cleavable
            base = calc_base(seq, [(m, site, "'"^i) for (i, m) in enumerate(linker.masses)])
        else
            seq[site:end] .+= another + linker.mass + MesMS.calc_mass(MesMS.Formula(H=2, O=1), tab_ele)
            base = calc_base(seq)
        end
        ions = [calc_ion(i, MesMS.calc_mass(t.Δ, tab_ele), z, t.type, "$(pep_sym)$(i.var)$(t.sym)", t.color; sym_abbr="$(i.var)$(t.sym)") for (t, zs) in types for z in zs for i in base[t.part]]
        ions = [match_peak(spec, ion, ε; abu) for ion in ions]
        x = maximum(sites) - site
        return [(; ion..., loc_=ion.loc + x + offset, shape) for ion in ions]
    end
    return ions
end

end
