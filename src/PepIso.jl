module PepIso

import GLPK
import Graphs
import JuMP
import ProgressMeter: @showprogress

import ..MesMS

prefilter(ion, spec, ε, V, mode=:mono) = begin
    f = x -> (x > 0) && !isempty(MesMS.query_ε(spec, MesMS.ipv_mz(ion, x, V), ε))
    if mode == :mono
        return f(1) && f(2)
    elseif mode == :max
        i = argmax(MesMS.ipv_w(ion, V))
        return f(i) && (f(i + 1) || f(i - 1))
    else
        @error "unknown mode: $(mode)"
        return false
    end
end

exclude(ions, xs, ms, τ, ε, V) = begin
    map(zip(ions, xs, ms)) do (i, x, m)
        (x <= 0.0 || m <= 0.0) && return false
        for (i_, x_) in zip(ions, xs)
            if i_.z % i.z == 0 && x < x_ * τ
                MesMS.in_moe(i.mz, i_.mz, ε) && i_ !== i && return false
                MesMS.in_moe(i.mz, MesMS.ipv_mz(i_, 2, V), ε) && return false
                MesMS.in_moe(i.mz, MesMS.ipv_mz(i_, 3, V), ε) && return false
            end
        end
        return true
    end
end

build_constraints_lp(ions, spec, ε, V) = begin
    items = [(; i, j, mz=MesMS.ipv_mz(ion, j, V)) for (i, ion) in enumerate(ions) for j in eachindex(MesMS.ipv_w(ion, V))]
    cs = map(p -> (; p.mz, p.inten, slots=empty(items)), spec)
    for i in items
        l, r = searchsortedfirst(spec, (1 - ε) * i.mz), searchsortedlast(spec, (1 + ε) * i.mz)
        εs = map(p -> abs(i.mz - p.mz), spec[l:r])
        l <= r && push!(cs[argmin(εs)+l-1].slots, i)
        l > r && push!(cs, (; i.mz, inten=0.0, slots=[i]))
    end
    return filter(c -> !isempty(c.slots), cs)
end

solve_lp(ions, cs, V) = begin
    model = JuMP.Model(GLPK.Optimizer)
    JuMP.set_optimizer_attribute(model, "msg_lev", GLPK.GLP_MSG_OFF)
    JuMP.@variable(model, x[eachindex(ions)] >= 0)
    JuMP.@variable(model, δ[eachindex(cs)] >= 0) # the constraint is redundant but faster
    JuMP.@objective(model, Min, sum(δ))
    s = map(c -> sum(s -> x[s.i] * MesMS.ipv_w(ions[s.i], V)[s.j], c.slots), cs)
    i = map(c -> c.inten, cs)
    JuMP.@constraint(model, i .- s .<= δ)
    JuMP.@constraint(model, s .- i .<= δ)
    JuMP.optimize!(model)
    return (xs=JuMP.value.(x), δs=i .- JuMP.value.(s))
end

calc_match_lp(xs, cs, δs) = begin
    map(eachindex(xs)) do i
        s = map(c -> any(s -> s.i == i, c.slots), cs)
        return 1 - sum(abs, δs[s], init=0.0) / sum(c -> c.inten, cs[s], init=1.0e-16)
    end
end

calc_mono_error_lp(xs, cs, δs) = begin
    map(eachindex(xs)) do i
        s = map(c -> any(s -> s.i == i, c.slots), cs)
        for (c, δ) in zip(cs[s], δs[s])
            any(s -> s.j == 1, c.slots) && return abs(δ) / c.inten
        end
        return Inf
    end
end

evaluate_lp(ions, spec, ε, V) = begin
    cs = build_constraints_lp(ions, spec, ε, V)
    xs, δs = solve_lp(ions, cs, V)
    ms = calc_match_lp(xs, cs, δs)
    es = calc_mono_error_lp(xs, cs, δs)
    return xs, ms, es
end

build_constraints_rlp(ions, spec, ε, V) = begin
    items = [(; i, j, mz=MesMS.ipv_mz(ion, j, V)) for (i, ion) in enumerate(ions) for j in eachindex(MesMS.ipv_w(ion, V))]
    cs = map(p -> (; p.mz, p.inten, slots=[]), spec)
    tab = [[] for _ in ions]
    for item in items
        l, r = searchsortedfirst(spec, (1 - ε) * item.mz), searchsortedlast(spec, (1 + ε) * item.mz)
        for (k, c) in enumerate(cs[l:r])
            push!(c.slots, (; item..., k))
        end
        l > r && push!(cs, (; item.mz, inten=0.0, slots=[(; item..., k=1),]))
        push!(tab[item.i], l <= r ? Vector(l:r) : [length(cs)])
    end
    return cs, tab
end

solve_rlp(ions, cs, tab, V) = begin
    model = JuMP.Model(GLPK.Optimizer)
    JuMP.@variable(model, x[eachindex(ions)] >= 0)
    JuMP.@variable(model, u[i=eachindex(tab), j=eachindex(tab[i]), k=eachindex(tab[i][j])] >= 0)
    JuMP.@variable(model, δ[eachindex(cs)] >= 0) # the constraint is redundant but faster
    JuMP.@objective(model, Min, sum(δ))
    s = map(c -> sum(s -> u[s.i, s.j, s.k] * MesMS.ipv_w(ions[s.i], V)[s.j], c.slots, init=0.0), cs)
    i = map(c -> c.inten, cs)
    JuMP.@constraint(model, i .- s .<= δ)
    JuMP.@constraint(model, s .- i .<= δ)
    for i in eachindex(tab)
        for j in eachindex(tab[i])
            JuMP.@constraint(model, sum(u[i, j, k] for k in eachindex(tab[i][j])) == x[i])
        end
    end
    JuMP.optimize!(model)
    return (xs=JuMP.value.(x), δs=i .- JuMP.value.(s), us=JuMP.value.(u))
end

calc_match_rlp(xs, cs, δs, us, tab) = begin
    map(enumerate(xs)) do (i, x)
        inten = 0.0
        delta = 0.0
        for (j, ps) in enumerate(tab[i])
            ws = [us[i, j, k] / max(x, 1.0e-16)  for k in eachindex(ps)]
            inten += sum(map(p -> cs[p].inten, ps) .* ws)
            delta += sum(map(p -> abs(δs[p]), ps) .* ws)
        end
        return 1 - delta / max(inten, 1.0e-16)
    end
end

evaluate_rlp(ions, spec, ε, V) = begin
    cs, tab = build_constraints_rlp(ions, spec, ε, V)
    xs, δs, us = solve_rlp(ions, cs, tab, V)
    ms = calc_match_rlp(xs, cs, δs, us, tab)
    return xs, ms
end

_deisotope(f, ions, spec, τ_max, ε, V) = begin
    xs = ms = es = nothing
    τs = [τ_max / 4, τ_max / 2, τ_max]
    τ = popfirst!(τs)
    while true
        xs, ms, es = f(ions, spec, ε, V)
        s = exclude(ions, xs, ms, τ, ε, V)
        while all(s)
            isempty(τs) && @goto done
            τ = popfirst!(τs)
            s = exclude(ions, xs, ms, τ, ε, V)
        end
        ions = ions[s]
    end
    @label done
    return [(; i.mz, i.z, x, m, e) for (i, x, m, e) in zip(ions, xs, ms, es)]
end

deisotope(f, ions, spec, τ_max, ε, V; split=false) = begin
    if split
        scores = map(zip(split_ions(ions, spec, ε, V)...)) do (ions_sub, spec_sub)
            return _deisotope(f, ions_sub, spec_sub, τ_max, ε, V)
        end
        return reduce(vcat, scores; init=typeof((; mz=0.0, z=0, x=0.0, m=0.0))[])
    else
        return _deisotope(f, ions, spec, τ_max, ε, V)
    end
end

deisotope(ions, spec, τ_max, ε, V; split=false) = deisotope(evaluate_lp, ions, spec, τ_max, ε, V; split)

split_ions(ions, spec, ε, V) = begin
    items = [(; i, n, mz=MesMS.ipv_mz(ion, n, V)) for (i, ion) in enumerate(ions) for n in eachindex(MesMS.ipv_w(ion, V))]
    cs = map(p -> (; p.mz, p.inten, slots=empty(items)), spec)
    for i in items
        l, r = searchsortedfirst(spec, (1 - ε) * i.mz), searchsortedlast(spec, (1 + ε) * i.mz)
        εs = map(p -> abs(i.mz - p.mz), spec[l:r])
        l <= r && push!(cs[argmin(εs)+l-1].slots, i)
    end
    cs = filter(c -> !isempty(c.slots), cs)
    g = Graphs.SimpleGraph(length(ions))
    for c in cs
        a = c.slots[begin]
        for b in c.slots[begin+1:end]
            Graphs.add_edge!(g, a.i, b.i)
        end
    end
    coms = Graphs.connected_components(g)
    tab = zeros(Int, length(ions))
    for (i, com) in enumerate(coms)
        tab[com] .= i
    end
    slices = map(_ -> MesMS.Peak[], coms)
    for c in cs
        push!(slices[tab[c.slots[begin].i]], MesMS.Peak(c.mz, c.inten))
    end
    return map(idxs -> ions[idxs], coms), slices
end

group_ions(I, gap, ε) = begin
    G = []
    z_max = [maximum(map(i -> i.z, ions)) for ions in I if !isempty(ions)] |> maximum
    z_min = [minimum(map(i -> i.z, ions)) for ions in I if !isempty(ions)] |> minimum
    for z in z_min:z_max
        @info "grouping (charge state: $(z))"
        tmp = []
        gs = ones(Int, length(tmp))
        @showprogress for ions in I
            s = gs .> gap
            append!(G, tmp[s])
            tmp, gs = tmp[.!s], gs[.!s]
            for ion in filter(i -> i.z == z, ions)
                grouped = false
                for (j, t) in enumerate(tmp)
                    if MesMS.in_moe(ion.mz, t[end].mz, ε)
                        push!(t, ion)
                        grouped = true
                        gs[j] = 0
                        break
                    end
                end
                if !grouped
                    push!(tmp, [ion])
                    push!(gs, 0)
                end
            end
            gs .+= 1
        end
        append!(G, tmp)
    end
    return G
end

end
