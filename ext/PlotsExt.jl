module PlotsExt

using Statistics
using Printf

import MesMS
import Plots

_plot_seq!(x, y, seq, mods, ions, font, vspace; frag_error=true) = begin
    colors = [:black for _ in seq]
    foreach(m -> colors[m[2]] = :red, mods)
    for (i, aa) in enumerate(seq)
        Plots.annotate!((x + i - 0.5, y, (string(aa), font, colors[i])))
    end
    n = Dict()
    for i in filter(i -> i.peak > 0, ions)
        n[(i.part, i.loc)] = get(n, (i.part, i.loc), 0) + 1
        if i.part == :l
            Plots.annotate!((x + i.loc - 0.1, y - vspace / 3, ("\$\\lrcorner\$", i.color, font)))
            Plots.annotate!((x + i.loc - 0.5, y - (n[(i.part, i.loc)] / 2 + 0.5) * vspace, (i.text, font ÷ 3, i.color, :left)))
        elseif i.part == :r
            Plots.annotate!((x + i.loc + 0.1, y + vspace / 3, ("\$\\ulcorner\$", i.color, font)))
            Plots.annotate!((x + i.loc, y + (n[(i.part, i.loc)] / 2 + 0.6) * vspace, (i.text, font ÷ 3, i.color, :left)))
        end
    end
    if frag_error
        errors = [i.error for i in ions]
        text = @sprintf("\$\\mu=%+0.2f\$ ppm\n\$\\sigma=\\pm %0.2f\$ ppm", mean(errors), std(errors))
        Plots.annotate!((x + length(seq) + 0.5, y, (text, font ÷ 3, :left)))
    end
end

MesMS.Plot.seq!(p::Plots.Plot, seq, mods, ions; n=32, font=12, vspace=1, frame=:none, frag_error=true) = begin
    seq = collect(seq)
    width = length(seq)
    font = floor(Int, n * font / max(n, width))
    vspace = vspace * font / 12
    _plot_seq!(0, 0, seq, mods, ions, font, vspace; frag_error)
    margin = max(4, (n - width) / 2)
    return Plots.plot!(p; xlim=(-margin, width + margin), ylim=(-2, 2), frame)
end

MesMS.Plot.seq_xl!(p::Plots.Plot, seqs, modss, sites, ionss; n=32, font=12, vspace=1, frame=:none, frag_error=true) = begin
    seqs = collect.(seqs)
    site_max = maximum(sites)
    xs = [site_max - site for site in sites] .+ [0.25, -0.25]
    ys = [2, -2]
    width = maximum(((seq, Δx),) -> length(seq) + Δx, zip(seqs, xs))
    font = floor(Int, n * font / max(n, width))
    vspace = vspace * font / 12
    Plots.annotate!((site_max - 0.5, 0, ("/", font)))
    for (x, y, seq, mods, ions) in zip(xs, ys, seqs, modss, ionss)
        _plot_seq!(x, y, seq, mods, ions, font, vspace; frag_error)
    end
    margin = max(4, (n - width) / 2)
    return Plots.plot!(p; xlim=(-margin, width + margin), ylim=(-3, 4), frame)
end

MesMS.Plot.spec!(p::Plots.Plot, spec, ions; vspace=5.0, hspace=1/20, font=4, linewidth=0.5, linecolor=:grey, linealpha=0.6, frame=:box) = begin
    xs = [p.mz for p in spec]
    xs_ = (xs .- xs[begin]) ./ max((xs[end] - xs[begin]), 1e-16)
    ys = [p.inten for p in spec] ./ (maximum(p -> p.inten, spec; init=1e-16) / 100)
    Plots.plot!(p, xs, ys; st=:sticks, label=nothing, linewidth, linecolor, linealpha)
    ions = sort(filter(i -> i.peak > 0, ions); by=i -> ys[i.peak])
    for i in ions
        Plots.plot!(p, [xs[i.peak]], [ys[i.peak]]; st=:sticks, label=nothing, linewidth, linecolor=i.color, linealpha)
    end
    ys_ = zeros(length(spec))
    for i in ions
        ys_[i.peak] = max(maximum(ys_[MesMS.argquery_δ(xs_, xs_[i.peak], hspace / 2)]) + vspace, ys[i.peak])
        Plots.annotate!((xs[i.peak], ys_[i.peak], (i.text, font, i.color, :bottom)))
    end
    return Plots.plot!(p; ylim=(0, 114), frame)
end

MesMS.Plot.error!(p::Plots.Plot, ions, ε; x=:mz_exp, y=:error, label=nothing, markersize=3, markeralpha=0.6, xformatter=_ -> "", frame=:box) = begin
    ions = filter(i -> i.peak > 0, ions)
    xs = [getproperty(i, x) for i in ions]
    ys = [getproperty(i, y) for i in ions]
    markershape = [:shape in propertynames(i) ? i.shape : :circle for i in ions]
    markercolor = [:color in propertynames(i) ? i.color : :black for i in ions]
    Plots.plot!(p, xs, ys; st=:scatter, label, markershape, markercolor, markersize, markeralpha, markerstrokewidth=0, markerstrokealpha=0.0)
    return Plots.plot!(p; ylim=(-ε * 1.0e6, ε * 1.0e6), xformatter, frame)
end

MesMS.Plots.psm(spec, seq, mods, ε, ions; dup=true) = begin
    p_seq = MesMS.Plot.seq!(Plots.plot(), seq, mods, dup ? MesMS.Plot.check_dup(ions) : ions)
    p_seq_err = MesMS.Plot.error!(Plots.plot(), ions, ε; x=:loc)
    Plots.xlims!(p_seq_err, Plots.xlims(p_seq))
    p_spec = MesMS.Plot.spec!(Plots.plot(), spec, ions)
    p_spec_err = MesMS.Plot.error!(Plots.plot(), ions, ε; x=:mz_exp)
    Plots.xlims!(p_spec_err, Plots.xlims(p_spec))
    return Plots.plot(p_seq, p_seq_err, p_spec, p_spec_err; layout=Plots.grid(4, 1, heights=[0.2, 0.15, 0.5, 0.15]))
end

MesMS.Plots.psm(spec, seq, mods, ε, tab_ele, tab_aa, tab_mod; types=[(MesMS.Plot.ion_b, 1:2), (MesMS.Plot.ion_y, 1:2)], dup=true, abu=false) = begin
    ions = MesMS.Plot.build_ions(spec, seq, mods, ε, tab_ele, tab_aa, tab_mod; types, abu)
    return MesMS.Plots.psm(spec, seq, mods, ε, ions; dup)
end

MesMS.Plots.psm_xl(spec, seqs, modss, sites, ε, ions; dup=true) = begin
    p_seq = MesMS.Plot.seq_xl!(Plots.plot(), seqs, modss, sites, dup ? MesMS.Plot.check_dup(ions...) : ions)
    ions = vcat(ions...)
    p_seq_err = MesMS.Plot.error!(Plots.plot(), ions, ε; x=:loc_)
    Plots.xlims!(p_seq_err, Plots.xlims(p_seq))
    p_spec = MesMS.Plot.spec!(Plots.plot(), spec, ions)
    p_spec_err = MesMS.Plot.error!(Plots.plot(), ions, ε; x=:mz_exp)
    Plots.xlims!(p_spec_err, Plots.xlims(p_spec))
    return Plots.plot(p_seq, p_seq_err, p_spec, p_spec_err; layout=Plots.grid(4, 1, heights=[0.3, 0.15, 0.4, 0.15]))
end

MesMS.Plots.psm_xl(spec, seqs, modss, linker, sites, ε, tab_ele, tab_aa, tab_mod; types=[(MesMS.Plot.ion_b, 1:3), (MesMS.Plot.ion_y, 1:3)], dup=true, abu=false) = begin
    ions = MesMS.Plot.build_ions_xl(spec, seqs, modss, linker, sites, ε, tab_ele, tab_aa, tab_mod; types, abu)
    return MesMS.Plots.psm_xl(spec, seqs, modss, sites, ε, ions; dup)
end

MesMS.Plot.peak!(p::Plots.Plot, spec, mzs, base, ε) = begin
    xs = [p.mz for p in spec]
    ys = [p.inten for p in spec]
    ys = ys ./ maximum(ys; init=1.0e-16) * 100
    base = base ./ maximum(base; init=1.0e-16) * 100
    Plots.plot!(p, xs, ys; st=:sticks, label=nothing, linealpha=0.5, linecolor=:grey)
    ions = [MesMS.Plot.match_peak(spec, (; idx, mz), ε) for (idx, mz) in enumerate(mzs)]
    xs = [i.peak > 0 ? xs[i.peak] : i.mz for i in ions]
    ys = [i.peak > 0 ? ys[i.peak] : 0.0 for i in ions]
    Plots.plot!(p, xs, ys; st=:sticks, label=nothing, linecolor=:green)
    Plots.plot!(p, xs[begin:begin]; st=:vline, label=nothing, linealpha=0.5, linecolor=:red, linestyle=:dot)
    Plots.plot!(p, xs, base; st=:scatter, label=nothing, markeralpha=0.5, markerstrokewidth=0)
    return Plots.plot!(p; ylims=(0, 114), grid=:none)
end

MesMS.Plot.peak_rt!(p::Plots.Plot, mzs, specs, ε; n=nothing, spec=nothing) = begin
    for (i, mz) in reverse(collect(enumerate(mzs)))
        xs = [m.retention_time for m in specs]
        ys = zeros(length(specs)) .+ mz
        ions = [MesMS.Plot.match_peak(s.peaks, (; mz), ε; abu=true) for s in specs]
        zs = [i.peak > 0 ? s.peaks[i.peak].inten : 0.0 for (i, s) in zip(ions, specs)]
        style = i == n ? (linewidth=1, linealpha=1) : (linewidth=0.5, linealpha=0.6)
        Plots.plot!(p, xs, ys, zs; label=nothing, style...)
    end
    if !isnothing(spec)
        ys = fill(spec.retention_time, length(mzs))
        ions = [MesMS.Plot.match_peak(spec.peaks, (; mz), ε; abu=true) for mz in mzs]
        zs = [i.peak > 0 ? spec.peaks[i.peak].inten : 0.0 for i in ions]
        Plots.plot!(p, ys, mzs, zs; st=:scatter, label=nothing, markeralpha=1, markercolor=:red, markerstrokewidth=0)
    end
    return Plots.plot!(p; camera=(60, 15), showaxis=false, tickfontsize=2, grid=:none)
end

end
