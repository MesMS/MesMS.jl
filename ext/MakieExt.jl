module MakieExt

using Statistics
using Printf

import Makie
import MesMS

__init__() = begin
    Makie.update_theme!(
        font = "Helvetica Neue",
        Axis = (
            backgroundcolor = :grey90,
            topspinevisible = false,
            bottomspinevisible = false,
            leftspinevisible = false,
            rightspinevisible = false,
            xgridcolor = :white,
            ygridcolor = :white,
        ),
    )
end

MesMS.Plot.pie!(ax::Makie.Axis, xs, labels; colors=Makie.cgrad(:Set2_8)[eachindex(xs)], text_radius=1.3, lim=1.6) = begin
    Makie.pie!(ax, xs, color=colors, inner_radius=0.5, strokecolor=:white, strokewidth=4)
    s = sum(xs)
    for i in eachindex(xs)
        θ = (sum(xs[1:i]) - xs[i] / 2) / s * 2π
        f = xs[i] / s * 100
        Makie.text!(ax, text_radius * cos(θ), text_radius * sin(θ);
            text=@sprintf("%s\n%s (%.2f%%)", labels[i], xs[i], f), align=(:center, :center),
        )
        Makie.text!(ax, 0, 0; text="total: $(s)", align=(:center, :center))
    end
    ax.aspect = 1
    ax.limits = (-1, 1, -1, 1) .* lim
    ax.backgroundcolor = :white
    Makie.hidedecorations!(ax)
    Makie.hidespines!(ax)
    return ax
end

_plot_seq!(ax::Makie.Axis, x, y, seq, mods, ions, fontsize, vspace; frag_error=true) = begin
    colors = [:black for _ in seq]
    foreach(m -> colors[m[2]] = :red, mods)
    for (i, aa) in enumerate(seq)
        Makie.text!(ax, x + i - 0.5, y - vspace/3; text=string(aa), fontsize, color=colors[i], align=(:center, :center))
    end
    n = Dict()
    for i in filter(i -> i.peak > 0, ions)
        n[(i.part, i.loc)] = get(n, (i.part, i.loc), 0) + 1
        if i.part == :l
            Makie.text!(ax, x + i.loc + 0.1, y - 3/2 * vspace; text="◞", i.color, fontsize=fontsize, align=(:right, :bottom))
            Makie.text!(ax, x + i.loc, y - (n[(i.part, i.loc)] + 3) / 2 * vspace; text=Makie.L"%$(i.tex_abbr)", fontsize=fontsize/3, i.color, align=(:right, :bottom))
        elseif i.part == :r
            Makie.text!(ax, x + i.loc - 0.1, y + 3/2 * vspace; text="◜", i.color, fontsize=fontsize, align=(:left, :top))
            Makie.text!(ax, x + i.loc, y + (n[(i.part, i.loc)] + 3) / 2 * vspace; text=Makie.L"%$(i.tex_abbr)", fontsize=fontsize/3, i.color, align=(:left, :top))
        end
    end
    if frag_error
        errors = [i.error for i in ions]
        s1, s2 = @sprintf("%+0.2f", mean(errors)), @sprintf("%0.2f", std(errors))
        Makie.text!(ax, x + length(seq) + 0.5, y; text=Makie.L"$μ=%$(s1)$ ppm", fontsize=fontsize / 3, align=(:left, :bottom))
        Makie.text!(ax, x + length(seq) + 0.5, y; text=Makie.L"$σ=±%$(s2)$ ppm", fontsize=fontsize / 3, align=(:left, :top))
    end
end

MesMS.Plot.seq!(ax::Makie.Axis, seq, mods, ions; n=32, fontsize=24, vspace=1, frag_error=true) = begin
    seq = collect(seq)
    width = length(seq) + frag_error * 4
    fontsize *= n / max(n, width)
    n = max(n, width)
    vspace = vspace * fontsize / 24
    _plot_seq!(ax, 0, 0, seq, mods, ions, fontsize, vspace; frag_error)
    Makie.scatter!([-(n-width)/2, -(n-width)/2, (n+width)/2, (n+width)/2], [-4, 4, 4, -4]; color=(:black, 0.0))
    ax.backgroundcolor = :white
    Makie.hidedecorations!(ax)
    Makie.hidespines!(ax)
    return ax
end

MesMS.Plot.seq_xl!(ax::Makie.Axis, seqs, modss, sites, ionss; n=32, fontsize=24, vspace=1, frag_error=true) = begin
    seqs = collect.(seqs)
    site_max = maximum(sites)
    xs = [site_max - site for site in sites] .+ [0.25, -0.25]
    ys = [2, -2]
    width = maximum(length.(seqs) .+ xs) + frag_error * 4
    fontsize *= n / max(n, width)
    n = max(n, width)
    vspace = vspace * fontsize / 24
    Makie.text!(ax, site_max - 0.5, -vspace/3; text="/", fontsize, align=(:center, :center))
    for (x, y, seq, mods, ions) in zip(xs, ys, seqs, modss, ionss)
        _plot_seq!(ax, x, y, seq, mods, ions, fontsize, vspace; frag_error)
    end
    Makie.scatter!([-(n-width)/2, -(n-width)/2, (n+width)/2, (n+width)/2], [-4, 4, 4, -4]; color=(:black, 0.0))
    ax.backgroundcolor = :white
    Makie.hidedecorations!(ax)
    Makie.hidespines!(ax)
    return ax
end

MesMS.Plot.spec!(ax::Makie.Axis, spec, ions; vspace=5.0, hspace=1/20, fontsize=8, linewidth=2, linecolor=:grey, linealpha=0.6) = begin
    xs = [p.mz for p in spec]
    xs_ = (xs .- xs[begin]) ./ max((xs[end] - xs[begin]), 1e-16)
    ys = [p.inten for p in spec] ./ (maximum(p -> p.inten, spec; init=1e-16) / 100)
    ions = sort(filter(i -> i.peak > 0, ions); by=i -> ys[i.peak])
    f = trues(length(xs))
    foreach(i -> f[i.peak] = false, ions)
    Makie.rangebars!(ax, xs[f], 0, ys[f]; linewidth=linewidth/2, color=(linecolor, linealpha))
    for i in ions
        Makie.rangebars!(ax, [xs[i.peak]], 0, [ys[i.peak]]; linewidth, i.color)
    end
    ys_ = zeros(length(spec))
    for i in ions
        ys_[i.peak] = max(maximum(ys_[MesMS.argquery_δ(xs_, xs_[i.peak], hspace / 2)]) + vspace, ys[i.peak])
        Makie.text!(ax, xs[i.peak], ys_[i.peak]; text=Makie.L"%$(i.tex)", fontsize, i.color, align=(:center, :bottom))
    end
    Makie.ylims!(ax, 0, 114)
    return ax
end

MesMS.Plot.error!(ax::Makie.Axis, ions, ε; x=:mz_exp, y=:error, markersize=8, markeralpha=0.6) = begin
    ions = filter(i -> i.peak > 0, ions)
    xs = [getproperty(i, x) for i in ions]
    ys = [getproperty(i, y) for i in ions]
    marker = [:shape in propertynames(i) ? i.shape : :circle for i in ions]
    color = [:color in propertynames(i) ? (i.color, markeralpha) : (:black, markeralpha) for i in ions]
    if !isempty(ions)
        Makie.scatter!(ax, xs, ys; marker, color, markersize)
    end
    Makie.ylims!(ax, -ε * 1.0e6 * 1.19, ε * 1.0e6 * 1.19)
    return ax
end

MesMS.Plot.psm!(fig::Makie.Figure, spec, seq, mods, ε, ions; dup=true) = begin
    ax_seq = MesMS.Plot.seq!(Makie.Axis(fig[1, 1]), seq, mods, dup ? MesMS.Plot.check_dup(ions) : ions)
    ax_seq_err = MesMS.Plot.error!(Makie.Axis(fig[2, 1]), ions, ε; x=:loc)
    Makie.hidexdecorations!(ax_seq_err; grid=false)
    Makie.linkxaxes!(ax_seq_err, ax_seq)
    ax_spec = MesMS.Plot.spec!(Makie.Axis(fig[3, 1]), spec, ions)
    ax_spec_err = MesMS.Plot.error!(Makie.Axis(fig[4, 1]), ions, ε; x=:mz_exp)
    Makie.linkxaxes!(ax_spec_err, ax_spec)
    Makie.hidexdecorations!(ax_spec; grid=false)
    Makie.rowsize!(fig.layout, 1, Makie.Auto(2))
    Makie.rowsize!(fig.layout, 2, Makie.Auto(1.5))
    Makie.rowsize!(fig.layout, 3, Makie.Auto(5))
    Makie.rowsize!(fig.layout, 4, Makie.Auto(1.5))
    return fig
end

MesMS.Plot.psm!(fig::Makie.Figure, spec, seq, mods, ε, tab_ele, tab_aa, tab_mod; types=[(MesMS.Plot.ion_b, 1:2), (MesMS.Plot.ion_y, 1:2)], dup=true, abu=false) = begin
    ions = MesMS.Plot.build_ions(spec, seq, mods, ε, tab_ele, tab_aa, tab_mod; types, abu)
    return MesMS.Plot.psm!(fig, spec, seq, mods, ε, ions; dup)
end

MesMS.Plot.psm_xl!(fig::Makie.Figure, spec, seqs, modss, sites, ε, ions; dup=true) = begin
    ax_seq = MesMS.Plot.seq_xl!(Makie.Axis(fig[1, 1]), seqs, modss, sites, dup ? MesMS.Plot.check_dup(ions...) : ions)
    ions = vcat(ions...)
    ax_seq_err = MesMS.Plot.error!(Makie.Axis(fig[2, 1]), ions, ε; x=:loc_)
    Makie.hidexdecorations!(ax_seq_err; grid=false)
    Makie.linkxaxes!(ax_seq_err, ax_seq)
    ax_spec = MesMS.Plot.spec!(Makie.Axis(fig[3, 1]), spec, ions)
    ax_spec_err = MesMS.Plot.error!(Makie.Axis(fig[4, 1]), ions, ε; x=:mz_exp)
    Makie.linkxaxes!(ax_spec_err, ax_spec)
    Makie.hidexdecorations!(ax_spec; grid=false)
    Makie.rowsize!(fig.layout, 1, Makie.Auto(4))
    Makie.rowsize!(fig.layout, 2, Makie.Auto(1.5))
    Makie.rowsize!(fig.layout, 3, Makie.Auto(4))
    Makie.rowsize!(fig.layout, 4, Makie.Auto(1.5))
    return fig
end

MesMS.Plot.psm_xl!(fig::Makie.Figure, spec, seqs, modss, linker, sites, ε, tab_ele, tab_aa, tab_mod; types=[(MesMS.Plot.ion_b, 1:3), (MesMS.Plot.ion_y, 1:3)], dup=true, abu=false) = begin
    ions = MesMS.Plot.build_ions_xl(spec, seqs, modss, linker, sites, ε, tab_ele, tab_aa, tab_mod; types, abu)
    return MesMS.Plot.psm_xl!(fig, spec, seqs, modss, sites, ε, ions; dup)
end

end
