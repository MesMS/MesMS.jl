module PlotlyExt

using Printf
using Statistics

import MesMS
using PlotlyBase

_plot_seq!(ls, x, y, seq, mods, ions, font; frag_error=true) = begin
    colors = [:black for _ in seq]
    foreach(m -> colors[m[2]] = :red, mods)
    for (i, aa) in enumerate(seq)
        push!(ls, scatter(x=[x + i - 0.5], y=[y], mode="text", name="", text=string(aa), textposition="middle center", textfont=attr(size=font*2, color=colors[i])))
    end
    n = Dict()
    for i in filter(i -> i.peak > 0, ions)
        n[(i.part, i.loc)] = get(n, (i.part, i.loc), 0) + 1
        if i.part == :l
            push!(ls, scatter(x=[x + i.loc], y=[y - (n[(i.part, i.loc)] / 2 + 0.8)], mode="text", name="", text="┛", hovertext=i.text, textposition="top left", textfont=attr(size=font, color=i.color)))
        elseif i.part == :r
            push!(ls, scatter(x=[x + i.loc], y=[y + (n[(i.part, i.loc)] / 2 + 1.6)], mode="text", name="", text="┏", hovertext=i.text, textposition="bottom right", textfont=attr(size=font, color=i.color)))
        end
    end
    if frag_error
        errors = [i.error for i in ions]
        text = @sprintf("%+0.2f ± %0.2f ppm", mean(errors), std(errors))
        push!(ls, scatter(x=[x + length(seq) + 2], y=[y], mode="text", name="", text=text, textposition="left", textfont=attr(size=font/2)))
    end
end

MesMS.Plotly.seq(seq, mod, ions; font=18, frag_error=true) = begin
    seq = collect(seq)
    ls = AbstractTrace[]
    _plot_seq!(ls, 0, 0, seq, mod, ions, font; frag_error)
    return Plot(ls, Layout(; showlegend=false, yaxis=attr(showticklabels=false, range=(-6, 6)), xaxis=attr(showticklabels=false), height=300))
end

MesMS.Plotly.seq_xl(seqs, modss, sites, ionss; font=18, frag_error=true) = begin
    seqs = collect.(seqs)
    site_max = maximum(sites)
    xs = [site_max - site for site in sites] .+ [0.25, -0.25]
    ys = [2, -2]
    ls = [scatter(x=[site_max - 0.5], y=[0], mode="text", name="", text="/", textfont=attr(size=font*2))]
    for (x, y, seq, mods, ions) in zip(xs, ys, seqs, modss, ionss)
        _plot_seq!(ls, x, y, seq, mods, ions, font; frag_error)
    end
    return Plot(ls, Layout(; showlegend=false, yaxis=attr(showticklabels=false, range=(-6, 6)), xaxis=attr(showticklabels=false), height=300))
end

MesMS.Plotly.spec(spec, ions) = begin
    ls = map(p -> scatter(x=[p.mz, p.mz], y=[0, p.inten], mode="lines", line_color="black", name=""), spec)
    for i in ions
        push!(ls, scatter(
            x=[spec[i.peak].mz, spec[i.peak].mz], y=[0, spec[i.peak].inten],
            mode="lines+text", line_color=i.color, name="",
            text=["", i.text], hovertext=["", i.text], textposition="top", textfont=attr(color=i.color)
        ))
    end
    p1 = Plot(ls, Layout(; xaxis_title="m/z", yaxis_title="abundance"))
    xs = [i.mz for i in ions]
    ys = [i.error for i in ions]
    cs = [i.color for i in ions]
    ns = [i.text for i in ions]
    ls = [scatter(x=xs, y=ys, mode="markers", name="", line_color=cs, text=ns)]
    p2 = Plot(ls, Layout(; xaxis_title="m/z", yaxis_title="error", yaxis=attr(range=(-1.2, 1.2) .* maximum(abs, ys; init=0.0))))
    p = [p1; p2]
    relayout!(p, Layout(Subplots(rows=2, cols=1, row_heights=[1, 3], shared_xaxes=true)), showlegend=false)
    return p
end

end
