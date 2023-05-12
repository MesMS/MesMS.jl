module PyPlotExt

using Statistics
using Printf

import MesMS
import PyCall: PyCall, @py_str

__init__() = begin
    MesMS.PyPlot.mpl = PyCall.pyimport("matplotlib")
    MesMS.PyPlot.plt = PyCall.pyimport("matplotlib.pyplot")
    MesMS.PyPlot.plt_venn = PyCall.pyimport("matplotlib_venn")

    MesMS.PyPlot.mpl.style.use("ggplot")
    MesMS.PyPlot.mpl.rc("text", usetex=true)
    MesMS.PyPlot.mpl.rc("font", family="Helvetica Neue")
    MesMS.PyPlot.mpl.rc("figure", dpi=360)
end

MesMS.Plot.draw_title!(ax::PyCall.PyObject, title) = begin
    ax.set_title(title)
    ax.axis("off")
    return ax
end

MesMS.Plot.draw_index!(ax::PyCall.PyObject, text; x=-0.1, y=0.9, hide_axis=true) = begin
    ax.text(x, y, text, size=24, transform=ax.transAxes)
    if hide_axis
        ax.axis("off")
    end
    return ax
end

MesMS.Plot.draw_xlabel!(ax::PyCall.PyObject, label) = begin
    ax.set_xlabel(label)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(false) 
    return ax
end

MesMS.Plot.draw_ylabel!(ax::PyCall.PyObject, label) = begin
    ax.set_ylabel(label)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(false) 
    return ax
end

MesMS.Plot.venn2!(ax::PyCall.PyObject, a_b, b_a, ab, label_a, label_b; colors=("darkcyan", "red")) = begin
    v = MesMS.PyPlot.plt_venn.venn2((a_b, b_a, ab), set_labels=(label_a, label_b), set_colors=colors, ax=ax)
    for text in v.set_labels
        isnothing(text) && continue
        text.set_fontsize(12)
    end
    for text in v.subset_labels
        isnothing(text) && continue
        text.set_fontsize(12)
    end
    v.subset_labels[1].set_va("bottom")
    v.subset_labels[2].set_va("bottom")
    v.subset_labels[3].set_va("top")
    return ax
end

MesMS.Plot.venn3!(ax::PyCall.PyObject, a, b, c, labels; colors=("darkcyan", "red", "skyblue")) = begin
    v = MesMS.PyPlot.plt_venn.venn3([py"set($a)", py"set($b)", py"set($c)"], set_labels=labels, set_colors=colors, ax=ax)
    for text in v.set_labels
        isnothing(text) && continue
        text.set_fontsize(12)
    end
    for text in v.subset_labels
        isnothing(text) && continue
        text.set_fontsize(12)
    end
    return ax
end

_plot_seq!(ax, x, y, seq, mods, ions, font, vspace; frag_error=true) = begin
    colors = [:black for _ in seq]
    foreach(m -> colors[m[2]] = :red, mods)
    for (i, aa) in enumerate(seq)
        ax.text(x + i - 0.5, y, string(aa); fontsize=font, color=colors[i], ha=:center)
    end
    n = Dict()
    for i in filter(i -> i.peak > 0, ions)
        n[(i.part, i.loc)] = get(n, (i.part, i.loc), 0) + 1
        if i.part == :l
            ax.text(x + i.loc - 0.1, y - vspace / 3, "]"; color=i.color, fontsize=font, ha=:center)
            ax.text(x + i.loc - 0.1, y - (n[(i.part, i.loc)] / 2 + 0.5) * vspace, i.text; fontsize=font ÷ 3, color=i.color, ha=:right)
        elseif i.part == :r
            ax.text(x + i.loc + 0.1, y + vspace / 3, "["; color=i.color, fontsize=font, ha=:center)
            ax.text(x + i.loc + 0.1, y + (n[(i.part, i.loc)] / 2 + 0.6) * vspace, i.text; fontsize=font ÷ 3, color=i.color, ha=:left)
        end
    end
    if frag_error
        errors = [i.error for i in ions]
        text = @sprintf("\$\\mu=%+0.2f\$ ppm\n\$\\sigma=\\pm %0.2f\$ ppm", mean(errors), std(errors))
        ax.text(x + length(seq) + 0.5, y, text; fontsize=font ÷ 2, ha="left")
    end
    return ax
end

MesMS.Plot.seq!(ax::PyCall.PyObject, seq, mods, ions; n=24, font=16, vspace=1, frag_error=true) = begin
    seq = collect(seq)
    width = length(seq)
    font = floor(Int, n * font / max(n, width))
    vspace = vspace * font / 12
    _plot_seq!(ax, 0, 0, seq, mods, ions, font, vspace; frag_error)
    margin = max(4, (n - width) / 2)
    ax.set_xlim(-margin, width + margin)
    ax.set_ylim(-2, 2)
    ax.axis("off")
    return ax
end

MesMS.Plot.spec!(ax::PyCall.PyObject, spec, ions; vspace=5.0, hspace=1/20, font=8, linewidth=1, color=:grey, alpha=0.6) = begin
    xs = [p.mz for p in spec]
    xs_ = (xs .- xs[begin]) ./ max((xs[end] - xs[begin]), 1e-16)
    ys = [p.inten for p in spec] ./ (maximum(p -> p.inten, spec; init=1e-16) / 100)
    ax.vlines(xs, 0, ys; linewidth, color, alpha)
    ions = sort(filter(i -> i.peak > 0, ions); by=i -> ys[i.peak])
    for i in ions
        ax.vlines([xs[i.peak]], 0, [ys[i.peak]]; linewidth, color=i.color, alpha)
    end
    ys_ = zeros(length(spec))
    for i in ions
        ys_[i.peak] = max(maximum(ys_[MesMS.argquery_δ(xs_, xs_[i.peak], hspace / 2)]) + vspace, ys[i.peak])
        ax.text(xs[i.peak], ys_[i.peak], i.text; fontsize=font, i.color, va=:bottom, ha=:center)
    end
    ax.set_xlabel("m/z (Th)")
    ax.set_ylim(0, 114)
    ax.set_ylabel("intensity (\\%)")
    return ax
end

MesMS.Plot.error!(ax::PyCall.PyObject, ions, ε; x=:mz_exp, y=:error, alpha=0.6) = begin
    ions = filter(i -> i.peak > 0, ions)
    xs = [getproperty(i, x) for i in ions]
    ys = [getproperty(i, y) for i in ions]
    color = [:color in propertynames(i) ? i.color : :black for i in ions]
    ax.scatter(xs, ys; color, alpha)
    ax.set_xlabel("m/z (Th)")
    ax.set_ylabel("error (ppm)")
    ax.set_ylim(-ε * 1.0e6, ε * 1.0e6)
    return ax
end

end
