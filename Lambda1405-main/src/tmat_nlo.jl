

using LinearAlgebra
using NLsolve
using DataFrames
using CSV


include("./potential.jl")


tmat_nlo = function tmat_nlo!(w, mch, m_meson::Tuple, decons::Vector, a::Vector, b::Dict, d::Dict; n=2, ch="ch11", rs="rs11", μ_tuple=(630e0,))
    vlo = VLO!(w, mch, decons, n=n)

    # Potential at next to leading order
    mπ, mK = m_meson
    vnlo = VNLO!(w, mπ, mK, mch, b, d, decons, n=n)
    v = vlo + vnlo

    id = diagm(ones(n))

    if length(μ_tuple) != 1
        gl = diagm([Gdr!(w, mchi, a, μ=μ) for (mchi, a, μ) in zip(mch[1:n], a[1:n], μ_tuple[1:n])])
    else
        gl = diagm([Gdr!(w, mchi, a, μ=μ_tuple[1]) for (mchi, a) in zip(mch[1:n], a[1:n])])
    end

    s2 = [2im * mch[i][2] * qcm!(w, mch[i][1], mch[i][2]) / (4π * w) for i in 1:n]
    # sheets
    if rs == "rs11" || rs == "rs1111"
        nothing
    elseif rs == "rs21" || rs == "rs2111"
        gl[1, 1] += s2[1]
    elseif rs == "rs22" || rs == "rs2211"
        gl[1, 1] += s2[1]
        gl[2, 2] += s2[2]
    elseif rs == "rs2221"
        gl[1, 1] += s2[1]
        gl[2, 2] += s2[2]
        gl[3, 3] += s2[3]
    end
    t = inv(id - v * gl) * v

    if ch == "ch11"
        return t[1, 1]
    elseif ch == "ch22"
        return t[2, 2]
    elseif ch == "all"
        return t
    end
end

tdet_NLO = function tdet_NLO!(w, mch, m_meson::Tuple, decons::Vector, a::Vector, b::Dict, d::Dict; n=2, rs="rs11", μ_tuple=(630e0,))
    vlo = VLO!(w, mch, decons, n=n)

    # Potential at next to leading order
    mπ, mK = m_meson
    vnlo = VNLO!(w, mπ, mK, mch, b, d, decons, n=n)
    v = vlo + vnlo

    id = diagm(ones(n))

    if length(μ_tuple) != 1
        gl = diagm([Gdr!(w, mchi, a, μ=μ) for (mchi, a, μ) in zip(mch[1:n], a[1:n], μ_tuple[1:n])])
    else
        gl = diagm([Gdr!(w, mchi, a, μ=μ_tuple[1]) for (mchi, a) in zip(mch[1:n], a[1:n])])
    end

    s2 = [2im * mch[i][2] * qcm!(w, mch[i][1], mch[i][2]) / (4π * w) for i in 1:n]
    # sheets
    if rs == "rs11" || rs == "rs1111"
        nothing
    elseif rs == "rs21" || rs == "rs2111"
        gl[1, 1] += s2[1]
    elseif rs == "rs22" || rs == "rs2211"
        gl[1, 1] += s2[1]
        gl[2, 2] += s2[2]
    elseif rs == "rs2221"
        gl[1, 1] += s2[1]
        gl[2, 2] += s2[2]
        gl[3, 3] += s2[3]
    end

    return det(id - v * gl)
end

plot_tdet_NLO = function plot_tdet_NLO!(ax, rew, imw, params::Dict, a::Vector, b::Dict, d::Dict; n=2, rs="rs21", μ_tuple=(630e0,))
    nn = length(rew)
    abstdet = zeros(nn, nn)
    decons = params[:decons_vec]
    mch = params[:mch]
    m_meson = params[:m_meson]

    for i in 1:nn
        for j in 1:nn
            abstdet[i, j] = abs(tdet_NLO(rew[j] + imw[i] * 1im, mch, m_meson, decons, a, b, d, n=n, rs=rs, μ_tuple=μ_tuple))
        end
    end

    ax.set_ylim(minimum(imw), maximum(imw))
    ax.set_xlim(minimum(rew), maximum(rew))
    ax.set_xlabel(L"$Re[E_{cm}]$ [MeV]")
    ax.set_ylabel(L"$Im[E_{cm}]$ [MeV]")
    ax.set_title(rs)
    # Plot thrshold
    ax.vlines(params[:threshold][2], minimum(imw), maximum(imw), color=:gray, linewidth=2, label=L"$m_{\bar{K}N}$")
    ax.vlines(params[:threshold][1], minimum(imw), maximum(imw), color="k", linewidth=2, linestyles=:dashed, label=L"$m_{\pi\Sigma}$")

    l = range(minimum(abstdet), maximum(abstdet), 13)
    ax.contour(rew, imw, abstdet, colors=["#000", "#000"], linestyles="--", levels=l, linewidths=0.9)
    ax.contourf(rew, imw, abstdet, levels=l, cmap="jet")
    ax.legend()
end

"""Find the pole"""
pole_NLO = function pole_NLO!(init_x::Vector, params::Dict, a, b, d; n=2, rs="rs21", μ_tuple=(630e0,))
    mch = params[:mch]
    decons = params[:decons_vec]
    m_meson = params[:m_meson]
    sol = nlsolve(x -> cmplx!(tdet_NLO(x[1] + x[2] * 1e0im, mch, m_meson, decons, a, b, d, n=n, rs=rs, μ_tuple=μ_tuple)), init_x)

    if converged(sol) == true
        return @. round(sol.zero, digits=2)
    else
        println("The routine was not converged.")
        return [NaN, NaN]
    end
end