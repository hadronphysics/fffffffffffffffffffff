

using LinearAlgebra
using NLsolve
using DataFrames
using CSV
#import Plots


include("./loop_fun.jl")
include("./potential.jl")

function tmat!(w, mch, decons, a::Vector; n=2, ch="ch11", rs="rs21", μ_tuple=(630e0,))
    v = VLO!(w, mch, decons, n=n)

    id = diagm(ones(n) )

    if length(μ_tuple) != 1
        gl = diagm([Gdr!(w, mchi, a, μ=μ) for (mchi, a, μ) in zip(mch[1:n], a[1:n], μ_tuple[1:n] ) ] )
    else
        gl = diagm([Gdr!(w, mchi, a, μ=μ_tuple[1]) for (mchi, a) in zip(mch[1:n], a[1:n]) ] )
    end

    s2 = [2im * mch[i][2] * qcm!(w, mch[i][1], mch[i][2]) / (4π*w) for i in 1:n]
    # sheets
    if rs == "rs11" || rs == "rs1111"
        nothing
    elseif rs == "rs21" || rs == "rs2111"
        gl[1, 1] += s2[1]
    elseif rs == "rs22" || rs == "rs2211"
        gl[1, 1] += s2[1]
        gl[2, 2] += s2[2]
    elseif rs == "rs12"
        gl[2, 2] = s2[2]
    elseif rs == "rs2221"
        gl[1, 1] += s2[1]
        gl[2, 2] += s2[2]
        gl[3, 3] += s2[3]
    end

    t = inv(id - v*gl) * v

    if ch == "ch11"
        return t[1, 1]
    elseif ch == "ch22"
        return t[2, 2]
    elseif ch == "ch33"
        return t[3, 3]
    elseif ch == "ch44"
        return t[4, 4]
    elseif ch == "all"
        return t
    end

end

tdet = function tdet!(w, mch, decons, a::Vector; n=2, rs="rs21", μ_tuple=(630e0,))
    v = VLO!(w, mch, decons, n=n)

    id = diagm(ones(n) )

    if length(μ_tuple) != 1
        gl = diagm([Gdr!(w, mchi, a, μ=μ) for (mchi, a, μ) in zip(mch[1:n], a[1:n], μ_tuple[1:n] ) ] )
    else
        gl = diagm([Gdr!(w, mchi, a, μ=μ_tuple[1]) for (mchi, a) in zip(mch[1:n], a[1:n]) ] )
    end

    s2 = [2im * mch[i][2] * qcm!(w, mch[i][1], mch[i][2]) / (4π*w) for i in 1:n]
    # sheets
    if rs == "rs11" || rs == "rs1111"
        nothing
    elseif rs == "rs21" || rs == "rs2111"
        gl[1, 1] += s2[1]
    elseif rs == "rs22" || rs == "rs2211"
        gl[1, 1] += s2[1]
        gl[2, 2] += s2[2]
    elseif rs == "rs12"
        gl[2, 2] = s2[2]
    elseif rs == "rs2221"
        gl[1, 1] += s2[1]
        gl[2, 2] += s2[2]
        gl[3, 3] += s2[3]
    end
    return det(id - v*gl)
end

"""Plot tdet in an complex energy plane to verify the position of pole"""
plot_tdet_LO = function plot_tdet_LO!(ax, rew, imw, params::Dict, a::Vector; n=2, rs="rs21", μ_tuple=(630e0,))
    nn = length(rew)
    abstdet = zeros(nn, nn)
    decons = params[:decons_vec]
    mch = params[:mch]

    for i in 1:nn
        for j in 1:nn
            abstdet[i, j] = abs(tdet!(rew[j]+imw[i]*1im, mch, decons, a, n=n, rs=rs, μ_tuple=μ_tuple) )
        end
    end
    
    ax.set_ylim(minimum(imw), maximum(imw) )
    ax.set_xlim(minimum(rew), maximum(rew) )
    ax.set_xlabel(L"$Re[E_{cm}]$ [MeV]")
    ax.set_ylabel(L"$Im[E_{cm}]$ [MeV]")
    #ax.set_title(rs)
    # Plot thrshold
    ax.vlines(params[:threshold][2], minimum(imw), maximum(imw), color=:gray, linewidth=2, label=L"$m_{\bar{K}N}$") 
    ax.vlines(params[:threshold][1], minimum(imw), maximum(imw), color="k", linewidth=2, linestyles=:dashed, label=L"$m_{\pi\Sigma}$")  

    l = range(minimum(abstdet), maximum(abstdet), 13)
    ax.contour(rew, imw, abstdet, colors=["#000", "#000"], linestyles="--", levels=l, linewidths=1)
    cf = ax.contourf(rew, imw, abstdet, levels=l, cmap=:jet) #RdYlBu
    ax.legend(title=(rs=="rs21" ? L"RS[+-]" : L"RS[+---]") )
    return cf
end

"""Find the pole"""
pole_LO = function pole_LO!(init_x::Vector, params::Dict, a::Vector; n=2, rs="rs21", μ_tuple::Tuple=(630e0,), ftol=1e-8, xtol=0)
    mch = params[:mch]
    decons = params[:decons_vec]
    sol = nlsolve(x -> cmplx!(tdet(x[1]+x[2]*1e0im, mch, decons, a, n=n, rs=rs, μ_tuple=μ_tuple) ), init_x, ftol=ftol, xtol=xtol)
    if converged(sol) == true
        return sol.zero#@. round(sol.zero, digits=2)
    else
        #println("The routine was not converged.")
        return ["nan", "nan"]#[NaN, NaN]#@. round(sol.zero, digits=2)
    end
end

function pole_LO_RS_2ch!(init_x::Vector, params::Dict, a::Vector; ftol=1e-8, xtol=0)
    RS = ["rs11", "rs21", "rs12", "rs22"]
    p = []
    for rs in RS
        pp = pole_LO!(init_x, params, a, rs=rs, n=2, ftol=ftol, xtol=xtol)
        if pp == ["nan", "nan"]
            pp = [NaN, NaN]
        end
            append!(p, pp[1]+pp[2]*1im)
    end
    df = DataFrame(RS=RS, pole=p)
    return df
end
