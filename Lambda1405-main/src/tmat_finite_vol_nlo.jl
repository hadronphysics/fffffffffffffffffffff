

using LinearAlgebra
using NLsolve
using DataFrames


include("./potential.jl")
include("./loop_fun.jl")


function tdet_finite_vol_nlo!(w, params::Dict, a, b, d, N; cutoff=2500e0, μ=630e0)
    mch = params[:mch]
    decons = params[:decons_vec]
    mπ, mK = params[:m_meson]
    L = params[:L]

    vlo = VLO!(w, mch, decons)
    vnlo = VNLO!(w, mπ, mK, mch, b, d, decons)
    v = real(vlo + vnlo)
    n = 2
    id = diagm(ones(n) )
    gl = @fastmath diagm([Gtilde!(w, mchi, a, L, N, μ=μ, cutoff=cutoff) for (mchi, a) in zip(mch[1:n], a[1:n]) ] )

    return REAL!(det(id - v*gl) )
end
#=
function tdet_finite_vol_nlo!(w, mch, decons, m_meson, L, a, b, d, N; cutoff=2500e0, μ=630e0)
    mπ, mK = m_meson
    vlo = VLO!(w, mch, decons)
    vnlo = VNLO!(w, mπ, mK, mch, b, d, decons)
    v = real(vlo + vnlo)

    n = 2
    id = diagm(ones(n) )
    gl = @fastmath diagm([Gtilde!(w, mchi, a, L, N, μ=μ, cutoff=cutoff) for (mchi, a) in zip(mch[1:n], a[1:n]) ] )

    return REAL!(det(id - v*gl) )
end
=#

function plot_tdet_finite_vol_nlo!(ax, w, params, a, b, d, N)
    tdet = [tdet_finite_vol_nlo!(ww, params, a, b, d, N) for ww in w]
    ax.hlines(0, minimum(w), maximum(w), color="r")
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlim(minimum(w), maximum(w) )
    ax.set_xlabel(L"$E_{cm}$ [MeV]")
    ax.grid()
    ax.scatter(w, tdet, s=12)
end

function bisection_pole_finite_vol_nlo!(params, a, b, d, N)
    x0 = bisection_internal_sol_nlo!(x -> tdet_finite_vol_nlo!(x, params, a, b, d, N), 1300, 1700, 1).x0
    x0 = N == [0, 1, 1] ? x0[1:3] : x0[1:4]
    return x0
end


function bisection_internal_sol_nlo!(fun, xx0, xx1, step)
    x0 = xx0
    inter = []
    inter_updated = []
    bis_sol_vec = []
    sol_vec::Vector{Any} = []
    for x in xx0:step:xx1
        xx = @fastmath fun(x0) * fun(x)
        if xx < 0
            append!(inter, [(x0, x)] )
            x0 = x
        else
            x0 = x
        end
    end
    for inter in inter
        sol = @fastmath find_zero(fun, inter, Bisection() )
        if abs(fun(sol) ) < 1e-7
            append!(inter_updated, [inter] )
            append!(bis_sol_vec, sol)
            sl = nlsolve(x -> fun(x[1]), [sol*1e0])
            if converged(sl) == true
                append!(sol_vec, sl.zero)
            else
                println("The routine was not converged.")
                append!(sol_vec, sl.zero)
            end
        end
    end
    return DataFrame(bis_inter=inter_updated, bis_x0=bis_sol_vec, x0=sol_vec)
end

function pole_finite_vol_nlo!(init_x, params, a, b, d, N)
    #mch, decons, m_meson, L = params[:mch], params[:decons_vec], params[:m_meson], params[:L]
    
    sol = nlsolve(x -> tdet_finite_vol_nlo!(x[1], params, a, b, d, N), [init_x*1e0] )
    if converged(sol) == true
        return sol.zero[1]
    else
        println("The routine was not converged.")
        return sol.zero[1]
    end
end