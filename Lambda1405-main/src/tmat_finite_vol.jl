

using LinearAlgebra
using NLsolve

include("./loop_fun.jl")
include("./potential.jl")



function tdet_finite_vol!(w, mch, decons, a, L, N; n=2, μ=630e0, cutoff=2500e0)
    v = VLO!(w, mch, decons, n=n)
    v = real(v)
    id = diagm(ones(n) )
    gl = @fastmath diagm([Gtilde!(w, mchi, a, L, N, μ=μ, cutoff=cutoff) for (mchi, a) in zip(mch[1:n], a[1:n]) ] )

    return REAL!(det(id - v*gl) )
end

function tdet_finite_vol_params!(w, params, a, N; n=2, μ=630e0, cutoff=2500e0)
    mch = params[:mch]
    L = params[:L]
    decons = params[:decons_vec]
    return tdet_finite_vol!(w, mch, decons, a, L, N, n=n, μ=μ, cutoff=cutoff)
end

function plot_tdet_finite_vol!(ax, w, params::Dict, a, N; n=2, μ=630e0, cutoff=2500e0)

    decons = params[:decons_vec]
    mch = params[:mch]
    L = params[:L]
    
    tdet = [tdet_finite_vol!(ww, mch, decons, a, L, N, n=n, μ=μ, cutoff=cutoff) for ww in w]
    ax.hlines(0, minimum(w), maximum(w), color="r")
    

    ax.set_ylim(-0.5, 0.5)
    ax.set_xlim(minimum(w), maximum(w) )
    ax.set_xlabel(L"$E_{cm}$ [MeV]")
    ax.grid()
    ax.scatter(w, tdet, s=12)
end


function pole_finite_vol!(init_x, params, a, N)
    sol = nlsolve(x -> tdet_finite_vol_params!(x[1], params, a, N), [init_x*1e0] )
    if converged(sol) == true
        return sol.zero[1]
    else
        println("The routine was not converged.")
        return sol.zero[1]
    end
end

function pole_finite_vol_all_boost!(init_x, params, a, Nv)
    p = DataFrame([])
    len = 4
    columns = ["[0, 0, 0]", "[0, 0, 1]", "[0,1,1]", "[1,1,1]"]
    for (i, xx, N) in zip(1:4, init_x, Nv)
        pp = (i==4 ? bisection_internal_sol!(x -> tdet_finite_vol_params!(x, params, a, N), 1350, 1550, 0.5).x0 : [pole_finite_vol!(x, params, a, N) for x in xx] )
        if length(pp) != len
            append!(pp, NaN)
        end
        p[!, columns[i] ] = pp
    end
    return p
end

function bisection_pole_finite_vol!(params, a, N)
    x0 = bisection_internal_sol!(x -> tdet_finite_vol_params!(x, params, a, N), 1300, 1700, 1).x0
    x0 = N == [0, 1, 1] ? x0[1:3] : x0[1:4]
    return x0
end


function bisection_internal_sol!(fun, a, b, step)
    x0 = a
    inter = []
    inter_updated = []
    bis_sol_vec = []
    sol_vec::Vector{Any} = []
    for x in a:step:b
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