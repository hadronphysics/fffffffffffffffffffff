

using LegendrePolynomials

include("./formula.jl")
# 🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟 Formula in finite volume 🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟

"""Total momentum of two-particle system"""
P = 
function P!(L, N::Vector)
    return (2π) / L * N
end

P0 = 
function P0!(w, L, N)
    s = w*w
    p = P(L, N)
    PP = p' * p
    return sqrt(s + PP)
end


q = function q!(L, n::Vector)
    return (2π) / L * n
end

"""Energy of particle i in center mass system. For particle j, change the position of them"""
qstar0 = function qstar0!(w, mi, mj)
    s = w*w
    return (s + mi^2 - mj^2) / (2w)
end

"""Three momentum at rest momentum"""
qstar = function qstar!(w, mi, mj, L, n::Vector, N::Vector{Int})
    qP = q(L, n)' * P(L, N)
    PP = P(L, N)' * P(L, N)

    if N == [0, 0, 0]
        return q(L, n)
    else
        return q(L, n) + ((w/P0(w, L, N) -1) * qP/PP - qstar0(w, mi, mj) / P0(w, L, N) ) * P(L, N)
    end
end

function qqon!(w, m1, m2, l1, l2, L, n, N)
    lt = l1 + l2
    if rem(lt, 2) == 0
        solu = 1e0
    else
        PP = P!(L, N)' * P!(L, N)
        ss = (sqrt(m1^2 + 4*π^2/L^2 * (n' * n) )+ sqrt(m2^4 + 4*π^2/L^2*((N-n)' * (N-n) )) )^2 - PP
        qon = qcm!(sqrt(ss), m1, m2)

        if abs(qon) == 0.0
            solu = 1e0
        else
            qst = qstar!(w, m1, m2, L, n, N)
            qst_sq = qst' * qst
            solu = sqrt(qst_sq) / qon
        end
    end
    return solu
end

function sph_harm!(l, m, ϕ, θ)
    if m>=0
        return Plm(cos(θ), l, m) * exp(1im*m*ϕ)
    else
        return ((-1)^m * Plm(cos(θ), l, -m) ) * exp(1im*m*ϕ)
    end
end

