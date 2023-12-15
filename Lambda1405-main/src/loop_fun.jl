

using QuadGK

include("./formula.jl")
include("./formula1.jl")


# Propagator at infinite volume
"""DR scheme"""
function Gdr!(w, mch_i, a; μ=630e0)
    d = mch_i[2]^2 - mch_i[1]^2
    s = w^2
    q = qcm!(w, mch_i[1], mch_i[2])

    g = 2 / (16 * π^2) * mch_i[2] * (a + 2log(mch_i[2] / μ) + (-d + s) / (2s) * 2log(mch_i[1] / mch_i[2]) + q / w * (log(s - d + 2q * w) + log(s + d + 2q * w) - log(-s + d + 2q * w) - log(-s - d + 2q * w)))

    #if real(w) > sum(mch_i)
    #    return g + 2im * mch_i[2] * qcm!(w, mch_i...) / (4π*w)
    #else
    return g
    #end

end


"""analysis formula"""
function Gcut!(w, mch_i; cutoff=2500e0)
    qmax = cutoff
    s = w * w
    m1, m2 = mch_i # m1: mass of meson, m2: mass of baryon

    if qmax == 0
        solu = 0
    else
        solu = (-(((-m1^2 + m2^2) * log(complex(m1^2 / m2^2))) / s) + (2 * (-m1^2 + m2^2) * (log(complex(1 + xsqrt(complex(1 + m1^2 / qmax^2)))) - log(complex(1 + xsqrt(complex(1 + m2^2 / qmax^2)))))) / s - 2 * log(complex((1 + xsqrt(complex(1 + m1^2 / qmax^2)))) * (1 + xsqrt(complex(1 + m2^2 / qmax^2)))) + log(complex((m1^2 * m2^2) / qmax^4)) + (xsqrt(complex((-(m1 - m2)^2 + s) * (-(m1 + m2)^2 + s))) * (-log(complex(-m1^2 + m2^2 - s + xsqrt(complex(1 + m1^2 / qmax^2)) * xsqrt(complex((-(m1 - m2)^2 + s) * (-(m1 + m2)^2 + s))))) + log(complex(m1^2 - m2^2 + s + xsqrt(complex(1 + m1^2 / qmax^2)) * xsqrt(complex((-(m1 - m2)^2 + s) * (-(m1 + m2)^2 + s))))) - log(complex(m1^2 - m2^2 - s + xsqrt(complex(1 + m2^2 / qmax^2)) * xsqrt(complex((-(m1 - m2)^2 + s) * (-(m1 + m2)^2 + s))))) + log(complex(-m1^2 + m2^2 + s + xsqrt(complex(1 + m2^2 / qmax^2)) * xsqrt(complex((-(m1 - m2)^2 + s) * (-(m1 + m2)^2 + s))))))) / s) / (32.0 * π^2)
    end

    g = solu * 2m2
    return g
end

ω = function ω!(m, q)
    return sqrt(m^2 + q^2)
end
#=
I = function I!(w, m, M, q, EPS=1e-8)
    Em = ω(m, q)
    EM = ω(M, q)

    return (Em + EM) / (2 * Em * EM) * (1 / (w^2 - (Em + EM)^2 + EPS * 1im))
end
=#

"""Form factor"""
function form_factor!(qmax, q)
    return (1300 + qmax)^8 / ((1300 + qmax)^8 + q^9)
end

"""Function I"""
I = function I!(w, m, M, q; qmax=2500e0, EPS=0)
    E_meson = ω!(m, q)
    E_baryon = ω!(M, q)
    return 1 / (2 * E_meson * E_baryon ) * (1 / (w - E_meson - E_baryon + EPS*1im) ) * form_factor!(qmax, q)
end

function quad_integrate!(fun, a, b)
    inte = quadgk(fun, a, b)
    return inte[1]
end

function Gcut_integral!(w, mchi; qmax=630e0, EPS=10)
    m, M = mchi
    return (M)/(2π)^2 * quad_integrate!(q -> I!(w, m, M, q, qmax=qmax, EPS=EPS) * q^2, 0, qmax)
end


function Gtilde!(w, mchi::Tuple, a, L, N; cutoff=2500e0, μ=630e0)
    m, M = mchi

    qmax = cutoff

    g = 0e0

    nmax::Int = 0
    i::Int = 0

    (qst1_sq, qst2_sq, qst3_sq) = (0e0, 0e0, 0e0)

    while qst1_sq < qmax^2 && qst2_sq < qmax^2 && qst3_sq < qmax^2
        i += 1
        qstar1 = qstar(w, m, M, L, [i, 0, 0], N)
        qstar2 = qstar(w, m, M, L, [0, 0, i], N)
        qstar3 = qstar(w, m, M, L, [0, i, 0], N)

        (qst1_sq, qst2_sq, qst3_sq) = (qstar1' * qstar1, qstar2' * qstar2, qstar3' * qstar3)
    end
    nmax = i

    for n1 in 0:nmax
        for n2 in 1:nmax
            for n3 in 1:nmax
                nn = [n1, n2, n3]
                qst = qstar(w, m, M, L, nn, N)
                qst_sq = qst' * qst
                if qst_sq < qmax^2
                    g += I(w, m, M, sqrt(qst_sq)) * 8
                end
            end
        end
    end

    for n1 in -nmax:nmax
        for n2 in -nmax:nmax
            nn = [n1, n2, 0]
            qst = qstar(w, m, M, L, nn, N)
            qst_sq = qst' * qst

            if qst_sq < qmax^2
                g += 2 * I(w, m, M, sqrt(qst_sq))
            end
        end
    end

    for n3 in 1:nmax
        qst = qstar(w, m, M, L, [0, 0, n3], N)
        qst_sq = qst' * qst
        if qst_sq < qmax^2
            g += 2 * I(w, m, M, sqrt(qst_sq))
        end
    end

    gtilde = (M) / L^3 * (w / P0(w, L, N)) * g

    return REAL!(real(Gdr!(w, mchi, a, μ=μ) + (gtilde - Gcut!(w, mchi, cutoff=cutoff))) )
end


function Gsummixed!(w, mchi, l1, l2, mm1, mm2, L, N; cutoff=2500e0)
    m, M = mchi
    qmax = cutoff

    g = 0e0

    nmax::Int = 0
    i::Int = 0

    (qst1_sq, qst2_sq, qst3_sq) = (0e0, 0e0, 0e0)

    while qst1_sq < qmax^2 && qst2_sq < qmax^2 && qst3_sq < qmax^2
        i += 1
        qstar1 = qstar(w, m, M, L, [i, 0, 0], N)
        qstar2 = qstar(w, m, M, L, [0, 0, i], N)
        qstar3 = qstar(w, m, M, L, [0, i, 0], N)

        (qst1_sq, qst2_sq, qst3_sq) = (qstar1' * qstar1, qstar2' * qstar2, qstar3' * qstar3)
    end
    nmax = i

    for n1 in -nmax:nmax
        for n2 in -nmax:nmax
            for n3 in -nmax:nmax
                nn = [n1, n2, n3]
                qst = qstar!(w, m, M, L, nn, N)
                qst_sq = qst' * qst
                if qst_sq <= qmax^2
                    qst_norm = @. 1 / sqrt(qst_sq) * qst
                    theta = acos(qst_norm[3])
                    phi = atan(qst_norm[2], qst_norm[1])
                    if all(N .== 0)
                        if any(qst .!= 0)
                            angl = ((4π) * conj(sph_harm!(l1, mm1, phi, theta)) * sph_harm!(l2, mm2, phi, theta))
                            g += angl * qqon!(w, m, M, l1, l2, L, nn, N) * I!(w, m, M, sqrt(qst_sq))
                        elseif l1 == 0 && l2 == 0
                            angl = 1e0
                            g += angl * qqon!(w, m, M, l1, l2, L, nn, N) * I!(w, m, M, sqrt(qst_sq))
                        end
                    else
                        angl = ((4π) * conj(sph_harm!(l1, mm1, phi, theta)) * sph_harm!(l2, mm2, phi, theta))
                        g += angl * qqon!(w, m, M, l1, l2, L, nn, N) * I!(w, m, M, sqrt(qst_sq))
                    end
                end
            end
        end
    end
    return (1 / L^3) * g * (w / P0!(w, L, N)) * (2M)
end

function Gtildesummixed!(w, mchi, l1, l2, mm1, mm2, a, L, N; cutoff=2500e0, μ=630e0)
    gtilde = Gdr!(w, mchi, a, μ=μ) + (Gsummixed!(w, mchi, l1, l2, mm1, mm2, L, N; cutoff=cutoff) - Gcut!(w, mchi, cutoff=cutoff))
    return gtilde
end

function Gmat!(w, mch, a, L, N; cutoff=2500e0, μ=630e0, n=2)
    mch1, mch2 = mch[1:n]
    a1, a2 = a[1:n]

    gmat = zeros(ComplexF64, (8, 8))

    g1i = [Gtildesummixed!(w, mch1, 0, 0, 0, 0, a1, L, N), 0,
        Gtildesummixed!(w, mch1, 0, 1, 0, -1, a1, L, N), 0,
        Gtildesummixed!(w, mch1, 0, 1, 0, 0, a1, L, N), 0,
        Gtildesummixed!(w, mch1, 0, 1, 0, 1, a1, L, N), 0]

    g2i = [0, Gtildesummixed!(w, mch2, 0, 0, 0, 0, a2, L, N),
        0, Gtildesummixed!(w, mch2, 0, 1, 0, -1, a2, L, N),
        0, Gtildesummixed!(w, mch2, 0, 1, 0, 0, a2, L, N),
        0, Gtildesummixed!(w, mch2, 0, 1, 0, 1, a2, L, N)]
    
    g3i = [Gtildesummixed!(w, mch1, 1, 0, -1, 0, a1, L, N), 0,
    Gtildesummixed!(w, mch1, 1, 1, -1, -1, a1, L, N), 0,
    Gtildesummixed!(w, mch1, 1, 1, -1, 0, a1, L, N), 0,
    Gtildesummixed!(w, mch1, 1, 1, -1, 1, a1, L, N), 0]

    g4i = [0, Gtildesummixed!(w, mch2, 1, 0, -1, 0, a2, L, N),
        0, Gtildesummixed!(w, mch2, 1, 1, -1, -1, a2, L, N),
        0, Gtildesummixed!(w, mch2, 1, 1, -1, 0, a2, L, N),
        0, Gtildesummixed!(w, mch2, 1, 1, -1, 1, a2, L, N)]

    g5i = [Gtildesummixed!(w, mch1, 1, 0, 0, 0, a1, L, N), 0,
    Gtildesummixed!(w, mch1, 1, 1, 0, -1, a1, L, N), 0,
    Gtildesummixed!(w, mch1, 1, 1, 0, 0, a1, L, N), 0,
    Gtildesummixed!(w, mch1, 1, 1, 0, 1, a1, L, N), 0]

    g6i = [0, Gtildesummixed!(w, mch2, 1, 0, 0, 0, a2, L, N),
        0, Gtildesummixed!(w, mch2, 1, 1, 0, -1, a2, L, N),
        0, Gtildesummixed!(w, mch2, 1, 1, 0, 0, a2, L, N),
        0, Gtildesummixed!(w, mch2, 1, 1, 0, 1, a2, L, N)]

    g7i = [Gtildesummixed!(w, mch1, 1, 0, 1, 0, a1, L, N), 0,
    Gtildesummixed!(w, mch1, 1, 1, 1, -1, a1, L, N), 0,
    Gtildesummixed!(w, mch1, 1, 1, 1, 0, a1, L, N), 0,
    Gtildesummixed!(w, mch1, 1, 1, 1, 1, a1, L, N), 0]

    g8i = [0, Gtildesummixed!(w, mch2, 1, 0, 1, 0, a2, L, N),
        0, Gtildesummixed!(w, mch2, 1, 1, 1, -1, a2, L, N),
        0, Gtildesummixed!(w, mch2, 1, 1, 1, 0, a2, L, N),
        0, Gtildesummixed!(w, mch2, 1, 1, 1, 1, a2, L, N)]

    g = [g1i, g2i , g3i ,g4i, g5i, g6i, g7i, g8i]
    for (i, g) in enumerate(g)
        [gmat[i, j] = cmplx_num!(g[j] ) for j in 1:8]
    end
    return real(gmat)
end
