

include("./formula.jl")


function nor!(w, mchi)
    m, M = mchi
    return sqrt((M + Ecm!(w, M, m)) / (2M) + 0im)
end



function VLO!(w, mch, decons; n=2)
    Cij = [4 -sqrt(3 / 2) 0 sqrt(3 / 2)
        -sqrt(3 / 2) 3 3/sqrt(2) 0
        0 3/sqrt(2) 0 -3/sqrt(2)
        sqrt(3 / 2) 0 -3/sqrt(2) 3]
    f = decons

    v = zeros(ComplexF64, n, n)
    if typeof(f) == Vector{Float64}
        for i in 1:n
            for j in 1:n
                v[i, j] = Cij[i, j] * (1 / f[i] / f[j]) * (2w - mch[i][2] - mch[j][2]) * nor!(w, mch[i]) * nor!(w, mch[j])
            end
        end
    elseif typeof(f) == Matrix{Float64}
        for i in 1:n
            for j in 1:n
                v[i, j] = Cij[i, j] * (1 / f[i, j] / f[i, j]) * (2w - mch[i][2] - mch[j][2]) * nor!(w, mch[i]) * nor!(w, mch[j])
            end
        end
    end
    return -1 / 4 * v
end


"""Potential at next to leading order"""
function VNLO!(w, mπ, mK, mch, b::Dict, d::Dict, decons::Vector; n=2)
    # Dij
    b0, bD, bF = b[:b0], b[:bD], b[:bF]
    μ1 = sqrt(mK^2 + mπ^2)
    μ2 = sqrt(5mK^2 - 3mπ^2)
    μ3 = sqrt(4mK^2 - mπ^2)
    μ4 = sqrt(16mK^2 - 7mπ^2)

    D11 = 4 * (b0 + bD) * mπ^2
    D12 = -sqrt(3 / 2) * (bD - bF) * μ1^2
    D13 = -(4bD * mπ^2) / sqrt(3)
    D14 = sqrt(3 / 2) * (bD + bF) * μ1^2

    D21 = D12
    D22 = 2 * (2b0 + 3bD + bF) * mK^2
    D23 = (bD + 3bF) * μ2^2 / (3sqrt(2))
    D24 = 0e0

    D31, D32 = D13, D23
    D33 = (4 / 9) * (3b0 * μ3^2 + bD * μ4^2)
    D34 = -(bD - 3bF) * μ2^2 / (3sqrt(2))
    
    D41, D42, D43 = D14, D24, D34
    D44 = 2 * (2b0 + 3bD - bF) * mK^2

    D = [D11 D12 D13 D14; D21 D22 D23 D24; D31 D32 D33 D34; D41 D42 D43 D44]

    #=D = [4*(b0+bD)*mπ^2 -sqrt(3 / 2)*(bD-bF)*μ1^2 -(4bD * mπ^2)/sqrt(3) sqrt(3 / 2)*(bD+bF)*μ1^2;
        -sqrt(3 / 2)*(bD-bF)*μ1^2 2*(2b0+3bD+bF)*mK^2 (bD+3bF)*μ2^2/(3sqrt(2)) 0;
        -(4bD * mπ^2)/sqrt(3) (bD+3bF)*μ2^2/(3sqrt(2)) (4/9)*(3b0*μ3^2+bD*μ4^2) -(bD - 3bF)*μ2^2/(3sqrt(2));
        sqrt(3 / 2)*(bD+bF)*μ1^2 0 -(bD - 3bF)*μ2^2/(3sqrt(2)) 2*(2b0+3bD-bF)*mK^2]=#

    # Lij
    d1, d2, d3, d4 = d[:d1], d[:d2], d[:d3], d[:d4]

    L11 = -4d2+4d3+2d4
    L12 = sqrt(3 / 2)*(d1+d2-2d3)
    L13 = -sqrt(3)d3
    L14 = sqrt(3 / 2)*(d1-d2+2d3)

    L21= L12
    L22 = d1+3d2+2*(d3+d4)
    L23 = (d1-3d2+2d3)/sqrt(2)
    L24 = (6d2-3d3)

    L31, L32 = L13, L23
    L33 = 2 * (d3 + d4)
    L34 = (d1+3d2-2d3)/sqrt(2)

    L41, L42, L43 = L14, L24, L34
    L44 = -d1 + 3d2 + 2 * (d3 + d4)

    L = [L11 L12 L13 L14; L21 L22 L23 L24; L31 L32 L33 L34; L41 L42 L43 L44]
    #=
    L = [-4d2+4d3+2d4 sqrt(3 / 2)*(d1+d2-2d3) -sqrt(3)d3 sqrt(3 / 2)*(d1-d2+2d3);
        sqrt(3 / 2)*(d1+d2-2d3) d1+3d2+2*(d3+d4) (d1-3d2+2d3)/sqrt(2) (6d2-3d3);
        -sqrt(3)*d3 (d1-3d2+2d3)/sqrt(2) 2*(d3+d4) (d1+3d2-2d3)/sqrt(2);
        sqrt(3 / 2)*(d1-d2+2d3) 6d2-3d3 (d1+3d2-2d3)/sqrt(2) -d1+3d2+2*(d3+d4)]
    =#

    v = zeros(ComplexF64, n, n)
    for i in 1:n
        for j in 1:n

            ωi = ωω!(w, mch[i]) # energy of initial meson sqrt(qcm!(w, mch[i]...)^2 + mch[i][1]^2)#
            ωj = ωω!(w, mch[j]) # energy of finial meson sqrt(qcm!(w, mch[j]...)^2 + mch[j][1]^2)#

            qi_sq = qcm!(w, mch[i]...)^2
            qj_sq = qcm!(w, mch[j]...)^2

            EBi = EB!(w, mch[i])
            EBj = EB!(w, mch[j])

            Mi = mch[i][2]
            Mj = mch[j][2]

            v[i, j] = (D[i, j] - 2 * L[i, j] * (ωi * ωj + qi_sq * qj_sq / (3 * (Mi + EBi) * (Mj + EBj)) ) ) * nor!(w, mch[i]) * nor!(w, mch[j]) * (1 / (decons[i] * decons[j]))
        end
    end

    return v
end

"""Potential at next to leading order"""
function VNLO_old!(w, mπ, mK, mch, b::Dict, d::Dict, decons::Vector; n=2)
    # Dij
    b0, bD, bF = b[:b0], b[:bD], b[:bF]
    μ1 = sqrt(mK^2 + mπ^2)
    μ2 = sqrt(5mK^2 - 3mπ^2)
    μ3 = sqrt(4mK^2 - mπ^2)
    μ4 = sqrt(16mK^2 - 7mπ^2)

    D = [4*(b0+bD)*mπ^2 -sqrt(3 / 2)*(bD-bF)*μ1^2 -(4bD * mπ^2)/sqrt(3) sqrt(3 / 2)*(bD+bF)*μ1^2;
        -sqrt(3 / 2)*(bD-bF)*μ1^2 2*(2b0+3bD+bF)*mK^2 (bD+3bF)*μ2^2/(3sqrt(2)) 0;
        -(4bD * mπ^2)/sqrt(3) (bD+3bF)*μ2^2/(3sqrt(2)) (4/9)*(3b0*μ3^2+bD*μ4^2) -(bD - 3bF)*μ2^2/(3sqrt(2));
        sqrt(3 / 2)*(bD+bF)*μ1^2 0 -(bD - 3bF)*μ2^2/(3sqrt(2)) 2*(2b0+3bD-bF)*mK^2]

    # Lij
    d1, d2, d3, d4 = d[:d1], d[:d2], d[:d3], d[:d4]
    L = [-4d2+4d3+2d4 sqrt(3 / 2)*(d1+d2-2d3) -sqrt(3)d3 sqrt(3 / 2)*(d1-d2+2d3);
        sqrt(3 / 2)*(d1+d2-2d3) d1+3d2+2*(d3+d4) (d1-3d2+2d3)/sqrt(2) (6d2-3d3);
        -sqrt(3)*d3 (d1-3d2+2d3)/sqrt(2) 2*(d3+d4) (d1+3d2-2d3)/sqrt(2);
        sqrt(3 / 2)*(d1-d2+2d3) 6d2-3d3 (d1+3d2-2d3)/sqrt(2) -d1+3d2+2*(d3+d4)]

    v = zeros(ComplexF64, n, n)
    for i in 1:n
        for j in 1:n

            ωi = ωω!(w, mch[i]) # energy of initial meson sqrt(qcm!(w, mch[i]...)^2 + mch[i][1]^2)#
            ωj = ωω!(w, mch[j]) # energy of finial meson sqrt(qcm!(w, mch[j]...)^2 + mch[j][1]^2)#

            qi_sq = qcm!(w, mch[i]...)^2
            qj_sq = qcm!(w, mch[j]...)^2

            EBi = EB!(w, mch[i])
            EBj = EB!(w, mch[j])

            Mi = mch[i][2]
            Mj = mch[j][2]

            v[i, j] = (D[i, j] - 2 * L[i, j] * (ωi * ωj)) * nor!(w, mch[i]) * nor!(w, mch[j]) * (1 / (decons[i] * decons[j]))
        end
    end

    return v
end