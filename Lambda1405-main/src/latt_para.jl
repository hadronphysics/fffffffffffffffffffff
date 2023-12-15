

lattice_params = function lattice_params!()
    hbarc =  197.326980
    L = 4.05 / hbarc # MeV^(-1)
    mN, mΛ, mΣ, mXi = 979.77677431, 1132.83767033, 1193.93733555, 1295.16
    mπ, mK = 203.65515961, 486.36580442
    meta = 551.1
    mch = [(mπ, mΣ), (mK, mN), (meta, mΛ), (mK, mXi)]
    fπ, fK = 131.73835979 / sqrt(2), 153.06089602 / sqrt(2)
    feta = 119.21657697632926
    latt = Dict(
        :mch => mch,
        :decons_vec => [fπ, fK, feta, fK],
        :threshold => [sum(mch[i]) for i in 1:4],
        :m_meson => (mπ, mK),
        :threshold => [sum(mch[i])  for i in 1:4],
        :L => L
    )
    return latt
end