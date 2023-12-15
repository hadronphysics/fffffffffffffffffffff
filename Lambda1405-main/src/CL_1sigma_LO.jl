

using DataFrames
using CSV


include("./tmat.jl")
include("./tmat_finite_vol.jl")



function CL_1σ_LO!(cost_fun, fitted_para, para_err; nrand=2)
    a1, a2, fϕ, f = [rand(lower:0.01:upper, nrand) for (lower, upper) in zip(fitted_para-para_err, fitted_para+para_err) ]

    UP = 1.5#4.88
    bestredχ2 = cost_fun(fitted_para)
    
    df = DataFrame([[] [] [] [] []], [:a1, :a2, :fϕ, :f, :redχ2])
    CSV.write("./out/CL1sigmaLO.csv", df)
    
    para = zeros((1, 5) )
    for (aa1, aa2, fϕi, fi) in zip(a1, a2, fϕ, f)
        ppi = [aa1, aa2, fϕi, fi]

        redχ2 = cost_fun(ppi)

        if redχ2 <= bestredχ2 + UP
            para[1, 5] = redχ2
            for i in 1:4
                para[1, i] = ppi[i]
            end
            df = DataFrame(para, :auto)
            CSV.write("./out/CL1sigmaLO.csv", df, append=true)
        end
    end
    return DataFrame(CSV.File("./out/CL1sigmaLO.csv") )
end

function Λ1405_1σ!(para1σ, params, latt_pole, best_fitted_para::DataFrame)
    columns = [:lv1, :lv2, :lv3, :lv4]
    Gdf = DataFrame([[] [] [] []], columns)
    [CSV.write(path, Gdf) for path in ["./out/G01sigma.csv", "./out/G11sigma.csv", "./out/G21sigma.csv", "./out/G31sigma.csv"] ]

    a1vec, a2vec, fϕvec, fvec = [], [], [], []
    repH, impH, repL, impL = [], [], [], []
    params_copy = copy(params) # para used in the loop

    # pole with best fitted parameters
    # params is the lattice setup, we change the decay constants using the strategy 1

    best_params = copy(params)

    ba1, ba2, bfϕ, bf = best_fitted_para.para
    best_params[:decons_vec] = [bfϕ, bfϕ] .* bf

    bpH = pole_LO!([1466e0, 10e0], best_params, [ba1, ba2])
    bpL = pole_LO!([1397e0, 40e0], best_params, [ba1, ba2])


    Nv = [[0,0,0], [0,0,1],[0,1,1],[1,1,1]]
    ppf1, ppf2, ppf3, ppf4 = zeros((1, 4) ), zeros((1, 4) ), zeros((1, 4) ), zeros((1, 4) )

    for (a1, a2, fϕ, f) in zip(para1σ.a1, para1σ.a2, para1σ.fϕ, para1σ.f)

        params_copy[:decons_vec] = [fϕ, fϕ] .* f


        pH = pole_LO!([1466e0, 10e0], params_copy, [a1, a2])
        pL = pole_LO!([1397e0, 40e0], params_copy, [a1, a2])
        
        if pH == ["nan", "nan"] || pL == ["nan", "nan"]
            continue
        else
            pH = abs.(pH)
            pL = abs.(pL)
        end

        if 1365 <= pL[1] && pL[1] <= 1419 && 1423 <= pH[1] && pH[1] <= 1487 && pL[2] <= 80
            append!(a1vec, a1)
            append!(a2vec, a2)
            append!(fϕvec, fϕ)
            append!(fvec, f)

            append!(repH, pH[1])
            append!(impH, pH[2])

            append!(repL, pL[1])
            append!(impL, pL[2])

            pf = pole_finite_vol_all_boost!(latt_pole, params_copy, [a1, a2], Nv)
            for i in 1:4
                ppf1[1, i] = pf[i, 1]
                ppf2[1, i] = pf[i, 2]
                ppf3[1, i] = pf[i, 3]
                ppf4[1, i] = pf[i, 4]
            end

            dfv = [DataFrame(ppf, :auto) for ppf in [ppf1, ppf2, ppf3, ppf4] ]

            [CSV.write(path, dfvv, append=true) for (path, dfvv) in zip(["./out/G01sigma.csv", "./out/G11sigma.csv", "./out/G21sigma.csv", "./out/G31sigma.csv"], dfv) ]

        else
            continue
        end
    end
    #impH, impL = abs.(impH), abs.(impL)
    
    repH_upper = maximum(repH) - bpH[1]
    repH_lower = bpH[1] - minimum(repH)
    
    impH_upper = maximum(impH) - bpH[2]
    impH_lower = bpH[2] - minimum(impH)

    repL_upper = maximum(repL) - bpL[1]
    repL_lower = bpL[1] - minimum(repL)
    impL_upper = maximum(impL) - bpL[2]
    impL_lower = bpL[2] - minimum(impL)

    para_df = DataFrame(a1=a1vec, a2=a2vec, fϕ=fϕvec, f=fvec)

    pole_df = DataFrame(repH=repH, impH=impH, repL=repL, impL=impL)
    pole_df_err = DataFrame(pH=bpH, pH_lower_upper=[(repH_lower, repH_upper), (impH_lower, impH_upper)], pL=bpL, pL_lower_upper=[(repL_lower, repL_upper), (impL_lower, impL_upper)])

    CSV.write("./out/resonable_1σ_fitted_params.csv", para_df)
    CSV.write("./out/Lambda1405_1σ.csv", pole_df)
    CSV.write("./out/Lambda1405_error_band.csv", pole_df_err)
    return nothing
end