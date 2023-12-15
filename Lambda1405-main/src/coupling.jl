

using QuadGK

include("./tmat.jl")
include("../gauss_quad/gaussian_quad.jl")



function residue_LO!(pole, params, a; n=2, ch="ch11")
    w = pole[1] - pole[2] * 1im
    mch = params[:mch]
    decons = params[:decons_vec]
    EPS = 1e-5
    d = 1
    last_res = 0e0
    delta = 1e-1
    while abs(delta) > EPS
        d *= 0.9
        res = abs((d - d * 1im) * tmat!(w + d - d * 1im, mch, decons, a, n=n, ch=ch))
        delta = res - last_res
        last_res = res
    end
    #d = 1e-9#1 / abs(tmat!(w, mch, decons, a,n=n, ch=ch) )
    #res = (d-d*1im) * tmat!(w+d-d*1im, mch, decons, a, n=n, ch=ch)
    return last_res, d
end

function coupling_LO!(pole, params, a; n=2, ch="ch11", return_r=false)
    w = pole[1] - pole[2] * 1im
    mch = params[:mch]
    decons = params[:decons_vec]

    EPS = 1e-8
    delta = 1e-1
    last_inte = 0e0
    r = 1
    while abs(delta) > EPS
        r *= 0.1
        #inte = (quadgk(x -> real(tmat!(w + r * exp(x * 1im), mch, decons, a, n=n, ch=ch) * exp(x * 1im) ), 0, 2π)[1] + quadgk((x -> imag(tmat!(w + r * exp(x * 1im), mch, decons, a, n=n, ch=ch) * exp(x * 1im) ), 0, 2π) ) * 1im )[1] * (r / (2π) )
        inte = @fastmath gaussian_quad!(x -> tmat!(w + r * exp(x * 1im), mch, decons, a, n=n, ch=ch) * exp(x * 1im), 0, 2π) * (r / (2π))
        delta = last_inte - abs(inte)
        last_inte = abs(inte)
    end
    if return_r == true
        return last_inte, r
    else
        return last_inte
    end
end

"""Calculation of coupling square for all channels"""
function df_coupling_LO!(pole, params, a; n=2)
    chv = n == 2 ? ["ch11", "ch22"] : ["ch11", "ch22", "ch33", "ch44"]
    columns = n == 2 ? [:gs11, :gs22] : [:gs11, :gs22, :gs33, :gs44]
    gs = zeros((1, n))
    for (i, ch) in enumerate(chv)
        gs[1, i] = coupling_LO!(pole, params, a, n=n, ch=ch)
    end
    df = DataFrame(gs, columns)
    return df
end


function coupling_1sigma_LO!(params, fitted_para_1sigma; n=2, a34=[-2e0, -2e0])

    if n == 2
        path = ["./out/coupling/pH_coupling_1sigma_2ch_LO.csv", "./out/coupling/pL_coupling_1sigma_2ch_LO.csv"]
        dfpH = DataFrame([[] [] [] []], [:repH, :impH, :gsqr11, :gsqr22])
        dfpL = DataFrame([[] [] [] []], [:repL, :impL, :gsqr11, :gsqr22])
        CSV.write(path[1], dfpH)
        CSV.write(path[2], dfpL)
        gspLmat = zeros((1, 4))
        gspHmat = zeros((1, 4))
    elseif n == 4
        path = ["./out/coupling/pH_coupling_1sigma_4ch_LO.csv", "./out/coupling/pL_coupling_1sigma_4ch_LO.csv"]
        dfpH = DataFrame([[] [] [] [] [] []], [:repH, :impH, :gsqr11, :gsqr22, :gsqr33, :gsqr44])
        dfpL = DataFrame([[] [] [] [] [] []], [:repL, :impL, :gsqr11, :gsqr22, :gsqr33, :gsqr44])
        CSV.write(path[1], dfpH)
        CSV.write(path[2], dfpL)
        gspLmat = zeros((1, 6))
        gspHmat = zeros((1, 6))
    end
    params_copy = copy(params)
    for (a1, a2, fϕ, f) in zip(fitted_para_1sigma.a1, fitted_para_1sigma.a2, fitted_para_1sigma.fϕ, fitted_para_1sigma.f)
        params_copy[:decons_vec] = [fϕ, fϕ, fϕ, fϕ] .* f
        #println(par)
        a = [a1, a2, a34...]

        pH = pole_LO!([1466e0, 10e0], params_copy, a, n=n)
        pL = pole_LO!([1390e0, 40e0], params_copy, a, n=n)
        if pH == ["nan", "nan"] || pL == ["nan", "nan"] || pL[1] < 1365 || pH[1] < 1423
            continue
        else
            pH = abs.(pH)
            pL = abs.(pL)
        end
        gsL = df_coupling_LO!(pL, params_copy, a, n=n)
        gsH = df_coupling_LO!(pH, params_copy, a, n=n)
        gsL, gsH = ([gsL[1, :]...], [gsH[1, :]...])
        #println([eachrow(gsL)...])
        condi = 30
        if all(gsL .<= condi) && all(gsH .<= condi)
            gsL, gsH = sqrt.(gsL), sqrt.(gsH)
            [gspLmat[1, i] = val for (i, val) in enumerate(append!([], pL, gsL))]
            [gspHmat[1, i] = val for (i, val) in enumerate(append!([], pH, gsH))]

            dfpL = DataFrame(gspLmat, :auto)
            dfpH = DataFrame(gspHmat, :auto)
            CSV.write(path[1], dfpH, append=true)
            CSV.write(path[2], dfpL, append=true)
        end


    end

    return nothing
end