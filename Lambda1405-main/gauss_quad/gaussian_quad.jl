

#include("./xw.jl")
#include("./x.jl")
using FastGaussQuadrature

function gaussian_quad!(fcn, a, b; npt::Int=100)
    EPS = 1e-7
    Δ = 1e-2
    last_inte = 0e0
    i::Int = npt

    while abs(Δ) >= EPS
        inte = 0
        xw = gausslegendre(i)#xw[i]
        x, w = xw
        xp, wp = (b - a)/2 .* x .+ (a+b)/2, (b-a)/2 .* w

        for (ii, xxp) in enumerate(xp)
            inte += wp[ii] * fcn(xxp)
        end
        Δ = inte - last_inte
        last_inte = inte
        i += 10
        if length(x) == 900
            println("Reaching the maximum order of Gaussian quadrature")
            break
        end
    end
    return last_inte
end


