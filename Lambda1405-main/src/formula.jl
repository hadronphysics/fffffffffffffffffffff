"""redefine sqrt so that its cut is along the positive x axis"""
function xsqrt(x)
    imag(x) >=0 ? sqrt(x+0im) : -sqrt(x-0im)
end

function kallen!(x, y, z)
    return x^2 + y^2 + z^2 - 2x*y - 2y*z - 2z*x
end

"""memuntum in center mass"""
function qcm!(w, m1, m2)
    return xsqrt(kallen!(w^2, m1^2, m2^2)+0im) / (2w)
end

"""Energy in center mass"""
function Ecm!(w, m1, m2)
    return (w^2 + m1^2 - m2^2) / (2w)
end

function ωω!(w, mchi)
    return Ecm!(w, mchi...)
end

function EB!(w, mchi)
    m, M = mchi
    return Ecm!(w, M, m)
end

function cmplx!(xx)
    EPS = 1e-9
    re, im = real(xx), imag(xx)

    if abs(re) <= EPS && abs(im) <= EPS
        return (0e0, 0e0)
    elseif abs(im) <= EPS
        return (re, 0e0)
    else
        return (re, im)
    end
end

function cmplx_num!(xx)
    EPS = 1e-9
    re, im = real(xx), imag(xx)

    if abs(re) <= EPS && abs(im) <= EPS
        return 0e0 + 0e0im
    elseif abs(im) <= EPS
        return re + 0e0im
    else
        return re + im * 1.0im
    end
end

function REAL!(xx)
    EPS = 1e-10
    if abs(xx) <= EPS
        return 0e0
    else
        return xx
    end
end