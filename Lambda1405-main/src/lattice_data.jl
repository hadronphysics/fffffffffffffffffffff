

using CSV
using DataFrames



_2lv0 = [6.75, 7.104761904761905] .* 203.65515961
_2lv0_err = [6.783333333333333, 7.145238095238096].* 203.65515961

_2lv1 = [6.783333333333333, 7.130952380952381].* 203.65515961
_2lv1_err = [6.81904761904762, 7.192857142857143].* 203.65515961

_2lv2 = [6.823809523809524, 7.121428571428572].* 203.65515961
_2lv2_err = [6.859523809523809, 7.180952380952381].* 203.65515961

_2lv3 = [6.814285714285714, 7.147619047619048].* 203.65515961
_2lv3_err = [6.85, 7.204761904761905].* 203.65515961

_2lv0_err = _2lv0_err - _2lv0
_2lv1_err = _2lv1_err - _2lv1
_2lv2_err = _2lv2_err - _2lv2
_2lv3_err = _2lv3_err - _2lv3


_2lv = [_2lv0, _2lv1, _2lv2, _2lv3]
_2lv_err = [_2lv0_err, _2lv1_err, _2lv2_err, _2lv3_err]

# All data
p0 = [6.75, 7.104761904761905, 7.340506329113924, 7.689873417721519] .* 203.65515961
p0_err = [6.783333333333333, 7.145238095238096, 7.40126582278481, 7.7449367088607595].* 203.65515961

p1 = [6.783333333333333, 7.130952380952381, 7.368987341772152, 7.4525316455696204].* 203.65515961
p1_err = [6.81904761904762, 7.192857142857143, 7.422151898734177, 7.507594936708861].* 203.65515961

p2 = [6.823809523809524, 7.121428571428572, 7.3955696202531644].* 203.65515961
p2_err = [6.859523809523809, 7.180952380952381, 7.4525316455696204].* 203.65515961

p3 = [6.814285714285714, 7.147619047619048, 7.441139240506329, 7.492405063291139].* 203.65515961
p3_err = [6.85, 7.204761904761905, 7.498101265822784, 7.55126582278481].* 203.65515961

p0_err = p0_err - p0
p1_err = p1_err - p1
p2_err = p2_err - p2
p3_err = p3_err - p3

latt_pole = [p0, p1, p2, p3]
latt_pole_err = [p0_err, p1_err, p2_err, p3_err]


function plot_latt_ener_2level!(ax, latt_pole, latt_pole_err, threshold)
    
    [ax.errorbar([i, i], latt_pole[i], yerr=latt_pole_err[i], fmt="s", c="g") for i in 1:4]
    #ax.errorbar([2, 2, 2], latt_pole[2], yerr=latt_pole_err[2], fmt="s", c="g")
    ax.hlines(threshold[1], -0.1, 4.2, label=L"$m_{\pi\Sigma}$", linestyle="--")
    ax.hlines(threshold[2], -0.1, 4.2, label=L"$m_{\bar{K}N}$")
    ax.legend()
    ax.set(xlim=(0.9, 4.2), xlabel="boost", ylabel=L"$E_{cm}$ [MeV]")
    ax.set_xticks([1, 2, 3, 4], [L"$[0,0,0]$", L"$[0,0,1]$", L"$[0,1,1]$", L"$[1,1,1]$"])
    
end

function plot_latt_ener!(ax, latt_pole, latt_pole_err, threshold)
    
    [ax.errorbar([i, i, i, i], latt_pole[i], yerr=latt_pole_err[i], fmt="o", c="g", label=(i==1 ? "Lattice data" : false)) for i in [1, 2, 4]]
    ax.errorbar([3, 3, 3], latt_pole[3], yerr=latt_pole_err[3], fmt="o", c="g")
    
    #ax.errorbar([2, 2, 2], latt_pole[2], yerr=latt_pole_err[2], fmt="s", c="g")
    ax.hlines(threshold[1], -0.1, 4.2, label=L"$m_{\pi\Sigma}$", linestyle="--")
    ax.hlines(threshold[2], -0.1, 4.2, label=L"$m_{\bar{K}N}$")
    ax.legend()
    ax.set(xlim=(0.9, 4.2), xlabel="Boost", ylabel=L"$E_{cm}$ [MeV]")
    ax.set_xticks([1, 2, 3, 4], [L"$[0,0,0]$", L"$[0,0,1]$", L"$[0,1,1]$", L"$[1,1,1]$"])
    
end

function plot_theory_ener_2level!(ax, pole)
    ds = 0.1
    [ax.scatter([i, i] .+ ds, pole[i], color="#e6091c", marker="^", label=(i==1 ? "Lattice data" : false)) for i in 1:4]
end


function plot_theory_ener!(ax, pole)
    ds = 0.1
    [ax.scatter([i, i, i, i] .+ ds, pole[i], color="#e6091c", marker="^") for i in 1:4]
end

path = ["./out/G0lo_errbar.csv", "./out/G1lo_errbar.csv", "./out/G2lo_errbar.csv", "./out/G3lo_errbar.csv"]

function plot_lv_errbar!(ax, path::Vector=path)
    ds = 0.1
    
    for (i, p) in enumerate(path)
        ii = i==3 ? [i, i, i] : [i, i, i, i]
        df = DataFrame(CSV.File(p) )
        ax.errorbar(ii .+ ds, df.pole, yerr=[df.lower, df.upper], color="#e6091c", fmt="s", label=(i==1 ? "Fitted energy level" : false))
    end
    ax.legend()
end