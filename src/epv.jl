module epv

push!(LOAD_PATH, @__DIR__)
#include("param.jl")
#include("commonblocks.jl")
using .param: EP, EE, CSPEED
using .commonblocks: specie


# Quantities:
# beta - speed in units of c (speed of light)
# gamma - lorentz factor
# v - velocity
# rp - ??? some sort of mass quantity (in GeV/c^2?)
# e -
#
# rnm - ???
# rnz - ???

function e2p(e::Float64, spe::specie)::Float64
    if spe.rnm > 0.5
        return sqrt(2.0*EP*e + e*e) * rnm / rnz
    else
        return sqrt(2.0*EE*e + e*e)
    end
end

function rp2e(rp::Float64, spe::specie)::Float64
    if spe.rnm > 0.5
        rp1 = rp*rnz/rnm
        return hypot(EP, rp1) - EP
    else
        rp1 = rp
        return hypot(EE, rp1) - EE
    end
end

function rp2beta(rp::Float64, spe::specie)::Float64
    if spe.rnm > 0.5
        rp1 = rp * rnz / rnm
        return rp1 / hypot(EP, rp1)
    else
        rp1 = rp
        return rp1 / hypot(EE, rp1)
    end
end

function e2beta(e::Float64, spe::specie)::Float64
    rp = e2p(e, spe)
    return rp2beta(rp, spe)
end

"Convert β (= v/c) to the Lorentz factor γ"
function beta2gamma(beta)
    return 1 / sqrt(1 - beta^2)
end

function e2gamma(e::Float64, spe::specie)::Float64
    beta = e2beta(e, spe)
    return 1 / sqrt(1 - beta^2)
end

function e2v(e)
    return e2beta(e, spe) * CSPEED
end

function v2p(v, spe)
    return beta2p(v / CSPEED, spe)
end

function beta2p(beta::Float64, spe::specie)::Float64
    if spe.rnm > 0.5
        return beta*EP/sqrt(1.0-beta^2)*rnm/rnz
    else
        return beta*EE/sqrt(1.0-beta^2)
    end
end
end
