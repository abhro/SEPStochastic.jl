module fb
OMGP1::Float64 = 574745     #in 1/min
OMGE1::Float64 = 1055307413 #in 1/min

function fs0(
        tacc::Float64, xpk::Vector{Float64}, bm::Float64, u1::Float64,
        # densw in cm^-3, fs in cm^-3 (GeV/c)^-3,
        densw::Float64,
        ob::Float64, amach::Float64, rsh::Float64,
        # vth thermal speed of protons with downstream Maxwellian temperature
        vth::Float64,
        pinj::Float64, pc::Float64)

    length(xpk) == 6 || throw(DomainError("xpk must be a 6-vector"))

    # USES common block acc
    #   fampb is amplification of upstream magnetic field
    #   ratkp is parallel to perpendicular diff ceof ratio

    local fs0::Float64

    # USES common block specie

    sgm::Float64
    facob::Float64
    pth::Float64
    e1c::Float64
    f0::Float64
    omegae::Float64

    p1::Float64
    omegap::Float64
    p1inj::Float64

    ob2 = atan(rsh * sin(ob) / cos(ob))
    dob = amach / sqrt(40 * (1 + rsh^2))
    obinj = max(cos(ob), 0.25d0) #min(ob2,ob2-dob)
    etinj = 2.5 / obinj #injection at 3*u1
    #   shock slope
    #rpw = rsh * (amach - 1) / (amach + sqrt(rsh))
    #if rpw < 1.00001
    #   rpw = 1.00001
    #end
    sgm = 3*rsh/(rsh-1)
    #skp = sgm > 5.0 ? sgm/2 - 1 : 1.5

    #skp = 30.0

    #   obliuity factor for normal diffusion coeffiecient
    facob = sqrt(ratkp)
    facob = facob*cos(ob)^2 + sin(ob)^2/facob
    if rnm > 0.5
        pth = EP * vth / CSPEED #*sqrt((skp-1.5)/skp) #to kappa
        #   assume all species have the same thermal speed
        #    or temperature proportional to mass (Jacco Vink et al 2015,A&A)
        pth = rnm / rnz * pth #in rigidity
        p1 = EP * u1 / CSPEED
        omegap = OMGP1 * bm * fampb
        p1inj = etinj * p1
        pinj = rnm / rnz * p1inj
        e1c = hypot(EP, p1inj) + p1*p1/EP*(1-1/rsh)*tacc*omegap/facob*rnz/rnm
        pc = sqrt(e1c*e1c - EP*EP) * rnm / rnz
    else
        # same thermal speed as protons
        # no thermalization between e and p
        pth = EE * vth / CSPEED

        p1 = EE * u1 / CSPEED
        omegae = OMGE1*bm*fampb
        p1inj = etinj*p1
        pinj = p1inj
        e1c = hypot(EE, p1inj) + p1*p1/EE*(1-1/rsh)*tacc*omegae/facob/100
        pc = sqrt(e1c*e1c-EE*EE)
    end

    #ecr = (ecrmax + ecrmin - (ecrmax-ecrmin) * tanh((ob-ob0)/obw)) / 2

    #   only inject particles reached shock acceleration equilibrium
    if pc > xpk[4]
        #   matching distribution at pinj
        f0 = densw*rsh/(π*pth*pth)^1.5 / exp((pinj/pth)^2) #Maxwellian
        #f0 = densw * rsh / (PI*pth*pth)**1.5 * &
        #   exp(log_gamma(skp+1) - log_gamma(skp-0.5)) / skp**1.5 / &
        #   (1 + (pinj/pth)**2 / skp)**(skp+1)
        fs0 = f0 * sgm * u1 * (1 - 1/rsh) / 3 / (xpk[4]/pinj)^sgm
        if isnan(fs0)
            #write(*,*) fs0, densw, rsh, pth, pinj, pc, sgm
            fs0 = 0.0
        end
    else
        fs0 = 0.0
    end

    return fs0
end


function fb0(torg::Float64, rpb::Vector{Float64})
    length(rpb) == 5 || throw(DomainError("rpb must be a 5-vector"))
    # USES common block fbcnst
    # USES common block srcmod

    #rp = rpb[4]
    #if torg < 0 || torg > deltat
    #  fb0 = 0.0
    #else
    #  if nfbconst == 0
    #    tmodel = torg + tmodel0
    #    fbtimes = exp(-tc/tmodel - tmodel/tl) / tmodel
    #    ad = cos(rpb[2]) * trgtfs[2] + sin(rpb[2]) * trgtfs[1] *
    #        (cos(rpb[3]) * trgtfs[4] + sin(rpb[3]) * trgtfs[3])
    #    if ad < cos(scanw)
    #       fbtimes = 0.0
    #    end
    #  else
    #    fbtimes = 1.0
    #  end
    #    fb0 = rp2e(rp)^rk / (rp^2) * fbtimes
    #end

    fb0 = 0.0

    return fb0
end

#function getdeltat()
#end

#function getfbcnst(rk0, deltat0, tc0, tl0, tmodel00)
#end

end
