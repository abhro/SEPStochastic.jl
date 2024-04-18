module sim3d_utils
function f0mod(r::Vector{Float64}, pa::Float64)
    #  function modification factor for distribution
    #    f=f0mod*g -- g is the new function to be solve the stochastic

    length(r) == 3 || throw(DomainError("r must be a 3-vector"))

    df0 = zeros(Float64, 3)
    ddf0 = zeros(Float64, 3, 3)

    # rb0 is from fbcnst common block
    rs = zeros(Float64, 3)
    rs[1] = rb0
    rs[2] = r[2]
    rs[3] = r[3] + (r[1] - rs[1]) / vom # vom from bpark common block

    sintheta = sin(rs[2])
    if sintheta == 0
        sintheta = 1e-27
    end
    costheta = cos(rs[2])
    sinphi = sin(rs[3])
    cosphi = cos(rs[3])

    # trgtfs from srcmod common block
    cd = costheta * trgtfs[2] + sintheta * trgtfs[1] * (cosphi * trgtfs[4] + sinphi * trgtfs[3])
    dcddt = -sintheta * trgtfs[2] + costheta * trgtfs[1] * (cosphi * trgtfs[4] + sinphi * trgtfs[3])
    dcddp = -sintheta * trgtfs[1] * (sinphi * trgtfs[4] - cosphi * trgtfs[3])
    ddcddtt = -cd
    ddcddtp = -costheta * trgtfs[1] * (sinphi * trgtfs[4] - cosphi * trgtfs[3])
    ddcddpp = -sintheta * trgtfs[1] * (cosphi * trgtfs[4] + sinphi * trgtfs[3])
    ad = acos(cd)
    sd = sin(ad)
    daddt = -1 / sd * dcddt
    daddp = -1 / sd * dcddp
    ddaddtt = -cd / sd * daddt * daddt - 1 / sd * ddcddtt
    ddaddtp = -cd / sd * daddp * daddt - 1 / sd * ddcddtp
    ddaddpp = -cd / sd * daddp * daddp - 1 / sd * ddcddpp

    #f0 = (1+pa/sp) * exp(-ad*ad/ap/ap/2)
    f0 = log(1+pa/sp) - ad*ad/ap/ap/2   #used as exp(f0)

    df0dmu = 1/(sp+pa)   #normalized by f0
    ddf0dmu2 = 0 #gp*(gp-1)/(sp+pa)**2

    #  df0(1-3) = derivative in sphere coordinates /f0
    #    first order
    if sd == 0
        df0[2] = 0.0
        df0[3] = 0.0
    else
        df0[2] = -ad / ap / ap * daddt
        df0[3] = -ad / ap / ap * daddp
    end
    df0[1] = df0(3)/vom
    #    second order
    if sd == 0
      if ad == 0
        ddf0[2,2] = df0[2] * df0[2] - 1 / ap / ap
        ddf0[2,3] = df0[2] * df0[3] + 0
        ddf0[3,3] = df0[3] * df0[3]
      else
        sdp = 1e-10
        ddf0[2,2] = df0[2] * df0[2] + PI / ap / ap / sdp
        ddf0[2,3] = df0[2] * df0[3] + 0
        ddf0[3,3] = df0[3] * df0[3] + PI / ap / ap / sdp * trgtfs[1]^2
      end
    else
      #ddsd = (sd-ad*cd)/sd/sd/ap/ap
      #ddf0[2,2] = df0[2]*df0[2] - ad*cd/ap/ap/sd + ddsd*dcddt*dcddt
      #ddf0[2,3] = df0[2]*df0[3] + ad/ap/ap/sd*ddcddtp + ddsd*dcddt*dcddp
      #ddf0[3,3] = df0[3]*df0[3] + ad/ap/ap/sd*ddcddpp + ddsd*dcddp*dcddp
      ddf0[2,2] = -ad/ap/ap*ddaddtt + df0[2]*df0[2]*(1-ap*ap/ad/ad)
      ddf0[2,3] = -ad/ap/ap*ddaddtp + df0[2]*df0[3]*(1-ap*ap/ad/ad)
      ddf0[3,3] = -ad/ap/ap*ddaddpp + df0[3]*df0[3]*(1-ap*ap/ad/ad)
    end
    #ddf0[1,1] = ddf0[3,3] / vom^2
    ddf0[1,1] = -(ad*ddaddpp+daddp*daddp)/ap/ap/vom/vom + df0[1]^2
    #ddf0(1,2) = ddf0(2,3)/vom
    ddf0[1,2] = -(ad*ddaddtp+daddt*daddp)/ap/ap/vom + df0[1]*df0[2]
    #ddf0(1,3) = ddf0(3,3)/vom
    ddf0[1,3] = -(ad*ddaddpp+daddp*daddp)/ap/ap/vom + df0[1]*df0[3]
    #   gradient-f0 /f0 in spherical
    df0[2] = df0[2] / r[1]
    df0[3] = df0[3] / (r[1]*sintheta)
    #   grdaient-gradient-f0 /f0
    ddf0[2,1] = ddf0[1,2]/r[1] - df0[2]/r[1]^2
    ddf0[3,1] = (ddf0[1,3] - df0(3)/r(1))/r(1)/sintheta
    ddf0[2,2] = (ddf0[2,2]/r[1] + df0[1])/r[1]
    ddf0[3,2] = (ddf0[2,3]-costheta/sintheta*df0[3])/r[1]/r[1]/sintheta
    ddf0[3,3] = (ddf0[3,3]/r[1]/sintheta + df0[1]*sintheta + df0[2]/r[1]*costheta)/r[1]/sintheta
    ddf0[1,2] = ddf0[2,1]
    ddf0[1,3] = ddf0[3,1]
    ddf0[2,3] = ddf0[3,2]

    return (f0=f0, df0=df0, ddf0=ddf0, df0dmu=df0dmu, ddf0dmu2=ddf0dmu2)
end

function solarwindtemp(r)

    #  empirical model in the corona from Withbroe (ApJ 325,442,1988) 10.1086/166015

    length(r) == 3 || throw(DomainError("r must be a 3-vector"))

    if r[1] > 10.0
        temp = 1.4e6 * (10/r[1])^1.3333
    end

    if r[1] <= 10.0 && r[1] >= 1.125
        temp = 1.4e6
    end

    if r[1] < 1.125
      temp = (1.0e5^3.5 + 25.33965e6 * (1/1.0287 - 1/r[1])) ^ (1/3.5)
    end

    return temp
end

function compress(amach::Float64, smach::Float64, ob::Float64)
    rsh = zeros(Float64, 3)
    root = zeros(ComplexF64, 3)

    pbeta = (amach/smach)^2
    sintheta2 = sin(ob)^2
    amach2 = amach * amach
    amach4 = amach2 * amach2
    amach6 = amach4 * amach2

    a3 = (-2*pbeta + pbeta*sintheta2*4 - pbeta*sintheta2*sintheta2*2
          - amach2*sintheta2
          - amach2*gamma + amach2 + amach2*gamma*sintheta2)
    a2 = -amach2 * (
      -gamma - pbeta*4 + gamma*sintheta2 + pbeta*sintheta2*4 + sintheta2
      -amach2*gamma*2 + amach2*2 + amach2*gamma*sintheta2 - 1)
    a1 = -amach4 * (
      gamma*2 + pbeta*2 - gamma*sintheta2 - sintheta2*2 + amach2*gamma - amach2 + 2)
    a0 = amach6 * gamma + amach6

    rsh0 = -a2 / (3*a3)
    b1 = (3*a3*a1 - a2*a2) / (3*a3*a3)
    b0 = (2*a2*a2*a2-9*a3*a2*a1+27*a3*a3*a0) / (27*a3*a3*a3)
    del = 4*b1*b1*b1 + 27*b0*b0


    if del < 0
        pr3 = sqrt(Complex(-b1/3))
        wth = 3*b0/(2*b1) * sqrt(Complex(-3/b1))
        wth = acos(wth) / 3
        for k = 0:2
            root[k+1] = 2*pr3*cos(wth-2*PI*k/3) + rsh0
            if abs(aimag(root(k+1))) < 1.0d-10
                rsh[k+1] = real(root[k+1])
            else
                rsh[k+1] = 1
            end
        end
    else
        if b1 < 0
            wth = -3*abs(b0) / (2*b1)*sqrt(Complex(-3/b1))
            wth = acosh(wth)/3
            #wth = log(wth + sqrt(wth + 1) * sqrt(wth - 1)) / 3
            root[1] = -2*abs(b0)/b0*(-b1/3)^0.5*cosh(wth) + rsh0
        else
            wth = 3 * b0/(2*b1)*(3/b1)^0.5
            #wth = asinh(wth)/3
            wth = log(wth + sqrt(wth * wth + 1)) / 3
            root[1] = -2*(b1/3)^0.5*sinh(wth) + rsh0
        end
        root[2] = 1
        root[3] = 1
        rsh[:] = real(root(1:3))
    end
    return rsh
end

function split(iseed, rpb, ck, fs, t, nsplvl, dnsk, bv0, flx, dflx, walk3d)
end
end
