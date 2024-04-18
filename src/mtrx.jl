module mtrx

using LinearAlgebra: cross, dot, norm, normalize, tr
using StaticArrays: @MVector, MVector

push!(LOAD_PATH, @__FILE__)
include("param.jl")
using .param: N_R, N_θ, N_φ


"""calculate matrix for xyz' ellipsoid coordinate Kwon to spheric"""
function mxptr(gm::Float64)::Matrix{Float64}
    sing = sin(gm)
    cosg = cos(gm)

    return [0.0 0.0 1.0; cosg -sing 0.0; sing cosg 0.0]
end

"""
calculate matrix for xyz' ellipsoid coordinate Kwon to spheric
"""
function dmxptr(gm::Float64, dgm::Float64)::Matrix{Float64}
    sing = sin(gm)
    cosg = cos(gm)
    return [
        0.0       0.0       0.0
        -sing*dgm -cosg*dgm 0.0
        cosg*dgm  -sing*dgm 0.0
    ]
end


"""calculate martix from magnetic to polar spheric coordinates"""
function mbtr(
        uax1::Vector{Float64},
        uax2::Vector{Float64},
        uax3::Vector{Float64}
    )::Matrix{Float64}

    return [uax1 uax2 uax3]
end

"""calculate martix from polar spheric to xyz coordinates"""
function mrtx(sinθ::Float64, cosθ::Float64,
              sinφ::Float64, cosφ::Float64)::Matrix{Float64}

    return [sinθ*cosφ  cosθ*cosφ  -sinφ;
            sinθ*sinφ  cosθ*sinφ   cosφ;
            cosθ         -sinθ      0.0]
end

"""calculate martix from polar spheric to xyz coordinates"""
function dmrtx(
        sinθ::Float64, cosθ::Float64,
        sinφ::Float64, cosφ::Float64,
        dθ::Float64, dφ::Float64
    )::Matrix{Float64}

    dr2x = [
            cosθ*dθ*cosφ - sinθ*sinφ*dφ -sinθ*dθ*cosφ - cosθ*sinφ*dφ -cosφ*dφ;
            cosθ*dθ*sinφ + sinθ*cosφ*dφ -sinθ*dθ*sinφ + cosθ*cosφ*dφ -sinφ*dφ;
            -sinθ*dθ                    -cosθ*dθ                      0.0
           ]
    return dr2x
end

"""calculate solar wind velocity in the corotating frame
and its gradient in sphereical coordinate system
"""
function solarwind1(r::Vector{Float64})::Vector{Float64}
    # TODO get the common block variables
    if r[1] < 2.5
        vpl = zeros(Float64, 3)
    else
        vpl = [
               vsw / (1 + k4ok2/r[1]^2 + k6ok2/r[1]^4),
               0.0,
               -omega * (r[1]-2.5) * sin(r[2]) #rss=2.5Rs for pfss
              ]
    end
    vpl
end


"""Trilinear interpolation"""
function trilinear(
        #  fc the value of f at the corner of cubic box of side 1
        fc::AbstractArray{Float64, 3},
        #  location inside the cube (0<=x<=1) or outside x<0 x>1
        x::Vector{Float64}
    )::Float64

    #  interpolate f
    f = fc[1,1,1] * (1-x[1]) * (1-x[2]) * (1-x[3]) +
        fc[2,1,1] *    x[1]  * (1-x[2]) * (1-x[3]) +
        fc[1,2,1] * (1-x[1]) *    x[2]  * (1-x[3]) +
        fc[1,1,2] * (1-x[1]) * (1-x[2]) *    x[3]  +
        fc[2,1,2] *    x[1]  * (1-x[2]) *    x[3]  +
        fc[1,2,2] * (1-x[1]) *    x[2]  *    x[3]  +
        fc[2,2,1] *    x[1]  *    x[2]  * (1-x[3]) +
        fc[2,2,2] *    x[1]  *    x[2]  *    x[3]

    f
end

function trilineardif(
        #  fc the value of f at the corner of cubic box of side 1
        fc::AbstractArray{Float64, 3},
        #  location inside the cube (0<=x<=1) or outside x<0 x>1
        x::Vector{Float64}
    )::MVector{3, Float64}

    df = @MVector [
                   (fc[2,1,1] - fc[1,1,1]) * (1-x[2]) * (1-x[3]) +
                   (fc[2,2,1] - fc[1,2,1]) *    x[2]  * (1-x[3]) +
                   (fc[2,1,2] - fc[1,1,2]) * (1-x[2]) *    x[3]  +
                   (fc[2,2,2] - fc[1,2,2]) *    x[2]  *    x[3],

                   (fc[1,2,1] - fc[1,1,1]) * (1-x[1]) * (1-x[3]) +
                   (fc[2,2,1] - fc[2,1,1]) *    x[1]  * (1-x[3]) +
                   (fc[1,2,2] - fc[1,1,2]) * (1-x[1]) *    x[3]  +
                   (fc[2,2,2] - fc[2,1,2]) *    x[1]  *    x[3],

                   (fc[1,1,2] - fc[1,1,1]) * (1-x[1]) * (1-x[2]) +
                   (fc[2,1,2] - fc[2,1,1]) *    x[1]  * (1-x[2]) +
                   (fc[1,2,2] - fc[1,2,1]) * (1-x[1]) *    x[2]  +
                   (fc[2,2,2] - fc[2,2,1]) *    x[1]  *    x[2]
                  ]
    df
end

function vfunc(t, xpk, dxpkdt, du, gxw2, gxw3, bv0, densw, vpl, gper, b1s)
    length(xpk) == 6 || throw(DomainError("xpk must be a 6 element vector"))

    r = zeros(3)
    r[1] = norm(xpk)
    r[2] = acos(xpk[3]/r[1])
    r[3] = atan(xpk[2], xpk[1])

    p = xpk[4]
    beta = rp2beta(p)
    pa = xpk[5]

    sinθ = sin(r[2])
    if sinθ == 0
        sinθ = 1e-27
    end
    cosθ = cos(r[2])
    sinφ = sin(r[3])
    cosφ = cos(r[3])

    #drvbmag(r, bv, bm, cvtu, gbmag, bxgb2, dbbds, pol, b1s, gb1s)
    (b, bmag, cvtu, gbmag, bxgb2, dbbds, pol, b1rs, gb1rs) = drvbmg(r)

    bv0 = bv
    lpb = pol

    # local coordinates
    #   unit outward field line direction: uax1
    uax1 = lpb * bv / vm
    #   unit perpendicular directions
    uax2::Vector{Float64} = [0, -lpb*bv[3], lpb*bv[2]]
    um = hypot(uax2[2], uax2[3])
    if um == 0
        uax2[2] = 1
        uax2[3] = 0
    else
        uax2[2] = uax2[2]/um
        uax2[3] = uax2[3]/um
    end
    uax3 = normalize(cross(uax1, uax2))

    # focusing along outward field line
    dbbds = lpb * dbbds

    # calculate matrix from polar spheric to xyz coordinates
    r2x = mrtx(sinθ, cosθ, sinφ, cosφ)
    # calculate matrix from magnetic to polar spheric coordinates
    b2r = mbtr(uax1, uaax2, uax3)
    # calculate matrix from magnetic to xyz coordinates
    b2x = r2x * b2r # matrix multiplication

    # diffusion tensor in magnetic coordinates
    cofm(r, p, pa, beta, bv, bm, cvtu, gbmag, dbbds, b1s, gb1s, gb, dgr)

    #  calculate diffusion coeficients
    #     in polar spheric coordinates
    gper = gb[2] # gb[2] = gb[3]
    gr = zeros(3, 3)
    for i = 1:3
        for j = 1:3
            gr[i,j] = sum(@. b2r[i,:] * gb[:] * b2r[j,:])
        end
    end

    #   1. parallel particle speed + solar wind speed
    vpcl = beta * CSPEED
    #   2. drift velocity in spherical coordinates
    dcs = 1 / norm(bxgb2)
    sinmu = 1 - pa*pa
    sinmu = sinmu < 0 ? 0 : sqrt(sinmu)
    rgmu2 = rg1 * p / bm * sinmu * 1.414 # sqrt(2),linear; 2,step
    if dcs > rgmu2 # regular drift
      culpar = dot(bv(1:3), cvtu(1:3))/bm/bm*bv
      culper = cvtu - culpar
      vd = rg1 * p * vpcl / bm * (pa*pa*culper/bm + (1-pa*pa)/2/bm*culpar +
         (1+pa*pa)/2*bxgb2)
    else # current sheet drift (square delta function
      vd = -vpcl * sinmu / 2.828 * bxgb2 * dcs
    end

    if rnz < 0
        vd = -vd
    end

    #   3. solar wind speed
    #   gvpl::Matrix{Float64} # 3x3
    (vpl, gvpl, densw) = solarwind(r)

    #   4 artificial drift
    # df0::Vector{Flaot64} # 3-vector
    f0mod(r, pa, f0, df0, ddf0, df0dmu, ddf0dmu2)
    avr = 2 * gr * df0
    avx = r2x * avr

    #   total
    for i = 1:3
      vx[i] = dot(r2x[i,1:3], vpl[1:3])
      dxpkdt[i] = -b2x[i,1]*vpcl*pa - vx[i] - dot(r2x[i,:], vd) - dot(r2x[i,:], vd)
      #  add dg and artifical drift
      dxpkdt[i] = dxpkdt[i] + avx[i] + dot(r2x[i,:], dgr[1:3])
    end

    #  get pitch angle diffusion coefficient
    cofdu(p, pa, beta, bm, du, ddu)
    if r[1] > 20.0
        cosp = bv[1] / bm
        du = du * cosp * cosp / facip
        ddu = ddu * cosp * cosp / facip
    end

    #   pitch angle drift term
    #   focusing
    dpadt = -vpcl * dbbds * (1 - pa^2) / 2.0

    #   cooling
    bbgv = 0.0
    for i = 1:3
        for j = 1:3
            bbgv = bv[i]/bm*gvpl[i,j]*bv[j]/bm + bbgv
        end
    end
    divv = tr(gvpl)
    #divv1 = 2*vpl[1]/r[1] + gvpl[1,1]  # for test
    dpadt = dpadt + pa * (1-pa^2) / 2 * (divv - 3*bbgv)
    #   derivative of pitch-angle diffusion
    #   artificial due to modification of function
    adpadt = 2.0 * du * df0dmu
    dxpkdt[5] = -dpadt + ddu + adpadt

    #   momentum loss
    dxpkdt[4] = ((1-pa^2) / 2 * (divv-bbgv) + pa*pa*bbgv) * p

    #   killing factor
    fk = (-dpadt+ddu)*df0dmu + du*ddf0dmu2
    fk += (
      gr[1,1]*ddf0[1,1] + gr[1,2]*ddf0[2,1] + gr[1,3]*ddf0[3,1] +
      gr[2,1]*ddf0[1,2] + gr[2,2]*ddf0[2,2] + gr[2,3]*ddf0[3,2] +
      gr[3,1]*ddf0[1,3] + gr[3,2]*ddf0[2,3] + gr[3,3]*ddf0[3,3] )
    fk += dot(dgr, df0) - dot(vd,  df0) - dot(vpl, df0)
    dxpkdt[6] = fk

    #  diffusion vector applied to dw2 and dw3
    gxw2[:] = b2x[1:3,2] * sqrt(2 * gb[2])
    gxw3[:] = b2x[1:3,3] * sqrt(2 * gb[3])

    #write(54,"(f14.7,13(1pe14.5))") t, xpk, dxpkdt,bm
    #call flush(54)

    # TODO: output
end

function solarwind(
        r::Vector{Float64},
    )::Tuple{Vector{Float64}, Matrix{Float64}, Float64}
    #   calculate solar wind velocity in the corotating frame
    #     and its gradient
    #     in spherical coordinate system


    #     use leBalnc(1998) model
    # TODO port common block


    len(r) == 3 || throw(DomainError("r must be a 3-vector"))

    sinθ = sin(r[2])
    sinθp = sinθ == 0 ? 1e-20 : sinθ
    cosθ = cos(r[2])
    rpl = 1 + k40k2 / r[1]^2 + k6ok2 / r[1]^4
    if r[1] < 2.5
        vpl = zeros(Float64, 3)
        gvpl = zeros(Float64, 3, 3)
    else
        vpl::Vector{Float64} = [
                                vsw/rpl,
                                0,
                                -omega*(r[1]-2.5)*sinθ #rss=2.5Rs for pfss
                               ]
        r1 = r[1]
        gvpl11 = vpl[1] * (2*k4ok2/r1^3 + 4*k6ok2/r1^5) /
                          (1 + k4ok2/r1^2 + k6ok2/r1*4)
        gvpl::Matrix{Float64} = [
                                 gvpl11         0.0             -omega*sinθ;
                                 0.0        vpl[1]/r[1]        -omega*(r1-2.5)*cosθ/r1;
                                 -vpl[3]/r1  -cosθ/sinθp*vpl[3]/r1  vpl[1]/r1
                                ]
    end
    densw = densw0 * rpl / r[1]^2
    return (vpl, gvpl, densw)
end

function solarwind1(r::Vector{Float64})::Vector{Float64}
    #   calculate solar wind velocity in the corotating frame
    #     and its gradient in sphereical coordinate system

    len(r) == 3 || throw(DomainError("r must be a 3-vector"))
    #     use leBalnc(1998) model
    # TODO get common variables
    if r[1] < 2.5
        vpl = zeros(Float64, 3)
    else
        vpl = [
               vsw / (1 + k4ok2/r[1]^2 + k6ok2/r[1]^4),
               0.0,
               -omega * (r[1]-2.5) * sin(r[2]) # rss=2.5Rs for pfss
              ]
    end

    return vpl
end

function drvbmag(r1::Vector{Float64})
    # the following are output variables #  XXX pack into a struct?
    # b::Vector{Float64} # 3-vector
    # bmag::Float64
    # cvtu::Vector{Float64} # 3-vector
    # gbmag::Vector{Float64} # 3-vector
    # bxgb2::Vector{Float64} # 3-vector
    # dbbds::Float64
    # pol::Float64
    # b1rs::Float64
    # gb1rs::Vector{Float64} # 3-vector
    len(r1) == 3 || throw(DomainError("r1 must be a 3-vector"))

    # TODO get common block variables

    RSS = 2.5

    r = r1
    if r[1] > rss
        tpsw = ((r[1] - k4ok2/r[1] - k6ok2/r[1]^3/3) &
        - (rss - k4ok2/rss - k6ok2/rss^3/3) ) / vsw
    else
        tpsw = 0.0
    end
    r[2] = acos(cos(r[2]))
    r[3] = r[3] + omega * tpsw
    r[3] = atan(sin(r[3]), cos(r[3]))
    if r[3] < 0.0
        r[3] = r[3] + 2π
    end
    if r[1] > RSS
        r[1] = RSS
    end
    if r[1] < 1.0
        r[1] = 1.0
    end
    sinθ = sin(r[2])
    if sinθ == 0.0
        sinθ = 1e-20
    end

    #  find the grid cell
    irr = floor((r[1]-1) / (RSS-1) * N_R)
    if irr >= N_R
        irr = N_R - 1
    end

    itheta = floor(r[2] / π * N_θ)
    if itheta >= N_θ
        itheta = N_θ - 1
    end

    iphi = floor(r[3] / 2π * N_φ)
    if iphi >= N_φ
        iphi = N_φ - 1
    end

    #  relative displacement from lower grids
    px = [(r[1]-1.0) / (RSS-1.0) * N_R - irr,
          r[2] / π * N_θ - itheta,
          r[3] / 2π* N_φ - iphi]

    fc = zeros(Float64, 2, 2, 2)
    for m = 1:3
      fc[1,1,1] = bgrid[irr,   itheta,   iphi,   m]
      fc[2,1,1] = bgrid[irr+1, itheta,   iphi,   m]
      fc[1,2,1] = bgrid[irr,   itheta+1, iphi,   m]
      fc[2,2,1] = bgrid[irr+1, itheta+1, iphi,   m]
      fc[1,1,2] = bgrid[irr,   itheta,   iphi+1, m]
      fc[2,1,2] = bgrid[irr+1, itheta,   iphi+1, m]
      fc[1,2,2] = bgrid[irr,   itheta+1, iphi+1, m]
      fc[2,2,2] = bgrid[irr+1, itheta+1, iphi+1, m]
      f = trilinear(fc, px)
      b(m) = f
      # fc[1,1,1] = cvgrid[irr,     ith,     iph,     m]
      # fc[2,1,1] = cvgrid[irr + 1, ith,     iph,     m]
      # fc[1,2,1] = cvgrid[irr,     ith + 1, iph,     m]
      # fc[2,2,1] = cvgrid[irr + 1, ith + 1, iph,     m]
      # fc[1,1,2] = cvgrid[irr,     ith,     iph + 1, m]
      # fc[2,1,2] = cvgrid[irr + 1, ith,     iph + 1, m]
      # fc[1,2,2] = cvgrid[irr,     ith + 1, iph + 1, m]
      # fc[2,2,2] = cvgrid[irr + 1, ith + 1, iph + 1, m]
      # f = trilinear(fc,px)
      # cvtu[m] = f
      fc[1,1,1] = gbgrid[irr,   itheta,   iphi,   m]
      fc[2,1,1] = gbgrid[irr+1, itheta,   iphi,   m]
      fc[1,2,1] = gbgrid[irr,   itheta+1, iphi,   m]
      fc[2,2,1] = gbgrid[irr+1, itheta+1, iphi,   m]
      fc[1,1,2] = gbgrid[irr,   itheta,   iphi+1, m]
      fc[2,1,2] = gbgrid[irr+1, itheta,   iphi+1, m]
      fc[1,2,2] = gbgrid[irr,   itheta+1, iphi+1, m]
      fc[2,2,2] = gbgrid[irr+1, itheta+1, iphi+1, m]
      f = trilinear(fc,px)
      gbmag[m] = f
    end


    for i1 = 0:1
        for i2 = 0:1
            for i3 = 0:1
                bx  =  bgrid[irr+i1, itheta+i2, iphi+i3, 1]
                by  =  bgrid[irr+i1, itheta+i2, iphi+i3, 2]
                bz  =  bgrid[irr+i1, itheta+i2, iphi+i3, 3]

                bm  = hypot(bx, by, bz)

                #gbx = gbgrid[irr+i1, itheta+i2, iphi+i3, 1]
                gby = gbgrid[irr+i1, itheta+i2, iphi+i3, 2]
                gbz = gbgrid[irr+i1, itheta+i2, iphi+i3, 3]
                fc[i1+1,i2+1,i3+1] = (by*gbz-bz*gby) / bm^2
            end
        end
    end
    bxgb2[1] = trilinear(fc, px)

    for i1 = 0:1
        for i2 = 0:1
            for i3 = 0:1
                bx = bgrid[irr+i1, itheta+i2, iphi+i3, 1]
                by = bgrid[irr+i1, itheta+i2, iphi+i3, 2]
                bz = bgrid[irr+i1, itheta+i2, iphi+i3, 3]

                bm = hypot(bx, by, bz)

                gbx = gbgrid[irr+i1, itheta+i2, iphi+i3, 1]
                #gby = gbgrid[irr+i1, itheta+i2, iphi+i3, 2]
                gbz = gbgrid[irr+i1, itheta+i2, iphi+i3, 3]
                fc[i1+1,i2+1,i3+1] = (bz*gbx-bx*gbz)/bm^2
            end
        end
    end
    bxgb2[2] = trilinear(fc, px)

    for i1 = 0:1
        for i2 = 0:1
            for i3 = 0:1
                bx = bgrid[irr+i1, itheta+i2, iphi+i3, 1]
                by = bgrid[irr+i1, itheta+i2, iphi+i3, 2]
                bz = bgrid[irr+i1, itheta+i2, iphi+i3, 3]

                bm = hypot(bx, by, bz)

                gbx = gbgrid[irr+i1, itheta+i2, iphi+i3, 1]
                gby = gbgrid[irr+i1, itheta+i2, iphi+i3, 2]
                #gbz = gbgrid[irr+i1, itheta+i2, iphi+i3, 3]
                fc[i1+1,i2+1,i3+1] = (bx*gby-by*gbx)/bm^2
            end
        end
    end
    bxgb2[3] = trilinear(fc,px)


    fc[1,1,1] = abs(b1rsgrid[irr,   itheta,   iphi  ])
    fc[2,1,1] = abs(b1rsgrid[irr+1, itheta,   iphi  ])
    fc[1,2,1] = abs(b1rsgrid[irr,   itheta+1, iphi  ])
    fc[2,2,1] = abs(b1rsgrid[irr+1, itheta+1, iphi  ])
    fc[1,1,2] = abs(b1rsgrid[irr,   itheta,   iphi+1])
    fc[2,1,2] = abs(b1rsgrid[irr+1, itheta,   iphi+1])
    fc[1,2,2] = abs(b1rsgrid[irr,   itheta+1, iphi+1])
    fc[2,2,2] = abs(b1rsgrid[irr+1, itheta+1, iphi+1])
    f = trilinear(fc, px)
    df = trilineardif(fc, px)
    b1rs = f
    gb1rs = [df[1] / ((rss-1)/N_R),
             df[2] / (π/N_θ)/r[1],
             df[3] / (2π/N_φ)/r[1]/sinθ]

    (nplus, nminus) = (0, 0)
    (bplus, bminus) = (0.0, 0.0)
    for m3 = 0:1
        for m2 = 0:1
            for m1 = 0:1
                if b1rsgrid[irr+m1, itheta+m2, iphi+m3] > 0
                    nplus += 1
                    bplus[1] += bgrid[irr+m1,itheta+m2,iphi+m3,1]
                    bplus[2] += bgrid[irr+m1,itheta+m2,iphi+m3,2]
                    bplus[3] += bgrid[irr+m1,itheta+m2,iphi+m3,3]
                else
                    nminus += 1
                    bminus[1] += bgrid[irr+m1,itheta+m2,iphi+m3,1]
                    bminus[2] += bgrid[irr+m1,itheta+m2,iphi+m3,2]
                    bminus[3] += bgrid[irr+m1,itheta+m2,iphi+m3,3]
                end
            end
        end
    end
    if nplus > 0 && nminus > 0
        bplus = bplus / nplus
        bminus = bminus / nminus
        if dot(bplus, b) > dot(bminus, b)
            pol = 1.0
        else
            pol = -1.0
        end
    else
      if nplus == 0
          pol = -1.0
      end
      if nminus == 0
          pol = 1.0
      end
    end

    bmag = norm(b)
    dbbds = dot(b, gbmag) / bmag^2
    cvtu = zeros(Float64, 3) # potential field
    if r1[1] >= rss
      rsovr = rss / r1[1]
      plbs = b[1] >= 0 ? 1 : -1
      b[1] = b[1] * rsovr^2
      vpl = solarwind1(r1)
      tanp = vpl[3] / vpl[1]
      oneptan2 = 1 + tanp^2
      sroneptan2 = sqrt(oneptan2)
      b[3] = b[1] * tanp
      bmag = norm(b)
      dvr1 = vsw * (2*k4ok2/r1[1]^3 + 4*k6ok2/r1[1]^5) / &
        (1 + k4ok2/r1[1]^2 + k6ok2/r1[1]^4)^2
      dvr1 = dvr1/vpl[1]
      dvf1 = -omega * sinθ / vpl[3]
      dvf2 = -omega * (r1[1]-rss) * cos(r1[2]) / r1[1] / vpl[3]
      gbmags = gbmag
      gbmag[1] = -2*bmag/r1[1] + bmag*tanp^2/oneptan2*(dvf1-dvr1)
      gbmag[2] = gbmags[2] * rsovr^3 * sroneptan2 + bmag * tanp^2 / oneptan2 * dvf2
      gbmag[3] = gbmags[3] * rsovr^3 * sroneptan2
      bxgb2s = bxgb2
      bxgb2[1] = - bxgb2s[3]*rsovr*tanp/sroneptan2 - b[1]*tanp^3/oneptan2*dvf2/bmag
      bxgb2[2] = bxgb2[2]*rsovr/sroneptan2 - 2*b(1)/r1(1)*tanp/bmag + b(1)*tanp^3/oneptan2*(dvf1-dvr1)/bmag
      bxgb2[3] = bxgb2[3]*rsovr/sroneptan2 + b[1]*tanp^2/oneptan2*dvf2/bmag
      cvtu[1] = b[3] * (dvf2 + cos(r1[2])/sinθ/r1[1]) + tanp * rsovr^3 * plbs * gbmags[2]
      cvtu[2] = plbs * gbmags[3] * rsovr^3 + b[3] * (1/r1[1] - dvf1 + dvr1)
      cvtu[3] = -plbs * gbmags[2] * rsovr^3
      dbbds = dot(b, gbmag) / bmag^2
      b1rs = b1rs / rsovr^2 / sroneptan2
      gb1rs[1] = 2 * b1rs / r1[1] - b1rs * tanp^2 / oneptan2 * (dvf1-dvr1)
      gb1rs[2] = gb1rs[2] / rsovr / sroneptan2 - b1rs * tanp^2 / oneptan2 * dvf2
      gb1rs[3] = gb1rs[3] / rsovr / sroneptan2
    end

    return (
        b=b,
        bmag=bmag,
        cvtu=cvtu,
        gbmag=gbmag,
        bxgb2=bxgb2,
        dbbds=dbbds,
        pol=pol,
        b1rs=b1rs,
        gb1rs=gb1rs
    )
end # drvbmag
end
