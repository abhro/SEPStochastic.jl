module cme_cross
using OffsetArrays: zeros

using .mtrx: dmxptr

function preparecme
    # XXX using cmesk common block
    # XXX using dir common block
    # XXX using bpark common block
    # XXX using tmprm common block
    # XXX using vsksw common block
end

# TODO copy newer version from fortran
function inorout(tsh::Float64, x::Vector{Float64}, dnsk::Float64)

    dnsk0 = dnsk # make copy

    if tsh < tska(nsk)
      jt = locate(tska,nsk,tsh)
      if jt == 0
         jt = 1
      end
      #  maximum shock radial distance for quick search
      drskc = (pska(jt+1,1)-pska(jt,1))/(tska(jt+1)-tska(jt))
      rskc = pska(jt,1)+(tsh-tska(jt))*drskc
      dask[3] = (pska(jt+1,6)-pska(jt,6))/(tska(jt+1)-tska(jt))
      ask[3] = pska(jt,6)+(tsh-tska(jt))*dask(3)
      #vskf = vskf0 !only for fort.58; useless
      #rskf = pska(jt,1)+pska(jt,6)!only for fort.58; useless
      dthskc = (pska(jt+1,2)-pska(jt,2))/(tska(jt+1)-tska(jt))
      θ_skc = pska(jt,2)+(tsh-tska(jt))*dthskc
      dφskc = (pska(jt+1,3)-pska(jt,3))
      if (dφskc >  PI)
        dφskc = dφskc - 2π
      end
      if (dφskc < -PI)
        dφskc = dφskc + 2π
      end
      dφskc = dφskc / (tska(jt+1) - tska(jt))
      φskc = pska(jt,3)+(tsh-tska(jt))*dφskc
    else
      jt = nsk-1
      vcme0 = vskf0
      if tsh < tauc2_0
        rskf = pska(nsk,1)+pska(nsk,6)+vcme0*(tsh-tska[nsk]) #LC
        vskf = vcme0
      else
        if tauc2_0>tska(nsk)
          rskf = 3.0/2.0 * (vcme0-vsw)*tauc2*(((tsh-tska(1))/tauc2)^(2.0/3.0) -1.0e0) +
            vsw*tauc2*(((tsh-tska(1))/tauc2)-1.0d0) +
            pska(nsk,1)+pska(nsk,6)+vcme0*(tauc2_0-tska(nsk)) #LC
        else
          rskf = 3.0/2.0 * (vcme0-vsw)*tauc2*(((tsh-tska(1))/tauc2)^(2.0/3.0)
            -((tska(nsk)-tska(1))/tauc2)^(2.0/3.0)) + vsw*(tsh-tska(nsk))
            + pska(nsk,1)+pska(nsk,6)
        end
        vskf = (vcme0-vsw)*(((tsh-tska(1))/tauc2)^(-1.0d0/3.0d0)) + vsw
      end

      rhclf = pska(nsk,6) / (pska(nsk,1) + pska(nsk,6))
      rskc = (1-rhclf) * rskf
      drskc = (1-rhclf) * vskf
      if (pska(nsk,1) - pska(nsk-1,1)) < 1.e-3
        drskc = 0.0
        rskc = pska(nsk,1)
      end
      ask(3) = rhclf * rskf
      dask(3) = rhclf * vskf
    end
    rskmax = rskc + ask(3)
    rskmin = rskc - ask(3)
    ra = norm2(x(1:3))
    #write(58,"(f14.8,3(1pe12.4))") tsh,ra,rskf,vskf
    if (rskmax < ra*.99) || (rskmin > ra*1.01)
      dnsk = 0.01*ra  # use step size as a minimum
      if dnsk0 * dnsk >= 0.0
        vsk = 1e-6  #arbitrarily assigned
        vnx = 0.57735026918962576450e0  #arbitrarily assigned
        return
      end
    end

    if tsh < tska(nsk)
      dthskc = (pska(jt+1,2)-pska(jt,2)) / (tska(jt+1)-tska(jt))
      θ_skc = pska(jt,2)+(tsh-tska(jt))*dthskc
      dφskc = pska(jt+1,3) - pska(jt,3)
      if (dφskc >  PI)
        dφskc = dφskc - 2π
      end
      if (dφskc < -PI)
        dφskc = dφskc + 2π
      end
      dφskc = dφskc/(tska(jt+1)-tska(jt))
      φskc = pska(jt,3)+(tsh-tska(jt))*dφskc
      dgmsk = (pska(jt+1,7)-pska(jt,7))/(tska(jt+1)-tska(jt))
      gmsk = pska(jt,7)+(tsh-tska(jt))*dgmsk

      dask(1) = (pska(jt+1,4)-pska(jt,4))/(tska(jt+1)-tska(jt))
      ask(1) = pska(jt,4)+(tsh-tska(jt))*dask(1)
      dask(2) = (pska(jt+1,5)-pska(jt,5))/(tska(jt+1)-tska(jt))
      ask(2) = pska(jt,5)+(tsh-tska(jt))*dask(2)

      dtzskmin = (pska(jt+1,8)-pska(jt,8)) / (tska(jt+1)-tska(jt))
      tzskmin = pska(jt,8)+(tsh-tska(jt))*dtzskmin
    else
      tsh1 = tsh - tska(nsk)
      jt1 = floor(tsh1/30)
      dthskc = (pexsk(jt1+1,2)-pexsk(jt1,2))/30
      θ_skc = pexsk(jt1,2)+(tsh1-30*jt1)*dthskc
      dφskc = pexsk(jt1+1,3)-pexsk(jt1,3)
      if (dφskc >  PI)
        dφskc = dφskc - 2π
      end
      if (dφskc < -PI)
        dφskc = dφskc + 2π
      end
      dφskc = dφskc / 30
      φskc = pexsk(jt1,3) + (tsh1-30*jt1) * dφskc
      dgmsk = (pexsk(jt1+1,7)-pexsk(jt1,7)) / 30
      gmsk = pexsk(jt1,7) + (tsh1-30*jt1)*dgmsk
      rac = pska(nsk,4) / pska(nsk,6)
      ask(1) = rhclf * rskf * rac
      dask(1) = rhclf * vskf * rac
      rbc = pska(nsk,5) / pska(nsk,6)
      ask(2) = rhclf * rskf * rbc
      dask(2) = rhclf * vskf * rbc
      dtzskmin = (pexsk(jt1+1,8)-pexsk(jt1,8)) / 30
      tzskmin = pexsk(jt1,8) + (tsh1-30*jt1) * dtzskmin
    end
    tzskmin = cos(tzskmin)


    #tzskmin = -0.866
    # write(*,*) rskc, drskc, ask, dask
    #  θ_skc is colatitude
    sinθsk = sin(θ_skc)
    cosθsk = cos(θ_skc)
    sinφsk = sin(φskc)
    cosφsk = cos(φskc)
    mrtx(sinθsk, cosθsk, sinφsk, cosφsk, r2xsk)
    #  r2xsk matrix to convert r-the-φ of shock ellipsoid center to xyz of
    #    of HEEQ. r = z θ = x, φ = y
    mxptr(gmsk, xp2r)

    #xp2x = matmul(r2xsk, xp2r)
    for i = 1:3
      xp2x(i,1:3) = r2xsk(i,1)*xp2r(1,1:3) + r2xsk(i,2)*xp2r(2,1:3) + &
          r2xsk(i,3)*xp2r(3,1:3)
    end

    xskc(1) = rskc * sinθsk * cosφsk
    xskc(2) = rskc * sinθsk * sinφsk
    xskc(3) = rskc * cosθsk

    xc(1:3) = x(1:3) - xskc(1:3)
    xpc(1:3) = xc(1)*xp2x(1,1:3)+xc(2)*xp2x(2,1:3)+xc(3)*xp2x(3,1:3)

    #  distance to the shock
    grf1(1:3) = xpc(1:3)/ask(1:3)
    rho = norm2(grf1(1:3))
    grf2(1:3) = grf1(1:3) / ask(1:3)
    grf3(1:3) = grf1(1:3) * grf1(1:3)
    grf2m = norm2(grf2(1:3))
    dnsk = (rho-1) / grf2m

    if dnsk0*dnsk < 0. # crossed shock
      if grf1(3) < tzskmin #No shock on the backside
        vsk = 1d-6
        vnx = 0.57735026918962576450d0
        return
      end
      #   shock normal
      vnxp = grf2/grf2m
      vnx = matmul(xp2x, vnxp)
      #  calculate shock velocity (normal)
      dmrtx(sinθsk, cosθsk, sinφsk, cosφsk, dthskc, dφskc, dr2xsk)
      dxp2r = dmxptr(gmsk, dgmsk)
      for i = 1:3
        for j = 1:3
          dxp2x(i,j) = dot_product(dr2xsk(i,1:3), xp2r(1:3,j)) + &
                       dot_product(r2xsk(i,1:3), dxp2r(1:3,j))
        end
      end
      dot1 = dot_product(grf2, vnxp)
      dot2 = dot_product(grf3, dask(:)/ask(:))
      dxskc[1] = drskc * sinθsk * cosφsk &
               + rskc * cosθsk * dthskc * cosφsk &
               - rskc * sinθsk * sinφsk * dφskc
      dxskc[2] = drskc * sinθsk * sinφsk &
               + rskc * cosθsk * dthskc * sinφsk &
               + rskc * sinθsk * cosφsk * dφskc
      dxskc[3] = drskc * cosθsk - rskc * sinθsk * dthskc
      for i = 1:3
        dxpskc[i] = dot_product(dxskc, xp2x(:,i))
      end
      dot3 = dot_product(grf2, dxpskc)
      for i = 1:3
        drxc[i] = dot_product(xc, dxp2x(1:3,i))
      end
      dot4 = dot_product(grf2(1:3), drxc(1:3))
      vsk = (dot2 + dot3 - dot4) / dot1
    else
      vsk = 1.e-6
      vnx = 0.57735026918962576450e0
    end
    #write(60,"(f14.6,11(1pe12.4))") tsh, ra, rskc, vsk, vnx, ask, rskf, vskf

    return (dnsk=dnsk, vsk=vsk, vnx=vnx)
end

function locate(xx, n, x, j)
end
end