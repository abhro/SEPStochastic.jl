#!/usr/bin/env julia

using Dates: now, datetime2unix
using LinearAlgebra: LowerTriangular, norm, dot
using OffsetArrays: OffsetMatrix, zeros
using SpecialFunctions: loggamma
using LegendrePolynomials: Plm, collectPlm! # associated Legendre polynomial
using Comonicon: @main

push!(LOAD_PATH, @__DIR__)
include("file_op.jl")
using .file_op: read_shtc, write_maggrid
include("param.jl")
using .param: N_R, N_θ, N_φ


# number of spherical harmonic transform coefficients in ℓ
const NS = 40

# source surface radius, in solar radius units
const RSS::Float64 = 2.5

@doc raw"""Vector of lower-triangular matrices? Where each matrix is of the form
```math
\begin{bmatrix}
P_0^0    & \cdot    & \cdot    & \cdots & \cdot       & \cdot       & \cdot  \\
P_1^0    & P_1^1    & \cdot    & \cdots & \cdot       & \cdot       & \cdot  \\
P_2^0    & P_2^1    & P_2^2    & \cdots & \cdot       & \cdot       & \cdot  \\
\vdots   & \vdots   & \vdots   & \ddots & \vdots      & \vdots      & \vdots \\
P_{38}^0 & P_{38}^1 & P_{38}^2 & \cdots & P_{38}^{38} & \cdot       & \cdot  \\
P_{39}^0 & P_{39}^1 & P_{39}^2 & \cdots & P_{39}^{38} & P_{39}^{39} & \cdot  \\
P_{40}^0 & P_{40}^1 & P_{40}^2 & \cdots & P_{40}^{38} & P_{40}^{39} & P_{40}^{40}
\end{bmatrix}
```
where each ``P_ℓ^m`` is a matrix evaluated at an ``x`` given by the vector's
index, and each ``⋅`` signifies a zero entry in the matrix.
"""
#const ALP_TABLE = zeros(Float64, 0:N_θ, 0:NS, 0:NS)
const ALP_TABLE = zeros(Float64, 0:N_θ, 0:NS+2, 0:NS+2)
## first derivative
const DALP_TABLE = zeros(Float64, 0:N_θ, 0:NS, 0:NS)
## second derivative
const D2ALP_TABLE = zeros(Float64, 0:N_θ, 0:NS, 0:NS)

# for Schmidt semi-normalization
const FACTORIAL_TABLE = LowerTriangular(zeros(Float64, NS, NS))

@doc raw"""
    magfield(rvec, j, g, h) ->

Calculate the magnetic potential with the formula
```math
\Phi(\mathbf{r}) = \Phi(r, \theta, \phi) = R_0 \sum_{\ell=0}^m \sum_{m=0}^\ell
    P_\ell^m(\cos\theta)
    (g_{\ell m} \cos m\phi + h_{\ell m} \sin m \phi)
    \left[
        \left(\frac{R_0}{r}\right)^{\ell+1}
        -
        \left(\frac{R_0}{R_{ss}}\right)^{\ell+1}
        \left(\frac{r}{R_{ss}}\right)
    \right]
    \left[
        \ell + 1 + \ell \left(\frac{R_0}{R_{ss}}^{2\ell+1}\right)
    \right]^{-1}
```
where instead of going to infinity, we take ℓ to go up to 40 instead.

The magnetic field can be calculated as
```math
\mathbf{B}(r, \theta, \phi) = -\boldsymbol{\nabla} \Phi(r, \theta, \phi)
```
where each of the components are given by
```math
\begin{align*}
B_r(r, \theta, \phi) &= -\frac{\partial \Phi}{\partial r} \\
&= R_0 \sum_{\ell=0}\infty \sum_{m=0}^\ell
    P_\ell^m(\cos\theta)
    (g_{\ell m} \cos m\phi + h_{\ell m} \sin m\phi)
    \left[
        \frac{\ell + 1}{r} \left(\frac{R_0}{r}\right)^{\ell + 1}
        +
        \frac{\ell}{r} \left(\frac{R_0}{R_{ss}}\right)^{\ell+1}
        \left(\frac{r}{R_ss}\right)^\ell
    \right]
    \left[
        \ell + 1 + \ell \left(\frac{R_0}{R_ss}\right)^{2\ell + 1}
    \right]^{-1}
\\[2ex]
B_\theta(r, \theta, \phi) &= -\frac{1}{r} \frac{\partial Phi}{\partial \theta} \\
&= R_0 \sum_{\ell=0}\infty \sum_{m=0}^\ell
    \left.\frac{d P_\ell^m(x)}{dx}\right|_{x = \cos\theta}
    \sin\theta
    (g_{\ell m} \cos m\phi + h_{\ell m} \sin m\phi)
    \frac{1}{r}
    \left[
        \left(\frac{R_0}{r}\right)^{\ell+1}
        -
        \left(\frac{R_0}{R_{ss}}\right)^{\ell+1}
        \left(\frac{r}{R_{ss}}\right)
    \right]
    \left[
        \ell + 1 + \ell \left(\frac{R_0}{R_ss}\right)^{2\ell + 1}
    \right]^{-1}
\\[2ex]
B_\phi(r, \theta, \phi) &= -\frac{1}{r\sin\theta} \frac{\partial \Phi}{\partial \phi} \\
&= - R_0 \sum_{\ell=0}\infty \sum_{m=0}^\ell
    P_\ell^m(\cos\theta)
    \frac{m}{r\sin\theta}
    (g_{\ell m} \sin m\phi - h_{\ell m} \cos m\phi)
    \left[
        \left(\frac{R_0}{r}\right)^{\ell+1}
        -
        \left(\frac{R_0}{R_{ss}}\right)^{\ell+1}
        \left(\frac{r}{R_{ss}}\right)
    \right]
    \left[
        \ell + 1 + \ell \left(\frac{R_0}{R_ss}\right)^{2\ell + 1}
    \right]^{-1}
\end{align}
```
"""
function magfield(
        rvec::Vector{Float64},    # 3-vector, rvec = [r, θ, φ]
        j::Integer, # index of θ, or, equivalently, θ in degrees
        #b::Vector{Float64},     # 3-vector
        #bmag::Float64,
        #cvtu::Vector{Float64},  # 3-vector
        #∇normb::Vector{Float64}, # 3-vector
        #dbbds::Float64,
        g::OffsetMatrix{Float64},     # indices = (0:NS, 0:NS)
        h::OffsetMatrix{Float64}      # indices = (0:NS, 0:NS)
    )::@NamedTuple{b::Vector{Float64}, # mag field at (r, θ, φ)
                   bmag::Float64,      # norm of mag field at (r, θ, φ)
                   cvtu::Vector{Float64}, # curvature??
                   ∇normb::Vector{Float64},
                   dbbds::Float64}      #

    length(rvec) == 3 || throw(DomainError("rvec must be a 3-vector"))

    # PHYSICAL THEORY
    # Each ℓ,m term of the potential can be written as a product of three
    # functions depending uniquely on r, θ, φ each.
    # THEY ARE
    # F(r) = [(R_0/r)^(ℓ+1) - (R_0/R_ss)^(ℓ+1) (r/R_ss)^ℓ] / [ℓ+1 + ℓ(R_0/R_ss)^(2ℓ+1)]
    # G(θ) = P_ℓ^m(cos θ)
    # H(φ) = g_(ℓm) cos(mφ) + h_(ℓm) sin(mφ)

    # the derivatives are needed for various other functions. They are
    #
    # dG/dθ = - sinθ dP_ℓ^m(cosθ)/d(cosθ)
    # d²G/dθ² = sin²θ d²P_ℓ^m(cosθ)/d(cosθ)^2 - cosθ dP_ℓ^m(cosθ)/d(cosθ)
    #
    # dH/dφ = m [-g_(ℓm) sin(mφ) + h_(ℓm) cos(mφ)]
    # d²H/dφ² = -m² [g_(ℓm) cos(mφ) + h_(ℓm) sin(mφ)] = -m² H
    #
    # dF/dr = [-(ℓ+1) R_0^(ℓ+1)/r^(ℓ+2) - ℓ (R_0/R_ss)^(ℓ+1) r^(ℓ-1)/R_ss^ℓ]
    #         /
    #         [ℓ + 1 + ℓ (R_0/R_ss)^(2ℓ+1)]
    #       = [-(ℓ+1)/r (R_0/r)^(ℓ+1) - ℓ/r (R_0/R_ss)^(ℓ+1) (r/R_ss)^ℓ]
    #         /
    #         [ℓ + 1 + ℓ (R_0/R_ss)^(2ℓ+1)]
    # d²F/dr² = [(ℓ+1)(ℓ+2) R_0^(ℓ+1)/r^(ℓ+3) - ℓ(ℓ-1) (R_0/R_ss)^(ℓ+1) r^(ℓ-2)/R_ss^ℓ]
    #           /
    #           [ℓ + 1 + ℓ (R_0/R_ss)^(2ℓ+1)]
    # d²F/dr² = [(ℓ+1)(ℓ+2)/r² (R_0/r)^(ℓ+1) - ℓ(ℓ-1)/r² (R_0/R_ss)^(ℓ+1) (r/R_ss)^ℓ]
    #           /
    #           [ℓ + 1 + ℓ (R_0/R_ss)^(2ℓ+1)]

    # r and RSS are in units of solar radius (R_0)
    (r, θ, φ) = rvec

    # position vector outside the source surface. calculate it at
    # source surface for now and make adjustments later.
    if r > RSS
        r = RSS
    end

    cosθ = cos(θ)
    sinθ = sin(θ)
    if sinθ <= 1.745e-4
        sinθ = 1.745e-4
        cosθ = sqrt(1.0 - sinθ^2)
        if θ > 1.57
            cosθ = -cosθ
        end
    end

    alp = view(ALP_TABLE, j, :, :)

    b = zeros(Float64, 3)


    # Formulation of the Jacobian
    # The Jacobian for an R^3 is a 3^3 matrix.
    # For B, it is
    # JB(r, θ, φ) = [
    #       -∂_r^2 Φ              -∂_{θr}Φ / r            - ∂_{φr}Φ / (r sinθ)
    #       -∂_{rθ} Φ/r           -∂_θ^2 Φ / r^2          - ∂_{φr}Φ / (r^2 sinθ)
    #       -∂_{rφ} Φ/(r sinθ)    -∂_{θr}Φ / (r^2 sinθ)   - ∂_{φr}Φ / (r^2 sin^2 θ)
    # ]
    ∇b = zeros(Float64, 3, 3) # Jacobian

    # hold the b vector at each ℓ and m term,
    # if that makes sense (it's the summand)
    blm = zeros(Float64, 3)

    for ℓ = 1:NS
        r0rl2 = r^(-2-ℓ)
        rrs2l1 = (r / RSS)^(2ℓ+1)
        r0rs2l1 = ℓ + 1 + ℓ * RSS^(-2ℓ-1)
        for m = 0:ℓ
            # convert to Schmidt (semi-)normalized Associated Legendre Polynomial
            # XXX this should ideally be taken care of by init_aplm!, but we need
            # the raw values for now to calculate dalp and d2alp
            # TODO reformulate dalp and d2alp under normalization
            if m == 0
                schmt = 1
                dalp = alp[ℓ,1]
                d2alp = (-(ℓ+1)*ℓ*alp[ℓ,0] + alp[ℓ,2]) / 2.0
            else
                schmt = (-1)^m * FACTORIAL_TABLE[ℓ,m] # condon shortley phase
                dalp = (-(ℓ+m)*(ℓ-m+1)*alp[ℓ,m-1] + alp[ℓ,m+1]) / 2.0
                if m == 1
                    d2alp = (alp[ℓ,3] - (3*ℓ^2 - 3*ℓ - 2) * alp[ℓ,1]) / 4.0
                else
                    d2alp = (
                             alp[ℓ,m+2]
                             - ((ℓ+m+1)*(ℓ-m)+(ℓ+m)*(ℓ-m+1)) * alp[ℓ,m]
                             + (ℓ+m)*(ℓ-m+1)*(ℓ+m-1)*(ℓ-m+2) * alp[ℓ,m-2]
                            ) / 4.0
                end
            end

            cosmφ = cos(m*φ)
            sinmφ = sin(m*φ)

            ghmφ = g[ℓ,m]*cosmφ + h[ℓ,m]*sinmφ

            # derivative with respect to φ
            dghmφ = m * (-g[ℓ,m]*sinmφ + h[ℓ,m]*cosmφ)

            # r component: B_r =
            blm[1] =  r0rl2 * ghmφ  * schmt * alp[ℓ,m] * (ℓ+1+ℓ*rrs2l1) / r0rs2l1
            # θ component: B_θ =
            blm[2] = -r0rl2 * ghmφ  * schmt * dalp * (1-rrs2l1) / r0rs2l1
            # φ component: B_φ = - F(r)?
            # XXX why is this always 0?
            blm[3] = -r0rl2 * dghmφ * schmt * alp[ℓ,m] / sinθ * (1-rrs2l1) / r0rs2l1

            @. b += blm

            if blm[1] != 0
                ∇b[1,1] += -blm[1] * ((ℓ+2)*(ℓ+1)-ℓ*(ℓ-1)*rrs2l1) / (ℓ+1+ℓ*rrs2l1) / r
                ∇b[2,1] += blm[1] / alp[ℓ,m] * dalp
                ∇b[3,1] += blm[1] / ghmφ * dghmφ
            end
            if blm[2] != 0
                ∇b[1,2] += -blm[2] * (ℓ+2+(ℓ-1)*rrs2l1) / (1-rrs2l1) / r
                ∇b[2,2] += blm[2] / dalp * d2alp
                ∇b[3,2] += blm[2] / ghmφ * dghmφ
            end
            if blm[3] != 0
                ∇b[1,3] += -blm[3] * (ℓ+2+(ℓ-1)*rrs2l1) / (1-rrs2l1) / r
                ∇b[2,3] += blm[3] / alp[ℓ,m] * dalp - blm[3] * cosθ / sinθ
                ∇b[3,3] += blm[3] / dghmφ * m * m * ghmφ
            end
        end
    end

    local bmag = norm(b)

    local ∇normb = [dot(b, ∇b[1,:]) / bmag,
                    dot(b, ∇b[2,:]) / (bmag * r),
                    dot(b, ∇b[3,:]) / (bmag * r * sinθ)]
    if any(isnan.(∇normb))
        @error("∇normb is nan",
               r, θ/π, φ/π,
               bmag, # turns out to be 0, and the cause of all problems
               j,
               sinθ, cosθ
               #∇normb[1], ∇normb[2], ∇normb[3],
               #b[1], b[2], b[3]
              )
    end
    local dbbds = dot(b, ∇normb) / bmag^2

    local cvtu = zeros(Float64, 3)
    cvtu[1] = (
               b[1] * ∇b[1,1]
               + b[2] * ∇b[2,1] / r
               + b[3] * ∇b[3,1] / r/sinθ
               - (b[2]^2 + b[3]^2) / r
              )
    cvtu[2] = (
               b[1] * ∇b[1,2]
               + b[2] * ∇b[2,2]/r
               + b[3] * ∇b[3,2]/r/sinθ
               + (b[2]*b[1] - b[3]^2 * cosθ/sinθ) / r
              )
    cvtu[3] = (
               b[1] * ∇b[1,3]
               + b[2] * ∇b[2,3]/r
               + b[3] * ∇b[3,3]/r/sinθ
               + (b[3]*b[1] + b[3]*b[2]*cosθ/sinθ) / r
              )
    @. cvtu /= bmag

    cvtu .= (cvtu .- dbbds.*b) ./ bmag  # curvature perpendicular to B

    if rvec[1] > RSS
        # inverse square law once outside the source surface
        b[1] *= (RSS / rvec[1])^2

        bmag = norm(b)
        #cvtu[:] .= 0
        ∇normb[1] = -2 * bmag / rvec[1]
        ∇normb[2] *= (RSS/rvec[1])^3
        ∇normb[3] *= (RSS/rvec[1])^3
        dbbds = b[1] * ∇normb[1] / bmag^2
    end

    return (
            b=b,
            bmag=bmag,
            cvtu=cvtu,
            ∇normb=∇normb,
            dbbds=dbbds
           )
end

"""Populate the ALP_TABLE arrray for various cosθ, ℓ, and m"""
function init_aplm!(plm_table)
    #θs = Array(0:N_THETA) * π/180
    Threads.@threads for j = 0:N_θ
        θ = j * π/180
        cosθ = cos(θ)
        sinθ = sin(θ)
        if sinθ <= 1.745e-4
            sinθ = 1.745e-4
            cosθ = sqrt(1.0 - sinθ^2)
            if θ > π/2
                cosθ = -cosθ
            end
        end
        #  Associated Legendre Polynomial (really actually a triangular array)
        alp = zeros(Float64, 0:NS+2, 0:NS+2)
        for m = 0:NS+2
            collectPlm!(@view(alp[m:end,:]),
                        cosθ,
                        m=m, lmin=m, lmax=NS+2,
                        norm=Val(:standard),
                        csphase=false)
        end
        #for ℓ = 0:NS+2
        #    #alp[ℓ,0:ℓ] = Plm.(cosθ, ℓ, 0:ℓ, norm=Val(:schmidt))
        #    #alp[ℓ,0:ℓ] = Plm.(cosθ, ℓ, 0:ℓ, norm=Val(:schmidtquasi))
        #    alp[ℓ,0:ℓ] = Plm.(cosθ, ℓ, 0:ℓ, norm=Val(:standard), csphase=false)
        #    #for m = 0:ℓ
        #    #    p = Plm(cosθ, ℓ, m, norm=Val(:schmidt))
        #    #    alp[ℓ,m] = p
        #    #    if m == 0
        #    #        dalp[ℓ,m] = alp[ℓ,1]
        #    #    else
        #    #    end
        #    #end
        #end
        plm_table[j,:,:] .= alp
    end
end

function init_fac_table!(factorial_table)
    ns = size(factorial_table)[1]

    # we don't need the factor for ℓ = 0 (which only has m = 0) or any m = 0
    for ℓ = 1:NS
        for m = 1:ℓ
            # say, S = factorial_table
            # then, S_{ℓ,m} = sqrt(2 * (ℓ-m)! / (ℓ+m)!)
            # but using gamma directly would cause an overflow
            # LaTeX: S_{\ell,m} = \sqrt{2 \frac{(\ell - m)!}{(\ell + m)!}}
            factorial_table[ℓ,m] = sqrt(
                2 * exp(
                        # the actual ratio part (turns into subtraction for ln)
                        loggamma(ℓ-m+1) - loggamma(ℓ+m+1)
                       )
           )
        end
    end
end

"""
Create a magnetic field grid out of spherical harmonic transform coefficients.

# Options

- `--shtc-file`: path to NSO GONG SHTC file
- `--maggrid-file`: path to write maggrid to
- `--output-format`: Format to write maggrid file in (HDF5, unformatted, formatted)
"""
@main function maggrid(;shtc_file, maggrid_file, output_format)

    #chunk = 1

    g, h = read_shtc(shtc_file, NS)
    @debug("Finished loading shtc table")

    init_fac_table!(FACTORIAL_TABLE)
    @debug("Finished loading factorial table")

    bgrid = zeros(0:N_R,0:N_θ,0:N_φ,3)
    gbgrid = zeros(0:N_R,0:N_θ,0:N_φ,3)

    r = zeros(3)

    output_format = lowercase(output_format)
    write_format =
        if output_format == "hdf5"
            Val(:HDF5)
        elseif output_format == "unformatted"
            Val(:unformatted)
        else
            Val(:formatted)
        end

    init_aplm!(ALP_TABLE)
    @debug("Finished loading associated legendre polynomials")

    #println("Using pseudo-threading")
    #$OMP PARALLEL NUM_THREADS(5) DEFAULT(private) SHARED(chunk,bgrid,gbgrid,g,h,srfctrtl)
    #id = omp_get_thread_num()
    starttime = datetime2unix(now())
    @debug("Thread start time", starttime)

    flush(stdout)
    flush(stderr)
    #$OMP DO SCHEDULE(STATIC,chunk)
    Threads.@threads for i = 0:N_R
        r[1] = 1.0 + i/100.0
        for j = 0:N_θ
            r[2] = j * π/180.0
            for k = 0:N_φ
                r[3] = k * π/180.0

                #call magfield(r, b, bmag, cvtu, ∇normb, dbbds, srfctrtl, g, h)
                res = magfield(r, j, g, h)
                #dbbds = res.dbbds

                #  curvature not needed
                #cvgrid[i,j,k,:] = res.cvtu
                bgrid[i,j,k,:] = res.b
                gbgrid[i,j,k,:] = res.∇normb
            end
        end
        @info("Finished iteration i = $i")
        flush(stdout)
    end
    #!$OMP END PARALLEL
    @debug("Writing maggrid values")
    write_maggrid(maggrid_file, bgrid, gbgrid, write_format)

    stoptime = datetime2unix(now())
    @info("elapsed time:", stoptime - starttime)
end
