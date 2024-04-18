
using Printf: @printf
using LinearAlgebra: norm
using OffsetArrays: OffsetArray, zeros
using FortranFiles: FortranFile
using Comonicon: @main
using Dates: now, datetime2unix

push!(LOAD_PATH, @__FILE__)
include("param.jl")
include("mtrx.jl")
include("file_op.jl")
include("rksolvers.jl") # TODO do this better

using .param: N_R, N_θ, N_φ
using .mtrx: trilinear
using .file_op: read_maggrid, write_b1rs

"""
Velocity as function of position, time, and magnetic field
"""
function vfunc(
        _::Float64, # input variable (time), not necessary
        r::Vector{Float64},
        pol::Float64,
        magfieldgrid
    )::Vector{Float64}

    length(r) == 3 || throw(DomainError("Invalid length of (r = $r), should be 3"))

    @debug("Calling magfield with r = $r")
    b = magfield(r, magfieldgrid)
    @debug("Got back b = $b")

    bmag = norm(b)
    if bmag == 0
        return zeros(Float64, 3)
    end
    v = -pol * b / norm(b)
    @debug("Got v = $v")
    v[2] /= r[1]
    sinθ = sin(r[2])
    if sinθ == 0
        sinθ = 1.0e-6
    end
    v[3] /= r[1] * sinθ

    v
end # function


function magfield(
        r::Vector{Float64},
        magfieldgrid::OffsetArray{Float64,4}
    )::Vector{Float64}

    length(r) == 3 || throw(DomainError("Invalid length of (r = $r), should be 3"))

    RS = 2.5
    @debug("magfield called with r = $r")

    r = copy(r)
    r[2] = acos(cos(r[2]))
    r[3] = atan(sin(r[3]), cos(r[3]))
    if r[3] < 0.0
        r[3] = r[3] + 2π
    end

    #  find the grid cell
    irr = floor(Int64, (r[1]-1) / (RS-1) * N_R)
    if irr >= N_R
        irr = N_R - 1
    end
    #irr = max(0, irr) # don't let it go back into the sun

    iθ = floor(Int64, r[2] / π * N_θ)
    if iθ >= N_θ
        iθ = N_θ - 1
    end

    iφ = floor(Int64, r[3] / 2π * N_φ)
    if iφ >= N_φ
        iφ = N_φ - 1
    end

    if irr < 0
        return magfieldgrid[0,iθ,iφ,:]
    end

    #  relative displacement from lower grids
    px = [(r[1]-1.0) / (RS-1.0) * N_R - irr,
          r[2] / π * N_θ - iθ,
          r[3] / 2π * N_φ - iφ]

    φc = zeros(Float64, 2, 2, 2)
    b = zeros(Float64, 3)
    for m = 1:3
        φc[1,1,1] = magfieldgrid[irr,   iθ,   iφ,   m]
        φc[2,1,1] = magfieldgrid[irr+1, iθ,   iφ,   m]
        φc[1,2,1] = magfieldgrid[irr,   iθ+1, iφ,   m]
        φc[2,2,1] = magfieldgrid[irr+1, iθ+1, iφ,   m]
        φc[1,1,2] = magfieldgrid[irr,   iθ,   iφ+1, m]
        φc[2,1,2] = magfieldgrid[irr+1, iθ,   iφ+1, m]
        φc[1,2,2] = magfieldgrid[irr,   iθ+1, iφ+1, m]
        φc[2,2,2] = magfieldgrid[irr+1, iθ+1, iφ+1, m]
        b[m] = trilinear(φc, px)
        @debug("At m = $m, b = $b")
    end

    return b
end # function

@main function mapb2s(; maggrid_file, b1rs_file, maggrid_file_format, b1rs_file_format)
    b1rs = zeros(Float64, 0:N_R, 0:N_θ, 0:N_φ, 2)
    #v = zeros(3)
    r0 = zeros(Float64, 3)
    r1 = zeros(Float64, 3)
    r = zeros(Float64, 3)
    b = zeros(Float64, 3)
    #sinθ::Float64
    rmin = zeros(Float64, 3)
    rmin1 = zeros(Float64, 3)
    map = zeros(Int64, 0:N_R, 0:N_θ, 0:N_φ)

    #n1::Int64

    maggrid_file_format = lowercase(maggrid_file_format)
    read_format =
        if maggrid_file_format == "hdf5"
            Val(:HDF5)
        elseif maggrid_file_format == "unformatted"
            Val(:unformatted)
        else
            Val(:formatted)
        end

    b1rs_file_format = lowercase(b1rs_file_format)
    write_format =
        if b1rs_file_format == "hdf5"
            Val(:HDF5)
        elseif b1rs_file_format == "unformatted"
            Val(:unformatted)
        else
            Val(:formatted)
        end

    t = 0.0
    Δt = 0.005

    #  load magnetic field grid
    @info("Reading maggrid file")
    maggrid_reading_start_time = datetime2unix(now())
    magfieldgrid, gbgrid = read_maggrid(maggrid_file, read_format)
    maggrid_reading_end_time = datetime2unix(now())
    @info("Reading maggrid file took $(maggrid_reading_end_time - maggrid_reading_start_time) s")

    @info("Grid sizes", N_R, N_θ, N_φ)

    #!$OMP PARALLEL NUM_THREADS(20) DEFAULT(firstprivate) SHARED(magfieldgrid,b1rs,map)
    #!$OMP DO SCHEDULE(static,1)
    #  trace field lines from surface surface back to 1 Rs

    Threads.@threads for k = 0:N_R
        r0[1] = 1.000_010 + k*0.010
        for i = 0:N_θ
            r0[2] = i * π/180
            for j = 0:N_φ
                r0[3] = j * π / 180
                r = r0
                @debug("Getting magfield for", k, i, j)
                b = magfield(r, magfieldgrid)
                @debug("Got magfield b = $b")
                pol = sign(b[1])
                n = 0

                bmag = norm(b)
                b1rs[k, i, j, 2] = bmag
                if k == 1 && i == 1 && j == 1
                    println("b = $b")
                end
                rmin[1] = 100_000.0
                @debug("Stepping into while loop")
                while r[1] > 1.000_009
                    r[2] = acos(cos(r[2]))
                    r[3] = atan(sin(r[3]), cos(r[3]))
                    if r[3] < 0
                        r[3] += 2π
                    end
                    r1 = r
                    b = magfield(r, magfieldgrid)
                    bmag = norm(b)
                    #v = -pol*b/bmag
                    #v[2] /= r1[1]
                    #sinθ = sin(r1[2])
                    #if sinθ == 0
                    #    sinθ = 1.0e-6
                    #end
                    #v[3] /= r1[1] * sinθ
                    Δt = 0.001
                    @debug("Calling rk4")
                    #r = rk4(r1, v, t, Δt, vfunc, odefun_params=(pol, magfieldgrid))
                    r = rk4(vfunc, t, r1, Δt, odefun_params=(pol, magfieldgrid))
                    n += 1

                    if r[2] < 0
                        r[2] = -r[2]
                        r[3] = r[3] + pi
                        if r[3] > 2π
                            r[3] -= 2π
                        end
                    end

                    if r[1] < rmin[1]
                        rmin = r
                    end
                    if r[1] - rmin[1] > 2.5 || n > 10000 # wrong direction?
                        @debug("Wrong direction, reversing polarity")
                        pol = -pol # try back with opposite polarity
                        n1 = 0
                        r = r0
                        rmin1[1] = 100000.0
                        while r[1] > 1.000009
                            r[2] = acos(cos(r[2]))
                            r[3] = atan(sin(r[3]), cos(r[3]))
                            if r[3] < 0
                                r[3] += 2π
                            end
                            r1 = r
                            b = magfield(r, magfieldgrid)
                            bmag = norm(b)
                            #v = @. -pol * b / bmag
                            #v[2] /= r1[1]
                            #sinθ = sin(r1[2])
                            #if sinθ == 0
                            #    sinθ = 1.0e-6
                            #end
                            #v[3] /= r1[1] * sinθ
                            Δt = 0.001
                            #r = rk4(r1, v, t, Δt, vfunc, odefun_params=(pol, magfieldgrid))
                            r = rk4(vfunc, t, r1, Δt, odefun_params=(pol, magfieldgrid))
                            n1 += 1
                            if r[2] < 0
                                r[2] = -r[2]
                                r[3] = r[3] + π
                                if r[3] > 2π
                                    r[3] -= 2π
                                end
                            end
                            if r[1] < rmin1[1]
                                rmin1 = r
                            end
                            if (r[1] - rmin1[1]) > 2.5 || n1 > 10000 # double open field or stuck
                                if rmin[1] < rmin1[1]
                                    pol = -pol
                                    r = rmin # use rmin to remap
                                else
                                    pol = pol
                                    r = rmin1
                                end
                                lr = floor((r[1] - 1.0) / 0.01)
                                lθ = floor(r[2] / π * 180)
                                lϕ = floor(r[3] / π * 180)
                                if (((lr * (N_θ + 1) + lθ) * (N_φ + 1) + lϕ) <
                                    ((k * (N_θ + 1) + i) * (N_φ + 1) + j))
                                    map[k, i, j] = (lr * (N_θ + 1) + lθ) * (N_φ + 1) + lϕ
                                    bmag = 0
                                else
                                    println("New x point at $lr $lθ $lϕ")
                                    map[k, i, j] = ((lr - 1) * (N_θ + 1) + lθ) * (N_φ + 1) + lϕ
                                    bmag = 0
                                end # if
                                break
                            end # if
                        end # for
                        break
                    end # if
                end # while

                b1rs[k, i, j, 1] = pol * bmag
                map[k, i, j] *= pol
            end # for
        end # for
    end # for

    #$OMP BARRIER
    #$OMP END PARALLEL


    write_b1rs(b1rs_file, b1rs, map, write_format)

end # main
