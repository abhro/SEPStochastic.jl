module file_op
using FortranFiles: FortranFile
using OffsetArrays: OffsetArray, OffsetMatrix, zeros
using DelimitedFiles: readdlm, writedlm
using HDF5

push!(LOAD_PATH, @__DIR__)
include("param.jl")
using .param: N_R, N_θ, N_φ, bgrid, gbgrid, b1rsgrid, R_DISCRETE, θ_DISCRETE, φ_DISCRETE


"populate g, h"
function read_shtc(shtc_filename, n)::Tuple{OffsetMatrix{Float64},OffsetMatrix{Float64}}
    g = zeros(Float64, 0:n, 0:n)
    h = zeros(Float64, 0:n, 0:n)

    infile = open(shtc_filename)

    for l = 0:n
        for m = 0:l
            line = readline(infile)
            arr = readdlm(IOBuffer(line))
            g[l,m] = arr[3]
            h[l,m] = arr[4]
        end
    end

    close(infile)

    return (g, h)
end

write_maggrid(maggrid_filename, bgrid, gbgrid) =
    write_maggrid(maggrid_filename, bgrid, gbgrid, Val(:formatted))
function write_maggrid(
        maggrid_filename,
        bgrid::OffsetArray{Float64,4},
        gbgrid::OffsetArray{Float64,4},
        ::Val{:formatted}
    )

    outfile = open(maggrid_filename, "w")
    for i = 0:N_R
        r = R_DISCRETE[i]
        for j = 0:N_θ
            θ = θ_DISCRETE[j]
            for k = 0:N_φ
                φ = φ_DISCRETE[k]
                position = [r, θ, φ]
                writedlm(outfile,
                         vcat(position, bgrid[i,j,k,:], gbgrid[i,j,k,:])')
            end
        end
    end
    close(outfile)
end
function write_maggrid(
        maggrid_filename,
        bgrid::OffsetArray{Float64,4},
        gbgrid::OffsetArray{Float64,4},
        ::Val{:unformatted}
    )

    outfile = FortranFile(maggrid_filename, "w")
    for i = 0:N_R
        r = R_DISCRETE[i]
        for j = 0:N_θ
            θ = θ_DISCRETE[j]
            for k = 0:N_φ
                φ = φ_DISCRETE[k]
                position = [r, θ, φ]
                #vec2write = [bgrid[i,j,k,:], gbgrid[i,j,k,:]]'
                write(outfile,
                      position,
                      bgrid[i,j,k,:],
                      #cvgrid[i,j,k,:],
                      gbgrid[i,j,k,:])
            end
        end
    end
    close(outfile)
end
function write_maggrid(
        maggrid_filename,
        bgrid::OffsetArray{Float64,4},
        gbgrid::OffsetArray{Float64,4},
        ::Val{:HDF5}
    )

    outfile = h5open(maggrid_filename, "w")

    # the three vectors are in reverse order because the first one has index
    # changes the fastest
    positions = Iterators.product(φ_DISCRETE, θ_DISCRETE, R_DISCRETE) |>
                    collect |>
                    vec .|>
                    reverse # now each position is (φ, θ, r) so flip it

    write(outfile, "position", positions)
    write(outfile, "bgrid", bgrid)
    write(outfile, "gbgrid", gbgrid)

    close(outfile)
end

read_maggrid(maggrid_filename) = read_maggrid(maggrid_filename, Val(:formatted))
function read_maggrid(
        maggrid_filename,
        ::Val{:formatted}
    )::Tuple{OffsetArray{Float64,4}, OffsetArray{Float64,4}}

    magfieldgrid = zeros(Float64, 0:N_R, 0:N_θ, 0:N_φ, 3)
    gbgrid = zeros(Float64, 0:N_R, 0:N_θ, 0:N_φ, 3)

    maggrid_file = open(maggrid_filename)

    for i = 0:N_R
        for j = 0:N_θ
            for k = 0:N_φ
                line = readline(maggrid_file)
                arr = readdlm(IOBuffer(line))
                #position = arr[1:3]
                magfieldgrid[i,j,k,:] = arr[4:6]
                gbgrid[i,j,k,:] = arr[7:9]
            end
        end
    end

    close(maggrid_file)

    return (magfieldgrid, gbgrid)
end
function read_maggrid(
        maggrid_filename, ::Val{:unformatted}
    )::Tuple{OffsetArray{Float64,4}, OffsetArray{Float64,4}}

    magfieldgrid = zeros(Float64, 0:N_R, 0:N_θ, 0:N_φ, 3)
    gbgrid = zeros(Float64, 0:N_R, 0:N_θ, 0:N_φ, 3)
    position = zeros(3)

    maggrid_file = FortranFile(maggrid_filename)

    for i = 0:N_R
        for j = 0:N_θ
            for k = 0:N_φ
                read(maggrid_file,
                     position, magfieldgrid[i,j,k,:], gbgrid[i,j,k,:])
            end
        end
    end

    close(maggrid_file)

    return (magfieldgrid, gbgrid)
end # function
function read_maggrid(
        maggrid_filename, ::Val{:HDF5}
    )::Tuple{OffsetArray{Float64,4}, OffsetArray{Float64,4}}

    maggrid_file = h5open(maggrid_filename)

    magfieldgrid = read(maggrid_file, "bgrid")
    gbgrid = read(maggrid_file, "gbgrid")

    magfieldgrid = OffsetArray(magfieldgrid, 0:N_R, 0:N_θ, 0:N_φ, 1:3)
    gbgrid = OffsetArray(gbgrid, 0:N_R, 0:N_θ, 0:N_φ, 1:3)

    close(maggrid_file)
    return (magfieldgrid, gbgrid)
end


read_b1rs() = read_b1rs(b1rs_filename, Val(:unformatted))
function read_b1rs(b1rs_filename, ::Val{:unformatted})
    b1rsgrid = zeros(Float64, 0:N_R, 0:N_θ, 0:N_φ)
    b1rsfile = FortranFile(b1rs_filename)

    for i = 0:N_R
        for j = 0:N_θ
            for k = 0:N_φ
                read(b1rsfile, b1rsgrid[i,j,k])
            end
        end
    end

    close(b1rsfile)

    return b1rsgrid
end
function read_b1rs(b1rs_filename, ::Val{:formatted})
    b1rsgrid = zeros(Float64, 0:N_R, 0:N_θ, 0:N_φ)
    b1rsfile = open(b1rs_filename)

    for i = 0:N_R
        for j = 0:N_θ
            for k = 0:N_φ
                b1rsgrid[i,j,k] = read(b1rsfile, Float64)
            end
        end
    end

    close(b1rsfile)

    return b1rsgrid
end

write_b1rs(b1rs_filename, b1rs::OffsetArray{Float64,4},
           map::OffsetArray{Float64,3}) =
    write_b1rs(b1rs_filename, b1rs, map, Val(:formatted))
function write_b1rs(b1rs_filename, b1rs::OffsetArray{Float64,4},
                    map::OffsetArray{Float64,3}, ::Val{:formatted})
    b1rs_outfile = open(b1rs_filename)
    for k = 0:N_R
        for i = 0:N_θ
            for j = 0:N_φ
                if b1rs[k,i,j,1] == 0.0
                    nmap = abs(map[k,i,j])
                    lφ = mod(nmap, N_φ+1)
                    nmap = floor(nmap / (N_φ+1))
                    lθ = mod(nmap, N_θ+1)
                    lr = floor(nmap / (N_θ+1))
                    b1rs[k,i,j,1] = sign(map[k,i,j]) * abs(b1rs[lr,lθ,lφ,1])
                end

                write(b1rs_outfile, b1rs[k,i,j,1]/b1rs[k,i,j,2])
            end
        end
    end
    close(b1rs_outfile)
end
function write_b1rs(b1rs_filename, b1rs::OffsetArray{Float64,4},
                    map::OffsetArray{Float64,3}, ::Val{:unformatted})

    b1rs_outfile = FortranFile(b1rs_filename)
    for k = 0:N_R
        for i = 0:N_θ
            for j = 0:N_φ
                if b1rs[k,i,j,1] == 0.0
                    nmap = abs(map[k,i,j])
                    lφ = mod(nmap, N_φ+1)
                    nmap = floor(nmap / (N_φ+1))
                    lθ = mod(nmap, N_θ+1)
                    lr = floor(nmap / (N_θ+1))
                    b1rs[k,i,j,1] = sign(map[k,i,j]) * abs(b1rs[lr,lθ,lφ,1])
                end

                write(b1rs_outfile, b1rs[k,i,j,1]/b1rs[k,i,j,2])
            end
        end
    end
    close(b1rs_outfile)
end

function prepareptcl(ptcl_config_file::AbstractString, dir::AbstractString)

    np::Int32
    nf::Int32
    tf = zeros(NFMAX)
    rf = zeros(3, NFMAX)
    ef = zeros(NFMAX)
    rmuf = zeros(NFMAX)

    sp0::Float64
    gp::Float64
    ap::Float64
    trgtf = zeros(4)
    trgtfs = zeros(4)
    sp::Float64
    scanw::Float64
    h0::Float64

    master = parsefile(ptcl_config_file)
    s = specie(rnz = master["rnz"], rnm = master["rnm"])

    sp0 = master["artificial-drift"]["sp0"]
    gp  = master["artificial-drift"]["gp"]
    ap  = master["artificial-drift"]["ap"]

    tdl = master["tdl"]

    #ptcl_file = FortranFile(dir * "loadptcl.dat")
    ptcl_file = open(dir * "loadptcl.dat")
    readline(ptcl_file) # sp0, gp, ap header line
    readline(ptcl_file, sp0, gp, ap)
    readline(ptcl_file) # header line for tdl
    readline(ptcl_file, tdl)
    readline(ptcl_file) # another header line (nz, nm)
    readline(ptcl_file, rnz, rnm)
    readline(ptcl_file) # you guessed it, header (nf, np)
    readline(ptcl_file, nf, np)
    readline(ptcl_file) # headers yay (ime, ...)
    for i = 1:nf
        readline(ptcl_file, tf[i], rf[1:3,i], ef[i], rmuf[i])
    end

    close(ptcl_file)


    return (
            specie=s,
           )
end
end
