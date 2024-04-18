#!/usr/bin/env/julia
using FortranFiles
using Printf: @printf
using DelimitedFiles: readdlm

#   build histogram

# km
RSS = 2.5
RS = 6.96e5

# raw km/s converted to Rs/min
VSW = 400.0 / RS * 60

# angular frequency of sun's rotation
OMEGA = 2π / (27.27 * 24 * 60)
K4OK2 = 12.4242
K6OK2 = 242.4242

al = zeros(8)
hist = zeros(91, 46)
hs = zeros(91, 46)
xmu::Float64 = 0.0
#aa:: = zeros(Int32, 0)
nn:: int


#f = open("spec.dat_e00")
#do i = 1:10
#    read(f)
#end
#close(f)

n = 0
for i = 1:1_569_000
    try
        linebuf = IOBuffer(readline())
        parsed = readdlm(linebuf)
        al = parsed[1:8]
        nn = parsed[end]
    catch
        break
    end

    if al[1] == -1000.00
        break
    end
    ra = al[3]
    tpsw = (
            (ra - rss)
            - rss * log(ra/rss)
            + K4OK2 * (1/rss/2 - 1/ra + rss/ra^2/2)
            + K6OK2 * (1/rss^3/12 - 1/ra^3/3 + rss/ra^4/4)
    )
    tpsw = tpsw / vsw
    th = al[4]
    ph = al[5] + OMEGA*tpsw
    ph = atan(sin(ph), cos(ph))
    println(th, ph)
    if ph < 0
        ph = ph + 2π
    end
    n += 1
    if ra >= 2.50
        xmu = cos(th)
        iy = floor(th / 2π * 90) + 1
        ix = floor(ph / 2π * 90) + 1
        hist[ix,iy] = hist[ix,iy] + 1
    end
end

#99 continue
#  smoothing
for i = 1:45
    for j = 1:90
        if hist[j,i] < 10 # use 5x5 box smoothing
            ns = 0
            nz = 0
            for ii in range(i-2, i+2)
                for jj in range(j-2, j+2)
                    if ii > 0 && ii < 46
                        jjj = jj
                        if jj < 1
                            jjj = 90 + jj
                        end
                        if jj > 90
                            jjj = jj - 90
                        end
                        ns += 1
                        if hist[jjj,ii] != 0
                            nz += 1
                        end
                        hs[j,i] += hist[jjj,ii]
                    end
                end
            end
        end

        if 3*nz > ns
            hs[j,i] = hs[j,i] / ns
        else
            hs[j,i] = hist[j,i]
        end

        if hist[j,i] > 10 && hist[j,i] < 100 # use 3x3 box
            ns = 0
            for ii in [i-1, i, i+1]
                for jj in [j-1, j, j+1]
                    if ii > 0 && ii < 46
                        jjj: int = jj
                        if jj < 1
                            jjj = 90 + jj
                        end
                        if jj > 90
                            jjj = jj - 90
                        end
                        ns = ns + 1
                        if hist[jjj,ii] != 0
                            nz += 1
                        end
                        hs[j,i] += hist[jjj,ii]
                    end
                end
            end
        end

        if 3*nz > ns
            hs[j,i] = hs[j,i] / ns
        else
            hs[j,i] = hist[j,i]
        end

        if hist[j,i] >= 100
            hs[j,i] = hist[j,i]
        end
    end
end

hist = hs

for i = 1:45
    for j = 1:90
        th = (i-0.5) / 90 * 2π
        ph = (j-0.5) / 90 * 360
        @printf("%8.3f %8.3f %12.4e %6.0f",
                th*360/2π, ph, hist[j,i]/n/sin(th), hist[j,i])
    end
end
