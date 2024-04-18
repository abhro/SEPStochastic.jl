#!/usr/bin/env julia
"""
This program populates the seed file used by sim3d, specified in
\$SEEDS_FILENAME (usually seeds.dat)
"""

using Random: MersenneTwister, rand
using Comonicon: @main

push!(LOAD_PATH, @__FILE__)
include("datetime_utils.jl")
using .datetime_utils: seconds_of_year

"""
Generate seeds for use in simulation.

# Options

- `--seeds-count`: Number of seeds to generate (default 200).
- `--seeds-filename`: File to write the seeds to (hyphen is standard output).
"""
@main function seedgen(; seeds_count::Int64=200, seeds_filename::AbstractString)
    iseed = 1 * seconds_of_year()
    generator = MersenneTwister(iseed)

    nran::Vector{Int32} = []

    seeds_outfile = seeds_filename == "-" ? stdout : open(seeds_filename, "w")
    for i = 1:seeds_count
        local cur_ran
        while true
            cur_ran = -rand(generator, Int32)
            if cur_ran ∉ nran
                break
            end
        end
        push!(nran, cur_ran)
        println(seeds_outfile, cur_ran)
    end
    close(seeds_outfile)
end
