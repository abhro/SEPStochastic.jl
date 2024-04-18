# this file is only useful for building documentation, and has no bearing on
# running the code of the project

using Documenter

push!(LOAD_PATH, "../src/")

makedocs(sitename="SEPStochastic.jl documentation", remotes=nothing)
