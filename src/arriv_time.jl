push!(LOAD_PATH, @__DIR__)
using epv: e2v

function main()

    rnz = 2
    rnm = 4

    println("e = ?(GeV)")
    e = parse(Float64, readline())
    v = e2v(e)
    println("v = $v AU/day")
    println("arriving time = $(1.0/v) days")
end

main()
