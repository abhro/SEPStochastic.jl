function main()
    na = 0
    for i = 1:6
        read(stdin, Any)
    end
    fl::Float64 = 0
    fl2::Float64 = 0
    t = 0

    ar = zeros(Float64, 9)
    while t > -10
        t = read(stdin, Float64)
        read(stdin, ar)
        if !isnan(ar[7])
            na += 1
            fl += ar[7]
            fl2 += ar[7]^2
        end
    end

    println(fl/na)
    println(sqrt(fl2)/na)
end

main()
