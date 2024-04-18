module commonblocks
using OffsetArrays: OffsetMatrix

@kwdef struct specie
    rnz::Float64
    rnm::Float64
end

@kwdef struct bpark
    densw0::Float64
    vsw::Float64
    k4ok2::Float64
    k6ok2::Float64
    omega::Float64
    b1au::Float64
    vom::Float64
    facip::Float64
end

@kwdef struct srcmod
    sp::Float64
    sp0::Float64
    gp::Float64
    ap::Float64
    trgtfs::Vector{Float64} # 4-vector
    scanw::Float64
    h0::Float64
end

@kwdef struct fbcnst
    rb0::Float64
    rmax::Float64
    rk::Float64
    deltat::Float64
    tc::Float64
    tl::Float64
    tmodel0::Float64
    nfbconst::Int32
end

@kwdef struct vsksw
    vsksw::Float64
    tauf::Float64
    tauc1_0::Float64
    tauc2::Float64
    tauc2_0::Float64
    vcme0kmPs::Float64
    vcme0::Float64
end

@kwdef struct cmesk
    pexsk::OffsetMatrix{Float64} # = zeros(Float64, 0:144, 8)
    vskf0::Float64
    acsk::Float64
    tska::Vector{Float64} # = zeros(Float64, 20)
    pska::Matrix{Float64} # = zeros(Float64, 20, 8)
    nsk::Int32
end

@kwdef struct tmprm
    t0org::Float64
    te::Float64
    tdl::Float64
    dmapjul::Float64
    tcme0::Float64
end

@kwdef struct acc
    ecrmax::Float64
    ecrmin::Float64
    ob0::Float64
    obw::Float64
    etinj::Float64
    fampb::Float64
    ratkp::Float64
end

@kwdef struct ldptcl
    tf::Vector{Float64} # = zeros(NFMAX)
    rf::Matrix{Float64} # = zeros(3, NFMAX)
    ef::Vector{Float64} # = zeros(NFMAX)
    rmuf::Vector{Float64} # = zeros(NFMAX)
    np::Int32
    nf::Int32
end
end
