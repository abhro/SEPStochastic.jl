module param

using OffsetArrays: zeros, OffsetArray

# physical constants

const TAU = 2π

const CSPEED::Float64 = 25.8441774    # speed of light in R_sun/min
const QoMSI::Float64 = 9.57883376e7   # proton charge-to-mass ratio in coulomb/kg
const EP::Float64 = 0.938             # proton rest energy in GeV (mass in GeV/c^2)
const EE::Float64 = 0.5109989183e-3   # electron rest energy in GeV (mass in GeV/c^2)
const GAMMA::Float64 = 5 / 3          # heat capacity ratio


# program parameters

const NSEEDMAX::UInt64 = 2001
const NFMAX::UInt64 = 200

const NFL::UInt64 = 21
const NSTS::UInt64 = 25
const NNDS::UInt64 = 30

const NBASE::UInt64 = 40
const NSPMAX::UInt64 = 20

# for PFSS
const N_R::UInt32 = 150
const N_θ::UInt32 = 180
const N_φ::UInt32 = 360

# discretization
#    take 151 steps, each of length 0.01 R_⊙, from R_⊙ to 2.5 R_⊙ (= R_ss)
const R_DISCRETE = OffsetArray(range(1.0, 2.5, length=N_R+1), 0:N_R)
#    θ from 0 to π radians, in 181 (=180+1) steps. note that the index of each
#    element also corresponds to its degree value.
const θ_DISCRETE = OffsetArray(range(0, π, length=N_θ+1), 0:N_θ)
#    φ from 0 to 2π radians, in 361 (=360+1) steps. note that the index of each
#    element also corresponds to its degree value.
const φ_DISCRETE = OffsetArray(range(0, 2π, length=N_φ+1), 0:N_φ)

bgrid = zeros(0:N_R, 0:N_θ, 0:N_φ, 3)
#cvgrid = zeros(0:N_R, 0:N_θ, 0:N_PHI, 3)
gbgrid = zeros(0:N_R, 0:N_θ, 0:N_φ, 3)
b1rsgrid = zeros(0:N_R, 0:N_θ, 0:N_φ)

const EPSILON::Vector{Float64} = [0.04, 0.01, 0.04, 0.003]

# conversion factors
const rad_to_deg::Float64 = 180 / π
const deg_to_rad::Float64 = π / 180

end # module
