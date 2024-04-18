#!/usr/bin/env julia

push!(LOAD_PATH, @__FILE__)
using .datetime_utils: caldate

function main()
  #  Source injection at shocks
  #  time backward simulation to calculate fluxes at locations in IP
  #  pfss magnetic field model
  #  calculation is done in corotation reference frame
  #  pitch angle with outward magnetic field line
  #  pitch angle diffusion (symmetric D_{\mu\mu})
  #  perpendicular diffusion added

  NM1 = 16
  NMXID = 40

  rpb = zeros(5)
  rp0 = zeros(5)
  rp0org = zeros(5)
  r0 = zeros(3)
  rb = zeros(3)
  x0 = zeros(6)


  # USES common block specie
  # USES common block nsucmin
  # USES common block npmax
  # USES common block fbcnst
  # USES common block ndpdt
  # USES common block svsp
  # USES common block tmprm
  # USES common block ldptcl
  # USES common block dir
  # USES common block nseeds
  # USES common block bpark
  # USES common block srcmod
  # USES common block ndmumu

  df0 = zeros(3)
  ddf0 = zeros(3,3)
  bv0 = zeros(3)
  cvtu = zeros(3)
  gbmag = zeros(3)
  bxgb2 = zeros(3)
  gb1s = zeros(3)
  vnr = zeros(3)
  vnx = zeros(3)
  densw0 = 66410.0 #332820.0
  k4ok2 = 12.42420
  k6ok2 = 242.42420
end

main()
