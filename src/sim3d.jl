using Comonicon: @main

push!(LOAD_PATH, @__FILE__)

const NM1 = 16
const NMXID = 40

using file_op: read_b1rs

@main function sim3d(; maggrid_file, b1rs_file, ptcl_file)

  rpb = zeros(5)
  rp0 = zeros(5)
  rp0org = zeros(5)
  r0 = zeros(3)
  rb = zeros(3)
  x0 = zeros(6)

  # USES common block specie
  # USES common block sptm
  # USES common block nsucmin
  # USES common block npmax
  # USES common block fbcnst
  # USES common block ndpdt
  # USES common block svsp
  # USES common block tmprm
  # USES common block ldptcl
  # USES common block dir
  # USES common block rankstr
  # USES common block nseeds
  # USES common block bpark
  # USES common block srcmod
  # USES common block vsksw
  # USES common block ndmumu


  densw0 = 166410.0 #332820.e0
  k4ok2 = 12.4242e0
  k6ok2 = 242.4242e0

  nodes = 1 # default, will be read in readparm from 'dir.dat'
  if nodes > NMXID
      error("nodes > NMXID")
  end
  chunk = 1

  # USES environment variables PARAM_OUTDIR_PATH PARAM_NODES
  read_param(nodes)

  # USES environment variables MAGGRID_FILE
  (bgrid, gbgrid) = read_maggrid(maggrid_file)
  # USES environment variables B1RS_FILE
  b1rsgrid = read_b1rs(b1rs_file)

  # USES environment variables PTCL_FILE
  prepareptcl()
end

main()
