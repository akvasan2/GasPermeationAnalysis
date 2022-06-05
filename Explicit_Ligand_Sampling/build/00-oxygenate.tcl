# Load scripts
source oxygenate.tcl

# Load system
set mol0 [mol new apo.psf]
mol addfile apo.pdb

oxygenate $mol0 "02-sys"
