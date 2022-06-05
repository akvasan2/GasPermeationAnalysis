set mol0 [mol new 02-sys.psf]
mol addfile 02-sys.pdb

set mol1 [mol new 04-O2.psf]
mol addfile 04-O2.pdb

set sel0 [atomselect $mol0 all]
set sel1 [atomselect $mol1 all]

mol fromsels "$sel0 $sel1"
set saver [atomselect top all]

$saver writepsf 06-final.psf
$saver writepdb 06-final.pdb
