mol load psf 06-final.psf pdb 06-final.pdb
set sel [atomselect top "protein"]
set all [atomselect top all]
$all set beta 0
$sel set beta 1
$all writepdb const.pdb
