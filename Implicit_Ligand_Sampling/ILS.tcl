# Implicit Ligand Sampling run script
# ===================================

# Running ILS calculation:

# vmd -dispdev text -e <name_of_this_file>

# You will need VMD 1.8.7. or higher.

# ILStools plugin of the same or higher version than the one
# that generated it needs to be available for execution.
if {[catch {package require ilstools 1.5} pkgver]} {
    vmdcon -err "This script requires the ilstools plugin v1.5."
    exit 1
}
package require ilstools
# Change the input parameters below to your liking.
# The filenames used in this script are relative to the directory
# for which it was generated but you can of course change then.

# If you have a CUDA enabled Nvidia GPU VMD will use the GPU for
# the computation. Since the all the GPU resources will then be
# used for ILS your graphical display will freeze up, so don't be
# surprised. After finishing each frame the display will briefly
# be updated and freeze again.

# Comment this line out to prevent the use of the CUDA implementation:
set PWD "/Scr/akvasan2-new/AQP5/T41F_Sim"
set env(VMDCUDAILS) 1
# Set the device id of the GPU to be used with the CUDA ILS implementation:
set env(VMDILSCUDADEVICE) 0

# You might want to do a quick test with 1 frame first to see if
# the syntax is correct and to determine the approximate runtime
# per frame.

# Adjustable parameters:
# ----------------------

# First and last frames to process

#set first 0           THIS PORTION WAS COPIED BELOW FOR LATER SO WE COULD
#set last  25          GET THE LAST FROM BY ASKING VMD AFTER LOADING THE TRAJ.

# Resolution of final downsampled map in Angstrom
set res    1.0

# Subsampling of each dimension during computation
# i.e. each gridpoint of the final map will actually
# be downsampled from subres^3 points.
set subres 3

# Control of the angular spacing of probe orientation vectors,
# i.e. the number of probe conformers generated.
#
#   1: use 1 orientation only
#   2: use 6 orientations (vertices of octahedron)
#   3: use 8 orientations (vertices of hexahedron)
#   4: use 12 orientations (faces of dodecahedron)
#   5: use 20 orientations (vertices of dodecahedron)
#   6: use 32 orientations (faces+vert. of dodecahedron)
#  >6: geodesic subdivisions of icosahedral faces
#      with frequency 1, 2, ...
#
#   For each orientation a number of rotamers will be
#   generated. The angular spacing of the rotations
#   around the orientation vectors is chosen to be about
#   the same as the angular spacing of the orientation
#   vector itself.
#   If the probe ha at least one symmetry axis then the
#   rotations around the orientation vectors are reduced
#   accordingly. If there is an infinite oder axis (linear
#   molecule) the rotation will be omitted.
#   In case there is an additional perpendicular C2 axis
#   the half of the orientations will be ignored so that
#   there are no antiparallel pairs.
#
#   Probes with tetrahedral symmetry:
#   Here conf denotes the number of rotamers for each of
#   the 8 orientations defined by the vertices of the
#   tetrahedron and its dual tetrahedron.
set orient   6

# Cutoff energy above which the occupancy is regarded zero
# For GPUs energies of more than 87 always correspond to
# floating point values of zero for the occupancy. Hence
# there is no point going higher than that.
set maxen  50

# Temperature of the MD simualtion
set T   310

# Nonbonded interaction cutoff
set cutoff 12.0

# The minmax box defining the free energy map
# (two opposite corners of the grid)
set minmax {{-6 -6 -35} {6 6 35}}

if {0} {
# The DX file containing the free energy map
set dxfile 02-CenPoreILSCO2OUT.dx
}

# -------------------------------------------------------

# Set up the probe

# WARNING: The probe parameters have only been verified to
#          reproduce experimental solvation free energies for
#          xenon, oxygen, nitric oxide and carbon monoxide probes.
#          Use parameters for other probes as a starting point
#          for optimization.
#
#          For xenon, VDW radius = 2.24 & VDW episilon = -0.494
#          was used by Jordi Cohen, "Imagine the Migration Pathways
#          for O2, CO, No and Xe Inside Myoglobin," Biophysical Journal,
#          2006.

set pmol [mol new ilsprobe_carbondioxide.xyz]
set psel [atomselect $pmol all]
set c [atomselect $pmol "name C"]
$c set radius     2.1;  # VDW radius
$c set occupancy -0.11; # VDW epsilon
set o [atomselect $pmol "name O"]
$o set radius     1.7;  # VDW radius
$o set occupancy -0.12; # VDW epsilon

# -------------------------------------------------------

# This script depends on the ilstools package
package require ilstools

# Load molecule and trajectory
set molid [mol new $PWD/build/apo.psf]
mol addfile $PWD/build/apo.pdb
mol addfile "$PWD/Wrapping/Aligned.everyframe.dcd" step 1 waitfor -1

# First and last frames to process
set first 0
set last [molinfo top get numframes]

# Align
set aller [atomselect $molid all]
set prot_targ [atomselect $molid "name CA" frame 0]
set prot_from [atomselect $molid "name CA"]

for {set i0 0} {$i0 < $last} {incr i0} {
    $aller frame $i0
    $prot_from frame $i0
    $aller move [measure fit $prot_from $prot_targ]
}


# Set the radius and occupancy field for each atom to the
# VDW rmin and epsilon values from the force field
ILStools::readcharmmparams {toppar/par_all36m_prot.prm toppar/par_all36_lipid.prm toppar/toppar_water_ions.str toppar/toppar_ions_won.ils.str toppar/par_all36_na.prm toppar/par_all36_carb.prm toppar/par_all36_cgenff.prm toppar/par_interface.prm toppar/toppar_all36_moreions.str toppar/toppar_all36_nano_lig.str toppar/toppar_all36_nano_lig_patch.str toppar/toppar_all36_synthetic_polymer.str toppar/toppar_all36_synthetic_polymer_patch.str toppar/toppar_all36_polymer_solvent.str toppar/toppar_dum_noble_gases.str toppar/toppar_all36_prot_arg0.str toppar/toppar_all36_prot_c36m_d_aminoacids.str toppar/toppar_all36_prot_fluoro_alkanes.str toppar/toppar_all36_prot_heme.str toppar/toppar_all36_prot_na_combined.str toppar/toppar_all36_prot_retinol.str toppar/toppar_all36_prot_model.str toppar/toppar_all36_prot_modify_res.str toppar/toppar_all36_na_nad_ppi.str toppar/toppar_all36_na_rna_modified.str toppar/toppar_all36_lipid_sphingo.str toppar/toppar_all36_lipid_archaeal.str toppar/toppar_all36_lipid_bacterial.str toppar/toppar_all36_lipid_cardiolipin.str toppar/toppar_all36_lipid_cholesterol.str toppar/toppar_all36_lipid_dag.str toppar/toppar_all36_lipid_inositol.str toppar/toppar_all36_lipid_lnp.str toppar/toppar_all36_lipid_lps.str toppar/toppar_all36_lipid_mycobacterial.str toppar/toppar_all36_lipid_miscellaneous.str toppar/toppar_all36_lipid_model.str toppar/toppar_all36_lipid_prot.str toppar/toppar_all36_lipid_tag.str toppar/toppar_all36_lipid_yeast.str toppar/toppar_all36_lipid_hmmm.str toppar/toppar_all36_lipid_detergent.str toppar/toppar_all36_lipid_ether.str toppar/toppar_all36_carb_glycolipid.str toppar/toppar_all36_carb_glycopeptide.str toppar/toppar_all36_carb_imlab.str toppar/toppar_all36_label_spin.str toppar/toppar_all36_label_fluorophore.str}
ILStools::assigncharmmparams $molid

# Run ILS
set blocks 1
set ils_div [expr {$last/$blocks}]
for {set ils_it 0} {$ils_it < $blocks} {incr ils_it} {
    set first_frm [expr {$ils_it*$ils_div}]
    set last_frm  [expr {$first_frm+$ils_div}]

    # The DX file containing the free energy map
    #set dxfile [format 02-CenPoreILSCO2OUT$first_frm-$last_frm.dx]
    set dxfile "02-CenPoreILSCO2OUT$first_frm-$last_frm.dx"
    if {1} {
    volmap ils $molid $minmax -cutoff $cutoff \
        -res $res -subres $subres -probesel $psel -orient $orient \
        -maxenergy $maxen -T $T -first $first_frm \
        -last $last_frm \
        -o $dxfile
    }
}
# Quit VMD when done with ILS calculation
quit
