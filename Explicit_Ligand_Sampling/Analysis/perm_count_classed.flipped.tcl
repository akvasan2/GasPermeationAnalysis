################################################################################
# perm_count_classed.tcl by Eric Joon Shinn
################################################################################
# This script calculates the flow of particle across the membrane through the
# CP of AQP5.
################################################################################
#
#draw cylinder {0 0 -30} {0 0 30} radius 10

# Define the CP space
#                                
#             ______ ny             
#    .............................
#   /    nx/ (AQP5) /px         /|
#  ............................------- pz
#  |        ~~~~~~ py         | /
#  ............................------- nz

#
#   (0,3)       (0,3)      (0,3)
# +++++++++++++++++++++++++++++++   Progress value & origin value:
#           |          |                  (P, O)
#           |          |
#     0     |   (1,2)  |     0
#           |          |
#           |          |
# +++++++++++++++++++++++++++++++
#    +1    +1    +1    +1   +1
#
#
# Let's start by tracking downward flow
#
# frame i-1: Initiatialize. Set pre position on for above, +2
#
# frame i  : check position, compare to previous position
#            If particle moved from above and outside into CP, +1
#            If particle moved from CP to below, 0. Increment flow.
#            If particle moved from CP to out of CP
#            If particle moved from below into CP, no change.
#
#  Progress = occupancy ;  Ordinate = beta

set PWD "/Scr/akvasan2-new/AQP5/ELS/ELS_CO2_100"

set out_no_CPfile [open "Perm_noCP.flipped.dat" w]
set out_CPfile   [open "Perm_CP.flipped.dat"   w]

puts $out_no_CPfile "Frame, Flow count, Resid"
puts $out_CPfile "Frame, Flow count, Resid"

set bound_z_top  25
set bound_z_mid  0
set bound_z_bot -25
set cp_rad       30

if {1} {
set mol0 [mol new $PWD/Wrapping/not_water.psf];# add PSF
mol addfile $PWD/Wrapping/wrapped_no_water.dcd first 0 last -1 waitfor -1;# add trajectories
}

###### Flip every frame by 180 degrees ##########

set nf [molinfo $mol0 get numframes]

set all [atomselect top "all"]

for {set f 0} {$f<$nf} {incr f} {

        $all frame $f
        $all move [trans x 180]

}

# Initialize

set nf [molinfo top get numframes]
set aller [atomselect top all]
$aller set beta 0
set total_net_flux 0

set comp_3 [atomselect $mol0 "(resname CO2) and (name C)  and (z>$bound_z_top)" frame 0]
set comp_2 [atomselect $mol0 "(resname CO2) and (name C)  and (x**2+y**2<$cp_rad**2 and z<=$bound_z_top and z>$bound_z_bot)" frame 0]
set comp_1 [atomselect $mol0 "(resname CO2) and (name C)  and (x**2+y**2>=$cp_rad**2 and z<=$bound_z_top and z>$bound_z_bot)" frame 0]
set comp_0 [atomselect $mol0 "(resname CO2) and (name C)  and (z<=$bound_z_bot)" frame 0]

$comp_3 set beta 3
$comp_2 set beta 2
$comp_1 set beta 1
$comp_0 set beta 0

# remove ??
set aller [atomselect $mol0 "(resname CO2) and (name C)"]
$aller set occupancy 0
# remove ??

set comp_3_chk [atomselect $mol0 "(resname CO2) and (name C) and (z>$bound_z_top)" frame 0]
$comp_3_chk set occupancy 0

set comp_2_chk [atomselect $mol0 "(resname CO2) and (name C) and (x**2+y**2<$cp_rad**2 and z<=$bound_z_top and z>$bound_z_bot) and ((beta==3 and z>=$bound_z_mid) or (beta==2 and occupancy==1))" frame 0]
$comp_2_chk set occupancy 1

set comp_CP_chk [atomselect $mol0 "occupancy==1 and (same residue as (x^2+y^2)<121 and z>-10 and z<10)" frame 0]
$comp_CP_chk set occupancy 2

set comp_2_chk_no_prog [atomselect $mol0 "(resname CO2) and (name C) and (x**2+y**2<$cp_rad**2 and z<=$bound_z_top and z>$bound_z_bot) and (not ((beta==3 and z>=$bound_z_mid) or (beta==2 and occupancy>0)))" frame 0]
$comp_2_chk_no_prog set occupancy 0

set comp_1_chk [atomselect $mol0 "(resname CO2) and (name C) and (x**2+y**2>=$cp_rad**2 and z<=$bound_z_top and z>$bound_z_bot)" frame 0]
$comp_1_chk set occupancy 0

set comp_0_noCP_chk [ atomselect $mol0 "(resname CO2) and (name C) and (z<=$bound_z_bot) and (occupancy==1)" frame 0]
$comp_0_noCP_chk set occupancy 0

set comp_0_CP_chk [ atomselect $mol0 "(resname CO2) and (name C) and (z<=$bound_z_bot) and (occupancy==2)" frame 0]
$comp_0_CP_chk set occupancy 0

set comp_0_chk_no_prog [atomselect $mol0 "(resname CO2) and (name C) and (z<=$bound_z_bot) and occupancy==0" frame 0]
$comp_0_chk_no_prog set occupancy 0

set flow_cnt_CP 0
set flow_cnt_noCP 0

# Iterate

for {set i0 1} {$i0 < $nf} {incr i0} {
#           ^ we start at second frame, 1
    if {[file exist "stop"]} {
        error "Script stopped by user"
    }

    $comp_3_chk         frame $i0
    $comp_2_chk         frame $i0
    $comp_CP_chk        frame $i0
    $comp_2_chk_no_prog frame $i0
    $comp_1_chk         frame $i0
    $comp_0_noCP_chk    frame $i0
    $comp_0_CP_chk      frame $i0
    $comp_0_chk_no_prog frame $i0

    $comp_3_chk         update
    $comp_2_chk         update
    $comp_CP_chk        update
    $comp_2_chk_no_prog update
    $comp_1_chk         update
    $comp_0_noCP_chk    update
    $comp_0_CP_chk      update 
    $comp_0_chk_no_prog update

    $comp_3_chk         set occupancy 0
    $comp_2_chk         set occupancy 1
    $comp_CP_chk        set occupancy 2
    $comp_2_chk_no_prog set occupancy 0
    $comp_1_chk         set occupancy 0
    $comp_0_noCP_chk    set occupancy 0
    $comp_0_CP_chk      set occupancy 0
    $comp_0_chk_no_prog set occupancy 0

    ##### Flow count Not CP ######
    if {[llength [ $comp_0_noCP_chk get index] ] > 0} {
        set flow_cnt_noCP [expr {$flow_cnt_noCP+[llength [$comp_0_noCP_chk get index]]}]
        puts $out_no_CPfile "$i0, $flow_cnt_noCP, [lsort -unique [$comp_0_noCP_chk get resid]]"
    }

    ##### Flow count CP ######
    if {[llength [ $comp_0_CP_chk get index] ] > 0} {
        set flow_cnt_CP [expr {$flow_cnt_CP+[llength [$comp_0_CP_chk get index]]}]
        puts $out_CPfile "$i0, $flow_cnt_CP, [lsort -unique [$comp_0_CP_chk get resid]]"
    }
    

    $comp_3 frame $i0
    $comp_2 frame $i0
    $comp_1 frame $i0
    $comp_0 frame $i0

    $comp_3 update
    $comp_2 update
    $comp_1 update
    $comp_0 update

    $comp_3 set beta 3
    $comp_2 set beta 2
    $comp_1 set beta 1
    $comp_0 set beta 0

}

close $out_no_CPfile
close $out_CPfile
