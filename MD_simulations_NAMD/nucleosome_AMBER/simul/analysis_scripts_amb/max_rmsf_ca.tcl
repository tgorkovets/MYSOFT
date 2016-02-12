# Here we check maximum deviation of CA H3 core atoms
# in order to get reasonable constant for constraints in future simulations
# and also check it.

mol load psf ../analysis_data/only_nucl_init.psf

mol addfile ../analysis_data/only_nucl_init.pdb waitfor all
mol addfile ../analysis_data/md_nucl.dcd waitfor all

set alpha_core_ca "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"
 
set h3 [atomselect top "(segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) and name CA"]

#set h3_0 [atomselect top "(segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) and name CA" frame 0]


 set nframes [expr  [molinfo top get numframes] - 1 ]


 set outfile [open ../analysis_data/max_rmsf.dat w]	
# #set sel1 [atomselect top "segname CHI and resid '$r1'" frame $i]


# }
# puts -nonewline $outfile "\n"
for { set i 1 } { $i<=$nframes } { incr i } {
#set time [expr 0.1 * $i]
#puts -nonewline $outfile [format "%.3f" "$time"]
#puts  [format "Time %.3f" "$time"]
# for { set r1 -73 } { $r1<=73 } { incr r1 } {
# set sel1 [atomselect top "protein and noh" frame $i]
# set r2 [expr $r1 * (-1)]
# set sel2 [atomselect top "((segname CHI and resid '$r1') or (segname CHJ and resid '$r2')) and noh" frame $i]
# set contacts [lindex [measure contacts  3.0 $sel1 $sel2] 1]

foreach ind [$h3 get index] {

set sel0 [atomselect top "index $ind" frame 0]
set seli [atomselect top "index $ind" frame $i]

set c0 [lindex [$sel0 get {x y z}] 0]
set c1 [lindex [$seli get {x y z}] 0]

set d [veclength [vecsub $c0 $c1]]
#puts "KUKU"
puts $outfile "$d"

$sel0 delete
$seli delete
}



# set nc [llength $contacts]
# #puts "$contacts"
# #puts [format "Check: %d" "$delta"]
# puts -nonewline $outfile "\t$nc"

# $sel2 delete
# $sel1 delete
# #hbonds -sel1 $sel1 -sel2 $sel2 -dist 3.0 -ang 30 -frames $i:$i -writefile yes -outfile ../analysis_data/dna_hbonds.dat
# #-type unique -writefile yes -outfile ../analysis_data/dna_hbonds.dat

# }
#puts $outfile "\t$ralpha_core_ca\t$rh3\t$rh4\t$rh2a\t$rh2b\t$rh3_1\t$rh3_2\t$rh4_1\t$rh4_2\t$rh2a_1\t$rh2a_2\t$rh2b_1\t$rh2b_2\t$rdnai\t$rdnaj"
}

 close $outfile 

#hbonds -sel1 $chi -sel2 $chj -dist 3.0 -ang 30 -type unique -writefile yes -outfile ../analysis_data/dna_hbonds.dat

exit


