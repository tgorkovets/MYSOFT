#####################
#Default script to do some visualization in VMD
#and a rotation movie
#Launch as vmd -e view.tcl
#####################

####Setup useful funcions #####

#Easy text annotations
source add_text_layer.tcl 
source view_change_render.tcl

proc pause {{message "Hit Enter to continue ==> "}} {
    puts -nonewline $message
    flush stdout
    gets stdin
}

#Use this to stop-debug script
# pause;
# exit

###############
####Define view movements variables and view size
###############
set dispw 1000
set disph 750
set scale 1.75
set transx 0.3
set transy -0.40
set transz 0

display resize $dispw $disph

##########################
#Set general display parameters
#########################
display rendermode GLSL
axes location off

display antialias on

display eyesep       0.960000
display focallength  6.000000
display height       6.000000
display distance     -2.000000
display projection   Orthographic
display nearclip set 0.500000
display farclip  set 10.000000
display depthcue off

color Display Background white
color scale method RWB


# Some tachyon rendering tweaks
display ambientoclussion on
display shadows on
display antialias on
display aoambient 0.75
display aodirect 0.3

##############
#####Color adjustment
##############

color change rgb 24 0.15 0.25 0.93
color change rgb 30 0.81 0.18 0.18
color Highlight Nonback 6
color Highlight Nucback 2
color change rgb  0 0.1 0.2 0.7 ;# blue
color change rgb  1 0.7 0.2 0.1 ;# red
color change rgb  3 0.7 0.4 0.0 ;# orange
color change rgb  4 0.8 0.7 0.1 ;# yellow
color change rgb  7 0.1 0.7 0.2 ;# green
color change rgb 10 0.1 0.7 0.8 ;# cyan
color change rgb 11 0.6 0.1 0.6 ;# purple
color change rgb 23 0.550000011920929 0.7099999785423279 0.9800000190734863
color change rgb 24 0.15 0.25 0.93
color change rgb 30 0.81 0.18 0.18
####

###########
####Structure loading
###########
# mol load psf ../analysis_data/only_nucl_init.psf
# mol addfile ../analysis_data/only_nucl_init.pdb
# mol load psf ../analysis_data/only_nucl_init.psf
# mol addfile ../analysis_data/only_nucl_init.pdb
# mol ssrecalc top
# mol addfile ../analysis_data/md_nucl.dcd step 10 waitfor all


mol load pdb her2_active.pdb
#Remove defauls representation
mol delrep 0 top
#or redefine its style
# mol modstyle 0 0 NewCartoon 0.840000 20.000000 2.630000 0

mol load pdb her2_active_L869R.pdb
#Remove defauls representation
mol delrep 0 top
#or redefine its style
# mol modstyle 0 0 NewCartoon 0.840000 20.000000 2.630000 0





#####
####RMSD manipulation of structures
####

set sel0 [atomselect  0 "name CA" frame 0]	 
set sel0_all [atomselect  0 "all" frame 0]	 

set sel1 [atomselect 1 "name CA" frame 0]	 
set sel1_all [atomselect  1 "all" frame 0]	

$sel1_all move [measure fit $sel1 $sel0]


#########
#####Representations
########

###Add coloring by parameters
color Segname CHA blue3
color Segname CHB green
color Segname CHC yellow2
color Segname CHD red3
color Segname CHE blue3
color Segname CHF green
color Segname CHG yellow2
color Segname CHH red3

###Add representations
###								thickness resolution aspect_ratio
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color ColorID 24
mol selection {all}
mol material AOEdgy
mol addrep top
mol selupdate 0 top 0
mol colupdate 0 top 0

mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color ColorID 7
mol selection {(resid > 863) and (resid < 885)}
mol material AOEdgy
mol addrep top
mol selupdate 0 top 0
mol colupdate 0 top 0


mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color ColorID 9
mol selection {(resid > 761) and (resid < 776)}
mol material AOEdgy
mol addrep top
mol selupdate 0 top 0
mol colupdate 0 top 0

#					scale bond_radius sphere_resolution bond_resolution
mol representation CPK 2.0 2.5 12 10
mol color ColorID 1
mol selection {(resid 869)}
mol material AOEdgy
mol addrep top
mol selupdate 0 top 0
mol colupdate 0 top 0

#					scale bond_radius sphere_resolution bond_resolution
mol representation CPK 2.0 2.5 12 10
mol color ColorID 2
mol selection {(resid 869)}
mol material AOEdgy
mol addrep 0
mol selupdate 0 top 0
mol colupdate 0 top 0

mol representation CPK 2.0 2.5 12 10
mol color ColorID 4
mol selection {(same resid as ((within 4 of resid 869) and sidechain)) and not resid 869 868 870}
mol material AOEdgy
mol addrep top
mol selupdate 0 top 0
mol colupdate 0 top 0

#########Examples###############
# #  H3 H3 Dimer representation
# mol representation NewCartoon 0.840000 20.000000 2.630000 0
# mol color Segname
# mol selection {segname CHA CHE}
# mol material AOShiny
# mol addrep top
# mol selupdate 0 top 0
# mol colupdate 0 top 0

# # H4 H4 Dimer representation
# mol representation NewCartoon 0.840000 20.000000 2.630000 0
# mol color Segname
# mol selection {segname CHB CHF}
# mol material AOShiny
# mol addrep top
# mol selupdate 1 top 0
# mol colupdate 1 top 0

# # H2a H2a Dimer representation
# mol representation NewCartoon 0.840000 20.000000 2.630000 0
# mol color Segname
# mol selection {segname CHC CHG}
# mol material AOShiny
# mol addrep top
# mol selupdate 2 top 0
# mol colupdate 2 top 0

# # H2b H2b Dimer representation
# mol representation NewCartoon 0.840000 20.000000 2.630000 0
# mol color Segname
# mol selection {segname CHD CHH}
# mol material AOShiny
# mol addrep top
# mol selupdate 3 top 0
# mol colupdate 3 top 0

# # Dna sugar-phosphate backbone representation
# mol representation NewCartoon 1.090000 10.000000 2.070000 1
# mol color ColorID 6
# mol selection {nucleic and backbone}
# mol material AOEdgy
# mol addrep top
# mol selupdate 4 top 0
# mol colupdate 4 top 0

# # Dna nucleobases representation
# mol representation Licorice 0.300000 10.000000 10.000000
# mol color ColorID 6
# mol selection { nucleic and noh and not name P O1P O2P O3' O5' C5' and resname GUA CYT}
# mol material AOEdgy
# mol addrep top
# mol selupdate 5 top 0
# mol colupdate 5 top 0

# mol representation Licorice 0.300000 10.000000 10.000000
# mol color ColorID 3
# mol selection {nucleic and noh and not name P O1P O2P O3' O5' C5' and resname THY ADE}
# mol material AOEdgy
# mol addrep top
# mol selupdate 6 top 0
# mol colupdate 6 top 0

# #Key argininges
# mol representation VDW
# mol color Orange
# mol selection {(segname CHA CHE and resid 83 63 49) or (segname CHB CHF and resid 45) or (segname CHC CHG and resid 42 77) or (segname CHD CHH and resid 30)}
# mol material AOShiny
# mol addrep top
# mol selupdate 7 top 0
# mol colupdate 7 top 0
#################END Examples###############

#Apply scene scaling and translations
scale by $scale
translate by $transx $transy $transz
display update ui



#########
###Text annotations
#########

#set number of main molecules, so that add_text_layer knows where real molecules end
set molnum [expr [molinfo num]]

#Set parameters for text position and spacing
set txtx [expr -($dispw/2) + 20 ]
set txty [expr ($disph/2) - 20 ]
# set step of txt lines
set txtstep  [expr $disph/30 ]
set txtlncount 0

# Add heading
add_text_layer HEADNAME
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "HER2 active conformation (3pp0)" size 1.5 thickness 3
incr txtlncount 2

add_text_layer HISTONES

set matnum 22

#Set array of annotation labels and colors
set text(0) "Activation loop"
set text(1) "Alpha-helix"
set text(2) "L869"
set text(3) "L869R"
set text(4) "Interacting residues"
set color_text(0) 7
set color_text(1) 9
set color_text(2) 2
set color_text(3) 1
set color_text(4) 4


#Display them
for {set r 0} {$r < 5} {incr r 1} {
draw color $color_text($r)
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "$text($r)" size 1.5 thickness 3
incr txtlncount 1
}

display update

#return top to displayed molecule
mol top 0

########
####Setting manual view
#######
#This block makes use of view_change_render.tcl.
#Open VMD orient molecules as you like
#save_vp 1
#write_vps vps.tcl
#Uncomment the next two line
source vps.tcl
retrieve_vp 1


#####
###Rendering and rotation, trajectory play
####

#Check for dat dir
exec mkdir -p dat

for {set i 0} {$i <= 360} {incr i 1} {

rotate y by -1

##########Example for playing MD trajectory
# animate goto $i
# mol delete top
# add_text_layer TIME
# draw color 0
# set time [format "Time: %5.1f ns" [expr $i * 5.0]]
# draw text " $txtx [expr $txty-(27*$txtstep)] 0 " $time size 1.5 thickness 3
################

display update
render snapshot dat/$i.dat.tga
# render TachyonInternal dat/$i.dat.tga

}

#This is how one can make a move out of set of images
# ffmpeg -i dat/%d.dat.tga -s 1000x750 -q:v 0 -pix_fmt yuv420p  movie.mov
# ffmpeg -i dat/%d.dat.tga -s 1000x750 -q:v 0 -pix_fmt yuv420p  movie.wmv

