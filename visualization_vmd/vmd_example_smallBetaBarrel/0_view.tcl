#####################
#Default script to do some visualization in VMD
#and a rotation movie
#Launch as vmd -e view.tcl
#####################

####Setup useful funcions #####

##Rendering option
## Tachyon internal
## AA on, AA 1.00, AOdirect 0.6

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

mol new chainH.pdb
#Remove defauls representation
mol delrep 0 top
#or redefine its style
# mol modstyle 0 0 NewCartoon 0.840000 20.000000 2.630000 0

#########
#####Representations
########
set col [list "blue3" "green" "yellow2" "red3" "gray" "orange" "cyan" "purple"]

set sel [atomselect 0 "protein or nucleic"] 
set chains [lsort -ascii -unique [$sel get chain]] 


set i 0
foreach c $chains  {

color Chain $c [lindex $col $i]
incr i
if {$i>7} {set i 0}
}
###Add representations
###								thickness resolution aspect_ratio
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Structure
mol selection {all}
mol material AOEdgy
mol addrep top
mol selupdate 0 top 0
mol colupdate 0 top 0



#Apply scene scaling and translations
# scale by $scale
# translate by $transx $transy $transz
# display update ui



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
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Structure 1kq2" size 1.5 thickness 3
incr txtlncount 2

add_text_layer CHAINS



set matnum 22

set i 0
foreach c $chains  {
draw color [lindex $col $i]
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Chain $c" size 1.5 thickness 3
incr txtlncount 1
incr i
if {$i>7} {set i 0}
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
# source vps.tcl
# retrieve_vp 1


#####
###Rendering and rotation, trajectory play
####

Check for dat dir
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

exec rm -rf output/1kq2_A.mov
#This is how one can make a move out of set of images
exec ffmpeg -i dat/%d.dat.tga -s 1000x750 -q:v 0 -pix_fmt yuv420p  output/1kq2_A.mov
exec rm -r dat
# ffmpeg -i dat/%d.dat.tga -s 1000x750 -q:v 0 -pix_fmt yuv420p  movie.wmv
exit
