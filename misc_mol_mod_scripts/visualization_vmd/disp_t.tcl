# vmd tcl procedure: display simulation time with a thermometer like bar.
# 
# $Id: display_time_bar.tcl,v 1.1 2003/07/16 13:29:12 wwwadmin Exp $
# Time-stamp: <akohlmey 13.04.2003 15:16:19 timburton.bochum>
#
# Copyright (c) 2003 by <Axel.Kohlmeyer@theochem.ruhr-uni-bochum.de>
#
# TODO: find a smart way to position the bar at the side of the
#       screen regardless of the number and size of molecules loaded.
#

# arguments:
#  mult:   length of a timestep  (default 1.0)
#  tdelta: increment of the tics (default 50)
#  unit:   unit of a timestep    (default ps)
#  tmol:   molecule id where the dial is added (default: current top)

proc do_time {args} { disp_t }
trace variable vmd_frame(0) w do_time

proc disp_t {{tmol top}} {
    global disp_t_molid

    display update off
    if {![string compare $tmol top]} {
        set tmol [molinfo top]
    }

    set oldtop [molinfo top]
    if {![info exists disp_t_molid]} {
        set disp_t_molid [mol new]
        mol rename $disp_t_molid "time bar $tmol"
        mol fix $disp_t_molid
    }
    set bar $disp_t_molid
    mol top $oldtop

    # first delete the old bar (if any)
    graphics $bar delete all


    #########################################
    # step size
    set tstep [molinfo $tmol get frame]
    if {$tstep < 1} {set tstep 1}
    set tmax  [molinfo $tmol get numframes]
    if {$tmax < 1} {set tmax 1}

    # draw caption
    graphics $bar color white


    set ttime [expr $tstep * 1 ]
    set pt [format "%3d" $ttime]
    graphics $bar text {-0.7 -0.7 0} "Time $pt ps" size 1 thickness 1


    if { $tstep < 500 } { set temp 0
    } else  { set temp [expr ($tstep - 500) * 0.3] }
    if { $tstep > 1500 }  { set temp 300 }

    set ptemp [format "%3.0f" $temp]
    graphics $bar text {-0.7 -0.55 0} "Temp $ptemp K" size 1 thickness 1

    display update on
    # restore top molecule (just in case we messed it up).
    mol top $oldtop
}

############################################################
# Local Variables:
# mode: tcl
# time-stamp-format: "%u %02d.%02m.%y %02H:%02M:%02S %s"
# End:
############################################################
