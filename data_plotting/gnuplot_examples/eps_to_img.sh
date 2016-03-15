#!/bin/bash
#Note you can modify and save eps file with Adobe Illustartor perfectly!!!!!!!
#Possible gs devices  http://www.gnu.org/software/ghostscript/devices.html
#script converts eps to image files using GhostView
#scale resolution and image size simulaneously
#gnuplot plot_to_eps.scr;
#gs -dNOPAUSE -sDEVICE=tiff24nc -sOutputFile=test.tiff -q  -dBATCH -g130x90 -r20x20 test.eps;

#GS offers no color lzw compression, so png is the best variant
#otherwise output to ppmraw and use libtiff or Netpbm utils (Netpbm is better)
gs -dNOPAUSE -sDEVICE=png16m -sOutputFile=test.png -q  -dBATCH -g650x450 -r100x100 test.eps;
#gs -dNOPAUSE -sDEVICE=ppmraw -sOutputFile=test.ppm -q  -dBATCH -g650x450 -r100x100 test.eps;
#ppm2tiff test.ppm test.tif
#pnmtotiff -lzw msd_x_long.ppm > msd_x_long.tif



#mysend test.tiff;
