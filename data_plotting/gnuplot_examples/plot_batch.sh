#!/bin/bash
#script converts eps to tiff using GhostView
#scale resolution and image size simulaneously
#gnuplot plot_to_eps.scr;
#gs -dNOPAUSE -sDEVICE=tiff24nc -sOutputFile=test.tiff -q  -dBATCH -g130x90 -r20x20 test.eps;
LIST=$( ls . | grep .scr );
for n in $LIST;
do
echo "Processing $n ...";
gnuplot $n;
done
echo "Done!!!"
