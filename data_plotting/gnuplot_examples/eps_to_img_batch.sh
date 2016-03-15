#!/bin/bash
#script converts eps to tiff using GhostView
#scale resolution and image size simulaneously
#gnuplot plot_to_eps.scr;
#gs -dNOPAUSE -sDEVICE=tiff24nc -sOutputFile=test.tiff -q  -dBATCH -g130x90 -r20x20 test.eps;
LIST=$( ls . | grep .eps | sed -e 's/.eps//g');
for n in $LIST;
do
echo "Processing $n ...";
#gs -dNOPAUSE -sDEVICE=tiff24nc -sOutputFile=../tiff/$n.tiff -q  -dBATCH -g650x450 -r100x100 $n.eps;
#gs -dNOPAUSE -sDEVICE=png16m -sOutputFile=../tiff/$n.png -q  -dBATCH -g1300x900 -r200x200 $n.eps;
gs -dNOPAUSE -sDEVICE=ppmraw -sOutputFile=temp.ppm -q  -dBATCH  -g1300x900 -r200x200 $n.eps;
pnmcrop -margin=5 temp.ppm > temp2.ppm;
pnmtotiff -lzw temp2.ppm > ../img/$n.tif;
rm temp.ppm temp2.ppm;
#gs -dNOPAUSE -sDEVICE=x11  -q  -dBATCH -g1300x900 -r200x200 $n.eps;
done
echo "Done!!!"
