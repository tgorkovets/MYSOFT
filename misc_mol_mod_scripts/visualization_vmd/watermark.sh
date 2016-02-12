#!/bin/bash
#script adds a watermark to your .mov video using ffmpeg libvau filters
echo -n "Задайте расширение файла видео без точки :"
read type
echo -n "Файл содержащий водяной знак: "
read wmfile
LIST=$(ls . | grep .$type | sed -e "s/.$type//g");
for n in $LIST;
do
echo "Processing $n ...";
ffmpeg -i $n.mov -sameq -vf "movie=$wmfile [watermark]; [in][watermark] overlay=main_w-overlay_w-0:0 [out]" $n\_wm.mov;
done
echo "Done!"
