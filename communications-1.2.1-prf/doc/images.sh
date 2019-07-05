#!/bin/sh

images="awgn eyediagram scatterplot"
formats="eps pdf png"
script=commsimages.m

printf "IMAGES ="
for fmt in $formats; do
  for img in $images; do
    printf " $img.$fmt"
  done
done
printf "\n"

for fmt in $formats; do
  fmtupcase=`echo $fmt | awk '{print toupper($0)}'`
  printf "IMAGES_$fmtupcase ="
  for img in $images; do
    printf " $img.$fmt"
  done
  printf "\n"
done

for fmt in $formats; do
  for img in $images; do
    printf "$img.$fmt: commsimages.m\n"
    printf "\t\$(OCTAVE) --path=\$(CURDIR)/../inst -f -q -H --eval \"commsimages ('$img', '$fmt');\"\n"
  done
done
