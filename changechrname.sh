#!/bin/sh
for file in mapped/*.bam
do
filename=`echo $file | awk -F"/" '{ print $NF}' | cut -d "." -f1 `
samtools view -H mapped/$file | sed -e 's/SN:([0-9XY])/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - mapped/$file > mapped/${filename}_chr.bam
done
