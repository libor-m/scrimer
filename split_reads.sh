#!/bin/bash

# split fastq files for all .bcfiles found in current directory
# using fastx tools' barcode splitter, 
# input names are based on names of the .bcfile
# output to split/
#

# the splitter
SPL=~/data/sw_wog/fastx/bin/fastx_barcode_splitter.pl

for bcfile in *.bcfile
do
  sample=$(basename "$bcfile" .bcfile)
  echo "Procesing $sample"
  cat $sample.fastq | $SPL --bcfile $bcfile --prefix split/$sample_ --suffix ".fastq" --bol 
done
