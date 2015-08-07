#!/bin/bash


inputdata="testdata.txt"

######################################
# select which version you want to run:
# R01   - rtype = 1 #
# JTH07 - rtype = 2 #
#####################
 rtype=1 #
##########
if    [ $rtype -eq 1 ]; then
    outfile="tctv_R01_results.out"
elif  [ $ptype -eq 2 ]; then
    outfile="tctv_JTH07_results.out"
fi  
######################################
echo output results to ..................:-  $outfile

./EXAMPLE.exe  $inputdata $outfile



### run program ##
ipython  plot_tctv.py \
$outfile \
