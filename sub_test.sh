#!/bin/sh


for((j=0;j<4;j++))

do

bsub -ql   ./data.sh 21 $j

done

for((j=0;j<7;j++))

do

bsub -ql   ./data.sh 23 $j

done
