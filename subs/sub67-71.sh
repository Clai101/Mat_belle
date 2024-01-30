#!/bin/sh


for((j=0;j<10;j++))

do

bsub -ql   ./data.sh 43 $j

done

for((j=0;j<10;j++))

do

bsub -ql   ./data.sh 53 $j

done

for((j=0;j<7;j++))

do

bsub -ql   ./data.sh 67 $j

done

for((j=0;j<14;j++))

do

bsub -ql   ./data.sh 69 $j

done

for((j=0;j<22;j++))

do

bsub -ql   ./data.sh 71 $j

done

