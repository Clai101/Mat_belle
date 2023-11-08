#!/bin/sh



for((j=0;j<29;j++))

do 

 bsub -ql   ./mc.sh 7 $j

done

for((j=0;j<13;j++))

do

 bsub -ql   ./mc.sh 9 $j

done

for((j=0;j<14;j++))

do

 bsub -ql   ./mc.sh 11 $j

done

for((j=0;j<17;j++))

do

 bsub -ql   ./mc.sh 13 $j

done

for((j=0;j<15;j++))

do

 bsub -ql  ./mc.sh 15 $j

done

for((j=0;j<11;j++))

do

 bsub -ql   ./mc.sh 17 $j

done

for((j=0;j<18;j++))

do

 bsub -ql   ./mc.sh 19 $j

done













