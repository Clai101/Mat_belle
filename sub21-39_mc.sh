#!/bin/sh


for((j=0;j<4;j++))

do

bsub -ql   ./mc.sh 21 $j

done

for((j=0;j<7;j++))

do

bsub -ql   ./mc.sh 23 $j

done

for((j=0;j<22;j++))

do

bsub -ql   ./mc.sh 25 $j

done

for((j=0;j<17;j++))

do

bsub -ql   ./mc.sh 27 $j

done

for((j=0;j<18;j++))

do

bsub -ql   ./mc.sh 31 $j

done

for((j=0;j<9;j++))

do

bsub -ql   ./mc.sh 33 $j

done

for((j=0;j<7;j++))

do

bsub -ql   ./mc.sh 35 $j

done

for((j=0;j<20;j++))

do

bsub -ql   ./mc.sh 37 $j

done

for((j=0;j<14;j++))

do

bsub -ql   ./mc.sh 39 $j

done
