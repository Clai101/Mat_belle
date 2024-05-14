#!/bin/sh



for((j=0;j<29;j++))

do 

 bsub -ql   ./data.sh 7 $j

done

for((j=0;j<13;j++))

do

 bsub -ql   ./data.sh 9 $j

done

for((j=0;j<14;j++))

do

 bsub -ql   ./data.sh 11 $j

done

for((j=0;j<17;j++))

do

 bsub -ql   ./data.sh 13 $j

done

for((j=0;j<15;j++))

do

 bsub -ql  ./data.sh 15 $j

done

for((j=0;j<11;j++))

do

 bsub -ql   ./data.sh 17 $j

done

for((j=0;j<18;j++))

do

 bsub -ql   ./data.sh 19 $j

done

for((j=0;j<4;j++))

do

bsub -ql   ./data.sh 21 $j

done

for((j=0;j<7;j++))

do

bsub -ql   ./data.sh 23 $j

done

for((j=0;j<22;j++))

do

bsub -ql   ./data.sh 25 $j

done

for((j=0;j<17;j++))

do

bsub -ql   ./data.sh 27 $j

done

for((j=0;j<18;j++))

do

bsub -ql   ./data.sh 31 $j

done

for((j=0;j<9;j++))

do

bsub -ql   ./data.sh 33 $j

done

for((j=0;j<7;j++))

do

bsub -ql   ./data.sh 35 $j

done

for((j=0;j<20;j++))

do

bsub -ql   ./data.sh 37 $j

done

for((j=0;j<14;j++))

do

bsub -ql   ./data.sh 39 $j

done


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











