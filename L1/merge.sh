n=0
for((q=1;q<2000000;q+=100000))
   do
   q1="`expr $q + 49999`"
   q2="`expr $q1 + 1`"
   q3="`expr $q2 + 49999`"
   ~/code/L1/L1_merg ~/L1_data/L1_${q}_${q1}.dat ~/L1_data/L1_${q2}_${q3}.dat temp_${n}.dat
   n="`expr $n + 1`"
   done
l=0
for((m=0;m<20;m+=2))
   do
   m1="`expr $m + 1`"
   ~/code/L1/L1_merg temp_${m}.dat temp_${m1}.dat temp1_${l}.dat
   l="`expr $l + 1`"
   done
l=0
for((m=0;m<10;m+=2))
   do
   m1="`expr $m + 1`"
   ~/code/L1/L1_merg temp1_${m}.dat temp1_${m1}.dat temp2_${l}.dat
   l="`expr $l + 1`"
   done
l=0
for((m=0;m<4;m+=2))
   do
   m1="`expr $m + 1`"
   ~/code/L1/L1_merg temp2_${m}.dat temp2_${m1}.dat temp3_${l}.dat
   l="`expr $l + 1`"
   done
l=0
for((m=0;m<2;m+=2))
   do
   m1="`expr $m + 1`"
   ~/code/L1/L1_merg temp3_${m}.dat temp3_${m1}.dat temp4.dat
   l="`expr $l + 1`"
   done
~/code/L1/L1_merg temp4.dat temp2_4.dat temp.dat


