#!/bin/bash
for st in 14 5000 26000 236000
do
   echo "/projects/Zeros_of_Lfunctions/zzeros/zeros1/zeros_${st}.dat"
done
for((st=446000;st<11155646000;st+=2100000))
do 
   echo "/projects/Zeros_of_Lfunctions/zzeros/zeros1/zeros_${st}.dat"
done
for((st=11155646000;st<17203646000;st+=2100000))
do 
   echo "/projects/Zeros_of_Lfunctions/zzeros/zeros2A/zeros_${st}.dat"
done
for((st=17203646000;st<20950046000;st+=2100000))
do 
   echo "/projects/Zeros_of_Lfunctions/zzeros/zeros4A/zeros_${st}.dat"
done
for((st=20950046000;st<29988446000;st+=2100000))
do 
   echo "/projects/Zeros_of_Lfunctions/zzeros/zeros6A/zeros_${st}.dat"
done
for((st=29988446000;st<=30607946000;st+=2100000))
do 
   echo "/projects/Zeros_of_Lfunctions/zzeros/zeros8A/zeros_${st}.dat"
done
