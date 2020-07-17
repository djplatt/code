#!/bin/bash
zd="${HOME}/data/first_billion"
for((om=$1;om<$2;om++))
do
   ${HOME}/code/chris/M_sum 200 ${zd}/zeros_14.dat ${om} 1000 1.64e19 > ${HOME}/data/chris/temp_${3}.out
   echo "${om} $(grep "test on" ${HOME}/data/chris/temp_${3}.out | awk '// {print $16}')" >> ${HOME}/data/chris/core_${3}.out
done

