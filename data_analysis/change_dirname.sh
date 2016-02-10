#!/bin/csh

for name in `ls -l /scratch/madorisi/results/polymers/ | awk '$7==20 || $7==19' | awk '{print $9}' | awk -F'/' '{print $1}' | awk -F'_' '{print $2}'`; 
do
    new_name=$(printf 'avjmp_%.3e\n' $name)
    mv 'avjmp_'$name $new_name
done
