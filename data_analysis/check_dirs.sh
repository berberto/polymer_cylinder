#!/bins/bash

for dir in $1/avjmp_*; 
do 
    echo $dir, `ls -l $dir | awk '$5==0' | wc -l` ; 
done
