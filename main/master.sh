#!/bin/bash


make

list=$(cat output/genf_z/pars.dat)

for i in $list; do
	./drift $i
#	echo $i
done 
