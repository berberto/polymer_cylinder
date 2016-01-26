#!/bin/bash
#TO BE specified only in case of usage of a given queue (now commented)
##PBS -q reserved2

#define the job array
##PBS -t 0-19

#defne the number of nodes (and cores)
#PBS -l nodes=1:ppn=20

#define the walltime
##PBS -l walltime=12:00:00
#PBS -l walltime=00:01:00

#define the error output file as JOBNAME.e[JOBID] || JOBNAME.o[JOBID]
#PBS -j oe

#send me an e-mail when job ends
#PBS -m e

#Define the JOB name
#PBS -N cylinder

#IT HAS TO BE HERE
#PBS -T flush_cache
#PBS
#

# set up environment to import pre-compiled executables from AK.
#module load intel

#set the current directory as the working directory (where the job is submitted from)
#by default the job will start on the $HOME directory on the remote node
PBS_O_WORK=/scratch/apezzotta/control/main

#echo $PBS_O_WORK

# Read the file in parameter and fill the array named "array"
getArray() {
    array=() # Create array
    while IFS= read -r line # Read a line
    do
        array+=("$line") # Append line to the array
    done < "$1"
}

getArray "/scratch/apezzotta/control/main/input"

#for e in "${array[@]}"
#do
#    echo "$e"
#done


cd $PBS_O_WORK

module load openmpi-x86_64

module load gsl/1.16/gnu/4.9.2

#cd ../devel/random && make && cd ../../main && make

echo "${array[$PBS_ARRAYID]}"

python2.7 run.py ${array[$PBS_ARRAYID]}

#mpirun -n 20 cylinder ${array[$PBS_ARRAYID]} 10 
