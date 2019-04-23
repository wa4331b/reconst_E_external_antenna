#!/usr/bin/env bash

#
# BATCH request script
#
#PBS -S /usr/local/bin/tcsh
#PBS -q lx -b a
#PBS -N omi -m e -M scummyman0326@gmail.com -l elapstim_req=02:00:00
#

cd $PBS_O_WORKDIR
pwd
# limit coredumpsize 0
# setenv F_PROGINF DETAIL
# setenv OMP_NUM_THREADS 24
date

setenv DIR $PBS_O_WORKDIR

# now the real deal: the MLFMA computation
# mpirun -ppn 24 -np 360 ./code/MoM/mpi_mlfma --simudir ${SIMU_DIR}
python makeMat.py
python computeEobs.py

date
