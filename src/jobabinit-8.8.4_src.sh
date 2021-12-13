#!/bin/bash
#SBATCH -J  jobabinit-8.8.4_src
#SBATCH --ntasks 
#SBATCH --partition=ccq
#SBATCH --qos=ccq
#SBATCH --constraint=skylake
#SBATCH -o out.%J
#SBATCH -e err.%J
#SBATCH -t 00:00:00
#SBATCH --mail-type=ALL

module load intel/compiler/2017-4
module load intel/mpi/2017-4
module load intel/mkl/2017-4

mpirun -np   /mnt/home/cdreyer/apps/vlfrc_abinit-8.8.4/abinit-8.8.4/src/98_main/abinit <  >& abinit.out
rm *_xo_*DEN*
rm *_xo_*WF*
rm *_xo_*POT*
rm *_x_*
