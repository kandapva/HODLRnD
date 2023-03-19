#!/bin/bash
#PBS -o logfile.log
#PBS -e errorfile.err
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=8
tpdir=`echo $PBS_JOBID | cut -f 1 -d .`
tempdir=$HOME/scratch/job$tpdir
mkdir -p $tempdir
cd $tempdir
cp -R $PBS_O_WORKDIR/* .

export OMP_NUM_THREADS=8
export OMP_PLACES=cores
./generate_rhs.app 30 > rhs_30.txt

mv ../job$tpdir $PBS_O_WORKDIR/.

