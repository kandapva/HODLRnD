#!/bin/bash
#PBS -o logfile.log
#PBS -e errorfile.err
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=1
tpdir=`echo $PBS_JOBID | cut -f 1 -d .`
tempdir=$HOME/scratch/job$tpdir
mkdir -p $tempdir
cd $tempdir
cp -R $PBS_O_WORKDIR/* .

export OMP_NUM_THREADS=8
export OMP_PLACES=cores

./hodlrnd_dp_0.app 10 > hnd_0_10.txt
./hodlrnd_dp_0.app 15 > hnd_0_15.txt
./hodlrnd_dp_0.app 18 > hnd_0_18.txt
./hodlrnd_dp_0.app 27 > hnd_0_27.txt
./hodlrnd_dp_0.app 32 > hnd_0_32.txt
./hodlrnd_dp_0.app 22 > hnd_0_22.txt

mv ../job$tpdir $PBS_O_WORKDIR/.

