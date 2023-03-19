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

./hodlrnd_dp_1.app 10 > hnd_1_10.txt
./hodlrnd_dp_1.app 15 > hnd_1_15.txt
./hodlrnd_dp_1.app 18 > hnd_1_18.txt
./hodlrnd_dp_1.app 27 > hnd_1_27.txt
./hodlrnd_dp_1.app 32 > hnd_1_32.txt
./hodlrnd_dp_1.app 22 > hnd_1_22.txt

mv ../job$tpdir $PBS_O_WORKDIR/.

