#!/bin/bash
#PBS -o logfile.log
#PBS -e errorfile.err
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=8
tpdir=`echo $PBS_JOBID | cut -f 1 -d .`
tempdir=$HOME/scratch/job$tpdir
mkdir -p $tempdir
cd $tempdir
cp -R $PBS_O_WORKDIR/* .

export OMP_NUM_THREADS=8
export OMP_PLACES=cores

./hodlrnd_dp_3.app 10 3 > hnd_3_10.txt
./hodlrnd_dp_3.app 15 3 > hnd_3_15.txt
./hodlrnd_dp_3.app 18 3 > hnd_3_18.txt
./hodlrnd_dp_3.app 27 3 > hnd_3_27.txt
./hodlrnd_dp_3.app 30 3 > hnd_3_30.txt
./hodlrnd_dp_3.app 22 3 > hnd_3_22.txt

mv ../job$tpdir $PBS_O_WORKDIR/.

