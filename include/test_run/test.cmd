#!/bin/bash
#PBS -o logfile.log
#PBS -e errorfile.err
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1
tpdir=`echo $PBS_JOBID | cut -f 1 -d .`
tempdir=$HOME/scratch/job$tpdir
mkdir -p $tempdir
cd $tempdir
cp -R $PBS_O_WORKDIR/* .

./hd4dim 8  > hdd_0_8.txt
./hd4dim 10 > hdd_0_10.txt
./hd4dim 12 > hdd_0_12.txt
./hd4dim 14 > hdd_0_14.txt

mv ../job$tpdir $PBS_O_WORKDIR/.
