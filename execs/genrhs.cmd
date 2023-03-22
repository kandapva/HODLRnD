#!/bin/bash
#PBS -o logfile.log
#PBS -e errorfile.err
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8
tpdir=`echo $PBS_JOBID | cut -f 1 -d .`
tempdir=$HOME/scratch/job$tpdir
mkdir -p $tempdir
cd $tempdir
cp -R $PBS_O_WORKDIR/generate_rhs_1overR2.app .

export OMP_NUM_THREADS=8
export OMP_PLACES=cores

echo "OneOverR2" >> /lfs/usrhome/phd/ma16d300/HODLRnD/output/naive_mat_vec.txt
./generate_rhs_1overR2.app 10 >> /lfs/usrhome/phd/ma16d300/HODLRnD/output/naive_mat_vec.txt
./generate_rhs_1overR2.app 15 >> /lfs/usrhome/phd/ma16d300/HODLRnD/output/naive_mat_vec.txt
./generate_rhs_1overR2.app 18 >> /lfs/usrhome/phd/ma16d300/HODLRnD/output/naive_mat_vec.txt
./generate_rhs_1overR2.app 22 >> /lfs/usrhome/phd/ma16d300/HODLRnD/output/naive_mat_vec.txt
./generate_rhs_1overR2.app 27 >> /lfs/usrhome/phd/ma16d300/HODLRnD/output/naive_mat_vec.txt
./generate_rhs_1overR2.app 30 >> /lfs/usrhome/phd/ma16d300/HODLRnD/output/naive_mat_vec.txt

#mv ../job$tpdir $PBS_O_WORKDIR/.

