#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -p XXX
#SBATCH -A PROJECTNO
#SBATCH -o out.txt
#SBATCH -e err.txt
#SBATCH -t 38:15:00
#SBATCH -J JOBNAME
#SBATCH --mail-user=EMAILUSER
#SBATCH --mail-type=END

export NCORES=$(($SLURM_JOB_NUM_NODES*16))

cp *.py $TMPDIR
cd $TMPDIR


while sleep 1800; do
    # This will be executed once per every 3600 seconds
    cp -r run_*/ $SLURM_SUBMIT_DIR
done &     # The &-sign after the done-keyword places 
           # the while-loop in a sub-shell in the background
LOOPPID=$! # Save the PID of the subshell running the loop


mpirun python kmc_master_parallel.py 300. 2E3 1E3 1E-1

# Copy the files back after run:

kill $LOOPPID

cp -r run_*/ $SLURM_SUBMIT_DIR

