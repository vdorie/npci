#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -l mem=1GB
#PBS -N jobname
#PBS -M vjd4@nyu.edu
#PBS -m e
#PBS -j oe
#PBS -o jobs/jobname.o

module load r/intel/3.2.0
RUNDIR=$HOME/gpcomp

export METHOD=method
export SETTING=setting
export START=start
export END=end
cd $RUNDIR
R --no-save < runJob.R
