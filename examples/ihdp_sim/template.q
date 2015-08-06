#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -l mem=1GB
#PBS -N _JOBNAME_
#PBS -M vjd4@nyu.edu
#PBS -m e
#PBS -j oe
#PBS -o jobs/_JOBNAME_.o

module load r/intel/3.2.0
RUNDIR=$HOME/gpcomp

export METHOD=_METHOD_
export OVERLAP=_OVERLAP_
export COVARIATES=_COVARIATES_
export VERBOSE=_VERBOSE_
export START=_START_
export END=_END_
cd $RUNDIR
R --no-save < runJob.R
