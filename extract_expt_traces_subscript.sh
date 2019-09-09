#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --mem=80GB
#SBATCH --job-name=extractTraces
#SBATCH --mail-type=END
#SBATCH --array=1-50
#SBATCH -o /gpfs/home/stetlb01/logs/extract_%A_%a.log
$SBATCH -e /gpfs/home/stetlb01/logs/extract_%A_%a.log

module purge
module load matlab/R2018a

tabfile=$1
basedir=$2

nlines=$(wc -l < "$tabfile")
nlines=$(echo $nlines | awk '{print $1-1}')


if [ $SLURM_ARRAY_TASK_ID -gt nlines ]
then
    echo "Array task ID greater than needed. Exiting."
    exit
fi





{
  echo $tif
  echo running command "addpath('/gpfs/home/stetlb01/scripts');addpath('/gpfs/home/stetlb01/Holography_Analysis');addpath('/gpfs/home/stetlb01/JG_Functions');run_extract_from_expt_table('$tabfile','$basedir',$SLURM_ARRAY_TASK_ID);exit"
  matlab -nodisplay -r "addpath('/gpfs/home/stetlb01/scripts');addpath('/gpfs/home/stetlb01/Holography_Analysis');addpath('/gpfs/home/stetlb01/JG_Functions');run_extract_from_expt_table('$tabfile','$basedir',$SLURM_ARRAY_TASK_ID);exit"
     }  2>&1

exit
