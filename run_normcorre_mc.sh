#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=5:00:00
#SBATCH --mem=60GB
#SBATCH --job-name=motionCorrection
#SBATCH --mail-type=END
#SBATCH --array=0-100
#SBATCH -o /gpfs/home/stetlb01/logs/mc_%A_%a.log
$SBATCH -e /gpfs/home/stetlb01/logs/mc_%A_%a.log

module purge
module load matlab/R2018a

## vid_dir (1st arg): video and template file directory
## tifs (2nd arg, mapped into list form): text file containing video file names (just the file names, no parent directories!)
## template (3rd arg): alignment template (file name only)
## redchannel (4th arg): 0 for no red channel in the input, 1 for red channel in the input (dropped in the output in all cases)
## replace (5th arg): 1 to replace existing output file, 0 otherwise

vid_dir=$1
mapfile -t tifs < $2
template=$3
redchannel=$4
replace=$5

cd $vid_dir

if [ $SLURM_ARRAY_TASK_ID -ge  ${#tifs[@]}  ]
then
    echo "Array task ID greater than needed. Exiting."
    exit
fi


tif=${tifs[$SLURM_ARRAY_TASK_ID]}

{
  echo $tif
  matlab -nodisplay -r "addpath('/gpfs/home/stetlb01/scripts');normcorremotioncorrection('$tif','$template',$redchannel,$replace);exit"
     }  2>&1

exit
