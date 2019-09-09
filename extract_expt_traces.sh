#!/bin/bash

tabfile=$1
basedir=$2

nlines=$(wc -l < "$tabfile")
nlines=$(echo $nlines | awk '{print $1-1}')

sbatch --array=1-$nlines $(pwd)/extract_expt_traces_subscript.sh $tabfile $basedir

