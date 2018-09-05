#!/bin/bash

format='+%Y/%m/%d-%H:%M:%S'

date $format

job_num=$(($SLURM_ARRAY_TASK_ID))
output_file=$output_dir/$job_num.root
filelist=$job_num.list
sed -ie 'H;${x;s/\n/,/g;s/^,//;p;};d' $filelist

echo "loading " $hadesroot
source hadesroot
echo "executing $executable  $job_num.list $output_file -1"
$executable  $job_num.list $output_file -1
rm $filelist

echo JOB FINISHED!
date $format