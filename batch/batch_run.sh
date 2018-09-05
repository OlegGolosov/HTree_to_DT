#!/bin/bash

format='+%Y/%m/%d-%H:%M:%S'

date $format

job_num=$(($SLURM_ARRAY_TASK_ID))
output_file=$output_dir/$job_num.root
filelist=$lists_dir/$job_num.list
#sed -ie 'H;${x;s/\n/,/g;s/^,//;p;};d' $filelist

while read line; do    
    input_files=$input_files","$line    
done < $filelist

echo "loading " $hadesroot
source $hadesroot
echo "executing $executable  $input_files $output_file $n_events"
$executable  $input_files $output_file $n_events 
#rm $filelist

echo JOB FINISHED!
date $format
