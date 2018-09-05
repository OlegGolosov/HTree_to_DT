#!/bin/bash

current_dir=$(pwd)
n_events=-1

partition=main
[ "$partition" == "debug" ] && time=0:20:00
[ "$partition" == "main" ] && time=8:00:00

current_time=$(date +"%b_%d_%k_%M" | tr -d ' ')
file_list=${1}
configuration=$(basename "$file_list")
output_dir=${current_dir}/../output/${configuration}/${current_time}
log_dir=${output_dir}/log
lists_dir=${output_dir}/lists

hadesroot=/cvmfs/hades.gsi.de/install/6.12.06/hydra2-4.9w/defall.sh 

mkdir -p $output_dir
mkdir -p $lists_dir
mkdir -p $log_dir

csplit -s -f "$lists_dir/" -b %1d.list $file_list /01./ {*}
rm $lists_dir/0.list

if [ "$#" -ne "2" ]
then
    n_runs=$(ls $lists_dir/*.list | wc -l)
else 
    n_runs=${2}
fi

job_range=1-$n_runs

cp -r ${current_dir}/../src ${output_dir}/src
executable=${current_dir}/../build/HTree_to_DT

echo configuration=$configuration
echo executable=$executable
echo output_dir=$output_dir
echo lists_dir=$lists_dir
echo log_dir=$log_dir
echo n_runs=$n_runs
echo job_range=$job_range

sbatch -J HTree_to_DT -p $partition -t $time -a $job_range -e ${log_dir}/%A_%a.e -o ${log_dir}/%A_%a.o --export=executable=$executable,output_dir=$output_dir,lists_dir=$lists_dir,hadesroot=$hadesroot,n_events=$n_events batch_run.sh
