#!/bin/bash

day=108
user=$(whoami)   
currentDir=$(pwd)
n_events=10

partition=main
partition=debug
[ "$partition" == "debug" ] && time=0:20:00
[ "$partition" == "main" ] && time=8:00:00

currentTime=$(date +"%b_%d_%k_%M" | tr -d ' ')
file_list=${1}
configuration=$(basename "$file_list")
outputdir=${currentDir}/../output/${currentTime}/${configuration}
logdir=${currentDir}/../log/${currentTime}/${configuration}

hadesroot=/cvmfs/hades.gsi.de/install/6.12.06/hydra2-4.9w/defall.sh 

if [ "$#" -ne "2" ]
then
    echo "All runs from $filelist will be submitted in "
    n_runs=100000
else
    n_runs=${2}
fi

if [ ! -d $outputdir ]
then
   echo "===> CREATE OUTPUTDIR : $outputdir"
   mkdir -p $outputdir
else
   echo "===> USE OUTPUTDIR : $outputdir"
fi

if [ ! -d $logdir ]
then
   echo "===> CREATE LOGDIR : $logdir"
   mkdir -p $logdir
else
   echo "===> USE LOGDIR : $logdir"
fi

cd $outputdir;
csplit -f list_ -b %1d.list $file_list /_1/ {*}
rm list_0
n_runs=$(ls *.list | wc -l)
jobRange=1-$n_runs

cp -r ${currentDir}/../src ${currentDir}/../output/${currentTime}/${configuration}/src
executable=${currentDir}/../build/HTree_to_DT

sbatch -J HTree_to_DT -p $partition -t $time -a $jobRange -e ${logdir}/%A_%a.e -o ${logdir}/%A_%a.o --export=executable=$executable,outputdir=$outputdir,hadesroot=$hadesroot,n_events=$n_events run_conv.sh