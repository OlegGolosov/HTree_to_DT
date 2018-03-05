#!/bin/bash

day=108
user=$(whoami)   
currentDir=$(pwd)

jobscript=${currentDir}/jobScript.sh
currentTime=$(date +"%b_%d_%k_%M" | tr -d ' ')
configuration=${1}
filelist=${currentDir}/../lists/${configuration}
outputdir=${currentDir}/../output/${currentTime}/${configuration}
logdir=${currentDir}/../log/${currentTime}/${configuration}

. /cvmfs/hades.gsi.de/install/5.34.34/hydra2-4.9u/defall.sh 

if [ "$#" -ne "2" ]
then
    echo "All runs from $filelist will be submitted in "
    echo "3..."
    sleep 1
    echo "2..."
    sleep 1
    echo "1..."
    sleep 1
    maxNumberOfRuns=100000
else
    maxNumberOfRuns=${2}
fi

bash ./make.sh

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

cp -r ${currentDir}/../src ${currentDir}/../output/${currentTime}/${configuration}/src
cp ${currentDir}/../exeMakeTree ${currentDir}/../output/${currentTime}/${configuration}


env=/cvmfs/hades.gsi.de/install/5.34.34/hydra2-4.9u/defall.sh   # environment script
exe=${currentDir}/../output/${currentTime}/${configuration}/exeMakeTree             # executable
numberOfRuns=0   # counter for part number

prevrunnumber=$(cat $filelist | sort | sed -r 's/.*\/be//' | sed -r 's/.hld.+//' | rev | cut -c 3- | rev | head -1)
runlist=" "
for file in $(cat $filelist | sort)
do
    runnumber=$(echo $file | sed -r 's/.*\/be//' | sed -r 's/.hld.+//' | rev | cut -c 3- | rev)
    if [ $prevrunnumber -eq $runnumber ]
    then
        if [ "$runlist" = " " ]
        then
            runlist=${file}
        else
            runlist="${runlist},${file}"
        fi
    else
        log_err=${logdir}/${prevrunnumber}.err
        log_out=${logdir}/${prevrunnumber}.out
        command="sbatch --error=${log_err} --output=${log_out} ${jobscript} ${env} ${exe} ${runlist} ${outputdir}/tree_${prevrunnumber}.root -1"
        echo " "
        echo "sbatch"
        echo "--error=${log_err} --output=${log_out}"
        echo "${jobscript} ${env} ${exe}"
        echo "${runlist}"
        echo "${outputdir}/tree_${prevrunnumber}.root -1"
        $command
        ((numberOfRuns+=1))
        if [ ${numberOfRuns} -eq ${maxNumberOfRuns} ]
        then
            echo " "
            echo "${numberOfRuns} jobs have been submitted"
            break
        fi
        runlist=${file}
    fi
    prevrunnumber=$runnumber
done

echo " "
echo "${numberOfRuns} jobs have been submitted"
