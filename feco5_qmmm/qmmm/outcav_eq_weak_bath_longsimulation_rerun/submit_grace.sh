#!/bin/bash

#SBATCH --job-name=aimd_FeCO5_16
#SBATCH --time=28-00:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=16
#SBATCH --partition=general
#SBATCH -p pi_hammes_schiffer

source ~/.bashrc

PATTERN='<step>80000</step>'

CHECKPOINT=simu_1.checkpoint
if grep -q $PATTERN $CHECKPOINT; then
    echo "found checkpoint finished"
    echo "Skip this job and quit"
else
    echo "not found checkpoint finished"
    a=$(wc -c < $CHECKPOINT)
    if [ ! -f "$CHECKPOINT" ] || [ $a -le 1 ]; then
       echo "Performing eq simulation for 1 ps"
       i-pi input_1.xml &> log_ipi_1
    else
       echo "Continuing the simulation of 1 ps"
       i-pi $CHECKPOINT &> log_ipi_1
    fi
fi

cp RESTART $CHECKPOINT
