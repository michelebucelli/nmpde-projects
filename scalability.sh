#!/bin/bash

MAX_PROCESSES=10

mkdir -p runs

for (( CORES=1; CORES<=$MAX_PROCESSES; CORES+=1 )) do
        cp -f run.orig runs/${CORES}.sub
        sed -i "s/CORES/$CORES/g" runs/${CORES}.sub
        qsub runs/${CORES}.sub
    done
