#!/bin/bash

echo " $(date +%m%d-%H%M).out "
#./../build/bin/Density.exe > $(date +%m%d-%H%M).out

nsys profile --output $(date +%m%d-%H%M) mpirun -np 2 ./../build/bin/Density.exe > $(date +%m%d-%H%M)
