#!/bin/bash

# first number is number of cores
# second number is J 
mpiexec -machinefile lab_ma_file -np 16 ./main_exe 32
