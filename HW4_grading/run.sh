#!/bin/bash

# Define the number of times to run the command
num_runs=70

# Run the command in a loop
for ((i=1; i<=$num_runs; i++))
do
    echo "Running HW4_grading.sh - Attempt $i"
    bash HW4_grading.sh
    # echo "----------------------------------------"
done

echo "Script completed after $num_runs runs"
