#!/bin/bash

# Check if at least one argument is provided
# if [ $# -eq 0 ]; then
#     echo "Error: Missing argument. Usage: $0 <argument>"
#     exit 1
# fi

# arg1=$1

# if [ $arg1 -eq 1 ]; then
#     cd ~/CS6135-PDA-VLSI-hw4/
#     tar -zcvf CS6135_HW4_112065504.tar.gz HW4/
#     cd ~/CS6135-PDA-VLSI-hw4/HW4_grading/student/112065504
#     rm -rf *
#     cd ~/CS6135-PDA-VLSI-hw4/
#     mv CS6135_HW4_112065504.tar.gz HW4_grading/student/112065504
# elif [ $arg1 -eq 2 ]; then
#     cd ~/CS6135-PDA-VLSI-hw4/
#     tar -zcvf CS6135_HW4_112065504.tar.gz HW4/
#     cd ~/CS6135-PDA-VLSI-hw4/HW4_gd2/student/112065504
#     rm -rf *
#     cd ~/CS6135-PDA-VLSI-hw4/
#     mv CS6135_HW4_112065504.tar.gz HW4_gd2/student/112065504
# else
#     echo "Invalid arg: $arg1"
# fi

cd ~/CS6135-PDA-VLSI-hw4/
tar -zcvf CS6135_HW4_112065504.tar.gz HW4/
cd ~/CS6135-PDA-VLSI-hw4/HW4_grading/student/112065504
rm -rf *
cd ~/CS6135-PDA-VLSI-hw4/
mv CS6135_HW4_112065504.tar.gz HW4_grading/student/112065504