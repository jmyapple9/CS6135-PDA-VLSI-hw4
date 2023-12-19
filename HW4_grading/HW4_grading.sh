#!/bin/bash

# Set the variables
homeworkName=HW4
binaryName=hw4
timeLimit=600
testcasePool="public1 public2 public3 public4 public5 public6"
# testcasePool="public1 public2 public3"

# Set the directories and paths
root=$(pwd)
outputDir=$root/output
studentDir=$root/student
testcaseDir=$root/testcase
verifyBin=$root/verifier/verify
chmod +x "$verifyBin"
csvFile=$root/${homeworkName}_grade.csv

# Define the execute function for running tests
function execute() {
	if [[ ! -f "$binaryName" ]]; then
		echo "N/A"
	else
		local log=$(timeout $timeLimit time -p ./$binaryName "$testcaseDir/$1/$1.aux" "$outputDir/$1.gp.pl" 2>&1 >/dev/null)
		if [[ $log =~ "real " && $log =~ "user " && $log =~ "sys " ]]; then
			echo "$(echo "$log" | grep real | tail -1 | cut -d ' ' -f 2)"
		else
			echo "TLE"
		fi
	fi
}

# Define the verify function for result verification
function verify() {
	if [[ ! -s "$outputDir/$1.gp.pl" ]]; then
		echo "The output file of $1 is not found or is empty."
	else
		local log=$("$verifyBin" "$testcaseDir/$1/$1.aux" "$outputDir/$1.gp.pl")
		if [[ $log =~ "[Error]" ]]; then
			echo "There is an error in the output results of $1 ($(echo "$log" | grep "[Error]" | tail -n 1))."
		elif [[ $log =~ "[Success]" ]]; then
			local wirelength=$(echo "$log" | grep "HPWL:" | tail -n 1)
			if [[ $wirelength =~ HPWL:\ ([0-9]+) ]]; then
				wirelength="${BASH_REMATCH[1]}"
				echo "success $wirelength"
			else
				echo "Something goes wrong when verifying the output results of $1."
			fi
		else
			echo "Something goes wrong when verifying the output results of $1."
		fi
	fi
}

echo "|------------------------------------------------|"
echo "|                                                |"
echo "|    This script is used for PDA $homeworkName grading.    |"
echo "|                                                |"
echo "|------------------------------------------------|"

# Set CSV title
csvTitle="student id,tar.gz,file structure,readme,makefile,make clean,make"
for testcase in $testcasePool; do
	csvTitle="$csvTitle,$testcase wirelength,$testcase runtime"
done
echo "$csvTitle,status" >"$csvFile"

# Traverse the student directory for grading
cd "$studentDir"
for studentId in *; do
	if [[ -d "$studentId" ]]; then
		cd "$studentId"
		printf "grading on %s:\n" "$studentId"
		csvContent="$studentId"

		# Check various requirements
		correctTarGz=no
		correctStruct=no
		haveReadme=no
		haveMakefile=no
		correctMakeClean=no
		correctMake=no
		if [[ -f "CS6135_${homeworkName}_$studentId.tar.gz" ]]; then
			correctTarGz=yes
			find . -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +
			tar -zxf "CS6135_${homeworkName}_$studentId.tar.gz"
			if [[ -d "$homeworkName/src" && -d "$homeworkName/bin" ]]; then
				correctStruct=yes
				cd $homeworkName/src
				if [[ $(find . -maxdepth 1 -type f -iname README* | wc -l) -eq 1 ]]; then
					haveReadme=yes
				fi
				if [[ $(find . -maxdepth 1 -type f -iname Makefile | wc -l) -eq 1 ]]; then
					haveMakefile=yes
					if [[ ! -f "../bin/$binaryName" ]]; then
						touch "../bin/$binaryName"
					fi
					make clean >/dev/null
					if [[ ! -f "../bin/$binaryName" ]]; then
						correctMakeClean=yes
					else
						rm -f ../bin/$binaryName
					fi
					make >/dev/null
					if [[ -f "../bin/$binaryName" ]]; then
						correctMake=yes
					fi
				fi
				cd ../bin
			fi
		fi
		printf " checking item          | status\n"
		echo "------------------------|--------"
		printf " correct tar.gz         | %s\n" $correctTarGz
		printf " correct file structure | %s\n" $correctStruct
		printf " have README            | %s\n" $haveReadme
		printf " have Makefile          | %s\n" $haveMakefile
		printf " correct make clean     | %s\n" $correctMakeClean
		printf " correct make           | %s\n" $correctMake
		printf "\n"
		csvContent="$csvContent,$correctTarGz,$correctStruct,$haveReadme"
		csvContent="$csvContent,$haveMakefile,$correctMakeClean,$correctMake"

		# Verify the program by each testcase
		rm -rf "$outputDir"/*
		printf "%10s | %10s | %10s | %s\n" testcase wirelength runtime status
		echo "-----------|------------|------------|--------"
		statusList=""
		for testcase in $testcasePool; do
			printf "%10s | " $testcase
			wirelength=N/A
			runtime=$(execute $testcase)
			if [[ $runtime == "N/A" ]]; then
				status="The executable file ($binaryName) is not found."
				statusList=" $status"
			elif [[ $runtime == "TLE" ]]; then
				status="Time out while testing $testcase."
				statusList="$statusList $status"
			else
				status=$(verify $testcase)
				if [[ $status =~ "success " ]]; then
					IFS=' ' read -r status wirelength <<<"$status"
				else
					runtime=N/A
					statusList="$statusList $status"
				fi
			fi
			printf "%10s | %10s | %s\n" $wirelength $runtime "$status"
			csvContent="$csvContent,$wirelength,$runtime"
		done
		printf "\n"
		echo "$csvContent,$statusList" >>"$csvFile"

		cd "$studentDir"
	fi
done

echo "|-----------------------------------------------------|"
echo "|                                                     |"
echo "|    Successfully generate grades to ${homeworkName}_grade.csv    |"
echo "|                                                     |"
echo "|-----------------------------------------------------|"
