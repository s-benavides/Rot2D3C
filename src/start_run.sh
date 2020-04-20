#!/bin/sh
set -eu   ## Stop on errors and on undefined variables
read -p "Enter folder names you wish to start, separated by 'space' : " input

for i in ${input[@]}
do
   echo ""
   echo "Working on : "$i    # or do whatever with individual element of the array
   echo ""
	read -n 1 -p "See status? Press any key"
	vi $i/run/status.inp
	read -n 1 -p "See parameter? Press any key"
	vi $i/run/parameter.inp
	read -n 1 -p "See jobscript? Press any key"
	vi $i/run/jobscript.sh
	read -n 1 -p "Submit job?"
	cd $i/run
	sbatch jobscript.sh
	cd ../../
done

echo ""
echo done
