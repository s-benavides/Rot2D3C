#!/bin/sh
set -eu   ## Stop on errors and on undefined variables
read -p "Enter folder names you wish to restart, separated by 'space' : " input

for i in ${input[@]}
do
   echo ""
   echo "Working on : "$i    # or do whatever with individual element of the array
   echo ""
	while IFS=" " read -r value1 remainder
		do
  		lastout=$value1
		timef=$remainder
		done < "$i/run/time_field.txt"
	echo "Last output: "$lastout", removing old files and copying new ones to 'ins' folder"
	rm $i/ins/*.out
	rsync -aP $i/outs/*.$lastout.out $i/ins/
	echo "Done"
	echo ""
	echo "Working on : "$i    # or do whatever with individual element of the array
	echo ""
	echo $lastout $timef 
	read -n 1 -p "Continue?"
	vi $i/run/status.inp
	read -n 1 -p "Submit job?" job
	cd $i/run
	sbatch jobscript.sh
	cd ../../
done

echo ""
echo done
