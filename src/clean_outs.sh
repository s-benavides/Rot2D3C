#!/bin/sh

####
# Script for cleaning outputs, to reduce storage issues.
####

echo run?
read runname

tail ./$runname/run/time_field.txt

echo latest output?
read lb
echo delete down to what output?
read upto

numlist () {
PYTHON_ARG="$lb" PYTHON_ARG2="$upto" PYTHON_ARG3="$runname" python - <<END
import os
la = int(os.environ['PYTHON_ARG'])
upto = int(os.environ['PYTHON_ARG2'])
run = str(os.environ['PYTHON_ARG3'])
while la>(upto+10):
	s="rm ./"+run+"/outs/*{%03d..%03d}.out;" % (la-10,la-2)
	print(s)
	la = la-10
END
}

TEST=$(numlist $lb $upto)
echo $TEST
eval $TEST

####################################

##!/bin/sh
##echo Latest input?
##read last
##while ((last>0));
##do
##	echo {"%03d\n" $(($last-10)).."%03d\n" $(($last-2))}
##	last=$((last-10))
##done
