#!/bin/sh
echo new run name?
read newrun
echo copy from where?
read oldrun

mkdir $newrun

cd $newrun

mkdir run outs ins

cd ../

rsync -aP --exclude='*.txt' $oldrun/run/* $newrun/run/

cd $newrun/run/

ls

./clean_new_run.sh

echo done
