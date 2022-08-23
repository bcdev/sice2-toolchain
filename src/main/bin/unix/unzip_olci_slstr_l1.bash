#!/bin/bash

year=$1
month=$2
startDay=$3
endDay=$4

subdir=olci_slstr_l1

numDaysInMonth=`cal $month $year | awk 'NF {DAYS = $NF}; END {print DAYS}'`
olddir=${PWD}
for day in `seq -w $startDay $endDay`; do
  cd ${olddir}/${subdir}/$year/$month/$day
  #for i in `ls S3A_OL*.zip | sed -r 's/^.+\///'`; do
  for i in `ls S3A_OL*.zip`; do
    echo "unzip $i"
    unzip $i
    echo "rm -Rf $i"
    rm -Rf $i
  done
  for i in `ls S3A_SL*.zip`; do
    echo "unzip $i"
    unzip $i
    echo "rm -Rf $i"
    rm -Rf $i
  done
done
cd $olddir

