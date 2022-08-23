#!/bin/bash

startDay=$1
endDay=$2

subdir=MANIFEST_SLSTR_202108

for year in {2021..2021}; do
  for month in {08..08}; do
    numDaysInMonth=`cal $month $year | awk 'NF {DAYS = $NF}; END {print DAYS}'`
    #for day in `seq -w 02 ${numDaysInMonth}`; do
    for day in `seq -w $startDay $endDay`; do
      for i in `ls ${PWD}/${subdir}/S3A_SL_1_RBT____$year$month$day* | sed -r 's/^.+\///'`; do
        echo $i
        lftp -c "open ftp://ftp.ceda.ac.uk --user \"odanne\" --password \"dAU>]aQB{z|B\"; cd neodc/sentinel3a/data/SLSTR/L1_RBT/$year/$month/$day; mget ${i}*NT_004.zip; bye"
        mv ${PWD}/${subdir}/*.zip ${PWD}/${subdir}/data
      done
    done
  done
done
