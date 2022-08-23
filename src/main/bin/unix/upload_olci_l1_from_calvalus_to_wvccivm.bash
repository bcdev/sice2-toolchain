#!/bin/bash

year=$1
month=$2
startDay=$3
endDay=$4

subdir=MANIFEST_OLCI_202108

numDaysInMonth=`cal $month $year | awk 'NF {DAYS = $NF}; END {print DAYS}'`
#for day in `seq -w 02 ${numDaysInMonth}`; do
for day in `seq -w $startDay $endDay`; do
  for i in `ls ${PWD}/${subdir}/S3A_OL_1_EFR____$year$month$day* | sed -r 's/^.+\///'`; do
    l1name=${i:0:31}
    echo "sshpass -p "tasskaff" scp /calvalus/eodata/OLCI_FR_L1/v2/2021/08/$day/${l1name}*.SEN3.zip olafd@wvcci-vm:/home/olafd/sice2/workspace/geus-sice2/SICE2/olci_slstr_l1/$year/$month/$day"
    sshpass -p "tasskaff" scp /calvalus/eodata/OLCI_FR_L1/v2/2021/08/$day/${l1name}*.SEN3.zip olafd@wvcci-vm:/home/olafd/sice2/workspace/geus-sice2/SICE2/olci_slstr_l1/$year/$month/$day
  done
done

