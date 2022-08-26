#!/usr/bin/env bash

# Wrapper for running SICE pipeline

echo off

#set -o errexit 
#set -o nounset 
#set -o pipefail 
#set -x 

red='\033[0;31m' 
orange='\033[0;33m' 
green='\033[0;32m' 
nc='\033[0m' # No Color 
log_info() { echo -e "${green}[$(date --iso-8601=seconds)] [INFO] ${*}${nc}"; } 
log_warn() { echo -e "${orange}[$(date --iso-8601=seconds)] [WARN] ${*}${nc}"; } 
log_err() { echo -e "${red}[$(date --iso-8601=seconds)] [ERR] ${*}${nc}" 1>&2; }

# wvcci-vm
preprocessed_root=${SICE2_ROOT}/olci_slstr_preprocessed
mosaic_root=${SICE2_ROOT}/olci_slstr_mosaics
sice_results_root=${SICE2_ROOT}/olci_slstr_sice_results

# Slope correction
slopey=false 

xml_file=${SICE2_ROOT}/S3.xml 
LD_LIBRARY_PATH=. # SNAP requirement 

error=false
startDay=$1
endDay=$2

for year in {2021..2021}; do
  for month in {08..08}; do
    numDaysInMonth=`cal $month $year | awk 'NF {DAYS = $NF}; END {print DAYS}'`
    #for day in `seq -w 01 ${numDaysInMonth}`; do
    for day in `seq -w $startDay $endDay`; do
      thedate=${year}-${month}-${day}
      if [[ -d "${mosaic_root}/${thedate}" ]] && [[ -e "${mosaic_root}/${thedate}/conc.tif" ]]; then
        log_warn "${mosaic_root}/${thedate} already exists, date skipped"
        continue
      fi

      mkdir -p ${preprocessed_root}/${thedate}
      mkdir -p ${mosaic_root}/${thedate}
      mkdir -p ${sice_results_root}/${thedate}

      # SNAP: Reproject, calculate reflectance, extract bands, etc.
      # e.g. ./S3_proc.sh -i ./olci_slstr_l1/2021/08/01 -o ./olci_slstr_preprocessed/2021-08-01 -X ./S3.xml || error=true
      conda activate base
      echo "./S3_proc.sh -i ${SEN3_local}/${year}/${month}/${day} -o ${preprocessed_root}/${thedate} -X ${xml_file} -t || error=true"
      ./S3_proc.sh -i ${SEN3_local}/${year}/${month}/${day} -o ${preprocessed_root}/${thedate} -X ${xml_file} -t || error=true

      # Run the Simple Cloud Detection Algorithm (SCDA)
      # Works with SICE env only
      # --> run in advance: conda activate SICE
      conda activate SICE
      echo "python ./SCDA.py ${preprocessed_root}/${thedate} || error=true"
      python ./SCDA.py ${preprocessed_root}/${thedate} || error=true

      # Mosaic
      # Works with base env only
      # --> run in advance: conda activate base
      conda activate base
      echo "./dm.sh ${thedate} ${preprocessed_root}/${thedate} ${mosaic_root} || error=true"
      ./dm.sh ${thedate} ${preprocessed_root}/${thedate} ${mosaic_root} || error=true

      conda activate SICE
      if [ "$slopey" = true ]; then
        # Run the slopey correction
        echo "python ./get_ITOAR.py ${mosaic_root}/${thedate}/ $(pwd)/ArcticDEM/ || error=true"
        #python ./get_ITOAR.py ${mosaic_root}/${thedate}/ "$(pwd)"/ArcticDEM/ || error=true
      fi

      # SICE
      echo "python ./sice.py ${mosaic_root}/${thedate} || error=true"
      python ./sice.py ${mosaic_root}/${thedate} || error=true

      # gpt py_sice2_op -PtifInputDir=./olci_slstr_preprocessed/2021-08-01/20210801T133724 -PtifOutputDir=./olci_slstr_preprocessed/2021-08-01/20210801T133724/snap_output \ 
      #                 -f NetCDF4-BEAM -t ./delete_me.nc
      # Works with SICE env only
      # --> run in advance: conda activate SICE
      echo "gpt py_sice2_op -PtifInputDir=${mosaic_root}/${thedate} -PtifOutputDir=${sice_results_root}/${thedate} -f NetCDF4-BEAM -t ./delete_me.nc"
      gpt py_sice2_op -PtifInputDir=${mosaic_root}/${thedate} -PtifOutputDir=${sice_results_root}/${thedate} -f NetCDF4-BEAM -t ./delete_me.nc

      if [ "$error" = true ]; then
        echo "Processing of ${thedate} failed, please check logs."
      fi
    done
  done 
done
