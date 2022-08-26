#!/usr/bin/env bash 

conda activate SICE

#./test_redir.bash
#./S3_proc_test.sh
./S3_proc_test.sh -i ./olci_slstr_l1/2021/08/02 -o ./olci_slstr_preprocessed/2021-08-02 -X ./S3.xml || error=true
