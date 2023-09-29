#!/bin/bash

#define relative path to raw sequences
oe_read_ids='/home/jdubos2/mtstp/data/oe_aligned_reads'
dpl_read_ids='/home/jdubos2/mtstp/data/dpl_aligned_reads'
outdir='/home/jdubos2/mtstp/data/consolidated_aligned_read_ids'

#get input files
oe_files=`ls $oe_read_ids`

#get list of completed runs, incase run kills for some reason
completed_runs=`ls $outdir | cut -f1 -d_ | uniq`

for file in ${oe_files[@]}; do

    #define file handle
    prefix=$(basename ${file})
    prefix=`echo $prefix | cut -f1 -d_`

    #check to see if run has been completed
    run_complete=0
    for run in ${completed_runs[@]}; do
        if [ "$prefix" = "$run" ]; then
            run_complete=1
        fi
    done

    if [ $run_complete = 0 ]; then

        file_name=`echo $prefix`_reads_mapped.txt
        
        echo "Working on $prefix"

        cat $oe_read_ids/$file_name $dpl_read_ids/$file_name > $outdir/$file_name

    fi
done
        