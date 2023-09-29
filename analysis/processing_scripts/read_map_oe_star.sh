#!/bin/bash
#SBATCH --partition=3day-long

#define input
directory='/home/jdubos2/mtstp/data/raw_sequences'
genome_index='/home/jdubos2/mtstp/data/oe_genome/star_index'
outdir='/home/jdubos2/mtstp/data/oe_star_alignment_results'

#get input files
forward_files=`ls $directory/*[ci]i*_1.fq.gz`
#get list of completed runs, incase run kills for some reason
completed_runs=`ls $outdir | cut -f1 -d_ | uniq`

for forward_file in ${forward_files[@]}; do

    #define file handle
    prefix=$(basename ${forward_file})
    prefix=`echo $prefix | cut -f1 -d_`

    #check to see if run has been completed
    run_complete=0
    for run in ${completed_runs[@]}; do
        if [ "$prefix" = "$run" ]; then
            run_complete=1
        fi
    done

    if [ $run_complete = 0 ]; then

        #define reverse file
        reverse_file=`echo $forward_file | sed 's/_1/_2/g'`
        prefix=`echo $prefix`_

        echo "Working on $prefix"

        #run STAR
        STAR --runMode alignReads \
        --outSAMtype BAM Unsorted \
        --runThreadN 20 \
        --readFilesCommand zcat \
        --genomeDir $genome_index \
        --outFileNamePrefix $outdir/$prefix \
        --readFilesIn $forward_file $reverse_file \
        --quantMode GeneCounts

        echo "Completed $prefix"
    fi
done