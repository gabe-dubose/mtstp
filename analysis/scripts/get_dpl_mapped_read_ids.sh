#!/bin/bash
#SBATCH --partition=3day-long

directory='/home/jdubos2/mtstp/data/star_alignment_results'
outdir='/home/jdubos2/mtstp/data/dpl_aligned_reads'

#get input files
bam_files=`ls $directory/*.bam`

#get list of completed runs, incase run kills for some reason
completed_runs=`ls $outdir | cut -f1 -d_ | uniq`

for bam_file in ${bam_files[@]}; do

    #define file handle
    prefix=$(basename ${bam_file})
    prefix=`echo $prefix | cut -f1 -d_`

    #check to see if run has been completed
    run_complete=0
    for run in ${completed_runs[@]}; do
        if [ "$prefix" = "$run" ]; then
            run_complete=1
        fi
    done

    if [ $run_complete = 0 ]; then

        outfile=`echo $prefix`_reads_mapped.txt

        echo "Working on $prefix"

        #get mapped reads
        samtools view -F 4 $bam_file | cut -f1 | uniq > $outdir/$outfile
    
        echo "Completed $prefix"
    fi
done