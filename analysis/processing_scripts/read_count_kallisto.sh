#!/bin/bash
#SBATCH --partition=3day-long

#define relative path to raw sequences
directory='/home/jdubos2/mtstp/data/raw_sequences'
genome_index='/home/jdubos2/mtstp/data/dpl_genome/kallisto_index.idx'
outdir='/home/jdubos2/mtstp/data/dpl_kallisto_quantifications'
gtf='/home/jdubos2/mtstp/data/dpl_genome/genes.gtf'

#get input files
forward_files=`ls $directory/*_1.fq.gz`
#get list of completed runs, incase run kills for some reason
completed_runs=`ls $outdir | cut -f1 -d_ | uniq`

for forward_file in ${forward_files[@]}; do

    #make sample directory
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
        #make outdir
        mkdir $outdir/$prefix
        #define reverse file
        reverse_file=`echo $forward_file | sed 's/_1/_2/g'`
        
        kallisto quant \
        -i $genome_index \
        -o $outdir/$prefix \
        -b 100 \
        -t 5 \
        -g $gtf \
        $forward_file $reverse_file

    fi
done