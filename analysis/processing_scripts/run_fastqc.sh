#!/bin/bash

#define relative path to raw sequences
directory='../data/raw_sequences'
outdir='../data/fastqc_reports'

#get input files
files=`ls $directory`

#run fastqc on files
for file in ${files[@]}; do
    echo "Running FastQC on: $directory/$file"
    fastqc $directory/$file -o $outdir
    echo "Completed: $directory/$file"
    echo "Saved report to: $outdir"
done