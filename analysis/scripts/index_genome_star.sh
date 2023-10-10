#!/bin/bash

#define genome fasta path
genome_fasta='/home/jdubos2/mtstp/data/dpl_genome/genome.fasta'
#define genes file
genes_gtf='/home/jdubos2/mtstp/data/dpl_genome/genes.gtf'
#define output dir
index_dir='/home/jdubos2/mtstp/data/dpl_genome/star_index'

STAR -- runThreadN 4 \
--runMode genomeGenerate \
--genomeDir $index_dir \
--genomeFastaFiles $genome_fasta \
--sjdbGTFfile $genes_gtf \
--genomeSAindexNbases 12 \
--sjdbOverhang 149