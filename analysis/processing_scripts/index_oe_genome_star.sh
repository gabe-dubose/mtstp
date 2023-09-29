#!/bin/bash

#define genome fasta path
genome_fasta='/home/jdubos2/mtstp/data/oe_genome/genome.fasta'
#define genes file
genes_gff='/home/jdubos2/mtstp/data/oe_genome/genes.gff'
#define output dir
index_dir='/home/jdubos2/mtstp/data/oe_genome/star_index'

STAR -- runThreadN 4 \
--runMode genomeGenerate \
--genomeDir $index_dir \
--genomeFastaFiles $genome_fasta \
--sjdbGTFfile $genes_gff \
--sjdbGTFtagExonParentTranscript Parent\
--genomeSAindexNbases 10 \
--sjdbOverhang 149