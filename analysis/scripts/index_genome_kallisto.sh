#!/bin/bash

index_path='/home/jdubos2/mtstp/data/dpl_genome/kallisto_index.idx'
genes_path='/home/jdubos2/mtstp/data/dpl_genome/cds.fasta'

kallisto index -i $index_path $genes_path
