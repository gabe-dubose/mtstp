#!/bin/bash

psiblast -query /home/gabe/Desktop/mtstp/analysis/data/dplex_protein_db/protein_sequences.faa -db /home/gabe/Desktop/mtstp/analysis/data/dplex_protein_db/protein_sequences.faa -out /home/gabe/Desktop/mtstp/analysis/data/dplex_self_blast_results/psiblast_results_5.tsv -num_iterations 5 -outfmt "6 qseqid sseqid pident length evalue bitscore qlen slen"