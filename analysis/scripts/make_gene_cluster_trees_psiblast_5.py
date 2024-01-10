#!/usr/bin/env python3

from gnat import misc_utils
from gnat import phylogenetics


#load clusters and partition each cluster into an individual sequence file
clusters = '/home/gabe/Desktop/mtstp/data/intermediate_data/gene_cluster_diversity_analysis/psiblast_id-30_e-neg-10_cov-1_sequence_clusters.json'
sequences = '/home/gabe/Desktop/mtstp/analysis/data/dplex_protein_db/coding_sequences.fasta'
outdir = '/home/gabe/Desktop/mtstp/analysis/data/gene_family_analysis/gene_clusters_sequences/psiblast_inferred'

misc_utils.partition_sequence_clusters(clusters, sequences, outdir)

#define files and directories
directory = '/home/gabe/Desktop/mtstp/analysis/data/gene_family_analysis/gene_clusters_sequences/psiblast_inferred'
extension = 'fasta'
alignment_dir = '/home/gabe/Desktop/mtstp/analysis/data/gene_family_analysis/gene_clusters_sequences/psiblast_inferred_alignments'
phylogeny_dir = '/home/gabe/Desktop/mtstp/analysis/data/gene_family_analysis/gene_clusters_sequences/psiblast_inferred_phylogenies'

phylogenetics.batch_alignments(directory, extension, alignment_dir)
phylogenetics.batch_phylogenies(alignment_dir, extension, phylogeny_dir)