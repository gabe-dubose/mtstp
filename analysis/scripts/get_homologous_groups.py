#!/usr/bin/env python3

import json
from collections import defaultdict

def combine_groups(groups):
    # Create a defaultdict to store groups with common elements
    combined_groups = defaultdict(list)

    # Iterate through each group
    for group_key, group_values in groups.items():
        # Check against existing combined groups
        found_in_combined = False
        for combined_key, combined_values in combined_groups.items():
            if any(val in combined_values for val in group_values):
                # If a common element is found, combine groups
                combined_groups[combined_key].extend(group_values)
                found_in_combined = True
                break

        if not found_in_combined:
            # If no common element is found, create a new combined group
            combined_groups[group_key] = group_values

    return combined_groups

def assemble_groups(pairwise_hits):
    #store which ids have been assigned to a group
    group_cache = {}
    #dictionary to store groups
    homologous_groups = {}
    #iterate through each entry
    group = 0
    for entry in pairwise_hits:
        #check if either pair has been assined to a group
        if entry[0] not in group_cache and entry[1] not in group_cache:
            #if not, make new group
            group_id = f"group_{group}"
            #add to cache
            group_cache[entry[0]] = group_id
            group_cache[entry[1]] = group_id
            #add to output dictionary
            homologous_groups[group_id] = [entry[0], entry[1]]
            #add to group counter
            group += 1
        #if either pair has been recorded
        elif entry[0] in group_cache or entry[1] in group_cache:
            #get which group(s) they were assined to
            if entry[0] in group_cache:
                group_id = group_cache[entry[0]]
            elif entry[1] in group_cache:
                group_id = group_cache[entry[1]]
            #add to cache
            group_cache[entry[0]] = group_id
            group_cache[entry[1]] = group_id
            #add to output dictionary
            if entry[0] not in homologous_groups[group_id]:
                homologous_groups[group_id].append(entry[0])
            if entry[1] not in homologous_groups[group_id]:
                homologous_groups[group_id].append(entry[1])

    #combine groups with overlap
    homologous_groups_collapsed = combine_groups(homologous_groups)
    #remove redundant values
    for group in homologous_groups_collapsed:
        homologous_groups_collapsed[group] = list(set(homologous_groups_collapsed[group]))

    return homologous_groups_collapsed

#-outfmt "6 qseqid sseqid pident length evalue bitscore qlen slen"
def get_sequence_clusters(file, evalue_cutoff, alignment_coverage_cutoff, percent_identity_cutoff):
    #initialize list to store homologous genes
    homologous_genes = []

    #read file
    with open(file, 'r') as infile:
        lines = infile.readlines()
    
    #iterate through lines
    for line in lines:
        try:
            #separate fields
            line = line.split('\t')
            query_seq_id = line[0]
            subject_seq_id = line[1]
            percent_identity = float(line[2])
            alignment_length = int(line[3])
            evalue = float(line[4])
            query_length = int(line[6])

            #calculate alignment coverage
            query_coverge = alignment_length / query_length

            #populate homologous genes list
            #don't include same identical sequence
            if query_seq_id != subject_seq_id:
                #filter at evalue cutoff
                if evalue <= evalue_cutoff:
                    #filter at alignment length
                    if query_coverge >= alignment_coverage_cutoff:
                        #filter at percent identity cutoff
                        if percent_identity > percent_identity_cutoff:
                            #if query sequence has not been recorded, add pair to list
                            pair = [subject_seq_id, query_seq_id]
                            homologous_genes.append(pair)

        except:
            pass

    homologous_groups = dict(assemble_groups(homologous_genes))

    return homologous_groups

#get gene families based on psiblast results
gene_clusters = get_sequence_clusters(file='/home/gabe/Desktop/mtstp/data/intermediate_data/gene_cluster_diversity_analysis/dplex_self_blast_results/psiblast_results_5.tsv',
                                                evalue_cutoff=1*10**-10,
                                                alignment_coverage_cutoff=1,
                                                percent_identity_cutoff=30)

#write to file
clusters_json = json.dumps(gene_clusters, indent=4)

with open("/home/gabe/Desktop/mtstp/data/intermediate_data/gene_cluster_diversity_analysis/psiblast_id-30_e-neg-10_cov-1_sequence_clusters.json", "w") as outfile:
    outfile.write(clusters_json)

#get gene families based on less stringent sequence similarity
gene_clusters = get_sequence_clusters(file='/home/gabe/Desktop/mtstp/data/intermediate_data/gene_cluster_diversity_analysis/dplex_self_blast_results/psiblast_results_5.tsv',
                                                evalue_cutoff=1*10**-5,
                                                alignment_coverage_cutoff=0.7,
                                                percent_identity_cutoff=20)

#write to file
clusters_json = json.dumps(gene_clusters, indent=4)

with open("/home/gabe/Desktop/mtstp/data/intermediate_data/gene_cluster_diversity_analysis/psiblast_id-20_e-neg-5_cov-0.7_sequence_clusters.json", "w") as outfile:
    outfile.write(clusters_json)