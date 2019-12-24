#!/usr/bin/env python
# coding: utf-8

from Bio import SeqIO
import sys

query_fasta = sys.argv[2]  # fasta file of protein sequences used in BLAST
blast = query_fasta = sys.argv[3]  # r'D:\Геномный_центр\files\HaCENH3_vs_orf_1KP_CENH3.fasta'
with open(blast) as inFile, open('selectedORF_CENH3.faa', 'w') as outFile:
    ind = SeqIO.index(query_fasta, 'fasta')
    species_list = []
    for line in inFile:
        sp = line.split('\t')
        seq = ind[sp[1]]
        seq.id = seq.description.split(' ')[-1]
        seq.description = ''
        species, qstart, qend = seq.id, int(sp[7]), int(sp[8])
        if abs(qstart - qend) > 143 * 0.8 and species not in species_list:
            species_list.append(species)

            SeqIO.write(seq, outFile, 'fasta')
    print(f"Number of sequences selected: {len(species_list)}")