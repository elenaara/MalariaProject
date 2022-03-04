#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Name: remove_scaffolds.py

Created on Tue Mar  1 10:11:58 2022

@author: Elena Aramendía

This script removes scaffolds containing bird genes from the Ht fasta file.
Created for the Malaria Project on BINP29.

Modules:
    argparse
    

Example usage:
    python remove_scaffolds.py -i Ht2.fna -s scaffolds.txt -o only_ht.fna
"""
#%% Argument parser
import argparse 
parser = argparse.ArgumentParser()
usage = 'This script is created to retrieve scaffolds that contain host genes, identified by previously comparing the BLAST-hits with the swissprot database.'

parser.add_argument(                                                  # fasta file
    '-i',
    type = argparse.FileType('r'),
    dest = 'fasta_file',
    required = True,
    help = 'Input fasta file'
    )

parser.add_argument(                                                  # scaffolds to remove
    '-s',
    type = argparse.FileType('r'),
    dest = 'scaff_file',
    required = True,
    help = 'List of scaffolds to remove, each scaffold name should be in one line'
    )



parser.add_argument(                                                  # output file
    '-o',
    type = argparse.FileType('w'),
    dest = 'outfile',
    required = False,
    default = 'output',
    help = 'Output file'
    )



args = parser.parse_args() 

#%% test files
#scaff_file=open('scaffolds2.txt', 'r')
#fasta_file=open('/clean_H_tartak/H_tartakovsky_clean30gc.fa','r')

## BIRD SCAFFOLDS FILE
# Get scaffolds to remove in a list
bird_scaffolds = []
for line in args.scaff_file:
    bird_scaffolds.append(line.rstrip())
    
## FASTA FILE
# Parse to file and remove scaffolds present in the list
# If the scaffold is NOT in the list, id and sequence are kept and later printed

fasta_dict = {}
length_dict = {}
for line in args.fasta_file:
    if line.startswith('>'):
        scaff = line.rstrip().split()[0][1:]
        
        if scaff not in bird_scaffolds:
            fasta_dict[scaff] = ''
            
    else:
        if scaff not in bird_scaffolds:
            fasta_dict[scaff] += line
    
for key in fasta_dict:
    print('>' + key, file=args.outfile)
    print(fasta_dict[key], file=args.outfile)


