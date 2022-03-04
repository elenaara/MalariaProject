#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Name: genomic_gc.py

Created on Tue Mar  1 11:39:05 2022

@author: Elena Aramendia

This script calculates total GC-content of a fasta file.
It counts GCs and ATs (does not take into account other symbols, i. e Ns)
and calculates GC content expressed in percentage.

Modules:
    sys
    
Usage:
    genomic_gc.py fasta_file.fna

"""


#%% arguments
import sys

try:
    fastafile = sys.argv[1]

except:
    print("Run program like: genomic_gc.py fasta_file.fna")
    sys.exit()

#%%
gc = 0
at = 0
GCcontent = 0
with open(fastafile, "r") as fasta:
    for line in fasta:
        if line.startswith('>'):
            pass
        else:
            at += line.count('A') + line.count('T')
            gc += line.count('C') + line.count('G')
    GCcontent = gc*100 / (gc + at)
    print('Genomic GC content: {:.2f}%'.format(GCcontent))


