#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Name: clean_scaffold.py

Date: Feb 28 11:25:02 2022

@author: Elena Aramendia

This script takes a fasta file and checks length and GC-content of each sequence,
filtering the sequences according to this parameters.
Minimum length and max GC-content (in %) can be entered as optional parameters,
if not specified the default values are 30% GC-content (sequences with a higher
GC content will be removed) and 3000 minimum length (sequences shorter than 3000
nucleotides will be removed).
Created for the Malaria Project on BINP29.
   
Modules:
    argparse
    
Functions:
    - get_fasta: function that takes a fasta file and returns a dictionary with header
    as the key and sequence as the ??, gets multiple-line sequences in a single line.

USAGE:
    python clean_scaffold.py -i fasta_file.fa [-len length in nucleotides] [-gc % GC-content] [-o output_file.fa]
    
Example use:
    python clean_scaffold.py -i fasta_file.fa [-len length in nucleotides] [-gc % GC-content] [-o output_file.fa]
"""
#%% Functions
# Define function to make single-line sequences from a fasta file, in case they are in multiple lines
def get_fasta(fasta_file):
    # Define lists of sequences, sequence fragments and headers
    seq = []
    seq_fragments = []
    header_list = []
    fasta_dict = {}
    index = 0 # counter
    fasta = False # Check that is a valid FASTA file
    for line in fasta_file:
    # Search for the header line
        if line.startswith('>'):
            fasta = True
            header = line.rstrip() # Remove newline character
            header_list.append(header)
            #If this is NOT the first sequence
            if seq_fragments:
                #Add sequence to a list of sequences
                curr_seq = ''.join(seq_fragments)
                seq.append(curr_seq)
            seq_fragments = []
        else:
            # Found more of existing sequence
            curr_seq = line.rstrip() # remove new line character
            seq_fragments.append(curr_seq) # Add fragment to current sequence   
    if fasta == False:
        raise Exception('Not a valid format. Please try again with a FASTA file.')
    # Appending the last fasta sequence 
    if seq_fragments:
        #if the file is not empty
        curr_seq = ''.join(seq_fragments)
        seq.append(curr_seq)
    for i in header_list: # Create dictionary with all the sequences
        fasta_dict[i] = seq[index]
        index += 1
    return fasta_dict
#%% Argument parser
import argparse 
parser = argparse.ArgumentParser()
usage = 'This program takes an input fasta file and removes sequences above a certain GC-content and with a certain length'

parser.add_argument(                                                  # fasta file
    '-i',
    type = argparse.FileType('r'),
    dest = 'fasta_file',
    required = True,
    help = 'Input fasta file'
    )

parser.add_argument(                                                  # minimun length
    '-len',
    type = int,
    dest = 'min_scaff',
    required = False,
    default = 3000,
    help = 'Minimun length for the sequences'
    )

parser.add_argument(                                                  # gc content
    '-gc',
    type = float,
    dest = 'max_GC',
    required = False,
    default = 30,
    help = 'Cutoff for GC-content (in %), sequences with GC% above this will be removed'
    )

parser.add_argument(                                                  # output file
    '-o',
    type = argparse.FileType('w'),
    dest = 'output_file',
    required = False,
    default = 'clean_seq.fa',
    help = 'Output file'
    )



args = parser.parse_args() 

#%%

fasta = get_fasta(args.fasta_file)

for key in fasta:    # for each header (key in the dictionary)     
        header = key.rstrip() # Store header, strip newline
        scaff = key.split('  ')[0][1:]
        length = key.split('=')[1] # length of the scaffold
        seq = fasta[key] 
        useq = seq.upper() # make sequence uppercase
        if len(useq) >= int(args.min_scaff): # if the length is okay, count gc content and print to outfile
                # count gc content
                ATcount = useq.count('A') + useq.count('T')
                GCcount = useq.count('C') + useq.count('G')
                gc = GCcount*100 / (GCcount+ATcount)
                if gc <= float(args.max_GC):
                    print('>' + scaff, file=args.output_file)
                    print(useq, file=args.output_file)

 

            