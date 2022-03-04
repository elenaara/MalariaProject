# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 21:02:56 2022

@author: earam
"""

#%% Argument parser
import argparse
parser = argparse.ArgumentParser()
usage = 'This program takes a set of busco ids and a path to the location of the alignment files. It then concatenates all sequences from the same species into one sequence.'
parser.add_argument(                                                  # distance matrix
    '-i',
    type = argparse.FileType('r'),
    dest = 'ids',
    required = True,
    help = 'File containing BUSCO ids'
    )
parser.add_argument(                                                  # title
    '-d',
    type = str,
    dest = 'dir',
    required = True,
    help = 'Path to directory containing fasta files'
    )
parser.add_argument(                                                  # output file
    '-o',
    dest = 'output',
    required = False,
    default = 'concatenated _sequences.faa',
    help = 'Output file'
    )


args = parser.parse_args() 

#%%
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

#%%
# Open busco ids file
# For each id, find alignment file
# FOr each alignment file, add to dictionary sequence
seq_dict = {} # Dictionary for Species and Sequences (concatenate)

for name in args.ids: # Loop through BUSCO ids
    name = name.rstrip()
# We will loop through the alignment file, named after BUSCO ids
    if args.dir.endswith('/'): # check how the directory is written and add name of the file we are looking for
        file = args.dir + name + '_aligned.faa'
    else: 
        file = args.dir + '/' + name + '_aligned.faa'
    with open(file, 'r') as in1:
        fasta = get_fasta(in1) # get sequences in one line
    for key in fasta: # Loop through sequences
        sp = key.split()[0][1:]
        seq = fasta[key]
        if sp not in seq_dict: # Initialize dict entry if it is the first sequence
            seq_dict[sp] = seq
        else: # If it is not the first sequence add to dictionary in this species
            seq_dict[sp] += seq

# Print seqs to stdout
for key in seq_dict:
    print('>'+key)
    print(seq_dict[key])
        
        