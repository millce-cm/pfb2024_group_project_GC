#!/usr/bin/env python

########################################### Set up our main script ########################################### 
import os, sys, re, argparse
import pandas as pd

parser = argparse.ArgumentParser(description= 'This is a script which has two optional arguments: file_1 and file_2. If you include at least one file, you will get basic statistics for that sequence. If you include both, you will get statistics for the whole genome and coding DNS')
parser.add_argument('-f1', '--file_1', help ='optional: Path to FASTA file for the whole genome')
parser.add_argument('-f2', '--file_2', help ='optional: Path to FASTA file for the coding sequence')
args = parser.parse_args()
file1 = args.file_1
file2 = args.file_2

########################################### These are the functions used in the script ########################################### 

# Function that a) splits a dna string into triplets and b) translates the triplets into amino acids
def amino_acid_translation (dna):  #This function reads in a dna string and translates it into a list of aminoacids

    dna_upper = dna.upper() # Converts everything into upper letters
    my_list = re.findall(r'(\S{3})', dna_upper) # Splits the string into a list, each element being a codon 
    i = 0 #This is needed for iterating through my_list

    while i < len(my_list): # The translation is based on the "Standard DNA codon Table" from Wikipedia (https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables)
        if re.search(r'TT[TC]', my_list[i]):
            my_list[i] = 'Phenylalanine'
        if re.search(r'TT[AG]', my_list[i]):
            my_list[i] = 'Leucine'
        if re.search(r'CT[TCAG]', my_list[i]):
            my_list[i] = 'Leucine'
        if re.search(r'AT[TCA]', my_list[i]):
            my_list[i] = 'Isoleucine'
        if re.search(r'ATG', my_list[i]):
            my_list[i] = 'Methionine'     
        if re.search(r'GT[TCAG]', my_list[i]):
            my_list[i] = 'Valine'
        if re.search(r'TC[TCAG]', my_list[i]):
            my_list[i] = 'Serine'
        if re.search(r'CC[TCAG]', my_list[i]):
            my_list[i] = 'Proline'
        if re.search(r'AC[TCAG]', my_list[i]):
            my_list[i] = 'Threonine'
        if re.search(r'GC[TCAG]', my_list[i]):
            my_list[i] = 'Alanine'
        if re.search(r'TA[TC]', my_list[i]):
            my_list[i] = 'Tyrosine'
        if re.search(r'CA[TC]', my_list[i]):
            my_list[i] = 'Histidine'
        if re.search(r'CA[AG]', my_list[i]):
            my_list[i] = 'Glutamine'
        if re.search(r'AA[TC]', my_list[i]):
            my_list[i] = 'Asparagine'
        if re.search(r'AA[AG]', my_list[i]):
            my_list[i] = 'Lysine'
        if re.search(r'GA[TC]', my_list[i]):
            my_list[i] = 'Aspartic acid'
        if re.search(r'GA[AG]', my_list[i]):
            my_list[i] = 'Glutamic acid'
        if re.search(r'TG[TC]', my_list[i]):
            my_list[i] = 'Cysteine'
        if re.search(r'TGG', my_list[i]):
            my_list[i] = 'Tryptophan'
        if re.search(r'CG[TCAG]', my_list[i]):
            my_list[i] = 'Arginine'
        if re.search(r'AG[TC]', my_list[i]):
            my_list[i] = 'Serine'
        if re.search(r'AG[AG]', my_list[i]):
            my_list[i] = 'Arginine'
        if re.search(r'GG[TCAG]', my_list[i]):
            my_list[i] = 'Glycine'
        if re.search(r'TA[AG]', my_list[i]):
            my_list[i] = 'Stop'
        if re.search(r'TGA', my_list[i]):
            my_list[i] = 'Stop'
        if re.search(r'N+', my_list[i]):
            my_list[i] = 'Unknown'
        i += 1

    return(my_list) #This returns a list of amino acids 

# Function that counts the amino acids in a list and returns a dictionary with the name of the amino acid as a key, and the count as an element
def amino_acid_counting(aa_list):
    counts = {}
    for i in aa_list:
        counts[i] = counts.get(i, 0) + 1
    return(counts)

# Function that determines a)the length, b) the GC count, and c) the GC content of a string and returns a dictionary with the operation as a key, and the number as an element (used for whole DNA and coding DNA)
def basic_stats(dna):
    stats = {}
    stats['length_all'] = len(dna)
    stats['length_non'] = (dna.count('G') + dna.count('C') + dna.count('A') + dna.count('T'))
    stats['GC_count'] = dna.count('G') + dna.count('C')
    
    if stats['length_all'] > 0:
        stats['GC_content_all'] = stats['GC_count'] / stats['length_all']
    else:
        stats['GC_content_all'] = 0
    
    if stats['length_non'] > 0:
        stats['GC_content_non'] = stats['GC_count'] / stats['length_non']
    else:
        stats['GC_content_non'] = 0

    return stats


########################################### This should be executed if the user provides the whole sequence ########################################### 
if file1:
    
    chromosomes_full = {} 
    string = ''
    i=0  # line counter
    with open (file1, 'r') as fh:  #This is opening my Test File, we need to adapt this later to read in the actual file
        for line in fh: #Go through the line file by file
            i+=1 
            print(f'line {i}')
            line = line.rstrip().upper()  #Get rid of \n at the end of the line
            if line.startswith('>'): #These should be used as our keys
                if re.search(r'^>([1-9]|1[0-9]|2[0-3]|[XY])', line):
                    match_1 = re.search(r'^>([1-9]|1[0-9]|2[0-3]|[XY])', line)
                    key = match_1.group(1)
                    chromosomes_full[key] = {} 
                    string = ''
                else:
                    key = 'scaffold'
            elif key != 'scaffold':
                string += line
        stats = (basic_stats(string))
        chromosomes_full[key] = stats

    df_full= pd.DataFrame.from_dict(chromosomes_full , orient='index')
    print(df_full)
    df_full.to_csv('output_whole.csv', index=False) 


########################################### This should be executed if the user provides the coding DNS ########################################### 
if file2:
############## Tim's Code for reformatting the input file

############## This should use the output file from Tim and run it 

    chromosomes = {} #This is our outer dictionary, the structure should be, chromosomes as keys, elements should be an inner dictionary
    chromosomes_1 = {} #This is for merging the stats and amino acid analyses later - this is the one for amino acids
    chromosomes_2 = {} #This is for merging the stats and amino acid analyses later - this is the one for stats


with open (file2, 'r') as fh:  #This is opening my Test File, we need to adapt this later to read in the actual file
    for line in fh: #Go through the line file by file
        line = line.rstrip().upper()  #Get rid of \n at the end of the line
        if line.startswith('>'): #These should be used as our keys
            line = line.strip('>')
            chromosomes[line] = {} #This allows us to add the line as key to dictionary, again if we end up with a differnt header we need to adapt this
            key = line
            string = ''
            counter = 1
        else:
           string += line
        amino_acid_count = (amino_acid_counting(amino_acid_translation(string)))
        stats = (basic_stats(string))
        chromosomes_1[key] = amino_acid_count
        chromosomes_2[key] = stats
    #print(chromosomes_1)
    #print(chromosomes_2)


    # Combine the two dictionaries into one which we can save and use for further analyses
    chromosomes = {k: {**v, **chromosomes_1[k]} for k, v in chromosomes_2.items()}
    df = pd.DataFrame.from_dict(chromosomes , orient='index')
    print(df)
    df.to_csv('output_coding_dns.csv', index=False) 
