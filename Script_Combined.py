#!/usr/bin/env python

########################################### Set up our main script ########################################### 
import os, sys, re, argparse
import pandas as pd
import pickle

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
    
    if stats['length_all'] > 0: #Note: this is probably unnecessary because length is always >) if we have an input, kept it in here so it does not break if we have an empty one
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

    def run_whole(): #Note, this whole section did not run before putting it in a function. Noone udnerstand why it runs now, but it does. 
        
        chromosomes_full = {} #This is our outer dictionary
        
        this_seq = ''
        i=0  # line counter
        with open (file1, 'r') as fh:  

            for i, line in enumerate(fh): #Go through the file line by line, not printing this anymore because it takes lot of time
                #print(f'{i}:{line}')
                line = line.rstrip()  #Get rid of \n at the end of the line

                if line.startswith('>'): #All lines that start with >
                    print('in header line')
                    if this_seq: #We are moving this up to ensure that our operations are done once on the whole concatenated string
                        print('we have a seq for statistics') 
                        stats = (basic_stats(this_seq)) #This calls our stats function from above
                        chromosomes_full[key] = stats #The output from our stats function is a dictionary
                        print('we have created our dictionary')
                        this_seq = '' #Setting sequence back to emoty for our next key
                    if not 'scaffold' in line and not line.startswith('>MT'): #We are in a header, that does not include scaffold or MT
                        key = line.split(' ')[0].lstrip('>') #We are creating out key from the header line, by extracting everything before the first  empty character and then stripping the >
                        print(f'working on {key}')
                    else:
                        key = ''#This is if we are in a header line, that is either a scaffold or starts with MT
                    

                elif key:  # We are in the sequence, adding all lines into our this_seq, which we will use for our statistics
                    this_seq += line.upper()
                
        if this_seq: #This is needed to also include the stats for the last key - sequence pair
            stats = (basic_stats(this_seq))
            chromosomes_full[key] = stats
        
        print(chromosomes_full) #We kept this print in here to make sure it ran

        pickle.dump(chromosomes_full, open('output_whole_genome.p', 'wb')) #Saves our dictionary to later inport it

        df_full= pd.DataFrame.from_dict(chromosomes_full , orient='index') #Converst our dictionar to a dataframe and saves the dataframe
        df_full.to_csv('output_whole_genome.csv', index=True) 

    if __name__ == "__main__":
        run_whole()


########################################### This should be executed if the user provides the coding DNS ########################################### 
if file2:

############## Tim's Code for reformatting the input file

############## This should use the output file from Tim and run it 

#Because this runs, we kept it outside a function. However, we could move this into a function as well. 

    chromosomes = {} #This is our finl outer dictionary. 
    chromosomes_1 = {} 
    chromosomes_2 = {}
    this_seq = ''

    with open (file2, 'r') as fh:  #This should be adapted once we incorporate Tim's code into this script
        for line in fh: #Go through the line file by file
            line = line.rstrip()  #Get rid of \n at the end of the line
            if line.startswith('>'): #Find lines that start with >
                if this_seq: #See above - we are doing this to efficiently calculate our stats once on the whole sequence for one chromosome
                    amino_acid_count = (amino_acid_counting(amino_acid_translation(this_seq))) #This is calling the amino acid translation and counting functions
                    stats = (basic_stats(this_seq)) #This is calling the stats function from above. 
                    chromosomes_1[key] = amino_acid_count #This is creating an outer dictionary that only has the results from the aminoacid count as inner dictionary 
                    chromosomes_2[key] = stats #This is creating an outer dictionary that only has the results from the basic count as inner dictionary 
                    this_seq = '' #This is resetting our sequence for the next key
                if not 'scaffold' in line and not line.startswith('>MT'): #This is creating our keys, should not be needed, but double cheks Tim's script
                    key = line.split(' ')[0].lstrip('>')
                    print(f'working on {key}')
                else:
                    key = ''
                
            elif key:  # We are in our sequence and concatenate all lines to the sequence
                this_seq += line.upper()
    if this_seq: #Again, this is needed for the last key
        amino_acid_count = (amino_acid_counting(amino_acid_translation(this_seq)))
        stats = (basic_stats(this_seq))
        chromosomes_1[key] = amino_acid_count
        chromosomes_2[key] = stats

    #Here we combine our two helper dictionaries (which have the same keys) into our final dictionary
    chromosomes = {k: {**v, **chromosomes_1[k]} for k, v in chromosomes_2.items()}
    
    print(chromosomes)
    
    pickle.dump(chromosomes, open('output_coding_genome.p', 'wb'))
    df = pd.DataFrame.from_dict(chromosomes , orient='index')
    df.to_csv('output_coding_genome.csv', index=True) 


