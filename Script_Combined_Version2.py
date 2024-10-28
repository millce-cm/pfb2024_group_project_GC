#!/usr/bin/env python

########################################### Set up of the script ########################################### 
import os, sys, re, argparse
import pandas as pd
import pickle

parser = argparse.ArgumentParser(description= 'This is a script that has two optional arguments: file_1 and file_2. If you provide a whole genome (file_1), you will receive gene length, GC count, and GC content for each chromosome as output. If you provide the coding sequence, you will receive length, GC count, GC content, amino acid count, and amino acid content for each chromosome as output.')
parser.add_argument('-f1', '--file_1', help ='optional: Path to FASTA file for the whole genome')
parser.add_argument('-f2', '--file_2', help ='optional: Path to FASTA file for the coding sequence')
args = parser.parse_args()
file1 = args.file_1
file2 = args.file_2

########################################### Functions ########################################### 

# Function that determines a) the DNA length (with/without 'n's), b) the GC count, and c) the GC content (with/without taking 'n's into account) of a string and returns a dictionary with the name of the operation as a key and the number as an element (used for whole DNA and coding DNA).

def basic_stats(dna):
    stats = {} #Final dictionary to save the results. 
    stats['length_all'] = len(dna) #DNA length
    stats['length_non'] = (dna.count('G') + dna.count('C') + dna.count('A') + dna.count('T')) #DNA length without 'n's.
    stats['GC_count'] = dna.count('G') + dna.count('C') #GC count. 
    
    if stats['length_all'] > 0: 
        stats['GC_content_all'] = stats['GC_count'] / stats['length_all'] #GC content (displayed as 0.x; with taking 'n's into account)
    else:
        stats['GC_content_all'] = 0
    
    if stats['length_non'] > 0:
        stats['GC_content_non'] = stats['GC_count'] / stats['length_non'] #GC content (displayed as 0.xx; without taking 'n's into account)
    else:
        stats['GC_content_non'] = 0

    return stats

# Function that a) splits a DNA string into triplets and b) translates the triplets into amino acids.

def amino_acid_translation (dna):  

    dna_upper = dna.upper() # Converts everything into upper letters. This is how we decided to deal with soft maked regions. 
    my_list = re.findall(r'(\S{3})', dna_upper) # Splits the string into a list, each element being a codon.
    i = 0 #This is needed for iterating through my_list

    while i < len(my_list): # The translation is based on the "Standard DNA codon Table" from Wikipedia (https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables). 
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
        if re.search(r'N+', my_list[i]): #We included this for all instances where we find a N in a triplet. 
            my_list[i] = 'Unknown'
        i += 1

    return(my_list) #This returns a list of amino acids. 

# Function that a) counts the amino acids in a list and b) calculates the content of each amino acid, returning a dictionary with the name of the amino acid as a key and the count and content as elements. 

def amino_acid_counting(aa_list):
    counts = {}
    for i in aa_list:
        if i not in counts:
            counts[i] = 0
        counts[i] += 1

    total = sum(counts.values())
    
    output = {}
    for i in counts:
        output[f'{i}_count'] = counts[i]
        output[f'{i}_percentage'] = (counts[i] / total) 
    return(output)


########################################### This should be executed if the user provides no input files ########################################### 

if file1 == None and file2 == None:
    print('You did not provide an input file. Usage: Please provide either file_1 (whole genome), file_2 (coding genome), or both file_1 and file_2 (whole genome and coding genome).')

########################################### This should be executed if the user provides the whole sequence (file_1) ########################################### 

if file1:

    print('You provided a whole genome. The output will include gene length, GC count, and GC content for each chromosome. The output will be in the form of a double dictionary. It will be saved as a pickle file and a CSV file for further analyses.')

    def run_whole(): #Note, this whole section did not run before putting it in a function. Noone understand why it runs now, but it does. 
        
        chromosomes_full = {} #This is the outer dictionary. The keys should be the chromosomes we are looping through. 
        
        this_seq = ''
        i=0  # line counter

        with open (file1, 'r') as fh:  

            for i, line in enumerate(fh): #Go through the file line by line, not printing this anymore because it takes a lot of time but can do that for debugging. 
                #print(f'{i}:{line}')
                line = line.rstrip()  #Get rid of \n at the end of the line. 

                if line.startswith('>'): #Select all lines that start with '>' as headers. 
                    print('in header line') #We kept this in here, because it makes it easier to follow the progress of the script. 
                    if this_seq: #We moved this up to ensure that our operations are done once on the whole concatenated this_seq string. 
                        print('we have a seq for statistics') #We kept this in here, because it makes it easier to follow the progress of the script.
                        stats = (basic_stats(this_seq)) #This calls our stats function from above. 
                        chromosomes_full[key] = stats #The output from our stats function is a dictionary. 
                        print('we have created our dictionary') #We kept this in here, because it makes it easier to follow the progress of the script.
                        this_seq = '' #Setting sequence back to empty for our next key
                    if not 'scaffold' in line and not line.startswith('>MT'): #We are in a header, that does not include scaffold or MT. This is relevant for the human. 
                    #if not 'scaffold' in line and not line.startswith('>MT') and not line.startswith('>VNFC'): #We used this for the whale. VNFC stands for: Variable Number of FLN (Fibrinogen-like) Copies, which might be relevant for adaption to living in water. 
                        key = line.split(' ')[0].lstrip('>') #We are creating our key from the header line, by extracting everything before the first empty character and then stripping the '>'.
                        print(f'working on {key}') #We kept this in here, because it makes it easier to follow the progress of the script.
                    else:
                        key = ''#This is if we are in a header line, that is either a scaffold or starts with MT. (Or anything else specificed by the user). 
                    

                elif key:  # We are in the sequence, adding all lines into our this_seq, which we will use for our statistics.
                    this_seq += line.upper()
                
        if this_seq: #This is needed to also include the stats for the last key - sequence pair.
            stats = (basic_stats(this_seq))
            chromosomes_full[key] = stats
        
        print(chromosomes_full) #We kept this in here, because it makes it easier to follow the progress of the script.

        #Please specify the organism if you change it in your output names. If you do not change the output name, it will overwrite the files created previously. 
        pickle.dump(chromosomes_full, open('output_whole_genome_human.p', 'wb')) #Saves our dictionary to later inport it as a pickle. 

        df_full= pd.DataFrame.from_dict(chromosomes_full , orient='index') #Converst our dictionary to a dataframe and saves the dataframe. 
        df_full.to_csv('output_whole_genome_human.csv', index=True) 

    if __name__ == "__main__":
        run_whole()


########################################### This should be executed if the user provides the coding DNS (file_2) ########################################### 
if file2:
    
    print('You provided a coding genome. The script will create a new file that combines the coding sequence of each gene for a chromosome. The output files will include genome length, GC count, and GC content, as well as the count and content of all amino acids per chromosome.')

############## Cleaning the code
    print('Please wait—the code is being cleaned up, and a file for your further analyses is being created.')

    # Please note that you should change the name of the output file if you do not want to overwrite the previous created fa file. However, it is important that your input for the next section matches thes output name provided here. 
    
    #This initalizes the dictionaries used in the code below. 
    chrom_dict = {}
    gene_id_dict = {}

    with open(file2, 'r') as fasta_in_fh, open('coding_sequence_per_chromosome_human.fa', 'w') as fasta_out_fh:
        
        current_gene_id = None  # Store the current gene ID
        header = ''
        
        for line in fasta_in_fh:
            
            if line.startswith('>'):

                if 'protein_coding' in line and 'scaffold' not in line: #Check for 'protein_coding' in the header line and exclude headers that have a 'scaffold' in it.
                    chrom_match = re.search(r'^>.*:([\dXY]\d?):.*gene:(ENSG\S*)\s', line) #Extract chromosome and gene ID from the header.This is the code we used for the human genome.
                    #chrom_match = re.search(r'^>.*:([\dXY]\d?):.*gene:(ENS\S*)\s', line) #Extract chromosome and gene ID from the header.This is the code we used for the mouse and whale genome. Please check if you use the genome of a differnt species. 
                    
                    if chrom_match:
                        chrom_id = chrom_match.group(1)
                        current_gene_id = chrom_match.group(2)

                        # Initialize if this gene ID hasn't been seen before
                        if current_gene_id not in gene_id_dict:
                            gene_id_dict[current_gene_id] = True  # Mark it as seen
                            chrom_dict.setdefault(chrom_id, '')  # Initialize if not present
                        else:
                            current_gene_id = None  # Reset to avoid duplicates
                else:
                    # Skip lines that don't contain "protein_coding"
                    current_gene_id = None
                    
            else:
                # Only append if we have a valid gene ID and it's a multiple of 3
                if current_gene_id is not None and len(line.strip()) % 3 == 0:  # Check if line length is a multiple of 3
                    chrom_dict[chrom_id] += line  # Append sequence lines

        # Write output
        for key in chrom_dict:
            fasta_out_fh.write(f">{key}\n{chrom_dict[key]}")

############## This runs with the output file from above

    print('Your file has been created. The statistical functions are being ran next.')

    # Please note that if you changed the name of the output file above, you need to change the input file here. 

    def run_sequence():
    
        chromosomes = {} #This is our final outer dictionary. 
        chromosomes_1 = {} #This is the dictionary for the amino acid results. If wanted, one could only use one of the two dictionaries. 
        chromosomes_2 = {} #This is the dictionary for the basi stats results. If wanted, one could only use one of the two dictionaries. 
        this_seq = ''

        with open ('coding_sequence_per_chromosome_human.fa', 'r') as fh:  

            for line in fh: #Go through the file line by line. 
                line = line.rstrip()  #Get rid of \n at the end of the line

                if line.startswith('>'): #Find lines that start with '>'.
                    print('in header line') #We kept this in here, because it makes it easier to follow the progress of the script. 

                    if this_seq: #See above - we are doing this to efficiently calculate our stats once on the whole sequence for one chromosome
                        print('we have a seq for statistics') #We kept this in here, because it makes it easier to follow the progress of the script.
                        amino_acid_count = (amino_acid_counting(amino_acid_translation(this_seq))) #This is calling the amino acid translation and counting functions. 
                        stats = (basic_stats(this_seq)) #This is calling the stats function. 
                        chromosomes_1[key] = amino_acid_count #This is creating an outer dictionary that only has the results from the aminoacid count as inner dictionary. 
                        chromosomes_2[key] = stats #This is creating an outer dictionary that only has the results from the basic count as inner dictionary.
                        this_seq = '' #This is resetting our sequence for the next key
                    
                    if not 'scaffold' in line  and not line.startswith('>MT'): #The way the first part of the code works, this should not be necessary. However, depending on the input files, this can be used for further excluding parts of the genome. 
                        key = line.split(' ')[0].lstrip('>')
                        print(f'working on {key}') #We kept this in here, because it makes it easier to follow the progress of the script.
                    else:
                        key = ''
                    
                elif key:  #We are in our sequence and concatenate all lines to this_seq. 
                    this_seq += line.upper()

        if this_seq: #Again, this is needed for the last key
            amino_acid_count = (amino_acid_counting(amino_acid_translation(this_seq)))
            stats = (basic_stats(this_seq))
            chromosomes_1[key] = amino_acid_count
            chromosomes_2[key] = stats

        #Here we combine our two helper dictionaries (which have the same keys) into our final dictionary. 
        chromosomes = {k: {**v, **chromosomes_1[k]} for k, v in chromosomes_2.items()}
        
        print(chromosomes) #We kept this in here, because it makes it easier to follow the progress of the script.
        
        #Please specify the organism if you change it in your output names. If you do not change the output name, it will overwrite the files created previously. 
        pickle.dump(chromosomes, open('output_coding_genome_human.p', 'wb')) #Saves our dictionary to later inport it as a pickle. 
        
        df = pd.DataFrame.from_dict(chromosomes , orient='index') #Converst our dictionary to a dataframe and saves the dataframe. 
        df.to_csv('output_coding_genome_human.csv', index=True) 
    
    if __name__ == "__main__":
        run_sequence()


