#!/usr/bin/env python3
import re
import pickle

fasta_dict = {
'chromosome':[],
'sequence':[],
'gc_percentage_ATCG':[],
'gc_percentage':[],
'g_count':[],
'c_count':[],
'a_count':[],
't_count':[],
'total_count_ATCG':[],
'total_count':[] 
}  #is initialized as an empty dictionary to hold the headers and sequences
list_seq=[]
header = ''      #initializes header as empty string ensuring variable exists before using it in the loop 
seq = ''         #initializes seq as an empty string ensuring variable exists before using it in the loop
str_seq =''
def gc_content(sequence):  #this line defines a function names gc_content that takes one argument, sequence. 
    """Calculate the GC content of a given DNA sequence."""
    g_count = sequence.count('G')  #line counts how many times nucleotide 'G' appears in sequence
    c_count = sequence.count('C')  #line counts how many times nucleotide 'C' appears in sequence
    a_count = sequence.count('A')  #line counts how many times nucleotide 'A' appears in sequence
    t_count = sequence.count('T')  #line counts how many times nucleotide 'T' appears in sequence
   
    total_count_ATCG = g_count + c_count + a_count + t_count  # Total count of A, T, C, G

    total_count = len(sequence)    #calculates the total number of nucleotides in the sequence

    if total_count == 0:  # Avoid division by zero (this checks if the sequence is empty)
        return 0.0, 0, 0, 0, 0, 0   # if empty, function return 0GC content and counts if no sequence
    else:
        return (g_count + c_count) / (total_count_ATCG) * 100, (g_count + c_count) / (total_count) * 100, g_count, c_count, a_count, t_count, total_count_ATCG, total_count  # this is GC content over ACTG alone and total seq lenght, and invividual counts of A,C,T,G 
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# full genome file : Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.tim

with open("Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.tim", 'r') as file:   #opens file for reading
    for line in file:
        line = line.strip()         # Remove leading/trailing whitespace
        if line.startswith('>'):    #if line starts with > its a header
            if line.startswith('>K') or line.startswith('>G') or line.startswith('>MT'):  # If header has >K,G,MT its marked as bad
                header = 'bad'
            else:                    #otherwise it uses a regex to extract chromosome number and append to dictionary[chromosome]
                match = re.search(r'^>([\dXY]\d?)', line)
                header = match.group(1)  # Update the header variable
                fasta_dict['chromosome'].append(header)
                if len(list_seq) != 0:          # if sequences are collected in list_seq, it joins them into single string
                    str_seq = str_seq.join(list_seq)
                    fasta_dict['sequence'].append(str_seq)      #appends single string into dictionary[sequence]
                    list_seq = []       #reset the list for the next sequence
                    #Store GC content and counts into dictionary using function gc_content (processing sequence lines)
                    gc_percentage_ATCG, gc_percentage, g_count, c_count, a_count, t_count, total_count_ATCG, total_count  = gc_content(str_seq)
                    fasta_dict['gc_percentage_ATCG'].append(gc_percentage_ATCG)
                    fasta_dict['gc_percentage'].append(gc_percentage)
                    fasta_dict['g_count'].append(g_count)
                    fasta_dict['c_count'].append(c_count)
                    fasta_dict['a_count'].append(a_count)
                    fasta_dict['t_count'].append(t_count)
                    fasta_dict['total_count_ATCG'].append(total_count_ATCG)
                    fasta_dict['total_count'].append(total_count)
                    str_seq =''
        elif header != 'bad':           #if header is not marked as bad, it processes line as sequence...
                line = line.upper()     #making all lowercase letter to uppercase
                seq = line              # Assigns the sequence line
                #print(seq)
                list_seq.append(seq)   #add the sequence to the list
    str_seq = str_seq.join(list_seq)            #after the loop, it joins any remaining sequences and adds to dictionary
    fasta_dict['sequence'].append(str_seq)
    #Store GC content and counts into dictionary using function gc_content (processing sequence lines) - for the final sequence
    gc_percentage_ATCG, gc_percentage, g_count, c_count, a_count, t_count, total_count_ATCG, total_count  = gc_content(str_seq)
    fasta_dict['gc_percentage_ATCG'].append(gc_percentage_ATCG)
    fasta_dict['gc_percentage'].append(gc_percentage)
    fasta_dict['g_count'].append(g_count)
    fasta_dict['c_count'].append(c_count)
    fasta_dict['a_count'].append(a_count)
    fasta_dict['t_count'].append(t_count)
    fasta_dict['total_count_ATCG'].append(total_count_ATCG)
    fasta_dict['total_count'].append(total_count)

fasta_dict['sequence']=[]   #resetting sequence list in the dictionary
#line 89-92 is a loop that adds the number 2 to fictionary[seq] for each count in g_count, serving as a placeholder to maintain dictionary structure.
i = 0
while i <len(fasta_dict['g_count']):
    fasta_dict['sequence'].append(2)
    i+= 1
print(len(fasta_dict['c_count']), fasta_dict['sequence'])
with open('fasta_dict.p', 'wb') as fp:          #this block opens a new file in binary write mode and saves output using pickle
     pickle.dump(fasta_dict, fp)
     print('dictionary saved successfully to file')     #import this into jupitorlab to make dataframe and plots

#this line prints a summary of collected data including headers, GC content and seq lenghts
print(f"Header: {fasta_dict['chromosome']}\nGC % ATCG: {fasta_dict['gc_percentage_ATCG']}%\nGC % with n: {fasta_dict['gc_percentage']}%\nSequence lenght only ATCG: {fasta_dict['total_count_ATCG']}\nSequence lenght with n: {fasta_dict['total_count']}\n")
