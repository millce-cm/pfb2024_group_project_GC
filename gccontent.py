#!/usr/bin/env python3

import re

fasta_dict = {}  #is initialized as an empty dictionary to hold the headers and sequences
header = ''      #initializes header as empty string ensuring variable exists before using it in the loop 
seq = ''         #initializes seq as an empty string ensuring variable exists before using it in the loop


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
    
with open("test2dna.py", 'r') as file:
    for line in file:
        line = line.strip()  # Remove leading/trailing whitespace
        line = line.upper() #making all lowercase letter to uppercase
        if line.startswith('>'):  # If there's a header, make it a key in dict
            header = line  # Update the header variable
            fasta_dict[header] = ''  # Initialize the header in the dictionary
        else:
            seq = line  # Get the sequence line
            fasta_dict[header] += seq  # Append the sequence to the corresponding header

# Now calculate GC content for each sequence in the dictionary
for key in fasta_dict:
    sequence = fasta_dict[key]
    gc_percentage_ATCG, gc_percentage, g_count, c_count, a_count, t_count, total_count_ATCG, total_count  = gc_content(sequence) #gets leght,  GC content, and base counts
    print(f"Header: {key}\nSequence: {sequence}\nSequence length with n: {total_count}\nSequence lenght only ATCG: {total_count_ATCG}\nGC Count:{g_count+c_count}\nGCContent_n: {gc_percentage:.2f}%\nGCContent_ATCG: {gc_percentage_ATCG:.2f}%")

