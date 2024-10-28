#!/usr/bin/env python3
import re
import sys
import pickle

## Obejctive: General descriptives
fasta = open(sys.argv[1])

def GC_content(fasta):
    seq_dict = {}
    header = ''
    for line in fasta:
       line = line.rstrip().upper()
       if line.startswith(">"):        
            print(f"This is your sequence header: {line}")
            if re.search(r"^>([\dXY]\d?)\s.*", line):
                trunc_head = re.search(r"^>([\dXY]\d?)\s.*", line)
                header = trunc_head.group(1)
                print(f"This is your current header that should in the the dictionary: {header}")
                if header not in seq_dict: 
                   seq_dict[header] = {'length_all': 0, 'length_non': 0,'GC_count' : 0, 'GC_content_all': 0.0, 'GC_content_non': 0.0,} 
            else:
                header = 'scaffold'
                continue 
       else:
            ### if the header is not scaffold and is already in the dictionary:
            if header != 'scaffold' and header in seq_dict:
                 seqlength = len(line)
                 seq_dict[header]['length_all'] += seqlength
             ##count only ATCG bases (no Ns) for length_non():
                 countbases = 0
                 for base in line:
                    if base == 'G':
                         countbases += 1
                    if base == 'C':
                         countbases += 1
                    if base == 'A':
                        countbases += 1
                    if base == 'T':
                        countbases += 1
                    else:
                         countbases +=0
                 seq_dict[header]['length_non'] += countbases
                 count = 0
                 for nt in line:
                    if nt == 'G':
                         count += 1
                    if nt == 'C':
                         count += 1
                    if nt == 'A':
                        count += 0
                    if nt == 'T':
                        count += 0
                    else:
                         count +=0
                 seq_dict[header]['GC_count'] += count
                
                
                #calculate the GC content:
                #GC_content = count/seqlength
                 if seq_dict[header]['length_all'] > 0:
                ##create the GC content for all bases in the file(ATCGN):
                    seq_dict[header]['GC_content_all'] = seq_dict[header]['GC_count'] / seq_dict[header]['length_all']
                 if seq_dict[header]['length_non'] > 0:
                ##create the GC content in only ACTGs (no 'n's):
                    seq_dict[header]['GC_content_non'] = seq_dict[header]['GC_count'] / seq_dict[header]['length_non']

            
            else:
             header = 'scaffold'
             continue 

              
    return(seq_dict)
whole_gen_dict = GC_content(fasta)
print(whole_gen_dict)
pickle.dump(whole_gen_dict, open("save.p", "wb"))     

                       
                
           
             
