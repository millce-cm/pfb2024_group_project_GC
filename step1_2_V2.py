#!/usr/bin/env python3
import re
import sys

## Obejctive: General descriptives
# Create a dictionary within a dictionary:
# Sequence/chromosome name as the Main Key:
    #For each Main Key, create a dictionary consisting of:
    #   {'MainKey1': {'length'=123, 'GC_count' = 123, 'GC_content' =}, 'MainKey2'... }

fasta = open(sys.argv[1])

def GC_content(fasta):
    seq_dict = {}
    header = ''
    for line in fasta:
       line = line.rstrip().upper()
       
       if line.startswith(">"):        
            print(f"This is your sequence header: {line}")
        #   trunc_head = re.search(r">([1-9]|1[0-9]|2[0-3]|[XY])", line)
            trunc_head = re.search(r"^>.*:([\dXY]\d?):.*", line)
            ### up to this point we've search through a line starting with ">"
            ### and extracted the regex defined above. 
            ### if there the regex expression was identified, then continue with below:
            if trunc_head != None:
                 if trunc_head.group(1) not in seq_dict: 
                   ###up to now, we have extracted the value of the regex (eg., '1', 'X'...)
                   ### if that value is not already in the dictonary, then we will do something
                   # starting with print...
                   print(f"This is your truncated header: {(trunc_head.group(1))}")
                   seq_dict[trunc_head.group(1)] = {'length_all': 0, 'length_non': 0,'GC_count' : 0, 'GC_content_all': 0.0, 'GC_content_non': 0.0,} 
                 ###if the regex expression IS already in there, then:
                 else:
                      ##why is this here?
                      header = trunc_head.group(1)
            ###otherwise it is a scaffold:
            else:
                header = 'scaffold'
                continue 
       else:
            ### if the header is not scaffold and is already in the dictionary:
            if header != 'scaffold' and header in seq_dict:
                 print(f"This is your sequence that should be all CAPITALS: {line}")
                 ##then calculate the seq length
                 seqlength = len(line)
                 print(f"This is your entire sequence length (ATCGN): {seqlength}")
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
                 print(f"This is your shortened sequence length (ATCG): {countbases}")
                 seq_dict[header]['length_non'] += countbases
           
               ##calcualte the GC count:
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
                 print(count)
                 seq_dict[header]['GC_count'] += count
                
                
                #calculate the GC content:
                #GC_content = count/seqlength
                 if seq_dict[header]['length_all'] > 0:
                ##create the GC content for all bases in the file(ATCGN):
                    seq_dict[header]['GC_content_all'] = seq_dict[header]['GC_count'] / seq_dict[header]['length_all']
                 if seq_dict[header]['length_non'] > 0:
                ##create the GC content in only ACTGs (no 'n's):
                    seq_dict[header]['GC_content_non'] = seq_dict[header]['GC_count'] / seq_dict[header]['length_non']
                #print(GC_content)
                 print(seq_dict[header]['GC_content_all'])
                 print(seq_dict[header]['GC_content_non'])
            
            else:
             header = 'scaffold'
             continue 

              
    return(seq_dict)
result = GC_content(fasta)
print(result)

#               
#             
#                  
#                   
#                 
#                     
#                  
#                       
#                     
#           
#               
#              

                       
                
           
             