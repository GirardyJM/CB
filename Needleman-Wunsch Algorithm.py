##
!pip install Bio
import os
import sys
import copy
import traceback

import re

import pandas as pd
import numpy as np

from Bio import SeqIO

#Needleman-Wunsch Algorithm
</h3>

##Construct similarity matrix for a pair of sequences
"""

def get_sim_mat(seq1, seq2):
    """
    Construct a similarity matrix for two sequences
    Args:
        seq1 [str]: sequence 1 (vertical)
        seq2 [str]: sequence 2 (horizontal)
    Returns:
        sim_mat [2D-list or 2D-array]: similarity matrix where 1=match, 0=mismatch
    """
    sim_mat=np.zeros((len(seq1),len(seq2)))
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            if seq1[i]==seq2[j]:
                sim_mat[i][j]=1
            else:
              sim_mat[i][j]=0
    return sim_mat


# Test case
seq1 = "AYCYNRCKCRBP"
seq2 = "ABCNYRQCLCRPM"

sim_mat_test = get_sim_mat(seq1, seq2)
sim_mat_test



"""### Update cell value"""

def update_cell_value(current_matrix, i, j):
	
    currentvalue= current_matrix[i,j]
    if i+1 <len(current_matrix) and j+1 <len(current_matrix[0]):
        diagonal=current_matrix[i+1,j+1]
    else:
        diagonal=0
    if i+1 <len(current_matrix) and j+2 < len(current_matrix[0]):
        row_gap=max(current_matrix[i+1,j+2:])
    else:
        row_gap=0
    if i + 2 < len(current_matrix) and j+1 < len(current_matrix[0]):
        col_gap = max(current_matrix[x, j + 1] for x in range(i + 2, len(current_matrix)))
    else:
        col_gap = 0

    return max(diagonal,row_gap,col_gap) + currentvalue






update_cell_value(sim_mat_test, 9, 10)




def get_sum_mat(sim_mat):

    for i in range(len(sim_mat)-1,-1,-1):
        for j in range(len(sim_mat[0])-1,-1,-1):
            sim_mat[i,j]=update_cell_value(sim_mat,i,j)
    return sim_mat


sum_mat_test = get_sum_mat(sim_mat_test)
sum_mat_test





"""##Find cell with the maximum value (to start tracing back)"""

def find_max_cell(sum_mat):

    max_value=0
    max_cell=(0,0)
    for i in range(len(sum_mat)-1,-1,-1):
                for j in range(len(sum_mat[0])-1,-1,-1):
                      if max_value < sum_mat[i,j]:
                        max_value=sum_mat[i,j]
                        max_cell=(i,j)
    return max_cell



max_cell = find_max_cell(sum_mat_test)
max_cell



def find_next_maxcell(sum_matrix, i, j):

    max_value=0
    max_cell=(i,j)
    if i + 1 < len(sum_matrix):
        for y in range(j + 1, len(sum_matrix[0])):
            if sum_matrix[i + 1, y] > max_value:
                max_value = sum_matrix[i + 1, y]
                max_cell = (i + 1, y)
            elif sum_matrix[i + 1, y] == max_value and max_cell == (None, None):
                max_cell = (i + 1, y)

    if j + 1 < len(sum_matrix[0]):
        for x in range(i + 1, len(sum_matrix)):
            if sum_matrix[x, j + 1] > max_value:
                max_value = sum_matrix[x, j + 1]
                max_cell = (x, j + 1)
            elif sum_matrix[x, j + 1] == max_value and max_cell == (None, None):
                max_cell = (x, j + 1)

    if i == len(sum_matrix) - 1 or j == len(sum_matrix[0]) - 1 or max_cell==(i,j):
        return (None, None)

    return max_cell



next_max_cell = find_next_maxcell(sum_mat_test, 0, 0)
next_max_cell





"""###Find the trace-back path"""

def trace_back(sum_mat):

  max_cell=find_max_cell(sum_mat)
  trace_list=[]
  trace_list.append(max_cell)
  while max_cell[0]!=None and max_cell[1]!=None:
       max_cell=find_next_maxcell(sum_mat,max_cell[0],max_cell[1])
       if max_cell[0]!=None and max_cell[1]!=None:
          trace_list.append(max_cell)

  return trace_list



trace_list = trace_back(sum_mat_test)
trace_list





"""### Get alignment"""

def get_alignment(seq1, seq2):

      
  similarity_mat=get_sim_mat(seq1,seq2)

     
  sim_sum_mat=get_sum_mat(similarity_mat)

  trace_list_align=trace_back(sim_sum_mat)
  print(trace_list_align)

      
  seq1_aligned = ""
  seq2_aligned = ""
  prev=(0,0)

  for (i, j) in trace_list_align:
        if i == 0 and j == 0:
            seq1_aligned += seq1[i]
            seq2_aligned += seq2[j]
            prev=(i,j)
        elif i > 0 and  j > 0 and seq1[i]==seq2[j]:
          for k in range(prev[0]+1,i):
            seq1_aligned += seq1[k]
            seq2_aligned += "-"
          seq1_aligned += seq1[i]
          seq2_aligned += seq2[j]
          prev=(i,j)

  if prev[0] < len(seq1) - 1:
        for k in range(prev[0] + 1, len(seq1)):
            seq1_aligned += seq1[k]
            seq2_aligned += "-"


  print("Aligned seq1:", seq1_aligned)
  print("Aligned seq2:", seq2_aligned)


seq1 = "ATTCGA"
seq2 = "ACG"

alignment = get_alignment(seq1, seq2)
alignment


def parse_sam_row(row, headers):
    """
    Parse a single row from alignment section in SAM file to a dictionary
    """
    fields = row.strip().split('\t')
    sam_dict = {}
    for i in range(len(headers)):
        sam_dict[headers[i]] = fields[i]
    return sam_dict

with open('test_sam_row.txt', 'r') as f:
    test_row = f.readlines()

headers = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]
parse_sam_row(test_row[0].strip(), headers)



"""### Parse SAM file"""

import pandas as pd

def parse_sam_file(fname, headers, int_cols):
    with open(fname, 'r') as file:
        lines = file.readlines()

    data = [parse_sam_row(line, headers) for line in lines if not line.startswith('@')]
    sam_df = pd.DataFrame.from_dict(data)
    for col in int_cols:
        sam_df[col] = sam_df[col].astype(int)

    return sam_df


headers = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]
int_cols = ["FLAG", "POS", "MAPQ", "TLEN"]
parsed_sam = parse_sam_file("Clone_Seq_Aligned.sam", headers, int_cols)
parsed_sam.head()





"""### (c) Decode FLAG"""

def decode_FLAG(flag):
   
    bin_str = bin(flag)[2:].zfill(12)
    return bin_str

decode_FLAG(99)



"""

As the binary string 000001100011, bit 1 is set to 1, indicating that the read is part of a pair, and Bit 2 is also set, showing that the read is mapped in a proper pair, meaning it aligns correctly with its mate according to the sequencing protocol. Bit 3 is not set, which means the read itself is mapped, and Bit 4 is also unset, indicating that the mate is mapped as well. Bit 5 is not set, signifying that the read is on the forward strand, while Bit 6 is set, showing that the mate is on the reverse strand. Lastly, Bit 7 is set, indicating that this read is the first in the pair. The remaining bits are not set.So the read is paired. It is mapped in a proper pair.The mate is on the reverse strand and the read is the first in the pair.

"""
