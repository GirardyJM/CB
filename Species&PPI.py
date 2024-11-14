###
import os
import sys
import pandas as pd

sys.path.append("./") 
 


def reverse_complement(seq):

    rev_comp=list(seq[::-1])
    i=0
    while i<len(rev_comp):
        if rev_comp[i] == 'A':
            rev_comp[i] = 'T'
        elif rev_comp[i] == 'T':
            rev_comp[i] = 'A'
        elif rev_comp[i] == 'C':
            rev_comp[i] = 'G'
        elif rev_comp[i] == 'G':
            rev_comp[i] = 'C'
        i+=1

    return ''.join(rev_comp)

seq_test = 'ACGTATAGGCTGACACGTAGAGATGGATGACCATAG'


seq_rc = reverse_complement(seq_test)

print('The original strand: {}'.format(seq_test))
print('The reverse complementary strand: {}'.format(seq_rc))

"""### (b) Translation

Test your translation function below with the testing cases
"""

def translate(seq: str, codon_dict: dict):

    aa_seq = ''
    # ---- YOUR CODE STARTS HERE ----
    start=False

    if seq.find('ATG')>0:

       for i in range(seq.find('ATG'), len(seq)-2, 3):
          codon=seq[i:i+3]

          if codon == 'ATG':
              start=True

          if start:
              if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                  break

              if codon in codon_dict:
                  aa_seq += codon_dict[codon] #will give me the amino acid correspondent in the codon dict
      # ----- YOUR CODE ENDS HERE -----

    return aa_seq

CODON_DICT  = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
                'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
                'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*', 'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
                'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
                'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
                'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
                'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
                'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}


# Test your function with the following sequences
seq1 = 'ACGTATAGGCTGACACGTAGAGATGGATGACCATAG'
seq2 = 'CTGACGTATAGGCTGACACGTAGAGGATACCATAGT'

print('Translated protein sequence for seq1: {}'.format(translate(seq1, CODON_DICT)))
print('Translated protein sequence for seq2: {}'.format(translate(seq2, CODON_DICT)))

"""### (c) Global Alignment

Test your implementation of `aln_score_mat` function here.
"""

score_dict = {'match': 1, 'mismatch': -1, 'gap': -2}

def aln_score_mat(seq1, seq2, score_dict):
    
    # ---- Initialize ----
    n_row = len(seq1) + 1
    n_col = len(seq2) + 1
    aln_score = [[0 for col in range(n_col)] for row in range(n_row)]  # initialize score matrix with 0

    gap = score_dict['gap']

    for i in range(1, n_row):
        aln_score[i][0] = aln_score[i - 1][0] + gap

    for j in range(1, n_col):
        aln_score[0][j] = aln_score[0][j - 1] + gap

    for i in range(1, n_row):
        for j in range(1, n_col):
            if seq1[i - 1] == seq2[j - 1]:
                diagonal = aln_score[i - 1][j - 1] + score_dict['match']
            else:
                diagonal = aln_score[i - 1][j - 1] + score_dict['mismatch']

            from_top = aln_score[i - 1][j] + gap
            from_left = aln_score[i][j - 1] + gap

            aln_score[i][j] = max(diagonal, from_top, from_left)

    # ----- YOUR CODE ENDS HERE -----

    return aln_score



seq1 = 'GATTACA'
seq2 = 'GTCGACGCA'

aln_score_mat(seq1, seq2, score_dict)



### (a) Data Analysis
"""

# Load Data
# ---- WRITE YOUR CODE HERE ----
df = pd.read_csv('/biogrid_ppi.csv')
df.head(5)

"""#### ID for protein pair

a function that applies to a single row in the dataframe that returns a string for PPI ID (`prot1:prot2` or `prot2:prot1`)

"""

def join_ppi(row):
    """
    Args:
        row: a single row in the dataframe
    Returns:
        A string that joins two protein IDs with ':'
    """
    # ---- YOUR CODE STARTS HERE ----
    if row['taxidA'] == 9606:
      return ':'.join([row['prot1'], row['prot2']])
    elif row['taxidB'] == 9606:
      return ':'.join([row['prot2'], row['prot1']])
    else:
      return ':'.join([row['prot1'], row['prot2']])




"""2) apply this function to each row of the dataframe (HINT: you may find `pd.apply()` helpful)"""

# Apply the function to every row and assign the results to a new column called 'ppi'

df['ppi'] = df.apply(join_ppi, axis=1)
df.head(5)

"""#### Summarize Publication Information

Please find the top-5 publications (identified by `pub_id`) with most PPIs detected
"""


toppubid=df.groupby('pub_id').size().sort_values(ascending=False).head(5)
top5=toppubid.index.tolist()
print(top5)



"""### (b) Merge Dataframe"""

# Load Data
# ---- WRITE YOUR CODE HERE ----
info=pd.read_csv('/species_info.csv')
info.head(5)



"""#### Merge

Merge the species information (`info`) to the PPI table (`df`) by taxonomy ID (`taxid_x`). If a taxonomy ID cannot be found in `info`, please assign 'others' to its `family` and `species` columns.

Please create a new dataframe called `df_merge` for the resulted table. Please keep all records in the original PPI table.
"""

df_merge = pd.merge(df, info, how='left', left_on='taxid_x', right_on='taxid')
df_merge['family'].fillna('others', inplace=True)
df_merge['species'].fillna('others', inplace=True)



"""#### Save merged data to file"""

df_merge.to_csv('/merge_final.csv', index=False)

