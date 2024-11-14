##
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statistics import mean  # python built-in module (no need to install)
!pip install biopython
from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices


import pandas as pd

columns = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

blast_df = pd.read_csv('blast_mouse2human.txt', sep='\t', names=columns)
blast_df.head()





"""#### Box-plot showing E-value distribution for each query sequence


"""

import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(12, 6))
sns.boxplot(x='qseqid', y='evalue', data=blast_df)
plt.xlabel('Query Sequence IDs', fontsize=12)
plt.ylabel('E-value', fontsize=12)
plt.title('E-value Distribution by Query Sequence', fontsize=14)
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()



"""#### top-hits"""

top_hitmoments = blast_df.groupby('qseqid')['evalue'].idxmin()
top_hit = blast_df.loc[top_hitmoments].reset_index(drop=True)

top_hit.head()



"""#### Load UniProt_Partial_MM.fasta and UniProt_HS.fasta with `SeqIO` module"""

from Bio import SeqIO

mouse_seqs = SeqIO.to_dict(SeqIO.parse('UniProt_Partial_MM.fasta', 'fasta'))
human_seqs = SeqIO.to_dict(SeqIO.parse('UniProt_HS.fasta', 'fasta'))

list(mouse_seqs.keys())[:5], list(human_seqs.keys())[:5] #just checking if it worked





"""#### Pairwise alignment with `pairwise2`"""

from Bio import pairwise2
from Bio.Align import substitution_matrices

matrix = substitution_matrices.load("BLOSUM62")
filtered_hits = top_hit[top_hit['length'] <= 200]
for index, row in filtered_hits.iterrows():
    query_id = row['qseqid']
    subject_id = row['sseqid']

    query_seq = str(mouse_seqs[query_id].seq)
    subject_seq = str(human_seqs[subject_id].seq)

    alignments = pairwise2.align.localds(query_seq, subject_seq, matrix, -20, -10)

print(f"Alignment for {query_id} vs {subject_id}:")
print(pairwise2.format_alignment(*alignments[0]))

"""<h3 style="background-color: #daebff; border-color: #bad5f6; border-left: 5px solid #bad5f6; padding: 1.5em; color: #6f89a9">
Next-Generation Sequencing Data Analysis
</h3>

#### Load FASTQ file (sample.fastq) and create dictionaries
"""

from Bio import SeqIO

label2seq = {}
label2qual = {}
for record in SeqIO.parse('sample.fastq', 'fastq'):
    label = record.id
    label2seq[label] = str(record.seq)
    label2qual[label] = record.letter_annotations["phred_quality"]

# Checking if its working
list(label2seq.items())[:5], list(label2qual.items())[:5]





"""#### GC content per sequence

import matplotlib.pyplot as plt
import seaborn as sns

# Step 1
def gc_content(sequence):
    g = sequence.count('G')
    c = sequence.count('C')
    a = sequence.count('A')
    t = sequence.count('T')
    # Calculate the GC content percentage
    return (g + c) / (a + t + g + c) * 100 if (a + t + g + c) > 0 else 0

# Step 2
gc_contents = [gc_content(seq) for seq in label2seq.values()]

# Step 3
plt.figure(figsize=(10, 6))
plt.hist(gc_contents, bins=30, edgecolor='black', alpha=0.7)
plt.xlabel('GC-content (%)', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.title('GC-content Distribution', fontsize=14)
plt.show()





"""#### (Per-sequence average quality


from statistics import mean
import seaborn as sns
import matplotlib.pyplot as plt

# Step 1
def parse_quality_string(quality_scores):
    return quality_scores

# Step 2
average_qualities = [mean(parse_quality_string(qual)) for qual in label2qual.values()]

# Step 3
plt.figure(figsize=(10, 6))
sns.kdeplot(average_qualities, fill=True)

plt.xlabel('Average Quality Score', fontsize=12)
plt.ylabel('Density', fontsize=12)
plt.title('Average Quality Score Distribution', fontsize=14)
plt.show()





"""<h3 style="background-color: #daebff; border-color: #bad5f6; border-left: 5px solid #bad5f6; padding: 1.5em; color: #6f89a9">



