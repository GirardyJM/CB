from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices
# Load fasta file with general file operations
with open('./data/seq.fa', 'r') as f:
    fasta_lines = f.readlines()

fasta_lines

from Bio import SeqIO

for record in SeqIO.parse("./data/seq.fa", "fasta"):
    print("id:", record.id)
    print("seq:", record.seq)

type(record)



record_dict = SeqIO.to_dict(SeqIO.parse("./data/seq.fa", "fasta"))
record_dict

record_dict.keys()

record_dict['Query1']






substitution_matrices.load()

substitution_matrices.load("BLOSUM62")


# NW Result Check
for a in pairwise2.align.globalms("GATTACA", "GTCGACGCA", match=1, mismatch=-1, open=-2, extend=-2):
    print(format_alignment(*a))

for a in pairwise2.align.localms("GATTACA", "GTCGACGCA", match=1, mismatch=-1, open=-2, extend=-2):
    print(format_alignment(*a))

"""### Global Alignment

"""

# Input Setup
s1 = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAM"
s2 = "MEERLVSPSAQMPLSQETFSKLLPENNVLSPLPLVSQAM"
matrix = substitution_matrices.load("BLOSUM62")
gap_open = -10
gap_extend = -0.5

type(matrix)

align_g = pairwise2.align.globalds(s1, s2, match_dict=matrix, open=gap_open, extend=gap_extend)

align_g



seqA, seqB, score, begin, end = align_g[0]
formated_alignment =  pairwise2.format_alignment(seqA, seqB, score, begin, end)
print(formated_alignment)



# This yields the same result as above
formated_alignment =  pairwise2.format_alignment(*align_g[0])
print(formated_alignment)

"""### Local Alignment

"""

# Setup
s3 = "LSPADKTNVKAA"
s4 = "PEEKSAV"
matrix = substitution_matrices.load("BLOSUM62")
gap_open = -10
gap_extend = -0.5

align_l = pairwise2.align.localds(s3, s4, match_dict=matrix, open=gap_open, extend=gap_extend)
align_l


format_alignment =  pairwise2.format_alignment(*align_l[0])
print(format_alignment)


format_alignment = pairwise2.format_alignment(*align_l[0], full_sequences=True)
print(format_alignment)

"""# BLAST

"""
import pandas as pd
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

fasta_string = open("./data/seq.fa").read()


result = NCBIWWW.qblast("blastn", "nt", fasta_string, expect=1e-50)

result



blast_record = NCBIXML.read(result)


blast_record



blast_record.alignments[0]


for i in range(5):
    blast_aln = blast_record.alignments[i]
    for hsp in blast_aln.hsps:
        print("**** Alignment {} ****".format(i+1))
        print("sequence:", blast_aln.title)
        print("length:", blast_aln.length)
        print("e value:", hsp.expect)
        print(hsp.query[0:50] + "...")
        print(hsp.match[0:50] + "...")
        print(hsp.sbjct[0:50] + "...")

"""### Load BLAST result

"""



# Here I've just added a parser to read in the standard BLAST output
def blast2df(file, header=("qseqid", "sseqid", "pident", "length", "mismatch",
                           "gapopen", "qstart", "qend", "sstart", "send",
                           "evalue", "bitscore")):

    return pd.read_csv(file, sep="\t", header=None, names=header)

df = blast2df("data/blastn_mouse2human.txt")
df

# So we can see that we have many potential hits
# Let's do a little bit of re-formatting and filtering
df = df.sort_values(["qseqid", "evalue"], ascending=[True, True])
df

# Now let's only consider the top hit for each entry
df = df.drop_duplicates("qseqid", keep="first")
df
