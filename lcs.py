#design a DNA sequence to fit the input network's topology, geometry, and connectivity:
import os
import re
import sys
import random

################################################################################
#FUNCTION: find all common substrings >=3 residues:

def acs(seq1,seq2):
    min_match = 4
    #initialize 2D list of seq1xseq2:
    seq_array = [[0]*len(seq1) for i in range(0,len(seq2))]
    for i in range(0,len(seq1)):
        for j in range(0,len(seq2)):
            if seq1[i]==seq2[j]:
                seq_array[j][i] = 1

    #for line in seq_array:
    #    print(line)

    #find all 1-diagonals
    diag_startend = []
    for i in range(0,len(seq1)):
        for j in range(0,len(seq2)):
            if seq_array[j][i] == 1:
                for k in range(0,len(seq1)+len(seq2)):
                    if j+k >= len(seq2) or i+k >= len(seq1):
                        break
                    if seq_array[j+k][i+k] != 1:
                        break
                    seq_array[j+k][i+k] = 2
                if k>=min_match:
                    diag_startend.append([[j,j+k],[i,i+k]])
              
    return diag_startend

################################################################################
#FUNCTION: find the reverse complement of a given DNA sequence:
def revcomp(sequence):
    sequence = sequence.replace('A','a').replace('T','t').replace('G','g').replace('C','c')
    sequence = sequence.replace('a','T').replace('t','A').replace('g','C').replace('c','G')
    sequence = sequence[::-1] #reverse
    return sequence

#FUNCTION: find the longest common DNA sequence with a complement present in the system:
def acs_DNA(seq1, seq2):

    #find the reverse complement of seq2:
    seq2 = revcomp(seq2)
    return acs(seq1,seq2)

################################################################################
#print(acs('abcdeXfgh','abcdefgh'))
import random

S1='AGGCAGTTGAGACGAACATTCCTAAGTCTGAAATTTATCACCCGCCATAGTAGACGTATCACC'
S2='CTTGCTACACGATTCAGACTTAGGAATGTTCGACCATGCGAGGGTCCAATACCGACGATTACAG'
S2=''.join(random.sample(S2,len(S2)))
S3='GGTGATAAAACGTGTAGCAAGCTGTAATCGAAGAGCATGCCACTACTATGGCG'
S4='CCTCGCATGACTCAACTGCCTGGTGATACGAGGCATGCTCTACCGGTATTGGAC'

for i in acs_DNA(S1,S2):
    print(i[0][1]-i[0][0])



















