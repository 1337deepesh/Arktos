#design a DNA sequence to fit the input network's topology, geometry, and connectivity:
import os
import re
import sys
import random

################################################################################
#FUNCTIONS: generate random DNA:
def randomDNA(length):
    sequence = ''.join([random.choice('ATGC') for i in range(0,length-2)])
    sequence = random.choice('GC')+sequence+random.choice('GC') #flank edges with GC for stronger bends
    return sequence

def GCrandomDNA(length):
    sequence = ''.join([random.choice('GCGCGCAT') for i in range(0,length-2)])
    sequence = random.choice('GC')+sequence+random.choice('GC') #flank edges with GC for stronger bends
    return sequence

################################################################################
#FUNCTION: find the reverse complement of a given DNA sequence:
def revcomp(sequence):
    sequence = sequence.replace('A','a').replace('T','t').replace('G','g').replace('C','c')
    sequence = sequence.replace('a','T').replace('t','A').replace('g','C').replace('c','G')
    sequence = sequence[::-1] #reverse
    return sequence

################################################################################
#FUNCTION: find all common substrings for 2 IP strings:
def acs(seq1,seq2):
    min_match = 4
    #initialize 2D list of seq1xseq2:
    seq_array = [[0]*len(seq1) for i in range(0,len(seq2))]
    for i in range(0,len(seq1)):
        for j in range(0,len(seq2)):
            if seq1[i]==seq2[j]:
                seq_array[j][i] = 1

    #for line in seq_array:
    #  print(line)

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

#FUNCTION: find the longest common DNA sequence with a complement present in the system:
def acs_DNA(seq1, seq2):

    #find the reverse complement of seq2:
    seq2 = revcomp(seq2)
    return acs(seq1,seq2)

################################################################################
#standard usage/error message:
if len(sys.argv) != 2:
    print("usage: "+sys.argv[0]+" <IP_graph.txt>")
    sys.exit()

#read network into python:
IP_graph = open(sys.argv[1]).readlines()

################################################################################
#data processing:

#extract edges from graph:
edges = {}
for line in IP_graph:
    if '?' in line[0]:
        edge = line.replace('?','').replace('\n','').split(',') #{A:B} edge format
        edges[edge[0]] = edge[1]

#find which strand each edge belongs to:
strand_edge = {}
for edge in edges:
    for line in IP_graph:
        if '!' in line[0] and edge in line:
            strand = line.replace('!','').split(':')[0] #extract name of DNA strand
            strand_edge[strand] = edge
            strand_edge[edge] = strand
            break

#extract node information from graph:
DNA_strands = {}

#format: {... strand:{order:[node1,hinge1,node2]
#                     node1:{partner:XXXX,len:X,GC:True/False,seq:XXXX},
#                     hinge1:'A',
#                     node2:{,partner:XXXX,len:X,GC:True/False,seq:XXXX},
#                     fullseq:XXXXXXXX, rules:['GC','ATGC','GCGCGCAT']}  ...}

for line in IP_graph:
    if '!' in line[0]:
        strand = line.replace('!','').split(':')[0] #extract name of DNA strand
        strand_nodes = re.findall(r'\(.*?\)',line.replace('!','').split(':')[1]) #extract data in braces
        strand_nodes = [i.replace('(','').replace(')','') for i in strand_nodes] #remove braces
        DNA_strands[strand] = {'order':[i.split(',')[0] for i in strand_nodes], 'fullseq':'', 'rules':[]} 
        for node in strand_nodes:
            GC = False
            if node.split(',')[1].isnumeric() == False: #if the node is a hinge:
                DNA_strands[strand][node.split(',')[0]] = node.split(',')[1]
                continue
            elif len(node.split(',')) == 2: #no GC
                [node, length] = node.split(',')
            elif len(node.split(',')) == 3: #GC given
                [node, length, GC] = node.split(',')
                GC = True
            DNA_strands[strand][node] = {'partner':edges[node], 
                                      'len':int(length), 'GC':GC, 'seq':''}
            #create a random DNA sequence for the above strand:
            partner = DNA_strands[strand][node]['partner']
            #if a complementary sequence DOES NOT exist:
            if strand_edge[partner] not in DNA_strands:
                if GC == False:
                    DNA_strands[strand][node]['seq'] = randomDNA(int(length))
                else:
                    DNA_strands[strand][node]['seq'] = GCrandomDNA(int(length))
            #if a complementary sequence already exists:
            else:
                if(DNA_strands[strand][node]['len'] != DNA_strands[strand_edge[partner]][partner]['len']):
                    print("error: nodes "+node+" and "+partner+" are of unequal lengths")
                DNA_strands[strand][node]['seq'] = revcomp(DNA_strands[strand_edge[partner]][partner]['seq'])

        #concatenate all sequences belonging to the same strand together:
        hinge_flag = False
        for node in DNA_strands[strand]['order']:
            if isinstance(DNA_strands[strand][node], dict): #if NOT a hinge
                DNA_strands[strand]['fullseq'] += DNA_strands[strand][node]['seq']
                #create rules for sequence design:
                if DNA_strands[strand][node]['GC'] == True:
                    for base in DNA_strands[strand][node]['seq']:
                        DNA_strands[strand]['rules'].append('GCGCGCAT')
                        if hinge_flag == True: #flank a hinge with GC residues on either end
                            DNA_strands[strand]['rules'][-1] = 'GC'
                            hinge_flag = False
                else:
                    for base in DNA_strands[strand][node]['seq']:
                        DNA_strands[strand]['rules'].append('ATGC')
                        if hinge_flag == True: #flank a hinge with GC residues on either end
                            DNA_strands[strand]['rules'][-1] = 'GC'
                            hinge_flag = False
            else: #if it's a hinge
                DNA_strands[strand]['fullseq'] += DNA_strands[strand][node]
                DNA_strands[strand]['rules'][-1] = 'GC' #flank a hinge with GC residues on either end
                for base in DNA_strands[strand][node]:
                    DNA_strands[strand]['rules'].append(base)
                hinge_flag = True


#check if all edges are symmetric (A->B and B->A):
sym_flag = False
for edge in edges:
    if edges[edge] not in edges:   
        print("error: edges "+edge+" and "+edges[edge]+" not symmetric")
        sys.exit()    
    if edges[edges[edge]] != edge:
        print("error: edges "+edge+" and "+edges[edge]+" not symmetric")
        sys.exit()

##############################################################################
#dictionary of all residue-residue pairings:

#calculate cumulative residue lengths of each node:
strand_len = {}
for strand in DNA_strands:

    node_len = {}
    counter = 0
    for node in DNA_strands[strand]['order']:
        if isinstance(DNA_strands[strand][node], dict): #if NOT a hinge
            node_len[strand+'|'+node] = [DNA_strands[strand][node]['len']+counter, DNA_strands[strand][node]['partner']]
            counter += DNA_strands[strand][node]['len']
        else:
            node_len[strand+'|'+node] = [len(DNA_strands[strand][node])+counter, None]
            counter += len(DNA_strands[strand][node])
    
        strand_len[strand] = node_len

#store residue lengths in this format: {'strand|node|base|rules':'strand|node|base|rules', ...}
res_list = []
for strand in strand_len:
    #print(strand_len[strand])
    counter = 0
    hinge_counter = 0
    hinge_flag = False
    for strandnode in strand_len[strand]:
        residue_count = strand_len[strand][strandnode][0]
        #print(residue_count)
        for res in range(counter, residue_count):
            res_strand = strand
            res_node = strandnode.split('|')[1]
            if isinstance(DNA_strands[res_strand][res_node], dict): #if NOT a hinge
                GC = DNA_strands[res_strand][res_node]['GC']
                if GC == True:
                    rules = 'GCGCGCAT'
                else:
                    rules = 'ATGC'
                if hinge_flag == True: #add GC residues after hinge:
                    rules = 'GC'
                hinge_flag = False
            else: #if it's a hinge
                if hinge_flag == False:
                    #add GC residues begore hinge:
                    prev = res_list[-1].split('|') 
                    prev = prev[0]+'|'+prev[1]+'|'+prev[2]+'|'+'GC'
                    res_list[-1] = prev
                    hinge_counter = 0
                    hinge_flag = True
                rules = DNA_strands[res_strand][res_node][hinge_counter]
                hinge_counter += 1    
            res_list.append(strandnode+'|'+str(res)+'|'+rules)
        
        counter = residue_count

#pair residues based on their node-pairings (edges):
for edge in edges:
    print(edge, edges[edge])

for res in res_list:
    print(res)

for edge1 in edges:
    for res in res_list:
        if edge1 in res:
            print(edge1, res)

'''
resres = {}
edge1_res = []
edge2_res = []
for edge1 in edges:
    for res in res_list:
        if edge1 in res:
            edge1_res.append(res)
    edge2 = edges[edge1]
    for res in res_list:
        if edge2 in res:
            edge2_res.append(res)
    edge2_res = edge2_res[::-1]
    #error message:
    if len(edge1_res) != len(edge2_res):
        print("error: edges "+edge1+" and "+edge2+" have unequal residue counts")
        sys.exit()

    for i in range(0,len(edge1_res)):
        resres[edge1_res[i]] = edge2_res[i]

for i in resres:
    print(i+'\t'+resres[i])
'''
    



##############################################################################

#sequence optimization:
match_string = ''
for i in DNA_strands:
    for j in DNA_strands:
            matches = acs_DNA(DNA_strands[i]['fullseq'],DNA_strands[j]['fullseq'])
            for k in matches:
                a = 1
                #print(k[0][1]-k[0][0])
#print(len(match_string))

##############################################################################
#output sequences:
OP_filename = sys.argv[1].split('.')
OP_filename = '.'.join(OP_filename[0:len(OP_filename)-1])+'.OP.fasta'
OP_file = open(OP_filename, 'w+')

for strand in DNA_strands:
    OP_file.write('>'+strand+'\n') #write fasta header/seqname
    concat_seq = ''
    for node in DNA_strands[strand]['order']: #write sequence
        if isinstance(DNA_strands[strand][node], str):
            concat_seq += DNA_strands[strand][node]+'|'
        else:
            concat_seq += DNA_strands[strand][node]['seq']+'|'
    concat_seq = concat_seq[:-1]# remove terminal '|'
    OP_file.write(concat_seq+'\n') #write sequence

##############################################################################
#NOTES:

#1962 nucleotides in total in this assembly
#4^5 = 1024
#4^6 = 4096


#print(acs('89gabcde','123456789abcde'))













