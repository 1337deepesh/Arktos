#design a DNA sequence to fit the input network's topology, geometry, and connectivity:
import os
import re
import sys
import copy
import math
import numpy
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
    #initialize 2D list of seq1 x seq2:
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

################################################################################
#FUNCTION: find all DNA complements in the system >=4, score them:
def acs_DNA(allseqs):

    seqscore = 0
    for i in range(0,len(allseqs)-1):
        for j in range(i+1,len(allseqs)):
            seq1 = allseqs[i]
            seq2 = revcomp(allseqs[j])
            acs_list = acs(seq1,seq2)
            #calculate seqscore. higher numbers indicate poorer specificity.
            for score in acs_list:
                my_score = score[0][1]-score[0][0]
                if my_score >= 4:
                    seqscore += math.pow(my_score,2)
    
    return seqscore

################################################################################
#FUNCTION: calculate baseline seqscore (expected from designed complementarity):
def baseline_seqscore(DNA_strands):
    
    seqscore = 0
    for strand in DNA_strands:
        for node in DNA_strands[strand]:
             if node in edges:
                seqscore += math.pow(len(DNA_strands[strand][node]),2)

    return seqscore/2 #divide to remove a->b b->a redundancy

################################################################################
#FUNCTION: extract DNA strands
def seq_extractor(DNA_strands):

    allseqs = []
    for strand in DNA_strands:
        fullseq = ''
        for node in DNA_strands[strand]:
            if node != 'fullseq' and node != 'order':
                residues = DNA_strands[strand][node]
                residues = numpy.sort(list(residues.keys()))
                for res in residues:
                    fullseq += DNA_strands[strand][node][res]['seq']
            DNA_strands[strand]['fullseq'] = fullseq
        allseqs.append(fullseq)

    return allseqs

################################################################################
#FUNCTION: mutate a single residue-pair anywhere along DNA_strands:
def mutate(DNA_strands):

    #pick a random strand, node, residue:
    strand = random.choice(list(DNA_strands.keys()))

    nodes = []
    for node in DNA_strands[strand]['order']:
        if node in edges:
            nodes.append(node)
    node = random.choice(nodes)
    residue = random.choice(list(DNA_strands[strand][node].keys()))
    #mutate a given residue:
    DNA_strands[strand][node][residue]['seq'] = random.choice(DNA_strands[strand][node][residue]['rule'])
    #also mutate the complement residue:
    mutant = DNA_strands[strand][node][residue]['seq']
    c = DNA_strands[strand][node][residue]['partner']
    DNA_strands[c['strand']][c['node']][c['res']]['seq'] = revcomp(mutant)

    return DNA_strands

################################################################################
#FUNCTION: cooling schedules for simualted annealing:
def cooling_scheduler(i, cycles, cooling_schedule, start_temperature, end_temperature):
    i = float(i)
    T0 = float(start_temperature)
    Tn = float(end_temperature)
    N = float(cycles)
    
    if cooling_schedule == 0:
        Ti = T0 -i*(T0-Tn)/N
    
    if cooling_schedule == 1:
        Ti = T0*(Tn/T0)**(i/N)
    
    if cooling_schedule == 2:
        A = (T0-Tn)*float(N+1)/N
        B = T0 -A
        Ti = A/(i+1) +B
    
    if cooling_schedule == 3:
        print("warning: cooling_schedule '3' does not work as described")
        print("switching to cooling_schedule '5'")
        cooling_schedule = 5
    
    if cooling_schedule == 4:
        Ti = (T0-Tn)/(1+math.exp(0.01*(i-N/2))) +Tn;
    
    if cooling_schedule == 5:
        Ti = 0.5*(T0 -Tn)*(1+math.cos(i*math.pi/N)) +Tn
    
    if cooling_schedule == 6:
        Ti = 0.5*(T0-Tn)*(1-math.tanh(i*10/N-5)) +Tn;
    
    if cooling_schedule == 7:
        Ti = (T0-Tn)/math.cosh(i*10/N) +Tn;
    
    if cooling_schedule == 8:
        A = (1/N)*math.log(T0/Tn)
        Ti = T0*math.exp(-A*i)

    if cooling_schedule == 9:
        A = (1/N**2)*math.log(T0/Tn)
        Ti = T0*math.exp(-A*i**2);
    
    return Ti

################################################################################
#FUNCTION: acceptance decision for simulated annealing:
def acceptance_decision(current_score, previous_score, temperature, probability_flag):
    if current_score < previous_score:
        return "accept"
    else:
        if previous_score == 0:
            previous_score = -0.001
        if current_score == 0:
            current_score = -0.001
        score_ratio = current_score/previous_score
        if probability_flag == True:
            acceptance_threshold = temperature*(  1/(1+1000*math.exp(-12*score_ratio))  )
        else:
            acceptance_threshold = temperature
        #acceptance_threshold = temperature
        if random.random() < acceptance_threshold:
            return "accept"
        else:
            return "reject"

################################################################################
#standard usage/error message:
if len(sys.argv) != 2:
    print("usage: "+sys.argv[0]+" <IP_graph.txt>")
    sys.exit()

#read network into python:
IP_graph = open(sys.argv[1]).readlines()

################################################################################
# ?: extract edges from graph:
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
            strand = line.replace('!','').split('→')[0] #extract name of DNA strand
            strand_edge[strand] = edge
            strand_edge[edge] = strand
            break

################################################################################
# $: extract variables for graph:
variables = {}
for line in IP_graph:
    if '$' in line[0]:
        variable = line.replace(' ', '').replace('\n','').split('=')
        variables[variable[0]] = variable[1]

################################################################################
# !: extract node information from graph:
DNA_strands = {}

for line in IP_graph:
    if '!' in line[0]:
        #extract name of DNA strand
        strand = line.replace(' ','').replace('!','').split('→')[0] 
        #extract data in braces
        nodes = re.findall(r'\(.*?\)',line.replace(' ','').replace('!','').split('→')[1])
        #remove braces
        nodes = [i.replace('(','').replace(')','') for i in nodes] 

        DNA_strands[strand] = {'order':[i.split(',')[0] for i in nodes], 'fullseq':''} 

        for node in nodes:
            [node, seq_length, rules] = node.split(',')
            seq_length = int(seq_length)
            
            #extract rules, build a sequence from them:
            rules = rules.replace('[','').replace(']','')
            #replace any variables with their values:
            for var in variables:
                if var in rules:
                    rules = rules.replace(var,variables[var])
            rules = rules.split(';')
            sequence = {}
            counter = 0
            #handle nonpositional rule first:
            for rule in rules:
                if ':' not in rule:
                    for i in range(0,seq_length):
                        sequence[i] = {'seq':random.choice(rule), 'rule':rule, 'partner':None}
            #overwrite sequence with positional rules:
            for rule in rules:
                if ':' in rule:
                    [position,ruleseq] = rule.split(':')
                    position = int(position)
                    if position<0:
                        position += seq_length
                    sequence[position] = {'seq':random.choice(ruleseq), 
                                           'rule':ruleseq, 'partner':None} 
            DNA_strands[strand][node] = sequence
            #compile full sequence:
            residues = numpy.sort(list(sequence.keys()))
            fullseq = ''
            for i in residues:
                fullseq += sequence[i]['seq']
            DNA_strands[strand]['fullseq'] = fullseq

#assign partners for every residue:
for strand in DNA_strands:
    for node in DNA_strands[strand]:
        if node != 'order' and node != 'fullseq':
            if node in edges:
                
                #calculate sequence length for node:
                node_length = DNA_strands[strand][node].keys()
                node_length = max(node_length)+1

                #find parther
                partner_node = edges[node]
                partner_strand = strand_edge[partner_node]

                #calculate sequence length for partner:
                print(partner_strand, partner_node)
                partner_length = DNA_strands[partner_strand][partner_node].keys()
                partner_length = max(partner_length)+1
                if node_length != partner_length:
                    print("error: nodes "+node+" and "+partner_node+
                          " are of unequal lengths")
                    sys.exit()
                #assign partners to every residue:
                for i in range(0,node_length):
                    DNA_strands[strand][node][i]['partner'] = {'strand':partner_strand,
                                                               'node':partner_node, 
                                                               'res':node_length-i-1}

##############################################################################
#enforce complementarity on partner:

for strand in DNA_strands:
    for node in DNA_strands[strand]:
        if node in edges:
            for residue in DNA_strands[strand][node]:
                
                #mutate a given residue:
                DNA_strands[strand][node][residue]['seq'] = random.choice(DNA_strands[strand][node][residue]['rule'])
                
                #also mutate the complement residue:
                mutant = DNA_strands[strand][node][residue]['seq']
                c = DNA_strands[strand][node][residue]['partner']
                DNA_strands[c['strand']][c['node']][c['res']]['seq'] = revcomp(mutant)
                
                #check for rule complementarity as well:
                rule1 = DNA_strands[c['strand']][c['node']][c['res']]['rule']
                rule2 = DNA_strands[strand][node][residue]['rule']
                if rule1 != rule2:
                    print("error: residues "+c['node']+"|"+str(c['res'])+" and "+node+"|"+str(residue)+" have different design rules.")
                    sys.exit()

##############################################################################
#count total number of mutatable residues:

mutate_count = 0
for strand in DNA_strands:
    for node in DNA_strands[strand]:
        if node in edges:
            mutate_count += len(DNA_strands[strand][node])

##############################################################################
#sequence optimization:
baseline = baseline_seqscore(DNA_strands)
cycles = int(mutate_count/2*10) #sample each residue 10 times
temperature = 1
cooling_schedule = 7
curr_score = float("inf")
prev_score = float("inf")
best_score = float("inf")
curr_state = copy.deepcopy(DNA_strands)
prev_state = {}
best_state = {}

#simple gradient minimization:
for i in range(0,cycles):

    #calculate global temperature:
    temperature = cooling_scheduler(i, cycles, cooling_schedule, 1, 1/cycles)
    #backup existing state:
    print("score: "+str(curr_score)+", "+str(i)+"/"+str(cycles))
    prev_score = curr_score
    prev_state = copy.deepcopy(curr_state)

    #mutate a single residue pair & calculate seqscore:
    curr_state = mutate(curr_state)
    curr_score = acs_DNA(seq_extractor(curr_state))-baseline

    decision = acceptance_decision(curr_score, prev_score, temperature, False)

    if decision == "accept":     
        #check if the current state occupies the lowest energy:
        if curr_score <= best_score:
            best_score = curr_score
            best_state = copy.deepcopy(curr_state)
    if decision == "reject":
        #revert to previous state:
        curr_score = prev_score
        curr_state = copy.deepcopy(prev_state)
    
##############################################################################
#output sequences:
DNA_strands = copy.deepcopy(best_state)
OP_filename = sys.argv[1].split('.')
OP_filename = '.'.join(OP_filename[0:len(OP_filename)-1])+'.OP.fasta'
OP_file = open(OP_filename, 'w+')

OP_file.write("#ArkScore: "+str(best_score)+"\n")

for strand in DNA_strands:
    OP_file.write('>'+strand+'\n') #write fasta header/seqname
    fullseq = ''
    for node in DNA_strands[strand]:
        if node != 'fullseq' and node != 'order':
            residues = DNA_strands[strand][node]
            residues = numpy.sort(list(residues.keys()))
            fullseq += '|'
            for res in residues:
                fullseq += DNA_strands[strand][node][res]['seq']
    
    OP_file.write(fullseq+'|'+'\n') #write sequence

OP_file.close()

##############################################################################
#NOTES:

#1962 nucleotides in total in this assembly
#4^5 = 1024
#4^6 = 4096


#print(acs('89gabcde','123456789abcde'))

#python3 arktos.py IP_small.txt|grep "⚫.α.1.1"









