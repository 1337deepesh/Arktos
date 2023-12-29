# Arktos: a simple tool for the design of polyhedral DNA nanostructures

Copyright (C) 2023 Deepesh Nagarajan (1337deepesh@gmail.com)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License <http://www.gnu.org/licenses/>
for more details.

## Usage

python3 arktos.py <X_graph.txt>
Expect 'X_graph.OP.txt' as the output

There are no dependencies. The script should run out of the box.

## Concept

Arktos attempts to design simple DNA polyhedra using positive and negative design. The user inputs the DNA polyhedra they intend to design as a graph (format described later). This format represents individual single stranded complementary regions as nodes. ssDNA-ssDNA interactions that form dsDNA are represented as edges, or connections.

Once an input is read, Arktos attempts to design a DNA sequence that is expected to fold into the polyhedra using positive and negative selection.

Positive design involves trying to design on-target matches (base-pair complementarity) for every matching strand. This step is trivial.

Negative design ensures that off-target matches that arise by chance in unpaired strands are elimininated. Arktos does this using a simple scoring function:

Let UM = mismatch (length >=4) in unpaired strands.

ArkScore = Σ(1)(N+NC2) (UM)^2

For a system of N strands, there are N+NC2 total strand pairings (including self pairings at off-target regions). A summation of the squares of all off-target matches of sequence lengths >= 4 across all N+NC2 strands is the ArkScore. A sequence length of >=4 was chosen as the minimum threshold for a mismatch since  mismatch lengths <= 3 only transiently bind and will not affect the foldability of our designed polyhedron. A square term is added to disproportionately penalize off-target matches of greater lengths.

The ArkScore will be zero or positive. By design, it can never be negative. A simulated annealing algorithm is used to minimize the ArkScore in order to maximize folding into the designed polyhedral structure. 

Designed sequences are outputted in the fasta format, but with '|' breaks between different complementary regions (nodes) within a strand.

## Input

The full input format described her can be found in benchmarking_graph.txt

Assume you want to design a very simple DNA nanostructure like the one below. It's  not a polyhedron yet but the principles and methods used to design it are identical.


                 A
               5'| |3'
                 | |
                γ| |γ
                 | |
                 | |
                 ○ ○
                /   \
               / /○\ \
             α/ /   \ \β
             / /α   β\ \
          3'/ /       \ \5'
             /5'     3'\ C 
             B            
        

There are three strands: A, B, and C. Each strand has 11 residues.
- 5 residues on strand A pair with 5 residues on strand B (α-α)
- 5 residues on strand B pair with 5 residues on strand C (β-β)
- 5 residues on strand C pair with 5 residues on strand A (γ-γ)

There is also a single unpaired residue on every strand that is neccessary to allow the strand to turn.

### Input rules:

The input format for strand A is described below. The input format for strands B and C will follow the same rules.

!A→(γ.A,5,$any) (hingeA,1,A)  (α.A,5,$any)

- !: indicates that the line gives information about a single DNA strand, and not about complementary region pairing (node pairing)
- A: name of strand
- →: denotes end of strand name and start of strand data
- (γ.A,5,$any): gives information about a complementary region
- γ.A: name of complementary region
- 5: length of complementary region
- $any: sequence design rules (see below)

Likewise, (hingeA,1,A) represents a single residue adenine hinge, and (α.A,5,$any) represents another complementary region (node).

### Sequence design rules:

Every complementary region (node) need a set of sequences constraints before design.

- (node_name, len, ATGC): samples all 4 nucleotides at 25% frequency for each. Used as the standard for most complementary regions (nodes).
- (node_name, len, GCGCGCAT): Used to design GC rich regions. The frequency of every letter within the string denotes the frequency it will be sampled at during design. For G occurs 3 times. C occurs 3 times. G+C = 6, string length = 8. 6/8*100 = GC sampling at 75% frequency
- $GCrich=GCGCGCAT: $ denotes a variable that can be used in place of a sequence. (node_name, len, GCGCGCAT) becomes (node_name, len, $GCrich)

### Complementary region (node) pairing rules:

?γ.A,γ.C

- ?: incidates thhat the line gives information about about complementary region pairing (node pairing), and not information about a single DNA strand.
- γ.A: the name of the first complementary region. α.A was defined as belonging to strand A in A→(γ.A,5,$any) ...
- γ.C: the name of the second complementary region. γ.C was defined as belonging to strand C in C→ ... (γ.C,5,$any)

?γ.C,γ.A
?γ.A,γ.C

Complementary regions need not be redundantly described, but if redundant inputs aere provided (as in benchmarking_graph.txt) the script will filter out the redundancy.

# Output format:

The full output format described her can be found in benchmarking_graph.OP.txt

#ArkScore: Arktos score of the designed polyhedron. Lower scores are preferable. Scores can never be negative.

> \>s1
>
> |ATCTCATTAC|A|GCTGTCTTCACCTGGAGGTG|A|TTGTGGCTCC|AGGGGTAGCC|A|GTTCCGCCGG|

Strands are outputted in the .fasta format. '|' breaks denote different complementary regions (nodes) within a strand. '|' can be removed via a simple pipe:

cat benchmarking_graph.OP.fasta |tr -d '|'

These sequences can now be used for primer synthesis and experimental validation.
