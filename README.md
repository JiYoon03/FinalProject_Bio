# FinalProject_Bio
explore the evolutionary relationships between the selected members of the two families, Felines (tiger, lion, leopard) and Canidae (coyote, wolf, fox). 

## Pairwise

Input: Two sequences in FASTA format, a substitution matrix (blosum-62.txt from lab 4), and a gap penalty(3).

>$ python3 pairwise.py -a data_pw/<data>.fasta -b data_pw/<data>.fasta -m blosum-62.txt -g -3 -o <output>.txt 

Steps:
1) Read Sequences & Matrix: Parse input sequences and substitution matrix.
2) Dynamic Programming Table: Compute scores for alignments using the Smith-Waterman algorithm.
3) Traceback: Identify optimal local alignments by backtracking through the DP table.

Output:
Alignment sequences with matching regions and gap information.
Maximum alignment score, identity percentage, and gap count.
