"""
CS364, Pairwise Sequence Alignment
Purpose:
Making pairwise sequence alignment algorithms
By applying the Smith-Waterman algorithm, we can assess the similarities and differences in amino acid sequences across species, which may provide insights into their evolutionary relationships.

Author: jiyoon
"""

#-------------------------------------------------------------------------------
# MAIN
#-------------------------------------------------------------------------------
import numpy as np
import argparse

def read_fasta(file_path):
    """
    Read a FASTA file and return the sequence.
    file_path: name of the file 
    return: sequence from readline in file expect header(first line) 
    """
    try:
        with open(file_path) as f:
            lines = f.readlines()
        sequence = ''.join(line.strip() for line in lines[1:]) 
        return sequence
    except FileNotFoundError:
        raise ValueError(f"File not found: {file_path}")
    except Exception as e:
        raise ValueError(f"An error occurred while reading {file_path}: {e}")

def read_blosum62(file_path):
    """
    Read the BLOSUM62 substitution matrix from a file.
    if it is empty send user error message
    
    file_path: name of the file (in lab4, blosum-62.txt as default)
    
    return: amino_acids, matrix
    """
    with open(file_path) as f:
        lines = f.readlines()
    
    if not lines:
        raise ValueError("The substitution matrix file is empty.")

    amino_acids = lines[0].strip().split()
    matrix = np.zeros((len(amino_acids), len(amino_acids)), dtype=int)
    
    for i, line in enumerate(lines[1:]):
        values = line.strip().split()
        if len(values) != len(amino_acids) + 1:
            raise ValueError(f"Line {i + 2} does not match expected number of columns.")
        matrix[i] = list(map(int, values[1:]))  

    return amino_acids, matrix

def get_score(a, b, amino_acids, substitution_matrix):
    """
    Calculate the substitution score for given amino acids.
    
    a: seq1[x]
    b: seq2[y]
    amino_acids:amino_acids that created by def read_blosum62()
    substitution_matrix:substitution_matrix that created by def read_blosum62()
    
    return:substitution_matrix[index of a][index of b]
    """
    # check a or b is gaps
    if a == '-' or b == '-':
        return 0  
    i = amino_acids.index(a)
    j = amino_acids.index(b)
    return substitution_matrix[i][j]

def smith_waterman(seq1, seq2, amino_acids, substitution_matrix, gap_penalty):
    '''
    Building DP table and calculate max score and  max position
    
    seq1: first sequence
    seq2: second sequence
    amino_acids: list of valid amino acid characters
    substitution_matrix: 2D array that is scoring matrix for substitutions
    gap_penalty: penalty score for gaps
    
    return:
    dp: dynamic programming table that is 2D array
    max_score: maximum alignment score
    max_positions: list of positions where max_score occurs
    '''
    len1, len2 = len(seq1), len(seq2)
    dp = np.zeros((len1 + 1, len2 + 1), dtype=int)
    
    max_score = 0
    max_positions = []    
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            match = dp[i - 1][j - 1] + get_score(seq1[i - 1], seq2[j - 1], amino_acids, substitution_matrix)
            delete = dp[i - 1][j] + gap_penalty
            insert = dp[i][j - 1] + gap_penalty
            dp[i][j] = max(0, match, delete, insert)
            if dp[i][j] > max_score:
                max_score = dp[i][j]
                max_positions = [(i, j)]
            elif dp[i][j] == max_score:
                max_positions.append((i, j))
    return dp, max_score, max_positions

def back_traceback(dp, seq1, seq2, max_positions, amino_acids, substitution_matrix, gap_penalty):
    """
    Trace back through the DP table to find the optimal alignments.
    
    Parameters:
    seq1: First sequence
    seq2: Second sequence
    max_positions: List of starting positions of maximum scores
    amino_acids: List of amino acids
    substitution_matrix: Substitution matrix
    gap_penalty: Gap penalty
    
    Returns:
    alignments: List of tuples containing the alignments and pos1 & pos2
    """
    alignments = []
    for start in max_positions:
        i, j = start
        align1, align2 = [], []
        pos1, pos2 = [], []  # Track positions in the original sequences

        while dp[i][j] > 0:
            # Move DIAGONAL
            # Original position in seq1
            # Original position in seq2
            if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + get_score(seq1[i - 1], seq2[j - 1], amino_acids, substitution_matrix):  
                align1.append(seq1[i - 1])
                align2.append(seq2[j - 1])
                pos1.append(i) 
                pos2.append(j)  
                i -= 1
                j -= 1
            # Move UP
            # Original position in seq1
            # No position in seq2
            elif i > 0 and dp[i][j] == dp[i - 1][j] + gap_penalty:  
                align1.append(seq1[i - 1])
                align2.append('-')
                pos1.append(i)  
                pos2.append(None)  
                i -= 1
            # Move LEFT
            # No position in seq1
            # Original position in seq2
            elif j > 0 and dp[i][j] == dp[i][j - 1] + gap_penalty:  
                align1.append('-')
                align2.append(seq2[j - 1])
                pos1.append(None)  
                pos2.append(j)  
                j -= 1

            if dp[i][j] == 0:
                break

        # Reverse to correct order
        alignments.append((''.join(reversed(align1)), ''.join(reversed(align2)), list(reversed(pos1)), list(reversed(pos2))))

    return alignments

def write_output(alignments, max_score, output_file):
    """
    Write the alignments, score, identity, and gap information to the output file.
    alignments: contatin the alignments that should be save in the output file
    max_score: max score of the alignment
    output_file: name of file that will contain the output
    """
    with open(output_file, 'w') as f:
        f.write(f"Alignment score: {max_score}\n")
        for align1, align2, pos1, pos2  in alignments:
            f.write(align1 + '\n')
            f.write(align2 + '\n')
            
            # Calculate identity and gaps
            total_length = len(align2)
            matches = sum(1 for a, b in zip(align1, align2) if a == b and a != '-' and b != '-')
            gaps = align1.count('-') + align2.count('-')
            identity = (matches / total_length) * 100 if total_length > 0 else 0
            f.write(f"Identity: {identity:.2f}%  Gaps: {gaps}/{total_length}\n")
            f.write("\n")
            

def main():
    parser = argparse.ArgumentParser(description="Smith-Waterman Local Sequence Alignment")
    parser.add_argument('-a', required=True, help="File containing the first sequence")
    parser.add_argument('-b', required=True, help="File containing the second sequence")
    parser.add_argument('-m', required=True, help="Substitution matrix file")
    parser.add_argument('-g', type=int, required=True, help="Gap penalty")
    parser.add_argument('-o', required=True, help="Output file for alignments")
    
    args = parser.parse_args()

    seq1 = read_fasta(args.a)
    seq2 = read_fasta(args.b)
    amino_acids, substitution_matrix = read_blosum62(args.m)
    
    dp, max_score, max_positions = smith_waterman(seq1, seq2, amino_acids, substitution_matrix, args.g)

    print(f"Alignment score: {max_score}\n")
    
    alignments = back_traceback(dp, seq1, seq2, max_positions,amino_acids,substitution_matrix, args.g)
    for align1, align2, pos1, pos2 in alignments:
        start1 = pos1[0] + 1 if pos1[0] is not None else ''
        start2 = pos2[0] + 1 if pos2[0] is not None else ''
        align1_clean = ''.join(a for a in align1 if a != '-')
        align2_clean = ''.join(b for b in align2 if b != '-')
        end1 = start1 + len(align1_clean) - 2
        end2 = start2 + len(align2_clean) - 2

        # Print the first alignment
        print(f"{' ' * 4}{start1} {align1} {end1}")
        # visual representation of matches
        match_line = ''.join('|' if a == b and a != '-' else '.' if a != '-' and b != '-' else ' ' for a, b in zip(align1, align2))
        print(f"{' ' * 6}{match_line}")
        # Print the second alignment
        print(f"{' ' * 4}{start2} {align2} {end2}")

        total_length = len(align2) 
        matches = sum(1 for a, b in zip(align1, align2) if a == b and a != '-' and b != '-')
        gaps = align1.count('-') + align2.count('-') 
        print(f"Identity {matches}/{total_length} Gaps {gaps}/{total_length}\n")

    write_output(alignments, max_score, args.o)



#-------------------------------------------------------------------------------
# HELPERS
#-------------------------------------------------------------------------------

if __name__ == '__main__':
    main()