"""
CS364, Lab 5 Phylogenetic Trees
Purpose: File to implement building Neighbor Joining phylogenetic tree. 
Author: Jishi Lin
Date: 10/28/24
"""
import pygraphviz as pgv
from copy import deepcopy
import numpy as np
import optparse
import os
import sys

#-------------------------------------------------------------------------------
# MAIN
#-------------------------------------------------------------------------------

def main():
    opts = parse_args()
    delta = opts.dissimilarity
    output_filename = opts.output

    # run Neighbor-Joining algorithm
    nj_tree = NeighborJoiningTree(delta)
    nj_tree.nj(output_filename)
    print(nj_tree)

#-------------------------------------------------------------------------------
# HELPERS
#-------------------------------------------------------------------------------

class NeighborJoiningTree:
    def __init__(self, delta):
        self.delta = self.read_dissimilarity(delta)
        # deepcopy reference: https://stackoverflow.com/questions/63532921/
        # deep-copy-of-list-in-python
        # save initial delta for display
        self.initial_delta = deepcopy(self.delta)
        self.tree = pgv.AGraph()
        self.nodes = list(self.delta.keys())
        self.edge_labels = {}
        self.iterations_info = []

    def nj(self, output_filename):
        """
        Contruct Neighbor-Joning phylogenetic tree with a top-down hierarchical 
        structure.
        Input: filename (str)
        """
        while len(self.nodes) > 2:
            # compute the Q matrix and find f, g minimizing Q(i, j)
            Q = self.calculate_q()
            f, g = self.min_q_pair(Q)
            n = len(self.nodes)

            # create a new node v and compute distances
            d_fg = self.delta[f][g]
            S_f = sum(self.delta[f][k] for k in self.nodes if k != f)
            S_g = sum(self.delta[g][k] for k in self.nodes if k != g)

            d_fv = 0.5 * d_fg + 1 / (2 * (n - 2)) * (S_f - S_g)
            d_gv = 0.5 * d_fg + 1 / (2 * (n - 2)) * (S_g - S_f)

            # create the new node v and add it to the distance matrix
            v = f"{f}-{g}"
            # initialize new entry for v in delta
            self.delta[v] = {}
            self.tree.add_node(v, label="", fillcolor="blue", style="filled",
                               fontsize=40)
            self.tree.add_edge(v, f, label=f"{d_fv:.2f}", len=d_fv/5, fontsize=15)
            self.tree.add_edge(v, g, label=f"{d_gv:.2f}", len=d_gv/5, fontsize=15)

            # calculate distances from v to all remaining nodes
            for i in self.nodes:
                if i != f and i != g:
                    d_iv = 0.5 * (self.delta[f][i] - d_fv) + 0.5 * (
                        self.delta[g][i] - d_gv)
                    if v not in self.delta:
                        self.delta[v] = {}
                    if i not in self.delta:
                        self.delta[i] = {}
                    self.delta[v][i] = d_iv
                    self.delta[i][v] = d_iv

            # update distance matrix: remove f and g, add v
            del self.delta[f]
            del self.delta[g]
            for k in list(self.delta.keys()):
                if f in self.delta.get(k, {}):
                    del self.delta[k][f]
                if g in self.delta.get(k, {}):
                    del self.delta[k][g]
            self.nodes.remove(f)
            self.nodes.remove(g)
            self.nodes.append(v)

            # store iteration details for printing
            self.iterations_info.append({"merge": (f, g), "matrix": 
                                         self.format_matrix()})

        # termination: connect the last 2 nodes
        remaining_nodes = list(self.nodes)
        if len(remaining_nodes) == 2:
            a, b = remaining_nodes
            d_ab = self.delta[a][b]
            self.tree.add_edge(a, b, label=f"{d_ab:.2f}", fontsize=15)

        # save tree as png
        self.tree.draw(output_filename, prog="neato")

    def calculate_q(self):
        """
        Compute the Q matrix for NJ.
        Return: Q (dictionary)
        """
        # total number of nodes in tree
        n = len(self.nodes)
        # S[i] for each node i
        # S[i] is sum of distance from i to all other nodes
        S = {i: sum(self.delta[i][j] for j in self.nodes if j != i) 
             for i in self.nodes}
        Q = {}

        # iterate each pair of nodes to compute Q[i][j]
        for i in self.nodes:
            Q[i] = {}
            for j in self.nodes:
                if i != j:
                    Q[i][j] = (n - 2) * self.delta[i][j] - S[i] - S[j]
        return Q

    def min_q_pair(self, Q):
        """
        Finds pair of nodes with minimum Q score from Q matrix.
        Input: Q (dictioanry)
        Return: tuple pair of nodes
        """
        # initialize to infinity
        min_q = float('inf')
        # intialize min pair
        min_pair = (None, None)
        # iterate each pair in Q matrix
        for i in Q:
            for j in Q[i]:
                # check if Q score for pair is smaller than current
                if Q[i][j] < min_q:
                    # update score and pairs
                    min_q = Q[i][j]
                    min_pair = (i, j)
        return min_pair

    def read_dissimilarity(self, filename):
        """
        Read and return dissimilarity map (delta) for file.
        Input: filename (str)
        Return: delta (dictionary)
        """
        delta = {}

        # read file
        with open(filename, 'r') as file:
            column_headers = file.readline().split()
            for row in file:
                row = row.split() # remove spaces
                # row labels
                row_headers = row[0]
                values = row[1:]
                # inner dictionary for each column
                inner = {}
                for i in range(len(column_headers)):
                    inner[column_headers[i]] = float(values[i])
                
                delta[row_headers] = inner

        return delta
        
    def format_matrix(self, initial=False):
        """
        Formats the delta matrix for printing.
        Input: boolean for original delta matrix printing
        Return: list of strings for header and rows
        """
        delta_matrix = self.initial_delta if initial else self.delta
        clusters = list(delta_matrix.keys())
        header = ["\t" + "\t".join(clusters)]
        rows = []
        for row in clusters:
            row_values = "\t".join(f"{delta_matrix[row].get(col, 0.00):>5.2f}" 
                                   for col in clusters)
            rows.append(f"{row:>3}\t{row_values}")
        return header + rows

        
    def __str__(self):
        output = ["\nWelcome to NJ phylogenetic tree builder!\n"]
        # initial input dissimilarity map
        output.append("\nInput dissimilarity map:\n")
        output += self.format_matrix(initial=True)
        
        # iteration informations
        for iteration, data in enumerate(self.iterations_info, start=1):
            output.append(f"\n--------------")
            output.append(f"Iteration {iteration}")
            output.append(f"--------------\n")
            output.append(f"Merging: {data['merge'][0]} {data['merge'][1]}\n")
            output += data["matrix"]

        return "\n".join(output)
    
# from previous labs
def parse_args():
    """
    Parse command line arguments.
    """
    parser = optparse.OptionParser(description='UPGMA phylogentic tree builder')
    
    parser.add_option('-d', '--dissimilarity', type='string',
                      help='File containing dissimilarity map')
    parser.add_option('-t', '--output', type='string',
                      help='Output filename for tree image')

    (opts, args) = parser.parse_args()

    mandatories = ['dissimilarity', 'output']
    for m in mandatories:
        if not opts.__dict__[m]:
            print('mandatory option ' + m + ' is missing\n')
            parser.print_help()
            sys.exit()

    return opts


if __name__ == '__main__':
    main()
