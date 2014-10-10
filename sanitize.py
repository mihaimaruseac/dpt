#!/usr/bin/env python
# -*- coding: utf8 -*-

import sys

import numpy
import scipy

class Args:
    def __init__(self, epsilon):
        self.epsilon = epsilon

def print_call(args):
    print "Called with\n\tepsilon = {:5.2f}".format(args.epsilon)

def usage():
    print sys.argv[0], ": Differentially private trajectory mining"
    print "Usage:"
    print "\t{0} epsilon".format(sys.argv[0])
    sys.exit(-1)

def convert_pairs_counts_to_matrix(pairs, nodes):
    """
    Converts a dictionary { (i, j): v } to numpy matrix m[i][j] = v.
    Assume 0 <= i,j < nodes.
    """
    ret = numpy.zeros([nodes + 1, nodes + 1])
    for k in pairs:
        ret[k] = pairs[k]
    return ret

def main():
    if len(sys.argv) != 2:
        usage()
    args = Args(float(sys.argv[1]))
    print_call(args)
    # TODO:  build graph
    nodes = 3
    lmax = 3
    graph = numpy.ones(nodes) - numpy.eye(nodes)
    print "Graph is:\n", graph
    transactions = [[1, 2], [1, 2, 3], [2, 3], [1, 3]]
    print "Transactions:\n", transactions
    real_pairs = {}
    for i in xrange(1, nodes+1):
        real_pairs[(0, i)] = real_pairs[(i, 0)] = 0
        for j in xrange(1, nodes+1):
            if graph[i-1][j-1]:
                real_pairs[(i, j)] = 0
    for t in transactions:
        t = [0] + t + [0] # add extra edge
        for edge in zip(t, t[1:]):
            if not real_pairs.has_key(edge):
                raise Exception("Invalid transaction found " + t)
            real_pairs[edge] += 1
    print "Real pair counts:\n", real_pairs
    print convert_pairs_counts_to_matrix(real_pairs, nodes)
    # add noise to pairs
    noisy_pairs = {}
    for k in real_pairs:
        noisy_pairs[k] = real_pairs[k] + numpy.random.laplace(scale=(lmax+1)/args.epsilon)
    noisy_pairs = {(0, 1): 3.03, (0, 2): 1.03, (0, 3): 0.27,
            (1, 0): -0.35, (1, 2): 0.93, (1, 3): -0.53,
            (2, 0): 1.03, (2, 1): 0, (2, 3): 2.01,
            (3, 0): 2.72, (3, 1): 0.18, (3, 2):-0.25}
    print "Noisy pair counts:\n", noisy_pairs
    print convert_pairs_counts_to_matrix(noisy_pairs, nodes)
    # compute approximation
    dim = nodes + 1
    equation_matrix = numpy.zeros([2 * dim, len(noisy_pairs)])
    b = numpy.zeros([2 * dim, 1])
    for i in xrange(nodes+1):
        s0 = s1 = 0
        for j in xrange(len(noisy_pairs)):
            p = noisy_pairs.keys()[j]
            if p[0] == i:
                equation_matrix[(2*i,j)] = 1
                s0 += noisy_pairs[p]
            if p[1] == i:
                equation_matrix[(2*i+1,j)] = 1
                s1 += noisy_pairs[p]
        s = (s0 + s1) / 2
        b[2*i] = s - s0
        b[2*i + 1] = s - s1
    print equation_matrix
    print b
    print numpy.linalg.matrix_rank(equation_matrix)
    sol = numpy.linalg.lstsq(equation_matrix, b)[0]
    print sol
    filtered_pairs = {}
    for j in xrange(len(noisy_pairs)):
        filtered_pairs[noisy_pairs.keys()[j]] = sol[(j,0)] + noisy_pairs.values()[j]
    print filtered_pairs
    print convert_pairs_counts_to_matrix(filtered_pairs, nodes)
    # get errors

if __name__ == '__main__':
    main()
