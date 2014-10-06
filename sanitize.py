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
    real_pairs_array = numpy.zeros([nodes + 1, nodes + 1])
    for k in real_pairs:
        real_pairs_array[k] = real_pairs[k]
    print real_pairs_array
    # add noise to pairs
    noisy_pairs = {}
    for k in real_pairs:
        noisy_pairs[k] = real_pairs[k] + numpy.random.laplace(scale=(lmax+1)/args.epsilon)
    noisy_pairs = {(0, 1): 3.03, (0, 2): 1.03, (0, 3): 0.27,
            (1, 0): -0.35, (1, 2): 0.93, (1, 3): -0.53,
            (2, 0): 1.03, (2, 1): 0, (2, 3): 2.01,
            (3, 0): 2.72, (3, 1): 0.18, (3, 2):-0.25}
    print "Noisy pair counts:\n", noisy_pairs
    noisy_pairs_array = numpy.zeros([nodes + 1, nodes + 1])
    for k in noisy_pairs:
        noisy_pairs_array[k] = noisy_pairs[k]
    print noisy_pairs_array
    # compute approximation
    dim = nodes + 1
    equation_matrix = numpy.zeros([2 * dim, dim * (dim - 1)])
    b = numpy.zeros([2 * dim, 1])
    for i in xrange(nodes+1):
        for j in xrange(nodes+1):
            pass
    #sums = {}
    #for i in xrange(nodes + 1):
    #    s1 = s2 = 0
    #    for k in noisy_pairs:
    #        if k[0] == i: s1 += noisy_pairs[k]
    #        if k[1] == i: s1 += noisy_pairs[k]
    #    sums[i] = (s1 + s2) / 2
    #equation_matrix = numpy.zeros([2 * len(sums), 2 * len(sums)])
    #b = numpy.zeros([2 * len(sums), 1])
    #for i in xrange(2 * len(sums)):
    #    j = i % len(sums)
    #    b[(i,0)] = sums[j]
    print equation_matrix
    print b
    # get errors

if __name__ == '__main__':
    main()
