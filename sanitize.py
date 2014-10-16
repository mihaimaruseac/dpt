#!/usr/bin/env python
# -*- coding: utf8 -*-

import sys

import numpy
import scipy

class Args:
    def __init__(self, epsilon):
        self.epsilon = epsilon

def usage():
    print sys.argv[0], ": Differentially private trajectory mining"
    print "Usage:"
    print "\t{0} epsilon".format(sys.argv[0])
    sys.exit(-1)

def print_call(args):
    print "Called with\n\tepsilon = {:5.2f}".format(args.epsilon)

def convert_pairs_counts_to_matrix(pairs, nodes):
    """
    Converts a dictionary { (i, j): v } to numpy matrix m[i][j] = v.
    Assume 0 <= i,j < nodes.
    """
    ret = numpy.zeros([nodes + 1, nodes + 1])
    for k in pairs:
        ret[k] = pairs[k]
    return ret

def get_real_counts(graph, nodes, transactions, lmax, size=2):
    """
    Returns the real counts for all subtransactions of length size. Counts
    only the possible subpaths in the given graph. Graph is given as a numpy
    matrix -- the adjacency matrix -- and transactions as a list of lists.
    """
    if size != 2:
        raise Exception("Cases for size > 2 not implemented yet!")
    # TODO: size, lmax
    ret = {}
    for i in xrange(1, nodes+1):
        ret[(0, i)] = ret[(i, 0)] = 0
        for j in xrange(1, nodes+1):
            if graph[i-1][j-1]:
                ret[(i, j)] = 0
    for t in transactions:
        t = [0] + t + [0] # add extra edge
        for edge in zip(t, t[1:]):
            if not ret.has_key(edge):
                raise Exception("Invalid transaction found " + t)
            ret[edge] += 1
    return ret

def add_noise(counts, sensitivity, epsilon):
    """
    Add Laplace noise to all items in counts dictionary.
    """
    ret = {}
    for k in counts:
        ret[k] = counts[k] + numpy.random.laplace(scale=sensitivity/epsilon)
    return ret

def postprocess(noisy_counts, nodes):
    """
    Post process the noisy_counts to ensure sum consistency across rows and
    columns of implied matrix.
    """
    dim = nodes + 1
    A = numpy.zeros([2 * dim, len(noisy_counts)])
    b = numpy.zeros([2 * dim, 1])
    keys = noisy_counts.keys()
    for i in xrange(nodes+1):
        s = [0, 0]
        for j in xrange(len(keys)):
            p = keys[j]
            for k in [0, 1]:
                if p[k] == i:
                    A[(2 * i + k, j)] = 1
                    s[k] += noisy_counts[p]
        avg_s = sum(s) / 2
        for k in [0, 1]:
            b[2 * i + k] = avg_s - s[k]
    assert numpy.linalg.matrix_rank(A) == 2 * dim - 1
    sol = numpy.linalg.lstsq(A, b)[0]

    ret =  {}
    for j in xrange(len(noisy_counts)):
        ret[noisy_counts.keys()[j]] = sol[(j,0)] + noisy_counts.values()[j]
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
    real_pairs = get_real_counts(graph, nodes, transactions, lmax)
    print "Real pair counts:\n", real_pairs
    print convert_pairs_counts_to_matrix(real_pairs, nodes)
    # add noise to pairs
    noisy_pairs = add_noise(real_pairs, lmax + 1, args.epsilon)
    print "Noisy pair counts:\n", noisy_pairs
    print convert_pairs_counts_to_matrix(noisy_pairs, nodes)
    # compute approximation
    filtered_pairs = postprocess(noisy_pairs, nodes)
    print filtered_pairs
    print convert_pairs_counts_to_matrix(filtered_pairs, nodes)
    # get errors

if __name__ == '__main__':
    main()
