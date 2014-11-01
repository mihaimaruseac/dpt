#!/usr/bin/env python
# -*- coding: utf8 -*-

import itertools
import numpy
import random
import scipy
import sys

def usage():
    print sys.argv[0], ": Differentially private trajectory mining"
    print "Usage:"
    print "\t{0} epsilon nodes lmax tmax size [seed]".format(sys.argv[0])
    sys.exit(-1)

def print_call(*args):
    print "Called with\n\
            \tepsilon = {:5.2f}\n\
            \tnodes = {:d}\n\
            \tlmax = {:d}\n\
            \ttmax = {:d}\n\
            \tsize = {:d}\n\
            \tseed = {:d}\n\
            ".format(*args)

def convert_pairs_counts_to_matrix(pairs, size):
    """
    Converts a dictionary { (i, j): v } to numpy matrix m[i][j] = v.
    Assume 0 <= i,j <= size.
    """
    ret = numpy.zeros([size + 1, size + 1])
    for k in pairs:
        ret[k] = pairs[k]
    return ret

def is_path_in_graph(graph, path):
    """
    Returns True if path is a valid path in graph. Graph is given by the
    adjancency matrix as a numpy matrix.
    """
    for i, j in zip(path, path[1:]):
        if not graph[(i - 1, j - 1)]:
            return False
    return True

def get_real_counts(graph, nodes, transactions, lmax, size=2):
    """
    Returns the real counts for all subtransactions of length size. Counts
    only the possible subpaths in the given graph. Graph is given as a numpy
    matrix -- the adjacency matrix -- and transactions as a list of lists.
    """
    ret = {}
    for fragment in itertools.product(xrange(1, nodes+1), repeat=size - 1):
        if is_path_in_graph(graph, fragment):
            ret[tuple([0] + list(fragment))] = 0
            ret[tuple(list(fragment) + [0])] = 0
    for path in itertools.product(xrange(1, nodes+1), repeat=size):
        if is_path_in_graph(graph, path):
            ret[path] = 0
    for t in transactions:
        t = [0] + t + [0] # add extra edges to start and end
        lists = [t[i:] for i in xrange(size)]
        for path in itertools.izip(*lists):
            if not ret.has_key(path):
                raise Exception("Invalid transaction found: {}".format(t))
            ret[path] += 1
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
        k = keys[j]
        ret[k] = sol[(j,0)] + noisy_counts[k]
    return ret

def env_setup(seed):
    numpy.set_printoptions(formatter={'float': '{: 0.3f}'.format})
    random.seed(seed)
    numpy.random.seed(seed)

def main(epsilon, nodes, lmax, tmax, size, seed):
    env_setup(seed)

    # Read graph
    graph = numpy.ones(nodes) - numpy.eye(nodes)
    print "Graph is:\n", graph

    # Read transactions
    frequencies = list(itertools.chain.from_iterable([
        [i] * i for i in xrange(1, nodes + 1)]))
    transactions = []
    for i in xrange(tmax):
        tlen = random.choice(xrange(3, lmax + 1))
        t = [random.choice(frequencies)]
        while len(t) < tlen:
            e = random.choice([n for n in frequencies if graph[(t[-1]-1, n-1)]])
            t.append(e)
        transactions.append(t)
    #print "Transactions:\n", transactions

    # Compute pair counts
    real_pairs = get_real_counts(graph, nodes, transactions, lmax, size)
    print "Real counts:"
    print real_pairs

    # Add noise to pairs
    noisy_pairs = add_noise(real_pairs, lmax + 1, epsilon)
    print "Noisy counts:"
    print noisy_pairs

    # Compute approximation
    filtered_pairs = postprocess(noisy_pairs, nodes)
    print "Final approximation:"
    print filtered_pairs

    # Get errors
    P = convert_pairs_counts_to_matrix(real_pairs, nodes)
    P_star = convert_pairs_counts_to_matrix(noisy_pairs, nodes)
    P_hat = convert_pairs_counts_to_matrix(filtered_pairs, nodes)
    fr = numpy.linalg.norm(P_hat - P, ord='fro')
    nr = numpy.linalg.norm(P_star - P, ord='fro')
    fn = numpy.linalg.norm(P_hat - P_star, ord='fro')
    print fr, nr, fn, fr <= nr + fn, fr <= nr
    fre = 0.0
    cnt = 0
    for i in xrange(nodes + 1):
        for j in xrange(nodes + 1):
            loc = (i, j)
            if P[loc] != 0:
                cnt += 1
                fre += abs(P_hat[loc] - P[loc])/P[loc]
    fre /= cnt
    nre = 0.0
    cnt = 0
    for i in xrange(nodes + 1):
        for j in xrange(nodes + 1):
            loc = (i, j)
            if P[loc] != 0:
                cnt += 1
                nre += abs(P_star[loc] - P[loc])/P[loc]
    nre /= cnt
    fne = 0.0
    cnt = 0
    for i in xrange(nodes + 1):
        for j in xrange(nodes + 1):
            loc = (i, j)
            if P_star[loc] != 0:
                cnt += 1
                fne += abs(P_hat[loc] - P_star[loc])/P_star[loc]
    fne /= cnt
    print fre, nre, fne, fre <= nre + fne, fre <= nre

if __name__ == '__main__':
    if len(sys.argv) not in [6, 7]:
        usage()
    try:
        epsilon = float(sys.argv[1])
        nodes = int(sys.argv[2])
        lmax = int(sys.argv[3])
        tmax = int(sys.argv[4])
        size = int(sys.argv[5])
        seed = 42 if len(sys.argv) == 6 else int(sys.argv[6])
    except Exception as e:
        print e
        usage()
    print_call(epsilon, nodes, lmax, tmax, size, seed)
    main(epsilon, nodes, lmax, tmax, size, seed)
