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
    pad = [0]
    for start in xrange(size):
        for fragment in itertools.product(xrange(1, nodes+1), repeat=size - start):
            pad = [0] * start
            if is_path_in_graph(graph, fragment):
                ret[tuple(pad + list(fragment))] = 0
                ret[tuple(list(fragment) + pad)] = 0
    for t in transactions:
        t = pad + t + pad # add extra edges to start and end
        lists = [t[i:] for i in xrange(size)]
        for path in itertools.izip(*lists):
            if not ret.has_key(path):
                raise Exception("Invalid transaction found: {} {}".format(t, path))
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
    keys = noisy_counts.keys()
    path_fragments = list(set(itertools.chain.from_iterable(
        [[k[1:], k[:-1]] for k in keys])))
    dim = len(path_fragments)
    A = numpy.zeros([2 * dim, len(keys)])
    b = numpy.zeros([2 * dim])
    for i in xrange(dim):
        p = path_fragments[i]
        s_from = sum([noisy_counts[k] for k in keys if k[1:]  == p])
        s_to   = sum([noisy_counts[k] for k in keys if k[:-1] == p])
        s = sum([s_from, s_to]) / 2
        A[2*i]   = [1 if k[1:]  == p else 0 for k in keys]
        A[2*i+1] = [1 if k[:-1] == p else 0 for k in keys]
        b[2*i]   = s - s_from
        b[2*i+1] = s - s_to
    print "rank = {}, cond_num = {}".format(numpy.linalg.matrix_rank(A),
            numpy.linalg.cond(A))
    assert numpy.linalg.matrix_rank(A) <= 2 * dim - 1
    sol, residuals, rank, s = numpy.linalg.lstsq(A, b)
    print "residuals", residuals
    print "rank", rank
    print "s", s

    ret =  {}
    for j in xrange(len(keys)):
        k = keys[j]
        ret[k] = noisy_counts[k] + sol[j]
    return ret

def compute_relative_error(term1, term2):
    err = 0.0
    cnt = 0
    for k in term1.keys():
        if term1[k] != 0:
            cnt += 1
            err += abs(term1[k] - term2[k])/abs(0.0 + term1[k])
    return err/cnt

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
    print "Transactions:\n", transactions

    # Compute pair counts
    real_pairs = get_real_counts(graph, nodes, transactions, lmax, size)
    print "Real counts:"
    print real_pairs
    print postprocess(real_pairs, nodes)

    # Add noise to pairs
    noisy_pairs = add_noise(real_pairs, lmax + 1, epsilon)
    print "Noisy counts:"
    print noisy_pairs

    # Compute approximation
    filtered_pairs = postprocess(noisy_pairs, nodes)
    print "Final approximation:"
    print filtered_pairs

    # Get errors
    fr = compute_relative_error(filtered_pairs, real_pairs)
    nr = compute_relative_error(noisy_pairs, real_pairs)
    fn = compute_relative_error(filtered_pairs, noisy_pairs)
    print fr, nr, fn, fr <= nr

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
