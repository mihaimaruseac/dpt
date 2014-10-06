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
    graph = numpy.ones(3) - numpy.eye(3)
    print "Graph is:\n", graph
    transactions = [(1, 2), (1, 2, 3), (2, 3), (1, 3)]
    print "Transactions:\n", transactions
    # add noise to pairs
    # compute approximation
    # get errors

if __name__ == '__main__':
    main()
