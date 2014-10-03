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

if __name__ == '__main__':
    main()
