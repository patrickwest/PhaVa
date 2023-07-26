#!/usr/bin/env python



'''
PhaVa- Utility functions and classes
'''

import sys
import os
#from collections import Counter
#from Bio.Seq import Seq
import operator

class IRDb:

    genome = {}
    genomeName = ""
    IRs = []
    genes = {}

    #def __init__(self, genome):
    #    self.genome = genome

class IR:

    chr = ""
    leftStart = 0
    leftStop = 0
    rightStart = 0
    rightStop = 0

    leftSeq = ""
    rightSeq = ""
    middleSeq = ""

    flankSize = 0
    fSeq = ""
    rSeq = ""

    forwardReads = {}
    reverseReads = {}
    ratio = {}

    geneOverlaps = []
    upstreamGene = "NA"
    upstreamStrand = "NA"
    upstreamDistance = 1000000000
    downstreamGene = "NA"
    downstreamStrand = "NA"
    downstreamDistance = 1000000000

    #forwardReads = 0
    #reverseReads = 0
    #ratio = 0

    def __init__(self, chr, leftStart, leftStop, rightStart, rightStop):
        self.chr = chr
        self.leftStart = leftStart
        self.leftStop = leftStop
        self.rightStart = rightStart
        self.rightStop = rightStop

