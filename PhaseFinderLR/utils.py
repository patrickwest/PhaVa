#!/usr/bin/env python

'''
MitoScan- Utility functions and classes
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

    #forwardReads = 0
    #reverseReads = 0
    #ratio = 0

    def __init__(self, chr, leftStart, leftStop, rightStart, rightStop):
        self.chr = chr
        self.leftStart = leftStart
        self.leftStop = leftStop
        self.rightStart = rightStart
        self.rightStop = rightStop

# class SingleCopyGene:
#
#     ident = ""
#     name = ""
#     start = ""
#     stop = ""
#     strand = ""
#     seq = ""
#     translation = ""
#     dClass = ""
#
#     def __init__(self, name, ident, seq, translation):
#         self.name = name
#         self.ident = ident
#         self.seq = seq
#         self.translation = translation
#
# class Mito:
#
#     scg = []
#     completeness = 0
#     name = ""
#     seq = ""
#     dClass = ""
#     code = 0
#
#     def __init__(self, name, seq):
#         self.name = name
#         self.seq = seq
#
#     def exportSCG(self, out_base_fn):
#         # uses append so delete existing files of same name first
#         prot_fn = out_base_fn + '.protein.scg.faa'
#         rrna_fn = out_base_fn + '.rrna.scg.fna'
#         #if os.path.exists(prot_fn):
#         #    try:
#         #        os.remove(prot_fn)
#         #    except:
#         #        print("Error while deleting file ", prot_fn)
#         #if os.path.exists(rrna_fn):
#         #    try:
#         #        os.remove(rrna_fn)
#         #    except:
#         #        print("Error while deleting file ", rrna_fn)
#
#         prot_out_fh = open(out_base_fn + '.protein.scg.faa', 'a+')
#         rrna_out_fh = open(out_base_fn + '.rrna.scg.fna', 'a+')
#         for gene in self.scg:
#             if (gene.name != 'LSU_rRNA_bacteria' and gene.name != 'SSU_rRNA_bacteria' and gene.name != 'tRNA'):
#                 #print(str(gene.ident) + " " + str(gene.name))
#                 prot_out_fh.write('>' + gene.ident + ' ' + gene.name + '\n')
#                 prot_out_fh.write(gene.translation + '\n')
#             elif (gene.name != 'tRNA'):
#                 rrna_out_fh.write('>' + gene.ident + ' ' + gene.name + '\n')
#                 rrna_out_fh.write(gene.seq + '\n')
#         prot_out_fh.close()
#         rrna_out_fh.close()
#
#     def classifyDomain(self):
#         counts = {}
#         for gene in self.scg:
#             # skip tRNAs for classification
#             if gene.name != 'tRNA':
#                 if gene.dClass not in counts:
#                     counts[gene.dClass] = 0
#                 counts[gene.dClass] += 1
#         if len(counts) > 0:
#             self.dClass = max(counts.items(), key=operator.itemgetter(1))[0]
