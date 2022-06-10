#!/usr/bin/env python

import logging
import PhaseFinderLR.utils
import PhaseFinderLR.locate
import os
from Bio.Seq import Seq

def main(args, irDb):

    #if args.irs is not None:
    #    irDb.genome = PhaseFinderLR.locate.parseGenome(args.fasta)
    #    irDb.genomeName = os.path.basename(args.fasta)
    #    irDb = parseIRs(args.irs, irDb)
    #    PhaseFinderLR.locate.exportIRs(irDb.IRs, args.dir)

    irDb = createInSilico(args, irDb)
    exportMockInvertons(irDb, args.dir)

    if (args.mockGenome):
        createMockIRGenome(args, irDb)

    return irDb


def createInSilico(args, irDb):

    IRs = irDb.IRs
    for ir in irDb.IRs:
        IRs[ir].flankSize = args.flankSize

        leftFlankStart = irDb.IRs[ir].leftStart - IRs[ir].flankSize
        if leftFlankStart < 0:
            leftFlankStart = 0
        rightFlankEnd = irDb.IRs[ir].rightStop + IRs[ir].flankSize

        leftFlank = irDb.genome[irDb.IRs[ir].chr][leftFlankStart:irDb.IRs[ir].leftStart]
        #rIR = str(Seq(irDb.IRs[ir].leftSeq + irDb.IRs[ir].middleSeq + irDb.IRs[ir].rightSeq).reverse_complement())
        rIR = str(irDb.IRs[ir].leftSeq + str(Seq(irDb.IRs[ir].middleSeq).reverse_complement()) + irDb.IRs[ir].rightSeq)
        rightFlank = irDb.genome[irDb.IRs[ir].chr][irDb.IRs[ir].rightStop:rightFlankEnd]

        IRs[ir].fSeq = leftFlank + irDb.IRs[ir].leftSeq + irDb.IRs[ir].middleSeq + irDb.IRs[ir].rightSeq + rightFlank
        IRs[ir].rSeq = leftFlank + rIR + rightFlank

    return irDb

def parseIRs(path, irDb):
    IRs = {}

    ir_fh = open(path)
    for line in ir_fh:
        sline = line.strip().split()
        id = sline[0] + ':' + str(sline[1]) + '-' + str(sline[2]) + '-' + str(sline[3]) + '-' + str(sline[4])
        IRs[id] = PhaseFinderLR.utils.IR(sline[0], int(sline[1]), int(sline[2]), int(sline[3]), int(sline[4]))
        IRs[id].leftSeq = sline[5]
        IRs[id].middleSeq = sline[6]
        IRs[id].rightSeq = sline[7]

    irDb.IRs = IRs
    return irDb

def createMockIRGenome(args, irDb):
    mockG = irDb.genome
    for ir in irDb.IRs:
        mockG[irDb.IRs[ir].chr] = mockG[irDb.IRs[ir].chr][0:irDb.IRs[ir].leftStop] + str(Seq(irDb.IRs[ir].middleSeq).reverse_complement()) + mockG[irDb.IRs[ir].chr][irDb.IRs[ir].rightStart:len(mockG[irDb.IRs[ir].chr])]

    out = open(args.dir + "/mockGenome.fasta", 'w')
    for chrom in mockG:
        out.write('>' + chrom + "\n")
        out.write(mockG[chrom] + "\n")

    out.close()


def exportMockInvertons(irDb, outpath):
    out = open(outpath + '/invertedSeqs.fasta', 'w')
    for ir in irDb.IRs:
        out.write('>' + ir + '_f\n')
        out.write(irDb.IRs[ir].fSeq + '\n')
        out.write('>' + ir + '_r\n')
        out.write(irDb.IRs[ir].rSeq + '\n')
