#!/usr/bin/env python

import logging
import PhaVa.utils
import PhaVa.locate
import os
import random
import warnings
from Bio.Seq import Seq
from Bio import SeqIO

def main(args, irDb):

    if args.irs is not None:
        irDb.genome = PhaVa.locate.parseGenome(args.fasta)
        irDb.genomeName = os.path.basename(args.fasta)
        irDb = parseIRs(args.irs, irDb)
        PhaVa.locate.exportIRs(irDb.IRs, args.dir)

    irDb = createInSilico(args, irDb)
    exportMockInvertons(irDb, args.dir)

    if args.genes is not None:
        if args.genesFormat == 'gff':
            genes = parseGFF(args.genes)
        else:
            genes = parseGBK(args.genes)

        findGeneOverlaps(genes, irDb)
        exportGeneOverlaps(irDb, args.dir)

    if (args.mockGenome is not None and args.mockGenome):
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

    irNames = []
    if args.mockNumber > 0:
        irs = random.sample(irDb.IRs.keys(), args.mockNumber)
        for ir in irs:
            mockG[irDb.IRs[ir].chr] = mockG[irDb.IRs[ir].chr][0:irDb.IRs[ir].leftStop] + str(Seq(irDb.IRs[ir].middleSeq).reverse_complement()) + mockG[irDb.IRs[ir].chr][irDb.IRs[ir].rightStart:len(mockG[irDb.IRs[ir].chr])]
            irNames.append(ir)
    else:
        for ir in irDb.IRs:
            mockG[irDb.IRs[ir].chr] = mockG[irDb.IRs[ir].chr][0:irDb.IRs[ir].leftStop] + str(Seq(irDb.IRs[ir].middleSeq).reverse_complement()) + mockG[irDb.IRs[ir].chr][irDb.IRs[ir].rightStart:len(mockG[irDb.IRs[ir].chr])]
            irNames.append(ir)

    out = open(args.dir + "/mockGenome.fasta", 'w')
    outTsv = open(args.dir + "/mockGenomeInvertedInvertons.tsv", 'w')
    for chrom in mockG:
        out.write('>' + chrom + "\n")
        out.write(mockG[chrom] + "\n")
    for name in irNames:
        outTsv.write(name + "\n")

    out.close()
    outTsv.close()


def exportMockInvertons(irDb, outpath):
    out = open(outpath + '/invertedSeqs.fasta', 'w')
    for ir in irDb.IRs:
        out.write('>' + ir + '_f\n')
        out.write(irDb.IRs[ir].fSeq + '\n')
        out.write('>' + ir + '_r\n')
        out.write(irDb.IRs[ir].rSeq + '\n')

# parse ncbi gbk file
def parseGBK(path):
    genes = {}
    for record in SeqIO.parse(path, "genbank"):
        for seq_feature in record.features:
            # edge case where genes overlapping genome start index look like they cover the whole genome
            if seq_feature.type == 'CDS' and ( int(seq_feature.location.end) - int(seq_feature.location.start) ) < 100000:
                # overlap detection assumes start coordinate comes before stop coordinate
                if (int(seq_feature.location.end) - int(seq_feature.location.start) < 0):
                    warnings.warn("Warning...........Gene stop coordinate is less than gene start coordinate at: " + str(seq_feature.qualifiers['locus_tag'][0]))
                if int(seq_feature.location.strand) == -1:
                    strand = '-'
                else:
                    strand = '+'
                genes[seq_feature.qualifiers['locus_tag'][0]] = [record.id, int(seq_feature.location.start), int(seq_feature.location.end), strand]

    return genes

# Parse a prodigal gff file
def parseGFF(path):
    genes = {}

    fh = open (path)
    for line in fh:
        if not line.startswith('#'):
            sline = line.strip().split()
            if sline[2] == 'CDS':
                strand = sline[6]
                id = sline[0] + '_' + sline[8].split(';')[0].split('_')[1]
                genes[id] = [sline[0], int(sline[3]), int(sline[4]), strand]

    return genes

def findGeneOverlaps(genes, irDb):

    for ir in irDb.IRs:
        #clear gene overlaps from any previous results
        irDb.IRs[ir].geneOverlaps = []

        for gene in genes:
            if genes[gene][0] == irDb.IRs[ir].chr:
                if irDb.IRs[ir].leftStop >= genes[gene][1] and irDb.IRs[ir].rightStart <= genes[gene][2]:
                    irDb.IRs[ir].geneOverlaps.append([gene, "intragenic"])
                elif irDb.IRs[ir].leftStop <= genes[gene][1] and irDb.IRs[ir].rightStart >= genes[gene][1]:
                    if genes[gene][3] == "+":
                        irDb.IRs[ir].geneOverlaps.append([gene, "partial overlap, start"])
                    else:
                        irDb.IRs[ir].geneOverlaps.append([gene, "partial overlap, stop"])
                elif irDb.IRs[ir].leftStop <= genes[gene][2] and irDb.IRs[ir].rightStart >= genes[gene][2]:
                    if genes[gene][3] == "+":
                        irDb.IRs[ir].geneOverlaps.append([gene, "partial overlap, stop"])
                    else:
                        irDb.IRs[ir].geneOverlaps.append([gene, "partial overlap, start"])

    for ir in irDb.IRs:
        if len(irDb.IRs[ir].geneOverlaps) == 0:
            irDb.IRs[ir].geneOverlaps.append(["","intergenic"])

def exportGeneOverlaps(irDb, outpath):
    out = open(outpath + '/geneOverlaps.tsv', 'w')
    for ir in irDb.IRs:
        out.write(ir + "\t")
        for overlap in irDb.IRs[ir].geneOverlaps:
            out.write(overlap[0] + "," + overlap[1] + ";")
        out.write("\n")
    out.close()
