#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from multiprocessing.pool import ThreadPool
import subprocess
from functools import partial
import logging
import re
import os
import PhaVa.utils

def main(args):
    processedLog = []
    irDb = PhaVa.utils.IRDb()

    pool = ThreadPool(args.cpus)
    pool.map(partial(runEinverted, contigsPath=args.fasta, wd=args.dir), [0,1])

    # ensure processes are finished before moving on
    pool.close()
    pool.join()

    irDb.genome = parseGenome(args.fasta)
    irDb.genomeName = os.path.basename(args.fasta)
    irDb.IRs = parseEinverted(args.dir)
    irDb.IRs = parseIRSeqs(irDb.IRs, irDb.genome)
    irDb.IRs = applyFilters(irDb.IRs)

    exportIRs(irDb.IRs, args.dir + '/data')

    return irDb


def runEinverted(index, contigsPath, wd):
    if index == 0:
        (out,err) = runEinverted51(contigsPath, wd)
        logging.info(err.strip())
    elif index == 1:
        (out,err) = runEinverted75(contigsPath, wd)
        logging.info(err.strip())

def runEinverted51(contigsPath, wd):
    basename = os.path.basename(contigsPath)
    outDir = wd + "/intermediate/"

    command = "einverted -maxrepeat 750  -gap 100 -threshold 51 -match 5 -mismatch -9 -outfile " + outDir + "einverted.51.outfile -outseq " + outDir + "einverted.51.outseq -sequence " + contigsPath
    logging.info(command)
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return (out, err)

def runEinverted75(contigsPath, wd):
    basename = os.path.basename(contigsPath)
    outDir = wd + "/intermediate/"

    command = "einverted -maxrepeat 750  -gap 100 -threshold 75 -match 5 -mismatch -15 -outfile " + outDir + "einverted.75.outfile -outseq " + outDir + "einverted.75.outseq -sequence " + contigsPath
    logging.info(command)
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return (out, err)

def parseGenome(fasta_fn):
    seqs = {}
    for record in SeqIO.parse(fasta_fn, "fasta"):
        seqs[record.id] = str(record.seq)
    return seqs

def parseEinverted(wd):
    #i51 = wd + "/intermediate/einverted.51.outseq"
    #i75 = wd + "/intermediate/einverted.75.outseq"
    i51 = wd + "/intermediate/einverted.51.outfile"
    i75 = wd + "/intermediate/einverted.75.outfile"

    #i51IRs = parseEinvertedSeq(i51)
    #i75IRs = parseEinvertedSeq(i75)
    i51IRs = parseEinvertedOutput(i51)
    i75IRs = parseEinvertedOutput(i75)

    # only include non-overlapping IRs
    # TODO: think about optimizing with sort
    IRs = {}
    for ir51 in i51IRs:
        overlap = False
        for ir75 in i75IRs:
            # check if on same chromosome
            if ir51[0] == ir75[0]:
                # check if overlapping:
                if (ir51[4] >= ir75[1] and ir75[4] >= ir51[1]):
                    overlap = True
        if not overlap:
            #IRs.append(ir51)
            IRs[ir51[0] + ':' + str(ir51[1]) + '-' + str(ir51[2]) + '-' + str(ir51[3]) + '-' + str(ir51[4])] = \
                PhaVa.utils.IR(ir51[0], ir51[1], ir51[2], ir51[3], ir51[4])

    for ir75 in i75IRs:
        #IRs.append(ir75)
        IRs[ir75[0] + ':' + str(ir75[1]) + '-' + str(ir75[2]) + '-' + str(ir75[3]) + '-' + str(ir75[4])] = \
            PhaVa.utils.IR(ir75[0], ir75[1], ir75[2], ir75[3], ir75[4])

    return IRs

# parse IR coordinates from einverted outseq file
def parseEinvertedSeq(path):
    IRs = []
    tempIR = []

    i = 1
    for record in SeqIO.parse(path, "fasta"):
        # chromosome and start/stop located in fasta header. Always in pairs, so group together every two fasta entries
        temp = record.id.split('_')
        if (i%2 != 0):
            tempIR.extend([temp[0], int(temp[1])-1, int(temp[2])])
        else:
            tempIR.extend([int(temp[1])-1, int(temp[2])])
            IRs.append(tempIR)
            tempIR = []
        i += 1
    return IRs

# parse IR coordinates from einverted outfile file
def parseEinvertedOutput(path):
    IRs = []

    fh = open(path)
    lines = []
    for line in fh:
        lines.append(line)
        if len(lines) >= 5:
            tempIR = []
            tempIR.append(lines[1].split(':')[0])
            fCoor = lines[2].split()
            tempIR.append(int(fCoor[0])-1)
            tempIR.append(int(fCoor[2]))
            rCoor = lines[4].split()
            tempIR.append(int(rCoor[2])-1)
            tempIR.append(int(rCoor[0]))
            IRs.append(tempIR)
            lines = []
    return IRs

# pull out inverted repeat seqeucnes and sequence located between them
def parseIRSeqs(IRs, seqs):
    for ir in IRs:
        IRs[ir].leftSeq = seqs[IRs[ir].chr] [IRs[ir].leftStart:IRs[ir].leftStop]
        IRs[ir].rightSeq = seqs[IRs[ir].chr] [IRs[ir].rightStart:IRs[ir].rightStop]
        IRs[ir].middleSeq = seqs[IRs[ir].chr] [IRs[ir].leftStop:IRs[ir].rightStart]

    return IRs

def applyFilters(IRs):
    toBeDel = []
    for ir in IRs:
        lgc = gc_fraction(IRs[ir].leftSeq, ambiguous='ignore') * 100
        rgc = gc_fraction(IRs[ir].rightSeq, ambiguous='ignore') * 100
        # GC filter
        #if (lgc <= 15 or rgc <= 15 or lgc >= 85 or rgc >= 85):
        #    IRs.remove(ir)
        # homopolymer filter
        #elif len(re.findall(r'([ACGT])\1{4,}', ir.leftSeq)) > 0 and len(re.findall(r'([ACGT])\1{4,}', ir.rightSeq)) > 0:
        #    IRs.remove(ir)
        # minimum distance between IRs filter
        if (len(IRs[ir].middleSeq) < 30):
            toBeDel.append(ir)

    for key in toBeDel:
        del IRs[key]

    return IRs

def exportIRs(IRs, outpath):
    out = open(outpath + '/IRs.tsv', 'w')
    out.write("IR_Chr\tLeftIRStart\tLeftIRStop\tRightIRStart\tRightIRStop\tLeftIRSequence\tInvertibleSequence\tRightIRSequence\n")
    for ir in IRs:
        out.write(IRs[ir].chr + '\t' + \
        str(IRs[ir].leftStart) + '\t' + \
        str(IRs[ir].leftStop) + '\t' + \
        str(IRs[ir].rightStart) + '\t' + \
        str(IRs[ir].rightStop) + '\t' + \
        IRs[ir].leftSeq + '\t' + \
        IRs[ir].middleSeq + '\t' + \
        IRs[ir].rightSeq + '\n')

    out.close()
