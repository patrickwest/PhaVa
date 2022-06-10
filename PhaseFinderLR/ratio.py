#!/usr/bin/env python

import logging
import PhaseFinderLR.utils
from multiprocessing.pool import ThreadPool
import subprocess
from functools import partial
from Bio.Seq import Seq
import os
import pysam
from os.path import exists

def main(args, irDb):

    # run minimap2
    # TODO: check if sam present in a smarter way
    if not exists(args.dir + '/intermediate/' + os.path.basename(args.fastq) + "_vs_" + irDb.genomeName + '.sam'):
        pool = ThreadPool(args.cpus)
        pool.map(partial(runMinimap, wd=args.dir, reads=args.fastq, threads=args.cpus, genomeName=irDb.genomeName), [0])
        # ensure processes are finished before moving on
        pool.close()
        pool.join()

        #index alignment so it can be used for pileup
        #pysam.sort("-o", args.dir + '/intermediate/' + os.path.basename(args.fastq) + "_vs_" + irDb.genomeName + '.bam', args.dir + '/intermediate/' + os.path.basename(args.fastq) + "_vs_" + irDb.genomeName + '.sam')
        #pysam.index(args.dir + '/intermediate/' + os.path.basename(args.fastq) + "_vs_" + irDb.genomeName + '.bam')

    # read in reads
    reads = parseSam(args.dir, args.fastq, irDb.genomeName, irDb, args.maxMismatch)

    # compute ratios
    reads = computeRatios(reads, irDb, os.path.basename(args.fastq), args.dir)

    # produce output
    reportInvertons(irDb, os.path.basename(args.fastq), args.dir)

    if (args.keepSam):
        reportMappingReads(args.dir, args.fastq, irDb.genomeName, reads)

    cleanup(args.dir, args.fastq, irDb.genomeName)

    return irDb

def runMinimap(index, wd, reads, genomeName, threads):
    command = "minimap2 -a --MD -o " + wd + '/intermediate/' + os.path.basename(reads) + "_vs_" + genomeName + '.sam -t ' + str(threads) + ' ' + wd + '/invertedSeqs.fasta ' + reads
    logging.info(command)
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return (out, err)

def parseSam(wd, reads, genomeName, irDb, maxMismatch):
    samfile = pysam.AlignmentFile(wd + '/intermediate/' + os.path.basename(reads) + "_vs_" + genomeName + '.sam', "r")

    reads = []
    for read in samfile:
        if (read.mapping_quality > 2):
            # mismatch density filter
            splithead = read.reference_name.rsplit('_', 1)
            ir = splithead[0]

            # get ir length
            if (irDb.IRs[ir].leftStart < irDb.IRs[ir].flankSize):
                start = irDb.IRs[ir].leftStart
            else:
                start = irDb.IRs[ir].flankSize
            stop = (irDb.IRs[ir].rightStop + irDb.IRs[ir].flankSize) - irDb.IRs[ir].leftStart
            cutoff = (stop - start) * maxMismatch

            num_mismatch = 0
            for pair in read.get_aligned_pairs(with_seq=True):
                if pair[1] != None and pair[1] > start and pair[1] < stop:
                    if pair[2].islower() or pair[0] == None:
                        num_mismatch += 1

            #print(str(cutoff) + "\t" + str(num_mismatch))

            if num_mismatch < cutoff:
                reads.append([splithead[0], splithead[1], read.reference_start, read.reference_end, read.query_name, read.mapping_quality])
    return reads

def computeRatios(reads, irDb, basename, wd):
    fdb = {}
    rdb = {}
    for read in reads:
        if read[1] == 'f':
            fdb[read[0] + '_' + read[4]] = read
        else:
            rdb[read[0] + '_' + read[4]] = read
            #print(read[0] + '_' + read[4])

    #rdb = filterMismatch(wd, rdb, basename, irDb)

    final_filt_reads = {}
    ratios = {}

    # clear read counts for this sample for (re)counting
    for ir in irDb.IRs:
        irDb.IRs[ir].forwardReads = {}
        irDb.IRs[ir].reverseReads = {}
        irDb.IRs[ir].ratio = {}
        irDb.IRs[ir].forwardReads[basename] = 0
        irDb.IRs[ir].reverseReads[basename] = 0
        irDb.IRs[ir].ratio[basename] = 0

    for key in fdb:
        if key not in rdb and fdb[key][0] in irDb.IRs:
            # filter out reads that dont map into flanking sequence
            length =  (irDb.IRs[fdb[key][0]].rightStop + irDb.IRs[fdb[key][0]].flankSize) - (irDb.IRs[fdb[key][0]].leftStart - irDb.IRs[fdb[key][0]].flankSize)
            if (int(fdb[key][2]) < irDb.IRs[fdb[key][0]].flankSize - 20 and int(fdb[key][3]) > (length - (irDb.IRs[fdb[key][0]].flankSize + 20))):
                irDb.IRs[fdb[key][0]].forwardReads[basename] += 1
                final_filt_reads[fdb[key][0] + '_f_' + str(fdb[key][4]) + '_' + str(fdb[key][2]) + '_' + str(fdb[key][3])] = fdb[key]

    for key in rdb:
        if key not in fdb and rdb[key][0] in irDb.IRs:
            # filter out reads that dont map into flanking sequence
            length =  (irDb.IRs[rdb[key][0]].rightStop + irDb.IRs[rdb[key][0]].flankSize) - (irDb.IRs[rdb[key][0]].leftStart - irDb.IRs[rdb[key][0]].flankSize)
            if (int(rdb[key][2]) < irDb.IRs[rdb[key][0]].flankSize - 20 and int(rdb[key][3]) > (length - (irDb.IRs[rdb[key][0]].flankSize + 20))):
                irDb.IRs[rdb[key][0]].reverseReads[basename] += 1
                final_filt_reads[rdb[key][0] + '_r_' + str(rdb[key][4]) + '_' + str(rdb[key][2]) + '_' + str(rdb[key][3])] = rdb[key]

    for ir in irDb.IRs:
        if irDb.IRs[ir].forwardReads[basename] == 0:
            irDb.IRs[ir].ratio[basename] = 10
        else:
            irDb.IRs[ir].ratio[basename] = irDb.IRs[ir].reverseReads[basename] / (irDb.IRs[ir].forwardReads[basename] + irDb.IRs[ir].reverseReads[basename])

    return final_filt_reads


def reportInvertons(irDb, basename, wd):

    out_fh = open(wd + '/' + basename + "_vs_" + irDb.genomeName + '_ratio.tsv', 'w')
    for ir in irDb.IRs:
        if irDb.IRs[ir].reverseReads[basename] >= 1:
            print(ir + "\t" + str(irDb.IRs[ir].forwardReads[basename]) + "\t" + str(irDb.IRs[ir].reverseReads[basename]) + "\t" + str(irDb.IRs[ir].ratio[basename]) + "\t" + basename)
            out_fh.write(ir + "\t" + str(irDb.IRs[ir].forwardReads[basename]) + "\t" + str(irDb.IRs[ir].reverseReads[basename]) + "\t" + str(irDb.IRs[ir].ratio[basename]) + "\t" + basename + "\n")

    out_fh.close()

def reportMappingReads(wd, readsFn, genomeName, reads):
    # open sam
    bam = pysam.AlignmentFile(wd + '/intermediate/' + os.path.basename(readsFn) + "_vs_" + genomeName + '.sam', "r")
    # open out same
    filtreads = pysam.AlignmentFile(wd + "/" + os.path.basename(readsFn) + "_vs_" + genomeName + '.filtered.bam', "wb", template=bam)

    # write reads that pass filters
    for read in bam.fetch():
        if (str(read.reference_name) + '_' + str(read.query_name) + '_' + str(read.reference_start) + '_' + str(read.reference_end)) in reads:
            filtreads.write(read)

def cleanup(wd, readsFn, genomeName):
    os.remove(wd + '/intermediate/' + os.path.basename(readsFn) + "_vs_" + genomeName + '.sam')
