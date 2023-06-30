#!/usr/bin/env python

import logging
from multiprocessing.pool import ThreadPool
import subprocess
from functools import partial
import os
import pysam
from os.path import exists
from collections import defaultdict


def main(args, irDb):
    samfile = args.dir + '/intermediate/' + os.path.basename(args.fastq) + "_vs_" + irDb.genomeName + '.sam'
    if args.short_reads:
        if not exists(samfile):
            checkBowtieIdx(args.dir)
            pool = ThreadPool(args.cpus)
            pool.map(partial(runBowtie, wd=args.dir, reads_forward=args.fastq, reads_reverse=args.fastq2,
                             threads=args.cpus, genomeName=irDb.genomeName), [0])
            # ensure processes are finished before moving on
            pool.close()
            pool.join()
        reads = parseSamPaired(args.dir, args.fastq, irDb.genomeName, args.maxMismatch)
        ratios = computeRatiosPaired(reads)
        reportInvertonsPaired(ratios, os.path.basename(args.fastq), args.dir, irDb)
    else:
        # run minimap2
        # TODO: check if sam present in a smarter way
        if not exists(samfile):
            pool = ThreadPool(args.cpus)
            pool.map(partial(runMinimap, wd=args.dir, reads=args.fastq,
                             threads=args.cpus, genomeName=irDb.genomeName), [0])
            # ensure processes are finished before moving on
            pool.close()
            pool.join()

        # read in reads
        reads = parseSam(args.dir, args.fastq, irDb.genomeName, irDb, args.maxMismatch)

        # compute ratios
        reads = computeRatios(reads, irDb, os.path.basename(args.fastq), args.dir)

        # produce output
        reportInvertons(irDb, os.path.basename(args.fastq), args.dir, args.reportAll)

        if args.keepSam:
            reportMappingReads(args.dir, args.fastq, irDb.genomeName, reads)

    cleanup(args.dir, args.fastq, irDb.genomeName)

    return irDb


def runMinimap(index, wd, reads, genomeName, threads):
    command = "minimap2 -a --MD -o " + wd + '/intermediate/' + os.path.basename(reads) + "_vs_" + \
              genomeName + '.sam -t ' + str(threads) + ' ' + wd + '/invertedSeqs.fasta ' + reads
    logging.info(command)
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return out, err


def checkBowtieIdx(wd):
    ref_data = wd + "/invertedSeqs.fasta"
    if not exists(ref_data + ".1.bt2"):
        command = "bowtie2-build " + ref_data + " " + ref_data
        logging.info("Building bowtie2 index")
        logging.info(command)
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        streamdata = p.communicate()
        logging.info("Finished building bowtie2 index")
        return streamdata[0]
    else:
        return True


def runBowtie(index, wd, reads_forward, reads_reverse, genomeName, threads):
    # run bowtie
    samfile = wd + '/intermediate/' + os.path.basename(reads_forward) + '_vs_' + genomeName + ".sam"
    command = "bowtie2 -p " + str(threads) + " -x " + wd + "/invertedSeqs.fasta" + " -S " + samfile + \
              " -1 " + reads_forward + " -2 " + reads_reverse
    logging.info(command)
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return out, err


def parseSam(wd, reads, genomeName, irDb, maxMismatch):
    samfile = pysam.AlignmentFile(wd + '/intermediate/' + os.path.basename(reads) + "_vs_" + genomeName + '.sam', "r")

    reads = []
    for read in samfile:
        if (read.mapping_quality > 2):
            # mismatch density filter
            splithead = read.reference_name.rsplit('_', 1)
            ir = splithead[0]

            # get ir length
            if (int(irDb.IRs[ir].leftStart) < int(irDb.IRs[ir].flankSize)):
                # start = irDb.IRs[ir].leftStart
                start = irDb.IRs[ir].leftStop
            else:
                # start = irDb.IRs[ir].flankSize
                start = irDb.IRs[ir].flankSize + (irDb.IRs[ir].leftStop - irDb.IRs[ir].leftStart)
            # stop = (irDb.IRs[ir].rightStop + start) - irDb.IRs[ir].leftStart
            stop = (irDb.IRs[ir].rightStart + start) - irDb.IRs[ir].leftStop
            cutoff = (stop - start) * maxMismatch

            num_mismatch = 0
            num_flank_mismatch = 0
            flank_length = 0
            for pair in read.get_aligned_pairs(with_seq=True):
                if pair[1] != None:
                    if pair[1] > start and pair[1] < stop:
                        if pair[2].islower() or pair[0] == None:
                            num_mismatch += 1
                    else:
                        flank_length += 1
                        if pair[2].islower() or pair[0] == None:
                            num_flank_mismatch += 1

            # also check that flanking region doesnt have a high mismatch rate, but consider them separately to avoid masking a signal form one or the other
            flank_cutoff = flank_length * maxMismatch

            if num_mismatch < cutoff and num_flank_mismatch < flank_cutoff:
                reads.append([splithead[0], splithead[1], read.reference_start, read.reference_end,
                              read.query_name, read.mapping_quality])

    return reads


def check_aln_paired(read, maxMismatch):
    if read.mapping_quality >= 1:
        splithead = read.reference_name.rsplit('_', 1)
        ir = splithead[0]
        num_mismatch = 0
        num_flank_mismatch = 0
        flank_length = 0
        mid_length = 0

        ir_length = ir.split(':')[1].split('-')
        leftStart = int(ir_length[0])
        leftStop = int(ir_length[1])
        rightStart = int(ir_length[2])
        rightStop = int(ir_length[3])
        flankSize = 1000
        if int(leftStart) < int(flankSize):
            # start = irDb.IRs[ir].leftStart
            start = leftStop
        else:
            # start = irDb.IRs[ir].flankSize
            start = flankSize + (leftStop - leftStart)
        # stop = (irDb.IRs[ir].rightStop + start) - irDb.IRs[ir].leftStart
        stop = (rightStart + start) - leftStop
        cutoff = (stop - start) * maxMismatch
        for pair in read.get_aligned_pairs(with_seq=True):
            if pair[1] is not None:
                if pair[1] > start and pair[1] < stop:
                    mid_length += 1
                    if pair[2].islower() or pair[0] is None:
                        num_mismatch += 1
                else:
                    flank_length += 1
                    if pair[2].islower() or pair[0] is None:
                        num_flank_mismatch += 1

        flank_cutoff = flank_length * maxMismatch
        # check how many mismatched are in the mid-IR matching region
        if mid_length > 0:
            mid_mismatch_ratio = num_mismatch/mid_length
        else:
            mid_mismatch_ratio = 0

        if flank_length > 0:
            flank_mismatch_ratio = num_flank_mismatch/flank_length
        else:
            flank_mismatch_ratio = 0

        if num_mismatch < cutoff and num_flank_mismatch <= flank_cutoff:

            return [read.reference_name,
                    flank_length, mid_length, flank_mismatch_ratio, mid_mismatch_ratio,
                    read.reference_start, read.reference_end, read]
        else:
            return None
    else:
        return None


def parseSamPaired(wd, reads, genomeName, maxMismatch):
    samfile = pysam.AlignmentFile(wd + '/intermediate/' + os.path.basename(reads) + "_vs_" + genomeName + '.sam', "r")
    read_dict = defaultdict(lambda: [None, None])
    for read in samfile:
        qname = read.query_name
        if read.is_read1:
            read_dict[qname][0] = check_aln_paired(read, maxMismatch)
        else:
            read_dict[qname][1] = check_aln_paired(read, maxMismatch)
    return read_dict


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

    # TODO: combine duplicate code into one method
    for key in fdb:
        if key not in rdb and fdb[key][0] in irDb.IRs:
            # filter out reads that dont map into flanking sequence
            length = (irDb.IRs[fdb[key][0]].rightStop - irDb.IRs[fdb[key][0]].leftStart)
            # get start point of ir (its the raw start if start < flanksize)
            if (irDb.IRs[fdb[key][0]].leftStart < irDb.IRs[fdb[key][0]].flankSize):
                start = irDb.IRs[fdb[key][0]].leftStart
            else:
                start = irDb.IRs[fdb[key][0]].flankSize

            if (int(fdb[key][2]) < start - 20 and int(fdb[key][3]) > (length + (start + 20))):
                irDb.IRs[fdb[key][0]].forwardReads[basename] += 1
                final_filt_reads[fdb[key][0] + '_f_' + str(fdb[key][4]) + '_' + str(fdb[key][2]) + '_' + str(fdb[key][3])] = fdb[key]

    for key in rdb:
        if key not in fdb and rdb[key][0] in irDb.IRs:
            # filter out reads that dont map into flanking sequence
            length = (irDb.IRs[rdb[key][0]].rightStop - irDb.IRs[rdb[key][0]].leftStart)
            # get start point of ir (its the raw start if start < flanksize)
            if (irDb.IRs[rdb[key][0]].leftStart < irDb.IRs[rdb[key][0]].flankSize):
                start = irDb.IRs[rdb[key][0]].leftStart
            else:
                start = irDb.IRs[rdb[key][0]].flankSize

            if (int(rdb[key][2]) < start - 20 and int(rdb[key][3]) > (length + (start + 20))):
                irDb.IRs[rdb[key][0]].reverseReads[basename] += 1
                final_filt_reads[rdb[key][0] + '_r_' + str(rdb[key][4]) + '_' + str(rdb[key][2]) + '_' + str(rdb[key][3])] = rdb[key]

    for ir in irDb.IRs:
        if irDb.IRs[ir].forwardReads[basename] == 0:
            irDb.IRs[ir].ratio[basename] = 1
        else:
            irDb.IRs[ir].ratio[basename] = irDb.IRs[ir].reverseReads[basename] / (irDb.IRs[ir].forwardReads[basename] + irDb.IRs[ir].reverseReads[basename])

    return final_filt_reads


def computeRatiosPaired(read_dict):
    ir_counts = defaultdict(lambda: [0, 0, 0, 0])  # forward, forward_IR, reverse, reverse_IR
    for qname in read_dict.keys():
        reads = [x for x in read_dict[qname] if x is not None]
        if len(reads) == 0:
            continue
        # check for those overlapping the IR

        # also remove read pairs if the IR_matching region is too many mismatches
        # need to check this in the future!
        hit_ir = sum([x[2] if x[4] < 0.25 else 0 for x in reads])
        hit_flank = sum([x[1] if x[3] < 0.25 else 0 for x in reads])
        # remove read pairs that are either fully in the flanking region
        # or fully within the IR (they give no information)
        if hit_ir <= 20 or hit_flank <= 20:
            continue

        if len(reads) == 1:
            # for single reads (they are by definition overlapping, add 0.5 to the respective IR)
            ref = reads[0][0].rsplit('_', 1)

            if ref[1] == 'f':
                idx = 0
            else:
                idx = 2
            ir_counts[ref[0]][idx] = ir_counts[ref[0]][idx] + 0.5
            ir_counts[ref[0]][idx + 1] = ir_counts[ref[0]][idx + 1] + 0.5
        if len(reads) == 2:
            ref = list(set([x[0].rsplit('_', 1)[0] for x in reads]))
            # check that the reference is the same (Otherwise sth is wrong!)
            if len(ref) != 1:
                continue

            direction = [x[0].rsplit('_', 1)[1] for x in reads]

            # if the reference and direction is the same, check if the orientation is right
            if len(set(direction)) == 1:
                read1_direction = reads[0][7].is_reverse
                read2_direction = reads[1][7].is_reverse
                correct_direction = (read1_direction and not read2_direction) or (
                            read2_direction and not read1_direction)
                if correct_direction:
                    if direction[0] == 'f':
                        idx = 0
                    else:
                        idx = 2
                else:
                    if direction[0] == 'f':
                        idx = 2
                    else:
                        idx = 0
                ir_counts[ref[0]][idx] = ir_counts[list(ref)[0]][idx] + 1
                # check if a read overlaps the IR border
                if (reads[0][1] > 0 and reads[0][2] > 0) or (reads[1][1] > 0 and reads[1][2] > 0):
                    ir_counts[ref[0]][idx + 1] = ir_counts[ref[0]][idx + 1] + 1
            else:
                # if the reference is the same, but direction is different, what to do?
                # if one of the reads spans the IR, take this orientation
                read1_span = (reads[0][1] > 20 and reads[0][2] > 20)
                read2_span = (reads[1][1] > 20 and reads[1][2] > 20)
                span = False
                if read1_span or read2_span:
                    span = True

                    # if not, take the one mapping in the IR as reference,
                    #   pretend the other mapped on the same orientation flank (which it would)
                    #   check orientation and take the one that make sense
                ir_match = [x[2] for x in reads]
                which_read = ir_match.index(max(ir_match))
                read1_direction = reads[0][7].is_reverse
                read2_direction = reads[1][7].is_reverse
                correct_direction = (read1_direction and not read2_direction) or (
                        read2_direction and not read1_direction)
                if correct_direction:
                    if direction[which_read] == 'f':
                        idx = 0
                    else:
                        idx = 2
                else:
                    if direction[which_read] == 'f':
                        idx = 2
                    else:
                        idx = 0
                ir_counts[ref[0]][idx] = ir_counts[list(ref)[0]][idx] + 1
                # check if a read overlaps the IR border
                if span:
                    ir_counts[ref[0]][idx + 1] = ir_counts[ref[0]][idx + 1] + 1

    return ir_counts


def reportInvertons(irDb, basename, wd, reportAll):

    out_fh = open(wd + '/' + basename + "_vs_" + irDb.genomeName + '_ratio.tsv', 'w')
    out_fh.write("Inverton\tforwardReads\treverseReads\tratio\tread_file\n")
    for ir in irDb.IRs:
        if irDb.IRs[ir].reverseReads[basename] >= 1 or reportAll:
            outstring = ir + "\t" + str(irDb.IRs[ir].forwardReads[basename]) + "\t" + \
                  str(irDb.IRs[ir].reverseReads[basename]) + "\t" + \
                  str(irDb.IRs[ir].ratio[basename]) + "\t" + basename
            print(outstring)
            out_fh.write(outstring + "\n")

    out_fh.close()


def reportInvertonsPaired(ratios, basename, wd, irDb):
    out_fh = open(wd + '/' + basename + "_vs_" + irDb.genomeName + '_ratio.tsv', 'w')
    out_fh.write('InvertonID\tforwardReadPairs\treverseReadPairs\tforwardReadPairsIR'+
                 '\treverseReadPairsIR\tratio\tratioIR\treadfile\n')
    for inv in sorted(ratios.keys()):
        outstring = inv + '\t' + str(ratios[inv][0])  + '\t' + str(ratios[inv][2])  + '\t' + \
            str(ratios[inv][1])  + '\t' + str(ratios[inv][3])  + '\t' + \
            str(ratios[inv][2]/(ratios[inv][2] + ratios[inv][0])) + '\t' + \
            str(ratios[inv][3]/(ratios[inv][3] + ratios[inv][1])) + "\t" + basename
        print(outstring)
        out_fh.write(outstring + '\n')

    out_fh.close()


def reportMappingReads(wd, readsFn, genomeName, reads):
    # open sam
    bam = pysam.AlignmentFile(wd + '/intermediate/' + os.path.basename(readsFn) + "_vs_" + genomeName + '.sam', "r")
    # open out same
    filtreads = pysam.AlignmentFile(wd + "/" + os.path.basename(readsFn) + "_vs_" +
                                    genomeName + '.filtered.bam', "wb", template=bam)

    # write reads that pass filters
    for read in bam.fetch():
        if (str(read.reference_name) + '_' + str(read.query_name) + '_' + str(
                read.reference_start) + '_' + str(read.reference_end)) in reads:
            filtreads.write(read)


def cleanup(wd, readsFn, genomeName):
    os.remove(wd + '/intermediate/' + os.path.basename(readsFn) + "_vs_" + genomeName + '.sam')
