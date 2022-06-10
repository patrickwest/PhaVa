#!/usr/bin/env python

def main(args, irDb):

    reportRepeats(irDb, args.dir)





def reportRepeats(irDb, wd):

    outFh = open(wd + "/repeatSequences.tsv", 'w')

    outFh.write("IR\tLeft Repeat\tRight Repeat\n")
    for ir in irDb.IRs:
        outFh.write(ir + "\t" + irDb.IRs[ir].leftSeq + "\t" + irDb.IRs[ir].rightSeq + "\n")
