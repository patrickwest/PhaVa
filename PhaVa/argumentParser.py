#!/usr/bin/env python

'''
PhaVa- parse command-line arguemnts
'''

__author__ = "Patrick T West"
__license__ = "MIT"
__email__ = "ptwest10@gmail.com"
__status__ = "Development"

import argparse
import os
import sys

import PhaVa
from PhaVa.manager import Manager


def version():
    versionFile = open(os.path.join(PhaVa.__path__[0], 'VERSION'))
    return versionFile.read().strip()

VERSION = version()


class SmartFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


def printHelp():
    print('')
    print('                ,,,;;;| PhaVa v' + VERSION + ' |;;;,,,''')
    print('''\
  Patrick T West. MIT License. Bhatt Lab, Stanford. 2022 (last updated 2022)
  Choose one of the operations below for more detailed help.

  Example: phava variation_wf -h
  Commands:
    variation_wf            -> Combine locate, create, and ratio steps 
                               into one workflow for identifying phase variation

    locate                  -> Identify regions flanked by invertable repeats
    create                  -> Create in silico flipped versions of invertons
    ratio                   -> Map reads and compute ratio of forward vs reverse 
                               orientation invertons

    summarize               -> Report general statistics (requires running 
                               locate/create or variation_wf beforehand)
    
    test                    -> Test your installation
    ''')


def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=SmartFormatter)
    subparsers = parser.add_subparsers(help='Desired operation', dest='operation')

    parent_parser = argparse.ArgumentParser(add_help=False)

    ReqFlags = parent_parser.add_argument_group('REQUIRED PARAMETERS')
    ReqFlags.add_argument("-d","--dir", required=True, help="R|Directory where data and output are stored\
    \n*** USE THE SAME WORK DIRECTORY FOR ALL PHASEFINDERLR OPERATIONS ***")

    SysFlags = parent_parser.add_argument_group('SYSTEM PARAMETERS')
    SysFlags.add_argument("-t", "--cpus", help="Number of threads to use", default=1, type=int)
    SysFlags.add_argument("-l", "--log",
                          help="Should the logging info be output to stdout? Otherwise, it will be written to 'PhaVa.log'",
                          action='store_true')

    locate_parent = argparse.ArgumentParser(add_help=False)
    LocateFlags = locate_parent.add_argument_group('LOCATE PARAMETERS')
    LocateFlags.add_argument("-i", "--fasta", help="Name of input assembly file to be searched",
                         type=str)

    create_parent = argparse.ArgumentParser(add_help=False)
    CreateFlags = create_parent.add_argument_group('CREATE PARAMETERS')
    CreateFlags.add_argument("-f", "--flankSize", help="Size flanking size to include on either side of invertable regions (in bps)", default=1000,
                         type=int)
    #CreateFlags.add_argument("--geneOverlap", help="Report gene/inverton overlaps. Requires a list of features in genbank format from --genes"
    #                     action='store_true')
    CreateFlags.add_argument("--genes", help="List of gene features in ncbi genbank format, for detecting gene/inverton overlaps",
                         type=str)
    geneFormats = ["gff","gbff"]
    CreateFlags.add_argument("--genesFormat", help="File format of the list of gene features. Gff must be in prodigal gff format",
                         choices=geneFormats, default="gbff")
    CreateFlags.add_argument("--mockGenome", help="Create a mock genome where all putative IRs are flipped to opposite of the reference orientation",
                         action='store_true')
    CreateFlags.add_argument("--mockNumber", help="If creating a mockGenome, the number of invertons to invert. A value of 0 inverts all predicted inverton locations",
                         default=0, type=int)

    createReq_parent = argparse.ArgumentParser(add_help=False)
    CreateReqFlags = createReq_parent.add_argument_group('CREATE SPECIFIC PARAMETERS')
    CreateReqFlags.add_argument("-i", "--fasta", help="Name of input assembly file to be searched",
                         type=str)
    CreateReqFlags.add_argument("--irs", help="Table of identified invertable repeats (eg. if locate command was never run)",
                         type=str)

    ratio_parent = argparse.ArgumentParser(add_help=False)
    RatioFlags = ratio_parent.add_argument_group('RATIO PARAMETERS')
    RatioFlags.add_argument("-r", "--fastq", help="Name of the reads file to be used for mapping",
                         type=str)
    RatioFlags.add_argument("-m", "--maxMismatch", help="Maximum proportion of inverton sequence that can be mismatch before a read is removed", default=0.15,
                         type=float)
    RatioFlags.add_argument("--keepSam", help="Keep the sam file from the mapping",
                         action='store_true')
    RatioFlags.add_argument("--reportAll", help="Report mapping results for all putative invertons, regardless of outcome",
                         action='store_true')

    ratioReq_parent = argparse.ArgumentParser(add_help=False)
    RatioReqFlags = ratioReq_parent.add_argument_group('RATIO SPECIFIC PARAMETERS')
    RatioReqFlags.add_argument("--inv", help="Fasta file of forward and inverted repeats (ie. generated from Create command)",
                         type=str)

    # create subparsers
    test_parser = argparse.ArgumentParser(add_help=False)
    locate_parser = subparsers.add_parser("locate", formatter_class=SmartFormatter,\
                                        parents=[parent_parser, locate_parent])
    create_parser = subparsers.add_parser("create", formatter_class=SmartFormatter,\
                                        parents=[parent_parser, create_parent, createReq_parent])
    ratio_parser = subparsers.add_parser("ratio", formatter_class=SmartFormatter,\
                                        parents=[parent_parser, ratio_parent, ratioReq_parent])
    variation_wf_parser = subparsers.add_parser("variation_wf", formatter_class=SmartFormatter,\
                                        parents=[parent_parser, locate_parent, create_parent, ratio_parent])
    summarize_parser = subparsers.add_parser("summarize", formatter_class=SmartFormatter,\
                                        parents=[parent_parser])
    test_parser = subparsers.add_parser("test", formatter_class=SmartFormatter,
                                        parents=[test_parser])

    # Handle the situation where the user wants the raw help
    if len(args) == 0 or args[0] == '-h' or args[0] == '--help':
        printHelp()
        sys.exit(0)
    else:
        return parser.parse_args(args)
