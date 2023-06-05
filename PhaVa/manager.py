
#!/usr/bin/env python3

'''
Manager- takes input from argparse and calls correct modules
'''


__author__ = "Patrick T West"
__license__ = "MIT"
__email__ = "ptwest10@gmail.com"
__status__ = "Development"

import argparse
import logging
import pickle
import os
import sys
from os.path import exists

import PhaVa
import PhaVa.summarize
import PhaVa.utils
import PhaVa.locate
import PhaVa.create
import PhaVa.ratio
import PhaVa.test
from PhaVa.fileManager import WorkDirectory

def version():
    x = ""
    #versionFile = open(os.path.join(PhaseFinderLR.__path__[0], 'VERSION'))
    #return versionFile.read().strip()

#VERSION = version()

class Manager():
    def __init__(self):
        self.logger = logging.getLogger()

    def main(self, args):
        ''' Parse user options and call the correct pipeline'''
        if (args.operation == "test"):
            self.test_installation(args)
        else:
            #setup working directory and logging
            wd = WorkDirectory(args.dir)
            if (args.log):
                logging.basicConfig(stream=sys.stdout, level=logging.INFO)
            else:
                logging.basicConfig(filename=args.dir + "/Phava.log", level=logging.INFO)

            if(args.operation == "locate"):
                logging.info("------Beginning IR locating operation------")
                self.ir_locate_operation(args)
                logging.info("------Finished IR locating operation------")

            if(args.operation == "create"):
                logging.info("------Beginning IR create operation------")
                self.ir_create_operation(args)
                logging.info("------Finished IR create operation------")

            if(args.operation == "ratio"):
                logging.info("------Beginning IR ratio operation------")
                self.ir_ratio_operation(args)
                logging.info("------Finished IR ratio operation------")

            if(args.operation =="variation_wf"):
                logging.info("------Beginning variation workflow operation------")
                self.variation_workflow_operation(args)
                logging.info("------Finished variation workflow operation------")

            if(args.operation =="summarize"):
                self.summarize_operation(args)



    def ir_locate_operation(self, args):
        # annotate contigs with prodigal, cmscan, and hmmscan
        logging.info("------Beginning IR search------")
        irDb = PhaVa.locate.main(args)
        logging.info("------Finished IR search------")
        self.pickleDb(args, irDb)

    def ir_create_operation(self, args):
        logging.info("------Beginning mock IR creation------")
        irDb = self.unpickleDb(args)
        if irDb is None:
            irDb = PhaVa.utils.IRDb()
        irDb = PhaVa.create.main(args, irDb)
        logging.info("------Finished mock IR creation------")
        self.pickleDb(args, irDb)

    def ir_ratio_operation(self, args):
        logging.info("------Beginning IR ratio calculation------")
        irDb = self.unpickleDb(args)
        if irDb is not None:
            irDb = PhaVa.ratio.main(args, irDb)
            logging.info("------Finished IR ratio calculation------")
            # pickling causes issues with running multiple ratio commands from the same directy, a common use
            #self.pickleDb(args, irDb)

    def summarize_operation(self, args):
        logging.info("------Beginning summarize operation------")
        irDb = self.unpickleDb(args)
        if irDb is not None:
            PhaVa.summarize.main(args, irDb)
            logging.info("------Finished summarize operation------")

    def test_installation(self, args):
        PhaVa.test.main(args)

    def variation_workflow_operation(self, args):
        self.ir_locate_operation(args)
        self.ir_create_operation(args)
        self.ir_ratio_operation(args)

    def pickleDb(self, args, irDb):
        pickleOut = open(args.dir + "/irDb.pickle", 'wb')
        pickle.dump(irDb, pickleOut)
        pickleOut.close()
        logging.info("------Finished pickling IR database------")

    def unpickleDb(self, args):
        if exists(args.dir + "/irDb.pickle"):
            pickled_db = open(args.dir + "/irDb.pickle", 'rb')
            irDb = pickle.load(pickled_db)
            pickled_db.close()
            logging.info("------Finished unpickling IR database------")
            return irDb
        else:
            logging.info("------No pickled IR database------")
            return None
