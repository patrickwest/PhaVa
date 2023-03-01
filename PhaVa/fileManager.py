#!/usr/bin/env python3

import os
import argparse

class WorkDirectory(object):

    def __init__(self, location):
        self.location = os.path.abspath(location)
        self.initFileStructure()

    def initFileStructure(self):
        '''
        Make the top level file structure
        '''
        location = self.location

        if not os.path.exists(location):
            os.makedirs(location)

        loc = location + '/' + 'intermediate'
        if not os.path.exists(loc):
            os.makedirs(loc)
