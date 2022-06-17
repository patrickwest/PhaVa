#!/usr/bin/python

import pytest
import subprocess
import os
import filecmp
import shlex
from pathlib import Path

phaseFinderLR_rootDir = Path(__file__).parent.absolute()

phaseFinderWorkingDir = os.path.dirname(phaseFinderLR_rootDir) + "/tests/test_outputs"
refGenomePath = os.path.dirname(phaseFinderLR_rootDir) + "/tests/GCA_000011065.1_ASM1106v1_genomic.fna"
readsDirectory = os.path.dirname(phaseFinderLR_rootDir) + "/tests/simulated_reads"
comparisonOutputFile = os.path.dirname(phaseFinderLR_rootDir) + "/tests/test_outputs/comparison_inverton_hits.tsv"
phaseFinderLR = os.path.dirname(phaseFinderLR_rootDir) + "/bin/PhaseFinderLR"
filterFile = os.path.dirname(phaseFinderLR_rootDir) + "/tests/runFilter.sh"

def runPhaseFinderLRPipeline():
	# Running the phaseFinder locate process
	locate_command = 'python ' + phaseFinderLR + '  locate -d ' + phaseFinderWorkingDir + ' -t 16 -i ' + refGenomePath
	formatted_locate_command = shlex.split(locate_command)
	locate_process = subprocess.Popen(formatted_locate_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	locate_out, locate_err = locate_process.communicate()
	print('PhaseFinder Locate:')
	print(locate_out)
	print(locate_err)
	
	# Running the phaseFinder create process
	create_command = 'python ' + phaseFinderLR + '  create -d ' + phaseFinderWorkingDir + ' -t 16 -i ' + refGenomePath
	formatted_create_command = shlex.split(create_command)
	create_process = subprocess.Popen(formatted_create_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	create_out, create_err = create_process.communicate()
	print('PhaseFinder Create:')
	print(create_out)
	print(create_err)

	# Running the phaseFinder ratio process over all the reads in a given directory
	for filename in os.listdir(readsDirectory):
		f = os.path.join(readsDirectory, filename)
		if os.path.isfile(f):
			ratio_command = 'python ' + phaseFinderLR + '  ratio -d ' + phaseFinderWorkingDir + ' -t 16 -r ' + f
			formatted_ratio_command = shlex.split(ratio_command)
			ratio_process = subprocess.Popen(formatted_ratio_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
			ratio_out, ratio_err = ratio_process.communicate()
			print('PhaseFinder Ratio for ' + f + ':')
			print(ratio_out)
			print(ratio_err)

	# Combining all .tsv output files into a single .tsv file
	# Can be done by running the filter bash script previously written
	# Keeps only the invertons w/ >= 4 reverse hits and rev/all ratio > 0.009
	filter_command = 'bash ' + filterFile
	formatted_filter_command = shlex.split(filter_command)
	filter_process = subprocess.call(formatted_filter_command)

def test_comparePhaseFinderOutputs():
	runPhaseFinderLRPipeline()
	# Want to compare the just-generated output file to the previously manually-constructed output file
	assert filecmp.cmp('test_outputs/filtered_inverton_hits.tsv',  comparisonOutputFile, shallow=False) 
