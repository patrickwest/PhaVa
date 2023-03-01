#!/usr/bin/python

import pytest
import subprocess
import os
import filecmp
import shlex
from pathlib import Path

PhaVa_rootDir = Path(__file__).parent.absolute()

PhaVaWorkingDir = os.path.dirname(PhaVa_rootDir) + "/tests/test_outputs"
refGenomePath = os.path.dirname(PhaVa_rootDir) + "/tests/GCA_000011065.1_ASM1106v1_genomic.fna"
readsDirectory = os.path.dirname(PhaVa_rootDir) + "/tests/simulated_reads"
comparisonOutputFile = os.path.dirname(PhaVa_rootDir) + "/tests/test_outputs/comparison_inverton_hits.tsv"
PhaVa = os.path.dirname(PhaVa_rootDir) + "/bin/phava"
filterFile = os.path.dirname(PhaVa_rootDir) + "/tests/runFilter.sh"

def runPhaVaPipeline():
	# Running the PhaVa locate process
	locate_command = 'python ' + PhaVa + '  locate -d ' + PhaVaWorkingDir + ' -t 16 -i ' + refGenomePath
	formatted_locate_command = shlex.split(locate_command)
	locate_process = subprocess.Popen(formatted_locate_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	locate_out, locate_err = locate_process.communicate()
	print('PhaVa Locate:')
	print(locate_out)
	print(locate_err)
	
	# Running the PhaVa create process
	create_command = 'python ' + PhaVa + '  create -d ' + PhaVaWorkingDir + ' -t 16 -i ' + refGenomePath
	formatted_create_command = shlex.split(create_command)
	create_process = subprocess.Popen(formatted_create_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	create_out, create_err = create_process.communicate()
	print('PhaVa Create:')
	print(create_out)
	print(create_err)

	# Running the PhaVa ratio process over all the reads in a given directory
	for filename in os.listdir(readsDirectory):
		f = os.path.join(readsDirectory, filename)
		if os.path.isfile(f):
			ratio_command = 'python ' + PhaVa + '  ratio -d ' + PhaVaWorkingDir + ' -t 16 -r ' + f
			formatted_ratio_command = shlex.split(ratio_command)
			ratio_process = subprocess.Popen(formatted_ratio_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
			ratio_out, ratio_err = ratio_process.communicate()
			print('PhaVa Ratio for ' + f + ':')
			print(ratio_out)
			print(ratio_err)

	# Combining all .tsv output files into a single .tsv file
	# Can be done by running the filter bash script previously written
	# Keeps only the invertons w/ >= 4 reverse hits and rev/all ratio > 0.009
	filter_command = 'bash ' + filterFile
	formatted_filter_command = shlex.split(filter_command)
	filter_process = subprocess.call(formatted_filter_command)

def test_comparePhaVaOutputs():
	runPhaVaPipeline()
	# Want to compare the just-generated output file to the previously manually-constructed output file
	assert filecmp.cmp('test_outputs/filtered_inverton_hits.tsv',  comparisonOutputFile, shallow=False) 
