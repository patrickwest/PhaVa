#!/usr/bin/env python

import sys
import os
import subprocess
import errno
import shlex
import filecmp

# where is this script
path_phava = os.path.realpath(__file__)
path_array = path_phava.split("/")
relative_path = "/".join(path_array[0:-2])
relative_path = relative_path + "/"


# function to check if a tool is installed
def is_tool(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True


def is_tool_and_return0(name):
    try:
        devnull = open(os.devnull)
        popenCMD = shlex.split(name)
        child = subprocess.Popen(popenCMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        streamdata = child.communicate()
        rc = child.wait()
        if rc == 0:
            return True, str(streamdata[0], 'utf-8'), str(streamdata[1], 'utf-8')
        else:
            return False, ' ', ' '
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False, ' ', ' ' 
    return True, ' ', ' '


def main(args=None):
    sys.stderr.write(" ------------------------------------------------\n")
    sys.stderr.write("|                TEST PhaVa                      |\n")
    sys.stderr.write(" ------------------------------------------------\n\n")

    error_found = False

    # check for tools being installed
    sys.stderr.write("1-- Tools and versions:\n\n")

    # check python version
    sys.stderr.write("- python:   ")
    python_version = sys.version_info
    if python_version >= (3, 0, 0):
        sys.stderr.write("      success, found v" + str(python_version[0]) + "." +
                         str(python_version[1]) + "." + str(python_version[2]) + "\n")
    else:
        sys.stderr.write("      ERROR: found v " + str(python_version[0]) + "." +
                         str(python_version[1]) + "." + str(python_version[2]) + ". Required 3.0.0 (or higher)\n")
        error_found = True

    # check minimap
    sys.stderr.write("- minimap2:      ")
    minimap, minimap_out, minimap_err = is_tool_and_return0("minimap2 --version")
    if minimap:
        sys.stderr.write(" success, found v" + minimap_out)
    else:
        sys.stderr.write(" ERROR. Minimap2 is not installed\n")
        error_found = True

    # check samtools
    sys.stderr.write("- samtools:      ")
    samtools, samtools_out, samtools_err = is_tool_and_return0("samtools --version")
    if samtools:
        samtools_version = samtools_out.split('\n')[0].replace("samtools ", "")
        sys.stderr.write(" success, found v" + samtools_version + "\n")
    else:
        sys.stderr.write(" ERROR. Samtools is not installed\n")
        error_found = True

    # check einverted
    sys.stderr.write("- einverted:     ")
    einverted, einverted_out, einverted_err = is_tool_and_return0("einverted --version")
    if einverted:
       sys.stderr.write(" success, found v" + einverted_err.replace("EMBOSS:", ""))
    else:
        sys.stderr.write(" ERROR. EMBOSS is not installed\n")
        error_found = True

    # optionally, check bowtie
    sys.stderr.write("- bowtie2 (opt.):")
    bowtie, bowtie_out, bowtie_err = is_tool_and_return0("bowtie2 --version")
    if bowtie:
        bowtie_version = bowtie_out.split('\n')[0].split("version ")[1]
        sys.stderr.write(" success, found v" + bowtie_version + "\n")
    else:
        sys.stderr.write(" WARNING. Bowtie2 not found. The pipeline will not work for short reads!\n")

    # run pipeline
    if error_found:
        sys.stderr.write("### Skipping test pipeline because of errors\n")
        return 1

    sys.stderr.write("\n2-- Run test pipeline:\n\n")

    # set correct directories
    PhaVaWorkingDir = os.path.dirname(relative_path) + "/tests/test_outputs"
    refGenomePath = os.path.dirname(relative_path) + "/tests/test_small.fna"
    PhaVa = os.path.dirname(relative_path) + "/bin/phava"
    simulatedReads = os.path.dirname(relative_path) + "/tests/simulated_reads.fastq"
    expOutput = os.path.dirname(relative_path) + "/tests/expected_output.tsv"
    pipelineOutput = os.path.dirname(relative_path) + "/tests/test_outputs/simulated_reads.fastq_vs_test_small.fna_ratio.tsv"

    # locate
    locate_command = 'python ' + PhaVa + '  locate -d ' + PhaVaWorkingDir + ' -l -i ' + refGenomePath
    formatted_locate_command = shlex.split(locate_command)
    print('PhaVa Locate:')
    locate_process = subprocess.Popen(formatted_locate_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                      universal_newlines=True)
    locate_out, locate_err = locate_process.communicate()
    print(locate_out)
    print(locate_err)

    # create
    create_command = 'python ' + PhaVa + '  create -d ' + PhaVaWorkingDir + ' -l -i ' + refGenomePath
    formatted_create_command = shlex.split(create_command)
    print('PhaVa Create:')
    create_process = subprocess.Popen(formatted_create_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                      universal_newlines=True)
    create_out, create_err = create_process.communicate()
    print(create_out)
    print(create_err)

    # ratio computation
    ratio_command = 'python ' + PhaVa + '  ratio -d ' + PhaVaWorkingDir + ' -l -r ' + simulatedReads
    formatted_ratio_command = shlex.split(ratio_command)
    print('PhaVa Ratio:')
    ratio_process = subprocess.Popen(formatted_ratio_command, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, universal_newlines=True)
    ratio_out, ratio_err = ratio_process.communicate()
    print(ratio_out)
    print(ratio_err)

    # compare to expected output
    assert filecmp.cmp(expOutput, pipelineOutput, shallow=False)

