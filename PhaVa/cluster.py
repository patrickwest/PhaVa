#!/usr/bin/env python
import os
from os.path import exists
import pickle
import csv
import logging
import subprocess
import PhaVa.utils
from PhaVa.fileManager import WorkDirectory
from collections import defaultdict, Counter


def main(args):
    irDb = PhaVa.utils.IRDb()
    # check if one or many IRs need to be combined
    if ',' in args.dir:
        all_dirs = args.dir.split(',')
        new_outdir = args.new_dir
    else:
        all_dirs = [args.dir]
        if args.new_dir == './phava_out':
            new_outdir = args.dir
        else:
            new_outdir = args.new_dir
    wd = WorkDirectory(new_outdir)
    # load all the database pickles
    # make this a dictionary (and keep the name)
    if len(all_dirs) == 1:
        all_irs = {os.path.basename(all_dirs[0].rstrip('/')): unpickleDb(all_dirs[0])}
    else:
        all_irs = {os.path.basename(x.rstrip('/')): unpickleDb(x) for x in all_dirs}
    # create fasta with flanks and invertable regions
    flanks = {}
    invertibles = {}
    for db in all_irs.keys():
        IRs = all_irs[db].IRs
        for ir in IRs:
            IRs[ir].flankSize = args.flankSize
            leftFlankStart = IRs[ir].leftStart - IRs[ir].flankSize
            if leftFlankStart < 0:
                leftFlankStart = 0
            rightFlankEnd = IRs[ir].rightStop + IRs[ir].flankSize
            leftFlank = all_irs[db].genome[IRs[ir].chr][leftFlankStart:IRs[ir].leftStart]
            rightFlank = all_irs[db].genome[IRs[ir].chr][IRs[ir].rightStop:rightFlankEnd]

            flanks[db + ':' + ir] = leftFlank + rightFlank
            invertibles[db + ':' + ir] = IRs[ir].middleSeq
    with open(new_outdir + '/intermediate/flanks.fa', 'w+') as output_handle:
        for name in flanks.keys():
            output_handle.write('>' + name + '\n')
            output_handle.write(flanks[name] + '\n')
    with open(new_outdir + '/intermediate/invertibles.fa', 'w+') as output_handle:
        for name in invertibles.keys():
            output_handle.write('>' + name + '\n')
            output_handle.write(invertibles[name] + '\n')

    # call mmseqs to cluster everything
    if not os.path.exists(new_outdir + '/intermediate/tmp_mmseqs'):
        os.makedirs(new_outdir + '/intermediate/tmp_mmseqs')

    mmseq_loc = new_outdir + "/intermediate/tmp_mmseqs"
    run_mmseqs(new_outdir + "/intermediate/flanks.fa", mmseq_loc, 'flanks', args.pident, args.cpus)
    run_mmseqs(new_outdir + "/intermediate/invertibles.fa", mmseq_loc, 'invertibles', args.pident, args.cpus)

    # parse the clustering tsv file and assign new ids
    logging.info("------Creating new clustered database------")
    cluster_flanks = {}
    cluster_invs = {}
    with open(new_outdir + '/intermediate/flanks.tsv', 'r') as flanks_file:
        flanks_reader = csv.reader(flanks_file, delimiter='\t')
        for row in flanks_reader:
            cluster_flanks[row[1]] = row[0]
    with open(new_outdir + '/intermediate/invertibles.tsv', 'r') as inv_file:
        inv_reader = csv.reader(inv_file, delimiter='\t')
        for row in inv_reader:
            cluster_invs[row[1]] = row[0]
    df_cluster = defaultdict(list)
    for member in cluster_invs.keys():
        df_cluster[(cluster_flanks[member], cluster_invs[member])].append(member)

    cluster_counter = 0
    df_result = []
    IR_list = {}
    genome_list = {}
    for group in df_cluster.keys():
        new_id = f'id{cluster_counter}'
        for m in df_cluster[group]:
            df_result.append([new_id, group[0], group[1], m])
        sdir, inv_id = df_cluster[group][0].split(':', 1)
        ir = all_irs[sdir].IRs[inv_id]
        new_genome_name = sdir + ":" + ir.chr
        genome_list[new_genome_name] = all_irs[sdir].genome[ir.chr]
        ir.chr = new_genome_name
        IR_list[new_id] = ir
        cluster_counter += 1
    irDb.IRs = IR_list
    irDb.genome = genome_list

    # save the new clustering into new list and to file
    with open(new_outdir + '/clustering.tsv', 'w') as clustering_file:
        clustering_writer = csv.writer(clustering_file, delimiter='\t')
        clustering_writer.writerow(['new_cluster_id', 'cluster_flank', 'cluster_inv', 'inverton_id'])
        clustering_writer.writerows(df_result)

    # create new irDb
    irDb.genomeName = 'clustered'
    exportIRs(irDb.IRs, new_outdir)
    # give some stats about the clustering
    get_stats(df_result)
    return irDb, new_outdir


def exportIRs(IRs, outpath):
    out = open(outpath + '/IRs.tsv', 'w')
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


def unpickleDb(dir):
    if exists(dir + "/irDb.pickle"):
        pickled_db = open(dir + "/irDb.pickle", 'rb')
        irDb = pickle.load(pickled_db)
        pickled_db.close()
        logging.info("------Finished unpickling IR database------")
        return irDb
    else:
        logging.info("------No pickled IR database------")
        return None


def run_mmseqs(fasta_file, tmp_loc, id, pident, cpus):
    logging.info("------Starting to cluster the " +  id +  "------")
    # create database
    command = "mmseqs createdb " + fasta_file + " " + tmp_loc + "/" + id + " --dbtype 2 -v 3"
    logging.info(command)
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    # linclust
    command = "mmseqs linclust " + tmp_loc + "/" + id + " " + tmp_loc + "/" + id + "_clustered " + \
              tmp_loc + "/tmp" + " --min-seq-id " + str(pident) + \
              " --cov-mode 0 -c 0.8 -v 3 --cluster-mode 0 --threads " + str(cpus)
    logging.info(command)
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    # get tsv
    command = "mmseqs createtsv " + tmp_loc + "/" + id + " " + tmp_loc + "/" + id + " " + \
              tmp_loc + "/" + id + "_clustered " + fasta_file.rstrip('.fa') + ".tsv --threads " + str(cpus)
    logging.info(command)
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    logging.info("------Finished clustering the " + id + "------")


def get_stats(clustering_df):
    n_all = len(clustering_df)
    n_clusters = len(set([x[0] for x in clustering_df]))
    n_flanks = len(set([x[1] for x in clustering_df]))
    n_invertibles = len(set([x[2] for x in clustering_df]))
    # biggest cluster
    counts = Counter([x[0] for x in clustering_df])
    counts = sorted(counts.values())


    logging.info("------Clustering info:")
    logging.info("      Total number of invertons: " + str(n_all))
    logging.info("      Total number of resulting clusters: " + str(n_clusters))
    logging.info("      Number of unique flanking regions: " + str(n_flanks))
    logging.info("      Number of unique invertible regions: " + str(n_invertibles))
    logging.info("      The five biggest clusters contain: " + str(counts[-5:]))



