#!/usr/bin/env python
import os
from Bio import SeqIO
import argparse
from multiprocessing import Pool
import subprocess
import time
import sys

### Function definitons:
def log(txt, out):
    if LOGLEVEL > 0:
        print(txt)
    out.write(txt+'\n')
    out.flush()


def process_psl(sample, contigs, psl):
    """
    Iterate over a PSL file, pull out the higest scoring hit and do a number of other
    calculations. TODO: Break this into a couple of functions to clean up the logic a bit.
    """
    sequences = SeqIO.index("./finished_" + sample + "/" + contigs, "fasta")  # how should i avoid doing this?
    blatfinal = {}
    for line in open(psl, 'r'):
        line = line.strip().split()
        matches = int(line[0])
        strand = line[8]
        qname = line[9]
        qsize = int(line[10])
        qstart = int(line[11])+1
        qend = int(line[12])
        tname = line[13]
        if len(tname.split('_:_')) > 1:
            tnamesplit = tname.split('_:_')[1]
        elif len(tname.split('_:_')) == 1:
            tnamesplit = tname
        tsize = int(line[14])
        tstart = int(line[15])+1
        tend = int(line[16])
        blocksizes = [int(m) for m in line[18].split(",")[0:-1]]
        qstarts = [int(m)+1 for m in line[19].split(",")[0:-1]]
        tstarts = [int(m)+1 for m in line[20].split(",")[0:-1]]
        if tstart == 1 and qstart > 1:
            hsubsize = (qstart-1)
        elif tend == tsize and qend < qsize:
            hsubsize = (qsize-qend)
        else:
            hsubsize = 0
        #calculate percent id
        id = matches/(qsize-hsubsize)
        #calculate cigar
        cigar = []
        cigarstring = str()
        if strand == "+":
            sqstart = qstart
            sqend = qend
        elif strand == "-":
            sqstart = (qsize-qend)+1
            sqend = (qsize-qstart)+1
        if sqstart > 1:
            L = sqstart-1
            cigar.append(str(L)+"S")
        if len(blocksizes) == 1:
            cigar.append(str(blocksizes[0])+"M")
        elif len(blocksizes) > 1:
            for i in range(len(blocksizes[0:-1])):
                cigar.append(str(blocksizes[i])+"M")
                if qstarts[i] + blocksizes[i] != qstarts[i+1]:
                    L = qstarts[i+1] - (qstarts[i] + blocksizes[i])
                    cigar.append(str(L)+"I")
                if tstarts[i] + blocksizes[i] != tstarts[i+1]:
                    L = tstarts[i+1]-(tstarts[i]+blocksizes[i])
                    cigar.append(str(L)+"D")
            cigar.append(str(blocksizes[-1])+"M")
        if qsize > sqend:
            L = qsize-sqend
            cigar.append(str(L)+"S")
        for i in cigar:
            cigarstring = cigarstring+i
        #pull out sequence
        qseq = sequences[qname].seq
        #only keep if it hits the right target
        if qname.split('_:_')[1] == tnamesplit:
            #if hit doesn't match something in blatfinal already, append
            if qname not in blatfinal.keys():
                blatfinal.update({qname: [tname, tstart, tend, strand, cigarstring, qseq, id]})
            #if hit matches something in blatfinal, replace if id is higher
            elif qname in blatfinal.keys():
                if id > blatfinal[qname][6]:
                    blatfinal.update({qname: [tname, tstart, tend, strand, cigarstring, qseq, id]})
    return blatfinal


def run_blat(sample, targets, contigs, outf_s):
    """
    Given a sample name and targets file, this will run blat, sending all output to a log
    file, and returning the PSL file path.
    """
    print("processing " + sample)
    #this needs to be more flexible
    originalcontigs = list(SeqIO.parse("./finished_" + sample + "/" + contigs, "fasta"))
    if len(originalcontigs) == 0:
        print("./finished_" + sample + "/" + contigs + " file empty")
    elif len(originalcontigs) > 0:
        args = ["blat", "-noHead", targets, "./finished_" + sample + "/" + contigs,
                "./finished_" + sample + "/blat_" + sample + ".psl"]
        ret = subprocess.call(args, stdout=outf_s, stderr=outf_s)
        if ret != 0:
            print "Error running blat for Sample:", sample
            return(None)
        else:
            psl = "./finished_" + sample + "/blat_" + sample + ".psl"
            return(psl)


def check_status(results):
    """
    Checks the stats of a dictionary full of jobs and returns the number of jobs which
    have not finished.
    """
    unfinished = 0
    finished = 0
    for r in results:
        if not results[r].ready():
            unfinished += 1
        else:
            finished += 1
    return unfinished


def process_sample(sample, targets, contigs):
    """
    Wrapper function which handles processing for a single sample.
    """
    outf_s = open("./finished_" + sample + '/align_log.txt', 'w')
    #First run blat to generate PSL:
    psl = run_blat(sample, targets, contigs, outf_s)
    #TODO add appropriate behavior if output is None
    psl_processed = process_psl(sample, psl)

    with open(str("./finished_" + sample + "/processedcontigs.fa"), "w") as outfile:
        for contig_id, value in psl_processed.items():
            #write out relevant info in line
            rec = str(">"+contig_id+" "+value[0]+" "+str(value[1])+" "+str(value[2])+" "+value[3]+" "+str(value[4])+"\n"+value[5]+"\n")
            outfile.write(rec)


############### Main #######################
parser = argparse.ArgumentParser()
parser.add_argument('-t', action='store', dest='targets', default=None,
                    help='targets file for mapping (uses path specified in config file by default')
parser.add_argument('-p', action='store', dest='processes', default=7,
                    help='Number of processes to use', type=int)
parser.add_argument("-c", action='store', dest='contigs', default='contigs.fasta',
                    help='alternate names for contigs, i.e. if hets have been injected.')

args = parser.parse_args()
processes = int(args.processes)
targets = args.targets
contigs = args.contigs

#processes = 7
samples = []

config = open('./ARC_config.txt', 'r').read().strip().split('\n')

for line in config:
    if 'reference=' in line and targets is None:
        targets = line.strip().split('=')[1].strip()
    if not line.startswith('#'):
        if not line.startswith('Sample_ID'):
            if len(line.split()) > 1:
                samples.append(line.split()[0])

#Remove duplicate sample ids:
samples = list(set(samples))

#multiprocessing wrapper
p = Pool(processes=processes, maxtasksperchild=1)
results = {}

for sample in samples:
    results[sample] = p.apply_async(process_sample, (sample, targets, contigs, ))

allfinished = False
while not allfinished:
    time.sleep(5)
    np = check_status(results)
    if np == 0:
        allfinished = True
