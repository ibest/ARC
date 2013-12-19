#!/usr/bin/env python
"""
This script takes a finished ARC run, maps the reads which were recruited per-target against
the contigs which were assembled from those reads, and then calls varaints.

New strategy:
    iterate over samples and targets withing samples, pass sample + target + i to a
    multiprocessing call

1) Read in ARC_config.txt and get SAMPLE_ID

2) For each sample in SAMPLE_ID
    - for each set of contigs belonging to a target in finished_SAMPLE/contigs.fasta
        -extract contigs and PE + SE reads
        -map reads
        - SAM --> sort --> filter by mapq > 10 --> BAM
        - call variants to a TARGET_variants.vcf
    -combine all TARGET_variants.vcf into a SAMPLE_variants.vcf
    -inject hets into contigs
        -filter variants based on QUAL, DP (or DP4), perhaps other things
        -hets get encoded with ambiguity codes
        -homozygous other (1/1) calls get converted to "N" reflecting the lack of
            confidence in this position.
        -indels are ignored
        -Produce a new set of contigs with injected HETS

"""
from Bio import SeqIO
import os
import sys
import time
#import glob
from multiprocessing import Pool
#from Bio.SeqRecord import SeqRecord
from copy import deepcopy
import subprocess
import argparse


#Function definitions:
def log(txt, out):
    if VERBOSE:
        print(txt)
    out.write(txt+'\n')
    out.flush()


def sp_caller(args, out):
    try:
        log(' '.join(args), out)
        ret = subprocess.call(args, stdout=out, stderr=out)
        if ret != 0:
            log("Error running: ", ' '.join(args))
            raise Exception("Error running: ", ' '.join(args))
    except Exception as e:
        raise e


def make_vcf(i, j, target, PE, SE, caller='GATK'):
    """
    i is sample
    j is target

    """
    #output handle for stderr and stdout:
    with open('./make_vcf_temp/S%s_%s.log' % (i, j), 'a') as out:
        log("Starting processing for target %s" % target, out)
        #Build indexes:
        if os.path.exists("./make_vcf_temp/idx_%s_%s" % (i, j)):
            os.system("rm -rf ./make_vcf_temp/idx_%s_%s" % (i, j))
        os.system("mkdir ./make_vcf_temp/idx_%s_%s" % (i, j))
        log("Building index..", out)
        args = ['bowtie2-build', '-f', './make_vcf_temp/ref_%s_%s.fasta' % (i, j)]
        args += ['./make_vcf_temp/idx_%s_%s/idx' % (i, j)]
        sp_caller(args, out)
        #cmd = "bowtie2-build -f ./make_vcf_temp/ref_%s_%s.fasta" % (i, j)
        #cmd += " ./make_vcf_temp/idx_%s_%s/idx >>./make_vcf_temp/S%s_%s.log 2>&1" % (i, j, i, j)
        #log(cmd, out)
        #os.system(cmd)

        #Do mapping:
        log("Running bowtie..", out)
        ## Note that RG fields are included in order to make the output compatible with GATK
        args = ['bowtie2', '--seed', '42', '-p', '2', '-t', '-I', '0', '-X', '1500',
                '--rg-id', 'none', '--rg', 'PL:ILLUMINA', '--rg', 'SM:none', '-x',
                './make_vcf_temp/idx_%s_%s/idx' % (i, j), PE, SE, '-S',
                './make_vcf_temp/tmp_%s_%s.sam' % (i, j)]
        sp_caller(args, out)
        # cmd = "bowtie2 --seed 42 -p 2 -t -I 0 -X 1500 --rg-id none --rg PL:ILLUMINA --rg SM:none"
        # cmd += " -x ./make_vcf_temp/idx_%s_%s/idx" % (i, j)
        # cmd += PE + SE + " -S ./make_vcf_temp/tmp_%s_%s.sam >> ./make_vcf_temp/S%s_%s.log 2>&1" % (i, j, i, j)
        #log(' '.join((str(i), str(j), cmd)), out)
        #os.system(cmd)

        # Screen for low map q, sort and convert to bam
        log("Filtering and converting to bam", out)
        #Note, you can't use subprocess.call() with this command because stdout is getting
        # piped from one call to the next and subprocess.call() seems to mess this up:
        cmd = "samtools view -q 10 -bS ./make_vcf_temp/tmp_%s_%s.sam | samtools sort - ./make_vcf_temp/tmp_%s_%s" % (i, j, i, j)
        cmd += " >> ./make_vcf_temp/S%s_%s.log 2>&1" % (i, j)
        log(cmd, out)
        os.system(cmd)

        #Create index for GATK:
        args = ['samtools', 'index', './make_vcf_temp/tmp_%s_%s.bam' % (i, j)]
        sp_caller(args, out)
        args = ['samtools', 'faidx', "./make_vcf_temp/ref_%s_%s.fasta" % (i, j)]
        sp_caller(args, out)

        #Call variants:
        if caller == 'GATK':
            #Build dict for GATK:
            args = [
                'java', '-jar', PICARDPATH + 'CreateSequenceDictionary.jar',
                'REFERENCE=./make_vcf_temp/ref_%s_%s.fasta' % (i, j),
                'OUTPUT=./make_vcf_temp/ref_%s_%s.dict' % (i, j)
            ]
            sp_caller(args, out)

            #Call with GATK:
            args = [
                'java', '-Xmx6g', '-jar', GATKPATH + 'GenomeAnalysisTK.jar',
                '-T', 'HaplotypeCaller', '-R', './make_vcf_temp/ref_%s_%s.fasta' % (i, j),
                '-I', './make_vcf_temp/tmp_%s_%s.bam' % (i, j),
                '-o', './make_vcf_temp/tmp_%s_%s.vcf' % (i, j)]
            sp_caller(args, out)

        else:
            # Call variants:
            log("Calling variants", out)
            cmd = "samtools mpileup -D -u -f ./make_vcf_temp/ref_%s_%s.fasta ./make_vcf_temp/tmp_%s_%s.bam |" % (i, j, i, j)
            cmd += " bcftools view -vgc - > ./make_vcf_temp/tmp_%s_%s.vcf 2> ./make_vcf_temp/S%s_%s.log" % (i, j, i, j)
            log(cmd, out)
            os.system(cmd, out)

        # Cleanup:
        log("Cleaning up...", out)
        os.system("rm ./make_vcf_temp/PE1_%s_%s.fastq" % (i, j))
        os.system("rm ./make_vcf_temp/PE2_%s_%s.fastq" % (i, j))
        os.system("rm ./make_vcf_temp/SE_%s_%s.fastq" % (i, j))
        os.system("rm ./make_vcf_temp/tmp_%s_%s.sam" % (i, j))
        os.system("rm ./make_vcf_temp/tmp_%s_%s.bam" % (i, j))
        os.system("rm ./make_vcf_temp/tmp_%s_%s.bam.bai" % (i, j))
        os.system("rm -rf ./make_vcf_temp/idx_%s_%s" % (i, j))
        os.system("rm ./make_vcf_temp/ref_%s_%s.*" % (i, j))
        log("target %s complete" % target, out)


def check_status(results):
    """
    run through each of the lists in results, if all of the jobs are done return the key
    """
    for sample in results:
        unfinished = 0
        finished = 0
        for r in results[sample]:
            if not r[2].ready():
                unfinished += 1
            else:
                finished += 1
        print "Sample ", sample, "finished:", finished, " unfinished:", unfinished
        if unfinished == 0:
            return sample
    return None


def inject_variants(sample):
    """
    A function which takes a sample ID, assumes that the appropriate vcf and contig files
    will exist in the correct locations (based on ARC folder structure), and creates a new
    set of contig files containing ambiguity codes.
    Filtering for SNPS:
    1) Indels are ignored
    2) Hets are only injected if support is sufficient
        -supporting reads between 30 and 70%
    """
    vcf_file = "./make_vcf_temp/%s.vcf" % sample
    contig_file = './finished_%s/contigs.fasta' % sample
    out_file = './finished_%s/contigs_hets.fasta' % sample

    vcff = open(vcf_file, 'r')
    contigf = open(contig_file, 'r')
    outf = open(out_file, 'w')
    out = open('./make_vcf_temp/%s.log' % sample, 'w')

    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  S151_s.bam
    #Abiguity codes:
    ACs = {
        'AC': 'M',
        'CA': 'M',
        'AG': 'R',
        'GA': 'R',
        'AT': 'W',
        'TA': 'W',
        'CG': 'S',
        'GC': 'S',
        'CT': 'Y',
        'TC': 'Y',
        'GT': 'K',
        'TG': 'K'
    }

    #Parse vcf and load it into a dictionary:
    vcf = {}  # contig IDs are keys,
    # seqlens = {}  # Store sequence lengths

    filters = {}
    filters['minQual'] = 100
    filters['minDP4'] = .1
    filters['mincov'] = 11

    for l in vcff:
        if l[0] != '#':
            l = l.split('\t')
            if len(l) != 10:
                log("ERROR in sample", sample, "VCF record has", str(len(l)), "columns instead of 10.", out)
                log(' '.join(l), out)
                continue
            # variants which don't start with DP are INDELS
            contig = l[0]
            pos = l[1]
            ref = l[3]
            alt = l[4]
            qual = float(l[5])
            gt = l[9][0:3]

            #Currently we only handle HET calls (ref and alt have length of 1)
            if len(ref) == 1 and len(alt) == 1 and qual > filters['minQual']:
                #pull details out of info line and check minDP4
                if contig not in vcf:
                    vcf[contig] = {}
                vcf[contig][pos] = {'ref': ref,
                                    'alt': alt,
                                    'gt': gt
                                    }
                #print "SNP Inserted: Contig %s, pos %s, ref %s, alt %s, gt %s, r %s, a %s, qual %s" % (contig, pos, ref, alt, gt, r, a, float(l[5]))
                log("Variant passed filter: Contig %s, pos %s, ref %s, alt %s, gt %s, qual %s" % (contig, pos, ref, alt, gt, qual), out)

                # info = l[7].split(';')
                # r = a = 0
                # for i in info:
                #     i = i.split('=')
                #     if i[0] == 'DP4':
                #         r = float(int(i[1].split(',')[0]) + int(i[1].split(',')[1]))
                #         a = float(int(i[1].split(',')[2]) + int(i[1].split(',')[3]))
                # #Coverage filter:
                # if (r+a) < filters['mincov']:
                #     print "Filtered for mincoverage:", r+a
                #     print "Contig %s, pos %s, ref %s, alt %s, gt %s, r %s, a %s, qual %s" % (contig, pos, ref, alt, gt, r, a, float(l[5]))
                #     continue

                #het filter:
                # if gt == '0/1' and float(a)/(r + a) < filters['minDP4'] or float(r)/(r + a) < filters['minDP4']:
                #     print "Filtered for het and low coverage", float(a)/(r+a), float(r)/(r+a)
                #     print "Contig %s, pos %s, ref %s, alt %s, gt %s, r %s, a %s, qual %s" % (contig, pos, ref, alt, gt, r, a, float(l[5]))
                #     continue
                # if a/(r + a) < filters['minDP4'] or r/(r + a) < filters['minDP4'] or (r+a) < filters['mincov']:
                #     print "Failed filter:"
                #     print "Contig %s, pos %s, ref %s, alt %s, gt %s, r %s, a %s, qual %s" % (contig, pos, ref, alt, gt, r, a, float(l[5]))
                #     continue

                #finally insert this variant into the dictionary
    #Parse the contigs and apply variants:
    #outf = open(outf, 'w')
    # Note: VCF uses a 1-based index, so subtract 1 from index
    totalhets = 0
    totalhomozygous = 0
    for contig in SeqIO.parse(contigf, 'fasta'):
        contigm = contig.seq.tomutable()
        if contig.id in vcf:
            print contig.id
            variants = vcf[contig.id]
            print variants
            for pos in variants.keys():
                if variants[pos]['gt'] == '0/1':
                    # Het call
                    het = ''.join([variants[pos]['ref'], variants[pos]['alt']])
                    ac = ACs[het]
                    log(str(variants[pos]) + contigm[int(pos) - 1] + ac, out)
                    contigm[int(pos) - 1] = ac
                    totalhets += 1
                if variants[pos]['gt'] == '1/1':
                    #contigm[int(pos) - 1] = variants[pos]['alt']
                    contigm[int(pos) - 1] = 'N'
                    totalhomozygous += 1
        else:
            log("------" + sample + contig.name + ": No Variants pass filter--------", out)
        #All variants processed:
        new_seq = deepcopy(contig)
        new_seq.seq = contigm.toseq()
        #SeqRecord(seq=contigm.toseq(), id=contig.id, description=contig.description)
        SeqIO.write(new_seq, outf, 'fasta')
    outf.close()
    log(sample + "Total hets: %s " % totalhets + " Totalhomozygous:%s" % totalhomozygous, out)


## MAIN:

#Globals
#GATKPATH = '/bio/software/HT_sequence/GATK/'
#PICARDPATH = '/bio/software/HT_sequence/picard/picard-tools-1.103/'
#LOGLEVEL = 1

parser = argparse.ArgumentParser()
parser.add_argument('-g', action='store', dest='GATKPATH', default=None,
                    help='Path where GATK jar file is stored.')
parser.add_argument('-p', action='store', dest='PICARDPATH', default=None,
                    help='Path where Picard jar files are stored.')
parser.add_argument('-n', action='store', dest='processes', default=7,
                    help='Number of processes to use', type=int)
parser.add_argument('-c', action='store', dest='configf', default='ARC_config.txt',
                    help='ARC config file if different from ARC_config.txt.')
parser.add_argument('-v', action='store_true', help='Enable verbose output.')

#Deal with params:
args = parser.parse_args()
processes = args.processes
GATKPATH = args.GATKPATH
PICARDPATH = args.PICARDPATH
VERBOSE = args.v
CONFIGF = args.configf
#targets = args.targets

if GATKPATH is None or os.path.exists(GATKPATH + '/GenomeAnalysisTK.jar') is False:
    print 'Error, GATK jar file:' + str(GATKPATH) + ' /GenomeAnalysisTK.jar does not exist.'
    print 'Please indicate the path to GenomeAnalysisTK.jar with -g.'
    os.sys.exit()

if PICARDPATH is None or os.path.exists(PICARDPATH + '/CreateSequenceDictionary.jar') is False:
    print 'Error, Picard jar file:' + str(PICARDPATH) + ' /CreateSequenceDictionary.jar does not exist.'
    print 'Please indicate the path to CreateSequenceDictionary.jar with -p.'
    os.sys.exit()


###########
# Setup a working folder:
if os.path.exists("./make_vcf_temp"):
    print "make_vcf_temp already exists, please delete it before running this script"
    sys.exit()
os.mkdir("make_vcf_temp")

## Parse config to get a list of SAMPLE_IDs
samples = {}
format = None
if os.path.exists(CONFIGF):
    for line in open(CONFIGF):
        if line[0] != '#':
            line = line.strip().split()
            if line[0] != 'Sample_ID':
                if line[0] not in samples:
                    samples[line[0]] = {
                        'PE1': False,
                        'PE2': False,
                        'SE': False
                        }
                samples[line[0]][line[2]] = True
        elif line[0:2] == "# ":
            param = line.strip().split()[1].split("=")
            if param[0] == 'format':
                format = param[1]
else:
    print "Can't find config file:", CONFIGF
    sys.exit()

p = Pool(processes=processes, maxtasksperchild=1)
results = {}

lsamples = samples.keys()

for i, sample in enumerate(lsamples):
    #sample = lsamples[i]
    results[sample] = []
    print "Sample: ", sample

    contigsidx = SeqIO.index("./finished_%s/contigs.fasta" % sample, "fasta")

    #Get a list of targets, use set to get unique targets
    targets = list(set(map(lambda x: x.strip().split("_:_")[1], contigsidx.keys())))

    #Get contigs, reads, and do mapping and variant calling for each target:
    for j, target in enumerate(targets):
        #target = targets[j]
        #Get contigs, reads, and do mapping and variant calling for each target:
        #print "Target: ", target
        outf = open("./make_vcf_temp/ref_%s_%s.fasta" % (i, j), 'w')
        for contig in filter(lambda x: target == x.strip().split("_:_")[1], contigsidx.keys()):
            #print contig
            SeqIO.write(contigsidx[contig], outf, "fasta")
        outf.close()
        #Get reads for this target:
        if format == "fastq":
            A = 3
        elif format == "fasta":
            A = 1
        else:
            print "format %s not recognized, exiting" % format
        PE = SE = ""
        if samples[sample]['PE1'] and samples[sample]['PE2']:
            #For reasons that remain unclear to me, grep sometimes outputs a '--' line after the qual, causing major problems
            # Because bowtie2 views it as 2 additional quality characters and fails.
            cmd = "grep -A %s '%s$' ./finished_%s/PE1.%s > ./make_vcf_temp/PE1_%s_%s.%s" % (A, target, sample, format, i, j, format)
            #print cmd
            os.system(cmd)
            cmd = "grep -A %s '%s$' ./finished_%s/PE2.%s > ./make_vcf_temp/PE2_%s_%s.%s" % (A, target, sample, format, i, j, format)
            #print cmd
            os.system(cmd)
            PE = " -1 ./make_vcf_temp/PE1_%s_%s.%s -2 ./make_vcf_temp/PE2_%s_%s.%s" % (i, j, format, i, j, format)
        if samples[sample]['SE']:
            cmd = "grep -A %s '%s$' ./finished_%s/SE.%s > ./make_vcf_temp/SE_%s_%s.%s" % (A, target, sample, format, i, j, format)
            #print cmd
            os.system(cmd)
            SE = " -U ./make_vcf_temp/SE_%s_%s.%s" % (i, j, format)

        results[sample].append((i, j, p.apply_async(make_vcf, (i, j, target, PE, SE,))))
        #make_vcf(i,j,target,PE, SE)

#Check results, if anything is done combine the vcf
while len(results) > 0:
    finished_sample = check_status(results)
    if finished_sample is not None:
        #Sample is finished, combine and delete VCF files:
        contigs = []
        header = []
        variants = []
        for result in results[finished_sample]:
            i, j = result[0:2]
            #print "i:", i, " j:", j
            if os.path.exists("./make_vcf_temp/tmp_%s_%s.vcf" % (i, j)):
                for line in open("./make_vcf_temp/tmp_%s_%s.vcf" % (i, j), 'r'):
                    if j == 0 and line[0] == "#" and line[0:9] != "##contig=":
                        header.append(line)
                    elif line[0:9] == "##contig=":
                        contigs.append(line)
                    elif line[0] != '#':
                        variants.append(line)
            os.system("rm ./make_vcf_temp/tmp_%s_%s.vcf" % (i, j))
            os.system("rm ./make_vcf_temp/tmp_%s_%s.vcf" % (i, j))
        outf = open('./make_vcf_temp/%s.vcf' % finished_sample, 'w')
        outf.write("".join(header[0:4]))
        outf.write("".join(contigs))
        outf.write("".join(header[4:]))
        outf.write("".join(variants))
        outf.close()
        del results[finished_sample]
        inject_variants(finished_sample)

    else:
        time.sleep(1)


