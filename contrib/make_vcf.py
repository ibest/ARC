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
import glob
from multiprocessing import Pool
from Bio.SeqRecord import SeqRecord
from copy import deepcopy
import subprocess

#Globals
GATKPATH = '/bio/software/HT_sequence/GATK/GenomeAnalysisTK-2.4-9-g532efad/'
PICARDPATH = '/bio/software/HT_sequence/picard/picard-tools-1.103/'


def make_vcf(i, j, target, PE, SE, caller='GATK'):
    """
    i is sample
    j is target

    """
    #output handle for stderr and stdout:
    out = open('./make_vcf_temp/S%s_%s.log' % (i, j), 'a')
    #Build indexes:
    if os.path.exists("./make_vcf_temp/idx_%s_%s" % (i, j)):
        os.system("rm -rf ./make_vcf_temp/idx_%s_%s" % (i, j))
    os.system("mkdir ./make_vcf_temp/idx_%s_%s" % (i, j))
    #print "Building index.."
    cmd = "bowtie2-build -f ./make_vcf_temp/ref_%s_%s.fasta" % (i, j)
    cmd += " ./make_vcf_temp/idx_%s_%s/idx >>./make_vcf_temp/S%s_%s.log 2>&1" % (i, j, i, j)
    #print cmd
    os.system(cmd)

    #Do mapping:
    #print "Running bowtie.."
    ## Note that RG fields are included in order to make the output compatible with GATK
    cmd = "bowtie2 --seed 42 -p 2 -t -I 0 -X 1500 --rg-id none --rg PL:ILLUMINA --rg SM:none"
    cmd += " -x ./make_vcf_temp/idx_%s_%s/idx" % (i, j)
    cmd += PE + SE + " -S ./make_vcf_temp/tmp_%s_%s.sam >> ./make_vcf_temp/S%s_%s.log 2>&1" % (i, j, i, j)
    #print i, j, cmd
    os.system(cmd)

    # Screen for low map q, sort and convert to bam
    #print "Filtering and converting to bam"
    #Note, you can't use subprocess.call() with this command because stdout is getting
    # piped from one call to the next and subprocess.call() seems to mess this up:
    cmd = "samtools view -q 10 -bS ./make_vcf_temp/tmp_%s_%s.sam | samtools sort - ./make_vcf_temp/tmp_%s_%s" % (i, j, i, j)
    cmd += " >> ./make_vcf_temp/S%s_%s.log 2>&1" % (i, j)
    #print cmd
    os.system(cmd)

    #Create index for GATK:
    args = ['samtools', 'index', './make_vcf_temp/tmp_%s_%s.bam' % (i, j)]
    ret = subprocess.call(args, stdout=out, stderr=out)
    args = ['samtools', 'faidx', "./make_vcf_temp/ref_%s_%s.fasta" % (i, j)]
    ret = subprocess.call(args, stdout=out, stderr=out)


    #Call variants:
    if caller == 'GATK':
        #Build dict for GATK:
        args = [
            'java', '-jar', PICARDPATH + 'CreateSequenceDictionary.jar',
            'REFERENCE=', './make_vcf_temp/ref_%s_%s.fasta' % (i, j),
            'OUTPUT=', './make_vcf_temp/ref_%s_%s.dict' % (i, j)
        ]
        ret = subprocess.call(args, stdout=out, stderr=out)

        #Call with GATK:
        args = [
            'java', '-Xmx6g', '-jar', GATKPATH + 'GenomeAnalysisTK.jar', '-T', 'HaplotypeCaller', '-R', './make_vcf_temp/ref_%s_%s.fasta' % (i, j),
            '-I', './make_vcf_temp/tmp_%s_%s.bam' % (i, j),
            '-o', './make_vcf_temp/tmp_%s_%s.vcf' % (i, j)]
        ret = subprocess.call(args, stdout=out, stderr=out)
    else:
        #Call with samtools

# gtk <- "java -Xmx6g -jar /mnt/home/msettles/opt/src/GenomeAnalysisTK-2.7-4/GenomeAnalysisTK.jar"
# if (opt$generate_vcf){## Extract Unmapped Reads
#   vcf_out <- mclapply(bowtie, function(index){
#     try({
#       system(paste(gtk,
#                  "-T HaplotypeCaller",
#                  "-R",paste(index$target_path,"fasta",sep="."),
#                  "-I",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),
#                  "-o",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"vcf",sep=".")),
#                  ">",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"gtk","out",sep=".")),sep=" "));
#     })
#   },mc.cores=procs)
# }

    # Call variants:
    print "Calling variants"
    cmd = "samtools mpileup -D -u -f ./make_vcf_temp/ref_%s_%s.fasta ./make_vcf_temp/tmp_%s_%s.bam |" % (i, j, i, j)
    cmd += " bcftools view -vgc - > ./make_vcf_temp/tmp_%s_%s.vcf 2> ./make_vcf_temp/S%s_%s.log" % (i, j, i, j)
    #print cmd  1>> ./make_vcf_temp/S%s_%s.log
    os.system(cmd)

    # Cleanup:
    os.system("rm ./make_vcf_temp/PE1_%s_%s.fastq" % (i, j))
    os.system("rm ./make_vcf_temp/PE2_%s_%s.fastq" % (i, j))
    os.system("rm ./make_vcf_temp/SE_%s_%s.fastq" % (i, j))
    os.system("rm ./make_vcf_temp/tmp_%s_%s.sam" % (i, j))
    os.system("rm ./make_vcf_temp/tmp_%s_%s.bam" % (i, j))
    os.system("rm -rf ./make_vcf_temp/idx_%s_%s" % (i, j))
    os.system("rm ./make_vcf_temp/ref_%s_%s.*" % (i, j))


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
    filters['minQual'] = 15
    filters['minDP4'] = .1
    filters['mincov'] = 11

    for l in vcff:
        if l[0] != '#':
            l = l.split('\t')
            # variants which don't start with DP are INDELS
            if l[7][0:2] == 'DP' and float(l[5]) > filters['minQual']:
                contig = l[0]
                pos = l[1]
                ref = l[3]
                alt = l[4]

                gt = l[9][0:3]
                #pull details out of info line and check minDP4
                info = l[7].split(';')
                r = a = 0
                for i in info:
                    i = i.split('=')
                    if i[0] == 'DP4':
                        r = float(int(i[1].split(',')[0]) + int(i[1].split(',')[1]))
                        a = float(int(i[1].split(',')[2]) + int(i[1].split(',')[3]))
                #Coverage filter:
                if (r+a) < filters['mincov']:
                    print "Filtered for mincoverage:", r+a
                    print "Contig %s, pos %s, ref %s, alt %s, gt %s, r %s, a %s, qual %s" % (contig, pos, ref, alt, gt, r, a, float(l[5]))
                    continue

                #het filter:
                if gt == '0/1' and float(a)/(r + a) < filters['minDP4'] or float(r)/(r + a) < filters['minDP4']:
                    print "Filtered for het and low coverage", float(a)/(r+a), float(r)/(r+a)
                    print "Contig %s, pos %s, ref %s, alt %s, gt %s, r %s, a %s, qual %s" % (contig, pos, ref, alt, gt, r, a, float(l[5]))
                    continue
                # if a/(r + a) < filters['minDP4'] or r/(r + a) < filters['minDP4'] or (r+a) < filters['mincov']:
                #     print "Failed filter:"
                #     print "Contig %s, pos %s, ref %s, alt %s, gt %s, r %s, a %s, qual %s" % (contig, pos, ref, alt, gt, r, a, float(l[5]))
                #     continue

                #finally insert this variant into the dictionary
                if contig not in vcf:
                    vcf[contig] = {}
                vcf[contig][pos] = {'ref': ref,
                                    'alt': alt,
                                    'gt': gt
                                    }
                print "SNP Inserted: Contig %s, pos %s, ref %s, alt %s, gt %s, r %s, a %s, qual %s" % (contig, pos, ref, alt, gt, r, a, float(l[5]))

    #Parse the contigs and apply variants:
    #outf = open(outf, 'w')
    # Note: VCF uses a 1-based index, so subtract 1 from index
    totalhets = 0
    totalhomoozygous = 0
    for contig in SeqIO.parse(contigf, 'fasta'):
        contigm = contig.seq.tomutable()
        if contig.id in vcf:
            print contig.id
            variants = vcf[contig.id]
            print variants
            for pos in variants.keys():
                if variants[pos]['gt'] == '0/1':
                    # Het call
                    ac = ACs[''.join([variants[pos]['ref'], variants[pos]['alt']])]
                    print variants[pos], contigm[int(pos) - 1], ac
                    contigm[int(pos) - 1] = ac
                    totalhets += 1
                if variants[pos]['gt'] == '1/1':
                    #contigm[int(pos) - 1] = variants[pos]['alt']
                    contigm[int(pos) - 1] = 'N'
                    totalhomoozygous += 1
        else:
            print "------", sample, contig.name, ": No Variants pass filter--------"
        #All variants processed:
        new_seq = deepcopy(contig)
        new_seq.seq = contigm.toseq()
        #SeqRecord(seq=contigm.toseq(), id=contig.id, description=contig.description)
        SeqIO.write(new_seq, outf, 'fasta')
    outf.close()
    print sample, "Total hets: ", totalhets, " Totalhomoozygous:", totalhomoozygous


## MAIN:

############
## TODO: add some switches
configf = 'ARC_config.txt'


###########
# Setup a working folder:
if os.path.exists("./make_vcf_temp"):
    print "make_vcf_temp already exists, please delete it before running this script"
    sys.exit()
os.mkdir("make_vcf_temp")

## 1) Parse config to get a list of SAMPLE_IDs
samples = {}
format = None
if os.path.exists(configf):
    for line in open(configf):
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
    print "Can't find config file:", configf
    sys.exit()

p = Pool(processes=8, maxtasksperchild=1)
results = {}

lsamples = samples.keys()

for i in range(len(lsamples)):
    sample = lsamples[i]
    results[sample] = []
    print "Sample: ", sample

    contigsidx = SeqIO.index("./finished_%s/contigs.fasta" % sample, "fasta")

    #Get a list of targets, use set to get unique targets
    targets = list(set(map(lambda x: x.strip().split("_:_")[1], contigsidx.keys())))

    #Get contigs, reads, and do mapping and variant calling for each target:
    for j in range(len(targets)):
        target = targets[j]
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
            A = 2
        else:
            print "format %s not recognized, exiting" % format
        PE = SE = ""
        if samples[sample]['PE1'] and samples[sample]['PE2']:
            cmd = "grep -A %s '%s' ./finished_%s/PE1.%s > ./make_vcf_temp/PE1_%s_%s.%s" % (A, target, sample, format, i, j, format)
            #print cmd
            os.system(cmd)
            cmd = "grep -A %s '%s' ./finished_%s/PE2.%s > ./make_vcf_temp/PE2_%s_%s.%s" % (A, target, sample, format, i, j, format)
            #print cmd
            os.system(cmd)
            PE = " -1 ./make_vcf_temp/PE1_%s_%s.%s -2 ./make_vcf_temp/PE2_%s_%s.%s" % (i, j, format, i, j, format)
        if samples[sample]['SE']:
            cmd = "grep -A %s '%s' ./finished_%s/SE.%s > ./make_vcf_temp/SE_%s_%s.%s" % (A, target, sample, format, i, j, format)
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
            inf = open("./make_vcf_temp/tmp_%s_%s.vcf" % (i, j), 'r')
            for line in inf:
                if j == 0 and line[0] == "#" and line[0:9] != "##contig=":
                    header.append(line)
                elif line[0:9] == "##contig=":
                    contigs.append(line)
                elif line[0] != '#':
                    variants.append(line)
            #os.system("rm ./make_vcf_temp/tmp_%s_%s.vcf" % (i, j))
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


