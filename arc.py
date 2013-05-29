import glob
from Bio import SeqIO
from collections import Counter
import os
from multiprocessing import Pool
import time

#######                                                                                                                                                                                 
"""
Input:

Reads are split by barcode:  lane5_******_PE1_adp_qual.fastq and lane5_******_PE2_adp_qual.fastq
/mnt/lfs/dmorales/Capture/Raw_data/

Blat resulst are stored:
/mnt/lfs/dmorales/Capture/Blat:   geneblat_lane5_CAAAAG_PE2_adp_qual.psl, and   geneblat_lane5_CAAAAG_PE1_adp_qual.psl

As of Version 3: Current filename format is listed below, but the directory structure remains the same.

"""
#######
def process_barcode(barcode, fastq1, fastq2, fastq3, psl1, psl2, psl3):
    #This dictionary will use readID as key, and put geneID in the values in a list:
    os.makedirs(barcode, 0755)
    os.chdir(barcode)
    print "Processing PSL+fastq for %s" % barcode
    read_map = {}

    #open the first psl file:
    inf = open(psl1, 'r')

    #skip the first 5 header lines
    for i in range(5):
        inf.readline()

    print "Processing PSL1: %s" % barcode
    for line in inf:
        line = line.strip().split("\t")
        readid = line[9].split("#")[0]
        gene = line[13]
        if readid in read_map:
            if gene not in read_map[readid]:
                read_map[readid].append(gene)
        else:
            read_map[readid] = list()
            read_map[readid].append(gene)                                                                                                                                               

    #open the second psl file:
    inf = open(psl2, 'r')

    #skip the first 5 header lines
    for i in range(5):
        inf.readline()

    print "Processing PSL2: %s" % barcode
    for line in inf:
        line = line.strip().split("\t")
        readid = line[9].split("#")[0]
        gene = line[13]
        if readid in read_map:
            if gene not in read_map[readid]:
                read_map[readid].append(gene)
        else:
            read_map[readid] = list()
            read_map[readid].append(gene)

    #open the third psl file (the single end reads)
    inf = open(psl3, 'r')

    #skip the first five header lines
    for i in range(5):
        inf.readline()

    print "Processing PSL3: %s" % barcode
    for line in inf:
        line = line.strip().split("\t")
        readid = line[9].split("#")[0]
        gene = line[13]
        if readid in read_map:
            if gene not in read_map[readid]:
                read_map[readid].append(gene)
        else:
            read_map[readid] = list()
            read_map[readid].append(gene)


    #Parse fastq files, output reads to appropriate "target" files
    iter1 = SeqIO.parse(fastq1, 'fastq')
    iter2 = SeqIO.parse(fastq2, 'fastq')
    outf = {}
    i = 0

    all_genes = Counter()

    print "Writing reads for barcode: %s" % barcode
    try:
        while 1:
            read1 = iter1.next()
            key = read1.id.split("#")[0]
            #ks = key.split(":")
            #read1.id = read1.name = "_".join([ks[0],ks[2],ks[4],ks[5],ks[6]]) + ":0:0:0:0#0/1"

            read2 = iter2.next()
            #read2.id = "_".join([ks[0],ks[2],ks[4],ks[5],ks[6]]) + ":0:0:0:0#0/2"
            
            if key in read_map:
                genes = read_map[key]
            else:
                continue
                
            for gene in genes:
                if gene not in outf:  
                    #open new file
                    outf[gene] = [open("%s_PE1.fastq" % gene, 'w'), open("%s_PE2.fastq" % gene, 'w')]
                SeqIO.write(read1, outf[gene][0], "fastq")
                SeqIO.write(read2, outf[gene][1], "fastq")
                all_genes[gene] += 1
                
            i += 1
            if i % 10000 == 0:
                print "PE process analyzed %s records" % i
    except StopIteration:
        pass
    finally:
        for key in outf:
            outf[key][0].close()
            outf[key][1].close()

    print "Finished processing PE.  Starting SE processing."
    #iterate over reads in SE file and write out
    iter3 = SeqIO.parse(fastq3, 'fastq')
    outfse = {}
    i = 0
    all_genesse = Counter()
    
    try:
        while 1:                
            read3 = iter3.next()
            keyse = read3.id.split("#")[0]
            #ksse = keyse.split(":")
            #read3.id = read3.name = "_".join([ksse[0],ksse[2],ksse[4],ksse[5],ksse[6]]) + ":0:0:0:0"
            if keyse in read_map:
                genese = read_map[keyse]
            else:
                continue
            for gene in genese:
                if gene not in outfse:
                    outfse[gene] = [open("%s_SE.fastq" % gene, 'w')]
                SeqIO.write(read3, outfse[gene][0], "fastq")
                all_genesse[gene] += 1
            i += 1
            if i % 10000 == 0:
                print "SE process analyzed %s records" % i
    except StopIteration:
        pass
    finally:
        for keyse in outfse:
            outfse[keyse][0].close()
  
    os.chdir("..")
    
    ## combine all_genes and all_genesse into a counter (to ensure a single entry per gene
    combined_genes = Counter()
    for gene in all_genes:
      combined_genes[gene] += 1
    for gene in all_genesse:
      combined_genes[gene] += 1    
    return(combined_genes)


def runassembly(barcode, gene):
    os.chdir("./%s" % barcode)
    ## Do assembly:
    finished = os.path.realpath('./%s/454AllContigs.fna' % gene )
    i=0
    # There is a potential bug here, if PE doesn't exist and SE does, runAssembly will probably complain
    # (maybe fail?)
    while not os.path.isfile(finished):
        print "Starting assembly: %s:%s, iteration:%s" % (barcode,gene, i)
        os.system("runAssembly -noace -cpu 2 -force -o %s %s_PE1.fastq %s_PE2.fastq %s_SE.fastq &>> %s_run.log"  % (gene, gene, gene, gene, gene))
    print "Assembly: %s:%s complete on iteration %s" % (barcode,gene, i)
    return()

def check_status(results):
    unfinished = 0
    finished = 0
    for r in results:
        if not r[2].ready():
            print r[0], ":", r[1], ":", r[2].ready()
            unfinished += 1
        else:
            finished += 1
    print "------Process Status Report--------"
    print "Finished: %s, still processing: %s" % (finished, unfinished)
    print "-----------------------------------"
    return(unfinished)


########## PRINCIPAL EXECUTION ##########

#Brice's .psl format: 151_PE1.mt.fasta.psl

#Brice's Raw_data format: 151_PE1.fastq


#psl string processing for Brice's files.  Get list of barcodes.
psl1_list = glob.glob("./Blat/*_PE1*.psl")
barcodes = []
for psl in psl1_list:
    temp = psl.split("_")[0]
    barcodes.append(temp.split("/")[2])

results = []
pool = Pool(processes=16, maxtasksperchild=1)  # This pool is used for assemblies
pool2 = Pool(processes=1, maxtasksperchild=1)  # This pool is used for extracting reads from fastq
for barcode in barcodes:
    print "Processing barcode ", barcode
    psl1 = os.path.realpath("./Blat/%s_PE1.psl" % barcode)
    psl2 = os.path.realpath("./Blat/%s_PE2.psl" % barcode)
    psl3 = os.path.realpath("./Blat/%s_SE.psl" % barcode)
    fastq1 = os.path.realpath("./Raw_data/%s_seqycleaned_PE1.fastq" % barcode)
    fastq2 = os.path.realpath("./Raw_data/%s_seqycleaned_PE2.fastq" % barcode)
    fastq3 = os.path.realpath("./Raw_data/%s_seqycleaned_SE.fastq" % barcode)
    #split the file using pool2
    bc_Result = pool2.apply_async(process_barcode, (barcode, fastq1, fastq2, fastq3, psl1, psl2, psl3,))
    bc_Result.wait()  # This is a non-blocking wait.. meaning the other pool can continue to start workers
    all_genes = bc_Result.get()
    #Add assemblies to the pool of workers:
    i=0
    for gene in all_genes:
        results.append((barcode, gene, pool.apply_async(runassembly, (barcode, gene,))))
        i+=1
    print "Total genes: %s, Total jobs added: %s." % (len(all_genes), i)
    check_status(results)

### At this point all processes have been added to the job queue.. so we just have to wait for all of them to finish:
allfinished = False
while not allfinished:
    time.sleep(10)
    np = check_status(results)
    if np == 0:
        allfinished = True


"""
The following 5 processes never finished
------Process Status Report--------
Finished: 8153, still processing: 5
-----------------------------------
GAGTGG : gene32182 : False       <------ this actually finished properly
CACTCA : gene30526 : False       <------ this actually finished properly  
CACTCA : gene03608 : False       <------ this actually finished properly  
ATTCCT : gene15629 : False       <------ this actually finished properly  
TAGCTT : gene23736 : False       <------ this actually finished properly
"""