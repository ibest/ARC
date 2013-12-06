import os
import re
from Bio import SeqIO
import argparse
from multiprocessing import Pool
import subprocess
import time

config=open('./ARC_config.txt', 'r').read().strip().split('\n')

for line in config:
	if line.startswith('# reference'):
		targetpath=line.split('=')[1]

parser = argparse.ArgumentParser()
parser.add_argument('-t', action='store', dest='targets', default=targetpath,
                   help='targets file for mapping (uses path specified in config file by default')
parser.add_argument('-p', action='store', dest='processes', default=7,
                   help='Number of processes to use')
args = parser.parse_args()
targets=args.targets
processes = int(args.processes)

def preprocessing(blat):
	blattemp=[]
	blatfinal=[]
	#check blat list for incorrect hits (where target in name and target blat mapped to don't match) and ditch, else keep
	for line in blat:
		if len(line[13].split('_:_'))>1:
			if line[9].split('_:_')[1]==line[13].split('_:_')[1]:
				blattemp.append(line)
		elif len(line[13].split('_:_'))==1:
			if line[9].split('_:_')[1]==line[13]:
				blattemp.append(line)	
	#if multiple hits, this now takes only the occurrence with highest percent identity (matches that aren't repeats/(size of query-extension left or right most end of target))
	blatnames=[]
	for line in blattemp:
		blatnames.append(line[9])
	for j in list(set(blatnames)):
		entry=[]
		id=[]
		for line in blattemp:
			qstart=int(line[11])+1
			tstart=int(line[15])+1	
			tend=int(line[16])
			tsize=int(line[14])
			qsize=int(line[10])
			qend=int(line[12])
			if j==line[9]:
				entry.append(line)
				hsubsize=0
				if tstart==1 and qstart>1:
					hsubsize=hsubsize+(qstart-1)
				if tend==tsize and qend<qsize:
					hsubsize=hsubsize+(qsize-qend)
				id.append(float(line[0])/(float(line[10])-hsubsize))
		blatfinal.append(entry[id.index(max(id))])	
	return blatfinal

#cigar is relative to whichever strand the read mapped to (left to right for forward strand, right to left for reverse strand)
def makecigar(blatline):
	qsize=int(blatline[10])
	tsize=int(blatline[14])
	blocksizes=[int(m) for m in blatline[18].split(",")[0:-1]]
	qstarts=[int(m)+1 for m in blatline[19].split(",")[0:-1]]
	tstarts=[int(m)+1 for m in blatline[20].split(",")[0:-1]]
	qstart=int(blatline[11])+1
	tstart=int(blatline[15])+1	
	qend=int(blatline[12])
	tend=int(blatline[16])
	if blatline[8]=="+":
		sqstart=qstart
		sqend=qend
	elif blatline[8]=="-":
		sqstart=(qsize-qend)+1
		sqend=(qsize-qstart)+1	
	cigar=[]
	cigarstring=str()
	if sqstart > 1:
		L=sqstart-1
		cigar.append(str(L)+"S")
	if len(blocksizes)==1:
		cigar.append(str(blocksizes[0])+"M")
	elif len(blocksizes) > 1:
		for i in range(len(blocksizes[0:-1])):
			cigar.append(str(blocksizes[i])+"M")
			if qstarts[i]+blocksizes[i]!=qstarts[i+1]:
				L=qstarts[i+1]-(qstarts[i]+blocksizes[i])
				cigar.append(str(L)+"I")
			if tstarts[i]+blocksizes[i]!=tstarts[i+1]:
				L=tstarts[i+1]-(tstarts[i]+blocksizes[i])
				cigar.append(str(L)+"D")
		cigar.append(str(blocksizes[-1])+"M")
	if qsize > sqend:
		L=qsize-sqend
		cigar.append(str(L)+"S")
	for i in cigar:
		cigarstring=cigarstring+i
	return cigarstring

#start and end are relative to forward strand
def findstart(blatline):
	fqstart=int(blatline[11])+1
	ftstart=int(blatline[15])+1	
	fstart=ftstart-fqstart+1
	return fstart
	
def findend(blatline):
	fqend=int(blatline[12])
	ftend=int(blatline[16])
	fqsize=int(blatline[10])
	fend=ftend+(fqsize-fqend)
	return fend

def outputseqs(blatline, originalcontigs):
	for seq_record in originalcontigs:
		if seq_record.id==blatline[9]:
                       seq=seq_record.seq
			#if blatline[8]=="+":
			#	seq=seq_record.seq
			#elif blatline[8]=="-":
			#	seq=(seq_record.seq).reverse_complement()
	return seq		
		
def process_blats(kk):
	print("processing "+samples[kk])
	originalcontigs=list(SeqIO.parse("./finished_"+samples[kk]+"/contigs.fasta", "fasta"))#this needs to be more flexible
	if len(originalcontigs)==0:
		print("contigs.fasta file empty")
	elif len(originalcontigs)>0:
		os.system(str("blat -noHead "+targets+" ./finished_"+samples[kk]+"/contigs.fasta ./finished_"+samples[kk]+"/blat_"+samples[kk]+".psl"))
		f=open(str("./finished_"+samples[kk]+"/blat_"+samples[kk]+".psl")).readlines()
		outfile=open(str("./finished_"+samples[kk]+"/processedcontigs.fa"),"w")
		blat=[]
		for line in f:
			blat.append(line.strip().split("\t"))
		for line in preprocessing(blat):
			#print line[9]
			#print line[13]
			#print findstart(line)
			#print findend(line)
			#print line[8]
			#print makecigar(line)
			#print outputseqs(line, originalcontigs)

			outfile.write(str(">"+line[9]+" "+line[13]+" "+str(findstart(line))+" "+str(findend(line))+" "+line[8]+" "+makecigar(line)+"\n"+outputseqs(line, originalcontigs)+"\n"))
			#out = open('.
		outfile.close()

def check_status(results):
	unfinished = 0
	finished = 0
	for r in results:
		if not results[r].ready():
			unfinished += 1
		else:
			finished += 1
	print "Finished:", finished, " Unfinished:", unfinished
	return unfinished 

samples=[]

for line in config:
	if not line.startswith('#'):
		if not line.startswith('Sample_ID'):
			if len(line.split())>1:
				samples.append(line.split()[0])

samples=list(set(samples))

#multiprocessing wrapper
p = Pool(processes=processes, maxtasksperchild=1)
results = {}

for kk in range(len(samples)):
	results[samples[kk]] = p.apply_async(process_blats, (kk,))

print results.keys()
for r in results:
	print r


check_status(results)
allfinished = False
while not allfinished:
	time.sleep(5)
	np = check_status(results)
	if np == 0:
		allfinished = True

