#assembly.py
#modular python script to pass arguments to the assembler and manage output


from optparse import OptionParser
import os, re, math


params = {'PE1': 'PE1.fastq', 'PE2': 'PE2.fastq', 'SE': 'SE.fastq', 't': 5}

def run_spades(params):
	"""
	Several arguments can be passed to spades.py: -1 [PE1], -2 [PE2], -s [SE], -t [# of threads]
	"""
	
	


parser = OptionParser(usage = "usage: python %prog [-f FILENAME] [-o output_file_name]", version = "%prog 0.1")

parser.add_option("-t", "--threads", dest="threads", help="Number of threads to run SPADES on.")






parser.add_option("-f", "--file", dest="filename", help="Path to the *.codeml file to be parsed")
parser.add_option("-o", "--outfile", dest="outfile", help="Output File")
(options, args) = parser.parse_args()
if not options.filename and not options.outfile:
	parser.error("You must specify a codeml output file and a file to write results to.")



