#assembly_spades.py
#spades subprocess calls and exit status management

import os
import subprocess
from ARC import exceptions


params = {'t': 5, 'PE1': "phiX_S1_L001_R1_001_first15k.fastq", 'PE2': "phiX_S1_L001_R2_001_first15k.fastq", 'o': "spades_test"}

# 'PE1': "PE1.fastq", 'PE2': "PE2.fastq", 'SE': "SE.fastq", 

#'PE1': "phiX_S1_L001_R1_001_first15k.fastq"

print params

def RunSpades(params):

	"""
	Several arguments can be passed to spades.py: -1 [PE1], -2 [PE2], -s [SE], and -o [target_dir]
	"""
	#Check that required params are available	
	if not ('PE1' in params and 'PE2' in params) or ('SE' in params):
		raise exceptions.FatalException('Missing params in RunSpades.')

	#Check that the files actually exist
	if 'PE1' in params and 'PE2' in params and not(os.path.exists(params['PE1']) and os.path.exists(params['PE2'])):
		raise exceptioins.FatalException('Missing PE files in RunSpades.')
		
	if 'SE' in params and not(os.path.exists(params['SE'])):
		raise exceptioins.FatalException('Missing SE files in RunSpades.')
	
	#Build args for assembler call		
	args = ['spades.py']
	if 'PE1' in params and 'PE2' in params:
		args += ['-1', params['PE1'], '-2', params['PE2']]
		
	if 'SE' in params:
		args += ['-s', params['SE']]
	
	args += ['-o', params['target_dir', '-t', '1']]
	
	if 'verbose' in params:
		out = open(os.path.join(params['target_dir'], "assembly.log"), 'w')
		
	else:
		out = open(os.devnull, 'w')

	ret = subprocess.call(args, stderr=out, stdout=out)
	out.close()
	
	if ret != 0:
		raise exceptions.RerunnableError("Assembly failed")