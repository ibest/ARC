#!/usr/bin/env python
"""
A set of tests for each class/function.

These are meant to be run using the run_unit_tests.py in the test_data folder.

Note that some of these rely on previous tests having been run in order to have the necessary files setup
and available.
"""

import os
import sys
import pprint
import copy

os.system("rm -rf working_Sample1")
os.system("rm -rf working_Sample2")

lib_path = os.path.join(os.getcwd(), "../")
if lib_path not in sys.path:
    sys.path.insert(0, lib_path)

######### __init__.py tests ##########

from ARC.__init__ import read_config, setup


pp = pprint.PrettyPrinter(indent=2)
print "-------------Testing ARC.__init__ read_config()-------------"
config = read_config()
print "Config:"
pp.pprint(config)

print "\n\n-------------Testing ARC.__init__ setup()-------------"
setup(config)
pp.pprint(config)


######### Blat Mapper tests ##########
print "\n\n-------------Testing ARC.mapper MapperRunner()-------------"
from ARC.mapper import MapperRunner

#Set up params for calling the mapper:
blat_params = []
for sample in config['Samples']:
    s = config['Samples'][sample]
    params = {
        'reference': config['reference'],
        'numcycles': config['numcycles']
    }
    if 'PE1' in s and 'PE2' in s:
        params['PE1'] = s['PE1']
        params['PE2'] = s['PE2']
    if 'SE' in s:
        params['SE'] = s['SE']
    params['working_dir'] = s['working_dir']
    params['sample'] = sample
    params['mapper'] = config['mapper']
    params['verbose'] = True
    params['format'] = config['format']
    params['assembler'] = config['assembler']
    blat_params.append(params)

#Instantiate MapperRunner class objects:
blat_mapper1 = MapperRunner(blat_params[0])
blat_mapper2 = MapperRunner(blat_params[1])

#Test to_dict()
print "\n\n-------------Testing MapperRunner class objects: to_dict()-------------"
pp.pprint(blat_mapper1.to_dict())
pp.pprint(blat_mapper2.to_dict())

#Test start() which runs:
# start() --> run_blat() --> SAM_to_dict(), write_dict()
print "\n\n-------------Testing MapperRunner: start()-------------"
blat_mapper1.start()
blat_mapper2.start()










#from ARC import splitter
#from ARC import assembler
