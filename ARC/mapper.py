#!/usr/bin/env python

# Copyright 2013, Institute for Bioninformatics and Evolutionary Studies
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import time
from ARC import logger
import subprocess
import os
from ARC import exceptions


class MapperRunner:
    def __init__(self, params):
        self.params = params

    def start(self, params):
        #print "Running the mapper"
        if not('mapper' in params):
            raise exceptions.FatalException("mapper not defined in params")
        if params['mapper'] == 'bowtie2':
            self.run_bowtie2(params)
        if params['mapper'] == 'blat':
            self.run_blat(params)

    def run_bowtie2(self, params):
        """
        Expects params:
            sample - required
            target - required
            PE1 and PE2 or SE
        """
        #Check for necessary params:
        if not ('sample' in params and 'reference' in params and (('PE1' in params and 'PE2' in params) or 'SE' in params)):
            raise exceptions.FatalException('Missing params in run_bowtie2.')
        #Check for necessary files:
        if os.path.exists(params['reference']) is False:
            raise exceptions.FatalException("Missing reference file for mapping")
        if 'PE1' in params and 'PE2' in params:
            if not (os.path.exists(params['PE1']) and os.path.exists(params['PE2'])):
                raise exceptions.FatalException("One or both PE files can not be found for mapping.")
        if 'SE' in params:
            if not os.path.exists(params['SE']):
                raise exceptions.FatalException("SE file cannot be found.")

        #Make temporary working directory and idx directory
        try:
            working_dir = os.path.realpath('_'.join(['tmp', params['sample']]))
            os.mkdir(working_dir)
            idx_dir = os.path.realpath(os.path.join(working_dir, 'idx'))
            os.mkdir(idx_dir)
            params['working_dir'] = working_dir
        except Exception as exc:
            txt = "Error creating working directory for Sample: %s" % (params['sample']) + '\n\t' + str(exc)
            raise exceptions.FatalException(txt)

        #Check whether to log to temporary file, or default to os.devnull
        if ''

        base = os.path.join(idx_dir, 'idx')
        ret = subprocess.call(['bowtie2-build', '-f', params['target'], base], stdout=open(os.devnull, 'w'))
        if ret != 0:
            raise exceptions.FatalException("Error creating bowtie2 index for Sample: %s" % params['sample'])

        #Carry out mapping:
        args = ['bowtie2', '--local', '-x', base]
        if 'PE1' in params and 'PE2' in params:
            args += ['-1', params['PE1'], '-2', params['PE2']]
        if 'SE' in params:
            args += ['-U', params['SE']]
        args += ['-S', os.path.join(working_dir, 'mapping.sam')]
        ret = subprocess.call(args, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))




        txt = "bowtie2 -p %s --local -x "
        bowtie2 -p 15 --local -x ./idx/IsoTigs -1 ../01-Cleaned_Reads/Male_lane5/Male_PE1.fastq ../01-Cleaned_Reads/Male_lane5/Male_PE2.fastq -S male_IsoTig_hits.sam &



    def run_blat(self, params):
        "blat ./idx/454IsoTigs_renamed.fasta reads.txt -minScore=60 -minIdentity=95 -dots=1000000 -fastq Blat_IsoTig_hits.psl &"








    def SAM_to_dict(self, filename):
        """ Read a SAM file to a mapping dict and return it """
        try:
            inf = open(filename, 'r')
        except IOError:
            logger.error("Failed to open mapping dictionary %s." % filename)
            raise

        read_map = {}  # target:{read} dictionary of dictionaries
        i = 0
        startT = time.time()

        for l in inf:
            i += 1
            if l[0] != "@":
                l2 = l.strip().split()
                if l2[2] == "*":  # skip unmapped
                    continue
                readid = l2[0]
                target = l2[2].split("_")[0].split("=")[1]
                if target not in read_map:
                    read_map[target] = {}
                read_map[target][readid] = 1
        #Report total time:
        logger.info("Processed %s lines in %s seconds." % (i, time.time() - startT))
        return read_map

    def PSL_to_dict(self, filename):
        try:
            inf = open(filename, 'r')
        except Exception as inst:
            if type(inst) == IOError:
                logger.error("Failed to open mapping dictionary %s." % filename)
            raise inst
        read_map = {}
        i = 0
        startT = time.time()

        psl_header = False

        for l in inf:
            i += 1
            #Check for PSL header and skip 5 lines if it exists
            if i == 1 and l.split()[0] == 'psLayout':
                psl_header = True
            if psl_header and i <= 5:
                continue
            l2 = l.strip().split("\t")
            readid = l2[9]
            target = l2[13]
            if target not in read_map:
                read_map[target] = {}
            read_map[target][readid] = 1
        logger.info("Processed %s lines in %s seconds." % (i, time.time() - startT))
        return read_map

    def write_dict(self, filename, read_map):
        """ Write a mapping dictionary to a file. """
        startT = time.time()
        outf = open(filename, 'w')
        for k in read_map.keys():
            outf.write(k + '\t' + ",".join(read_map[k].keys()) + '\n')
        outf.close()
        logger.info("Wrote all values to txt in %s seconds" % (time.time() - startT))

    def read_dict(self, filename):
        """ Read a mapping dictionary from a file """
        startT = time.time()
        try:
            inf = open(filename, 'r')
        except Exception as inst:
            if type(inst) == IOError:
                logger.error("Failed to open mapping dictionary %s." % filename)
            raise inst
        new_map = {}
        for l in inf:
            l2 = l.split('\t')
            l3 = l2[1].strip().split(",")
            new_map[l2[0]] = {}
            for k in l3:
                new_map[l2[0]][k] = 1
        logger.info("Read all values to txt in %s seconds" % (time.time() - startT))
        return new_map