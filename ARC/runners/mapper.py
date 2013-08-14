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
import os
from ARC import logger
from ARC.runners import ProcessBase
from ARC.runners import Splitter
from ARC import FatalError


class Mapper(ProcessBase):
    """
    This calss handles mapping jobs, as well as converting map results into a text version of a dict.
    required params:
        PE1, PE2, SE, format, mapper, numcycles, reference, sample, verbose, working_dir
    params added:
        mapping_dict
    """
    def setup(self):
        required_params = 'sample' and 'reference' and 'working_dir' in self.params
        if not required_params:
            raise FatalError('Missing sample, reference or working directory params in mapper.')

        pe_run = 'PE1' and 'PE2' in self.params
        se_run = 'SE' in self.params
        if not (pe_run or se_run):
            raise FatalError('Missing assembly params in mapper.')

        pe_one_path = os.path.exists(self.params['PE1'])
        pe_two_path = os.path.exists(self.params['PE2'])
        if pe_run and not (pe_one_path and pe_two_path):
            raise FatalError('Missing PE files in assembler.')

        se_path = os.path.exists(self.params['SE'])
        if se_run and not se_path:
            raise FatalError('Missing SE file in assembler.')

        if os.path.exists(self.params['reference']) is False:
            raise FatalError("Missing reference file for mapping")

    def execute(self):
        if not('mapper' in self.params):
            raise FatalError("mapper not defined in params")

        if self.params['mapper'] == 'bowtie2':
            self.log("Sample: %s Running bowtie2." % self.params['sample'])
            self.run_bowtie2()
        elif self.params['mapper'] == 'blat':
            self.log("Sample: %s Running blat." % self.params['sample'])
            self.run_blat()

        #Mapping is done, run splitreads:
        self.submit(
            Splitter,
            procs=1,
            params=self.params.copy())

    def run_bowtie2(self):
        """
        Builds idx and runs bowtie2 -I 0 -X 1500 --local
        Expects params:
            sample, target, reference, working_dir, PE1 and PE2 and/or SE
        """
        #Make idx directory
        working_dir = self.params['working_dir']

        try:
            logger.debug("Running bowtie in %s" % (working_dir))
            idx_dir = os.path.realpath(os.path.join(working_dir, 'idx'))
            os.mkdir(idx_dir)
        except Exception as exc:
            txt = "Sample: %s Error creating working directory." % (
                self.params['sample']) + '\n\t' + str(exc)
            raise FatalError(txt)

        #Build index
        base = os.path.join(idx_dir, 'idx')

        args = [
            'bowtie2-build',
            '-f',
            self.params['reference'], base]
        self.shell(
            args,
            description="Bowtie2 build (Sample %s)" % (
                self.params['sample']),
            logfile='assembly.log',
            working_dir=working_dir,
            verbose=self.params['verbose'])

        #Do bowtie2 mapping:
        args = [
            'bowtie2',
            '-I', '0',
            '-X', '1500',
            '--local',
            '-p', str(self.params['mapping_procs']),
            '-x', base]
        if self.params['format'] == 'fasta':
            args += ['-f']
        if 'PE1' in self.params and 'PE2' in self.params:
            args += ['-1', self.params['PE1'], '-2', self.params['PE2']]
        if 'SE' in self.params:
            args += ['-U', self.params['SE']]
        args += ['-S', os.path.join(working_dir, 'mapping.sam')]

        self.shell(
            args,
            description="Bowtie2 mapping (Sample %s)" % (
                self.params['sample']),
            logfile='assembly.log',
            working_dir=working_dir,
            verbose=self.params['verbose'])

        #Extract the SAM to a dict
        self.params['mapping_dict'] = self.SAM_to_dict(
            os.path.join(working_dir, 'mapping.sam'))
        #clean up intermediary files:
        os.remove(os.path.join(working_dir, 'mapping.sam'))
        os.system("rm -rf %s" % idx_dir)

    def run_blat(self):
        #Blat doesn't need an index
        working_dir = self.params['working_dir']

        #Build a temporary txt file with all of the reads:
        allreads_outf = open(os.path.join(working_dir, 'reads.txt'), 'w')
        if 'PE1' in self.params and 'PE2' in self.params:
            allreads_outf.write(self.params['PE1'] + '\n')
            allreads_outf.write(self.params['PE2'] + '\n')
        if 'SE' in self.params:
            allreads_outf.write(self.params['SE'] + '\n')
        allreads_outf.close()

        #Do blat mapping
        args = [
            'blat',
            self.params['reference'],
            os.path.join(working_dir, 'reads.txt')]

        if 'format' in self.params and self.params['format'] == 'fastq':
            args.append('-fastq')
        if 'fastmap' in self.params:
            args.append('-fastMap')

        args.append(os.path.join(working_dir, 'mapping.psl'))

        self.shell(
            args,
            description="Blat mapping (Sample %s)" % (self.params['sample']),
            logfile='assembly.log',
            working_dir=working_dir,
            verbose=self.params['verbose'])

        #Extract the PSL to a dict
        self.params['mapping_dict'] = self.PSL_to_dict(
            os.path.join(working_dir, 'mapping.psl'))

    def SAM_to_dict(self, filename):
        """
            Read a SAM file to a mapping dict and return it
        """
        #Check for necessary files:
        if os.path.exists(filename) is False:
            raise FatalError("Missing SAM file")
        try:
            inf = open(filename, 'r')
        except Exception as exc:
            txt = "Failed to open SAM file %s" % filename
            txt += '\n\t' + str(exc)
            raise FatalError(txt)

        read_map = {}  # target:{read} dictionary of dictionaries
        i = 0
        startT = time.time()

        for l in inf:
            i += 1
            if l[0] != "@":  # skip header lines
                l2 = l.strip().split()
                if l2[2] == "*":  # skip unmapped
                    continue
                readid = l2[0].split("/")[0]
                target = l2[2]
                #handle references built using assembled contigs:
                if len(target.split("_:_")) > 1:
                    target = target.split("_:_")[1]
                if target not in read_map:
                    read_map[target] = {}
                read_map[target][readid] = 1

        #Report total time:
        self.log("Sample: %s, Processed %s lines in %s seconds." % (
            self.params['sample'], i, time.time() - startT))
        return read_map

    def PSL_to_dict(self, filename):
        """Process a PSL file to the dict format """
        try:
            inf = open(filename, 'r')
        except Exception as inst:
            if type(inst) == IOError:
                logger.error("Failed to open mapping dictionary %s." % (
                    filename))
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
            readid = l2[9].split("/")[0]  # remove unique part of PE reads
            target = l2[13]
            #handle references built using assembled contigs:
            if len(target.split("_:_")) > 1:
                target = target.split("_:_")[1]
            if target not in read_map:
                read_map[target] = {}
            read_map[target][readid] = 1

        self.log("Sample: %s, Processed %s lines in %s seconds." % (
            self.params['sample'], i, time.time() - startT))

        return read_map


