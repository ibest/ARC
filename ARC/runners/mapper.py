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
from collections import Counter
from copy import deepcopy
from Bio import SeqIO
from ARC import exceptions
from ARC import logger
from ARC.runners import BaseRunner
from ARC.runners import AssemblyRunner
from ARC.runners import Finisher


class MapperRunner(BaseRunner):
    """
    This calss handles mapping jobs, as well as converting map results into a text version of a dict.
    required params:
        PE1, PE2, SE, format, mapper, numcycles, reference, sample, verbose, working_dir
    params added:
        mapping_dict
    """
    def setup(self):
        if not 'sample' and 'reference' and 'working_dir' in self.params:
            raise exceptions.FatalError('Missing params in mapper.')

        pe_run = 'assembly_PE1' and 'assembly_PE2' in self.params
        se_run = 'assembly_SE' in self.params
        if not (pe_run or se_run):
            raise exceptions.FatalError('Missing params in mapper.')

        pe_one_path = os.path.exists(self.params['assembly_PE1'])
        pe_two_path = os.path.exists(self.params['assembly_PE2'])
        if pe_run and not (pe_one_path and pe_two_path):
            raise exceptions.FatalException('Missing PE files in assembler.')

        se_path = os.path.exists(self.params['assembly_SE'])
        if se_run and not se_path:
            raise exceptions.FatalException('Missing SE file in assembler.')

        if os.path.exists(self.params['reference']) is False:
            raise exceptions.FatalError("Missing reference file for mapping")

    def execute(self):
        if not('mapper' in self.params):
            raise exceptions.FatalError("mapper not defined in params")

        if self.params['mapper'] == 'bowtie2':
            logger.info("Sample: %s Running bowtie2." % self.params['sample'])
            self.run_bowtie2()
        elif self.params['mapper'] == 'blat':
            logger.info("Sample: %s Running blat." % self.params['sample'])
            self.run_blat()

        #Mapping is done, run splitreads:
        logger.info("Sample: %s Running splitreads." % self.params['sample'])
        self.splitreads()

    def run_bowtie2(self):
        """
        Builds idx and runs bowtie2 -I 0 -X 1500 --local
        Expects params:
            sample, target, reference, working_dir, PE1 and PE2 and/or SE
        """
        #Make idx directory
        try:
            working_dir = self.params['working_dir']
            idx_dir = os.path.realpath(os.path.join(working_dir, 'idx'))
            os.mkdir(idx_dir)
        except Exception as exc:
            txt = "Sample: %s Error creating working directory." % (
                self.params['sample']) + '\n\t' + str(exc)
            raise exceptions.FatalError(txt)

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
        n_bowtieprocs = int(round(max(float(self.params['nprocs']) / len(self.params['Samples']), 1)))
        #print "Number of bowtie2 procs:", n_bowtieprocs
        #args = ['nice', '-n', '19', 'bowtie2', '-I', '0', '-X', '1500',
        #    '--local', '-p', self.params['nprocs'], '-x', base]
        #args = ['nice', '-n', '19', 'bowtie2', '-I', '0', '-X', '1500',
        #   '--local', '-p', '1', '-x', base]
        args = [
            'bowtie2',
            '-I', '0',
            '-X', '1500',
            '--local',
            '-p', str(n_bowtieprocs),
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
            raise exceptions.FatalError("Missing SAM file")
        try:
            inf = open(filename, 'r')
        except Exception as exc:
            txt = "Failed to open SAM file %s" % filename
            txt += '\n\t' + str(exc)
            raise exceptions.FatalError(txt)

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
        logger.info("Sample: %s, Processed %s lines in %s seconds." % (
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

        logger.info("Sample: %s, Processed %s lines in %s seconds." % (
            self.params['sample'], i, time.time() - startT))

        return read_map

    def splitreads(self):
        """
            Split reads and then kick off assemblies once the reads are split
            for a target, use safe_targets for names
        """
        self.params['iteration'] += 1

        checker_params = deepcopy(self.params)
        checker_params['targets'] = {}
        iteration = self.params['iteration']

        if 'PE1' in self.params and 'PE2' in self.params:
            idx_PE1 = SeqIO.index_db(
                os.path.join(
                    self.params['working_dir'],
                    "PE1.idx"),
                key_function=lambda x: x.split("/")[0])
            idx_PE2 = SeqIO.index_db(
                os.path.join(
                    self.params['working_dir'],
                    "PE2.idx"),
                key_function=lambda x: x.split("/")[0])

        if 'SE' in self.params:
            idx_SE = SeqIO.index_db(
                os.path.join(
                    self.params['working_dir'],
                    "SE.idx"),
                key_function=lambda x: x.split("/")[0])

        if 'readcounts' not in self.params:
            checker_params['readcounts'] = {}

        assids = []
        for target in checker_params['mapping_dict']:
            startT = time.time()
            target_dir = os.path.join(
                checker_params['working_dir'],
                checker_params['safe_targets'][target])

            if target not in checker_params['readcounts']:
                checker_params['readcounts'][target] = Counter()
            os.mkdir(target_dir)
            assembly_params = deepcopy(self.params)
            assembly_params['target'] = target
            assembly_params['target_dir'] = target_dir
            reads = self.params['mapping_dict'][target]

            # track how many total reads were added for this cycle
            checker_params['readcounts'][target][iteration] = len(reads)
            SEs = PEs = 0

            if 'PE1' in self.params and 'PE2' in self.params:
                outf_PE1 = open(os.path.join(
                    target_dir,
                    "PE1." + self.params['format']), 'w')
                outf_PE2 = open(os.path.join(
                    target_dir,
                    "PE2." + self.params['format']), 'w')
                #idx_PE1 = self.params['indexes'][sample]['PE1']
                #idx_PE2 = self.params['indexes'][sample]['PE2']

            if 'SE' in self.params:
                outf_SE = open(os.path.join(
                    target_dir,
                    "SE." + self.params['format']), 'w')
                #idx_SE = self.params['indexes'][sample]['SE']

            for readID in reads:
                if 'PE1' in self.params and readID in idx_PE1:
                    read1 = idx_PE1[readID]
                    read2 = idx_PE2[readID]
                    new_readID = readID.replace(":", "_") + ":0:0:0:0#0/"
                    read1.id = read1.name = new_readID + "1"
                    read2.id = read2.name = new_readID + "2"
                    SeqIO.write(read1, outf_PE1, self.params['format'])
                    SeqIO.write(read2, outf_PE2, self.params['format'])
                    PEs += 1
                elif 'SE' in self.params and readID in idx_SE:
                    read1 = idx_SE[readID]
                    read1.id = read1.name = readID.replace(":", "_") + ":0:0:0:0#0/"
                    SeqIO.write(read1, outf_SE, self.params['format'])
                    SEs += 1
            if 'PE1' in self.params and 'PE2' in self.params:
                outf_PE1.close()
                outf_PE2.close()
            if 'SE' in self.params:
                outf_SE.close()

            # Properly handle the case where no reads ended up mapping for the
            # PE or SE inputs:
            if PEs > 0:
                assembly_params['assembly_PE1'] = os.path.join(
                    target_dir,
                    "PE1." + self.params['format'])
                assembly_params['assembly_PE2'] = os.path.join(
                    target_dir,
                    "PE2." + self.params['format'])

            if SEs > 0:
                assembly_params['assembly_SE'] = os.path.join(
                    target_dir, "SE." + self.params['format'])

            # All reads have been written at this point, submit the
            # assemblies
            msg = "Sample: %s target: %s iteration: %s Split %s" % (
                self.params['sample'],
                target,
                self.params['iteration'],
                len(reads))
            msg += "reads in %s seconds" % (time.time() - startT)
            logger.info(msg)

            # Only add an assembly job and AssemblyChecker target if there
            # are >0 reads:
            if PEs + SEs > 0:
                checker_params['targets'][target_dir] = False
                id = self.submit(
                    AssemblyRunner,
                    procs=1,  # This can now be changed from params!!!!!
                    params=assembly_params)
                assids.append(id)

        logger.info("------------------------------------")
        logger.info("Sample: %s Iteration %s of numcycles %s" % (
            checker_params['sample'],
            checker_params['iteration'],
            checker_params['numcycles']))
        logger.info("------------------------------------")
        if 'PE1' in self.params and 'PE2' in self.params:
            idx_PE1.close()
            idx_PE2.close()
        if 'SE' in self.params:
            idx_SE.close()

        # Kick off the finisher:
        self.submit(
            Finisher,
            procs=1,
            deps=assids,
            params=self.params)
