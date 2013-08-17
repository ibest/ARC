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
from Bio import SeqIO
from collections import Counter
from ARC import logger
from ARC.runners import ProcessBase
from ARC.runners import Assembler
from ARC.runners import Finisher
from ARC import FatalError


class Splitter(ProcessBase):
    """
        Deprecated in favor of doing the splitting of reads post mapping as
        part of Mapper. This class handles splitting reads into fastq
        files and launching assemblies. Once all assemblies are launched, add
        a job to check that the assemblies have finished.
    """
    def setup(self):
        required_params = 'sample' and 'reference' and 'working_dir' in self.params
        if not required_params:
            raise FatalError('Missing sample, reference or working directory params in mapper.')

    def execute(self):
        self.log("Sample: %s Running splitreads." % self.params['sample'])
        self.splitreads()

    def splitreads(self):
        """
            Split reads and then kick off assemblies once the reads are split
            for a target, use safe_targets for names
        """
        self.params['iteration'] += 1
        self.params['targets'] = {}

        iteration = self.params['iteration']

        if 'readcounts' not in self.params:
            self.params['readcounts'] = {}

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

        #if 'readcounts' not in self.params:
        #    checker_params['readcounts'] = {}

        assids = []
        for target in self.params['mapping_dict']:
            startT = time.time()
            target_dir = os.path.join(
                self.params['working_dir'],
                self.globals['safe_targets'][target])
            os.mkdir(target_dir)

            assembly_params = {}
            assembly_params['target'] = target
            assembly_params['target_dir'] = target_dir
            assembly_params['iteration'] = iteration
            assembly_params['sample'] = self.params['sample']
            assembly_params['targets'] = self.params['targets']

            if target not in self.params['readcounts']:
                self.params['readcounts'][target] = Counter()

            reads = self.params['mapping_dict'][target]
            self.params['readcounts'][target][iteration] = len(reads)

            SEs = PEs = 0

            if 'PE1' in self.params and 'PE2' in self.params:
                outf_PE1 = open(os.path.join(
                    target_dir,
                    "PE1." + self.globals['format']), 'w')
                outf_PE2 = open(os.path.join(
                    target_dir,
                    "PE2." + self.globals['format']), 'w')

            if 'SE' in self.params:
                outf_SE = open(os.path.join(
                    target_dir,
                    "SE." + self.globals['format']), 'w')

            for readID in reads:
                if 'PE1' in self.params and readID in idx_PE1:
                    read1 = idx_PE1[readID]
                    read2 = idx_PE2[readID]
                    new_readID = readID.replace(":", "_") + ":0:0:0:0#0/"
                    read1.id = read1.name = new_readID + "1"
                    read2.id = read2.name = new_readID + "2"
                    SeqIO.write(read1, outf_PE1, self.globals['format'])
                    SeqIO.write(read2, outf_PE2, self.globals['format'])
                    PEs += 1
                elif 'SE' in self.params and readID in idx_SE:
                    read1 = idx_SE[readID]
                    read1.id = read1.name = readID.replace(":", "_") + ":0:0:0:0#0/"
                    SeqIO.write(read1, outf_SE, self.globals['format'])
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
                    "PE1." + self.globals['format'])
                assembly_params['assembly_PE2'] = os.path.join(
                    target_dir,
                    "PE2." + self.globals['format'])

            if SEs > 0:
                assembly_params['assembly_SE'] = os.path.join(
                    target_dir, "SE." + self.globals['format'])

            # All reads have been written at this point, submit the
            # assemblies
            msg = "Sample: %s target: %s iteration: %s Split %s" % (
                self.params['sample'],
                target,
                self.params['iteration'],
                len(reads))
            msg += " reads in %s seconds" % (time.time() - startT)
            self.log(msg)

            # Only add an assembly job and AssemblyChecker target if there
            # are >0 reads:
            if PEs + SEs > 0:
                self.params['targets'][target_dir] = False
                logger.debug("Submitting new assembly")
                job = self.submit(
                    Assembler,
                    procs=self.globals['assembly_procs'],
                    params=assembly_params)
                assids.append(job.ident)

        self.log("------------------------------------")
        self.log("Sample: %s Iteration %s of numcycles %s" % (
            self.params['sample'],
            self.params['iteration'],
            self.globals['numcycles']))
        self.log("------------------------------------")

        finisher_params = {}
        finisher_params['reference'] = os.path.join(self.params['working_dir'], 'I%03d' % self.params['iteration'] + '_contigs.fasta')
        finisher_params['working_dir'] = self.params['working_dir']
        finisher_params['finished_dir'] = self.params['finished_dir']
        finisher_params['sample'] = self.params['sample']
        finisher_params['iteration'] = self.params['iteration']
        finisher_params['targets'] = self.params['targets']
        # These are passed through to every submission.
        finisher_params['readcounts'] = self.params['readcounts']
        finisher_params['mapping_dict'] = self.params['mapping_dict']

        if 'PE1' and 'PE2' in self.params:
            finisher_params['PE1'] = self.params['PE1']
            finisher_params['PE2'] = self.params['PE2']
            idx_PE1.close()
            idx_PE2.close()
        if 'SE' in self.params:
            finisher_params['SE'] = self.params['SE']
            idx_SE.close()

        # if 'PE1' in self.params and 'PE2' in self.params:
        #     idx_PE1.close()
        #     idx_PE2.close()
        # if 'SE' in self.params:
        #     idx_SE.close()

        # Kick off the finisher:
        self.submit(
            Finisher,
            procs=1,
            deps=assids,
            params=self.params)
