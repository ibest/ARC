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

import os
import time
from ARC.runners import BaseRunner
from ARC import exceptions


class AssemblyRunner(BaseRunner):
    """
    This class represents assembly jobs and handles running assemblies.
    required params:
        assembler, sample, target, PE1 and PE2 and/or SE, target_dir
    """
    def setup(self):
        pe_run = 'assembly_PE1' and 'assembly_PE2' in self.params
        se_run = 'assembly_SE' in self.params
        if not (pe_run or se_run):
            raise exceptions.FatalException('Missing self.params in assembler.')

        pe_one_path = os.path.exists(self.params['assembly_PE1'])
        pe_two_path = os.path.exists(self.params['assembly_PE2'])
        if pe_run and not (pe_one_path and pe_two_path):
            raise exceptions.FatalException('Missing PE files in assembler.')

        se_path = os.path.exists(self.params['assembly_SE'])
        if se_run and not se_path:
            raise exceptions.FatalException('Missing SE file in assembler.')

        self.target_dir = self.params['target_dir']

    def execute(self):
        if not('assembler' in self.params):
            raise exceptions.FatalException("assembler not defined in params")
        if self.params['map_against_reads'] and self.params['iteration'] == 1:
            self.RunMapAgainstReads()
        elif self.params['assembler'] == 'newbler':
            self.RunNewbler()
        elif self.params['assembler'] == 'spades':
            self.RunSpades()
        else:
            raise exceptions.FatalException(
                "Assembler %s isn't recognized." % self.params['assembler'])

    def RunMapAgainstReads(self):
        """
        A pseudo-assembler for cases where we don't actually assemble reads
        and instead just write them out as contigs.
        """
        outf = open(os.path.join(self.target_dir), 'finished'), 'w')
        outf.write("map_against_reads")
        outf.close()

    def RunNewbler(self):
        """
        Expects params keys:
            PE1 and PE2 and/or SE
            target_dir
            -urt
        """
        sample = self.params['sample']
        target = self.params['target']

        # Build args for newAssembly:
        self.info("Calling newAssembly for sample: %s target %s" % (sample, target))

        args = ['newAssembly', '-force', os.path.join(self.target_dir), 'assembly')]
        self.shell(
            args,
            description='Newbler newAssembly (sample: %s target %s)' % (sample, target),
            logfile='assembly.log',
            working_dir=self.target_dir),
            verbose=self.params['verbose'])

        if 'assembly_PE1' in self.params and 'assembly_PE2' in self.params:
            self.info("Calling addRun for sample: %s target %s" % (sample, target))

            args = ['addRun', os.path.join(self.target_dir), 'assembly')]
            args += [self.params['assembly_PE1']]
            self.shell(
                args,
                description='Newbler addRun PE1 (sample: %s target %s)' % (sample, target),
                logfile='assembly.log',
                working_dir=self.target_dir),
                verbose=self.params['verbose'])

            self.info("Calling addRun for sample: %s target %s" % (sample, target))

            args = ['addRun', os.path.join(self.target_dir), 'assembly')]
            args += [self.params['assembly_PE2']]
            self.shell(
                args,
                description='Newbler addRun PE2 (sample: %s target %s)' % (sample, target),
                logfile='assembly.log',
                working_dir=self.target_dir),
                verbose=self.params['verbose'])

        if 'assembly_SE' in self.params:
            self.info("Calling addRun for sample: %s target %s" % (sample, target))

            args = ['addRun', os.path.join(self.target_dir), 'assembly')]
            args += [self.params['assembly_SE']]
            self.shell(
                args,
                description='Newbler addRun SE (sample: %s target %s)' % (sample, target),
                logfile='assembly.log',
                working_dir=self.target_dir),
                verbose=self.params['verbose'])

        #Build args for runProject
        args = ['runProject', '-nobig', '-cpu', self.procs]
        if self.params['urt'] and self.params['iteration'] < self.params['numcycles']:
            #only run with the -urt switch when it isn't the final assembly
            args += ['-urt']
        args += [os.path.join(self.target_dir), 'assembly')]

        start = time.time()
        self.shell(
            args,
            description='Newbler Assembly (sample: %s target %s)' % (sample, target),
            logfile='assembly.log',
            working_dir=self.target_dir),
            verbose=self.params['verbose'],
            timeout=self.params['assemblytimeout'])

        self.info("Sample: %s target: %s iteration: %s Assembly finished in %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
        outf = open(os.path.join(self.target_dir), "finished"), 'w')
        outf.write("assembly_complete")
        outf.close()

    def RunSpades(self):
        """
        Several arguments can be passed to spades.py: -1 [PE1], -2 [PE2], -s [SE], and -o [target_dir]
        """
        sample = self.params['sample']
        target = self.params['target']

        #Build args for assembler call
        args = ['spades.py', '-t', '1']

        if self.params['format'] == 'fasta':
            # spades errors on read correction if the input isn't fastq
            args.append('--only-assembler')

        if 'assembly_PE1' in self.params and 'assembly_PE2' in self.params:
            args += [
                '-1', self.params['assembly_PE1'],
                '-2', self.params['assembly_PE2']]

        if 'assembly_SE' in self.params:
            args += ['-s', self.params['assembly_SE']]

        args += ['-o', os.path.join(self.target_dir), 'assembly')]

        self.info("Sample: %s target: %s Running spades assembler." % (sample, target))

        start = time.time()
        self.shell(
            args,
            description='Spades Assembly (sample: %s target %s)' % (sample, target),
            logfile='assembly.log',
            working_dir=self.target_dir),
            verbose=self.params['verbose'],
            timeout=self.params['assemblytimeout'])

        self.info("Sample: %s target: %s iteration: %s Assembly finished in %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
        outf = open(os.path.join(self.target_dir), "finished"), 'w')
        outf.write("assembly_complete")
        outf.close()
