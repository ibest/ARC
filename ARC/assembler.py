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
import subprocess
import os
from ARC import logger
from ARC import exceptions


class AssemblyRunner:
    """
    This class represents assembly jobs and handles running assemblies.
    required params:
        assembler, sample, target, PE1 and PE2 and/or SE, target_dir
    """
    def __init__(self, params):
        self.params = params

    def queue(self, ref_q):
        self.ref_q = ref_q

    def to_dict(self):
        return {'runner': self,
                'message': 'Assmelber for Sample: %s Target: %s' % (self.params['sample'], self.params['target']),
                'params': self.params}

    def start(self):
        #print "Running the mapper"
        if not('assembler' in self.params):
            raise exceptions.FatalException("assembler not defined in params")
        if self.params['assembler'] == 'newbler':
            self.run_newbler()
        elif self.params['assembler'] == 'spades':
            self.run_spades()
        else:
            raise exceptions.FatalException("Assembler %s isn't recognized." % self.params['assembler'])

    def RunNewbler(self, params):
        #Code for running newbler
        """
        Expects params keys:
            PE1 and PE2 and/or SE
            target_dir
            -urt
        """
        #Check for necessary params:
        if not (('PE1' in params and 'PE2' in params) or 'SE' in params):
            raise exceptions.FatalException('Missing params in RunNewbler.')

        #Check for necessary files:
        if 'PE1' in params and 'PE2' in params and not(os.path.exists(params['PE1']) and os.path.exists(params['PE2'])):
            raise exceptions.FatalException('Missing PE files in RunNewbler.')

        if 'SE' in params and not(os.path.exists(params['SE'])):
            raise exceptions.FatalException('Missing SE file in RunNewbler.')

        #Building the args
        args = ['runAssembly']
        args += ['-nobig', '-force', '-cpu', '1']
        if 'urt' in params:
            args += ['-urt']
        args += ['-o', params['target_dir']]
        if 'PE1' in params and 'PE2' in params:
            args += [params['PE1'], params['PE2']]
        if 'SE' in params:
            args += [params['SE']]
        if 'verbose' in params:
            out = open(os.path.join(params['target_dir'], "assembly.log"), 'w')
        else:
            out = open(os.devnull, 'w')

        ret = subprocess.call(args, stderr=out, stdout=out)
        out.close()
        if ret != 0:
            raise exceptions.RerunnableError("Newbler assembly failed")
        else:
            #Run finished without error
            outf = open(os.path.join(params['target_dir'], "assembly.log"), 'w')
            outf.write("1")
            outf.close()

    def RunSpades(self, params):
        """
        Several arguments can be passed to spades.py: -1 [PE1], -2 [PE2], -s [SE], and -o [target_dir]
        """
        #Check that required params are available
        if not (('PE1' in params and 'PE2' in params) or ('SE' in params)):
            raise exceptions.FatalException('Missing params in RunSpades.')

        #Check that the files actually exist
        if 'PE1' in params and 'PE2' in params and not(os.path.exists(params['PE1']) and os.path.exists(params['PE2'])):
            raise exceptions.FatalException('Missing PE files in RunSpades.')
        if 'SE' in params and not(os.path.exists(params['SE'])):
            raise exceptions.FatalException('Missing SE file in RunSpades.')

        #Build args for assembler call
        args = ['spades.py', '-t', '1']
        if params['format'] == 'fasta':
            args.append('--only-assembler')  # spades errors on read correction if the input isn't fastq
        if 'PE1' in params and 'PE2' in params:
            args += ['-1', params['PE1'], '-2', params['PE2']]
        if 'SE' in params:
            args += ['-s', params['SE']]
        args += ['-o', params['target_dir']]
        if 'verbose' in params:
            out = open(os.path.join(params['target_dir'], "assembly.log"), 'w')
        else:
            out = open(os.devnull, 'w')

        ret = subprocess.call(args, stderr=out, stdout=out)
        out.close()

        if ret != 0:
            raise exceptions.RerunnableError("Assembly failed")

    def queue(self, ref_q):
        self.ref_q = ref_q

# def run():
#     print "I'm running the assembler now"


# def cpu_intensive():
#     a, b = 0, 1
#     for i in range(100000):
#         a, b = b, a + b

