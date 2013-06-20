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

#import time
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
            logger.info("Running Newbler for sample: %s target: %s" % (self.params['sample'], self.params['target']))
            self.RunNewbler()
        elif self.params['assembler'] == 'spades':
            logger.info("Running Spades for sample: %s target: %s" % (self.params['sample'], self.params['target']))
            self.RunSpades()
        else:
            raise exceptions.FatalException("Assembler %s isn't recognized." % self.params['assembler'])

    def RunNewbler(self):
        #Code for running newbler
        """
        Expects params keys:
            PE1 and PE2 and/or SE
            target_dir
            -urt
        """
        #Check for necessary params:
        if not (('PE1' in self.params and 'PE2' in self.params) or 'SE' in self.params):
            raise exceptions.FatalException('Missing self.params in RunNewbler.')

        #Check for necessary files:
        if 'PE1' in self.params and 'PE2' in self.params and not(os.path.exists(self.params['PE1']) and os.path.exists(self.params['PE2'])):
            raise exceptions.FatalException('Missing PE files in RunNewbler.')

        if 'SE' in self.params and not(os.path.exists(self.params['SE'])):
            raise exceptions.FatalException('Missing SE file in RunNewbler.')

        #Building the args
        args = ['runAssembly']
        args += ['-nobig', '-force', '-cpu', '1']
        print self.params['iteration'],  self.params['numcycles']
        if 'urt' in self.params and self.params['iteration'] < self.params['numcycles']:
            #only run with the -urt switch when it isn't the final assembly
            args += ['-urt']
        args += ['-o', os.path.join(self.params['target_dir'], 'assembly')]
        if 'PE1' in self.params and 'PE2' in self.params:
            args += [self.params['PE1'], self.params['PE2']]
        if 'SE' in self.params:
            args += [self.params['SE']]
            out = open(os.path.join(self.params['target_dir'], "assembly.log"), 'w')
        else:
            out = open(os.devnull, 'w')

        logger.info("Calling newbler for sample: %s target %s" % (self.params['sample'], self.params['target']))
        logger.info(" ".join(args))
        ret = subprocess.call(args, stderr=out, stdout=out)
        out.close()
        if ret != 0:
            raise exceptions.RerunnableError("Newbler assembly failed")
        else:
            #Run finished without error
            outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
            outf.write("1")
            outf.close()

    def RunSpades(self):
        """
        Several arguments can be passed to spades.py: -1 [PE1], -2 [PE2], -s [SE], and -o [target_dir]
        """
        #Check that required params are available
        if not (('PE1' in self.params and 'PE2' in self.params) or ('SE' in self.params)):
            raise exceptions.FatalException('Missing self.params in RunSpades.')

        #Check that the files actually exist
        if 'PE1' in self.params and 'PE2' in self.params and not(os.path.exists(self.params['PE1']) and os.path.exists(self.params['PE2'])):
            raise exceptions.FatalException('Missing PE files in RunSpades.')
        if 'SE' in self.params and not(os.path.exists(self.params['SE'])):
            raise exceptions.FatalException('Missing SE file in RunSpades.')

        #Build args for assembler call
        args = ['spades.py', '-t', '1']
        if self.params['format'] == 'fasta':
            args.append('--only-assembler')  # spades errors on read correction if the input isn't fastq
        if 'PE1' in self.params and 'PE2' in self.params:
            args += ['-1', self.params['PE1'], '-2', self.params['PE2']]
        if 'SE' in self.params:
            args += ['-s', self.params['SE']]
        args += ['-o', os.path.join(self.params['target_dir'], 'assembly')]
        if 'verbose' in self.params:
            out = open(os.path.join(self.params['target_dir'], "assembly.log"), 'w')
        else:
            out = open(os.devnull, 'w')

        logger.info("Calling spades for sample: %s target %s" % (self.params['sample'], self.params['target']))
        logger.info(" ".join(args))
        ret = subprocess.call(args, stderr=out, stdout=out)
        out.close()

        if ret != 0:
            raise exceptions.RerunnableError("Assembly failed")
        else:
            #Run finished without error
            outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
            outf.write("1")
            outf.close()
