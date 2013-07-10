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
import errno
import os
import time
import signal
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
        if not('assembler' in self.params):
            raise exceptions.FatalException("assembler not defined in params")
        if self.params['map_against_reads'] and self.params['iteration'] == 1:
            #print "ASSEMBLER: %s %s MAP AGAINST READS" % (self.params['sample'], self.params['target'])
            self.RunMapAgainstReads()
        elif self.params['assembler'] == 'newbler':
            #logger.info("Running Newbler for sample: %s target: %s" % (self.params['sample'], self.params['target']))
            self.RunNewbler()
        elif self.params['assembler'] == 'spades':
            #logger.info("Running Spades for sample: %s target: %s" % (self.params['sample'], self.params['target']))
            self.RunSpades()
        else:
            raise exceptions.FatalException("Assembler %s isn't recognized." % self.params['assembler'])

    def RunMapAgainstReads(self):
        """
        A pseudo-assembler for cases where we don't actually assemble reads and instead just write them out as contigs.
        """
        #print "Creating finished file: " + os.path.join(self.params['target_dir'], 'finished')
        outf = open(os.path.join(self.params['target_dir'], 'finished'), 'w')
        outf.write("map_against_reads")
        outf.close()

    def kill_process_children(self, pid):
        """
            Base on code from:
            http://stackoverflow.com/questions/1191374/subprocess-with-timeout
            http://stackoverflow.com/questions/6553423/multiple-subprocesses-with-timeouts
        """
        p = subprocess.Popen('ps --no-headers -o pid --ppid %d' % pid, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        pids = [pid]
        pids.extend([int(q) for q in stdout.split()])
        print "Process had", len(pids), "children."
        try:
            for pid in pids:
                os.kill(pid, signal.SIGKILL)
                print "--->KILLED", pid
        except OSError:
            print "--->OSERROR"
            pass

    def RunNewbler(self):
        #Code for running newbler
        """
        Expects params keys:
            PE1 and PE2 and/or SE
            target_dir
            -urt
        """
        #Check for necessary params:
        if not (('assembly_PE1' in self.params and 'assembly_PE2' in self.params) or 'assembly_SE' in self.params):
            raise exceptions.FatalException('Missing self.params in RunNewbler.')

        #Check for necessary files:
        if 'assembly_PE1' in self.params and 'assembly_PE2' in self.params and not(os.path.exists(self.params['assembly_PE1']) or not(os.path.exists(self.params['assembly_PE2']))):
            raise exceptions.FatalException('Missing PE files in RunNewbler.')

        if 'assembly_SE' in self.params and not(os.path.exists(self.params['assembly_SE'])):
            raise exceptions.FatalException('Missing SE file in RunNewbler.')

        sample = self.params['sample']
        target = self.params['target']

        #Building the args
        #args = ['/bio/local/bin/runAssembly']
        args = ['runAssembly']
        args += ['-nobig', '-force', '-cpu', '1']
        if self.params['urt'] and self.params['iteration'] < self.params['numcycles']:
            #only run with the -urt switch when it isn't the final assembly
            args += ['-urt']
        args += ['-o', os.path.join(self.params['target_dir'], 'assembly')]
        if 'assembly_PE1' in self.params and 'assembly_PE2' in self.params:
            args += [self.params['assembly_PE1'], self.params['assembly_PE2']]
        if 'assembly_SE' in self.params:
            args += [self.params['assembly_SE']]
        if self.params['verbose']:
            out = open(os.path.join(self.params['target_dir'], "assembly.log"), 'w')
        else:
            out = open(os.devnull, 'w')

        logger.info("Calling newbler for sample: %s target %s" % (sample, target))
        logger.info(" ".join(args))
        killed = False
        failed = False
        try:
            #ret = subprocess.call(args, stdout=out, stderr=out)
            start = time.time()
            ret = subprocess.Popen(args, stdout=out, stderr=out)
            print "Assembly called"
            pid = ret.pid  # http://stackoverflow.com/questions/1191374/subprocess-with-timeout
            print "pid is", pid
            i = 0
            while ret.poll() is None:
                print "Assembly wait:", i
                i += 1
                time.sleep(.2)
                if time.time() - start > self.params['assemblytimeout']:
                    logger.warn("Sample: %s target: %s Killing assembly after %s seconds" % (sample, target, time.time() - start))
                    #print "os.waitpid..."
                    #vals = os.waitpid(pid, 0)
                    #print "os.waitpid returned", vals
                    print "Calling kill"
                    ret.kill()  # Newbler doesn't seem to actually respond to kill all that reliably
                    time.sleep(2)
                    print "Calling kill_process_children"
                    self.kill_process_children(pid)
                    killed = True
                    break
        except Exception as exc:
            txt = ("Sample: %s, Target: %s: Unhandeled error running Newbler assembly" % (self.params['sample'], self.params['target']))
            txt += '\n\t' + str(exc)
            logger.warn(txt)
            failed = True
            pass
        finally:
            out.close()

        #if ret != 0:
            #raise exceptions.RerunnableError("Newbler assembly failed.")
        if not killed and ret.poll() != 0:
            #raise exceptions.RerunnableError("Newbler assembly failed.")
            failed = True

        if failed:
            outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
            outf.write("assembly_failed")
            outf.close()
        if killed:
            outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
            outf.write("assembly_killed")
            outf.close()
        else:
            #Run finished without error
            logger.info("Sample: %s target: %s Assembly finished in %s seconds" % (sample, target, time.time() - start))
            outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
            outf.write("assembly_complete")
            outf.close()

    def RunSpades(self):
        """
        Several arguments can be passed to spades.py: -1 [PE1], -2 [PE2], -s [SE], and -o [target_dir]
        """
        #Check that required params are available
        if not (('assembly_PE1' in self.params and 'assembly_PE2' in self.params) or ('assembly_SE' in self.params)):
            raise exceptions.FatalException('Missing self.params in RunSpades.')

        #Check that the files actually exist
        if 'assembly_PE1' in self.params and 'assembly_PE2' in self.params and not(os.path.exists(self.params['assembly_PE1']) or not(os.path.exists(self.params['assembly_PE2']))):
            raise exceptions.FatalException('Missing PE files in RunSpades.')
        if 'assembly_SE' in self.params and not(os.path.exists(self.params['assembly_SE'])):
            raise exceptions.FatalException('Missing SE file in RunSpades.')

        sample = self.params['sample']
        target = self.params['target']

        #Build args for assembler call
        args = ['spades.py', '-t', '1']
        if self.params['format'] == 'fasta':
            args.append('--only-assembler')  # spades errors on read correction if the input isn't fastq
        if 'assembly_PE1' in self.params and 'assembly_PE2' in self.params:
            args += ['-1', self.params['assembly_PE1'], '-2', self.params['assembly_PE2']]
        if 'assembly_SE' in self.params:
            args += ['-s', self.params['assembly_SE']]
        args += ['-o', os.path.join(self.params['target_dir'], 'assembly')]
        if self.params['verbose']:
            out = open(os.path.join(self.params['target_dir'], "assembly.log"), 'w')
        else:
            out = open(os.devnull, 'w')

        logger.info("Sample: %s target: %s Running spades assembler." % (sample, target))
        logger.info(" ".join(args))
        killed = False
        failed = False
        start = time.time()
        try:
            #ret = subprocess.call(args, stderr=out, stdout=out)
            ret = subprocess.Popen(args, stdout=out, stderr=out)
            while ret.poll() is None:
                time.sleep(.1)
                if time.time() - start > self.params['assemblytimeout']:
                    ret.kill()
                    killed = True
                    logger.warn("Sample: %s target: %s Assembly killed after %s seconds." % (sample, target, time.time() - start))
                    break
        except Exception as exc:
            txt = ("Sample: %s, Target: %s: Unhandeled error running Spades assembly" % (sample, target))
            txt += '\n\t' + str(exc)
            logger.warn(txt)
            failed = True
            pass
        finally:
            out.close()

        if not killed and ret.poll() != 0:
            failed = True
        if failed:
            outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
            outf.write("assembly_failed")
            outf.close()
        if killed:
            outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
            outf.write("assembly_killed")
            outf.close()
        else:
            #Run finished without error
            logger.info("Sample: %s target: %s Assembly finished in %s seconds" % (sample, target, time.time() - start))
            outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
            outf.write("assembly_complete")
            outf.close()
