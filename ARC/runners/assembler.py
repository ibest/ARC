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

import subprocess
import os
import time
import signal
from ARC import logger
from ARC import exceptions
from ARC.runners import Base
import traceback
import sys


class Assembler(Base):
    """
    This class represents assembly jobs and handles running assemblies.
    required params:
        assembler, sample, target, PE1 and PE2 and/or SE, target_dir
    """

    def message(self):
        return 'Assembler for Sample: %s Target: %s' % (self.params['sample'], self.params['target'])

    def start(self):
        if not('assembler' in self.params):
            raise exceptions.FatalError("assembler not defined in params")
        if self.params['map_against_reads'] and self.params['iteration'] == 1:
            self.RunMapAgainstReads()
        elif self.params['assembler'] == 'newbler':
            self.RunNewbler()
        elif self.params['assembler'] == 'spades':
            self.RunSpades()
        else:
            raise exceptions.FatalError("Assembler %s isn't recognized." % self.params['assembler'])

    def RunMapAgainstReads(self):
        """
        A pseudo-assembler for cases where we don't actually assemble reads and instead just write them out as contigs.
        """
        #print "Creating finished file: " + os.path.join(self.params['target_dir'], 'finished')
        start = time.time()
        outf = open(os.path.join(self.params['target_dir'], 'finished'), 'w')
        outf.write("map_against_reads")
        sample = self.params['sample']
        target = self.params['target']
        logger.info("Sample: %s target: %s iteration: %s Assembly finished in %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
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
        #print "Process had", len(pids), "children."
        try:
            for pid in pids:
                os.kill(pid, signal.SIGKILL)
        except OSError:
            #print "--->OSERROR"
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
            raise exceptions.FatalError('Missing self.params in RunNewbler.')

        #Check for necessary files:
        if 'assembly_PE1' in self.params and 'assembly_PE2' in self.params and not(os.path.exists(self.params['assembly_PE1']) or not(os.path.exists(self.params['assembly_PE2']))):
            raise exceptions.FatalError('Missing PE files in RunNewbler.')

        if 'assembly_SE' in self.params and not(os.path.exists(self.params['assembly_SE'])):
            raise exceptions.FatalError('Missing SE file in RunNewbler.')

        sample = self.params['sample']
        target = self.params['target']
        killed = False
        failed = False

        #determine whether to pipe output to a file or /dev/null
        if self.params['verbose']:
            out = open(os.path.join(self.params['target_dir'], "assembly.log"), 'w')
        else:
            out = open(os.devnull, 'w')

        #Build args for newAssembly:
        if self.params['NewblerMap']:
            args = ['newMapping', '-force']
        else:
            args = ['newAssembly', '-force']
        if self.params['last_assembly'] and self.params['cdna']:
            #only run with cdna switch on the final assembly
            args += ['-cdna']
        args += [os.path.join(self.params['target_dir'], 'assembly')]
        logger.debug("Calling newAssembly for sample: %s target %s" % (sample, target))
        logger.info(" ".join(args))
        ret = subprocess.call(args, stdout=out, stderr=out)

        #PacBio
        if self.params['NewblerMap']:
            args = ['setRef', os.path.join(self.params['target_dir'], 'assembly'), self.params['reference']]
            logger.info(" ".join(args))
            ret = subprocess.call(args, stdout=out, stderr=out)

        #Build args for addRun:
        if 'assembly_PE1' in self.params and 'assembly_PE2' in self.params:
            args = ['addRun', os.path.join(self.params['target_dir'], 'assembly')]
            args += [self.params['assembly_PE1']]
            logger.debug("Calling addRun for sample: %s target %s" % (sample, target))
            logger.debug(" ".join(args))
            ret = subprocess.call(args, stdout=out, stderr=out)

            args = ['addRun', os.path.join(self.params['target_dir'], 'assembly')]
            args += [self.params['assembly_PE2']]
            logger.debug("Calling addRun for sample: %s target %s" % (sample, target))
            logger.debug(" ".join(args))
            ret = subprocess.call(args, stdout=out, stderr=out)
        if 'assembly_SE' in self.params:
            args = ['addRun', os.path.join(self.params['target_dir'], 'assembly')]
            args += [self.params['assembly_SE']]
            logger.debug("Calling addRun for sample: %s target %s" % (sample, target))
            logger.debug(" ".join(args))
            ret = subprocess.call(args, stdout=out, stderr=out)

        #Build args for runProject
        args = ['runProject']
        args += ['-cpu', '1']
        if self.params['last_assembly'] and self.params['cdna']:
            args += ['-noace']
        else:
            args += ['-nobig']
        if self.params['urt'] and not self.params['last_assembly']:
            #only run with the -urt switch when it isn't the final assembly
            args += ['-urt']
        if self.params['rip']:
            args += ['-rip']
        args += [os.path.join(self.params['target_dir'], 'assembly')]
        try:
            start = time.time()
            logger.debug("Calling runProject for sample: %s target %s" % (sample, target))
            logger.debug(" ".join(args))
            ret = subprocess.Popen(args, stdout=out, stderr=out)
            pid = ret.pid
            while ret.poll() is None:
                if time.time() - start > self.params['assemblytimeout']:
                    self.kill_process_children(pid)
                    logger.warn("Sample: %s target: %s iteration: %s Killing assembly after %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
                    killed = True
                    break
                time.sleep(.5)
        except Exception as exc:
            txt = "Sample: %s, Target: %s: Unhandeled error running Newbler assembly" % (self.params['sample'], self.params['target'])
            txt += '\n\t' + str(exc) + "".join(traceback.format_exception)
            logger.warn(txt)
            failed = True
            pass
        finally:
            out.close()

        #Sometimes newbler doesn't seem to exit completely:
        self.kill_process_children(pid)

        #if ret != 0:
            #raise exceptions.RerunnableError("Newbler assembly failed.")

        if not killed and ret.poll() != 0:
            #raise exceptions.RerunnableError("Newbler assembly failed.")
            failed = True

        if failed:
            logger.info("Sample: %s target: %s iteration: %s Assembly failed after %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
            outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
            outf.write("assembly_failed\t" + str(time.time() - start))
            outf.close()
        if killed:
            logger.info("Sample: %s target: %s iteration: %s Assembly killed after %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
            outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
            outf.write("assembly_killed\t" + str(time.time() - start))
            outf.close()
        else:
            #Run finished without error
            logger.info("Sample: %s target: %s iteration: %s Assembly finished in %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
            outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
            outf.write("assembly_complete\t" + str(time.time() - start))
            outf.close()

    def RunSpades(self):
        """
        Several arguments can be passed to spades.py: -1 [PE1], -2 [PE2], -s [SE], and -o [target_dir]
        """
        #Check that required params are available
        if not (('assembly_PE1' in self.params and 'assembly_PE2' in self.params) or ('assembly_SE' in self.params)):
            raise exceptions.FatalError('Missing self.params in RunSpades.')

        #Check that the files actually exist
        if 'assembly_PE1' in self.params and 'assembly_PE2' in self.params and not(os.path.exists(self.params['assembly_PE1']) or not(os.path.exists(self.params['assembly_PE2']))):
            raise exceptions.FatalError('Missing PE files in RunSpades.')
        if 'assembly_SE' in self.params and not(os.path.exists(self.params['assembly_SE'])):
            raise exceptions.FatalError('Missing SE file in RunSpades.')

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

        logger.debug("Sample: %s target: %s Running spades assembler." % (sample, target))
        logger.info(" ".join(args))
        killed = False
        failed = False
        start = time.time()
        try:
            #ret = subprocess.call(args, stderr=out, stdout=out)
            ret = subprocess.Popen(args, stdout=out, stderr=out)
            pid = ret.pid
            while ret.poll() is None:
                if time.time() - start > self.params['assemblytimeout']:
                    ret.kill()
                    killed = True
                    logger.warn("Sample: %s target: %s Assembly killed after %s seconds." % (sample, target, time.time() - start))
                    break
                time.sleep(.5)
        except Exception as exc:
            txt = ("Sample: %s, Target: %s: Unhandeled error running Spades assembly" % (sample, target))
            txt += '\n\t' + str(exc)
            logger.warn(txt)
            failed = True
            pass
        finally:
            out.close()

        #Ensure that assembler exits cleanly:
        self.kill_process_children(pid)

        if not killed and ret.poll() != 0:
            failed = True
        if failed:
            logger.info("Sample: %s target: %s iteration: %s Assembly failed after %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
            outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
            outf.write("assembly_failed")
            outf.close()
        elif killed:
            logger.info("Sample: %s target: %s iteration: %s Assembly killed after %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
            outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
            outf.write("assembly_killed")
            outf.close()
        else:
            #Run finished without error
            logger.info("Sample: %s target: %s iteration: %s Assembly finished in %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
            outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
            outf.write("assembly_complete")
            outf.close()
