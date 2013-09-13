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
#import errno
import os
import time
import signal
from ARC import logger
from ARC import FatalError
from ARC import SubprocessError
from ARC import Runner
from ARC.runners import Base
import sys


class Assembler(Base):
    """
    This class represents assembly jobs and handles running assemblies.
    required params:
        assembler, sample, target, PE1 and PE2 and/or SE, target_dir
    """
    def to_dict(self):
        return {'runner': self,
                'message': 'Assembler for Sample: %s Target: %s' % (self.params['sample'], self.params['target']),
                'params': self.params}

    def setup(self):
        pe_run = 'assembly_PE1' and 'assembly_PE2' in self.params
        se_run = 'assembly_SE' in self.params
        if not (pe_run or se_run):
            raise FatalError('Missing self.params in assembler.')

        if pe_run:
            pe_one_path = os.path.exists(self.params['assembly_PE1'])
            pe_two_path = os.path.exists(self.params['assembly_PE2'])
            if not (pe_one_path and pe_two_path):
                raise FatalError('Missing PE files in assembler.')

        if se_run:
            se_path = os.path.exists(self.params['assembly_SE'])
            if not se_path:
                raise FatalError('Missing SE file in assembler.')

        self.target_dir = self.params['target_dir']

    def execute(self):
        if not('assembler' in self.universals):
            raise FatalError("assembler not defined in params")
        if self.universals['map_against_reads'] and self.params['iteration'] == 1:
            self.RunMapAgainstReads()
        elif self.universals['assembler'] == 'newbler':
            self.RunNewbler()
        elif self.universals['assembler'] == 'spades':
            self.RunSpades()
        else:
            raise FatalError("Assembler %s isn't recognized." % self.universals['assembler'])

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
            raise FatalError('Missing self.params in RunNewbler.')

        #Check for necessary files:
        if 'assembly_PE1' in self.params and 'assembly_PE2' in self.params and not(os.path.exists(self.params['assembly_PE1']) or not(os.path.exists(self.params['assembly_PE2']))):
            raise FatalError('Missing PE files in RunNewbler.')

        if 'assembly_SE' in self.params and not(os.path.exists(self.params['assembly_SE'])):
            raise FatalError('Missing SE file in RunNewbler.')

        sample = self.params['sample']
        target = self.params['target']
        killed = False
        failed = False

        #determine whether to pipe output to a file or /dev/null
        # if self.universals['verbose']:
        #     out = open(os.path.join(self.params['target_dir'], "assembly.log"), 'w')
        # else:
        #     out = open(os.devnull, 'w')

        #Build args for newAssembly:
        args = [
            'newAssembly',
            '-force']

        if self.params['last_assembly'] and self.universals['cdna']:
            #only run with the -urt switch when it isn't the final assembly
            args += ['-cdna']

        args += [os.path.join(self.params['target_dir'], 'assembly')]

        # logger.debug("Calling newAssembly for sample: %s target %s" % (sample, target))
        # logger.info(" ".join(args))
        # ret = subprocess.call(args, stdout=out, stderr=out)

        self.shell(
            args,
            description='Newbler newAssembly (Sample: %s Target: %s)' % (
                sample, target),
            logfile='assembly.log',
            working_dir=self.target_dir,
            verbose=self.universals['verbose'])

        #Build args for addRun:
        if 'assembly_PE1' in self.params and 'assembly_PE2' in self.params:
            args = ['addRun', os.path.join(self.params['target_dir'], 'assembly')]
            args += [self.params['assembly_PE1']]
            # logger.debug("Calling addRun for sample: %s target %s" % (sample, target))
            # logger.debug(" ".join(args))
            # ret = subprocess.call(args, stdout=out, stderr=out)
            self.shell(
                args,
                description='Newbler addRun PE1 (Sample: %s Target: %s)' % (
                    sample, target),
                logfile='assembly.log',
                working_dir=self.target_dir,
                verbose=self.universals['verbose'])

            args = ['addRun', os.path.join(self.params['target_dir'], 'assembly')]
            args += [self.params['assembly_PE2']]
            # logger.debug("Calling addRun for sample: %s target %s" % (sample, target))
            # logger.debug(" ".join(args))
            # ret = subprocess.call(args, stdout=out, stderr=out)
            self.shell(
                args,
                description='Newbler addRun PE2 (Sample: %s Target: %s)' % (
                    sample, target),
                logfile='assembly.log',
                working_dir=self.target_dir,
                verbose=self.universals['verbose'])

        if 'assembly_SE' in self.params:
            args = ['addRun', os.path.join(self.params['target_dir'], 'assembly')]
            args += [self.params['assembly_SE']]
            # logger.debug("Calling addRun for sample: %s target %s" % (sample, target))
            # logger.debug(" ".join(args))
            # ret = subprocess.call(args, stdout=out, stderr=out)
            self.shell(
                args,
                description='Newbler addRun SE (Sample: %s Target: %s)' % (
                    sample, target),
                logfile='assembly.log',
                working_dir=self.target_dir,
                verbose=self.universals['verbose'])

        #Build args for runProject
        args = ['newbler']
        args += ['-cpu', '1']
        if self.params['last_assembly'] and self.universals['cdna']:
            args += ['-noace']
        else:
            args += ['-nobig']
        if self.universals['urt'] and self.params['iteration'] < self.params['numcycles']:
            #only run with the -urt switch when it isn't the final assembly
            args += ['-urt']
        if self.universals['rip']:
            args += ['-rip']
        args += [os.path.join(self.params['target_dir'], 'assembly')]

        start = time.time()
        self.shell(
            args,
            description='Newbler Assembly (Sample: %s Target: %s)' % (
                sample, target),
            logfile='assembly.log',
            working_dir=self.target_dir,
            verbose=self.universals['verbose'],
            timeout=self.universals['assemblytimeout'],
            kill_children=True,
            expected_exitcode=255)

        self.log("Sample: %s target: %s iteration: %s Assembly finished in %s seconds" % (sample, target, self.params['iteration'], time.time() - start))


        # try:
        #     start = time.time()
        #     logger.debug("Calling runProject for sample: %s target %s" % (sample, target))
        #     logger.debug(" ".join(args))
        #     ret = subprocess.Popen(args, stdout=out, stderr=out)
        #     pid = ret.pid
        #     while ret.poll() is None:
        #         if time.time() - start > self.universals['assemblytimeout']:
        #             self.kill_process_children(pid)
        #             logger.warn("Sample: %s target: %s iteration: %s Killing assembly after %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
        #             killed = True
        #             break
        #         time.sleep(.5)
        # except Exception as exc:
        #     txt = "Sample: %s, Target: %s: Unhandeled error running Newbler assembly" % (self.params['sample'], self.params['target'])
        #     txt += '\n\t' + str(exc) + "".join(traceback.format_exception)
        #     logger.warn(txt)
        #     failed = True
        #     pass
        # finally:
        #     out.close()

        # #Sometimes newbler doesn't seem to exit completely:
        # self.kill_process_children(pid)

        #if ret != 0:
            #raise exceptions.RerunnableError("Newbler assembly failed.")

        # if not killed and ret.poll() != 0:
        #     #raise exceptions.RerunnableError("Newbler assembly failed.")
        #     failed = True

        # if failed:
        #     logger.info("Sample: %s target: %s iteration: %s Assembly failed after %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
        #     outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
        #     outf.write("assembly_failed\t" + str(time.time() - start))
        #     outf.close()
        # if killed:
        #     logger.info("Sample: %s target: %s iteration: %s Assembly killed after %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
        #     outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
        #     outf.write("assembly_killed\t" + str(time.time() - start))
        #     outf.close()
        # else:
        #     #Run finished without error
        #     logger.info("Sample: %s target: %s iteration: %s Assembly finished in %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
        #     outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
        #     outf.write("assembly_complete\t" + str(time.time() - start))
        #     outf.close()

    def RunSpades(self):
        """
        Several arguments can be passed to spades.py: -1 [PE1], -2 [PE2], -s [SE], and -o [target_dir]
        """
        #Check that required params are available
        if not (('assembly_PE1' in self.params and 'assembly_PE2' in self.params) or ('assembly_SE' in self.params)):
            raise FatalError('Missing self.params in RunSpades.')

        #Check that the files actually exist
        if 'assembly_PE1' in self.params and 'assembly_PE2' in self.params and not(os.path.exists(self.params['assembly_PE1']) or not(os.path.exists(self.params['assembly_PE2']))):
            raise FatalError('Missing PE files in RunSpades.')
        if 'assembly_SE' in self.params and not(os.path.exists(self.params['assembly_SE'])):
            raise FatalError('Missing SE file in RunSpades.')

        sample = self.params['sample']
        target = self.params['target']

        #Build args for assembler call
        args = ['spades.py', '-t', '1']
        if self.universals['format'] == 'fasta':
            args.append('--only-assembler')  # spades errors on read correction if the input isn't fastq
        if 'assembly_PE1' in self.params and 'assembly_PE2' in self.params:
            args += ['-1', self.params['assembly_PE1'], '-2', self.params['assembly_PE2']]
        if 'assembly_SE' in self.params:
            args += ['-s', self.params['assembly_SE']]
        args += ['-o', os.path.join(self.params['target_dir'], 'assembly')]
        # if self.universals['verbose']:
        #     out = open(os.path.join(self.params['target_dir'], "assembly.log"), 'w')
        # else:
        #     out = open(os.devnull, 'w')

        # logger.debug("Sample: %s target: %s Running spades assembler." % (sample, target))
        # logger.info(" ".join(args))
        # killed = False
        # failed = False
        start = time.time()
        self.shell(
            args,
            description='Spades Assembly (Sample: %s Target: %s)' % (
                sample, target),
            logfile='assembly.log',
            working_dir=self.target_dir,
            verbose=self.universals['verbose'],
            timeout=self.universals['assemblytimeout'],
            callback_on_ok=self.output_on_ok,
            callback_on_error=self.output_on_error)

        self.log("Sample: %s target: %s iteration: %s Assembly finished in %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
        # try:
        #     #ret = subprocess.call(args, stderr=out, stdout=out)
        #     ret = subprocess.Popen(args, stdout=out, stderr=out)
        #     pid = ret.pid
        #     while ret.poll() is None:
        #         if time.time() - start > self.universals['assemblytimeout']:
        #             ret.kill()
        #             killed = True
        #             logger.warn("Sample: %s target: %s Assembly killed after %s seconds." % (sample, target, time.time() - start))
        #             break
        #         time.sleep(.5)
        # except Exception as exc:
        #     txt = ("Sample: %s, Target: %s: Unhandeled error running Spades assembly" % (sample, target))
        #     txt += '\n\t' + str(exc)
        #     logger.warn(txt)
        #     failed = True
        #     pass
        # finally:
        #     out.close()

        #Ensure that assembler exits cleanly:
        # self.kill_process_children(pid)

        # if not killed and ret.poll() != 0:
        #     failed = True
        # if failed:
        #     logger.info("Sample: %s target: %s iteration: %s Assembly failed after %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
        #     outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
        #     outf.write("assembly_failed")
        #     outf.close()
        # if killed:
        #     logger.info("Sample: %s target: %s iteration: %s Assembly killed after %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
        #     outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
        #     outf.write("assembly_killed")
        #     outf.close()
        # else:
        #     #Run finished without error
        #     logger.info("Sample: %s target: %s iteration: %s Assembly finished in %s seconds" % (sample, target, self.params['iteration'], time.time() - start))
        #     outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
        #     outf.write("assembly_complete")
        #     outf.close()

    def run_on_exit_ok(self):
        if self.universals['assembler'] == "newbler":
            # Newbler developers in their infinite wisdom decided to have
            # it return the exitcode 255 for everything... So we need to check
            # to see if the fna file actually exists and then raise an exception
            # if it does not. POS
            fna = os.path.join(self.target_dir, 'assembly', 'assembly', '454AllContigs.fna')
            if not os.path.exists(fna):
                raise SubprocessError("Newbler did not create the fna file.")

        outf = open(os.path.join(self.target_dir, "finished"), 'w')
        outf.write("assembly_complete")
        outf.close()

    def run_on_exit_error(self, retval):
        if retval == Runner.TIMEOUTERROR:
            outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
            outf.write("assembly_killed")
            outf.close()
        else:
            outf = open(os.path.join(self.params['target_dir'], "finished"), 'w')
            outf.write("assembly_failed")
            outf.close()
