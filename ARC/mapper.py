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
from Bio import SeqIO
from ARC import exceptions
from ARC import logger
#from ARC import AssemblyRunner
from ARC.assembler import AssemblyRunner
from ARC.assembly_checker import AssemblyChecker


class MapperRunner:
    """
    This calss handles mapping jobs, as well as converting map results into a text version of a dict.
    required params:
        PE1, PE2, SE, format, mapper, numcycles, reference, sample, verbose, working_dir
    params added:
        mapping_dict
    """
    def __init__(self, params):
        self.params = params

    def queue(self, ref_q):
        self.ref_q = ref_q

    def to_dict(self):
        return {'runner': self, 'message': 'Starting mapper for sample %s' % self.params['sample'], 'params': self.params}

    def start(self):
        if not('mapper' in self.params):
            raise exceptions.FatalError("mapper not defined in params")
        if self.params['mapper'] == 'bowtie2':
            logger.info("Running bowtie2 for %s" % self.params['sample'])
            self.run_bowtie2()
        if self.params['mapper'] == 'blat':
            logger.info("Running blat for %s" % self.params['sample'])
            self.run_blat()
        #Mapping is done, run splitreads:
        logger.info("Running splitreads for %s" % self.params['sample'])
        self.splitreads()

    def run_bowtie2(self):
        """
        Builds idx and runs bowtie2 -I 0 -X 1500 --local
        Expects params:
            sample, target, reference, working_dir, PE1 and PE2 and/or SE
        """
        #Check for necessary params:
        if not ('sample' in self.params and 'reference' in self.params and 'working_dir' in self.params and (('PE1' in self.params and 'PE2' in self.params) or 'SE' in self.params)):
            raise exceptions.FatalError('Missing params in run_bowtie2.')
        #Check for necessary files:
        if os.path.exists(self.params['reference']) is False:
            raise exceptions.FatalError("Missing reference file for mapping")
        if 'PE1' in self.params and 'PE2' in self.params:
            if not (os.path.exists(self.params['PE1']) and os.path.exists(self.params['PE2'])):
                raise exceptions.FatalError("One or both PE files can not be found for mapping.")
        if 'SE' in self.params:
            if not os.path.exists(self.params['SE']):
                raise exceptions.FatalError("SE file cannot be found.")

        #Make idx directory
        try:
            working_dir = self.params['working_dir']
            idx_dir = os.path.realpath(os.path.join(working_dir, 'idx'))
            os.mkdir(idx_dir)
        except Exception as exc:
            txt = "Error creating working directory for Sample: %s" % (self.params['sample']) + '\n\t' + str(exc)
            raise exceptions.FatalError(txt)

        #Check whether to log to temporary file, or default to os.devnull
        if 'verbose' in self.params:
            out = open(os.path.join(working_dir, "mapping_log.txt"), 'w')
        else:
            out = open(os.devnull, 'w')

        #Build index
        base = os.path.join(idx_dir, 'idx')
        ret = subprocess.call(['bowtie2-build', '-f', self.params['reference'], base], stdout=out, stderr=out)
        if ret != 0:
            raise exceptions.FatalError("Error creating bowtie2 index for Sample: %s" % self.params['sample'])

        #Do bowtie2 mapping:
        args = ['bowtie2', '-I', '0', '-X', '1500', '--local', '-x', base]
        if 'PE1' in self.params and 'PE2' in self.params:
            args += ['-1', self.params['PE1'], '-2', self.params['PE2']]
        if 'SE' in self.params:
            args += ['-U', self.params['SE']]
        args += ['-S', os.path.join(working_dir, 'mapping.sam')]
        ret = subprocess.call(args, stdout=out, stderr=out)
        out.close()
        if ret != 0:
            raise exceptions.FatalError("Error running bowtie2 mapping for Sample: %s" % self.params['sample'])

        #Extract the SAM to a dict
        self.params['mapping_dict'] = self.SAM_to_dict(os.path.join(working_dir, 'mapping.sam'))
        #clean up intermediary files:
        os.remove(os.path.join(working_dir, 'mapping.sam'))
        os.system("rm -rf %s" % idx_dir)

    def run_blat(self):
        #Check for necessary params:
        if not ('sample' in self.params and 'reference' in self.params and 'working_dir' in self.params and (('PE1' in self.params and 'PE2' in self.params) or 'SE' in self.params)):
            raise exceptions.FatalError('Missing self.params in run_bowtie2.')
        #Check for necessary files:
        if os.path.exists(self.params['reference']) is False:
            raise exceptions.FatalError("Missing reference file for mapping")
        if 'PE1' in self.params and 'PE2' in self.params:
            if not (os.path.exists(self.params['PE1']) and os.path.exists(self.params['PE2'])):
                raise exceptions.FatalError("One or both PE files can not be found for mapping.")
        if 'SE' in self.params:
            if not os.path.exists(self.params['SE']):
                raise exceptions.FatalError("SE file cannot be found.")

        #Blat doesn't need an index
        working_dir = self.params['working_dir']

        #Check whether to log to temporary file, or default to os.devnull
        if 'verbose' in self.params:
            out = open(os.path.join(working_dir, "mapping_log.txt"), 'w')
        else:
            out = open(os.devnull, 'w')

        #Build a temporary txt file with all of the reads:
        allreads_outf = open(os.path.join(working_dir, 'reads.txt'), 'w')
        if 'PE1' in self.params and 'PE2' in self.params:
            allreads_outf.write(self.params['PE1'] + '\n')
            allreads_outf.write(self.params['PE2'] + '\n')
        if 'SE' in self.params:
            allreads_outf.write(self.params['SE'] + '\n')
        allreads_outf.close()

        #Do blat mapping
        args = ['blat', self.params['reference'], os.path.join(working_dir, 'reads.txt')]
        if 'format' in self.params and self.params['format'] == 'fastq':
            args.append('-fastq')
        if 'fastmap' in self.params:
            args.append('-fastMap')
        args.append(os.path.join(working_dir, 'mapping.psl'))

        ret = subprocess.call(args, stdout=out, stderr=out)
        out.close()
        if ret != 0:
            raise exceptions.FatalError('Error running blat mapping for sample: %s \n\t %s' % (self.params['sample'], " ".join(args)))

        #Extract the PSL to a dict
        self.params['mapping_dict'] = self.PSL_to_dict(os.path.join(working_dir, 'mapping.psl'))
        #os.remove(os.path.join(working_dir, 'mapping.psl'))

    def SAM_to_dict(self, filename):
        """ Read a SAM file to a mapping dict and return it """
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
        logger.info("Processed %s lines in %s seconds." % (i, time.time() - startT))
        return read_map

    def PSL_to_dict(self, filename):
        """Process a PSL file to the dict format """
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
            readid = l2[9].split("/")[0]  # remove unique part of PE reads
            target = l2[13]
            #handle references built using assembled contigs:
            if len(target.split("_:_")) > 1:
                target = target.split("_:_")[1]
            if target not in read_map:
                read_map[target] = {}
            read_map[target][readid] = 1
        logger.info("Processed %s lines in %s seconds." % (i, time.time() - startT))
        return read_map

    # def write_dict(self, filename, read_map):
    #     """
    #     Write a mapping dictionary to a file.
    #     Experimentally not in use, it is faster just to keep the mapping_dict in memory, and probably
    #     doesn't sacrifice much from a functionality standpoint.
    #     """
    #     #startT = time.time()
    #     outf = open(filename, 'w')
    #     for k in read_map.keys():
    #         outf.write(k + '\t' + ",".join(read_map[k].keys()) + '\n')
    #     outf.close()
    #     #logger.info("Wrote all values to txt in %s seconds" % (time.time() - startT))

    # def read_dict(self):
    #     """ Read a mapping dictionary from a file """
    #     startT = time.time()
    #     try:
    #         inf = open(filename, 'r')
    #     except Exception as inst:
    #         if type(inst) == IOError:
    #             logger.error("Failed to open mapping dictionary %s." % filename)
    #         raise inst
    #     new_map = {}
    #     for l in inf:
    #         l2 = l.split('\t')
    #         l3 = l2[1].strip().split(",")
    #         new_map[l2[0]] = {}
    #         for k in l3:
    #             new_map[l2[0]][k] = 1
    #     logger.info("Read all values to txt in %s seconds" % (time.time() - startT))
    #     return new_map

    def splitreads(self):
        """ Split reads and then kick off assemblies once the reads are split for a target"""
        startT = time.time()
        checker_params = self.params
        checker_params['targets'] = {}
        if 'testing' in self.params:  # added for unit testing
            print "Running splitreads for %s" % self.params['sample']
            ars = []
        for target in self.params['mapping_dict'].keys():
            logger.info("Running splitreads for Sample: %s target: %s" % (self.params['sample'], target))
            target_dir = os.path.realpath(self.params['working_dir'] + "/" + target)
            checker_params['targets'][target_dir] = False
            os.mkdir(target_dir)
            assembly_params = {
                'reference': self.params['reference'],
                'working_dir': self.params['working_dir'],
                'sample': self.params['sample'],
                'assembler': self.params['assembler'],
                'format': self.params['format'],
                'verbose': self.params['verbose'],
                'target': target,
                'target_dir': target_dir
            }
            reads = self.params['mapping_dict'][target]
            if 'PE1' in self.params and 'PE2' in self.params:
                outf_PE1 = open(os.path.realpath(target_dir + "/PE1." + self.params['format']), 'w')
                outf_PE2 = open(os.path.realpath(target_dir + "/PE2." + self.params['format']), 'w')
                idx_PE1 = SeqIO.index_db(os.path.realpath(self.params['working_dir'] + "/PE1.idx"), key_function=lambda x: x.split("/")[0])
                idx_PE2 = SeqIO.index_db(os.path.realpath(self.params['working_dir'] + "/PE2.idx"), key_function=lambda x: x.split("/")[0])
                assembly_params['PE1'] = os.path.realpath(target_dir + "/PE1." + self.params['format'])
                assembly_params['PE2'] = os.path.realpath(target_dir + "/PE2." + self.params['format'])
            if 'SE' in self.params:
                outf_SE = open(os.path.realpath(target_dir + "/SE." + self.params['format']), 'w')
                idx_SE = SeqIO.index_db(os.path.realpath(self.params['working_dir'] + "/SE.idx"), key_function=lambda x: x.split("/")[0])
                assembly_params['SE'] = os.path.realpath(target_dir + "/SE." + self.params['format'])
            for readID in reads:
                if 'PE1' in self.params and readID in idx_PE1:
                    read1 = idx_PE1[readID]
                    read2 = idx_PE2[readID]
                    new_readID = readID.replace(":", "_") + ":0:0:0:0#0/"
                    read1.id = read1.name = new_readID + "1"
                    read2.id = read2.name = new_readID + "2"
                    SeqIO.write(read1, outf_PE1, self.params['format'])
                    SeqIO.write(read2, outf_PE2, self.params['format'])
                elif 'SE' in self.params and readID in idx_SE:
                    SeqIO.write(idx_SE[readID], outf_SE, self.params['format'])
            if 'PE1' in self.params:
                outf_PE1.close()
                outf_PE2.close()
            if 'SE' in self.params:
                outf_SE.close()
            #All reads have been written at this point, add an assembly to the queue:
            ar = AssemblyRunner(assembly_params)
            logger.info("Split %s reads for target %s in %s seconds" % (len(reads), target, time.time() - startT))

            if 'testing' in self.params:
                #This is only here to make testing somewhat possible
                ars.append(ar)
            else:
                #Add job to list:
                self.ref_q.put(ar.to_dict())

        #Kick off a job which checks if all assemblies are done, and if not adds a copy of itself to the job queue
        checker_params['iteration'] += 1
        print "------------------------------------"
        print "Iteration %s of numcycles %s" % (checker_params['iteration'], checker_params['numcycles'])
        print "------------------------------------"

        del checker_params['mapping_dict']
        checker = AssemblyChecker(checker_params)
        if 'testing' in self.params:
            ars.append(checker)
            return ars
        else:
            self.ref_q.put(checker.to_dict())
