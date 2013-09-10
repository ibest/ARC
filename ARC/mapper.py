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
from collections import Counter
from copy import deepcopy
from Bio import SeqIO
from ARC import exceptions
from ARC import logger
#from ARC import AssemblyRunner
from ARC.assembler import AssemblyRunner
from ARC.assembly_checker import AssemblyChecker
import traceback
import sys


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
        return {'runner': self, 'message': 'Sample: %s Starting mapper.' % self.params['sample'], 'params': self.params}

    def start(self):
        if not('mapper' in self.params):
            raise exceptions.FatalError("mapper not defined in params")
        if self.params['mapper'] == 'bowtie2':
            logger.info("Sample: %s Running bowtie2." % self.params['sample'])
            self.run_bowtie2()
        if self.params['mapper'] == 'blat':
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
            txt = "Sample: %s Error creating working directory." % (self.params['sample']) + '\n\t' + str(exc)
            raise exceptions.FatalError(txt)

        #Check whether to log to temporary file, or default to os.devnull
        if 'verbose' in self.params:
            out = open(os.path.join(working_dir, "mapping_log.txt"), 'w')
        else:
            out = open(os.devnull, 'w')

        #Build index
        base = os.path.join(idx_dir, 'idx')
        logger.info("Sample: %s Calling bowtie2-build." % self.params['sample'])
        logger.info(" ".join(['bowtie2-build', '-f', self.params['reference'], base]))
        try:
            ret = subprocess.call(['bowtie2-build', '-f', self.params['reference'], base], stdout=out, stderr=out)
        except Exception as exc:
            txt = ("Sample %s: Unhandeled error running bowtie2-build" % self.params['sample']) + '\n\t' + str(exc)
            out.close()  # make sure that out is closed before throwing exception
            raise exceptions.FatalError(txt)

        if ret != 0:
            out.close()
            raise exceptions.FatalError("Sample: %s Error creating bowtie2 index, check log file." % self.params['sample'])

        #Do bowtie2 mapping:
        n_bowtieprocs = int(round(max(float(self.params['nprocs'])/len(self.params['Samples']), 1)))
        #print "Number of bowtie2 procs:", n_bowtieprocs
        #args = ['nice', '-n', '19', 'bowtie2', '-I', '0', '-X', '1500', '--local', '-p', self.params['nprocs'], '-x', base]
        #args = ['nice', '-n', '19', 'bowtie2', '-I', '0', '-X', '1500', '--local', '-p', '1', '-x', base]
        #args = ['bowtie2', '-I', '0', '-X', '1500', '--local', '-p', str(n_bowtieprocs), '-x', base]
        args = ['bowtie2', '-I', '0', '-X', '1500', '--local', '-p', str(n_bowtieprocs), '-x', base]
        if self.params['bowtie2_k'] > 1:
            args += ['-k', self.params['bowtie2_k']]
        if self.params['format'] == 'fasta':
            args += ['-f']
        if 'PE1' in self.params and 'PE2' in self.params:
            args += ['-1', self.params['PE1'], '-2', self.params['PE2']]
        if 'SE' in self.params:
            args += ['-U', self.params['SE']]
        args += ['-S', os.path.join(working_dir, 'mapping.sam')]
        logger.info("Sample: %s Calling bowtie2 mapper" % self.params['sample'])
        logger.info(" ".join(args))

        try:
            ret = subprocess.call(args, stdout=out, stderr=out)
            out.close()
        except Exception as exc:
            txt = ("Sample %s: Unhandeled error running bowtie2 mapping" % self.params['sample']) + '\n\t' + str(exc)
            raise exceptions.FatalError(txt)

        out.close()
        if ret != 0:
            raise exceptions.FatalError("Sample %s: Bowtie2 mapping returned an error, check log file." % self.params['sample'])

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
        if self.params['format'] == 'fastq':
            args.append('-fastq')
        if 'fastmap' in self.params:
            args.append('-fastMap')
        args.append(os.path.join(working_dir, 'mapping.psl'))

        logger.info("Sample: %s Calling blat mapper" % self.params['sample'])
        logger.debug(" ".join(args))
        try:
            ret = subprocess.call(args, stdout=out, stderr=out)
        except Exception as exc:
            txt = ("Sample %s: Unhandeled error running blat mapping, check log file." % self.params['sample']) + '\n\t' + str(exc)
            raise exceptions.FatalError(txt)
        finally:
            out.close()
        if ret != 0:
            raise exceptions.FatalError('Sample: %s Error running blat mapping, check log file. \n\t %s' % (self.params['sample'], " ".join(args)))

        #Extract the PSL to a dict
        self.params['mapping_dict'] = self.PSL_to_dict(os.path.join(working_dir, 'mapping.psl'))

        #Cleanup
        os.remove(os.path.join(working_dir, 'mapping.psl'))
        out.close()



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
        logger.info("Sample: %s, Processed %s lines from SAM in %s seconds." % (self.params['sample'], i, time.time() - startT))
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
        logger.info("Sample: %s, Processed %s lines from PSL in %s seconds." % (self.params['sample'], i, time.time() - startT))
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
        """ Split reads and then kick off assemblies once the reads are split for a target, use safe_targets for names"""
        try:
            self.params['iteration'] += 1
            #checker_params = deepcopy(self.params)
            checker_params = {}
            for k in self.params:
                checker_params[k] = self.params[k]
            del checker_params['mapping_dict']
            checker_params['targets'] = {}
            iteration = self.params['iteration']
            if 'PE1' in self.params and 'PE2' in self.params:
                idx_PE1 = SeqIO.index_db(os.path.join(self.params['working_dir'], "PE1.idx"), key_function=lambda x: x.split("/")[0])
                idx_PE2 = SeqIO.index_db(os.path.join(self.params['working_dir'], "PE2.idx"), key_function=lambda x: x.split("/")[0])
            if 'SE' in self.params:
                idx_SE = SeqIO.index_db(os.path.join(self.params['working_dir'], "SE.idx"), key_function=lambda x: x.split("/")[0])
            if 'readcounts' not in checker_params:
                checker_params['readcounts'] = {}
            for target in self.params['mapping_dict']:
                startT = time.time()
                #logger.info("Running splitreads for Sample: %s target: %s" % (self.params['sample'], target))
                target_dir = os.path.join(self.params['working_dir'], self.params['safe_targets'][target])
                if target not in checker_params['readcounts']:
                    checker_params['readcounts'][target] = Counter()
                if os.path.exists(target_dir):
                    os.system("rm -rf %s" % target_dir)
                os.mkdir(target_dir)

                reads = self.params['mapping_dict'][target]
                # track how many total reads were added for this cycle
                checker_params['readcounts'][target][iteration] = len(reads)
                statsf = open(os.path.join(self.params['finished_dir'], "mapping_stats.tsv"), 'a')
                statsf.write('\t'.join([self.params['sample'], target, str(iteration), str(len(reads))]) + '\n')
                statsf.close()

                SEs = PEs = 0
                if 'PE1' and 'PE2' in self.params:
                    outf_PE1 = open(os.path.join(target_dir, "PE1." + self.params['format']), 'w')
                    outf_PE2 = open(os.path.join(target_dir, "PE2." + self.params['format']), 'w')
                if 'SE' in self.params:
                    outf_SE = open(os.path.join(target_dir, "SE." + self.params['format']), 'w')
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

                #Build assembly job:
                assembly_params = {}
                assembly_params['target'] = target
                assembly_params['target_dir'] = target_dir
                assembly_params['iteration'] = iteration
                assembly_params['last_assembly'] = False
                assembler_keys = ['assembler', 'sample', 'verbose', 'format', 'assemblytimeout', 'map_against_reads', 'urt', 'numcycles', 'cdna', 'rip']
                for k in assembler_keys:
                    assembly_params[k] = self.params[k]
                cur_reads = checker_params['readcounts'][target][iteration]  # note that this is a counter, so no key errors can occur
                previous_reads = checker_params['readcounts'][target][iteration - 1]

                #Turn off URT in situations where this will be the last iteration due to readcounts:

                if cur_reads <= previous_reads and iteration > 2 or iteration >= self.params['numcycles']:
                    logger.info("Sample: %s target: %s iteration: %s Setting last_assembly to True" % (self.params['sample'], target, self.params['iteration']))
                    assembly_params['last_assembly'] = True

                #properly handle the case where no reads ended up mapping for the PE or SE inputs:
                if PEs > 0:
                    assembly_params['assembly_PE1'] = os.path.join(target_dir, "PE1." + self.params['format'])
                    assembly_params['assembly_PE2'] = os.path.join(target_dir, "PE2." + self.params['format'])
                if SEs > 0:
                    assembly_params['assembly_SE'] = os.path.join(target_dir, "SE." + self.params['format'])

                #All reads have been written at this point, add an assembly to the queue:
                logger.info("Sample: %s target: %s iteration: %s Split %s reads in %s seconds" % (self.params['sample'], target, self.params['iteration'], len(reads), time.time() - startT))

                #Only add an assembly job and AssemblyChecker target if is there are >0 reads:
                if PEs + SEs > 0:
                    checker_params['targets'][target_dir] = False
                    self.ref_q.put(AssemblyRunner(assembly_params).to_dict())

            logger.info("------------------------------------")
            logger.info("| Sample: %s Iteration %s of numcycles %s" % (checker_params['sample'], checker_params['iteration'], checker_params['numcycles']))
            logger.info("------------------------------------")
            if 'PE1' in self.params and 'PE2' in self.params:
                idx_PE1.close()
                idx_PE2.close()
                del idx_PE1
                del idx_PE2
            if 'SE' in self.params:
                idx_SE.close()
                del idx_SE

            #Kick off a job which checks if all assemblies are done, and if not adds a copy of itself to the job queue
            if len(checker_params['targets']) > 0:
                checker = AssemblyChecker(checker_params)
                self.ref_q.put(checker.to_dict())
            else:
                logger.info("Sample: %s No reads mapped, no more work to do." % checker_params['sample'])
        except:
            print "".join(traceback.format_exception(*sys.exc_info()))
            raise exceptions.FatalError("".join(traceback.format_exception(*sys.exc_info())))

