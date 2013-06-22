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
from copy import deepcopy
from Bio import SeqIO
from ARC import logger
#from ARC import exceptions


class Finisher:
    """
    Iterate through all targets, pull out the assembled contigs, rename them to:
    Sample_:_Target_:_ContigN
    If params['iteration'] >= params['numcycles'] create a ./finished_SampleN/contigs.fasta folder/file and write there
    else: write to working_dir/IN_contigs.fasta

    Output to finished contigs depends on:
    1) if params['iteration'] >= params['numcycles'] all output goes to the final contigs file
        reads also go to final contigs file
    2) if params['readcounts'][iteration] <= params['readcount'][iteration - 1] meaning no more reads were incorporated
        reads + contigs go to final output, remove from further mapping/assembly
    3) if reads were mapped but no contigs were generated, write the reads out as new targets for the next mapping cycle
    4) if params['readcounts'][iteration] / params['readcount'][iteration - 1] > max_incorportaion:
        This is probably repetitive sequence,
    5) if 'map_against_reads' in params and params['iteration'] == 1:
        On the first iteration, write all reads out as contigs
        (there is no expectation of an assembly having been done in this case)

    """

    def __init__(self, params):
        self.params = params

    def queue(self, ref_q):
        self.ref_q = ref_q

    def to_dict(self):
        return {'runner': self,
                'message': 'Finisher for Sample: %s' % self.params['sample'],
                'params': self.params}

    def start(self):
        logger.info("Starting finisher for sample: %s" % self.params['sample'])
        finished_dir = self.params['finished_dir']
        sample_finished = False
        map_against_reads = False
        targets_written = 0
        #Set up output for both finished and additional mapping outputs
        fin_outf = open(os.path.join(finished_dir, 'contigs.fasta'), 'w')
        remap_outf = open(os.path.join(self.params['working_dir'], 'I' + str(self.params['iteration']) + '_contigs.fasta'), 'w')
        #check whether the sample is globally finished
        if self.params['iteration'] >= self.params['numcycles']:
            sample_finished = True
        if self.params['map_against_reads'] and self.params['iteration'] == 1:
            logger.info("Sample %s: map_against_reads is set, writing all reads to contigs" % self.params['sample'])
            map_against_reads = True

        #loop over the current set of targets_folders
        for target_folder in self.params['targets']:
            safe_target = target_folder.split("/")[-1]  # get last element of path name
            target = self.params['safe_targets'][safe_target]
            logger.info("Starting finisher for sample: %s, target %s" % (self.params['sample'], target))
            if sample_finished:  # everything goes into the final file/folders.
                self.write_target(target, target_folder, outf=fin_outf, finished=True)
            elif map_against_reads:
                self.write_target(target, target_folder, outf=remap_outf, finished=False, map_against_reads=True)
                targets_written += 1
            else:
                iteration = self.params['iteration']
                cur_reads = self.params['readcounts'][target][iteration]  # note that this is a counter, so no key errors can occur
                previous_reads = self.params['readcounts'][target][iteration - 1]
                #Check read counts and retire target, or send it back for re-mapping depending on mapped reads
                if iteration > 1 and cur_reads != 0 and previous_reads != 0:
                    if cur_reads / previous_reads > self.params['max_incorporation']:
                        logger.info("Sample %s target %s identified as a repeat, no more mapping will be done" % (self.params['sample'], target))
                        self.write_target(target, target_folder, outf=fin_outf, finished=True)
                    elif cur_reads <= previous_reads:
                        logger.info("Sample %s target %s did not incorporate any more reads, no more mapping will be done" % (self.params['sample'], target))
                        self.write_target(target, target_folder, outf=fin_outf, finished=True)
                    else:
                        #nothing fancy is going on, just write the contigs out for remapping
                        self.write_target(target, target_folder, outf=remap_outf, finished=False)
                        targets_written += 1
                else:
                    #nothing fancy is going on, just write the contigs out for remapping
                    self.write_target(target, target_folder, outf=remap_outf, finished=False)
                    targets_written += 1
        if targets_written > 0:
            # Build a new mapper and put it on the queue
            from ARC.mapper import MapperRunner
            params = deepcopy(self.params)
            params['reference'] = os.path.join(self.params['working_dir'], 'I' + str(self.params['iteration']) + '_contigs.fasta')
            if 'PE1' in self.params and 'PE2' in self.params:
                params['PE1'] = self.params['PE1']
                params['PE2'] = self.params['PE2']
            if 'SE' in self.params:
                params['SE'] = self.params['SE']

            mapper = MapperRunner(params)
            self.ref_q.put(mapper.to_dict())
            logger.info("Added new mapper to queue: Sample %s iteration %s" % (self.params['sample'], self.params['iteration']))
        else:
            logger.info("MapperRunner not added to queue")

    def write_target(self, target, target_folder, outf, finished=False, map_against_reads=False):
        if not map_against_reads:
            if self.params['assembler'] == 'newbler':
                contigf = os.path.join(self.params['working_dir'], target_folder, "assembly", "454AllContigs.fna")
            elif self.params['assembler'] == 'spades':
                contigf = os.path.join(self.params['working_dir'], target_folder, "assembly", "contigs.fasta")
            i = 0
            contig_inf = open(contigf, 'r')
            for contig in SeqIO.parse(contig_inf, 'fasta'):
                i += 1
                contig.name = contig.id = self.params['sample'] + "_:_" + target + "_:_" + "Contig%03d" % i
                SeqIO.write(contig, outf, "fasta")
            contig_inf.close()
            logger.info("Finished with contigs for sample %s target %s" % (self.params['sample'], target))
        # either map_against_reads was passed in, or
        # no contigs were assembled and target isn't finished --> write reads as contigs
        if map_against_reads or (i == 0 and not finished):
            i = 0
            logger.info("Sample %s target %s: Writing reads as contigs for mapping" % (self.params['sample'], target))
            if 'PE1' in self.params and 'PE2' in self.params:
                inf_PE1n = os.path.join(target_folder, "PE1." + self.params['format'])
                inf_PE2n = os.path.join(target_folder, "PE2." + self.params['format'])
                if os.path.exists(inf_PE1n) and os.path.exists(inf_PE2n):
                    inf_PE1 = open(inf_PE1n, 'r')
                    inf_PE2 = open(inf_PE2n, 'r')
                    for r in SeqIO.parse(inf_PE1, self.params['format']):
                        i += 1
                        r.name = r.id = self.params['sample'] + "_:_" + target + "_:_" + "Contig%03d" % i
                        SeqIO.write(r, outf, "fasta")
                    for r in SeqIO.parse(inf_PE2, self.params['format']):
                        i += 1
                        r.name = r.id = self.params['sample'] + "_:_" + target + "_:_" + "Contig%03d" % i
                        SeqIO.write(r, outf, "fasta")
                    inf_PE1.close()
                    inf_PE2.close()

            if 'SE' in self.params:
                inf_SEn = os.path.join(target_folder, "SE." + self.params['format'])
                if os.path.exists(inf_SEn):
                    inf_SE = open(inf_SEn, 'r')
                for r in SeqIO.parse(inf_SE, self.params['format']):
                    i += 1
                    r.name = r.id = self.params['sample'] + "_:_" + target + "_:_" + "Contig%03d" % i
                    SeqIO.write(r, outf, "fasta")
                inf_SE.close()

        if finished:
            #Also write reads:
            if 'PE1' in self.params and 'PE2' in self.params:
                inf_PE1n = os.path.join(target_folder, "PE1." + self.params['format'])
                inf_PE2n = os.path.join(target_folder, "PE2." + self.params['format'])
                if os.path.exists(inf_PE1n) and os.path.exists(inf_PE2n):
                    inf_PE1 = open(inf_PE1n, 'r')
                    inf_PE2 = open(inf_PE2n, 'r')

                    outf_PE1 = open(os.path.join(self.params['finished_dir'], "PE1." + self.params['format']), 'a')
                    outf_PE2 = open(os.path.join(self.params['finished_dir'], "PE2." + self.params['format']), 'a')

                    for r in SeqIO.parse(inf_PE1, self.params['format']):
                        r.description = self.params['sample'] + "_:_" + target
                        SeqIO.write(r, outf_PE1, self.params['format'])
                    for r in SeqIO.parse(inf_PE2, self.params['format']):
                        r.description = self.params['sample'] + "_:_" + target
                        SeqIO.write(r, outf_PE2, self.params['format'])
                    outf_PE1.close()
                    outf_PE2.close()

            if 'SE' in self.params:
                inf_SEn = os.path.join(target_folder, "SE." + self.params['format'])
                if os.path.exists(inf_SEn):
                    inf_SE = open(inf_SEn, 'r')
                    outf_SE = open(os.path.join(self.params['finished_dir'], "SE." + self.params['format']), 'a')
                    for r in SeqIO.parse(inf_SE, self.params['format']):
                        r.description = self.params['sample'] + "_:_" + target
                        SeqIO.write(r, outf_SE, self.params['format'])
                    outf_SE.close()

        os.system("rm -rf %s" % target_folder)

    def old_start(self):
        """
        """
        logger.info("Starting Finisher for sample:%s iteration %s" % (self.params['sample'], self.params['iteration']))

        #Check whether we have reached max iterations:
        if self.params['iteration'] >= self.params['numcycles']:
            finished_dir = self.params['finished_dir']
            outf = open(os.path.join(finished_dir, 'contigs.fasta'), 'w')
            sample_finished = True
        else:
            outfn = os.path.join(self.params['working_dir'], 'I' + str(self.params['iteration']) + '_contigs.fasta')
            outf = open(outfn, 'w')
            sample_finished = False
        for target_folder in self.params['targets']:
            safe_target = target_folder.split("/")[-1]
            target = self.params['safe_targets'][safe_target]
            if self.params['assembler'] == 'newbler':
                contigf = os.path.join(self.params['working_dir'], target_folder, "assembly", "454AllContigs.fna")
            elif self.params['assembler'] == 'spades':
                contigf = os.path.join(self.params['working_dir'], target_folder, "assembly", "contigs.fasta")
            i = 1
            contig_inf = open(contigf, 'r')
            for contig in SeqIO.parse(contig_inf, 'fasta'):
                contig.name = contig.id = self.params['sample'] + "_:_" + target + "_:_" + "Contig%03d" % i
                SeqIO.write(contig, outf, "fasta")
                i += 1
            logger.info("Finished with contigs for sample %s target %s" % (self.params['sample'], target))
            contig_inf.close()
            os.system("rm -rf %s" % target_folder)
        outf.close()
        if sample_finished:
            # do some kind of suprious contig filtering etc
            logger.info("Assembly finished for sample: %s" % self.params['sample'])
            #write mapping statistics:
            rc_outf = open(os.path.join(finished_dir, 'readcounts.tsv'), 'w')
            targets = self.params['readcounts'].keys()
            rc_outf.write("iteration\t"+"\t".join(targets))

            #for i in range(self.params['numcycles']):

            #for t in targets:

            return
        if not sample_finished:
            # Build a new mapper and put it on the queue
            from ARC.mapper import MapperRunner
            params = deepcopy(self.params)
            params['reference'] = outfn
            if 'PE1' in self.params and 'PE2' in self.params:
                params['PE1'] = self.params['PE1']
                params['PE2'] = self.params['PE2']
            if 'SE' in self.params:
                params['SE'] = self.params['SE']

            mapper = MapperRunner(params)
            self.ref_q.put(mapper.to_dict())
            logger.info("Added new mapper to que: Sample %s iteration %s" % (self.params['sample'], self.params['iteration']))

