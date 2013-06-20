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
