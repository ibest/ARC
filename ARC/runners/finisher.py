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
from Bio import SeqIO
from ARC.runners import ProcessBase


class Finisher(ProcessBase):
    """
    Iterate through all targets, pull out the assembled contigs, rename
    them to:

        Sample_:_Target_:_ContigN
        Output to: finished_SAMPLE/contigs.fasta
                                  /PE1.fast(a/q)
                                  /PE2.fast(a/q)
                                  /SE.fast(a/q)
            OR
                   working_SAMPLE/IN_contigs.fasta

    Finisher knows that an assembly is finished because it finds a
    "finished" file inside of the assemblies temporary folder.  The
    finished file will contain a status which can be one of the following:
        assembly_complete : a normal assembly which completed without error
        assembly_failed   : The assembly failed.
                            Don't copy contigs.
                            Copy reads as contigs with name:
                                Sample_:_Target_:_ReadN
                            (note that if no new reads are mapped, this
                                should auto-terminate, so there should be
                                no need for a minreads type of control,
                                however this may generate spurious assemblies
                                which will need to be accounted for in the
                                cleanup module.)
        map_against_reads : No assembly was attempted.
                            Copy reads as contigs with name:
                                Sample_:_Target_:_ReadN
                            (the logic for when to do this will be handled
                            by the Assembler)

    Output to finished contigs depends on:
    1) If params['iteration'] >= params['numcycles'] all output goes to
        the final contigs file. In the case where no contigs were assembled
        (finished:assembly_failed or finished:map_against_reads) reads go to
        contigs file.
        Reads have the naming convention: Sample_:_Target_:_ReadN
    2) if params['readcounts'][iteration] <= params['readcount'][iteration - 1]
        meaning no more reads were incorporated reads + contigs go to final
        output, remove from further mapping/assembly. In the case where no
        contigs were assembled (finished:assembly_failed or
        finished:map_against_reads) reads go to contigs file.
        Reads have the naming convention: Sample_:_Target_:_ReadN

    3) if reads were mapped but no contigs were generated, write the reads
        out as new targets for the next mapping cycle
    4) if params['readcounts'][iteration] / params['readcount'][iteration - 1] > max_incorportaion:
        This is probably repetitive sequence,
    5) if 'map_against_reads' in params and params['iteration'] == 1:
        On the first iteration, write all reads out as contigs
        (there is no expectation of an assembly having been done in this case)

    """
    def execute(self):
        self.log("Sample: %s Starting finisher" % self.params['sample'])
        finished_dir = self.params['finished_dir']
        sample_finished = False
        targets_written = 0
        #Set up output for both finished and additional mapping outputs
        fin_outf = open(os.path.join(finished_dir, 'contigs.fasta'), 'a')
        remap_outf = open(os.path.join(self.params['working_dir'], 'I%03d' % self.params['iteration'] + '_contigs.fasta'), 'w')
        #check whether the sample is globally finished
        if self.params['iteration'] >= self.globals['numcycles']:
            sample_finished = True

        #loop over the current set of targets_folders
        for target_folder in self.params['targets']:
            target_map_against_reads = False
            safe_target = target_folder.split("/")[-1]
            target = self.globals['safe_targets'][safe_target]
            self.log("Sample: %s target: %s finishing target.." % (
                self.params['sample'], target))
            finishedf = open(os.path.join(target_folder, 'finished'), 'r')
            l = finishedf.readline().strip()
            finishedf.close()
            if l == 'assembly_failed' or l == 'map_against_reads':
                target_map_against_reads = True
            iteration = self.params['iteration']
            # note that this is a counter, so no key errors can occur
            cur_reads = self.params['readcounts'][target][iteration]
            previous_reads = self.params['readcounts'][target][iteration - 1]

            if l == 'assembly_killed':
                #only write out the reads, assembly won't have contigs
                self.write_target(
                    target,
                    target_folder,
                    outf=fin_outf,
                    finished=False,
                    map_against_reads=False,
                    killed=True)
            elif sample_finished:  # everything goes into the final file/folders.
                self.write_target(
                    target,
                    target_folder,
                    outf=fin_outf,
                    finished=True)
            elif target_map_against_reads and cur_reads > previous_reads and iteration < 3:
                # Only map against reads if we have improvement in mapping and
                # we haven't been mapping for multiple iterations
                self.write_target(
                    target,
                    target_folder,
                    outf=remap_outf,
                    finished=False,
                    map_against_reads=True)
                targets_written += 1
            else:
                # Check read counts and retire target, or send it back for
                # re-mapping depending on mapped reads
                if iteration > 1 and cur_reads != 0 and previous_reads != 0:
                    if cur_reads / previous_reads > self.globals['max_incorporation']:
                        self.log("Sample %s target %s hit a repeatitive region, no more mapping will be done" % (self.params['sample'], target))
                        self.write_target(target, target_folder, outf=fin_outf, finished=True)
                    elif cur_reads <= previous_reads and iteration > 3:
                        # Give the mapper a couple extra iterations in case
                        # the first mapping got a lot of reads which didn't
                        # assemble
                        self.log("Sample %s target %s did not incorporate any more reads, no more mapping will be done" % (self.params['sample'], target))
                        self.write_target(target, target_folder, outf=fin_outf, finished=True)
                    else:
                        # Nothing fancy is going on, just write the contigs
                        # out for remapping
                        self.write_target(target, target_folder, outf=remap_outf, finished=False)
                        targets_written += 1
                else:
                    # Nothing fancy is going on, just write the contigs out
                    # for remapping
                    self.write_target(target, target_folder, outf=remap_outf, finished=False)
                    targets_written += 1
        if targets_written > 0:
            # Build a new mapper and put it on the queue
            from ARC.runners import Mapper
            mapper_params = {}
            mapper_params['reference'] = os.path.join(self.params['working_dir'], 'I%03d' % self.params['iteration'] + '_contigs.fasta')
            mapper_params['working_dir'] = self.params['working_dir']
            mapper_params['finished_dir'] = self.params['finished_dir']
            mapper_params['sample'] = self.params['sample']
            mapper_params['iteration'] = self.params['iteration']
            mapper_params['readcounts'] = self.params['readcounts']
            if 'PE1' and 'PE2' in self.params:
                mapper_params['PE1'] = self.params['PE1']
                mapper_params['PE2'] = self.params['PE2']
            if 'SE' in self.params:
                mapper_params['SE'] = self.params['SE']

            # mapper_params = self.params.copy()
            # mapper_params['reference'] = os.path.join(self.params['working_dir'], 'I%03d' % self.params['iteration'] + '_contigs.fasta')
            # if 'PE1' and 'PE2' in self.params:
            #     mapper_params['PE1'] = self.params['PE1']
            #     mapper_params['PE2'] = self.params['PE2']
            # if 'SE' in self.params:
            #     mapper_params['SE'] = self.params['SE']

            self.submit(
                Mapper,
                procs=self.globals['mapping_procs'],
                params=mapper_params)

            self.log("Sample: %s Added new mapper to queue: iteration %s" % (self.params['sample'], self.params['iteration']))
        else:
            self.log("Sample: %s Mapper not added to queue. Work finished." % self.params['sample'])

    def write_target(self, target, target_folder, outf, finished=False, map_against_reads=False, killed=False):
        # either map_against_reads was passed in, or
        # no contigs were assembled and target isn't finished, or
        # assembler crashed and no contig file was created
        # --> write reads as contigs
        if map_against_reads is False and killed is False:
            if self.globals['assembler'] == 'newbler':
                contigf = os.path.join(self.params['working_dir'], target_folder, "assembly", "assembly", "454AllContigs.fna")
            elif self.globals['assembler'] == 'spades':
                contigf = os.path.join(self.params['working_dir'], target_folder, "assembly", "contigs.fasta")
            i = 0
            if os.path.exists(contigf):
                contig_inf = open(contigf, 'r')
                for contig in SeqIO.parse(contig_inf, 'fasta'):
                    i += 1
                    contig.name = contig.id = self.params['sample'] + "_:_" + target + "_:_" + "Contig%03d" % i
                    SeqIO.write(contig, outf, "fasta")
                contig_inf.close()
                self.log("Sample: %s target: %s Finished writing %s contigs " % (self.params['sample'], target, i))
            if i == 0 and finished is False:
                map_against_reads = True

        if map_against_reads:
            i = 0
            self.log("Sample: %s target %s: Writing reads as contigs." % (self.params['sample'], target))
            if 'PE1' and 'PE2' in self.params:
                inf_PE1n = os.path.join(target_folder, "PE1." + self.globals['format'])
                inf_PE2n = os.path.join(target_folder, "PE2." + self.globals['format'])
                if os.path.exists(inf_PE1n) and os.path.exists(inf_PE2n):
                    inf_PE1 = open(inf_PE1n, 'r')
                    inf_PE2 = open(inf_PE2n, 'r')
                    for r in SeqIO.parse(inf_PE1, self.globals['format']):
                        i += 1
                        r.name = r.id = self.params['sample'] + "_:_" + target + "_:_" + "Read%04d" % i
                        SeqIO.write(r, outf, "fasta")
                    for r in SeqIO.parse(inf_PE2, self.globals['format']):
                        i += 1
                        r.name = r.id = self.params['sample'] + "_:_" + target + "_:_" + "Read%04d" % i
                        SeqIO.write(r, outf, "fasta")
                    inf_PE1.close()
                    inf_PE2.close()

            if 'SE' in self.params:
                inf_SEn = os.path.join(target_folder, "SE." + self.globals['format'])
                if os.path.exists(inf_SEn):
                    inf_SE = open(inf_SEn, 'r')
                for r in SeqIO.parse(inf_SE, self.globals['format']):
                    i += 1
                    r.name = r.id = self.params['sample'] + "_:_" + target + "_:_" + "Read%04d" % i
                    SeqIO.write(r, outf, "fasta")
                inf_SE.close()

        if finished or killed:
            #Write reads:
            if 'PE1' and 'PE2' in self.params:
                inf_PE1n = os.path.join(target_folder, "PE1." + self.globals['format'])
                inf_PE2n = os.path.join(target_folder, "PE2." + self.globals['format'])
                if os.path.exists(inf_PE1n) and os.path.exists(inf_PE2n):
                    inf_PE1 = open(inf_PE1n, 'r')
                    inf_PE2 = open(inf_PE2n, 'r')

                    outf_PE1 = open(
                        os.path.join(
                            self.params['finished_dir'],
                            "PE1." + self.globals['format']
                        ), 'a')
                    outf_PE2 = open(
                        os.path.join(
                            self.params['finished_dir'],
                            "PE2." + self.globals['format']
                        ), 'a')

                    for r in SeqIO.parse(inf_PE1, self.globals['format']):
                        r.description = self.params['sample'] + "_:_" + target
                        SeqIO.write(r, outf_PE1, self.globals['format'])
                    for r in SeqIO.parse(inf_PE2, self.globals['format']):
                        r.description = self.params['sample'] + "_:_" + target
                        SeqIO.write(r, outf_PE2, self.globals['format'])
                    outf_PE1.close()
                    outf_PE2.close()

            if 'SE' in self.params:
                inf_SEn = os.path.join(target_folder, "SE." + self.globals['format'])
                if os.path.exists(inf_SEn):
                    inf_SE = open(inf_SEn, 'r')
                    outf_SE = open(os.path.join(self.params['finished_dir'], "SE." + self.globals['format']), 'a')
                    for r in SeqIO.parse(inf_SE, self.globals['format']):
                        r.description = self.params['sample'] + "_:_" + target
                        SeqIO.write(r, outf_SE, self.globals['format'])
                    outf_SE.close()
        #Cleanup temporary assembly, and reads:
        os.system("rm -rf %s" % target_folder)
