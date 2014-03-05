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
from Bio.Seq import Seq
from ARC import logger
from ARC import exceptions
from ARC.runners import Base
from collections import Counter
import traceback
import sys


class Finisher(Base):
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

    Finisher knows that an assembly is finished because it finds a "finished"
    file inside of the assemblies temporary folder. The finished file will
    contain a status which can be one of the following:

    assembly_complete : a normal assembly which completed without error
    assembly_failed   : The assembly failed.
                        Don't copy contigs.
                        Copy reads as contigs with name: Sample_:_Target_:_ReadN
                        (note that if no new reads are mapped, this should
                            auto-terminate, so there should be no need for a
                            minreads type of control, however this may generate
                            spurious assemblies which will need to be accounted
                            for in the cleanup module.)
    map_against_reads : No assembly was attempted.
                        Copy reads as contigs with name: Sample_:_Target_:_ReadN
                        (the logic for when to do this will be handled by the
                            Assembler)

    Output to finished contigs depends on:
    1) If params['iteration'] >= params['numcycles'] all output goes to the
        final contigs file.  In the case where no contigs were assembled
        (finished:assembly_failed or finished:map_against_reads) reads go to
        contigs file.  Reads have the naming convention:
        Sample_:_Target_:_ReadN
    2) If params['readcounts'][iteration] <= params['readcount'][iteration - 1]
        meaning no more reads were incorporated reads + contigs go to final
        output, remove from further mapping/assembly.  In the case where no
        contigs were assembled (finished:assembly_failed or
        finished:map_against_reads) reads go to contigs file. Reads have the
        naming convention: Sample_:_Target_:_ReadN

    3) If reads were mapped but no contigs were generated, write the reads out
        as new targets for the next mapping cycle.
    4) If params['readcounts'][iteration] / params['readcount'][iteration - 1]
        > max_incorportaion: This is probably repetitive sequence,
    5) If 'map_against_reads' in params and params['iteration'] == 1:
        On the first iteration, write all reads out as contigs
        (there is no expectation of an assembly having been done in this case)

    """
    def message(self):
        return 'Finisher for Sample: %s' % self.params['sample']

    def start(self):
        sample = self.params['sample']
        logger.info("Sample: %s Starting finisher" % self.params['sample'])
        finished_dir = self.params['finished_dir']
        sample_finished = False
        targets_written = 0
        iteration = self.params['iteration']

        #Set up output for both finished and additional mapping outputs
        fin_outf = open(os.path.join(finished_dir, 'contigs.fasta'), 'a')
        remap_outf = open(os.path.join(self.params['working_dir'], 'I%03d' % self.params['iteration'] + '_contigs.fasta'), 'w')

        #check whether the sample is globally finished
        if self.params['iteration'] >= self.params['numcycles']:
            sample_finished = True

        #loop over the current set of targets_folders
        for target_folder in self.params['targets']:
            #Extract target specific details:
            target_map_against_reads = False
            safe_target = target_folder.split("/")[-1]  # get last element of path name
            target = self.params['safe_targets'][safe_target]
            cur_reads = self.params['readcounts'][target][iteration]  # note that this is a counter, so no key errors can occur
            previous_reads = self.params['readcounts'][target][iteration - 1]

            #Get finished assembly status:
            with open(os.path.join(target_folder, 'finished'), 'r') as finishedf:
                l = finishedf.readline().strip().split()[0]

            logger.info("Sample: %s target: %s finishing target.." % (self.params['sample'], target))
            logger.info("Sample: %s target: %s iteration: %s Assembly reports status: %s." % (sample, target, self.params['iteration'], l))

            if l in ('assembly_failed', 'map_against_reads'):
                target_map_against_reads = True

            if l == 'assembly_killed':
                #only write out the reads, assembly won't have contigs
                self.write_target(target, target_folder, outf=fin_outf, finished=False, map_against_reads=False, killed=True)
            elif sample_finished:  # everything goes into the final file/folders.
                self.write_target(target, target_folder, outf=fin_outf, finished=True)
            elif target_map_against_reads and cur_reads > previous_reads and iteration < 3:
                #Only map against reads if we have improvement in mapping and we haven't been mapping for multiple iterations
                targets_written += self.write_target(target, target_folder, outf=remap_outf, finished=False, map_against_reads=True)
            else:
                #Check read counts and retire target, or send it back for re-mapping depending on mapped reads
                if iteration > 1 and cur_reads != 0 and previous_reads != 0:
                    if cur_reads / previous_reads > self.params['max_incorporation']:
                        logger.info("Sample %s target %s hit a repetitive region, no more mapping will be done" % (self.params['sample'], target))
                        self.write_target(target, target_folder, outf=fin_outf, finished=True)
                    elif cur_reads <= previous_reads and iteration > 2:
                        #Give the mapper a couple extra iterations in case the first mapping got a lot of reads which didn't assemble
                        logger.info("Sample %s target %s did not incorporate any more reads, no more mapping will be done" % (self.params['sample'], target))
                        self.write_target(target, target_folder, outf=fin_outf, finished=True)
                    else:
                        #nothing fancy is going on, just write the contigs out for remapping
                        targets_written += self.write_target(target, target_folder, outf=remap_outf, finished=False)
                else:
                    #nothing fancy is going on, just write the contigs out for remapping
                    targets_written += self.write_target(target, target_folder, outf=remap_outf, finished=False)

        fin_outf.close()
        remap_outf.close()

        if targets_written > 0:
            # Build a new mapper and put it on the queue
            from ARC.runners import Mapper
            mapper_params = {}
            for k in self.params:
                mapper_params[k] = self.params[k]
            del mapper_params['targets']
            mapper_params['reference'] = os.path.join(self.params['working_dir'], 'I%03d' % self.params['iteration'] + '_contigs.fasta')
            self.submit(Mapper.to_job(mapper_params))
            logger.info("Sample: %s Added new mapper to queue: iteration %s" % (self.params['sample'], self.params['iteration']))

        else:
            logger.info("Sample: %s Mapper not added to queue. Work finished." % self.params['sample'])

    #Functions for repeat masking:
    def num_unmers(self, seq, N):
        #Calculate the number of unique nmers in seq
        nmers = {}
        for i in range(len(seq) - (N - 1)):
            nmers[str(seq[i:i+N]).upper()] = True
        return(len(nmers))

    def mask_seq(self, seq, W=15, N=3):
        #Replace simple repeats with 'n' characters
        #This masks a window if the number of unique Nmers is <
        seq_copy = bytearray(seq)
        i = 0
        for i in range(len(seq) - (W - 1)):
            if self.num_unmers(seq[i:i + W], N) < 7:
                if self.params['mapper'] == 'blat':
                    seq_copy[i:i + W] = seq[i:i + W].lower()
                if self.params['mapper'] == 'bowtie2':
                    seq_copy[i:i + W] = 'n' * W
        return(seq_copy)

    def write_target(self, target, target_folder, outf, finished=False, map_against_reads=False, killed=False):
        # either map_against_reads was passed in, or
        # no contigs were assembled and target isn't finished, or
        # assembler crashed and no contig file was created
        # --> write reads as contigs
        num_contigs = 0  # store how many contigs were written out
        if map_against_reads is False and killed is False:
            if self.params['assembler'] == 'newbler' and not self.params['NewblerMap']:
                contigf = os.path.join(self.params['working_dir'], target_folder, "assembly", "assembly", "454AllContigs.fna")
            if self.params['assembler'] == 'newbler' and self.params['NewblerMap']:
                contigf = os.path.join(self.params['working_dir'], target_folder, "assembly", "mapping", "454AllContigs.fna")
            elif self.params['assembler'] == 'spades':
                contigf = os.path.join(self.params['working_dir'], target_folder, "assembly", "contigs.fasta")
            #add support for a special output if this is the final assembly and newbler -cdna was used:
            if finished and self.params['cdna'] and self.params['assembler'] == 'newbler':
                self.writeCDNAresults(target, target_folder, outf, contigf)
            elif os.path.exists(contigf):
                i = 0
                contig_inf = open(contigf, 'r')
                for contig in SeqIO.parse(contig_inf, 'fasta'):
                    i += 1
                    if finished:
                        contig.name = contig.id = self.params['sample'] + "_:_" + target + "_:_" + "Contig%03d" % i
                    else:
                        contig.name = contig.id = self.params['sample'] + "_:_" + target + "_:_" + "Unfinished%03d" % i
                    contig = contig.upper()
                    #Only mask repeats on intermediate iterations.
                    if self.params['maskrepeats'] and not finished:
                        contig.seq = Seq(str(self.mask_seq(contig.seq.tostring())))
                    #Bowtie2 crashes if a contig is all 'n' so only write it out if it isn't
                    if len(contig.seq) != contig.seq.count('n'):
                        SeqIO.write(contig, outf, "fasta")
                contig_inf.close()
                logger.info("Sample: %s target: %s iteration: %s Finished writing %s contigs " % (self.params['sample'], target, self.params['iteration'], i))
                num_contigs += i
                #if i == 0 and finished is False and self.params['iteration'] < 2:
                #    map_against_reads = True

        if map_against_reads:
            i = 0
            logger.info("Sample %s target %s: Writing reads as contigs." % (self.params['sample'], target))
            if 'PE1' in self.params and 'PE2' in self.params:
                inf_PE1n = os.path.join(target_folder, "PE1." + self.params['format'])
                inf_PE2n = os.path.join(target_folder, "PE2." + self.params['format'])
                if os.path.exists(inf_PE1n) and os.path.exists(inf_PE2n):
                    inf_PE1 = open(inf_PE1n, 'r')
                    inf_PE2 = open(inf_PE2n, 'r')
                    for r in SeqIO.parse(inf_PE1, self.params['format']):
                        i += 1
                        r.name = r.id = self.params['sample'] + "_:_" + target + "_:_" + "Read%04d" % i
                        SeqIO.write(r, outf, "fasta")
                    for r in SeqIO.parse(inf_PE2, self.params['format']):
                        i += 1
                        r.name = r.id = self.params['sample'] + "_:_" + target + "_:_" + "Read%04d" % i
                        SeqIO.write(r, outf, "fasta")
                    inf_PE1.close()
                    inf_PE2.close()

            if 'SE' in self.params:
                inf_SEn = os.path.join(target_folder, "SE." + self.params['format'])
                if os.path.exists(inf_SEn):
                    inf_SE = open(inf_SEn, 'r')
                for r in SeqIO.parse(inf_SE, self.params['format']):
                    i += 1
                    r.name = r.id = self.params['sample'] + "_:_" + target + "_:_" + "Read%04d" % i
                    SeqIO.write(r, outf, "fasta")
                inf_SE.close()
            num_contigs += i

        if finished or killed:
            #Write reads:
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
        # Finally a special case for situations where assembly of a target is killed, but contigs exist from
        # a previous assembly. Note that we only do this when not running in cDNA mode.
        if killed and self.params['iteration'] > 1 and not self.params['cdna']:
            #No contigs will be available, however contigs from the previous iteration will be present in
            # I00N_contigs.fasta, grab these and write them out instead
            logger.info("Sample: %s target: %s iteration: %s Writing contigs from previous iteration."
                        % (self.params['sample'], target, self.params['iteration']))
            contigf = os.path.join(self.params['working_dir'], 'I%03d' % (self.params['iteration'] - 1) + '_contigs.fasta')
            if os.path.exists(contigf):
                for contig in SeqIO.parse(contigf, 'fasta'):
                    if contig.id.split("_:_")[1] == target:
                        contig.name = contig.id = contig.id.replace("Unfinished", "Contig")
                        SeqIO.write(contig, outf, "fasta")
        #Cleanup temporary assembly, and reads:
        if not self.params['keepassemblies']:
            os.system("rm -rf %s" % target_folder)
        if finished or killed:
            return 0
        else:
            return num_contigs

    def writeCDNAresults(self, target, target_folder, outf, contigf):
        """
        This is ONLY called when a cDNA target is finished.

        When doing a cDNA type run, it is very useful to have both the following:
        1) All contigs that belong to a gene (isogroup)
            - It would be particularly good to re-orient these if they are in RC.
        2) Total number of reads assembled in each gene (isogroup)

        Additionally it would be excellent to some day also get the following:
        3) Transcript (isotig) structure
        4) Estimate of isotig specific reads.

        """
        if self.params['assembler'] == 'newbler':
            contigf = os.path.join(self.params['working_dir'], target_folder, "assembly", "assembly", "454AllContigs.fna")
            isotigsf = os.path.join(self.params['working_dir'], target_folder, "assembly", "assembly", "454IsotigsLayout.txt")
            readstatusf = os.path.join(self.params['working_dir'], target_folder, "assembly", "assembly", "454ReadStatus.txt")
        else:
            logger.info("WARNING writeCDNAresults called when assembler was not Newbler")
            return None
        if not (os.path.exists(contigf) and os.path.exists(isotigsf) and os.path.exists(readstatusf)):
            logger.info("CDNA WARNING MISSING FILE!! %s %s" % (target, self.params['sample']))
            logger.info(contigf, os.path.exists(contigf))
            logger.info(isotigsf, os.path.exists(isotigsf))
            logger.info(readstatusf, os.path.exists(readstatusf))
            return None
        #Storage data structures:
        isogroups = {}  # A dict of isogroups which each contain an in-order list of contigs
        readcounts = Counter()  # A dict of all contigs, these contain read counts (from ReadStatus)
        contig_orientation = {}
        contig_to_isogroup = {}
        contig_idx = SeqIO.index(contigf, "fasta")
        # Parse isotigsf:
        igroup = ""
        #print self.params['sample'], target, "Parsing isotigsf: %s" % isotigsf
        for l in open(isotigsf, 'r'):
            #Handle lines with only a '\n'
            if l == '\n':
                pass
            #Handle lines for isogroup:
            elif l[0:9] == '>isogroup':
                igroup = l.strip().split()[0].strip(">")
            #Handle lines containing all contigs:
            elif l.strip().split()[0] == 'Contig':
                l2 = l.strip().split()
                contigs = map(lambda x: "contig" + x, l2[2:-1])
                isogroups[igroup] = contigs
                for contig in contigs:
                    if contig not in contig_orientation:
                        contig_orientation[contig] = '+'
                        contig_to_isogroup[contig] = igroup
                    else:
                        raise exceptions.FatalError('Contig %s in %s more than once' % (contig, contigf))
            #Handle lines containing contig orientation info:
            elif l[0:6] == 'isotig':
                l2 = l[l.find(" ") + 1: l.rfind(" ") - 1]
                l3 = [l2[i:i+6] for i in range(0, len(l2), 6)]
                for i in range(len(l3)):
                    if l3[i][0] == '<':
                        # contig is in reverse orientation
                        contig = isogroups[igroup][i]
                        contig_orientation[contig] = '-'
        #print self.params['sample'], target, "Parsed isotigsf, contigs:", len(isogroups), "contig_to_isogroup", len(contig_to_isogroup), "contig_orientation", len(contig_orientation)
        #Now parse readstatus:
        inf = open(readstatusf, 'r')
        inf.readline()  # discard first line
        for l in inf:
            l2 = l.strip().split('\t')
            #Determine if this read was assembled
            if len(l2) == 8:
                contig = l2[2]
                # Note that there are some built in limits to the number of contigs that can be in an isogroup:
                # http://contig.wordpress.com/2010/08/31/running-newbler-de-novo-transcriptome-assembly-i/
                # These won't appear in the IsotigsLayout.txt, but ARE in the ReadStatus.txt file.
                if contig in contig_to_isogroup:
                    readcounts[contig_to_isogroup[contig]] += 1
                else:
                    readcounts['ExceedsThreshold'] += 1
        #print self.params['sample'], target, "Parse read status"

        #Finally, output all of this information appropriately:
        countsf = open(os.path.join(self.params['finished_dir'], "isogroup_read_counts.tsv"), 'a')
        sample = self.params['sample']
        #First write out readcounts: sample \t target \t isogroup \t readcount
        for isogroup in readcounts:
            countsf.write('\t'.join([sample, target, isogroup, str(readcounts[isogroup])]) + '\n')
        countsf.close()
        #print self.params['sample'], target, "Wrote readcounts"

        #Next write the contigs in proper order and orientation:
        ncontigs = 0
        nisogroups = 0
        for isogroup in isogroups:
            nisogroups += 1
            for contig in isogroups[isogroup]:
                ncontigs += 1
                seqrec = contig_idx[contig]
                #print self.params['sample'], target, seqrec
                if contig_orientation[contig] == '-':
                    seqrec.seq = seqrec.seq.reverse_complement()
                #print self.params['sample'], target, seqrec
                seqrec.name = seqrec.id = sample + "_:_" + target + "_:_" + isogroup + "|" + contig
                #print self.params['sample'], target, seqrec
                SeqIO.write(seqrec, outf, "fasta")
        ## TODO: add support for the ExceedsThreshold contigs
        logger.info("Sample: %s target: %s iteration: %s Finished writing %s contigs, %s isogroups " % (self.params['sample'],
                    target, self.params['iteration'], ncontigs, nisogroups))

