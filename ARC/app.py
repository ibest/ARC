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
import time
from Bio import SeqIO
from Bio.Seq import Seq
from ARC import Config
from ARC import logger
from ARC import FatalError
from ARC import Spawn
from ARC.functions import *


class App:
    def start(self, loglevel, configfile='ARC_config.txt'):
        try:
            logger.setup(loglevel=loglevel)

            logger.info("Reading config file...")
            config = Config(configfile)
            values = config.get()

            logger.info(
                "Setting up working directories and building indexes...")
            self.setup(values)

            spawn = Spawn(values)

            logger.info("Running ARC.")
            spawn.submit()
            spawn.run()

            logger.info("Cleaning up.")
            self.clean()

            return 0
        except FatalError as e:
            logger.error("A fatal error was encountered. \n\t%s" % str(e))
            return 1
        except (KeyboardInterrupt, SystemExit):
            self.clean()
            logger.error("%s unexpectedly terminated" % (__name__))
            return 1

    def setup(self, config):
        """
            Set up working folder for each sample. Also assign a "safe_target"
            name to each target so that folder creation works. This is a little
            bit tricky because if the user has targets with the _:_ seperator
            in the name it messes up the splitter and SAM_to_dict. This code is
            therefore written with the assumption that the user has put the _:_
            in the name purposely so that multiple entries in the reference
            fasta will be treated as a single target.
        """
        format = config['format']
        for sample in config['Samples']:
            s = config['Samples'][sample]
            working_dir = os.path.realpath(config['workingdirectory'] + '/working_' + sample)
            #working_dir = os.path.realpath('./working_' + sample)
            finished_dir = os.path.realpath('./finished_' + sample)
            config['Samples'][sample]['working_dir'] = working_dir
            config['Samples'][sample]['finished_dir'] = finished_dir
            if os.path.exists(working_dir):
                logger.info(
                    "WARNING working directory already exists for "
                    "sample %s, deleting old results if any." % (sample))
                os.system('rm -rf %s' % finished_dir)
                os.system('rm -rf %s/t__*' % working_dir)
                os.system('rm -rf %s/*.psl' % working_dir)
                os.system('rm %s/I*_contigs.fasta' % working_dir)
                if os.path.exists('%s/idx' % working_dir):
                    os.system('rm -rf %s/idx' % working_dir)
                os.mkdir(finished_dir)
            else:
                os.mkdir(working_dir)
                os.mkdir(finished_dir)

            # Create stats file:
            statsf = open(os.path.join(finished_dir, "mapping_stats.tsv"), 'w')
            statsf.write('\t'.join(
                ['Sample', 'Target', 'Iteration', 'Reads']) + '\n')
            statsf.close()

            # Create Target Summary Table
            tstf = open(os.path.join(finished_dir, "target_summary_table.tsv"), 'w')
            tstf.write('\t'.join(
                ['Sample', 'Target', 'RefLen', 'Status', 'Iteration', 'Reads', 'Contigs', 'ContigLength']) + '\n')
            tstf.close()

            #Create a stats file for cdna
            if config['cdna']:
                countsf = open(os.path.join(finished_dir, "isogroup_read_counts.tsv"), 'a')
                countsf.write('\t'.join(['Sample', 'Target', 'isogroup', 'readcount']) + '\n')
                countsf.close()

            # Build a separate index for each read file in the input, put them
            # in working_dir
            #Consider parallelizing this?
            start = time.time()
            if 'PE1' in s:
                if not os.path.exists(os.path.join(working_dir, "PE1.idx")):
                    print s['PE1']
                    p1 = SeqIO.index_db(
                        os.path.join(working_dir, "PE1.idx"),
                        s['PE1'],
                        format,
                        key_function=lambda x: x.split("/")[0])
            if 'PE2' in s:
                if not os.path.exists(os.path.join(working_dir, "PE2.idx")):
                    print s['PE2']
                    p2 = SeqIO.index_db(
                        os.path.join(working_dir, "PE2.idx"),
                        s['PE2'],
                        format,
                        key_function=lambda x: x.split("/")[0])
                    if len(p1) != len(p2):
                        logger.error("The number of reads in %s and %s do not match, "
                                     "check the config for errors" % (s['PE1'], s['PE2']))
            if 'SE' in s:
                if not os.path.exists(os.path.join(working_dir, "SE.idx")):
                    print s['SE']
                    SeqIO.index_db(
                        os.path.join(working_dir, "SE.idx"),
                        s['SE'],
                        format,
                        key_function=lambda x: x.split("/")[0])

            logger.info(
                "Sample: %s, indexed reads in %s seconds." % (
                    sample, time.time() - start))

            #Read through the references, mask them if necessary

            #mapper_params['reference'] = os.path.join(self.params['working_dir'], 'I%03d' % self.params['iteration'] + '_contigs.fasta')

        # Read through the reference, set up a set of safe names for the targets.
        # Also create the Target Summary Table which is indexed by original target name (following ARC conventions)
        # Also mask sequences and write them to a new set of output files
        #safe_targets is a two-way lookup, meaning it has both the safe target ID and the contig ID.
        summary_stats = {}
        safe_targets = {}
        new_refsf = {}
        for sample in config['Samples']:
            s = config['Samples'][sample]
            new_refsf[sample] = open(os.path.join(s['working_dir'], 'I000_contigs.fasta'), 'w')

        i = 0
        for t in SeqIO.parse(config['reference'], "fasta"):
            if len(t.name.split("_:_")) == 1:
                target = t.name
            else:
                target = t.name.split("_:_")[1]

            safe_targets[target] = "t__%06d" % i
            safe_targets["t__%06d" % i] = target
            i += 1
            if target not in summary_stats:
                summary_stats[target] = {'targetLength': len(t)}
            else:
                summary_stats[target]['targetLength'] = (summary_stats[target]['targetLength'] + len(t))

            #Write contigs:
            if config['maskrepeats']:
                #t.seq = Seq(str(mask_seq(t.seq.tostring(), config['mapper'])))
                t.seq = Seq(str(mask_seq(str(t.seq), config['mapper'])))
            #Bowtie2 crashes if a contig is all 'n' so only write it out if it isn't
            if len(t) != t.seq.count('n'):
                for outf in new_refsf.values():
                    SeqIO.write(t, outf, "fasta")
            else:
                writeTargetStats(finished_dir=s['finished_dir'],
                                 sample=sample,
                                 target=target,
                                 targetLength=summary_stats[target]['targetLength'],
                                 status='MaskedOut',
                                 iteration=0,
                                 readcount=0,
                                 num_contigs=0, contig_length=0)
                del summary_stats[target]

        config['safe_targets'] = safe_targets
        config['summary_stats'] = summary_stats

    def clean(self):
        pass
