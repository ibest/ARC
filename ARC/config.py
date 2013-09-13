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
import re
import subprocess
import time
from Bio import SeqIO
from subprocess import CalledProcessError
from ARC import logger
from ARC import exceptions


class Config:
    OPTIONS = {
        'numcycles': 1,
        'max_incorporation': 5,
        'bowtie2_k': 5,
        'format': None,
        'mapper': None,
        'assembler': None,
        'urt': False,
        'verbose': False,
        'map_against_reads': False,
        'assemblytimeout': 10,
        'cdna': False,
        'rip': False
    }

    UNIVERSALS = [
        'safe_targets'
        'format'
        'verbose'
        'assembler'
        'mapper'
        'urt'
        'numcycles'
        'max_incorporation'
        'assemblytimeout'
        'map_against_reads'
    ]

    FORMATS = ['fastq', 'fasta']

    ASSEMBLERS = {
        'newbler': ['runAssembly', 'newbler', 'addRun'],
        'spades': ['spades.py']
    }

    MAPPERS = {
        'blat': ['blat'],
        'bowtie2': ['bowtie2-build', 'bowtie2']
    }

    def __init__(self, filename):
        if os.path.exists(filename) is False:
            raise exceptions.FatalError(
                "Error, you must run ARC in a folder containing "
                "ARC_config.txt")
        self.filename = filename

        # Initialize config
        self.config = {}

        # Read config file, set the defaults, and check
        self.read()
        self.set_defaults()
        self.check()
        self.convert()

    def set_defaults(self):
        for key, value in self.OPTIONS.iteritems():
            if key not in self.config:
                if value is None:
                    raise exceptions.FatalError(
                        "Error, %s required but not specificed in "
                        "ARC_self.config.txt" % key)
                else:
                    logger.info(
                        "%s not specified in ARC_config.txt, defaulting to "
                        "%s" % (key, value))
                    self.config[key] = value
        # Anything listed below here is not expected to be in the config but
        # needs initialized
        self.config['iteration'] = 0

    def check_bins(self, bins):
        for bin in bins:
            try:
                subprocess.check_output(['which', bin])
            except CalledProcessError:
                raise exceptions.FatalError(
                    "Cannot find %s in path, or the 'which' "
                    "command is missing" % (bin))

    def read(self):
        infile = open(self.filename, 'r')

        # Read in comments and globals.  Treats '##' as comments and '#' as
        # global variables.
        while True:
            line = infile.readline()
            if not line:
                break

            line = line.strip()
            # Blank line
            if line == "":
                continue

            arr = line.split()

            if arr[0] == "#":
                cfg = arr[1].split('=')
                if len(cfg) != 2 or cfg[1] == "":
                    raise exceptions.FatalError(
                        "Error, parameters not specified correctly, please "
                        "use # name=value. Offending entry: \n\t%s" % arr[1]
                    )
                # Go ahead and convert the things that should be ints to ints
                key = cfg[0].strip()
                value = cfg[1].strip()
                if re.match(r"[0-9]+", value):
                    self.config[key] = int(value)
                elif value == 'True':
                    self.config[key] = True
                elif value == 'False':
                    self.config[key] = False
                else:
                    self.config[key] = value

            elif arr[0] == "##":
                pass
            else:
                # We just sucked in the header for the samples
                break

        # Now get the sample information
        self.config['Samples'] = {}
        while True:
            line = infile.readline()
            if not line:
                break

            line = line.strip()
            # Blank line
            if line == "" or line[0] == '#':
                continue

            arr = line.split()
            if len(arr) != 3:
                raise exceptions.FatalError(
                    "Error, sample description entry is not properly "
                    "formatted! Offending entry: %s" % line
                )

            sample_id = arr[0].strip()
            filename = arr[1].strip()
            filetype = arr[2].strip()

            if sample_id not in self.config['Samples']:
                self.config['Samples'][sample_id] = {}

            if filetype in self.config['Samples'][sample_id]:
                raise exceptions.FatalError(
                    "Error same FileType specified more than once "
                    "for sample_id %s." % sample_id
                )
            if not os.path.exists(filename):
                raise exceptions.FatalError(
                    "%s file indicated but not found: %s" % (
                    filetype, filename))
            else:
                self.config['Samples'][sample_id][filetype] = os.path.realpath(
                    filename)

    def check(self):
        # Check that the reference file exists
        if 'reference' in self.config:
            self.config['reference'] = os.path.realpath(
                self.config['reference']
            )
            if not os.path.exists(self.config['reference']):
                raise exceptions.FatalError(
                    "Error, cannot find reference %s" % (
                    self.config['reference']))
        else:
            raise exceptions.FatalError(
                'Error, reference not included in ARC_self.txt')

        # Check to see if the samples are valid
        if len(self.config['Samples']) > 0:
            for sample in self.config['Samples']:
                pe_one = 'PE1' in self.config['Samples'][sample]
                pe_two = 'PE2' in self.config['Samples'][sample]
                pe = pe_one and pe_two
                se = 'SE' in self.config['Samples'][sample]

                if not (pe or se):
                    raise exceptions.FatalError(
                        "Error you must specify PE files and/or a SE file for "
                        "each sample.")
        else:
            raise exceptions.FatalError(
                "Could not find samples in ARC_self.txt")

        if self.config['format'] not in self.FORMATS:
            raise exceptions.FatalError(
                "Error, file format not specificed in ARC_self.txt.")

        if self.config['mapper'] not in self.MAPPERS:
            raise exceptions.FatalError(
                "Error mapper must be either %s" % (
                self.MAPPERS.keys().join(',')))
        else:
            self.check_bins(self.MAPPERS[self.config['mapper']])

        if self.config['assembler'] not in self.ASSEMBLERS:
            raise exceptions.FatalError(
                "Error assembler must be either %s" % (
                self.ASSEMBLERS.keys().join(',')))
        else:
            self.check_bins(self.ASSEMBLERS[self.config['assembler']])

    def convert(self):
        # Convert minutes to seconds for assembly timeouts
        self.config['assemblytimeout'] *= 10

    def setup(self):
        """
            Set up working folder for each sample. Also assign a "safe_target"
            name to each target so that folder creation works. This is a little
            bit tricky because if the user has targets with the _:_ seperator
            in the name it messes up the splitter and SAM_to_dict. This code is
            therefore written with the assumption that the user has put the _:_
            in the name purposely so that multiple entries in the reference
            fasta will be treated as a single target.
        """
        format = self.config['format']
        for sample in self.config['Samples']:
            s = self.config['Samples'][sample]
            working_dir = os.path.realpath('./working_' + sample)
            finished_dir = os.path.realpath('./finished_' + sample)
            self.config['Samples'][sample]['working_dir'] = working_dir
            self.config['Samples'][sample]['finished_dir'] = finished_dir
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

            #Create a stats file for cdna
            if self.config['cdna']:
                countsf = open(os.path.join(
                    finished_dir,
                    "isogroup_read_counts.tsv"), 'a')
                countsf.write('\t'.join([
                    'Sample', 'Target', 'isogroup', 'readcount']) + '\n')
                countsf.close()

            # Build a separate index for each read file in the input, put them
            # in working_dir
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
                        logger.error(
                            "The number of reads in %s and %s do not match, "
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

            # Read through the reference, set up a set of safe names for the
            # targets:
            safe_targets = {}
            i = 0
            for t in SeqIO.parse(self.config['reference'], "fasta"):
                if len(t.name.split("_:_")) == 1:
                    safe_targets[t.name] = "t__%06d" % i
                    safe_targets["t__%06d" % i] = t.name
                    i += 1
                else:
                    target = t.name.split("_:_")[1]
                    safe_targets[target] = "t__%06d" % i
                    safe_targets["t__%06d" % i] = target
                    i += 1
            self.config['safe_targets'] = safe_targets

    def get(self):
        return self.config

    def params(self):
        return {k: self.config[k] for k in self.config if k not in self.UNIVERSALS}

    def universals(self):
        return {k: self.config[k] for k in self.config if k in self.UNIVERSALS}
