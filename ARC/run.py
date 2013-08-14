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
import subprocess
from subprocess import CalledProcessError
from ARC import logger
from ARC import exceptions
from ARC.runners import Mapper


class Run:
    def __init__(self, filename):
        if os.path.exists(filename) is False:
            raise exceptions.FatalError("Error, you must run ARC in a folder containing ARC_config.txt")
        self.filename = filename
        self.config = {}
        self.load()

    def load(self):
        """Read in ARC_config.txt and put it in a datastructure """
        self.config['Samples'] = {}

        inf = open(self.filename, 'r')
        header = True  # hack to skip the header line, maybe there is a more graceful way to do this?
        for line in inf:
            if len(line) > 2 and line[0:2] != '##':
                if line[0] == '#':
                    """ Handle global parameters """
                    line = line.strip().strip("# ")
                    line2 = line.strip().split("=")
                    if len(line2) != 2:
                        raise exceptions.FatalError("Error, parameters not specified correctly, "
                                                    "please use # name=value. Offending entry: \n\t%s" % line)
                    self.config[line2[0].strip()] = line2[1].strip()
                elif header is False:
                    """ Handle Sample information """
                    line2 = line.strip().split()
                    # Check that fields are formatted correctly:
                    if len(line2) != 3:
                        raise exceptions.FatalError("Error, sample description entry is not properly"
                                                    "formatted! Offending entry: %s" % line)
                    Sample_ID = line2[0].strip()
                    FileName = line2[1].strip()
                    FileType = line2[2].strip()
                    if Sample_ID not in self.config['Samples']:
                        self.config['Samples'][Sample_ID] = {}
                    if FileType in self.config['Samples'][Sample_ID]:
                        raise exceptions.FatalError("Error same FileType specified more than once for Sample_ID %s." % Sample_ID)
                    self.config['Samples'][Sample_ID][FileType] = FileName
                else:
                    header = False

        #Check that all files exist:
        if 'reference' in self.config:
            self.config['reference'] = os.path.realpath(self.config['reference'])
            if not os.path.exists(self.config['reference']):
                raise exceptions.FatalError("Error, cannot find reference %s" % self.config['reference'])
        else:
            raise exceptions.FatalError('Error, reference not included in ARC_config.txt')
        if len(self.config['Samples']) > 0:
            for Sample in self.config['Samples']:
                if not (('PE1' in self.config['Samples'][Sample] and 'PE2' in self.config['Samples'][Sample]) or 'SE' in self.config['Samples'][Sample]):
                    raise exceptions.FatalError("Error you must specify PE files and/or a SE file for each sample.")
                if 'PE1' in self.config['Samples'][Sample]:
                    self.config['Samples'][Sample]['PE1'] = os.path.realpath(self.config['Samples'][Sample]['PE1'])
                    if not os.path.exists(self.config['Samples'][Sample]['PE1']):
                        raise exceptions.FatalError("PE1 file indicated but not found: %s" % self.config['Samples'][Sample]['PE1'])
                if 'PE2' in self.config['Samples'][Sample]:
                    self.config['Samples'][Sample]['PE2'] = os.path.realpath(self.config['Samples'][Sample]['PE2'])
                    if not os.path.exists(self.config['Samples'][Sample]['PE2']):
                        raise exceptions.FatalError("PE2 file indicated but not found: %s" % self.config['Samples'][Sample]['PE2'])
                if 'SE' in self.config['Samples'][Sample]:
                    self.config['Samples'][Sample]['SE'] = os.path.realpath(self.config['Samples'][Sample]['SE'])
                    if not os.path.exists(self.config['Samples'][Sample]['SE']):
                        raise exceptions.FatalError("SE file indicated but not found: %s" % self.config['Samples'][Sample]['SE'])
        else:
            raise exceptions.FatalError("Could not find samples in ARC_config.txt")

        #Check that required parameters exist:
        if 'numcycles' not in self.config:
            logger.info("numcycles not specified in ARC_config.txt, defaulting to 1")
            self.config['numcycles'] = 1
        if 'max_incorportaion' not in self.config:
            self.config['max_incorporation'] = 5
        else:
            self.config['max_incorporation'] = int(self.config['max_incorportaion'])
        if 'format' not in self.config:
            raise exceptions.FatalError(
                "Error, file format not specificed in ARC_config.txt.")

        if self.config['format'] != 'fastq' and self.config['format'] != 'fasta':
            raise exceptions.FatalError(
                "Error, format is neither fastq or fasta")

        if 'mapper' not in self.config:
            raise exceptions.FatalError(
                "Error, mapper not specificed in ARC_config.txt")

        elif self.config['mapper'] != 'blat' and self.config['mapper'] != 'bowtie2':
            raise exceptions.FatalError(
                "Error mapper must be either 'blat' or 'bowtie2'")

        if 'assembler' not in self.config:
            raise exceptions.FatalError(
                "Error, assembler not specificed in ARC_config.txt")

        elif self.config['assembler'] != 'spades' and self.config['assembler'] != 'newbler':
            raise exceptions.FatalError(
                "Error assembler must be either 'spades' or 'newbler'")

        if 'urt' in self.config and self.config['urt'].lower() == 'true':
            self.config['urt'] = True
        else:
            self.config['urt'] = False
        if 'verbose' in self.config and self.config['verbose'].lower() == 'true':
            self.config['verbose'] = True
        else:
            self.config['verbose'] = False
        if 'map_against_reads' in self.config and self.config['map_against_reads'].lower() == 'true':
            self.config['map_against_reads'] = True
        else:
            self.config['map_against_reads'] = False
        # if 'minreads' not in self.config:
        #     self.config['minreads'] = 50
        # else:
        #     self.config['minreads'] = int(self.config['minreads'])
        if 'assemblytimeout' not in self.config:
            self.config['assemblytimeout'] = 10 * 60  # in seconds
            logger.warn("Defaulting to 10 minute timeout for assemblies")
        else:
            self.config['assemblytimeout'] = float(self.config['assemblytimeout']) * 60

        #Check that the mapper exists:
        if self.config['mapper'] == 'blat':
            try:
                subprocess.check_output(['which', 'blat'])
            except CalledProcessError:
                raise exceptions.FatalError("Cannot find BLAT mapper in path, or Linux 'which' command is missing")
        elif self.config['mapper'] == 'bowtie2':
            try:
                subprocess.check_output(['which', 'bowtie2-build'])
                subprocess.check_output(['which', 'bowtie2'])
            except CalledProcessError:
                raise exceptions.FatalError("Cannot find 'bowtie2-build' or bowtie2 in path, or Linux 'which' command is missing")

        if not 'mapping_procs' in self.config:
            self.config['mapping_procs'] = 1
            logger.info("Using the default 1 for number of processors dedicated for multicore mapping")
        else:
            self.config['mapping_procs'] = int(self.config['mapping_procs'])

        if not 'assembly_procs' in self.config:
            self.config['assembly_procs'] = 1
            logger.info("Using the default 1 for number of processors dedicated for multicore assembly")
        else:
            self.config['assembly_procs'] = int(self.config['assembly_procs'])

        #Check that the assembler exists:
        if self.config['assembler'] == 'spades':
            try:
                subprocess.check_output(['which', 'spades.py'])
            except CalledProcessError:
                raise exceptions.FatalError("Spades assembler specified, but cannot find spades.py")
        if self.config['assembler'] == 'newbler':
            try:
                subprocess.check_output(['which', 'runAssembly'])
            except CalledProcessError:
                raise exceptions.FatalError("Newbler assembler specified, but cannot find runAssembly")

        #Set internal defaults:
        self.config['iteration'] = 0
        self.config['numcycles'] = int(self.config['numcycles'])

        logger.debug("Config dictionary contains: %s" % (self.config))

    def setup(self):
        """
            Set up working folder for each sample. Also assign a "safe_target"
            name to each target so that folder creation works. This is a little
            bit tricky because if the user has targets with the _:_ seperator in
            the name it messes up the splitter and SAM_to_dict. This code is
            therefore written with the assumption that the user has put the _:_
            in the name purposely so that multiple entries in the reference fasta
            will be treated as a single target.
        """
        from Bio import SeqIO

        for sample in self.config['Samples']:
            s = self.config['Samples'][sample]
            working_dir = os.path.realpath('./working_' + sample)
            finished_dir = os.path.realpath('./finished_' + sample)
            self.config['Samples'][sample]['working_dir'] = working_dir
            self.config['Samples'][sample]['finished_dir'] = finished_dir
            if os.path.exists(working_dir):
                logger.info("WARNING working directory already exists for sample %s, deleting old results if any." % sample)
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

            format = self.config['format']
            """
                Build a separate index for each read file in the input, put
                them in working_dir
            """
            if 'PE1' in s:
                if not os.path.exists(os.path.join(working_dir, "/PE1.idx")):
                    SeqIO.index_db(
                        os.path.join(working_dir, "PE1.idx"),
                        s['PE1'],
                        format,
                        key_function=lambda x: x.split("/")[0])
            if 'PE2' in s:
                if not os.path.exists(os.path.join(working_dir, "PE2.idx")):
                    SeqIO.index_db(
                        os.path.join(working_dir, "PE2.idx"),
                        s['PE2'],
                        format,
                        key_function=lambda x: x.split("/")[0])
            if 'SE' in s:
                if not os.path.exists(os.path.join(working_dir, "SE.idx")):
                    SeqIO.index_db(
                        os.path.join(working_dir, "SE.idx"),
                        s['SE'],
                        format, key_function=lambda x: x.split("/")[0])

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

    def submit(self, batchqueue):
        for sample in self.config['Samples']:
            s = self.config['Samples'][sample]
            params = {}
            params['reference'] = self.config['reference']
            params['mapper'] = self.config['mapper']
            params['assembler'] = self.config['assembler']
            params['verbose'] = self.config['verbose']
            params['format'] = self.config['format']
            params['numcycles'] = self.config['numcycles']
            params['urt'] = self.config['urt']
            params['mapping_procs'] = self.config['mapping_procs']
            params['assembly_procs'] = self.config['assembly_procs']
            params['map_against_reads'] = self.config['map_against_reads']
            params['max_incorporation'] = self.config['max_incorporation']
            params['assemblytimeout'] = self.config['assemblytimeout']
            params['safe_targets'] = self.config['safe_targets']
            params['working_dir'] = s['working_dir']
            params['finished_dir'] = s['finished_dir']
            params['sample'] = sample
            params['iteration'] = 0

            if 'PE1' in s and 'PE2' in s:
                params['PE1'] = s['PE1']
                params['PE2'] = s['PE2']
            if 'SE' in s:
                params['SE'] = s['SE']

            logger.debug("Params in run submission: %s" % (params))

            job = batchqueue.submit(
                Mapper,
                procs=params['mapping_procs'],
                params=params)


