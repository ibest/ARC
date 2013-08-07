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

import logging
import os
import sys
import subprocess
from Bio import SeqIO
from subprocess import CalledProcessError
from ARC import logger
#from ARC import mapper
from ARC import spawn
from ARC import exceptions
#from ARC import config


def main():
    try:
        logger.setup(loglevel=logging.INFO)

        logger.info("Reading config file...")
        config = read_config()

        logger.info("Setting up working directories and building indexes...")
        setup(config)

        logger.info("Setting up multiprocessing...")
        run_spawner(config)

        #logger.info("Running mapping")

        clean()
        return 0
    except exceptions.FatalError as e:
        logger.error("A fatal error was encountered. \n\t%s" % str(e))
        return 1
    except (KeyboardInterrupt, SystemExit):
        logger.error("%s unexpectedly terminated" % (__name__))
        clean()
        return 1


def setup(config):
    """
        Set up working folder for each sample. Also assign a "safe_target" name to each target so that folder creation works.
        This is a little bit tricky because if the user has targets with the _:_ seperator in the name it messes up the splitter
        and SAM_to_dict. This code is therefore written with the assumption that the user has put the _:_ in the name purposely
        so that multiple entries in the reference fasta will be treated as a single target.
    """
    #format = config['format']
    for sample in config['Samples']:
        working_dir = os.path.realpath('./working_' + sample)
        finished_dir = os.path.realpath('./finished_' + sample)
        #Check wheter finished_dir already exists, and exit if it does
        if os.path.exists(finished_dir):
            raise exceptions.FatalError("%s already exists, please remove or rename this folder and re-run ARC." % finished_dir)
        config['Samples'][sample]['working_dir'] = working_dir
        config['Samples'][sample]['finished_dir'] = finished_dir
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
        stats_outf = open(os.path.join(working_dir, 'mapping_stats.tsv'), 'w')
        stats_outf.write('\t'.join(["Sample", "Target", "Iteration", "PE", "SE"]) + '\n')
        stats_outf.close()
        ## instead of building indexes, we'll just split them using one iteration
        ## Build a separate index for each read file in the input, put them in working_dir

        # if 'PE1' in s:
        #     if not os.path.exists(os.path.join(working_dir, "/PE1.idx")):
        #         SeqIO.index_db(os.path.join(working_dir, "PE1.idx"), s['PE1'], format, key_function=lambda x: x.split("/")[0])
        # if 'PE2' in s:
        #     if not os.path.exists(os.path.join(working_dir, "PE2.idx")):
        #         SeqIO.index_db(os.path.join(working_dir, "PE2.idx"), s['PE2'], format, key_function=lambda x: x.split("/")[0])
        # if 'SE' in s:
        #     if not os.path.exists(os.path.join(working_dir, "SE.idx")):
        #         SeqIO.index_db(os.path.join(working_dir, "SE.idx"), s['SE'], format, key_function=lambda x: x.split("/")[0])

        #Read through the reference, set up a set of safe names for the targets:
        safe_targets = {}
        i = 0
        for t in SeqIO.parse(config['reference'], "fasta"):
            if len(t.name.split("_:_")) == 1:
                safe_targets[t.name] = "t__%06d" % i
                safe_targets["t__%06d" % i] = t.name
                i += 1
            else:
                target = t.name.split("_:_")[1]
                safe_targets[target] = "t__%06d" % i
                safe_targets["t__%06d" % i] = target
                i += 1
        config['safe_targets'] = safe_targets


def read_config():
    """Read in ARC_config.txt and put it in a datastructure """
    config = {}
    config['Samples'] = {}
    if os.path.exists('ARC_config.txt') is False:
        raise exceptions.FatalError("Error, you must run ARC in a folder containing ARC_config.txt")
    inf = open("ARC_config.txt", 'r')
    header = True  # hack to skip the header line, maybe there is a more graceful way to do this?
    for line in inf:
        if len(line) > 2 and line[0:2] != '##':
            if line[0] == '#':
                # Handle global parameters
                line = line.strip().strip("# ")
                line2 = line.strip().split("=")
                if len(line2) != 2:
                    raise exceptions.FatalError("Error, parameters not specified correctly, "
                                                "please use # name=value. Offending entry: \n\t%s" % line)
                config[line2[0].strip()] = line2[1].strip()
            elif header is False:
                # Handle Sample information
                line2 = line.strip().split()
                # Check that fields are formatted correctly:
                if len(line2) != 3:
                    raise exceptions.FatalError("Error, sample description entry is not properly"
                                                "formatted! Offending entry: %s" % line)
                Sample_ID = line2[0].strip()
                FileName = line2[1].strip()
                FileType = line2[2].strip()
                if Sample_ID not in config['Samples']:
                    config['Samples'][Sample_ID] = {}
                if FileType in config['Samples'][Sample_ID]:
                    raise exceptions.FatalError("Error same FileType specified more than once for Sample_ID %s." % Sample_ID)
                config['Samples'][Sample_ID][FileType] = FileName
            else:
                header = False

    #Check that all files exist:
    if 'reference' in config:
        config['reference'] = os.path.realpath(config['reference'])
        if not os.path.exists(config['reference']):
            raise exceptions.FatalError("Error, cannot find reference %s" % config['reference'])
    else:
        raise exceptions.FatalError('Error, reference not included in ARC_config.txt')
    if len(config['Samples']) > 0:
        for Sample in config['Samples']:
            if not ('PE1' and 'PE2' in config['Samples'][Sample] or 'SE' in config['Samples'][Sample]):
                raise exceptions.FatalError("Error you must specify PE files and/or a SE file for each sample.")
            if 'PE1' in config['Samples'][Sample]:
                config['Samples'][Sample]['PE1'] = os.path.realpath(config['Samples'][Sample]['PE1'])
                if not os.path.exists(config['Samples'][Sample]['PE1']):
                    raise exceptions.FatalError("PE1 file indicated but not found: %s" % config['Samples'][Sample]['PE1'])
            if 'PE2' in config['Samples'][Sample]:
                config['Samples'][Sample]['PE2'] = os.path.realpath(config['Samples'][Sample]['PE2'])
                if not os.path.exists(config['Samples'][Sample]['PE2']):
                    raise exceptions.FatalError("PE2 file indicated but not found: %s" % config['Samples'][Sample]['PE2'])
            if 'SE' in config['Samples'][Sample]:
                config['Samples'][Sample]['SE'] = os.path.realpath(config['Samples'][Sample]['SE'])
                if not os.path.exists(config['Samples'][Sample]['SE']):
                    raise exceptions.FatalError("SE file indicated but not found: %s" % config['Samples'][Sample]['SE'])
    else:
        raise exceptions.FatalError("Could not find samples in ARC_config.txt")

    #Check that required parameters exist:
    if 'numcycles' not in config:
        logger.info("numcycles not specified in ARC_config.txt, defaulting to 1")
        config['numcycles'] = 1
    if 'max_incorportaion' not in config:
        config['max_incorporation'] = 5
    else:
        config['max_incorporation'] = int(config['max_incorportaion'])
    if 'format' not in config:
        raise exceptions.FatalError("Error, file format not specificed in ARC_config.txt.")
    if config['format'] != 'fastq' and config['format'] != 'fasta':
        raise exceptions.FatalError("Error, format is neither fastq or fasta")
    if 'mapper' not in config:
        raise exceptions.FatalError("Error, mapper not specificed in ARC_config.txt")
    elif config['mapper'] != 'blat' and config['mapper'] != 'bowtie2':
        raise exceptions.FatalError("Error mapper must be either 'blat' or 'bowtie2'")
    if 'assembler' not in config:
        raise exceptions.FatalError("Error, assembler not specificed in ARC_config.txt")
    elif config['assembler'] != 'spades' and config['assembler'] != 'newbler':
        raise exceptions.FatalError("Error assembler must be either 'spades' or 'newbler'")
    if 'urt' in config and config['urt'].lower() == 'true':
        config['urt'] = True
    else:
        config['urt'] = False
    if 'verbose' in config and config['verbose'].lower() == 'true':
        config['verbose'] = True
    else:
        config['verbose'] = False
    if 'map_against_reads' in config and config['map_against_reads'].lower() == 'true':
        config['map_against_reads'] = True
    else:
        config['map_against_reads'] = False
    # if 'minreads' not in config:
    #     config['minreads'] = 50
    # else:
    #     config['minreads'] = int(config['minreads'])
    if 'assemblytimeout' not in config:
        config['assemblytimeout'] = 10 * 60  # in seconds
        logger.warn("Defaulting to 10 minute timeout for assemblies")
    else:
        config['assemblytimeout'] = float(config['assemblytimeout']) * 60

    #Check that the mapper exists:
    if config['mapper'] == 'blat':
        try:
            subprocess.check_output(['which', 'blat'])
        except CalledProcessError:
            raise exceptions.FatalError("Cannot find BLAT mapper in path, or Linux 'which' command is missing")
    elif config['mapper'] == 'bowtie2':
        try:
            subprocess.check_output(['which', 'bowtie2-build'])
            subprocess.check_output(['which', 'bowtie2'])
        except CalledProcessError:
            raise exceptions.FatalError("Cannot find 'bowtie2-build' or bowtie2 in path, or Linux 'which' command is missing")

    #Check that the assembler exists:
    if config['assembler'] == 'spades':
        try:
            subprocess.check_output(['which', 'spades.py'])
        except CalledProcessError:
            raise exceptions.FatalError("Spades assembler specified, but cannot find spades.py")
    if config['assembler'] == 'newbler':
        try:
            subprocess.check_output(['which', 'runAssembly'])
        except CalledProcessError:
            raise exceptions.FatalError("Newbler assembler specified, but cannot find runAssembly")

    #Set internal defaults:
    config['iteration'] = 0
    config['numcycles'] = int(config['numcycles'])

    return config


def run_spawner(config):
    spawn.run(config)


def clean():
    """Clean up"""
