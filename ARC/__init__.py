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
from subprocess import CalledProcessError
from ARC import logger
#from ARC import mapper
from ARC import spawn
from ARC import exceptions
#from ARC import config


def main():
    try:
        logger.setup(loglevel=logging.DEBUG)

        logger.info("Reading config file...")
        config = read_config()

        logger.info("Setting up working directories...")
        setup(config)

        logger.info("Setting up multiprocessing...")
        run_spawner(config)

        logger.info("Running mapping")

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
    """ Set up working folder for each sample """
    for Sample in config['Samples']:
        working_dir = os.path.realpath('./working_' + Sample)
        os.mkdir(working_dir)


def read_config():
    config = {}
    config['Samples'] = {}
    if os.path.exists('ARC_config.txt') is False:
        raise exceptions.FatalError("Error, you must run ARC in a folder containing ARC_config.txt")
    inf = open("ARC_config.txt", 'r')
    header = True  # hack to skip the header line, maybe there is a more graceful way to do this?
    for line in inf:
        if len(line) > 2 and line[0:2] != '##':
            if line[0] == '#':
                """ Handle global parameters """
                line = line.strip().strip("# ")
                line2 = line.split("=")
                if len(line2) != 2:
                    raise exceptions.FatalError("Error, parameters not specified correctly, "
                                                "please use # name=value. Offending entry: \n\t%s" % line)
                config[line2[0]] = line2[1]
            elif header is False:
                """ Handle Sample information """
                line2 = line.strip().split('\t')
                # Check that fields are formatted correctly:
                if len(line2) != 3:
                    raise exceptions.FatalError("Error, sample description entry is not properly"
                                                "formatted! Offending entry: \n\t%s" % line)
                Sample_ID = line2[0]
                FileName = line2[1]
                FileType = line2[2]
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
            if not (('PE1' in config['Samples'][Sample] and 'PE2' in config['Samples'][Sample]) or 'SE' in config['Samples'][Sample]):
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

    #Check that required parameters exist:
    if 'numcycles' not in config:
        logger.info("numcycles not specified in ARC_config.txt, defaulting to 1")
        config['numcycles'] = 1
    if 'verbose' not in config:
        config['verbose'] = False
    if 'format' not in config:
        raise exceptions.FatalError("Error, file format not specificed in ARC_config.txt.")
    else:
        assert config['format'] == 'fastq' or config['format'] == 'fasta'

    return config


def run_spawner(config):
    spawn.run(config)


def clean():
    """Clean up"""
