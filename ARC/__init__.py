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
from ARC import logger
#from ARC import mapper
from ARC import spawn
from ARC import exceptions
#from ARC import config


def main():
    try:
        logger.setup(loglevel=logging.INFO)
        config = read_config()

        #Run modules:
        setup(config)
        run_spawner(config)
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
    """ Set up working folders, link reference """
    for Sample in config['Samples']:
        working_dir = os.path.realpath('_'.join(['working_', Sample]))
        os.mkdir(working_dir)
    #link in


def read_config():
    config = {}
    config['Samples'] = {}
    if os.path.exists('ARC_config.txt') is False:
        #raise exceptions.FatalException("Missing ")
        raise exceptions.FatalError("Error, you must run ARC in a folder containing ARC_config.txt")
    inf = open("ARC_config.txt", 'r')
    header = True
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
                if len(line2) != 4:
                    logger.error("Error, sample description entry is not properly formatted! Offending entry: \n\t%s" % line)
                    sys.exit()
                Sample_ID = line2[0]
                FileName = line2[1]
                FileType = line2[2]
                FileFormat = line2[3]
                if Sample_ID not in config['Samples']:
                    config['Samples'][Sample_ID] = {}
                if FileType in config['Samples'][Sample_ID]:
                    logger.error("Error same FileType specified more than once for Sample_ID %s." % Sample_ID)
                    sys.exit()
                config['Samples'][Sample_ID][FileType] = {}
                config['Samples'][Sample_ID][FileType]['FileName'] = FileName
                config['Samples'][Sample_ID][FileType]['FileFormat'] = FileFormat
            else:
                header = False

    #Check that all files exist:
    if 'reference' in config:
        config['reference'] = os.path.realpath(config['reference'])
        if not os.path.exists(config['reference']):
            logger.error("Error, cannot find reference %s" % config['reference'])
            sys.exit()
    else:
        logger.error('Error, reference not included in ARC_config.txt')
        sys.exit()
    if len(config['Samples']) > 0:
        for Sample in config['Samples']:
            if not (('PE1' in config['Samples'][Sample] and 'PE2' in config['Samples'][Sample]) or 'SE' in config['Samples'][Sample]):
                logger.error("Error you must specify PE files and/or a SE file for each sample.")
                sys.exit()
            if 'PE1' in config['Samples'][Sample]:
                config['Samples'][Sample]['PE1'] = os.path.realpath(config['Samples'][Sample]['PE1'])
                if not os.path.exists(config['Samples'][Sample]['PE1']):
                    logger.error("PE1 file indicated but not found: %s" % config['Samples'][Sample]['PE1'])
                    sys.exit()
            if 'PE2' in config['Samples'][Sample]:
                config['Samples'][Sample]['PE2'] = os.path.realpath(config['Samples'][Sample]['PE2'])
                if not os.path.exists(config['Samples'][Sample]['PE2']):
                    logger.error("PE2 file indicated but not found: %s" % config['Samples'][Sample]['PE2'])
                    sys.exit()
            if 'SE' in config['Samples'][Sample]:
                config['Samples'][Sample]['SE'] = os.path.realpath(config['Samples'][Sample]['SE'])
                if not os.path.exists(config['Samples'][Sample]['SE']):
                    logger.error("SE file indicated but not found: %s" % config['Samples'][Sample]['SE'])
                    sys.exit()

    return config


def run_mapper():
    # mapper.run()
    pass


def run_spawner(config):
    spawn.run(config)


def clean():
    """Clean up"""
