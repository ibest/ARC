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
        config = read_config()
        logger.setup(loglevel=logging.INFO)

        #Run modules:
        setup()
        run_spawner(config)
        clean()
    except exceptions.FatalError:
        logger.error("A fatal error was encountered.  See log for details")
    except (KeyboardInterrupt, SystemExit):
        logger.error("%s unexpectedly terminated" % (__name__))
        clean()


def setup():
    """Add setup"""


def read_config():
    config = {}
    config['Samples'] = {}
    if os.path.exists('ARC_config.txt') is False:
        #raise exceptions.FatalException("Missing ")
        print "Error, you must run ARC in a folder containing ARC_config.txt"
        sys.exit()
    inf = open("ARC_config.txt", 'r')
    for line in inf:
        if len(line) > 2 and line[0:2] != '##':
            if line[0] == '#':
                """ Handle global parameters """
                line = line.strip().strip("# ")
                line2 = line.split("=")
                if len(line2) != 2:
                    print "Error, parameters not specified correctly, please use # name=value. Offending entry: \n\t%s" % line
                    sys.exit()
                config[line2[0]] = line2[1]
            else:
                """ Handle Sample information """
                line2 = line.split('\t')
                # Check that fields are formatted correctly:
                if len(line2) != 4:
                    print "Error, sample description entry is not properly formatted! Offending entry: \n\t%s" % line
                    sys.exit()
                Sample_ID = line2[0]
                FileName = line2[1]
                FileType = line2[2]
                FileFormat = line2[3]
                if Sample_ID not in config['Samples']:
                    config['Samples'][Sample_ID] = {}
                if FileType in config['Samples'][Sample_ID]:
                    print "Error same FileType specified more than once for Sample_ID %s." % Sample_ID
                    sys.exit()
                config['Samples'][Sample_ID][FileType] = {}
                config['Samples'][Sample_ID][FileType]['FileName'] = FileName
                config['Samples'][Sample_ID][FileType]['FileName'] = FileFormat

    return config


def run_mapper():
    # mapper.run()
    pass


def run_spawner(config):
    spawn.run(config)


def clean():
    """Clean up"""
