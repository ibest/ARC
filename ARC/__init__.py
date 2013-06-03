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
from ARC import logger
#from ARC import mapper
from ARC import spawn
#from ARC import config


def main():
    try:
        config = read_config()
        logger.setup(loglevel=logging.INFO)

        #Run modules:
        setup()
        run_spawner(config)
        clean()
    except (KeyboardInterrupt, SystemExit):
        logger.error("%s unexpectedly terminated" % (__name__))
        clean()


def setup():
    """Add setup"""


def read_config():
    # config.read()
    return {}


def run_mapper():
    # mapper.run()
    pass


def run_spawner(config):
    spawn.run(config)


def clean():
    """Clean up"""
