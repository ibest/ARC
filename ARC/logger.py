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
import sys
import logging


def setup(logfile=None, loglevel=logging.INFO):
    #Setup a global logger:
    # global logger
    # How should we handle this gracefully for cases where each component is run independently?
    # How should get get log-level (command line switch), and what level should we default to?
    logger = logging.getLogger(__name__)
    # if config['logfile']:
    log_handler = logging.StreamHandler(sys.stdout)
    # else:
    #   log_handler = logging.FileHandler(config['logfile'])
    log_handler.setFormatter(logging.Formatter('[%(asctime)s %(levelname)s %(process)s] %(message)s'))
    log_handler.setLevel(loglevel)  # Here's where
    logger.addHandler(log_handler)
    logger.setLevel(loglevel)  # And here


def info(message):
    logger = logging.getLogger(__name__)
    logger.info("%s" % (message))


def error(message):
    logger = logging.getLogger(__name__)
    logger.error("%s" % (message))

def debug(message):
    logger = logging.getLogger(__name__)
    logger.debug("%s" % (message))
