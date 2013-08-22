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
import logging, multiprocessing


def setup(logfile=None, loglevel=logging.DEBUG):
    # Set up a multiprocessing logger to hopefully log from all N workers in a safe and simultaneous fashion:
    logger = multiprocessing.get_logger()
    logger.setLevel(loglevel)
    log_handler = logging.StreamHandler(sys.stdout)
    log_handler.setFormatter(logging.Formatter('[%(asctime)s %(levelname)s %(process)s] %(message)s'))
    log_handler.setLevel(loglevel)  # Here's where
    logger.addHandler(log_handler)


def info(message):
    logger = multiprocessing.get_logger()
    logger.info("%s" % (message))


def error(message):
    logger = multiprocessing.get_logger()
    logger.error("%s" % (message))


def debug(message):
    logger = multiprocessing.get_logger()
    logger.debug("%s" % (message))


def warn(message):
    logger = multiprocessing.get_logger()
    logger.warn("%s" % (message))
