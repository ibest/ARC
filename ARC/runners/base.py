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

from ARC import logger
from ARC import FatalError
from ARC import TimeoutError
from ARC import RerunnableError
from ARC import SubprocessError
from copy import deepcopy
import os
import time
import subprocess
import signal
import sys
import logging


class Base:
    def __init__(self, params):
        self.params = params

    #
    # Helper methods
    #
    def name(self):
        return self.__class__.__name__

    def message(self):
        return 'Starting %s' % self.name

    def queue(self, job_q):
        self.job_q = job_q

    def submit(self, job):
        self.job_q.put(job)

    @classmethod
    def to_job(obj, params):
        return {'runner': obj.__name__,
                'params': deepcopy(params)}

    #
    # Flow
    #
    def runner(self):
        self.setup()
        self.start()
        self.teardown()
        self.clean()

    def clean(self):
        del self.params
        del self.job_q
        self.params = None
        self.job_q = None

    def start(self):
        pass

    def setup(self):
        pass

    def teardown(self):
        pass

    def execute(self):
        pass

    #
    # Logging
    #
    def log(self, msg):
        if logger.level() == logging.DEBUG:
            name = self.name
        else:
            name = self.__class__.__name__
        logger.info("%-12s| %s" % (name, msg))

    def info(self, msg):
        if self.loglevel == logging.DEBUG:
            name = self.name
        else:
            name = self.__class__.__name__
        logger.info("%-12s| %s" % (name, msg))

    def debug(self, msg):
        if self.loglevel == logging.DEBUG:
            name = self.name
            logger.debug("%-12s| %s" % (name, msg))

    def warn(self, msg):
        if logger.level() == logging.DEBUG:
            name = self.name
        else:
            name = self.__class__.__name__
        logger.warn("%-12s| %s" % (name, msg))

    def error(self, msg):
        if logger.level() == logging.DEBUG:
            name = self.name
        else:
            name = self.__class__.__name__
        logger.error("%-12s| %s" % (name, msg))

    def exception(self, exc):
        logger.exception(exc)
