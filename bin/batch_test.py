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
import os
import logging
from random import randint
from optparse import OptionParser

lib_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../')
if lib_path not in sys.path:
    sys.path.insert(0, lib_path)

from ARC import logger
from ARC import Batch
from ARC import BatchQueues
from ARC.runners import TestRunner

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option(
        "-d", "--debug",
        action="store_true", dest="debug", default=False,
        help="Turn on debug output")

    (options, args) = parser.parse_args()

    if options.debug:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    logger.setup(loglevel=loglevel)

    bq = BatchQueues()
    b = Batch(bq, procs=4)

    ids = []
    for i in xrange(50):
        procs = randint(1, 4)
        params = {
            'value': randint(1, 10),
            'sleep': randint(1, 10),
            'num': i
        }
        if params['value'] > 9 and i > 3:
            d = -(randint(1, 3))
            job = bq.submit(
                TestRunner,
                procs=procs,
                deps=ids[d:],
                params=params)
        else:
            job = bq.submit(
                TestRunner,
                procs=procs,
                params=params)
        ids.append(job.ident)
    b.run()
