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
import time
from ARC import logger
from ARC import queue
from random import randint
from ARC import exceptions


class TestRunner:

    def __init__(self, params={}):
        self.params = params

    def to_dict(self):
        return {'runner': self, 'message': 'Sample Run', 'params': self.params}

    def start(self):
        logger.info("Running foo with %s" % (self.params['foo']))
        self.cpu_intensive()

        num = randint(0, 10)
        if num > 8:
            # logger.info("Not finished yet, adding another test job.")
            # newjob = TestRunner({'foo': num})
            # queue.add(self.ref_q, newjob.to_dict())
            return
        elif num == 5:
            # Need to create a new one to make sure the joinablequeue isn't
            # passed along.
            rerunjob = TestRunner(self.params)
            queue.add(self.ref_q, rerunjob.to_dict())
            raise exceptions.RerunnableError("Oops, need to rerun this one.")
        elif num == 1:
            try:
                # f = open("/foobar")
                pass
            except IOError as e:
                raise exceptions.FatalError(e.strerror)

    def cpu_intensive(self):
        a, b = 0, 1
        for i in range(100000):
            a, b = b, a + b
        time.sleep(randint(2, 9))

    def queue(self, ref_q):
        self.ref_q = ref_q
