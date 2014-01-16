# Copyright 2013, Institute for Bioinformatics and Evolutionary Studies
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
import time
import traceback
import os
import signal
from Queue import Empty
from multiprocessing import Process
from ARC import logger
from ARC import exceptions
import ARC.runners


class ProcessRunner(Process):
    def __init__(self, proc, q, status, stats, ppid):
        super(ProcessRunner, self).__init__()
        self.proc = proc
        self.q = q
        self.status = status
        self.stats = stats
        self.ppid = ppid

    def launch(self):
        # Block until there is an item on the queue
        item = self.q.get()

        # Run the job
        self.running()
        job = getattr(ARC.runners, item['runner'])(item['params'])
        logger.debug("[%s] Processing: %s" % (self.name, job.message()))
        job.queue(self.q)
        job.start()

        # Update stats
        self.update_jobstats(item['runner'])
        
        # Clean up
        job.clean()
        del job
        job = None
        del item
        item = None

        # Notify that the task has been completed
        self.q.task_done()

    def run(self):
        while True:
            try:
                self.waiting()
                self.launch()
                self.update_runstats()
            except exceptions.RerunnableError as e:
                logger.warn("[%s] A job needs to be rerun: %s" % (self.name, e))
                self.update_runstats(1)
            except exceptions.FatalError as e:
                logger.error("[%s] A fatal error occurred: %s" % (self.name, e))
                os.kill(self.ppid, signal.SIGINT)
            except (KeyboardInterrupt, SystemExit):
                logger.debug("Process interrupted")
                sys.exit()
            except Exception as e:
                ex_type, ex, tb = sys.exc_info()
                logger.error("\n".join(traceback.format_exception(ex_type, ex, tb)))
                logger.error("An unhandled exception occurred")
                os.kill(self.ppid, signal.SIGINT)

    def waiting(self):
        self.status[self.proc] = 1

    def is_waiting(self):
        self.status[self.proc] == 1

    def running(self):
        self.status[self.proc] = 2

    def update_runstats(self, result = 0):
        if result == 0:
            self.stats[0] += 1
        elif result == 1:
            self.stats[1] += 1

    def update_jobstats(self, jobtype):
        if jobtype == "Mapper":
            self.stats[2] += 1
        elif jobtype == "Assembler":
            self.stats[3] += 1
        elif jobtype == "AssemblyChecker":
            self.stats[4] += 1
        elif jobtype == "Finisher":
            self.stats[5] += 1
