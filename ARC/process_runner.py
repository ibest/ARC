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
import time
from Queue import Empty
from multiprocessing import Process
from ARC import logger
from ARC import exceptions
import ARC.runners
# from ARC.runners import Assembler
# from ARC.runners import AssemblyChecker
# from ARC.runners import Mapper
# from ARC.runners import Finisher
# from ARC.runners import Mapper


class ProcessRunner(Process):
    def __init__(self, job_q, result_q, finished, proc):
        super(ProcessRunner, self).__init__()
        self.job_q = job_q
        self.result_q = result_q
        self.finished = finished
        self.proc = proc
        #self.numjobs = 0
        self.retired = False

    def launch(self):
        # Moving this up into launch to see if running outside the scope of the
        # run function will release any reverence to the object and help the
        # garbage collecter clean up.
        item = self.job_q.get_nowait()
        # If we made it this far, we have found something on the
        # queue so we need to make sure we let the spawner know we
        # are not done prior to starting so spawner doesn't kill the
        # process
        self.not_done()

        job = getattr(ARC.runners, item['runner'])(item['params'])
        logger.debug("[%s] Processing: %s" % (self.name, job.message()))
        job.queue(self.job_q)
        job.start()
        job.clean()
        del job
        job = None
        del item
        item = None

    def run(self):
        """
        run() will initially sleep for 0.5 seconds, if an item is then found
        on the job_q, it will process items off of the job_q every .1 second
        until the job_q is empty, at which point it will get an Empty
        exception, and set the sleeptime to 5 seconds.
        """
        sleeptime = 0.5
        while True:
            try:
                time.sleep(sleeptime)
                if not self.retired:
                    self.launch()
                    self.result_q.put({"status": 0, "process": self.name})
                    sleeptime = 0

            except exceptions.RerunnableError as e:
                logger.warn("[%s] A job needs to be rerun: %s" % (self.name, e))
                self.result_q.put({"status": 1, "process": self.name})
            except exceptions.FatalError as e:
                logger.error("[%s] A fatal error occured: %s" % (self.name, e))
                self.result_q.put({"status": 2, "process": self.name})
            except Empty:
                # Since we aren't allowing the process to exit until the spawner
                # don't report the status if we are already done
                sleeptime = 5
                #print "got Empty exception, sleeptime", sleeptime
                if not self.is_done():
                    logger.debug("[%s] The queue is empty" % (self.name))
                    self.result_q.put({"status": 3, "process": self.name})
                    self.done()
            except (KeyboardInterrupt, SystemExit):
                logger.debug("Process interrupted")
                sys.exit()
            except Exception as e:
                logger.error("An unhandled exception occured")
                self.result_q.put({"status": 2, "process": self.name})
                raise e

    def done(self):
        self.finished[self.proc] = 1

    def not_done(self):
        self.finished[self.proc] = 0

    def is_done(self):
        if self.finished[self.proc] == 1:
            return True
        else:
            return False
