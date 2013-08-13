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
from Queue import Empty
from multiprocessing import Process
import time
#import os
from ARC import logger
from ARC import exceptions


class ProcessRunner(Process):
    def __init__(self, ref_q, result_q, finished, proc):
        super(ProcessRunner, self).__init__()
        self.ref_q = ref_q
        self.result_q = result_q
        self.finished = finished
        self.proc = proc
        self.numjobs = 0
        self.retired = False

    def run(self):
        """
        run() will initially sleep for .5 seconds, if an item is then found on the ref_q, it will process items off of the ref_q
        every .1 second until the ref_q is empty, at which point it will get an Empty exception, and set the sleeptime to 5 seconds.
        """
        sleeptime = .5
        while True:
            try:
                time.sleep(sleeptime)
                if not self.retired:
                    item = self.ref_q.get_nowait()
                    sleeptime = 0
                    #print "got Item", item['runner'], "sleep time", sleeptime
                    # If we made it this far, we have found something on the
                    # queue so we need to make sure we let the spawner know we
                    # are not done prior to starting so spawner doesn't kill the
                    # process
                    self.not_done()
                    # Begin the run
                    job = item['runner']
                    logger.debug("[%s] Processing: %s" % (self.name, item['message']))
                    job.queue(self.ref_q)
                    self.numjobs += 1
                    job.start()
                    #print "%s finished a job, total jobs %s" % (self.name, self.numjobs)
                    self.result_q.put({"status": 0, "process": self.name})
                if str(job.__class__) == 'ARC.mapper.MapperRunner':
                    logger.info(self.name + " got a MapperRunner job, asking to be retired after %s jobs" % self.numjobs)
                    self.result_q.put({"status": 4, "process": self.name})
                    sleeptime = 5
                    self.retired = True
                # if self.numjobs > 10 and not self.retired:
                #     #Ask for retirement
                #     logger.debug(self.name + " Asking to be retired after %s jobs" % self.numjobs)
                #     self.result_q.put({"status": 4, "process": self.name})
                #     sleeptime = 5
                #     self.retired = True
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
            # else:
                # self.not_done()

    def done(self):
        self.finished[self.proc] = 1

    def not_done(self):
        self.finished[self.proc] = 0

    def is_done(self):
        if self.finished[self.proc] == 1:
            return True
        else:
            return False
