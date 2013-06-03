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


class ProcessRunner(Process):
    def __init__(self, ref_q, result_q, finished, proc):
        super(ProcessRunner, self).__init__()
        self.ref_q = ref_q
        self.result_q = result_q
        self.finished = finished
        self.proc = proc

    def run(self):
        while True:
            try:
                item = self.ref_q.get_nowait()
                job = item['runner']
                logger.debug("[%s] Processing: %s" % (self.name, item['message']))
                job.queue(self.ref_q)
                job.start()
                self.result_q.put({"status": 0, "process": self.name})
            except exceptions.RerunnableError as e:
                logger.warn("[%s] A job needs to be rerun: %s" % (self.name, e))
                self.result_q.put({"status": 1, "process": self.name})
            except exceptions.FatalError as e:
                logger.error("[%s] A fatal error occured: %s" % (self.name, e))
                self.result_q.put({"status": 2, "process": self.name})
            except Empty:
                # Since we aren't allowing the process to exit until the spawner
                # don't report the status if we are already done
                if not self.is_done():
                    logger.debug("[%s] The queue is empty" % (self.name))
                    self.result_q.put({"status": 3, "process": self.name})
                    self.done()
            except (KeyboardInterrupt, SystemExit):
                logger.debug("Process interrupted")
                sys.exit()
            else:
                self.not_done()

    def done(self):
        self.finished[self.proc] = 1

    def not_done(self):
        self.finished[self.proc] = 0

    def is_done(self):
        if self.finished[self.proc] == 1:
            return True
        else:
            return False
