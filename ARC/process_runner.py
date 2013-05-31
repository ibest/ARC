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

from Queue import Empty
from multiprocessing import Process
from ARC import logger
from ARC import exceptions


class ProcessRunner(Process):
    def __init__(self, ref_q):
        super(ProcessRunner, self).__init__()
        self.ref_q = ref_q

    def run(self):
        while True:
            try:
                # logger.info("The queue currently contains %d jobs" % (self.ref_q.qsize()))
                item = self.ref_q.get_nowait()
                job = item['runner']
                logger.info("[%s] Processing: %s" % (self.name, item['message']))
                job.queue(self.ref_q)
                job.start()
            except Empty:
                return
            except exceptions.FatalError as e:
                logger.error("[%s] A fatal error occured: %s" % (self.name, e))
                raise
            except exceptions.RerunnableError as e:
                logger.error("[%s] An error occured: %s" % (self.name, e))
                self.ref_q.task_done()
            else:
                self.ref_q.task_done()
