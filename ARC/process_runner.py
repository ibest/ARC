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
from Queue import Empty
from multiprocessing import Process, Queue
from random import randint

class ProcessRunner(Process):
  def __init__(self,ref_q):
    super(ProcessRunner, self).__init__()
    self.ref_q = ref_q

  def run(self):
    while True:
      try:
        item = self.ref_q.get_nowait()
        job = item['runner']
        print "[%s] Processing: %s" % (self.name,item['name'])
        job.start()
        if not job.error:
          for next_job in job.next:
            self.ref_q.put(next_job)
        else:
          # Log the error ??Recover or die??
          print "[%s] ERROR"

      except Empty:
        print "Queue Empty!!!!!"
        return
      else:
        self.ref_q.task_done()