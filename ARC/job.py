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
import uuid
import time


class Job:
    IDLE = 0
    RUNNING = 1
    RERUN = 2
    COMPLETE = 3

    OK = 0
    FATALERROR = 255
    TIMEOUTERROR = 249
    RERUNERROR = 248
    UNKNOWNERROR = 247
    PROCESSERROR = 246

    def __init__(self, runner, **kwargs):
        # A unique ID for the job
        self.ident = str(uuid.uuid4())
        # The runner class that will be used for the job
        self.runner = runner
        # Keyword arguments for job configuration
        self.procs = kwargs.pop('procs', 1)
        self.priority = kwargs.pop('priority', 0)
        self.timeout = kwargs.pop('timeout', -1)
        self.params = kwargs.pop('params', {})
        self.deps = kwargs.pop('deps', [])
        # States and times
        self.exitcode = -1
        self.state = self.IDLE
        self.enquetime = time.time()
        self.starttime = 0
        self.endtime = 0
        self.executions = 0
        self.pid = None

    def __repr__(self):
        return self.ident

    def show(self):
        print "id: %s" % (self.ident)
        print "runner: %s" % (self.runner.__name__)
        print "procs: %d" % (self.procs)

    def reset(self):
        self.state = self.IDLE
        self.enquetime = time.time()
        self.starttime = 0
        self.endtime = 0
        self.executions += 1

    def clean(self):
        del(self.runner)
        del(self.params)
        del(self.deps)

    def setargs(self, **kwargs):
        self.procs = kwargs.pop('procs', self.procs)
        self.priority = kwargs.pop('priority', self.priority)
        self.timeout = kwargs.pop('timeout', self.timeout)
        self.params = kwargs.pop('params', self.params)
        self.deps = kwargs.pop('deps', self.deps)

    def setrerun(self):
        self.state = self.RERUN

    def setstart(self):
        self.state = self.RUNNING
        self.starttime = time.time()

    def setcomplete(self, exitcode):
        self.exitcode = exitcode
        self.state = self.COMPLETE
        self.endtime = time.time()
        self.clean()
