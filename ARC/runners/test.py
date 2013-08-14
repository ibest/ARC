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
from ARC.runners import ProcessBase
from random import randint


class Test(ProcessBase):
    def execute(self):
        self.log("Starting run with %d processors" % (self.procs))
        self.log("Received a value of %d" % (self.params['value']))
        self.log("Sleeping for %d seconds" % (self.params['sleep']))

        args = [
            ['sleep', str(self.params['sleep'])],
            # This one should cause dd to return a non zero value
            ["dd", "of=/dev/nulls", "if=/dev/urandom", "bs=1k", "count=10000"],
            ['sleep', str(self.params['sleep'])],
            ["dd", "of=/dev/null", "if=/dev/urandom", "bs=1k", "count=20000"],
            ['sleep', str(self.params['sleep'])],
            ["dd", "of=/dev/null", "if=/dev/urandom", "bs=1k", "count=30000"],
            ['sleep', str(self.params['sleep'])],
            ["dd", "of=/dev/null", "if=/dev/urandom", "bs=1k", "count=40000"],
            ['sleep', str(self.params['sleep'])],
            ["dd", "of=/dev/null", "if=/dev/urandom", "bs=1k", "count=50000"],
            # Popen should barf!
            ['foo', '-a 42', '/path/to/my/file']]

        self.shell(
            args[self.params['value']],
            description="Sample %d" % (self.params['num']),
            timeout=8)

        if self.params['value'] == 8:
            self.resubmit(
                params={
                    'sleep': 1,
                    'value': 1,
                    'num': self.params['num']})
            self.log("Resubmitting job")
        elif self.params['value'] == 10:
            newparams = {
                'value': randint(1, 10),
                'sleep': randint(1, 10),
                'num': self.params['num']}
            job = self.submit(
                Test,
                procs=randint(1, 4),
                params=newparams)
            self.log("Submitting new job %s" % (job.ident))

        self.log("Done")
