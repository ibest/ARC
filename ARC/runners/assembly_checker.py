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

import os
import time
from ARC.runner import Finisher
from ARC.runner import BaseRunner


class AssemblyChecker(BaseRunner):
    """
        Checks for "finished" files in each of the assembly folders. Set values
        of params['assemblies'] to True for all finished assemblies. If all
        assemblies are finished, kick off the finisher process.  Now that we
        have job dependency, this will probably be deprecated.

        required params:
            'targets': dictionary, keys are paths, values are boolean
    """
    def execute(self):
        """
            run through list of targets, check any that haven't finished
            already
        """
        sample = self.params['sample']
        completed = sum(self.params['targets'].values())

        msg = "Sample: %s AssemblyChecker started with " % (sample)
        msg += "%s of %s targets completed" % (
            sample, completed, len(self.params['targets']))
        self.info(msg)

        for target_folder in self.params['targets']:
            if not self.params['targets'][target_folder]:
                file = os.path.join(target_folder, 'finished')
                if os.path.exists(file):
                    self.params['targets'][target_folder] = True
                    self.info("%s exists" % file)
                    completed += 1

        # Now check whether all have finished, if not, add a new
        # AssemblyChecker to the queue
        if len(self.params['targets']) > sum(self.params['targets'].values()):
            # sleep 4 seconds before submitting
            time.sleep(5)
            self.submit(
                AssemblyChecker,
                procs=1,
                params=self.params)
            msg = "Sample: %s Assemblies not finished: " % (sample)
            msg += "%s of %s targets completed" % (
                completed, len(self.params['targets']))
            self.info(msg)
        else:
            self.submit(
                Finisher,
                procs=1,
                params=self.params)
            msg = "Sample: %s Assemblies finished: " % (sample)
            msg += "%s of %s targets completed" % (
                completed, len(self.params['targets']))
            self.info(msg)
