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
import sys
import traceback
# from copy import deepcopy
from ARC import exceptions
from ARC import logger
from ARC.runners import Base
from ARC.runners import Finisher
# from ARC.runner import AssemblyChecker


class AssemblyChecker(Base):
    """
    Checks for "finished" files in each of the assembly folders. Set values of params['assemblies'] to True for
    all finished assemblies. If all assemblies are finished, kick off the finisher process.
    required params:
        'targets': dictionary, keys are paths, values are boolean
    """

    def message(self):
        return 'Starting AssemblyChecker for sample %s' % self.params['sample']

    def start(self):
        """ run through list of targets, check any that haven't finished already """
        try:
            sample = self.params['sample']
            completed = sum(self.params['targets'].values())
            logger.info("Sample: %s AssemblyChecker started with %s of %s targets completed" % (sample, completed, len(self.params['targets'])))
            for target_folder in self.params['targets']:
                if not self.params['targets'][target_folder]:
                    f = os.path.join(target_folder, 'finished')
                    if os.path.exists(f):
                        self.params['targets'][target_folder] = True
                        logger.info("%s exists" % f)
                        completed += 1
            #Now check whether all have finished, if not, add a new AssemblyChecker to the queue
            if len(self.params['targets']) > sum(self.params['targets'].values()):
                #some jobs haven't completed yet
                checker_params = {}
                for k in self.params:
                    checker_params[k] = self.params[k]
                #checker_params = deepcopy(self.params)
                # checker = AssemblyChecker(checker_params)
                time.sleep(5)  # sleep 4 seconds before putting a checker back on the job_q
                self.submit(AssemblyChecker.to_job(checker_params))
                logger.info("Sample: %s Assemblies not finished: %s of %s targets completed" % (sample, completed, len(self.params['targets'])))
            else:
                params = {}
                for k in self.params:
                    params[k] = self.params[k]
                # params = deepcopy(self.params)
                # finisher = Finisher(params)
                self.submit(Finisher.to_job(params))
                logger.info("Sample: %s Assemblies finished: %s of %s targets completed" % (sample, completed, len(self.params['targets'])))
                time.sleep(1)  # test for occasionally not kicking off finisher properly.
        except:
            print "".join(traceback.format_exception(*sys.exc_info()))
            raise exceptions.FatalError("".join(traceback.format_exception(*sys.exc_info())))
