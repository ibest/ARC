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

import time
import multiprocessing
#import os
#from copy import deepcopy
from Queue import Empty
from ARC import ProcessRunner
from ARC import logger
from ARC import exceptions
from ARC.runners import Mapper


class Spawn:
    def __init__(self, config):
        self.workers = []
        self.config = config
        self.nprocs = int(config['nprocs'])
        self.q = multiprocessing.JoinableQueue()
        # Contains the state of the job
        # 0: Not set
        # 1: Waiting for jobs
        # 2: Running a job
        self.status = multiprocessing.Array('i', [0] * self.nprocs)
        # Contains stats for the run
        # [0]: Number of jobs returned ok
        # [1]: Number of jobs returned rerun
        # [2]: Number of Mapper jobs run
        # [3]: Number of Assembly jobs run
        # [4]: Number of Checker jobs run
        # [5]: Number of Finisher jobs run
        self.stats = multiprocessing.Array('i', [0] * self.nprocs)

    def submit(self):
        # Get the number of samples from the configuration
        logger.info("Submitting initial mapping runs.")

        for sample in self.config['Samples']:
            s = self.config['Samples'][sample]
            params = {}
            for k in self.config:
                params[k] = self.config[k]
            params['working_dir'] = s['working_dir']
            params['finished_dir'] = s['finished_dir']
            params['sample'] = sample

            if 'PE1' in s and 'PE2' in s:
                params['PE1'] = s['PE1']
                params['PE2'] = s['PE2']
            if 'SE' in s:
                params['SE'] = s['SE']

            # mapper = Mapper(params)
            self.q.put(Mapper.to_job(params))

    def run(self):
        logger.info("Starting...")
        logger.debug("Setting up workers.")

        for i in range(self.nprocs):
            worker = ProcessRunner(
                i,
                self.q,
                self.status,
                self.stats)
            self.workers.append(worker)
            worker.daemon = False
            worker.start()

        while True:
            try:
                self.q.join()
                # This shouldn't be needed but we will check just in case
                if self.all_workers_waiting():
                    logger.debug("Workers are all waiting and the queue is empty.  Exiting")
                    break
                else:
                    logger.debug("Workers are not in a waiting state.  Waiting for more.")
                    time.sleep(5)
                # Terminate the jobs

            except exceptions.FatalError:
                logger.error("A fatal error was encountered.")
                raise
            except (KeyboardInterrupt, SystemExit):
                logger.error("Terminating processes")
                self.kill_workers(workers)
                raise

        # Kill 'em all!
        self.killall()

        logger.info("-----")
        logger.info("%d processes returned ok." % (self.stats[0]))
        logger.info("%d processes had to be rerun." % (self.stats[1]))
        logger.info("-----")
        logger.info("%d Mapper jobs run." %(self.stats[2]))
        logger.info("%d Assembly jobs run." %(self.stats[3]))
        logger.info("%d Checker jobs run." %(self.stats[4]))
        logger.info("%d Finisher jobs run." %(self.stats[5]))

    def killall(self, workers):
        for worker in workers:
            logger.debug("Shutting down %s" % (worker.name))
            worker.terminate()
            worker.join()

    def all_workers_waiting(self):
        waiting = 0
        for worker in workers:
            if worker.is_waiting():
                waiting += 1

        return waiting == len(workers)
