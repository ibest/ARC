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
#import os
#from copy import deepcopy
from Queue import Empty
from ARC import ProcessQueue
from ARC import ProcessRunner
from ARC import logger
from ARC import exceptions
from ARC.runners import Mapper


class Spawn:
    def __init__(self, params, universals):
        self.nprocs = int(params['nprocs'])

        #Setup thread-safe shared objects
        # self.job_q = multiprocessing.Queue()
        # self.result_q = multiprocessing.Queue()
        # self.finished = multiprocessing.Array('i', [0] * self.nprocs)
        self.pq = ProcessQueue(self.nprocs)

        # Add the universals to the process queue
        self.pq.add_universals(universals)

        # Get the number of samples from the configuration
        logger.info("Submitting initial mapping runs.")

        for sample in params['Samples']:
            s = params['Samples'][sample]
            mapper_params = {}
            for k in params:
                mapper_params[k] = params[k]
            mapper_params['working_dir'] = s['working_dir']
            mapper_params['finished_dir'] = s['finished_dir']
            mapper_params['sample'] = sample

            if 'PE1' in s and 'PE2' in s:
                mapper_params['PE1'] = s['PE1']
                mapper_params['PE2'] = s['PE2']
            if 'SE' in s:
                mapper_params['SE'] = s['SE']

            mapper = Mapper(mapper_params)
            self.pq.job_q.put(mapper.to_dict())

    def run(self):
        logger.info("Starting...")
        logger.debug("Setting up workers.")
        workers = []
        for i in range(self.nprocs):
            worker = ProcessRunner(i, self.pq)
            worker.daemon = False
            workers.append(worker)
            worker.start()

        status_ok = 0
        status_rerun = 0
        status_timeout = 0
        sleeptime = 0.1
        while True:
            try:
                time.sleep(sleeptime)
                result = self.pq.result_q.get_nowait()
                #print "Spawner, setting sleeptime to %s" % sleeptime
                sleeptime = 0
                if result['status'] == 0:
                    logger.debug("Completed successfully %s" % (
                        result['process']))
                    status_ok += 1
                elif result['status'] == 1:
                    logger.debug("Rerunnable error returned from %s" % (
                        result['process']))
                    status_rerun += 1
                elif result['status'] == 2:
                    logger.error("Fatal error returned from %s" % (
                        result['process']))
                    logger.error("Terminating processes")
                    self.kill_workers(workers)
                    raise exceptions.FatalError("Unrecoverable error")
                elif result['status'] == 3:
                    logger.debug("Empty state returned from %s" % (
                        result['process']))
                elif result['status'] == 4:
                    logger.debug("%s worker needs to be retired" % (
                        result['process']))
                elif result['status'] == 5:
                    logger.debug("A subprocess timed out on %s" % (
                        result['process']))
                    status_timeout += 1
                else:
                    logger.error("Unknown state returned from %s" % (
                        result['process']))
                    self.kill_workers(workers)
            except (KeyboardInterrupt, SystemExit):
                logger.error("Terminating processes")
                self.kill_workers(workers)
                raise
            except Empty:
                sleeptime = 5
                #print "Spawn setting sleeptime to", sleeptime
                # In rare cases, the queue can be empty because a worker just
                # pulled a job, but hasn't yet gotten to uppdate the finished
                # status. This will cause not_done() to return False, causing
                # the loop to break. Adding a short sleep here allows workers
                # to update their status.
                time.sleep(1)
                if not self.not_done():
                    logger.debug(
                        "Results queue is empty and there are no active "
                        "processes.  Exiting")
                    break
                else:
                    logger.debug(
                        "Results queue is empty and there are still "
                        "active processes.  Waiting")

        logger.info("%d processes returned ok" % (status_ok))
        logger.info("%d processes had to be rerun" % (status_rerun))
        logger.info("%d processes timed out" % (status_timeout))
        self.kill_workers(workers)

    def kill_workers(self, workers):
        for worker in workers:
            logger.debug("Shutting down %s" % (worker.name))
            worker.terminate()

    def retire_worker(self, workers, process, pq):
        for i in range(len(workers)):
            worker = workers[i]
            if worker.name == process:
                logger.info("Terminating working %s" % worker.name)
                worker.terminate()
                worker = ProcessRunner(i, pq)
                worker.daemon = False
                workers[i] = worker
                worker.start()
                logger.info("Started new worker %s" % worker.name)

    def not_done(self):
        logger.debug("Checking to see if workers are done")
        done = 0
        for i in self.pq.finished:
            done += i
        logger.debug("Active Workers: %d" % (
            len(self.pq.finished) - done))
        return done < len(self.pq.finished)

    def not_empty(self, q):
        return not q.empty()
