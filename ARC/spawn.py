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
import os
from copy import deepcopy
from Queue import Empty
from process_runner import ProcessRunner
from ARC import logger
from ARC import exceptions
from ARC.mapper import MapperRunner


def run(config):
    logger.info("Starting...")

    # For test
    #nprocs = 10
    # Get the number of processors to use
    nprocs = int(config['nprocs'])

    ref_q = multiprocessing.Queue()
    result_q = multiprocessing.Queue()
    finished = multiprocessing.Array('i', [0]*nprocs)

    # from test import TestRunner
    # for i in range(50):
    #     s = TestRunner({'foo': i})
    #     ref_q.put(s.to_dict())

    # Get the number of samples from the configuration
    #'sample' in params and 'reference' in params and 'working_dir' in params and (('PE1' in params and 'PE2' in params) or 'SE' in params))
    for sample in config['Samples']:
        s = config['Samples'][sample]
        params = deepcopy(config)
        params['working_dir'] = s['working_dir']
        params['finished_dir'] = s['finished_dir']
        params['sample'] = sample

        if 'PE1' in s and 'PE2' in s:
            params['PE1'] = s['PE1']
            params['PE2'] = s['PE2']
        if 'SE' in s:
            params['SE'] = s['SE']

        mapper = MapperRunner(params)
        ref_q.put(mapper.to_dict())

    # samples = config['samples'] is dict
    # Get the target file from the configuration
    # target = config['target']
    # Load the queue with the intial mapping jobs
    #
    # for sample in samples
    #   mapper = Mapper(sample,target)
    #   ref_q.put(mapper.to_dict()) # {'runner': self, 'message': 'Map <Sample name>'}

    # Need signal handling for graceful exit

    logger.debug("Setting up workers.")
    workers = []
    for i in range(nprocs):
        worker = ProcessRunner(ref_q, result_q, finished, i)
        worker.daemon = False
        workers.append(worker)
        worker.start()

    status_ok = 0
    status_rerun = 0
    sleeptime = 0.1
    while True:
        try:
            time.sleep(sleeptime)
            result = result_q.get_nowait()
            #print "Spawner, setting sleeptime to %s" % sleeptime
            sleeptime = 0
            if result['status'] == 0:
                logger.debug("Completed successfully %s" % (result['process']))
                status_ok += 1
            elif result['status'] == 1:
                logger.debug("Rerunnable error returned from %s" % (result['process']))
                status_rerun += 1
            elif result['status'] == 2:
                logger.error("Fatal error returned from %s" % (result['process']))
                logger.error("Terminating processes")
                kill_workers(workers)
                raise exceptions.FatalError("Unrecoverable error")
            elif result['status'] == 3:
                logger.debug("Empty state returned from %s" % (result['process']))
            elif result['status'] == 4:
                logger.debug("%s worker needs to be retired" % result['process'])
                retire_worker(workers, result['process'], ref_q, result_q, finished)
            else:
                logger.error("Unknown state returned from %s" % (result['process']))
                kill_workers(workers)
        except (KeyboardInterrupt, SystemExit):
            logger.error("Terminating processes")
            kill_workers(workers)
            raise
        except Empty:
            sleeptime = 5
            #print "Spawn setting sleeptime to", sleeptime
            if not not_done(finished):
                logger.debug("Results queue is empty and there are no active processes.  Exiting")
                break
            else:
                logger.debug("Results queue is empty and there are still active processes.  Waiting")

    logger.info("%d processes returned ok" % (status_ok))
    logger.info("%d processes had to be rerun" % (status_rerun))
    kill_workers(workers)


def kill_workers(workers):
    for worker in workers:
        logger.debug("Shutting down %s" % (worker.name))
        worker.terminate()

def retire_worker(workers, process, ref_q, result_q, finished):
    for i in range(len(workers)):
        worker = workers[i]
        if worker.name == process:
            logger.info("Terminating working %s" % worker.name)
            worker.terminate()
            worker = ProcessRunner(ref_q, result_q, finished, i)
            worker.daemon = False
            workers[i] = worker
            worker.start()
            logger.info("Started new worker %s" % worker.name)

def not_done(finished):
    logger.debug("Checking to see if workers are done")
    done = 0
    for i in finished:
        done += i
    logger.debug("Active Workers: %d" % (len(finished) - done))
    return done < len(finished)


def not_empty(q):
    return not q.empty()
