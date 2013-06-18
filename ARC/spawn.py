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
        params = {
            'reference': config['reference'],
            'numcycles': config['numcycles'],
            'working_dir': os.path.realpath('./working_' + sample),
            'sample': sample,
            'mapper': config['mapper'],
            'assembler': config['assembler'],
            'format': config['format'],
            'verbose': config['verbose'],
            'iteration': config['iteration']
        }
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

    while True:
        try:
            time.sleep(0.5)
            result = result_q.get_nowait()
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
            else:
                logger.error("Unknown state returned from %s" % (result['process']))
                kill_workers(workers)
        except (KeyboardInterrupt, SystemExit):
            logger.error("Terminating processes")
            kill_workers(workers)
            raise
        except Empty:
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


def not_done(finished):
    logger.debug("Checking to see if workers are done")
    done = 0
    for i in finished:
        done += i
    logger.debug("Active Workers: %d" % (len(finished) - done))
    return done < len(finished)


def not_empty(q):
    return not q.empty()
