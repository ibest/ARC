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
import multiprocessing
from process_runner import ProcessRunner
from ARC import logger

def run(config={}):
  ref_q = multiprocessing.JoinableQueue()
  logger.info("Starting...")

  from test import TestRunner
  nprocs = 10
  for i in range(50):
    s = TestRunner()
    ref_q.put({'runner': s, 'message': 'Sample Run'})

  # Get the number of processors to use
  # nprocs = config['nprocs']
  # Get the number of samples from the configuration
  # samples = config['samples'] is dict
  # Get the target file from the configuration
  # target = config['target']
  # Load the queue with the intial mapping jobs
  #
  # for sample in samples
  #   mapper = Mapper(sample,target)
  #   ref_q.put(mapper.to_queue) # {'runner': self, 'name': 'Map <Sample name>'}

  # Need signal handling for graceful exit
  for i in range(nprocs):
    worker = ProcessRunner(ref_q)
    worker.start()

  ref_q.join()