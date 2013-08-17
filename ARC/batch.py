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
import resource
import multiprocessing
import gc
import logging
from copy import deepcopy
from multiprocessing import Manager
from multiprocessing import Lock
from ARC import Job
from ARC import logger


class Batch:
    def __init__(self, procs=multiprocessing.cpu_count()):
        self.bq = BatchQueues()
        self.runners = []
        self.procs = procs
        self.procs_available = procs
        self.processes = {}

    def queue(self):
        return self.bq

    def execute(self):
        jobs = self.bq.backfill(self.procs_available)
        for job in jobs:
            job.setstart()
            self.reserve(job)
            self.bq.execute(job)
            self.spawn(job)

        if jobs == []:
            return False
        else:
            return True

    def complete(self):
        jobs = []
        for job in self.bq.exec_queue_jobs():
            process = self.processes[job.ident]
            if not process.is_alive():
                job.setcomplete(process.exitcode)
                self.manage_exitcode(job)
                self.remove_process(job)
                # Find something different
                self.bq.complete(job)
                self.release(job)
                self.waiting = False
                jobs.append(job)

        if jobs == []:
            return False
        else:
            return True

    def manage_exitcode(self, job):
        code = job.exitcode

        if code == Job.OK:
            # update stats
            self.debug("Job finished with status OK")
        elif code == Job.FATALERROR:
            self.debug("A fatal error has occurred. Exiting")
            # Kill everything and exit
        elif code == Job.TIMEOUTERROR:
            self.debug("Timeout occured.  Dependent jobs are safe.")
            # update stats
        elif code == Job.RERUNERROR:
            self.warn("Deprecated. Use the resubmit method.")
            self.debug("%s - Requeuing" % (job.ident))
            job.reset()
            self.bq.rerun(job)
            # update stats
        elif code == Job.UNKNOWNERROR:
            self.debug("An unknown error occured.")
            # Kill all dependent jobs
            # update stats
        elif code == Job.PROCESSERROR:
            self.debug("The runner subprocess failed.")
            # self.bq.killall()
            # self.bq.drain()
        else:
            self.error("Unhandled return code %d from %s.  Killing all jobs." % (code, job.ident))
            self.killall()
            self.bq.drain()
            # update stats

    def spawn(self, job):
        process = self.create_process(job)
        process.setup()
        process.daemon = False
        self.debug("Starting runner[%s]" % (job.ident))
        process.start()

    def run(self):
        interval = 0
        sleeptime = 0.05
        while self.idle_notempty() or self.exec_notempty():
            # discard dead processes (active_children has the side effect of
            # joining finished processes)
            multiprocessing.active_children()
            gc.collect()

            executed = self.execute()
            completed = self.complete()
            interval += sleeptime
            time.sleep(sleeptime)
            #if executed or completed:
            #    self.stats()
            #if int(interval) % 2 == 0:
            #    self.resource_stats()
            #    interval = 0

    def reserve(self, job):
        self.procs_available -= job.procs
        self.debug("(Reserved) %d processors available" % (self.procs_available))

    def release(self, job):
        self.procs_available += job.procs
        self.debug("(Released) %d processors available" % (self.procs_available))

    def idle_notempty(self):
        return self.bq.num_idle() > 0

    def exec_notempty(self):
        return self.bq.num_execution() > 0

    def resource_stats(self):
        ru_self = resource.getrusage(resource.RUSAGE_SELF)
        ru_children = resource.getrusage(resource.RUSAGE_CHILDREN)
        memuse = (ru_self[2] + ru_children[2]) * resource.getpagesize() / 1024 / 1024
        self.log(str(ru_self))
        self.log(str(ru_children))
        self.log(str(memuse))

    def stats(self):
        self.debug("There are %d jobs in the idle queue" % (
            self.bq.num_idle()))
        self.debug("There are %d jobs in the exec queue" % (
            self.bq.num_execution()))
        self.debug("There are %d jobs in the comp queue" % (
            self.bq.num_complete()))

    def create_process(self, job):
        """
            Creates the process runner and stores it in the process hash
            identified by the unique identifier.  Returns the process runner
            object.

            :param job: The job that the process is associated with.
        """
        process = job.runner(
            job.ident,
            job.procs,
            job.params,
            self.bq)
        self.processes[job.ident] = process
        return process

    def remove_process(self, job):
        """
            Remove the process runner from the hash.

            :param job: The job that the process is associated with.
        """
        self.processes[job.ident].delete()
        del self.processes[job.ident]
        self.debug("processes: %d" % (len(self.processes.keys())))

    def killall(self):
        self.warn("KILLALL RECEIVED!  Terminating running jobs")
        for ident in self.processes:
            process = self.processes[ident]
            if process.is_alive():
                process.terminate()

    def log(self, msg):
        logger.info("Batch: %s" % (msg))

    def debug(self, msg):
        logger.debug("Batch: %s" % (msg))

    def warn(self, msg):
        logger.warn("Batch: %s" % (msg))

    def error(self, msg):
        logger.error("Batch: %s" % (msg))


class BatchQueues(object):
    def __init__(self):
        """
            Initialize the queue management system.  Creates multiprocessing
            SyncManager queues for storing idle, executing and completed
            processes.  Adds locking mechanisms and a hash to store the
            running process information.
        """
        self.mgr = Manager()
        self.idle_queue = self.mgr.list()
        self.exec_queue = self.mgr.list()
        self.comp_queue = self.mgr.list()
        self.jobs = self.mgr.dict()
        self.deps = self.mgr.list()
        self.lock = Lock()
        self.loglevel = logger.level()
        # Add globals for things that need to be passed around everywhere
        self.globals = self.mgr.dict()

    def add_global(self, key, value):
        self.globals[key] = value

    def backfill(self, procs):
        """
            Finds job(s) that will run on the available processors.  This
            is a firstfit backfill method and has shown to favor jobs
            requiring the least number of processors in a randomly distributed
            population of processor requirements. Removes the job from the
            idle queue and returns an array of jobs that will fit on the
            number of processors specified.

            :param procs: the number of processors available.
        """
        available = procs
        runjobs = []
        with self.lock:
            for index, jobid in enumerate(self.idle_queue):
                job = self.jobs[jobid]
                if job.procs <= available and self.runnable(job):
                    self.idle_queue.remove(jobid)
                    runjobs.append(job)
                    available -= job.procs
                if available <= 0:
                    break
        return runjobs

    def runnable(self, job):
        """
            Check for any outstanding dependancies.  Returns true or False

            :param job: the job to be checked
        """
        count = 0
        for ident in self.deps:
            if ident in job.deps:
                count += 1

        if count > 0:
            self.debug("%s waiting on %d dependancies" % (
                job.ident, count))
            return False
        else:
            self.debug("%s dependancies have been satisfied" % (job.ident))
            return True

    def submit(self, runner, **kwargs):
        """
            Submit a job to the idle queue

            :param runner: The process runner to be used
            :param kwargs: An argument hash containing the job details
        """
        # kwargs = deepcopy(kwargs)
        job = Job(runner, **kwargs)
        self.jobs[job.ident] = job
        # Take deps and put them in there own list so we don't need to
        # run through all completed
        for ident in job.deps:
            self.deps.append(ident)

        self.debug("Adding %s(%s)" % (runner.__name__, job.ident))
        with self.lock:
            self.idle_queue.append(job.ident)
        return job

    def execute(self, job):
        with self.lock:
            self.exec_queue.append(job.ident)

    def drain(self):
        with self.lock:
            for index, jobid in enumerate(self.idle_queue):
                self.idle_queue.pop(index)
                self.error("Removed waiting job from the queue")

    def complete(self, job):
        with self.lock:
            self.exec_queue.remove(job.ident)
            #self.comp_queue.append(job.ident)
            # Remove me from deps so my favorite process can run
            try:
                self.deps.remove(job.ident)
            except ValueError:
                # Do nothing.
                pass
            del(self.jobs[job.ident])
            self.debug("idle: %d | exec: %d | jobs: %d | deps: %d" % (
                len(self.idle_queue),
                len(self.exec_queue),
                len(self.jobs),
                len(self.deps)
            ))

    def exec_queue_jobs(self):
        jobs = []
        for ident in self.exec_queue:
            jobs.append(self.jobs[ident])
        return jobs

    def rerun(self, job):
        self.remove_process(job)
        with self.lock:
            self.idle_queue.append(job)

    def num_idle(self):
        with self.lock:
            return len(self.idle_queue)

    def num_execution(self):
        with self.lock:
            return len(self.exec_queue)

    def num_complete(self):
        with self.lock:
            return len(self.comp_queue)

    def log(self, msg):
        logger.info("BatchQueue: %s" % (msg))

    def debug(self, msg):
        if self.loglevel == logging.DEBUG:
            logger.debug("BatchQueue: %s" % (msg))

    def warn(self, msg):
        logger.warn("BatchQueue: %s" % (msg))

    def error(self, msg):
        logger.error("BatchQueue: %s" % (msg))
