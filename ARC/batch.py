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
from multiprocessing import Manager
from multiprocessing import Lock
from ARC import Job
from ARC import logger


class Batch:
    def __init__(self, bq, procs=multiprocessing.cpu_count()):
        self.bq = bq
        self.runners = []
        self.procs = procs
        self.procs_available = procs
        self.processes = []

    def execute(self):
        jobs = self.bq.find_backfill(self.procs_available)
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
        jobs = self.bq.find_completed()
        for job in jobs:
            if job.state == Job.RERUN:
                self.debug("%s - Requeuing" % (job.ident))
                job.reset()
                self.bq.rerun(job)
            else:
                process = self.bq.process(job)
                print "HIIIIIIIIII: %d" % process.exitcode
                job.setcomplete(process.exitcode)
                self.manage_exitcode(job)
                self.bq.complete(job)

            self.release(job)
            self.waiting = False

        if jobs == []:
            return False
        else:
            return True

    def manage_exitcode(self, job):
        code = job.exitcode

        if code == Job.OK:
            # update stats
            self.debug("OK")
        elif code == Job.FATALERROR:
            self.debug("A fatal error has occurred. Exiting")
            # Kill everything and exit
        elif code == Job.TIMEOUTERROR:
            self.debug("Timeout occured.  Killing dependent jobs.")
            # Kill all dependent jobs
            # update stats
        elif code == Job.RERUNERROR:
            self.warn("Deprecated. Use the resubmit method.")
            # Kill all dependent jobs
            # update stats
        elif code == Job.UNKNOWNERROR:
            self.debug("An unknown error occured.")
            # Kill all dependent jobs
            # update stats
        elif code == Job.PROCESSERROR:
            self.debug("The runner subprocess failed.")
            # Kill all dependent jobs
            # update stats
        else:
            self.error("Unhandled return code %d.  Killing dependent jobs." % (code))
            # Kill all dependent jobs
            # update stats

    def spawn(self, job):
        process = self.bq.process(job)
        process.setup()
        process.daemon = False
        # self.debug("Starting runner[%s]" % (job.ident))
        process.start()

    def run(self):
        self.stats()

        while self.idle_notempty() or self.exec_notempty():
            executed = self.execute()
            completed = self.complete()
            time.sleep(0.5)
            if executed or completed:
                self.stats()

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

    def stats(self):
        self.debug("There are %d jobs in the idle queue" % (
            self.bq.num_idle()))
        self.debug("There are %d jobs in the exec queue" % (
            self.bq.num_execution()))
        self.debug("There are %d jobs in the comp queue" % (
            self.bq.num_complete()))

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
        self.lock = Lock()
        self.processes = {}

    def index(self, ident, queue):
        """
            Find the index value of a job in the queue.  Does not perform
            any locking, so the calling function should lock this.  Returns
            an integer value representing the index of the job in the array
            or -1 if it does not exist.  This should probably raise an
            exception since -1 is a perfectly legal value when accessing
            an item in an array.

            :param ident: the unique identifier of the job (uuid4 format)
            :param queue: the SyncManager queue to search
        """
        for index, job in enumerate(queue):
            if job.ident == ident:
                return index
        return -1

    def pop(self, ident, queue):
        """
            Pops a job off of a queue.  Returns the job or None if not found.

            :param ident: the unique identifier of the job (uuid4 format)
            :param queue: the SyncManager queue to search
        """
        with self.lock:
            index = self.index(ident, queue)
            if index >= 0:
                return queue.pop(index)
            else:
                return None

    def get(self, ident, queue):
        """
            Get a job from the queue.  Returns a reference to the job or
            None if not found.

            :param ident: the unique identifier of the job (uuid4 format)
            :param queue: the SyncManager queue to search
        """
        with self.lock:
            index = self.index(ident, queue)
            if index >= 0:
                return queue[index]
            else:
                return None

    def find_backfill(self, procs):
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
        jobs = []
        with self.lock:
            for index, job in enumerate(self.idle_queue):
                if job.procs <= available and self.runnable(job):
                    jobs.append(self.idle_queue.pop(index))
                    available -= job.procs
                if available <= 0:
                    break
        return jobs

    def find_completed(self):
        """
            Find jobs that are no longer running.  Removes the jobs from
            the execution queue and returns an array will all the completed
            jobs.
        """
        completed = []
        with self.lock:
            for index, job in enumerate(self.exec_queue):
                if not self.process(job).is_alive():
                    completed.append(self.exec_queue.pop(index))
        return completed

    def runnable(self, job):
        """
            Check for any outstanding dependancies.  Returns true or False

            :param job: the job to be checked
        """
        count = 0
        for completed in self.comp_queue:
            if completed.ident in job.deps:
                count += 1
        completed_deps = len(job.deps) - count

        if completed_deps > 0:
            self.debug("%s waiting on %d dependancies" % (
                job.ident, completed_deps))
            self.debug(job.deps)
            return False
        else:
            self.debug("%s dependancies have been satisfied" % (job.ident))
            return True

    def create_process(self, job):
        """
            Creates the process runner and stores it in the process hash
            identified by the unique identifier.  Returns the process runner
            object.

            :param job: The job that the process is associated with.
        """
        self.processes[job.ident] = job.runner(
            job.ident,
            job.procs,
            job.params,
            self)
        return self.processes[job.ident]

    def remove_process(self, job):
        """
            Remove the process runner from the hash.

            :param job: The job that the process is associated with.
        """
        del self.processes[job.ident]

    def process(self, job):
        """
            Return a process runner from the hash.

            :param job: The job that the process is associated with.
        """
        return self.processes[job.ident]

    def submit(self, runner, **kwargs):
        """
            Submit a job to the idle queue

            :param runner: The process runner to be used
            :param kwargs: An argument hash containing the job details
        """
        job = Job(runner, **kwargs)
        self.debug("Adding %s(%s)" % (runner.__name__, job.ident))
        self.debug("Deps %s %s" % (job.ident, job.deps))
        with self.lock:
            self.idle_queue.append(job)
        return job

    def resubmit(self, ident, **kwargs):
        """
            Resubmit a running job.  Places the job with any modified
            details back into the rerun state and lets the former job run to
            completion.  For some reason the job must be popped off the queue,
            modified and placed back in order for the changes to stick.

            :param runner: The process runner to be used
            :param kwargs: An argument hash containing the job details
        """
        job = self.pop(ident, self.exec_queue)
        job.setrerun()
        job.setargs(**kwargs)
        self.debug("Resubmitting %s, %s" % (job.__class__.__name__, job.ident))
        with self.lock:
            self.exec_queue.append(job)

    def execute(self, job):
        self.create_process(job)
        with self.lock:
            self.exec_queue.append(job)

    def complete(self, job):
        # self.debug(self.comp_queue)
        self.remove_process(job)
        with self.lock:
            self.comp_queue.append(job)

    def rerun(self, job):
        self.remove_process(job)
        with self.lock:
            self.idle_queue.append(job)

    def list_idle(self):
        with self.lock:
            for i in self.idle_queue:
                print " - %s" % (i.ident)

    def list_execution(self):
        with self.lock:
            for i in self.exec_queue:
                print " - %s" % (i.ident)

    def list_complete(self):
        with self.lock:
            for i in self.comp_queue:
                print " - %s" % (i.ident)

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
        logger.debug("BatchQueue: %s" % (msg))

    def warn(self, msg):
        logger.warn("BatchQueue: %s" % (msg))

    def error(self, msg):
        logger.error("BatchQueue: %s" % (msg))
