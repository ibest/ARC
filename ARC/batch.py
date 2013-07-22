from multiprocessing import Manager
from multiprocessing import Lock
from multiprocessing import Process
from random import randint
import multiprocessing
import uuid
import time
import sys
import subprocess
import os
import signal


class Job:
    IDLE = 0
    RUNNING = 1
    RERUN = 2
    COMPLETE = 3

    def __init__(self, runner, **kwargs):
        # A unique ID for the job
        self.ident = str(uuid.uuid4())
        # The runner class that will be used for the job
        self.runner = runner
        # Keyword arguments for job configuration
        self.procs = kwargs.pop('procs', 1)
        self.priority = kwargs.pop('priority', 0)
        self.timeout = kwargs.pop('timeout', -1)
        self.params = kwargs.pop('params', {})
        self.deps = kwargs.pop('deps', [])
        # States and times
        self.exitcode = -1
        self.state = self.IDLE
        self.enquetime = time.time()
        self.starttime = 0
        self.endtime = 0
        self.executions = 0

    def __repr__(self):
        return self.ident

    def show(self):
        print "id: %s" % (self.ident)
        print "runner: %s" % (self.runner.__name__)
        print "procs: %d" % (self.procs)

    def reset(self):
        self.state = self.IDLE
        self.enquetime = time.time()
        self.starttime = 0
        self.endtime = 0
        self.executions += 1

    def setargs(self, **kwargs):
        self.procs = kwargs.pop('procs', self.procs)
        self.priority = kwargs.pop('priority', self.priority)
        self.timeout = kwargs.pop('timeout', self.timeout)
        self.params = kwargs.pop('params', self.params)
        self.deps = kwargs.pop('deps', self.deps)

    def setrerun(self):
        self.state = self.RERUN

    def setstart(self):
        self.state = self.RUNNING
        self.starttime = time.time()

    def setcomplete(self):
        self.state = self.COMPLETE
        self.endtime = time.time()


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
            self.log("%s waiting on %d dependancies" % (
                job.ident, completed_deps))
            self.log(job.deps)
            return False
        else:
            self.log("%s dependancies have been satisfied" % (job.ident))
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
        # print "Adding %s" % (job.ident)
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
        self.exec_queue.append(job)

    def execute(self, job):
        self.create_process(job)
        with self.lock:
            self.exec_queue.append(job)

    def complete(self, job):
        # self.log(self.comp_queue)
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
        print "BATCHQUEUE[MGR]: %s" % (msg)


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
            # self.log("[STATE] %s %d" % (job.ident, job.state))
            if job.state == Job.RERUN:
                # self.log("%s - Requeuing" % (job.ident))
                job.reset()
                self.bq.rerun(job)
            else:
                # self.log("%s - Complete" % (job.ident))
                job.setcomplete()
                self.bq.complete(job)
            self.release(job)
            self.waiting = False

        if jobs == []:
            return False
        else:
            return True

    def spawn(self, job):
        process = self.bq.process(job)
        process.setup()
        process.daemon = False
        # self.log("Starting runner[%s]" % (job.ident))
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
        # self.log("%d processors available" % (self.procs_available))

    def release(self, job):
        self.procs_available += job.procs
        # self.log("%d processors available" % (self.procs_available))

    def idle_notempty(self):
        return self.bq.num_idle() > 0

    def exec_notempty(self):
        return self.bq.num_execution() > 0

    def stats(self):
        self.log("There are %d jobs in the idle queue" % (
            self.bq.num_idle()))
        self.log("There are %d jobs in the exec queue" % (
            self.bq.num_execution()))
        self.log("There are %d jobs in the comp queue" % (
            self.bq.num_complete()))

    def log(self, msg):
        print "BATCH[MGR]: %s" % (msg)


class BaseRunner(Process):
    def __init__(self, ident, procs, params, bq):
        super(BaseRunner, self).__init__(
            name="%s[%s]" % (self.__class__.__name__, ident)
        )
        self.id = ident
        self.bq = bq
        self.procs = procs
        self.params = params

    def shell(self, args, **kwargs):
        logfile = kwargs.pop('logfile', 'log.txt')
        description = kwargs.pop('description', 'Shell')
        verbose = kwargs.pop('verbose', False)
        working_dir = kwargs.pop('working_dir', '.')
        kill_children = kwargs.pop('kill_children', True)
        timeout = kwargs.pop('timeout', 0)

        self.log("Running %s" % (" ".join(args)))

        if verbose:
            out = open(os.path.join(working_dir, '/', logfile), 'w')
        else:
            out = open(os.devnull, 'w')

        try:
            start = time.time()
            ret = subprocess.Popen(args, stdout=out, stderr=out)
            while ret.poll() is None:
                now = time.time()
                runtime = now - start
                if timeout > 0 and runtime > timeout:
                    ret.kill()
                    msg = "%s: " % (description)
                    msg += "Exceeded timeout. "
                    msg += "%s killed after %d seconds" % (args[0], timeout)
                    self.warn(msg)
                    break
                time.sleep(0.5)
        except Exception as exc:
            msg = "%s: " % (description)
            msg += "Unhandled python error running %s " % (" ".join(args))
            msg += "check log file.\n\t $ "
            msg += str(exc)
            # raise exceptions.FatalError(msg)
            raise Exception(msg)
        finally:
            out.close()
            if kill_children:
                self.kill_subprocess_children(ret.pid)

        if ret.returncode != 0:
            msg = "%s: " % (description)
            msg += "%s returned an error. " % (args[0])
            msg += "check log file.\n\t $ "
            msg += " ".join(args)
            # raise exceptions.FatalError(msg)
            raise Exception(msg)

    def kill_subprocess_children(self, pid):
        """
            Kill any remaining child processes that are left over from a shell
            command.

            Based on code from:
            http://stackoverflow.com/questions/1191374/subprocess-with-timeout
            http://stackoverflow.com/questions/6553423/multiple-subprocesses-with-timeouts

            :param pid: the pid of the parent process
        """
        p = subprocess.Popen(
            "ps --no-headers -o pid --ppid %d" % (pid),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        pids = [pid]
        pids.extend([int(q) for q in stdout.split()])
        try:
            for pid in pids:
                os.kill(pid, signal.SIGKILL)
        except OSError:
            #print "--->OSERROR"
            pass

    def submit(self, runner, **kwargs):
        return self.bq.submit(runner, **kwargs)

    def resubmit(self, **kwargs):
        self.bq.resubmit(self.id, **kwargs)

    def setup(self):
        pass

    def teardown(self):
        pass

    def run(self):
        self.setup()
        retval = self.execute()
        self.teardown()
        sys.exit(retval)

    def execute(self):
        pass

    def log(self, msg):
        print "%s: %s" % (self.name, msg)


class TestRunner(BaseRunner):
    def execute(self):
        self.log("Starting run with %d processors" % (self.procs))
        self.log("Received a value of %d" % (self.params['value']))
        self.log("Sleeping for %d seconds" % (self.params['sleep']))

        args = [
            ['sleep', str(self.params['sleep'])],
            ["dd", "of=/dev/null", "if=/dev/urandom", "bs=1k", "count=10000"],
            ['sleep', str(self.params['sleep'])],
            ["dd", "of=/dev/null", "if=/dev/urandom", "bs=1k", "count=20000"],
            ['sleep', str(self.params['sleep'])],
            ["dd", "of=/dev/null", "if=/dev/urandom", "bs=1k", "count=30000"],
            ['sleep', str(self.params['sleep'])],
            ["dd", "of=/dev/null", "if=/dev/urandom", "bs=1k", "count=40000"],
            ['sleep', str(self.params['sleep'])],
            ["dd", "of=/dev/null", "if=/dev/urandom", "bs=1k", "count=50000"],
            ['sleep', str(self.params['sleep'])]]
        self.shell(args[self.params['value']])

        if self.params['value'] == 9:
            self.resubmit(params={'sleep': 1, 'value': 1})
            self.log("Resubmitting job")
        elif self.params['value'] == 10:
            newparams = {
                'value': randint(1, 10),
                'sleep': randint(1, 10)}
            job = self.submit(
                TestRunner,
                procs=randint(1, 4),
                params=newparams)
            self.log("Submitting new job %s" % (job.ident))

        self.log("Done")


if __name__ == '__main__':
    bq = BatchQueues()
    b = Batch(bq)

    ids = []
    for i in xrange(50):
        procs = randint(1, 4)
        params = {
            'value': randint(1, 10),
            'sleep': randint(1, 10)
        }
        if randint(1, 10) > 9 and i > 3:
            d = -(randint(1, 3))
            id = bq.submit(
                TestRunner,
                procs=procs,
                deps=ids[d:],
                params=params)
        else:
            job = bq.submit(
                TestRunner,
                procs=procs,
                params=params)
        ids.append(job.ident)

    b.run()
