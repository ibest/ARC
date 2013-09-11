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
from ARC import logger
from ARC import FatalError
from ARC import TimeoutError
from ARC import RerunnableError
from ARC import SubprocessError
# from ARC import Job
import os
import time
import subprocess
import signal
import sys
import logging


class Base:
    # which queues do I need to pass?  also globals?
    def __init__(self, jobid, procs, params, bq):
        self.jobid = jobid
        self.bq = bq
        self.procs = procs
        self.params = params
        self.loglevel = logger.level()
        self.globals = bq.globals

    def delete(self):
        del self.jobid
        del self.bq
        del self.procs
        del self.params

    def run_on_exit_ok(self):
        pass

    def run_on_exit_error(self, retval):
        pass

    def shell(self, args, **kwargs):
        logfile = kwargs.pop('logfile', 'log.txt')
        description = kwargs.pop('description', 'Shell')
        verbose = kwargs.pop('verbose', False)
        working_dir = kwargs.pop('working_dir', '.')
        kill_children = kwargs.pop('kill_children', True)
        timeout = kwargs.pop('timeout', 0)

        # self.log("Running %s in %s" % (" ".join(args), working_dir))

        if verbose:
            path = os.path.join(working_dir, logfile)
            self.debug("Logging to %s" % (path))
            out = open(path, 'w')
        else:
            out = open(os.devnull, 'w')

        try:
            start = time.time()
            ret = subprocess.Popen(args, stdout=out, stderr=out)
            while ret.poll() is None:
                now = time.time()
                runtime = now - start
                if timeout > 0 and runtime >= timeout:
                    ret.kill()
                    msg = "%s: " % (description)
                    msg += "Exceeded timeout. "
                    msg += "%s killed after %d seconds" % (args[0], timeout)
                    raise TimeoutError(msg)
                time.sleep(0.1)
        except OSError as exc:
            msg = "Failed to run. \n\t$ %s\n\t! " % (" ".join(args))
            msg += str(exc)
            raise SubprocessError(msg)
        except Exception as exc:
            msg = "%s: " % (description)
            msg += "Unhandled python error running %s " % (" ".join(args))
            msg += "check log file.\n\t $ "
            msg += str(exc)
            raise Exception(msg)
        finally:
            out.close()
            if 'ret' in vars() and kill_children:
                self.kill_subprocess_children(ret.pid)

        if ret.returncode != 0:
            msg = "%s: " % (description)
            msg += "%s returned an error. " % (args[0])
            msg += "check log file.\n\t $ "
            msg += " ".join(args)
            raise SubprocessError(msg)

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
        self.bq.resubmit(self.jobid, **kwargs)

    def setup(self):
        pass

    def teardown(self):
        pass

    def run(self):
        try:
            self.setup()
            self.execute()
            retval = Job.OK
            self.run_on_exit_ok()
        except TimeoutError as exc:
            retval = Job.TIMEOUTERROR
            self.run_on_exit_error(retval)
            self.warn(exc.msg)
        except FatalError as exc:
            retval = Job.FATALERROR
            self.run_on_exit_error(retval)
            self.error(exc.msg)
        except RerunnableError as exc:
            retval = Job.RERUNERROR
            self.run_on_exit_error(retval)
            self.warn(exc.msg)
        except SubprocessError as exc:
            retval = Job.PROCESSERROR
            self.run_on_exit_error(retval)
            self.warn(exc.msg)
        except Exception as exc:
            retval = Job.UNKNOWNERROR
            self.run_on_exit_error(retval)
            self.exception(exc)
            self.error(exc)
        except (KeyboardInterrupt, SystemExit):
            retval = Job.UNKNOWNERROR
        finally:
            self.teardown()
            self.delete()

        self.debug("Exiting with exitcode %d." % (retval))
        sys.exit(retval)

    def execute(self):
        pass

    def log(self, msg):
        if logger.level() == logging.DEBUG:
            name = self.name
        else:
            name = self.__class__.__name__
        logger.info("%-12s| %s" % (name, msg))

    def info(self, msg):
        if self.loglevel == logging.DEBUG:
            name = self.name
        else:
            name = self.__class__.__name__
        logger.info("%-12s| %s" % (name, msg))

    def debug(self, msg):
        if self.loglevel == logging.DEBUG:
            name = self.name
            logger.debug("%-12s| %s" % (name, msg))

    def warn(self, msg):
        if logger.level() == logging.DEBUG:
            name = self.name
        else:
            name = self.__class__.__name__
        logger.warn("%-12s| %s" % (name, msg))

    def error(self, msg):
        if logger.level() == logging.DEBUG:
            name = self.name
        else:
            name = self.__class__.__name__
        logger.error("%-12s| %s" % (name, msg))

    def exception(self, exc):
        logger.exception(exc)
