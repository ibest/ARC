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
import sys
import logging
import multiprocessing


def setup(logfile=None, loglevel=logging.INFO):
    # Set up a multiprocessing logger to hopefully log from all N workers in a
    # safe and simultaneous fashion:
    logger = multiprocessing.get_logger()
    logger.setLevel(loglevel)
    log_handler = logging.StreamHandler(sys.stdout)
    log_handler.setFormatter(logging.Formatter(
        '[%(asctime)s %(levelname)s %(process)s] %(message)s'))
    log_handler.setLevel(loglevel)  # Here's where
    logger.addHandler(log_handler)


def info(message):
    logger = multiprocessing.get_logger()
    logger.info("%s" % (message))


def error(message):
    logger = multiprocessing.get_logger()
    logger.error("%s" % (message))


def debug(message):
    logger = multiprocessing.get_logger()
    logger.debug("%s" % (message))


def warn(message):
    logger = multiprocessing.get_logger()
    logger.warn("%s" % (message))



# From:
# http://stackoverflow.com/questions/641420/how-should-i-log-while-using-
# multiprocessing-in-python/894284#894284
# Modified for StreamHandler (instead of RotatingFileHandler)
#from logging.handlers import RotatingFileHandler
#from logging.handlers import StreamHandler
# import multiprocessing, threading, logging, sys, traceback

# class MultiProcessingLog(logging.Handler):
#     def __init__(self, name, mode, maxsize, rotate):
#         logging.Handler.__init__(self)

#         #self._handler = RotatingFileHandler(name, mode, maxsize, rotate)
#         self._handler = logging.StreamHandler(sys.stdout)
#         self.queue = multiprocessing.Queue(-1)

#         t = threading.Thread(target=self.receive)
#         t.daemon = True
#         t.start()

#     def setFormatter(self, fmt='[%(asctime)s %(levelname)s %(process)s] %(message)s'):
#         logging.Handler.setFormatter(self, fmt)
#         self._handler.setFormatter(fmt)

#     def receive(self):
#         while True:
#             try:
#                 record = self.queue.get()
#                 self._handler.emit(record)
#             except (KeyboardInterrupt, SystemExit):
#                 raise
#             except EOFError:
#                 break
#             except:
#                 traceback.print_exc(file=sys.stderr)

#     def send(self, s):
#         self.queue.put_nowait(s)

#     def _format_record(self, record):
#         # ensure that exc_info and args
#         # have been stringified.  Removes any chance of
#         # unpickleable things inside and possibly reduces
#         # message size sent over the pipe
#         if record.args:
#             record.msg = record.msg % record.args
#             record.args = None
#         if record.exc_info:
#             dummy = self.format(record)
#             record.exc_info = None

#         return record

#     def emit(self, record):
#         try:
#             s = self._format_record(record)
#             self.send(s)
#         except (KeyboardInterrupt, SystemExit):
#             raise
#         except:
#             self.handleError(record)

#     def close(self):
#         self._handler.close()
#         logging.Handler.close(self)
