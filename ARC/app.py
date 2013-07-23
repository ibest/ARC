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
import logging
from ARC.run import Run
from ARC import spawn
from ARC import logger
from ARC import exceptions

class App:
    def start(self):
        try:
            logger.setup(loglevel=logging.INFO)

            logger.info("Reading config file...")
            run = Run('ARC_config.txt')
            config = run.config()

            logger.info("Setting up working directories and building indexes...")
            run.setup()

            logger.info("Setting up multiprocessing...")
            spawn.run(config)
            self.clean()

            return 0
        except exceptions.FatalError as e:
            logger.error("A fatal error was encountered. \n\t%s" % str(e))
            return 1
        except (KeyboardInterrupt, SystemExit):
            logger.error("%s unexpectedly terminated" % (__name__))
            self.clean()
            return 1

    def clean():
        pass
