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
from ARC import Config
from ARC import logger
from ARC import FatalError
from ARC import Spawn


class App:
    def start(self, loglevel, configfile='ARC_config.txt'):
        try:
            logger.setup(loglevel=loglevel)

            logger.info("Reading config file...")
            config = Config(configfile)

            logger.info(
                "Setting up working directories and building indexes...")
            config.setup()

            params = config.params()
            # I really would like to name this 'globals' but...
            universals = config.universals()
            spawn = Spawn(params, universals)

            logger.info("Running ARC.")
            spawn.run()

            logger.info("Cleaning up.")
            self.clean()

            return 0
        except FatalError as e:
            logger.error("A fatal error was encountered. \n\t%s" % str(e))
            return 1
        except (KeyboardInterrupt, SystemExit):
            # if 'batch' in vars():
            #     batch.killall()
            self.clean()
            logger.error("%s unexpectedly terminated" % (__name__))
            return 1

    def clean(self):
        pass
