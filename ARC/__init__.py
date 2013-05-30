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

from ARC import mapper
from ARC import spawn
from ARC import config

def main():
  #Setup a global logger:
  # How should we handle this gracefully for cases where each component is run independently?
  # How should get get log-level (command line switch), and what level should we default to? 
  logging.basicConfig(filename='run.log', level=logging.INFO)

  #Run modules:
  setup()
  config = read_config()
  run_spawner(config)
  clean()

def setup():
  """Add setup"""

def read_config():
  # config.read()
  return {}

def run_mapper():
  mapper.run()

def run_spawner(config):
  spawn.run(config)

def clean():
  """Clean up"""

