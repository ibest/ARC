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
import logging

class MapperRunner:
    def __init__(self, target, sample)
        self.target = target
        self.sample = sample
        self.next = []
        self.error = None

    def start(self):
      print "Running the mapper"

    def SAM_to_dict(self, filename):
        """ Read a SAM file to a mapping dict and return it """
        try:
            inf = open(filename, 'r')
        except Exception as inst:
            if type(inst) == IOError:
                logging.error("Failed to open mapping dictionary %s." % filename)
            raise inst

        read_map = {}  # target:{read} dictionary of dictionaries
        i = 0
        startT = time.time()

        for l in inf:
            i += 1
            if l[0] != "@":
                l2 = l.strip().split()
                if l2[2] == "*":  # skip unmapped
                    continue
                readid = l2[0]
                target = l2[2].split("_")[0].split("=")[1]
                if target not in read_map: read_map[target] = {}
                read_map[target][readid] = 1
        #Report total time:
        logging.info("Processed %s lines in %s seconds." % (i, time.time() - startT ))
        return read_map

    def PSL_to_dict(self, filename):
        try:
            inf = open(filename, 'r')
        except Exception as inst:
            if type(inst) == IOError:
                logging.error("Failed to open mapping dictionary %s." % filename)
            raise inst
        read_map = {}
        i = 0
        startT = time.time()

        psl_header = False

        for l in inf:
            i += 1
            #Check for PSL header and skip 5 lines if it exists
            if i == 1 and l.split()[0] == 'psLayout':
                psl_header = True 
            if psl_header and i <= 5:
                continue
            l2 = l.strip().split("\t")
            readid = l2[9]
            target = l2[13]
            if target not in read_map:
                read_map[target] = {}
            read_map[target][readid] = 1
        logging.info("Processed %s lines in %s seconds." % (i, time.time() - startT))
        return read_map


    def write_dict(self, filename, read_map):
        """ Write a mapping dictionary to a file. """
        startT = time.time()
        outf = open(filename, 'w')
        for k in read_map.keys():
            outf.write(k + '\t' + ",".join(read_map[k].keys()) + '\n' )
        outf.close()
        logging.info("Wrote all values to txt in %s seconds" % (time.time() - startT ))


    def read_dict(self, filename):
        """ Read a mapping dictionary from a file """
        startT = time.time()
        try:
            inf = open(filename, 'r')
        except Exception as inst:
            if type(inst) == IOError:
                logging.error("Failed to open mapping dictionary %s." % filename)
            raise inst
        new_map = {}
        for l in inf:
            l2 = l.split('\t')
            l3 = l2[1].strip().split(",")
            new_map[l2[0]] = {}
            for k in l3:
                new_map[l2[0]][k] = 1
        logging.info("Read all values to txt in %s seconds" % (time.time() - startT ))
        return new_map
