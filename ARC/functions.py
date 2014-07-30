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

#Functions for repeat masking:
def num_unmers(seq, N):
    #Calculate the number of unique nmers in seq
    nmers = {}
    for i in range(len(seq) - (N - 1)):
        nmers[str(seq[i:i+N]).upper()] = True
    return(len(nmers))


def mask_seq(seq, mapper, W=15, N=3):
    #Replace simple repeats with 'n' characters
    #This masks a window if the number of unique Nmers is <
    seq_copy = bytearray(seq)
    i = 0
    for i in range(len(seq) - (W - 1)):
        if num_unmers(seq[i:i + W], N) < 7:
            if mapper == 'blat':
                seq_copy[i:i + W] = seq[i:i + W].lower()
            if mapper == 'bowtie2':
                seq_copy[i:i + W] = 'n' * W
    return(seq_copy)


# def setTargetStatus(target, targetLength, iteration, contigs, contig_length, params):
#     '''Call this each time a target is finished in order to track the status'''
#     params['summary_stats'][target]['status'] = status
#     params['summary_stats'][target]['Iteration'] = params['iteration']
#     params['summary_stats'][target]['Reads'] = params['readcounts'][target][params['iteration']]
#     params['summary_stats'][target]['Contigs'] = num_contigs
#     params['summary_stats'][target]['ContigLength'] = contig_length


def writeTargetStats(finished_dir, sample, target, targetLength, status, iteration, readcount, num_contigs, contig_length):
    '''Call this when a target is finished/killed/NoContigs/Repeat, write all outputs to a file'''
    #write out statistics:
    #if num_contigs == 0:
    #    status = "NoContigs"

    tstf = open(os.path.join(finished_dir, "target_summary_table.tsv"), 'a')
    tstf.write('\t'.join(
               [sample,
                target,
                str(targetLength),
                status,
                str(iteration),
                str(readcount),
                str(num_contigs),
                str(contig_length)]) + '\n')
    tstf.close()

