# ARC (Assembly by Reduced Complexity)

ARC is a pipeline which facilitates iterative, reference based assemblies with the intent of reducing bias in the resulting contigs as compared to a purely mapping based approach. The software is designed to work in situations where a whole-genome assembly is not the objective, or fails to produce good results. ARC decomplexifies the traditionally difficult problem of assembly by subsetting the reads into small, manageable pieces which can then be assembled quickly and efficiently. Applications include any type of targeted sequencing approach in which a set of targets is available to map against.

ARC has shown promise in:

* Assembly of bacterial plasmids
* Assembly of viral genomes
* Assembly of mitochondrial genomes using a distantly related reference
* Assembly of exome capture data
* Assembly of chloroplast genomes

ARC is designed to:

* Break large, complex problems into smaller manageable chunks
* Reduce memory footprint requirements (many assemblies should work on a desktop PC)
* Be highly scalable, running multiple jobs simultaneously in parallel
* Be easy to use, portable and simple to configure

The algorithm in a nutshell:

* map reads against a set of targets using BLAT or Bowtie2
* extract mapped reads
* assemble mapped reads into contigs using Roche/Newbler or Spades assemblers
* map reads against the newly formed contigs
* iterate until stopping conditions have been met


## Installation
### Prerequisites:

A mapper and assembler must be installed in your path, the following are currently supported.

* Mapper (either will work well):
    * Bowtie 2 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
    * Blat - http://genome.ucsc.edu/FAQ/FAQblat.html#blat3
        * Note that Blat currently only supports FASTA format as an input for the reads. We have contacted the author with a patch to add FASTQ support and hope that it will be incorporated in the release version soon. Until then we recommend using Bowtie 2 for aligning FASTQ reads.
* Assembler:
    * Roche/Newbler assembler - http://www.454.com/products/analysis-software/
    * Spades assembler - http://bioinf.spbau.ru/spades/

### Additional Requirements:
* Python 2.7.X - http://www.python.org/getit/ (not tested, but should work with 3.3X)
* Python Modules:
    * BioPython - http://biopython.org/wiki/Download (tested with v 1.6.0, but should work with any recent release)

Any combination of the mappers and assemblers above should work well, however Newbler and Bowtie 2 are typically faster in our tests.

### Get the source:
    $ git clone git://github.com/ibest/ARC.git

### Option 1: Install to the system path:

    $ python setup.py install

### Option 2: Install using virtualenv:
Move to the directory where you keep all of your python virtual environments and run the following commands:

    $ cd ~/pyenvs
    $ virtualenv arc
    $ source arc/bin/activate
    $ cd /path/to/arc/source
    $ python setup.py install

### Option 3: Run from the git-clone folder without installing:
    
    $ ./ARC/bin/ARC

### Testing your install:
    
    A small test assembly can be run using the included test data and the following commands:

    $ cd test_data
    $ ./runarc

    This will assemble PhiX reads into contigs for a pair of targets.

## Usage

* Combine reads into a maximum of 3 files per sample: PE1, PE2, SE (where PE stands for Paired End, and SE for Single End). These files can be fasta or fastq formatted.
* Ensure that the targets file is in fasta format and that all entries have unique names.
* Create a file named ARC_config.txt (see the files in test_data for an example). Put this file in a working directory on a drive with plenty of free space.
* Run ARC using the approach appropriate to the installation method selected.

### Outputs
ARC will create a set of folders corresponding to the samples included in the ARC_config.txt. These will be labeled working_* and finished_*.

The working_* folders contain the following files:
* In_contigs.fasta - these are the intermediate contig results for each iteration
* mapping_log.txt - the output produced by the mapper
* \*.idx files - Biopython indexes of the reads (these are kept to facilitate restarting without re-indexing the read files)

The finished_* folders contain the following files:
* contigs.fasta - this contains all contigs for the sample
    * Contig names have the following format: SampleID_:_TargetID_:_ContigN where there can be multiple contigs for one target
* PE1.fastq, PE2.fastq, SE.fastq - these files contain the reads which were mapped on the final iteration of ARC.
    * reads have a slightly modified name to ensure compatibility with Newbler and contain a Sample_:_Target field in the description.


## Configuring ARC
All configuration for ARC is stored in the ARC_config.txt file. Consult the example configuration file in test_data for exact formatting requirements. The following options are currently supported:
* reference - path to the reference fasta which contains targets
* numcycles - the maximum number of iterations ARC is allowed to run
* mapper - which mapper to use (blat/bowtie2)
* assembler - which assembler to use (newbler/spades)
* nprocs - the number of processors to use
* format - the read format (fasta/fastq)
* urt - only applies to Newbler, causes it to 'use read tips' for all but the final assembly
* map_against_reads - causes ARC to map against reads instead of contigs on the second iteration



## Tips and tricks for advanced users
* Starting with multiple sequences for a single Target:
    * Name the sequences so that they share a common element in the reference fasta file in the following way:
        * ID_:_Target1_:_seq1
        * ID_:_Target1_:_seq2
* Setting map_against_reads:
    * Only a small number of reads may map on the first iteration if your reference is very distantly related or coverage is low. In these situations, set map_against_reads 

[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/ccea24d058d3315f3610784acc00af67 "githalytics.com")](http://githalytics.com/ibest/ARC)