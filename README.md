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

## Usage

### Prerequisites:

* Either bowite2 or blat must be available in your path. 
    * Note that blat currently only supports FASTA format. We have contacted the author with a patch to add FASTQ support and hope that it will be incorporated soon.
* Combine reads into a maximum of 3 files per sample: PE1, PE2, SE (where PE stands for Paired End, and SE for Single End). These files can be fasta or fastq formatted.
* Ensure that the targets file is in fasta format and that all entries have unique names.
* Create a file named ARC_config.txt (see the files in test_data for an example). Put this file in a working directory on a drive with plenty of free space.
* Run ARC using the approach appropriate to the installation method selected.

### Outputs
ARC will create a set of folders corresponding to the samples included in the ARC_config.txt. These will be labeled working_* and finished_*.


## Testing your install:
    
    A small test assembly can be run using the included test data and the following commands:

    $ cd test_data
    $ ./runarc

    This will assemble PhiX reads into contigs for a pair of targets.

## Configuring ARC
All configuration for ARC is stored in the ARC_config.txt file. 



## Tips and tricks for advanced users
* Starting with multiple sequences for a single Target: ARC

[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/ccea24d058d3315f3610784acc00af67 "githalytics.com")](http://githalytics.com/ibest/ARC)