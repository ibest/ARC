## News

ARC Updates:
2013-10-21:
 * After a month with no problems on the code-refactor develop is being merged back into stable. Please report any issues.
 * Two new features will be implemented soon. These include
   1. If an assembly is killed after iteration 1, contigs assembled in the previous iteration will be copied to the finished folder from the previous iteration.
   2. A new parameter will be added allowing for sub-sampling of reads. This can help tremendously in situations of very high coverage when using the Newbler assembler.
 * Support for subsampling reads has been added. For example, add "subsample=0.8" to include only 80% of mapped reads in assemblies at each iteration.

2013-09-18:
 * A major code refactor was recently pushed to the "develop" branch. You can test out this branch by typing "git checkout develop" after cloning the ARC repository. Please help us by submitting any bugs you find. A big thanks to Rob Lyon for his hard work re-organizing and modifying the code for better readability/maintainability.
 * The develop branch now comes with a patch which adds Fastq support to the BLAT mapper. If you are not familiar with the Unix/Linux "patch" tool, please see the updated installation instructions for information on how to apply this patch.
 * Basic support for cDNA assemblies has been integrated into ARC. Enable cDNA assemblies using cdna=True in the ARC_config.txt. Currently this only works with the Roche/Newbler assembler. After the run, check the finished_* folders for "isogroup_read_counts.tsv" files. These numbers represent the numbers of reads incorporated into assembled transcripts within the isogroup, and can be used to calculate an RPKM-like measure of gene expression. Preliminary results using reads generated from ovine tissue samples show good clustering of treatment and tissue based on these methods.
 * We plan to release a tool-kit for ARC-based cDNA assembly, annotation, and expression analysis in the near future.
 * ARC was presented in a talk by Dr. Matt Settles at the 5th International Symposium on Animal Functional Genomics in Brazil.

2013-08-30:
Updates on ARC progress:
 * ARC just finished 1.3 million assemblies in ~80hrs on a 60 core server! This was accomplished using a dataset containing 52 samples and ~6500 targets.
 * Tackling this big dataset with ARC exposed some issues with memory usage and speed which have now been addressed. ARC should be faster than ever.
    * These improvements to memory as well as some improvements to logging etc have been added to the Stable (default) branch.
 * ARC now outputs mapping statistics on each iteration with details per-target. See "mapping_stats.tsv" in the finished_* folders.
 * Some preliminary R scripts for profiling and plotting memory usage have been added in the Extensions folder.
    * Run profilemem.R during your ARC run to collect data (this collections information on ALL processes for a user, so it can be misleading if you use your server for something else while running ARC).
    * After profilemem.R has collected some data, run plot_memprofile.R (or source it from within R) to see some plots.
 * Work is currently underway to add preliminary RNA-seq support to ARC, stay tuned.


2013-07-09:
New features added and some modifications to output:
 * A an improvement to the way ARC handles splitting read which results in a major speedup.
 * You can now set an assemblytimeout in ARC_config.txt. ARC will monitor assemblies, and if they run longer than assemblytimeout minutes they will be killed. Preliminary testing has shown that this works great for large projects where some of the targets may contain repeats or flanking regions with repeats which can cause the assembler to founder and block other assemblies from running.
  * assemblytimeout defaults to 10 minutes and accepts fractional values if you would like to limit assemblies to less than a minute.
  * No contigs are output for assemblies that are killed for running too long, however all of the reads which mapped to that contig are written into the finished_* reads files for further analysis if necessary. Additionally, previously generated contigs (if any) will remain in the working*/I*_contigs.fasta intermediate contig files.
  *  Previously assemblies which didn't result in any contigs had their reads written to the finished* contigs.fasta. We found that this cluttered up the contigs.fasta file and this behavior has been changed so that reads are now written to the PE1, PE2 and SE fastq/fasta files and no contigs are written.
 * Output has been mostly standardized so that all output starts with "Sample: SampleID target: TargetID" which makes it easy to grep the assembly log and get progress information for a single sample. It is recommended that you run ARC with stdout pipped to a log file to facilitate this: ARC > arc.log
 * Added an Extensions folder to hold ARC extensions, currently there is only one, check_status.R. Please contact us if you are interested in contributing an extension.


# ARC (Assembly by Reduced Complexity)

ARC is a pipeline which facilitates iterative, reference guided *de nono* assemblies with the intent of
1. Reducing time in analysis and increasing accuracy of results by only considering those reads which should assemble together
2. Reducing reference bias as compared to mapping based approaches.
The software is designed to work in situations where a whole-genome assembly is not the objective, but rather when the researcher wishes to assemble discreet 'targets' contained within next-generation shotgun sequence data. ARC decomplexifies the traditionally difficult problem of assembly by breaking the reads into small, manageable subsets which can then be assembled quickly and efficiently. Applications include those in which the researcher wishes to *de novo* assemble specific content and a set of semi-similar reference targets is available to initialize the assembly process.

ARC has shown promise in:

* Assembly of bacterial plasmids
* Assembly of viral genomes
* Assembly of mitochondrial genomes using a distantly related reference
* Assembly of exome capture data
* Assembly of chloroplast genomes
* Transcriptomes (Not tested yet)

ARC is designed to:

* Break large, complex problems into smaller manageable chunks
* Reduce memory footprint requirements (many assemblies should work on a desktop/labtop PC)
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
* For more information on installing these programs see the Wiki: https://github.com/ibest/ARC/wiki/Details-on-installation

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

### Patching BLAT to support FASTQ (Optional):
We have found that the BLAT mapper by Jim Kent is an excellent, fast, and efficient piece of software. Unfortunately it doesn't natively support the FASTQ format. Ilya Zhbannikov (author of SeqyClean https://bitbucket.org/izhbannikov/seqyclean) wrote a patch which adds FASTQ support to BLAT. This patch is distributed in the "contrib" folder (but only in the develop branch for now). Please ensure that you have the proper academic/nonprofit status, or acquire a license from Kent Informatics (http://www.kentinformatics.com/) and then use the following set of commands to download, patch, and install BLAT with FASTQ support.

    $ wget http://users.soe.ucsc.edu/~kent/src/blatSrc.zip
    $ unzip blatSrc.zip
    $ patch -p0 </path/to/ARC/contrib/blat+fastq_support.patch
    $ cd blatSrc
    $ export MACHTYPE=x86_64
    $ mkdir ~/bin
    $ mkdir ~/bin/x86_64
    $ make

    Executables for blat with FASTQ support will now be located in your ~/bin/x86_64 folder. Either add this folder to your path, or move the executable files to a folder already in your path.



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
    * Only a small number of reads may map on the first iteration if your reference is very distantly related or coverage is low. In these situations, set map_against_reads=True in the ARC_config.txt file. This will use all of the mapped reads as the new set of targets instead of using contigs. Please note that this feature is experimental and sometimes leads to incorporation of repetitive elements.

## Questions/Comments/Suggestions
Please send questions, comments and suggestions to Sam Hunter (shunter {at} gmail {dot} com).

[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/ccea24d058d3315f3610784acc00af67 "githalytics.com")](http://githalytics.com/ibest/ARC)
