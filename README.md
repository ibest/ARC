# ARC (Assembly by Reduced Complexity)

### For details about ARC installation, usage, and output, please visit the [ARC web page: http://ibest.github.io/ARC/](http://ibest.github.io/ARC/).


## News

### ARC Updates:

#### 2016-03-27:
* Updated version to v1.1.4-beta to reflect a number of modifications made over the last year+ these changes include:
    * Low map qual reads are now recruited. Previously reads with mapq < 10 were not included (Bowtie2 specific).
    * Added a double check to ensure that an exception would be raised if a read was somehow present in PE1 but absent from PE2 (this should only occur if there is an error in the input files, or reads are named in an unexpected way) (issue #50).
    * Added option to make Spades assemblies faster by using "--only-assembler" option in all but the last assembly (issue #42). This parameter is disabled by default, set only-assembler=True to enable. Preliminary testing indicates:
        * ~30% improved speed with slightly more recruited reads
        * it isn't known how this impacts assembly quality or convergence rate
    * Modified bowtie2 mapping to add the --no-unal parameter, resulting in a significant reduction in disk space requirements in some cases (should also reduce I/O overhead). This is enabled by default (issue #44).
    * Added error handling around indexing of fastq files so that if ARC is killed during indexing, it will clean up partial index files before exiting (this seems to be one of the most common sources of error reports) (issue #46).
    * fixed silly issue where config file name was not properly printed when missing reference or reads were detected in config file.
    * fixed silly bug caused by incorrectly joining a list of strings when invalid assembler was detected.
    * modified the way that read IDs are handled, previously both old and new style Illumina reads would work, but other formats were not supported. Now by setting "sra=True", reads which differ in the last character but don't have pair number delimited by "/" will still work (e.g. for reads extracted from the SRA with read IDs ending in *.1 and *.2). (Issue #47)
    * Identified an issue where reads were not being recruited properly. In an earlier modification various parts of the code were modified to enforce assigning the same read to no more than N targets (set by the bowtie2_k parameter). This had the unintended side effect of making it impossible to use ARC contigs as targets for a second ARC run without renaming them first. This modification adds a message in the log file indicating that reads were mapped but discarded (bug #49).
    * Fixed regular expressions for parameter conversion so that non-numeric parameters starting in a numeric digit are detected properly (bug #51).

#### 2014-09-02:
* Released version v1.1.3 to reflect the following tweaks, improvements and bug fixes:
    * ARC no longer produces a "deprecated method" warning with Biopython v1.64
    * ARC now supports high-specificity mapping on iteration 1, set sloppymapping=False
    * ARC now supports doing work in a different folder, set workingdirectory='/path/to/fast/storage'
    * ARC now generates new summary tables named target_summary_table.tsv
    * ARC now supports simple-repeat-masking on iteration 1 when maskrepeats=True is set
    * ARC default max_incorporation now defaults to 10 instead of 5.

#### 2014-08-21:
* Two exciting updates for ARC (in the develop branch):
    * ARC's behavior has been modified so that it will mask the reference sequences for simple repeats on iteration 1 if maskrepeats=True (the default).
    * A new summary table named "target_summary_table.tsv" is now generated which complements the information already reported in "mapping_stats.tsv".
        * This file is created for each sample, and contains target specific details for every target in the original set of reference sequences.
        * Targets can have one of the following statuses:
            * NoReads - no reads were recruited by this target
            * Finished - target was finished normally
            * Killed - Assembly time was longer than assemblytimeout
            * Repeat - target recruited too many reads and may be a repeat
* It was found that the default max_incorporation=5 repeat detection setting was too aggressive in cases where divergent references were used. The default for this setting has been adjusted to a more conservative max_incorporation=10.

#### 2014-07-28:
* ARC (develop branch) now has a "sloppymapping" parameter which can be used to turn off the low-specificity mapping on iteration 1. This works with both Bowtie2 and BLAT.

#### 2014-07-24:
* ARC has a [web page: http://ibest.github.io/ARC/](http://ibest.github.io/ARC/)
* ARC has a [Google Group: https://groups.google.com/forum/?hl=en#!forum/arc-assembly](https://groups.google.com/forum/?hl=en#!forum/arc-assembly)

#### 2014-03-12:
* Some major improvements were made over the last few months, these include:
    * A total refactor of the job handling mechanism which improved efficiency and speed, and simplified the code significantly.
    * A bug was identified in the code which would allow bowtie2 to map a read to more than one target even when bowtie_2k=1 was set (in certain circumstances). This should be resolved now, but more testing wouldn't hurt.
    * Fixed a race condition which sometimes caused ARC to crash on slow file systems.
    * Version number is ALWAYS printed out to help with debugging.
    * The latest version should by v1.1.1. Upgrading to this version is recommended.

#### 2014-01-10:
* Fixing the Git branches so that we have "master" and "develop" again.
* A variety of improvements and bug fixes which all seem to be working well.
* Prepping for a v1.1.0 release.
