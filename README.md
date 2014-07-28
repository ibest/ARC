# ARC (Assembly by Reduced Complexity)

### For details about ARC installation, usage, and output, please visit the [ARC web page: http://ibest.github.io/ARC/](http://ibest.github.io/ARC/).


## News

### ARC Updates:

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
