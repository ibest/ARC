# ARC (Assembly by Reduced Complexity)

TODO: Write description here
ARC is a pipeline which facilitates iterative, reference based assemblies with the intent of reducing bias in the resulting contigs as compared to a purely mapping based approach. ARC 

## Installation
### Get the source:
    $ git clone git://github.com/ibest/ARC.git

### Install using virtualenv:
Move to your directory that you keep all of your python virtual environments and run the following commands:

    $ cd ~/pyenvs
    $ virtualenv arc
    $ source arc/bin/activate
    $ cd /path/to/arc/source
    $ python setup.py install

## Usage

Pre-requisites:

    *Install either bowite2 or blat for mapping. Make sure that the executable for these is in your path.
    *Combine reads into a maximum of 3 files per sample: PE1, PE2, SE (where PE stands for Paired End, and SE for Single End). These files can be fasta or fastq formatted.
    *Ensure that the targets file is in fasta format and that all entries have unique names.

1) Create an ARC_config_fastq.txt file using the one in "test_data" as a template. Put this file in a working directory on a drive with plenty of free space.
2) run ARC: /path_to_arc/bin/ARC


## Running unit tests:

    TODO: Write usage instructions here

## Contributing

1. Fork it
2. Create your feature branch (`git checkout -b my-new-feature`)
3. Commit your changes (`git commit -am 'Add some feature'`)
4. Push to the branch (`git push origin my-new-feature`)
5. Create new Pull Request


[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/1c58b31862cc6c8fdb36f4ccf190dd52 "githalytics.com")](http://githalytics.com/ibest/ARC)

[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/23406c26a8e65285a4b6e64bec82ec27 "githalytics.com")](http://githalytics.com/ibest/ARC)