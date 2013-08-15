#!/bin/bash

BRANCH="stable"
if [ $# -eq 1 ] ; then
  BRANCH=$1
fi
BLATZIPBALL="http://users.soe.ucsc.edu/~kent/src/blatSrc${BLATVER}.zip"
BOWTIE2ZIPBALL="http://sourceforge.net/projects/bowtie-bio/files/latest/download?source=files"
if [ "$(uname -s)" == "Darwin" ] ; then
  SPADESTARBALL="http://spades.bioinf.spbau.ru/release2.5.0/SPAdes-2.5.0-Darwin.tar.gz"
else
  SPADESTARBALL="http://spades.bioinf.spbau.ru/release2.5.0/SPAdes-2.5.0-Linux.tar.gz"
fi

# Is pip installed
command -v pip > /dev/null 2>&1 || { 
  echo >&2 "Please install pip before installing.";
  exit 1
}

# Is virtualenv installed
command -v virtualenv > /dev/null 2>&1 || { 
  echo >&2 "Please install virtualenv before installing.";
  exit 1
}

# Is git installed
command -v git > /dev/null 2>&1 || { 
  echo >&2 "Please install git before installing.";
  exit 1
}

# Setup a new virtualenv
virtualenv arc
source arc/bin/activate

echo "Installing required python modules (numpy)"
python -c "import numpy" > /dev/null 2>&1 || {
  pip install numpy
}

echo "Installing required python modules (Biopython)"
python -c "import Bio" > /dev/null 2>&1 || {
  pip install Biopython
}

echo "Cloning arc from the $BRANCH branch."
pushd arc
mkdir share
mkdir src

pushd src # arc/src
if [ ! -f ./ARC/README.md ] ; then
  git clone https://github.com/ibest/ARC.git
else
  echo "ARC has already been cloned in src."
fi

pushd ARC # arc/src/ARC
CURR_BRANCH=$(git branch | sed -n '/\* /s///p')
if [ "$CURR_BRANCH" != "$BRANCH" ] ; then
  echo "Branch differs.  Checking out the requested branch."
  git checkout $BRANCH
  # Reinstall since we changed branches
  echo "Reinstalling."
  python setup.py install
else
  command -v ARC > /dev/null 2>&1 || {
    echo "Installing."
    python setup.py install
  }
fi
popd # arc/src

export MACHTYPE=$(uname -m)
# Grab the path to our virtualenv bin folder
export BINDIR=$(echo $PATH | cut -d':' -f1)

echo "Installing Blat"
command -v blat > /dev/null 2>&1 || {
  if [ ! -f src/blat.zip ] ; then
    wget $BLATZIPBALL -O blat.zip
    unzip blat.zip
    pushd blatSrc # arc/src/blatSrc
    make
    popd # arc/src
  fi
}

echo "Installing Bowtie2"
command -v bowtie2 > /dev/null 2>&1 || {
  if [ ! -f src/bowtie2.zip ] ; then
    wget $BOWTIE2ZIPBALL -O bowtie2.zip
    unzip bowtie2.zip
    pushd bowtie2-* # arc/src/bowtie2-*
    make
    install -c bowtie2 $BINDIR
    for file in bowtie2-*[a-z] ; do 
      install -c $file $BINDIR
    done
    popd # arc/src
  fi
}

echo "Installing Spades"
command -v spades > /dev/null 2>&1 || {
  if [ ! -f src/spades.zip ] ; then
    wget $SPADESTARBALL -O SPAdes.tar.gz
    tar zxvf SPAdes.tar.gz
    pushd SPAdes-* # arc/src/SPAdes-*
    install -c bin/bwa-spades $BINDIR
    install -c bin/spades $BINDIR
    install -c bin/hammer $BINDIR
    install -c bin/spades.py $BINDIR
    install -c bin/spades_init.py $BINDIR
    mv share/spades ../../share
    popd # arc/src
  fi
}

popd # arc

if [ ! -d data ] ; then
  cp -rf src/ARC/test_data data
fi



