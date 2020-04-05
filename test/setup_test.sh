#!/bin/bash

# The goal of this script is to download data files hosted an a github repo.
# The script has to work both on Linux and MacOS so we could not use some
# nice Linux-specific features of tar and instead had to manually extract
# necessary files via a temp folder.
set -ex

PROJECT_DIR=$(pwd)
mkdir -p ${PROJECT_DIR}

TMPDATADIR=`mktemp -d 2>/dev/null || mktemp -d -t 'TMPDATADIR'`

# check if tmp dir was created
if [[ ! "$TMPDATADIR" || ! -d "$TMPDATADIR" ]]; then
    echo "Could not create temp dir"
    exit 1
fi

# deletes the temp directory at exit
function cleanup {
    rm -rf "$TMPDATADIR"
    echo "Deleted temp directory $TMPDATADIR"
}

trap cleanup EXIT


cd ${TMPDATADIR}
#curl -LkSs https://api.github.com/repos/mirnylab/distiller-test-data/tarball | tar -zxf -
wget -O - https://api.github.com/repos/mirnylab/distiller-test-data/tarball | tar xvz
# cd to the first (and the only) folder that was extracted from the tarball
cd $(ls -d */|head -n 1)
mv -n ./genome ${PROJECT_DIR}
mv -n ./fastq ${PROJECT_DIR}

cd ${PROJECT_DIR}
