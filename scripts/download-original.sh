#!/bin/bash

set -e
set -o pipefail

# Download only original, non-subsetted files.
URL="https://d3b-openaccess-us-east-1-prd-pbta.s3.amazonaws.com/splicing-neoepitopes"
RELEASE="v4"

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# If RUN_LOCAL is used, the time-intensive steps are skipped because they cannot
# be run on a local computer -- the idea is that setting RUN_LOCAL=1 will allow for
# local testing running/testing of all other steps
RUN_LOCAL=${RUN_LOCAL:-0}

# Get base directory of project
cd ..
BASEDIR="$(pwd)"
cd -

# check if the release folder exists, if not, create a release folder
[ ! -d "$BASEDIR/data/$RELEASE/" ] && mkdir $BASEDIR/data/$RELEASE/

# The md5sum file provides our single point of truth for which files are in a release.
curl -k --create-dirs $URL/$RELEASE/original-md5sum.txt -o data/$RELEASE/original-md5sum.txt -z data/$RELEASE/original-md5sum.txt

# Consider the filenames in the md5sum file and the release notes
FILES=(`tr -s ' ' < $BASEDIR/data/$RELEASE/original-md5sum.txt | cut -d ' ' -f 2` release-notes.md)

for file in "${FILES[@]}"
do
  if [ ! -e "$BASEDIR/data/$RELEASE/$file" ]
  then
    echo "Downloading $file"
    curl --create-dirs $URL/$RELEASE/$file -o $BASEDIR/data/$RELEASE/$file
  fi
done

#check md5sum
cd $BASEDIR/data/$RELEASE
echo "Checking MD5 hashes..."
md5sum -c original-md5sum.txt
cd $BASEDIR

# Make symlinks in data/ to the files in the just downloaded release folder.
for file in "${FILES[@]}"
do
  ln -sfn $RELEASE/$file data/$file
done


