#!/bin/sh

# generate hdf5 files for all tank data asocated with an exper
# - first use dbfind to generate list of all datafiles
# - feed datafile list to tank2hdf5.py

if [ "$1" = "" ]; then
  echo "usage: $(basename $0) [-force] exper"
elif [ "$1" = "-force" ]; then
  force="-force"
  shift
else
  force=""
fi
dbfind $1 | xargs tank2hdf5.py $force
