#!/bin/sh

# exper2pyt -- pyt=PYthonTank
# For all datafiles associated with expers listed on the command
# line, this will generate a matching HDF5 file containing the
# data from the TDT tank assocated with the exper/file.
#
# 1. Use dbfind to generate list of all datafiles assocated with exper
# 2. feed datafile list to tank2hdf5.py one file at a time.
#


if [ "$1" = "" ]; then
  echo "usage: $(basename $0) [-force] exper"
elif [ "$1" = "-force" ]; then
  force="-force"
  shift
else
  force=""
fi
dbfind $1 | xargs tank2hdf5.py $force
