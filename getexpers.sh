#!/bin/sh

# print list of all expers contained in directories on
# the command line to stdout

find $* -name '*[0-9][0-9][0-9][0-9].*.[0-9][0-9][0-9]' | \
     awk -F/ '{print $NF}' | awk -F. '{print $1}'  | \
     grep -v '0000$' | sort -n | uniq
