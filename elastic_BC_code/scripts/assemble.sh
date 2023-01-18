#!/bin/bash

# Assemble the .fe file to be called by surface evolver from the geometry file
# (Created in Mathematica, build_bc_param.nb) and file containing functions
# (commands.txt)
# This script is intended to be called by a launching script, prototype:
# evolve.sh

# first arguement should be the name of the file to be written

set -e

#path of this script
MY_PATH=$(dirname "$0")

# insert the preamble at the pointer text
sed -e "/preamblehere/r $MY_PATH/../se_code/preamble.txt" "$MY_PATH/../geometry/bc_initial.geom" > "$1"

# add the commands to the end of the geom file
cat "$MY_PATH/../se_code/commands.txt" >> "$1"
