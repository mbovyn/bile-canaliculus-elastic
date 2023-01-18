#!/bin/bash

# Assemble .fe file from base commands and geometry, then run surface evolver
# with a set of modified commands

set -e

fname="bc.fe"

# path to the "elastic_BC_code" folder
path_to_code="../elastic_BC_code"
# assemble .fe file
"$path_to_code"/scripts/assemble.sh $fname

# SE code to run, custom commands in elastic_BC_code/se_code/commands.txt
# write it to the file
echo "showq" >> $fname
echo "setup" >> $fname
echo "gotorest" >> $fname
echo "gotodeformed" >> $fname
#echo "nickw:=(.25+(1.18-.25)/2)/3.6" >> $fname
#echo "setnick" >> $fname
#echo "nick" >> $fname

# run evolver on the file that was created
evolver $fname

#tips
#screen -ls to see all running screens
#they should exit when the program completes
#e.g. screen -r c1 to attach to that screen
#pkill SCREEN to kill all
#ctrl+a ctrl+d to detach
