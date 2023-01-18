#!/bin/bash

#create and run .fe files which sweep over pressure in the BC

set -e

lvals=(0 1.212928)
vcritvals=(3.975 4.85096)
fol_here=$PWD
outfile_fol=p-v
outfile_base=v

path_to_code="../../elastic_BC_code"

mkdir -p "$outfile_fol"

i=0

for v in "${lvals[@]}"; do

  outfile="$outfile_fol/$outfile_base-$i.fe"

  "$path_to_code"/scripts/assemble.sh "$outfile"
  echo "


//Commands for junction sweep

v_deformed:=${vcritvals[$i]}
setup
gotorest
fname:=\"$outfile_fol/$outfile_base-$i.txt\"
v_exp:=(v_deformed-rest_vol)/4
v_deformed:=v_deformed-(v_deformed-rest_vol)/2
bulkheadyoungs:=$v
gotodeformed
getmax
get_sums
print_to_file
expand

q 0
" >> "$outfile"

  #open each instance in a screen, don't attach
  screen -d -m -S v$i evolver "$outfile"

  i=$((i+1))

done

#screen -ls to see all running screens
#they should exit when the program completes
#e.g. screen -r c1 to attach to that screen
#pkill SCREEN to kill all
#ctrl+a ctrl+d to detach
