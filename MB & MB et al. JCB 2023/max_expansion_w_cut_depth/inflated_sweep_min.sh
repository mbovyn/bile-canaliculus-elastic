#!/bin/bash

set -e

lvals=(.01 .1 .316 1 2.15 4.64 10 21.5 46.4 100)
fol_here=$PWD
outfile_fol="$fol_here/inflated_sweep_min"
outfile_base=v

mkdir -p "$outfile_fol"

i=0

for v in "${lvals[@]}"; do

  outfile="$outfile_fol/$outfile_base-$i.fe"
  #echo "$outfile"
  ./../assemble.sh "$outfile"

  echo "


//Commands for inflated sweep

setup
gotorest
fname:=\"$outfile_fol/$outfile_base-$i.txt\"
bulkheadyoungs:=$v
set_rest
g500
also_ablate:=1
nick_ablate:=1
nickw:=.25/3.6
fname_expand:=\"$outfile_fol/$outfile_base-$i\_expand.txt\"
v_exp:=.5
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
