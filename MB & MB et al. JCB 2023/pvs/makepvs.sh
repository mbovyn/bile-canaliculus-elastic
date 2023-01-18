#!/bin/bash

set -e

lvals=(0 .1 1 2 3 4 5 6 7 8 9 10)
vcritvals=(3.975 4.0214669352014 4.78391370538611 5.11913735542822 5.30454512642403 5.46351480111755 5.59990843175743 5.69616687296858 5.79242531417974 5.8886837553909 5.98494219660205 6.08120063781321)
fol_here=$PWD
outfile_fol=p-v
outfile_base=v

mkdir -p "$outfile_fol"

i=0

for v in "${lvals[@]}"; do

  outfile="$outfile_fol/$outfile_base-$i.fe"

  ./../assemble.sh "$outfile"
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
