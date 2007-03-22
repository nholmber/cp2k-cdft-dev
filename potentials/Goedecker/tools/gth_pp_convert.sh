#!/bin/sh
# Convert the psp.par files in the whole potentials/Goedecker/build tree
  exe=$(pwd)/gth_pp_convert.x
  echo Scanning build/ tree ... this may take some time, please wait
  for xx in $(find ../build -name XX); do
    cd $(dirname $xx) >/dev/null
    echo New directory: $(pwd)
    $exe XX ../atom/atom.dat psp.par
#   cat INFO
    cd - >/dev/null
  done
