#!/bin/sh
# Generate the QS files in the whole potentials/Goedecker/build tree
  exe=$(pwd)/cpmd_to_qs.x
  echo Scanning build/ tree ... this may take some time, please wait
  for xx in $(find ../build -name psp.par); do
    cd $(dirname $xx)
    echo New directory: $(pwd)
    $exe XX ../atom/atom.dat QS
#   cat INFO
    cd -
  done
