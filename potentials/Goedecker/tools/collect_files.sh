#!/bin/ksh
# Collects the pseudopotential files XX (CPMD format) and QS (CP2K format)
# in the Goedecker/build tree
#
  cp2klibpath=../cp2k
  cpmdlibpath=../cpmd
  cd ../build
  echo Collecting new XX and QS files ...
  for xcfun in blyp bp hcth120 hcth407 pade pbe; do
    for xxfile in $(find $xcfun -name XX); do
      el=$(echo $xxfile | cut -d"/" -f2)
      typeset -L q=$(head -3 $xxfile | tail -1 | cut -d"=" -f2 | cut -d"." -f1)
      qsfile=$(dirname $xxfile)/QS
      if [[ -n $(diff $xxfile $cpmdlibpath/$xcfun/$el-q$q) ]]; then
        cp $xxfile $cpmdlibpath/$xcfun/$el-q$q
        echo "$xxfile was copied to $cpmdlibpath/$xcfun/$el-q$q"
        if [[ -f $qsfile ]]; then
          cp $qsfile $cp2klibpath/$xcfun/$el-q$q
          echo "$qsfile was copied to $cp2klibpath/$xcfun/$el-q$q"
        fi
      fi
    done
  done
