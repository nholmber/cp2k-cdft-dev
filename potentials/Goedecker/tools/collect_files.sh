#!/bin/sh
# Collects the pseudopotential files XX (CPMD format) and QS (CP2K format)
# in the Goedecker/build tree
#
  cp2klibpath=../cp2k
  cpmdlibpath=../cpmd
  cd ../build
  echo Collecting new XX and QS files ...
  for xcfun in blyp bp hcth120 hcth407 pade pbe olyp; do
    for xxfile in $(find $xcfun -name XX); do
      el=$(echo $xxfile | cut -d"/" -f2)
      q=$(head -3 $xxfile | tail -1 | cut -d"=" -f2 | cut -d"." -f1)
      q=$(echo $q)
      qsfile=$(dirname $xxfile)/QS
      cpmdlibfile=$(echo ${cpmdlibpath}/${xcfun}/${el}-q${q})
      cp2klibfile=$(echo ${cp2klibpath}/${xcfun}/${el}-q${q})
      if [[ -f $cpmdlibfile ]]; then
        if [[ -n $(diff $xxfile $cpmdlibfile) ]]; then
          cp $xxfile $cpmdlibfile
          echo "Changed file $xxfile was moved to $cpmdlibfile"
          mv $qsfile $cp2klibfile
          echo "Changed file $qsfile was moved to $cp2klibfile"
        fi
      else
        cp $xxfile $cpmdlibfile
        echo "New file $xxfile was moved to $cpmdlibfile"
        mv $qsfile $cp2klibfile
        echo "New file $qsfile was moved to $cp2klibfile"
      fi
    done
  done
