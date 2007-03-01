#!/bin/sh
# Collects the pseudopotential files XX (CPMD format),
# CP2K, and ABINIT in the Goedecker/build tree
#
  typeset -i line1=0 line2=0
  cp2klibpath=../cp2k
  cpmdlibpath=../cpmd
  abinitlibpath=../abinit
  texlibpath=../tex
  cd ../build
  echo Collecting new XX, CP2K, and ABINIT files ...
  for xcfun in blyp bp hcth120 hcth407 pade pbe olyp; do
    for xxfile in $(find $xcfun -name XX); do
      el=$(echo $xxfile | cut -d"/" -f2)
      q=$(head -3 $xxfile | tail -n 1 | cut -d"=" -f2 | cut -d"." -f1)
      q=$(echo $q)
      cp2kfile=$(dirname $xxfile)/CP2K
      cpmdfile=$(dirname $xxfile)/CPMD
      abinitfile=$(dirname $xxfile)/ABINIT
      texfile=$(dirname $xxfile)/TEXTAB
      cpmdlibfile=$(echo ${cpmdlibpath}/${xcfun}/${el}-q${q})
      cp2klibfile=$(echo ${cp2klibpath}/${xcfun}/${el}-q${q})
      abinitlibfile=$(echo ${abinitlibpath}/${xcfun}/${el}-q${q})
      texlibfile=$(echo ${texlibpath}/${xcfun}/${el}-q${q})
      line1=$(grep -n "&POTENTIAL" $xxfile | cut -f1 -d":")
      line2=$(wc -l $xxfile | cut -f1 -d" ")
      head -7 $xxfile >$cpmdfile
      cat $(dirname $xxfile)/INFO >>$cpmdfile
      tail -n $((line2 - line1 + 2)) $xxfile >>$cpmdfile
      if [[ -f $cpmdlibfile ]]; then
        if [[ -n $(diff $cpmdfile $cpmdlibfile) ]]; then
          mv $cpmdfile $cpmdlibfile
          echo "Changed file $cpmdfile was moved to $cpmdlibfile"
        fi
      else
        mv $cpmdfile $cpmdlibfile
        echo "New file $cpmdfile was moved to $cpmdlibfile"
      fi
      if [[ -f $cp2klibfile ]]; then
        if [[ -n $(diff $cp2kfile $cp2klibfile) ]]; then
          mv $cp2kfile $cp2klibfile
          echo "Changed file $cp2kfile was moved to $cp2klibfile"
        fi
      else
        mv $cp2kfile $cp2klibfile
        echo "New file $cp2kfile was moved to $cp2klibfile"
      fi
      if [[ -f $abinitlibfile ]]; then
        if [[ -n $(diff $abinitfile $abinitlibfile) ]]; then
          mv $abinitfile $abinitlibfile
          echo "Changed file $abinitfile was moved to $abinitlibfile"
        fi
      else
        mv $abinitfile $abinitlibfile
        echo "New file $abinitfile was moved to $abinitlibfile"
      fi
      if [[ -f $texlibfile ]]; then
        if [[ -n $(diff $texfile $texlibfile) ]]; then
          mv $texfile $texlibfile
          echo "Changed file $texfile was moved to $texlibfile"
        fi
      else
        mv $texfile $texlibfile
        echo "New file $texfile was moved to $texlibfile"
      fi
    done
  done
