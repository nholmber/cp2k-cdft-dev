#!/bin/sh
#
  texfile=$(pwd)/../tex/GTH_POTENTIALS.tex
  cat <<*** >${texfile}
\documentclass[10pt,a4paper,twoside]{article}
\usepackage{supertabular}
\setlength{\paperheight}{297 true mm}
\setlength{\paperwidth}{210 true mm}
\setlength{\textheight}{230 true mm}
\setlength{\textwidth}{155 true mm}
\setlength{\evensidemargin}{-0.4 true mm}
\setlength{\oddsidemargin}{4.6 true mm}
\setlength{\topmargin}{-5.4 true mm}
\setlength{\overfullrule}{2pt}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0pt}
\begin{document}
***
  for xcfun in blyp bp hcth120 hcth407 pade pbe pbesol olyp; do
    if [[ -d $(pwd)/../tex/${xcfun} ]]; then
      cd $(pwd)/../tex/${xcfun}
    else
      echo "ERROR: Directory $(pwd)/../tex/${xcfun} not found"
      exit 1
    fi
    XCFUN=$(echo ${xcfun} | tr [:lower:] [:upper:])
    cat <<*** >>${texfile}
$XCFUN:

\vspace{5mm}

\begin{supertabular}{lr@{\hspace{1cm}}r@{\hspace{1cm}}rrrr}
***
    for e in H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar\
             K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se\
             Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn\
             Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy\
             Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb\
             Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf\
             Es Fm Md No Lr; do
      for f in $(ls ${e}-q* 2>/dev/null); do
        cat ${f} >>${texfile}
      done
    done
    cat <<*** >>${texfile}
\end{supertabular}
\newpage
***
    cd ${OLDPWD}
  done
  cat <<*** >>${texfile}
\end{document}
***
  echo New LaTeX file created
