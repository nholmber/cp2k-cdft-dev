#!/bin/sh
# Purpose: Set the default GTH PPs in a GTH_POTENTIAL file
# Version: 1.0
# History: - Creation (12.08.2016, Matthias Krack)
#          - 

# Check file name argument
if [[ -n $1 ]]; then
   ppfile=$1
   if [[ ! -s ${ppfile} ]]; then
      echo "ERROR: Potential file ${ppfile} not found"
      exit 1
   fi
else
   echo "ERROR: No potential file specified"
   exit 1
fi

# Loop over all default GTH PPs
# Add a new default PP to the list below or change the default PP if needed
for label in \
 "H GTH-BLYP-q1" \
 "He GTH-BLYP-q2" \
 "Li GTH-BLYP-q3" \
 "Be GTH-BLYP-q4" \
 "B GTH-BLYP-q3" \
 "C GTH-BLYP-q4" \
 "N GTH-BLYP-q5" \
 "O GTH-BLYP-q6" \
 "F GTH-BLYP-q7" \
 "Ne GTH-BLYP-q8" \
 "Na GTH-BLYP-q9" \
 "Mg GTH-BLYP-q10" \
 "Al GTH-BLYP-q3" \
 "Si GTH-BLYP-q4" \
 "P GTH-BLYP-q5" \
 "S GTH-BLYP-q6" \
 "Cl GTH-BLYP-q7" \
 "Ar GTH-BLYP-q8" \
 "K GTH-BLYP-q9" \
 "Ca GTH-BLYP-q10" \
 "Sc GTH-BLYP-q11" \
 "Ti GTH-BLYP-q12" \
 "V GTH-BLYP-q13" \
 "Cr GTH-BLYP-q14" \
 "Mn GTH-BLYP-q15" \
 "Fe GTH-BLYP-q16" \
 "Co GTH-BLYP-q17" \
 "Ni GTH-BLYP-q18" \
 "Cu GTH-BLYP-q11" \
 "Zn GTH-BLYP-q12" \
 "Ga GTH-BLYP-q13" \
 "Ge GTH-BLYP-q4" \
 "As GTH-BLYP-q5" \
 "Se GTH-BLYP-q6" \
 "Br GTH-BLYP-q7" \
 "Kr GTH-BLYP-q8" \
 "Sr GTH-BLYP-q10" \
 "Y GTH-BLYP-q11" \
 "Zr GTH-BLYP-q12" \
 "Mo GTH-BLYP-q14" \
 "Ru GTH-BLYP-q16" \
 "Rh GTH-BLYP-q17" \
 "Pd GTH-BLYP-q18" \
 "Ag GTH-BLYP-q11" \
 "In GTH-BLYP-q13" \
 "Sb GTH-BLYP-q5" \
 "Te GTH-BLYP-q6" \
 "I GTH-BLYP-q7" \
 "Ba GTH-BLYP-q10" \
 "Ce GTH-BLYP-q12" \
 "Gd GTH-BLYP-q18" \
 "W GTH-BLYP-q14" \
 "Au GTH-BLYP-q11" \
 "Bi GTH-BLYP-q5" \
 "H GTH-BP-q1" \
 "He GTH-BP-q2" \
 "Li GTH-BP-q3" \
 "Be GTH-BP-q4" \
 "B GTH-BP-q3" \
 "C GTH-BP-q4" \
 "N GTH-BP-q5" \
 "O GTH-BP-q6" \
 "F GTH-BP-q7" \
 "Ne GTH-BP-q8" \
 "Na GTH-BP-q9" \
 "Mg GTH-BP-q10" \
 "Al GTH-BP-q3" \
 "Si GTH-BP-q4" \
 "P GTH-BP-q5" \
 "S GTH-BP-q6" \
 "Cl GTH-BP-q7" \
 "Ar GTH-BP-q8" \
 "K GTH-BP-q9" \
 "Ca GTH-BP-q10" \
 "Sc GTH-BP-q11" \
 "Ti GTH-BP-q12" \
 "V GTH-BP-q13" \
 "Cr GTH-BP-q14" \
 "Mn GTH-BP-q15" \
 "Fe GTH-BP-q16" \
 "Co GTH-BP-q17" \
 "Ni GTH-BP-q18" \
 "Cu GTH-BP-q11" \
 "Zn GTH-BP-q12" \
 "Ga GTH-BP-q13" \
 "Ge GTH-BP-q4" \
 "As GTH-BP-q5" \
 "Se GTH-BP-q6" \
 "Br GTH-BP-q7" \
 "Kr GTH-BP-q8" \
 "Zr GTH-BP-q12" \
 "Ru GTH-BP-q16" \
 "Te GTH-BP-q6" \
 "Cs GTH-BP-q9" \
 "H GTH-HCTH120-q1" \
 "C GTH-HCTH120-q4" \
 "N GTH-HCTH120-q5" \
 "O GTH-HCTH120-q6" \
 "F GTH-HCTH120-q7" \
 "P GTH-HCTH120-q5" \
 "Ar GTH-HCTH120-q8" \
 "H GTH-HCTH407-q1" \
 "C GTH-HCTH407-q4" \
 "N GTH-HCTH407-q5" \
 "O GTH-HCTH407-q6" \
 "H GTH-PADE-q1 GTH-LDA-q1" \
 "He GTH-PADE-q2 GTH-LDA-q2" \
 "Li GTH-PADE-q3 GTH-LDA-q3" \
 "Be GTH-PADE-q4 GTH-LDA-q4" \
 "B GTH-PADE-q3 GTH-LDA-q3" \
 "C GTH-PADE-q4 GTH-LDA-q4" \
 "N GTH-PADE-q5 GTH-LDA-q5" \
 "O GTH-PADE-q6 GTH-LDA-q6" \
 "F GTH-PADE-q7 GTH-LDA-q7" \
 "Ne GTH-PADE-q8 GTH-LDA-q8" \
 "Na GTH-PADE-q9 GTH-LDA-q9" \
 "Mg GTH-PADE-q10 GTH-LDA-q10" \
 "Al GTH-PADE-q3 GTH-LDA-q3" \
 "Si GTH-PADE-q4 GTH-LDA-q4" \
 "P GTH-PADE-q5 GTH-LDA-q5" \
 "S GTH-PADE-q6 GTH-LDA-q6" \
 "Cl GTH-PADE-q7 GTH-LDA-q7" \
 "Ar GTH-PADE-q8 GTH-LDA-q8" \
 "K GTH-PADE-q9 GTH-LDA-q9" \
 "Ca GTH-PADE-q10 GTH-LDA-q10" \
 "Sc GTH-PADE-q11 GTH-LDA-q11" \
 "Ti GTH-PADE-q12 GTH-LDA-q12" \
 "V GTH-PADE-q13 GTH-LDA-q13" \
 "Cr GTH-PADE-q14 GTH-LDA-q14" \
 "Mn GTH-PADE-q15 GTH-LDA-q15" \
 "Fe GTH-PADE-q16 GTH-LDA-q16" \
 "Co GTH-PADE-q17 GTH-LDA-q17" \
 "Ni GTH-PADE-q18 GTH-LDA-q18" \
 "Cu GTH-PADE-q11 GTH-LDA-q11" \
 "Zn GTH-PADE-q12 GTH-LDA-q12" \
 "Ga GTH-PADE-q13 GTH-LDA-q13" \
 "Ge GTH-PADE-q4 GTH-LDA-q4" \
 "As GTH-PADE-q5 GTH-LDA-q5" \
 "Se GTH-PADE-q6 GTH-LDA-q6" \
 "Br GTH-PADE-q7 GTH-LDA-q7" \
 "Kr GTH-PADE-q8 GTH-LDA-q8" \
 "Rb GTH-PADE-q9 GTH-LDA-q9" \
 "Sr GTH-PADE-q10 GTH-LDA-q10" \
 "Y GTH-PADE-q11 GTH-LDA-q11" \
 "Zr GTH-PADE-q12 GTH-LDA-q12" \
 "Nb GTH-PADE-q13 GTH-LDA-q13" \
 "Mo GTH-PADE-q14 GTH-LDA-q14" \
 "Tc GTH-PADE-q15 GTH-LDA-q15" \
 "Ru GTH-PADE-q16 GTH-LDA-q16" \
 "Rh GTH-PADE-q17 GTH-LDA-q17" \
 "Pd GTH-PADE-q18 GTH-LDA-q18" \
 "Ag GTH-PADE-q11 GTH-LDA-q11" \
 "Cd GTH-PADE-q12 GTH-LDA-q12" \
 "In GTH-PADE-q13 GTH-LDA-q13" \
 "Sn GTH-PADE-q4 GTH-LDA-q4" \
 "Sb GTH-PADE-q5 GTH-LDA-q5" \
 "Te GTH-PADE-q6 GTH-LDA-q6" \
 "I GTH-PADE-q7 GTH-LDA-q7" \
 "Xe GTH-PADE-q8 GTH-LDA-q8" \
 "Cs GTH-PADE-q9 GTH-LDA-q9" \
 "Ba GTH-PADE-q10 GTH-LDA-q10" \
 "La GTH-PADE-q11 GTH-LDA-q11" \
 "Ce GTH-PADE-q12 GTH-LDA-q12" \
 "Pr GTH-PADE-q13 GTH-LDA-q13" \
 "Nd GTH-PADE-q14 GTH-LDA-q14" \
 "Pm GTH-PADE-q15 GTH-LDA-q15" \
 "Sm GTH-PADE-q16 GTH-LDA-q16" \
 "Eu GTH-PADE-q17 GTH-LDA-q17" \
 "Gd GTH-PADE-q18 GTH-LDA-q18" \
 "Tb GTH-PADE-q19 GTH-LDA-q19" \
 "Dy GTH-PADE-q20 GTH-LDA-q20" \
 "Ho GTH-PADE-q21 GTH-LDA-q21" \
 "Er GTH-PADE-q22 GTH-LDA-q22" \
 "Tm GTH-PADE-q23 GTH-LDA-q23" \
 "Yb GTH-PADE-q24 GTH-LDA-q24" \
 "Lu GTH-PADE-q25 GTH-LDA-q25" \
 "Hf GTH-PADE-q12 GTH-LDA-q12" \
 "Ta GTH-PADE-q13 GTH-LDA-q13" \
 "W GTH-PADE-q14 GTH-LDA-q14" \
 "Re GTH-PADE-q15 GTH-LDA-q15" \
 "Os GTH-PADE-q16 GTH-LDA-q16" \
 "Ir GTH-PADE-q17 GTH-LDA-q17" \
 "Pt GTH-PADE-q18 GTH-LDA-q18" \
 "Au GTH-PADE-q11 GTH-LDA-q11" \
 "Hg GTH-PADE-q12 GTH-LDA-q12" \
 "Tl GTH-PADE-q13 GTH-LDA-q13" \
 "Pb GTH-PADE-q4 GTH-LDA-q4" \
 "Bi GTH-PADE-q5 GTH-LDA-q5" \
 "Po GTH-PADE-q6 GTH-LDA-q6" \
 "At GTH-PADE-q7 GTH-LDA-q7" \
 "Rn GTH-PADE-q8 GTH-LDA-q8" \
 "H GTH-PBE-q1" \
 "He GTH-PBE-q2" \
 "Li GTH-PBE-q3" \
 "Be GTH-PBE-q4" \
 "B GTH-PBE-q3" \
 "C GTH-PBE-q4" \
 "N GTH-PBE-q5" \
 "O GTH-PBE-q6" \
 "F GTH-PBE-q7" \
 "Ne GTH-PBE-q8" \
 "Na GTH-PBE-q9" \
 "Mg GTH-PBE-q10" \
 "Al GTH-PBE-q3" \
 "Si GTH-PBE-q4" \
 "P GTH-PBE-q5" \
 "S GTH-PBE-q6" \
 "Cl GTH-PBE-q7" \
 "Ar GTH-PBE-q8" \
 "K GTH-PBE-q9" \
 "Ca GTH-PBE-q10" \
 "Sc GTH-PBE-q11" \
 "Ti GTH-PBE-q12" \
 "V GTH-PBE-q13" \
 "Cr GTH-PBE-q14" \
 "Mn GTH-PBE-q15" \
 "Fe GTH-PBE-q16" \
 "Co GTH-PBE-q17" \
 "Ni GTH-PBE-q18" \
 "Cu GTH-PBE-q11" \
 "Zn GTH-PBE-q12" \
 "Ga GTH-PBE-q13" \
 "Ge GTH-PBE-q4" \
 "As GTH-PBE-q5" \
 "Se GTH-PBE-q6" \
 "Br GTH-PBE-q7" \
 "Kr GTH-PBE-q8" \
 "Rb GTH-PBE-q9" \
 "Sr GTH-PBE-q10" \
 "Y GTH-PBE-q11" \
 "Zr GTH-PBE-q12" \
 "Nb GTH-PBE-q13" \
 "Mo GTH-PBE-q14" \
 "Tc GTH-PBE-q15" \
 "Ru GTH-PBE-q16" \
 "Rh GTH-PBE-q17" \
 "Pd GTH-PBE-q18" \
 "Ag GTH-PBE-q11" \
 "Cd GTH-PBE-q12" \
 "In GTH-PBE-q13" \
 "Sn GTH-PBE-q4" \
 "Sb GTH-PBE-q5" \
 "Te GTH-PBE-q6" \
 "I GTH-PBE-q7" \
 "Xe GTH-PBE-q8" \
 "Cs GTH-PBE-q9" \
 "Ba GTH-PBE-q10" \
 "La GTH-PBE-q11" \
 "Ce GTH-PBE-q12" \
 "Pr GTH-PBE-q13" \
 "Nd GTH-PBE-q14" \
 "Pm GTH-PBE-q15" \
 "Sm GTH-PBE-q16" \
 "Eu GTH-PBE-q17" \
 "Gd GTH-PBE-q18" \
 "Tb GTH-PBE-q19" \
 "Dy GTH-PBE-q20" \
 "Ho GTH-PBE-q21" \
 "Er GTH-PBE-q22" \
 "Tm GTH-PBE-q23" \
 "Yb GTH-PBE-q24" \
 "Lu GTH-PBE-q25" \
 "Hf GTH-PBE-q12" \
 "Ta GTH-PBE-q13" \
 "W GTH-PBE-q14" \
 "Re GTH-PBE-q15" \
 "Os GTH-PBE-q16" \
 "Ir GTH-PBE-q17" \
 "Pt GTH-PBE-q18" \
 "Au GTH-PBE-q11" \
 "Hg GTH-PBE-q12" \
 "Tl GTH-PBE-q13" \
 "Pb GTH-PBE-q4" \
 "Bi GTH-PBE-q5" \
 "Po GTH-PBE-q6" \
 "At GTH-PBE-q7" \
 "Rn GTH-PBE-q8" \
 "B GTH-PBESol-q3" \
 "H GTH-OLYP-q1" \
 "O GTH-OLYP-q6"
do
   default_ppnames=""
   token_list=(${label})
   element=${token_list[0]}
   ppnames=${token_list[@]:1}
   for ppname in ${ppnames}
   do
      default_ppnames+=" "${ppname%-q[0-9]*}
   done
   sed -i "s/${element} \+${ppnames}.*/${element} ${ppnames}${default_ppnames}/" ${ppfile}
done

#EOF
