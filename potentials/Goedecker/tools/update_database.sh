#!/bin/sh
# Update the database for the Goedecker-Teter-Hutter pseudo potentials (GTH PP)
# in CPMD format and CP2K/Quickstep format

# Generate the QS files (CP2K/Quickstep format) by converting the XX files
# (CPMD format) in the build/ tree
  xx_to_qs.sh

# Copy all XX files to the database in ../cpmd and
# copy all QS files to the database in ../cp2k/
  collect_files.sh

# Remove all the QS files in the build/ tree
  echo Removing all the INFO files in the build/ tree ...
  find ../build -name INFO -exec rm {} \;
  echo Removing all the QS files in the build/ tree ...
  find ../build -name QS -exec rm {} \;
  echo Removing all the CPMD files in the build/ tree ...
  find ../build -name CPMD -exec rm {} \;
  echo Removing all the TEXTAB files in the build/ tree ...
  find ../build -name TEXTAB -exec rm {} \;

# Create a new POTENTIAL database file for CP2K/Quickstep w.r.t the
# updated database
  create_cp2k_potential_file.sh
  create_tex_file.sh

  echo GTH PP database update finished
