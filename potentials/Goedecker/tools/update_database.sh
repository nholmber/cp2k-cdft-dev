#!/bin/ksh
# Update the database for the Goedecker-Teter-Hutter pseudo potentials (GTH PP)
# in CPMD format and CP2K/Quickstep format

# Generate the QS files (CP2K/Quickstep format) by converting the XX files
# (CPMD format) in the build/ tree
  xx_to_qs.sh

# Copy all XX files to the database in ../cpmd and
# copy all QS files to the database in ../cp2k/
  collect_files.sh

# Remove all the QS files in the build/ tree
  echo Removing all the QS files in the build/ tree ...
  find ../build -name QS -exec rm {} \;

# Create a new POTENTIAL database file for CP2K/Quickstep w.r.t the
# updated database
  create_cp2k_potential_file.sh

  echo GTH PP database update finished
