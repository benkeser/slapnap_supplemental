#! /bin/bash

# make a file to be read by R with all inputs
[ -e bash_output/all_nabs.txt ] && rm bash_output/all_nabs.txt
touch bash_output/all_nabs.txt
for NAB in "2G12" "PG16" \
          "PG9" "PGT145" \
          "PGDM1400" "VRC26.08" \
          "VRC26.25" "PGT128" \
          "10-1074" "10-996" \
          "PGT121" "VRC38.01" \
          "PGT135" "DH270.1" \
          "DH270.5" "DH270.6" \
          "VRC29.03" "VRC01"  \
          "3BNC117" "VRC-PG04" \
          "NIH45-46" "VRC03" \
          "VRC-CH31" "CH01" \
          "HJ16" "VRC07" \
          "b12" "PGT151" \
          "VRC34.01" "8ANC195" \
          "35O22" "2F5" \
          "4E10" 
do
  DIR_NAB=$(echo $NAB | tr '[:upper:]' '[:lower:]' | sed 's/+/_/g' | sed 's/ //g')
  echo $DIR_NAB >> bash_output/all_nabs.txt
done
