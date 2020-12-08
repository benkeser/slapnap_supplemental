#! /bin/bash

# make a file to be read by R with all inputs
[ -e bash_output/all_nabs.txt ] && rm bash_output/all_nabs.txt
touch bash_output/all_nabs.txt
# flagged for follow-up 
# "10E8 + PG9 + PGT128 + VRC07" "10E8v4-5R-100cF + N6"
for NAB in "10E8 + PG9 + VRC07" "10E8 + 3BNC117 + PG9" \
          "10E8 + PG9 + PGT128" "10-1074 + 10E8" "10-1074 + 10E8 + 3BNC117" \
          "10-1074 + 10E8 + 3BNC117 + PG9" "10-1074 + 10E8 + PG9" "10-1074 + 3BNC117" \
          "10-1074 + 3BNC117 + PG9" "10-1074 + PG9" "10E8 + 3BNC117" "10E8 + PGT128 + VRC07" \
          "3BNC117 + PG9" "BG1 + BG18 + NC37" "PG9 + PGT128" \
          "PG9 + PGT128 + VRC07" "PG9 + VRC07" "PGT128 + VRC07"
do
  DIR_NAB=$(echo $NAB | tr '[:upper:]' '[:lower:]' | sed 's/+/_/g' | sed 's/ //g')
  echo $DIR_NAB >> bash_output/all_nabs.txt
done