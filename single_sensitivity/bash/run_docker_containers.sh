#! /bin/bash

# all the single nabs in CATNAP that were run by Hake and Pfeifer (2017) and/or Rawi et al. (2019)
ALL_NABS=( "2G12" "PG16" \
"PG9" "PGT145" \
"PGDM1400" "VRC26.08" \
"VRC26.25" "PGT128" \
"10-1074" "10-996" \
"PGT121" "VRC38.01" \
"PGT135" "DH270.1" \
"DH270.5" "DH270.6" \
"VRC29.03" "VRC01" \
"3BNC117" "VRC-PG04" \
"NIH45-46" "VRC03" \
"VRC-CH31" "CH01" \
"HJ16" "VRC07" \
"b12" "PGT151" \
"VRC34.01" "8ANC195" \
"35O22" "2F5" "4E10" )

# for each nab, run slapnap
for NAB in "${ALL_NABS[@]}";
do
  SLAPNAP_NAB=$(echo $NAB | sed 's/+/;/g' | sed 's/ //g')
  DIR_NAB=$(echo $NAB | tr '[:upper:]' '[:lower:]' | sed 's/+/_/g' | sed 's/ //g')
  DATA_NAB=$(echo $NAB | sed 's/ //g')
  mkdir -p ./docker_output/$DIR_NAB
  docker run \
    -d \
    -v "$(pwd)"/docker_output/$DIR_NAB/:/home/output/ \
    -e learners="rf;lasso;xgboost" \
    -e cvperf="TRUE" \
    -e cvtune="TRUE" \
    -e nab=$SLAPNAP_NAB \
    -e outcomes="sens" \
    -e sens_thresh="50" \
    -e return="report;data;learner" \
    -e nfolds="5" \
    slapnap
done
