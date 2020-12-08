#! /bin/bash

# all the nabs in CATNAP with > 125 evaluable combination sequences and individual sequences
ALL_NABS=("10E8 + PG9 + VRC07" "10E8 + 3BNC117 + PG9" \
          "10E8 + PG9 + PGT128" "10-1074 + 10E8" "10-1074 + 10E8 + 3BNC117" \
          "10-1074 + 3BNC117" \
          "10-1074 + 3BNC117 + PG9" "10-1074 + PG9" "10E8 + 3BNC117" \
          "10E8 + PGT128 + VRC07" "3BNC117 + PG9" \
          "BG1 + BG18 + NC37" "PG9 + PGT128" \
          "PG9 + PGT128 + VRC07" "PG9 + VRC07" "PGT128 + VRC07")

# for each nab, run slapnap once to create the data set that builds a combination
# learner based on the individual nabs; another that creates a data set that is used
# for out of bag prediction on the measured combination nabs
for NAB in ${ALL_NABS[@]}
do 
  SLAPNAP_NAB=$(echo $NAB | sed 's/+/;/g' | sed 's/ //g')
  DIR_NAB=$(echo $NAB | tr '[:upper:]' '[:lower:]' | sed 's/+/_/g' | sed 's/ //g')
  DATA_NAB=$(echo $NAB | sed 's/ //g')
  mkdir ./docker_output/$DIR_NAB
  mkdir ./docker_output/${DIR_NAB}_2
  docker run \
    -d \
    -v docker_output/$DIR_NAB:/home/output \
    -e learners="rf;lasso;xgboost" \
    -e cvperf="TRUE" \
    -e cvtune="TRUE" \
    -e nab=$SLAPNAP_NAB \
    -e outcomes="estsens;multsens" \
    -e sens_thresh="1" \
    -e var_thresh="0;4" \
    -e return="report;data;learner" \
    -e nfolds="5" \
    slapnap/slapnap

  docker run \
    -d \
    -v docker_output/${DIR_NAB}_2:/home/output \
    -e nab=$DATA_NAB \
    -e outcomes="sens" \
    -e sens_thresh="1" \
    -e return="data" \
    slapnap/slapnap
done