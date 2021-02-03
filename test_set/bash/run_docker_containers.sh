# Run 3BNC117 + 10-1074 on CATNAP, predicting quantitative and binary IC-50
ALL_THRESH=("2", "1")
NAB="10-1074 + 3BNC117"
SLAPNAP_NAB=$(echo $NAB | sed 's/+/;/g' | sed 's/ //g')
for THRESH in ${ALL_THRESH[@]}
do
    DIR_NAB="$(echo $NAB | tr '[:upper:]' '[:lower:]' | sed 's/+/_/g' | sed 's/ //g')_$(echo $THRESH)"
    if [ $THRESH == "2" ]; then
        outcomes="ic50;estsens;multsens"
    else
        outcomes="estsens;multsens"
    fi
    mkdir -p ./docker_output/$DIR_NAB
    sudo docker run \
        -d \
        -v "$(pwd)"/docker_output/$DIR_NAB:/home/output/ \
        -e nab=$SLAPNAP_NAB \
        -e outcomes=$outcomes \
        -e binary_outcomes="ic50" \
        -e learners="rf;lasso;xgboost" \
        -e sens_thresh=$THRESH \
        -e var_thresh="0;4" \
        -e nfolds="5" \
        -e cvtune="TRUE" \
        -e cvperf="TRUE" \
        -e return="report;data;learner" \
        slapnap/slapnap
done
