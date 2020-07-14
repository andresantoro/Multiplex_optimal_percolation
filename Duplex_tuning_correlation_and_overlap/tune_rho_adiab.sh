#!/bin/sh

####
##
## Create duplex networks in which the inter-layer correlation
## coefficient \rho is tuned adiabatically between the values RHO_MIN
## and RHO_MAX given as input, with the provided increase step
## RHO_STEP.
##
## N.B.: By "adiabatically" we mean that two consecutive duplexes in
## the sequence (i.e. two networks separated by RHO_STEP) are
## characterised by a very small difference in the pairing pattern.
##
##
##


LANG=en_US # This  is required for having  sequence with dot instead of comma e.g
# $ LANG=fr_FR seq 0. 0.1 0.2
# 0,0
# 0,1
# 0,2
# $ LANG=en_US seq 0. 0.1 0.2
# 0.0
# 0.1
# 0.2
if [ $# -le 6 ]; then
    echo "Usage: $0 <layer1> <layer2> <rho_min> <rho_step> <rho_max> <epsilon> <beta> <iteration_folder>" 
    exit 1
fi

L1="$1"
L2="$2"
RHO_MIN=$3
RHO_STEP=$4
RHO_MAX=$5
EPS=$6
BETA=$7


TUNE_RHO="./tune_rho"
RANK_NODES="rank_nodes.py"
DEG_SEQ="./deg_seq"
RELABEL_GRAPHS="relabel_graph.py"

##Computing the degree of each node and save it in a file
${DEG_SEQ} ${L1} > ${L1}_degs
##Computing the rank of each node according to the degree and save it in a file
python ${RANK_NODES} ${L1}_degs > ${L1}_rank

cur_L2=${L2}
val=0.95

# Starting from two uncorrelated systems,   increase the correlation to RHO_MAX (e.g. rho=1)
# Then I go back from the maximum of correlation to negative correlation
# in one way (e.g. from rho =1 to rho=-1. All this machinery because the process
# is not ergodic 

for rho in `seq -f "%.1f" ${RHO_MIN} ${RHO_STEP} ${RHO_MAX}`; do
    label_map=pairing_${L1}_${rho}.txt
    ##Computing the degree sequence and save it in the corresponding file
    ${DEG_SEQ} ${cur_L2} > ${cur_L2}_degs
    ##Computing the ranking of the nodes acorrding to the degree sequence
    ##and save it in the corresponding file
    python ${RANK_NODES} ${cur_L2}_degs > ${cur_L2}_rank
    ## When we are close to MP and MN correlations, it could be useful to change
    ## the EPS tollerance to a greater threshold (the threshold could be too small,
    ## hence it is not possible to obtain the duplex with that value of correlation)
    s_val=`echo "$rho > $val" | bc -l`
    echo $s_val
    if [ $s_val -eq 1 ]
    then
        echo "I'm changing $EPS, since $rho > $val"
        EPS=0.025
    fi
    #Tuning the values of inter-layer degree correlation by relabelling the nodes
    #in one of the two layers (based on the code of MAMMULT - https://github.com/KatolaZ/mammult) 
    ${TUNE_RHO} ${L1}_rank ${cur_L2}_rank ${rho} ${EPS} ${BETA} NAT > ${label_map}
    new_L2=${L2}_relabel_${rho}.txt
    python ${RELABEL_GRAPHS} ${cur_L2} ${label_map} > ${new_L2}
    cur_L2=${new_L2}
    cp ${L1} ${L1}_rho_${rho}
    cp ${new_L2} ${L2}_rho_${rho}
done
${DEG_SEQ} ${cur_L2} > ${cur_L2}_degs
python ${RANK_NODES} ${cur_L2}_degs > ${cur_L2}_rank
echo "CURRENT VALUE OF RHO: ${rho}"
# cp ${L1}_rho_${rho} ${L1}
# cp ${L2}_rho_${rho} ${L2}
rm pairing*




# #Doing the negative correlation
EPS=$6
val=-0.95
counter_steps=1
for rho in `seq -f "%.1f" 0.9 -${RHO_STEP} -${RHO_MAX}`; do
    if [ $counter_steps -eq 10 ]
    then 
        rho=0.0
    fi
    label_map=pairing_${L1}_${rho}.txt
    ${DEG_SEQ} ${cur_L2} > ${cur_L2}_degs
    python  ${RANK_NODES} ${cur_L2}_degs > ${cur_L2}_rank
    ##This for changing the EPS tollerance for the difficult case, a.k.a. when we are close to MP and MN correlations
    s_val=`echo "$rho < $val" | bc -l`
    if [ $s_val -eq 1 ]
    then
        echo "I'm changing $EPS, since $rho > $val"
        EPS=0.025
    fi
    ${TUNE_RHO} ${L1}_rank ${cur_L2}_rank ${rho} ${EPS} ${BETA} NAT > ${label_map}
    new_L2=${L2}_relabel_${rho}.txt
    python ${RELABEL_GRAPHS} ${cur_L2} ${label_map} > ${new_L2}
    cur_L2=${new_L2}
    cp ${L1} ${L1}_rho_${rho}
    cp ${new_L2} ${L2}_rho_${rho}
    counter_steps=$((counter_steps+1))

done
${DEG_SEQ} ${cur_L2} > ${cur_L2}_degs
python ${RANK_NODES} ${cur_L2}_degs > ${cur_L2}_rank


# Saving all the data in a folder
if [ $# -le 8 ];
then
    FOLDER_CREATION=$8
    DIRECTORY="./multiplex_rho_$FOLDER_CREATION"
    echo "$DIRECTORY"
    if [ ! -d "$DIRECTORY" ]; 
    then
        mkdir $DIRECTORY
    fi
    counter=0
    for rho in `seq -1 ${RHO_STEP} 1`; do
        if [ $counter -eq 10 ]
        then
            rho=0
        fi
        # echo "$rho $counter"
        ls -v ${L1}*rho*_${rho} ${L2}*rho*_${rho} > list_multiplex_$counter.txt
        cp ${L1}*rho*_${rho} $DIRECTORY
        cp ${L2}*rho*_${rho} $DIRECTORY
        cp list_multiplex_$counter.txt $DIRECTORY
        counter=$((counter+1))

    done

fi


rm pairing* list_multiplex_* ${L1}*relabel* ${L2}*relabel* *_rank *_degs
