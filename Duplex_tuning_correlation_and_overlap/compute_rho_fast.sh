######
## Code for computing the interlayer degree correlation of a duplex 
## network given the two layers as input 

#!/bin/sh
if [ $# -le 1 ]; then
    echo "Usage: $0 <layer1> <layer2>"
    exit 1
fi
layer1="$1"
layer2="$2"

DEG_SEQ="./deg_seq"
RANK_NODES="./rank_nodes.py"

${DEG_SEQ} ${layer1} > ${layer1}_degs
python ${RANK_NODES} ${layer1}_degs > ${layer1}_rank

${DEG_SEQ} ${layer2} > ${layer2}_degs
python ${RANK_NODES} ${layer2}_degs > ${layer2}_rank

python compute_rho.py ${layer1}_rank ${layer2}_rank
