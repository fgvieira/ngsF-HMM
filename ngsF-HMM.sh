#!/bin/bash


######### Functions ##########
in_array() {
    idx=""
    local CNT=0
    local hay needle=$1
    shift
    for hay; do
        CNT=$((CNT+1))
        if [[ $hay == $needle ]]; then
            idx=$CNT
            return 1
        fi
    done
    return 0
}


########## Variables ##########
TMP_DIR=/tmp
N_REP=20
ID=ngsF-HMM_$RANDOM


##############################
args=( $@ )

# find -o/-out/--out argument
in_array "--out" "${args[@]}"
if [[ $idx -eq 0 ]]; then
    in_array "-out" "${args[@]}"
fi
if [[ $idx -eq 0 ]]; then
    in_array "-o" "${args[@]}"
fi
if [[ $idx -eq 0 ]]; then
    echo "ERROR: could not find argument for output files (-o / --out)"
    exit -1
fi
OUT=${args[$idx]}


# Run each replicate
for REP in `seq -w 1 $N_REP`
do
    args[$idx]=$TMP_DIR/$ID.REP_$REP
    ${0%\.sh} ${args[@]}
done


# Find best replicate
BEST=`awk 'FNR==1{print FILENAME"\t"$1}' $TMP_DIR/$ID.REP_*.indF | sort -k 2,2gr | awk 'NR==1{sub(".indF","",$1); print $1}'`

# Get best replicate
mv $BEST.indF $OUT.indF
mv $BEST.ibd $OUT.ibd
mv $BEST.geno $OUT.geno
if [[ -s $BEST.log.gz ]]; then
    mv $BEST.log.gz $OUT.log.gz
fi

# Clean-up
rm -f $TMP_DIR/$ID.REP_*
