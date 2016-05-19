#!/bin/bash
shopt -s extglob



#################
### Variables ###
#################
N_REP=10
TMP_DIR=$HOME/scratch/ngsF-HMM



#################
### Functions ###
#################
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



#######################
### Check arguments ###
#######################
args=( $@ )
ID=ngsF-HMM_$RANDOM
mkdir -p $TMP_DIR

# find -S/-seed/--seed argument
in_array "--seed" "${args[@]}"
if [[ $idx -eq 0 ]]; then
    in_array "-seed" "${args[@]}"
fi
if [[ $idx -eq 0 ]]; then
    in_array "-S" "${args[@]}"
fi
if [[ $idx -ne 0 ]]; then
    RANDOM=${args[$idx]}
    SEED_idx=$idx
fi

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



##########################
### Run each replicate ###
##########################
for REP in `seq -w 1 $N_REP`
do
    if [[ $SEED_idx -gt 0 ]]; then
	args[$SEED_idx]=$RANDOM
    fi
    args[$idx]=$TMP_DIR/$ID.REP_$REP
    echo ${0%\.sh} ${args[@]}
done
exit


##########################
### Get best replicate ###
##########################
if [ -s $TMP_DIR/$ID.REP_*(0)1.indF ]; then
    # Find best replicate
    BEST=`awk 'FNR==1{print FILENAME"\t"$1}' $TMP_DIR/$ID.REP_*.indF | sort -k 2,2gr | awk 'NR==1{sub(".indF","",$1); print $1}'`

    if [ -s $BEST.indF ]; then
	# Get best replicate
	mv $BEST.indF $OUT.indF
	mv $BEST.ibd $OUT.ibd
	mv $BEST.geno $OUT.geno
	if [[ -s $BEST.log.gz ]]; then
	    mv $BEST.log.gz $OUT.log.gz
	fi
    else
	echo "ERROR: invalid BEST output files"
	exit -1
    fi
else
    echo "ERROR: no output files found"
    exit -1
fi	

# Clean-up
rm -f $TMP_DIR/$ID.REP_*
