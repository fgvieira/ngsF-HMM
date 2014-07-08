SIM_DATA=../../ngsTools/ngsSim/examples
ANGSD=../../ngsTools/angsd
NGSPOPGEN=../../ngsTools/ngsPopGen


##### Clean-up
rm -f testF*


N_IND=10
N_SITES=10000
FREQ=0.2
INDF=0.5
TRANS=0.01
SEED=12345

##### Simulate data from HMM
echo "========== Simulating data ==========" >&2
DEPTH=2
ERROR=0.01
Rscript ../R/ngsSim-HMM.R --n_ind $N_IND --n_sites $N_SITES --indF $INDF --freq $FREQ --trans $TRANS --depth $DEPTH --error $ERROR --seed $SEED --out testF-HMM.SIM >&2

TRANS=0.005-0.005

##### Infer F
for TYPE in TG GL GL_CG
do
    echo "========== $TYPE =========="
    
    FILE="testF-HMM.SIM.geno.gz"
    if [[ $TYPE == "GL" ]]; then
	FILE="testF-HMM.SIM.glf.gz --loglkl"
    fi
    if [[ $TYPE == "GL_CG" ]]; then
	FILE="testF-HMM.SIM.glf.gz --loglkl --call_geno"
    fi
    
    ID=TRUE
    ../ngsF-HMM -verbose 2 -n_threads 10 --seed $SEED --geno $FILE --n_ind $N_IND --n_sites $N_SITES --freq $FREQ --freq_fixed --trans $TRANS --trans_fixed --path sim.path.gz --path_fixed --out testF-HMM.$ID.$TYPE --log 1

    ID=BEST
    ../ngsF-HMM -verbose 2 -n_threads 10 --seed $SEED --geno $FILE --n_ind $N_IND --n_sites $N_SITES --freq $FREQ --trans $TRANS --path sim.path.gz --out testF-HMM.$ID.$TYPE --log 1

    ID=freq_fixed
    ../ngsF-HMM -verbose 2 -n_threads 10 --seed $SEED --geno $FILE --n_ind $N_IND --n_sites $N_SITES --freq $FREQ --freq_fixed --trans 0.1,0.1 --path 0 --out testF-HMM.$ID.$TYPE --log 1

    ID=trans_fixed
    ../ngsF-HMM -verbose 2 -n_threads 10 --seed $SEED --geno $FILE --n_ind $N_IND --n_sites $N_SITES --freq 0.1 --trans $TRANS --trans_fixed --path 0 --out testF-HMM.$ID.$TYPE --log 1

    ID=path_fixed
    ../ngsF-HMM -verbose 2 -n_threads 10 --seed $SEED --geno $FILE --n_ind $N_IND --n_sites $N_SITES --freq 0.1 --trans 0.1,0.1 --path sim.path.gz --path_fixed --out testF-HMM.$ID.$TYPE --log 1

    ID=normal
    ../ngsF-HMM -verbose 2 -n_threads 10 --seed $SEED --geno $FILE --n_ind $N_IND --n_sites $N_SITES --freq 0.1 --trans 0.1,0.1 --path 0 --out testF-HMM.$ID.$TYPE --log 1


    #echo "===== Plot ====="
    #for ID in TRUE BEST freq_fixed trans_fixed path_fixed normal
    #do
    #    Rscript ../R/ngsF-HMMplot.R --in_file testF-HMM.$ID.$TYPE.log.gz --n_ind $N_IND --n_sites $N_SITES --geno testF-HMM.SIM.geno.gz --path testF-HMM.SIM.path.gz --out testF-HMM.$ID.$TYPE
    #done
done >&2



##### Get genotype likelihoods
N_IND=20
$ANGSD/angsd -glf $SIM_DATA/testF.glf.gz -fai $SIM_DATA/testAF.ANC.fai -nInd $N_IND -doMajorMinor 1 -doGlf 2 -doMaf -1 -SNP_pval 1e-4 -out testF
$ANGSD/angsd -glf $SIM_DATA/testF.glf.gz -fai $SIM_DATA/testAF.ANC.fai -nInd $N_IND -doMajorMinor 1 -doGlf 3 -doMaf -1 -SNP_pval 1e-4 -out testF
gunzip testF.glf.gz



##### Estimate F
N_SITES=$((`zcat testF.beagle.gz | wc -l`-1))
../ngsF-HMM --verbose 2 -n_threads 10 --seed $SEED --geno testF.beagle.gz --lkl --n_ind $N_IND --n_sites $N_SITES --freq 0.1 --trans 0.1,0.1 --path 0 --out testF --log 1 >&2
../ngsF-HMM --verbose 2 -n_threads 10 --seed $SEED --geno testF.glf --loglkl --n_ind $N_IND --n_sites $N_SITES --freq 0.1 --trans 0.1,0.1 --path 0 --out testF_bin --log 1 >&2



##### Get genotypes' posterior probability with inbreeding prior
head -n $((N_IND+1)) testF.indF | tail -n $N_IND | cut -f 1 > /tmp/testF.indF
$ANGSD/angsd -glf $SIM_DATA/testF.glf.gz -fai $SIM_DATA/testAF.ANC.fai -nInd $N_IND -doMajorMinor 1 -doPost 1 -doMaf -1 -indF /tmp/testF.indF -doGeno 32 -doSaf 2 -anc $SIM_DATA/testAF.ANC.fas -out testF.indF



##### Calculate covariance matrix
gunzip testF.indF.geno.gz
$NGSPOPGEN/ngsCovar -probfile testF.indF.geno -outfile testF.indF.covar -nind $N_IND -nsites $N_SITES -call 0 -sfsfile testF.indF.saf -norm 0



##### Calculate population genetics statistics
$NGSPOPGEN/ngsStat -npop 1 -postfiles testF.indF.saf -nsites $N_SITES -iswin 1 -nind $N_IND -outfile testF.indF.stat -isfold 0 -islog 1 -block_size $N_SITES



##### SFS
# Calculating folded SFS
cat testF.indF.saf | hexdump -v -e "$((N_IND+1))/8 \"%.10g\t\"\"\n\"" | perl -na -e '$sum=0; $sum+=exp($_) for @F; next if($sum==0); for $i (0..$#F){$frq[$i]+=exp($F[$i])/$sum}; END{$tsum+=$_ for @frq; $_/=$tsum for @frq; print join("\t",@frq)."\n"}' > testF.indF.fold-saf_sum

# Calculating unfolded SFS
cat testF.indF.saf | hexdump -v -e "$((2*N_IND+1))/8 \"%.10g\t\"\"\n\"" | perl -na -e '$sum=0; $sum+=exp($_) for @F; next if($sum==0); for $i (0..$#F){$frq[$i]+=exp($F[$i])/$sum}; END{$tsum+=$_ for @frq; $_/=$tsum for @frq; print join("\t",@frq)."\n"}' > testF.indF.saf_sum





##### Check MD5
rm -f *.arg
md5sum testF* | sort -k 2,2 > /tmp/test.md5
if diff /tmp/test.md5 test.md5 > /dev/null
then
    echo "ngsF-HMM: All tests OK!"
else
    echo "ngsF-HMM: test(s) failed!"
fi
