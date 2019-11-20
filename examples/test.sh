SIM_DATA=../../ngsSim/examples
ANGSD=../../angsd
NGSPOPGEN=../../ngsPopGen


##### Clean-up
rm -f testF*


N_IND=10
N_SITES=10000
FREQ=0.2
SITE_POS=r
INDF=0.5
ALPHA=0.01
SEED=12345

##### Simulate data from HMM
echo "========== Simulating data ==========" >&2
DEPTH=2
ERROR=0.01

Rscript ../scripts/ngsF-HMMsim.R --n_ind $N_IND --n_sites $N_SITES --freq $FREQ --site_pos $SITE_POS --indF $INDF --alpha $ALPHA --depth $DEPTH --error $ERROR --seed $SEED --out testF-HMM.SIM >&2



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
    ../ngsF-HMM --verbose 2 --n_threads $N_IND --seed $SEED --geno $FILE --n_ind $N_IND --n_sites $N_SITES --pos testF-HMM.SIM.pos.gz --freq $FREQ --freq_est 0 --indF $INDF,$ALPHA --indF_fixed --out testF-HMM.$ID.$TYPE --log 1

    ID=BEST
    ../ngsF-HMM --verbose 2 --n_threads $N_IND --seed $SEED --geno $FILE --n_ind $N_IND --n_sites $N_SITES --pos testF-HMM.SIM.pos.gz --freq $FREQ --indF $INDF,$ALPHA --out testF-HMM.$ID.$TYPE --log 1

    ID=freq_fixed
    ../ngsF-HMM --verbose 2 --n_threads $N_IND --seed $SEED --geno $FILE --n_ind $N_IND --n_sites $N_SITES --pos testF-HMM.SIM.pos.gz --freq $FREQ --freq_est 0 --indF 0.1,0.2 --out testF-HMM.$ID.$TYPE --log 1

    ID=indF_fixed
    ../ngsF-HMM --verbose 2 --n_threads $N_IND --seed $SEED --geno $FILE --n_ind $N_IND --n_sites $N_SITES --pos testF-HMM.SIM.pos.gz --freq 0.1 --indF $INDF,$ALPHA --indF_fixed --out testF-HMM.$ID.$TYPE --log 1

    ID=normal
    ../ngsF-HMM --verbose 2 --n_threads $N_IND --seed $SEED --geno $FILE --n_ind $N_IND --n_sites $N_SITES --pos testF-HMM.SIM.pos.gz --freq 0.1 --indF 0.1,0.2 --out testF-HMM.$ID.$TYPE --log 1


    echo "===== Plot ====="
    for ID in TRUE BEST freq_fixed indF_fixed normal
    do
        Rscript ../scripts/ngsF-HMMplot.R --in_file testF-HMM.$ID.$TYPE.ibd --n_ind $N_IND --n_sites $N_SITES --geno testF-HMM.SIM.geno.gz --path testF-HMM.SIM.path.gz --pos testF-HMM.SIM.pos.gz --marg_prob --out testF-HMM.$ID.$TYPE.pdf
    done
done >&2



##### Get genotype likelihoods
N_IND=20
$ANGSD/angsd -glf $SIM_DATA/testF.glf.gz -fai $SIM_DATA/testAF.ANC.fas.fai -nInd $N_IND -doMajorMinor 1 -doGlf 2 -doMaf 1 -SNP_pval 1e-4 -out testF
$ANGSD/angsd -glf $SIM_DATA/testF.glf.gz -fai $SIM_DATA/testAF.ANC.fas.fai -nInd $N_IND -doMajorMinor 1 -doGlf 3 -doMaf 1 -SNP_pval 1e-4 -out testF
gunzip -f testF.glf.gz



##### Estimate F
zcat testF.beagle.gz | tail -n +2 | cut -f 1 | tr "_" "\t" > testF.pos
N_SITES=`cat testF.pos | wc -l`
../ngsF-HMM --verbose 2 -n_threads $N_IND --seed $SEED --geno testF.beagle.gz --lkl --n_ind $N_IND --n_sites $N_SITES --pos testF.pos --freq 0.1 --indF 0.1,0.2 --max_iters 20 --out testF     --log 1 >&2
../ngsF-HMM --verbose 2 -n_threads $N_IND --seed $SEED --geno testF.glf    --loglkl --n_ind $N_IND --n_sites $N_SITES --pos testF.pos --freq 0.1 --indF 0.1,0.2 --max_iters 20 --out testF_bin --log 1 >&2



##### Get genotypes' posterior probability with inbreeding prior
head -n $((N_IND+1)) testF.indF | tail -n $N_IND | cut -f 1 > /tmp/testF.indF
$ANGSD/angsd -glf $SIM_DATA/testF.glf.gz -fai $SIM_DATA/testAF.ANC.fas.fai -nInd $N_IND -doMajorMinor 1 -doPost 1 -doMaf 1 -indF /tmp/testF.indF -doGeno 32 -doSaf 2 -anc $SIM_DATA/testAF.ANC.fas -out testF.indF



##### Calculate covariance matrix
gunzip -f testF.indF.geno.gz testF.indF.saf.gz
$NGSPOPGEN/ngsCovar -probfile testF.indF.geno -nind $N_IND -nsites $N_SITES -call 0 -sfsfile testF.indF.saf -norm 0 -outfile testF.indF.covar



##### Calculate population genetics statistics
$NGSPOPGEN/ngsStat -npop 1 -postfiles testF.indF.saf -nsites $N_SITES -iswin 1 -nind $N_IND -block_size $N_SITES -outfile testF.indF.stat 



##### SFS
hexdump -v -s 8 -e "$((2*N_IND+1))/4 \"%.10g\t\"\"\n\"" testF.indF.saf | perl -na -e '$sum=0; $sum+=exp($_) for @F; next if($sum==0); for $i (0..$#F){$frq[$i]+=exp($F[$i])/$sum}; END{$tsum+=$_ for @frq; $_/=$tsum for @frq; print join("\t",@frq)."\n"}' > testF.indF.saf_sum



##### Check MD5
rm -f *.arg
TMP=`mktemp --suffix .ngsF-HMM`
md5sum testF* | sort -k 2,2 | fgrep -v '.pdf' | fgrep -v '.log.gz' > $TMP
if diff $TMP test.md5 > /dev/null
then
    echo "ngsF-HMM: All tests OK!"
else
    echo "ngsF-HMM: test(s) failed!"
fi
