####################################################
# path - path or states
# l,k -  given state element
# G - observed sequence
# lkl - observed GenoLkl
# b - given sequence element

# a - transition probabillities
# e - emission probabilities

# i - iterator
############################ FUNCTIONS ############################

get_IBD <- function(pos_dist, indF=0, a=0.1){
  state <- sample(0:1, 1, prob=c(1-indF,indF))
  for(s in 2:length(pos_dist))
    state[s] <- sample(0:1, 1, prob=calc_trans(state[s-1],pos_dist[s],indF,a))
  
  return(state)
}

calc_trans <- function(state, pos_dist, indF, a){
  X <- exp(-a*pos_dist)
  
  IBD_01 <- (1-X)*indF
  IBD_00 <- 1-IBD_01
  IBD_10 <- (1-X)*(1-indF)
  IBD_11 <- 1-IBD_10
  
  if(state == 0)
    prob <- c(IBD_00,IBD_01)
  else
    prob <- c(IBD_10,IBD_11)
  
  return(prob)
}

get_haplotype <- function(snp, freq){
  if(length(freq) == snp)
    p <- freq
  else
    p <- sample(freq, snp, replace=TRUE)
  
  prob <- t(matrix( c(p,1-p), nrow=length(p) ))
  apply(prob, 2, sample, x=c(1,0), size=1, replace=TRUE)
}

getLikes <- function(x,depth=5,error=0.01,norm=TRUE,loglikeR=FALSE){
  n <- length(x)
  depth <- rpois(n,depth)
  nA <- rbinom(n,depth,c(error,0.5,1-error)[x+1])
  res <- rbind(dbinom(nA,depth,error),dbinom(nA,depth,0.5),dbinom(nA,depth,1-error))
  
  # Normalize values in [0,1]
  if(norm){
    sum <- apply(res, 2, sum)
    res <- t(apply(res, 1, function(x){x/sum}))
  }
  
  # Standardise to log like ratios
  if(loglikeR){
    max <- apply(res, 2, max)
    res <- t(apply(res, 1, function(x){x/max}))
  }
  
  return(log(res))
}

#####  Parse command-line arguments
library(optparse)
option_list <- list(make_option(c("-n", "--n_ind"), action="store", type="integer", default=10, help="Number of individuals [%default]"),
                    make_option(c("-l", "--n_sites"), action="store", type="integer", default=1000, help="Number of independent sites [%default]"),
                    make_option(c("-F", "--indF"), action="store", type="character", default="0", help="Per-individual inbreeding coefficients [%default]"),
                    make_option(c("-f", "--freq"), action="store", type="character", default="0.1", help="Allele frequencies [%default]"),
                    make_option(c("-x", "--pos_dist"), action="store", type="character", default="1", help="Distance between sites [%default]"),
                    make_option(c("-t", "--trans"), action="store", type="character", default="0.01", help="Transition probabilities [%default]"),
                    make_option(c("-d", "--depth"), action="store", type="character", default="5", help="Sequencing depth [%default]"),
                    make_option(c("-e", "--error"), action="store", type="numeric", default=0.01, help="Error rate [%default]"),
                    make_option(c("-s", "--seed"), action="store", type="integer", default=12345, help="Seed for random number generator [%default]"),
                    make_option(c("-o", "--out"), action="store", type="character", default="sim", help="Output prefix [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

############################ Parsing input arguments ############################
set.seed(opt$seed)

### indF
indF = c()
cat("====> Parsing per-individual indF...",fill=TRUE)
if(file.exists(opt$indF)){
  cat("==> Reading values from file:",opt$indF,fill=TRUE)
  indF <- as.numeric(readLines(opt$indF))
  if(opt$n_ind != length(indF)){
    cat("ERROR: number of individuals and indF file do not match", fill=TRUE)
    quit("no",-1)
  }
}else if(opt$indF == "r"){
  cat("==> Using random values",fill=TRUE)
  indF <- runif(opt$n_ind)
}else{
  cat("==> Setting to fixed value:",as.numeric(opt$indF),fill=TRUE)
  for (i in 1:opt$n_ind)
    indF[i] <- as.numeric(opt$indF)
}



### freq
freq = c()
cat("====> Parsing allele frequencies...",fill=TRUE)
if(file.exists(opt$freq)){
  cat("==> Reading values from file:",opt$freq,fill=TRUE)
  freq <- as.numeric(readLines(opt$freq))
  if(opt$n_sites != length(freq)){
    cat("ERROR: number of sites and FREQ file do not match", fill=TRUE)
    quit("no",-1)
  }
}else if(opt$freq == "r"){
  cat("==> Using random values",fill=TRUE)
  freq <- runif(opt$n_sites)
}else{
  cat("==> Setting to fixed value:",as.numeric(opt$freq),fill=TRUE)
  for (s in 1:opt$n_sites)
    freq[s] <- as.numeric(opt$freq)
}



### trans
trans = c()
cat("====> Parsing per-individual transition probability...",fill=TRUE)
if(file.exists(opt$trans)){
  cat("==> Reading values from file:",opt$trans,fill=TRUE)
  trans <- as.numeric(readLines(opt$trans))
  if(opt$n_ind != length(trans)){
    cat("ERROR: number of individuals and transitions file do not match", fill=TRUE)
    quit("no",-1)
  }
}else if(opt$trans == "r"){
  cat("==> Using random values",fill=TRUE)
  trans <- runif(opt$n_ind)
}else{
  cat("==> Setting to fixed value:",as.numeric(opt$trans),fill=TRUE)
  for (i in 1:opt$n_ind)
    trans[i] <- as.numeric(opt$trans)
}



### depth
depth = c()
cat("====> Parsing per-individual depth...",fill=TRUE)
if(file.exists(opt$depth)){
  cat("==> Reading values from file:",opt$depth,fill=TRUE)
  depth <- as.numeric(readLines(opt$depth))
  if(opt$n_ind != length(depth)){
    cat("ERROR: number of individuals and DEPTH file do not match", fill=TRUE)
    quit("no",-1)
  }
}else if(opt$depth == "r"){
  cat("==> Using random values",fill=TRUE)
  depth <- runif(opt$n_ind)
}else{
  cat("==> Setting to fixed value:",as.numeric(opt$depth),fill=TRUE)
  for (i in 1:opt$n_ind)
    depth[i] <- as.numeric(opt$depth)
}



### Parse distance between sites
pos_dist = c()
cat("====> Parsing distance between sites...",fill=TRUE)
if(file.exists(opt$pos_dist)){
  cat("==> Reading values from file:",opt$pos_dist,fill=TRUE)
  pos_dist <- as.numeric(readLines(opt$pos_dist))
  if(opt$n_sites != length(pos_dist)){
    cat("ERROR: number of sites and DIST file do not match", fill=TRUE)
    quit("no",-1)
  }
}else if(opt$pos_dist == "r"){
  avg_dist <- 1e6
  cat("==> Setting to normally distributed with mean:",avg_dist,fill=TRUE)
  pos_dist <- as.integer(rnorm(opt$n_sites,mean=avg_dist,sd=avg_dist/10))

  # Print sites' positions
  dist_out <- paste(opt$out,"dist.gz",sep=".")
  if(!is.na(file.info(dist_out)[,"size"])){
    warning("WARN: DIST file already exists. Skipping...");
  } else {
    fh <- gzfile(dist_out, "w", compression=9)
    write.table(paste("chrS",cumsum(as.numeric(pos_dist)),sep="\t"), fh, quote=FALSE, sep="\n", row.names=FALSE, col.names=FALSE)
    close(fh)
  }
}else{
  cat("==> Setting to fixed value:",as.numeric(opt$pos_dist),fill=TRUE)
  for (s in 1:opt$n_sites)
    pos_dist[s] <- as.numeric(opt$pos_dist)
}
# Convert distance between positions to Mb
pos_dist <- pos_dist / 1e6;



############################ Generate data ############################

############ Generating TRUE (hidden) path
cat("====> Generating true path...",fill=TRUE)
path = list();
for (i in 1:opt$n_ind)
  path[[i]] = get_IBD(pos_dist, indF[i], trans[i])


# Print
path_out <- paste(opt$out,"path.gz",sep=".")
if(!is.na(file.info(path_out)[,"size"])){
  warning("WARN: PATH file already exists. Skipping...");
} else {
  fh <- gzfile(path_out, "w", compression=9)
  write.table(t(as.data.frame(path)), fh, quote=FALSE, sep="", row.names=FALSE, col.names=FALSE)
  close(fh)
}



############ Generating KNOWN (observed) sequence (G)
cat("====> Generating true genotypes...",fill=TRUE)
geno = list();
for (i in 1:opt$n_ind){
  hap1 <- get_haplotype(opt$n_sites, freq)
  hap2 <- get_haplotype(opt$n_sites, freq)
  # Get IBD haplotypes
  hap1[path[[i]]==1] <- hap2[path[[i]]==1]
  # Get genotypes
  geno[[i]] = hap1+hap2
}

# Print
seq_out <- paste(opt$out,"geno.gz",sep=".")
if(!is.na(file.info(seq_out)[,"size"])){
  warning("WARN: SEQ file already exists. Skipping...");
} else {
  fh <- gzfile(seq_out, "w", compression=9)
  write.table(geno, fh, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  close(fh)
}



############ Calculate genotype likelihoods
cat("====> Calculating genotype likelihoods...",fill=TRUE)
geno_lkl <- geno
for (i in 1:opt$n_ind)
  geno_lkl[[i]] = getLikes(geno[[i]],depth[i],opt$error,T,F)

# Print
lkl_out <- paste(opt$out,"glf.gz",sep=".")
if(!is.na(file.info(lkl_out)[,"size"])){
  warning("WARN: LKL file already exists. Skipping...");
} else {
  fh <- gzfile(lkl_out, "w", compression=9)
  for (s in 1:opt$n_sites) {
    for (i in 1:opt$n_ind)
      writeLines(as.character(geno_lkl[[i]][,s]), fh, "\t")
    write("", fh)
  }
  close(fh)
}
