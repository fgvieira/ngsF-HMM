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

getLikes<-function(x,depth=5,error=0.01,norm=TRUE,loglikeR=FALSE){
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
                    make_option(c("-t", "--trans"), action="store", type="numeric", default=0.01, help="Transition probabilities [%default]"),
                    make_option(c("-d", "--depth"), action="store", type="character", default="5", help="Sequencing depth [%default]"),
                    make_option(c("-e", "--error"), action="store", type="numeric", default=0.01, help="Error rate [%default]"),
                    make_option(c("-s", "--seed"), action="store", type="integer", default=12345, help="Seed for random number generator [%default]"),
                    make_option(c("-o", "--out"), action="store", type="character", default="sim", help="Output prefix [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

############################ Parsing input arguments ############################

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
  for (i in 1:opt$n_ind) indF[i] <- as.numeric(opt$indF)
}



### freq
freq = c()
cat("====> Parsing initial allele frequencies...",fill=TRUE)
if(file.exists(opt$freq)){
  cat("==> Reading values from file:",opt$freq,fill=TRUE)
  freq <- as.numeric(readLines(opt$freq))
  if(opt$n_sites != length(freq)){
    cat("ERROR: number of sites and freq file do not match", fill=TRUE)
    quit("no",-1)
  }
}else if(opt$freq == "r"){
  cat("==> Using random values",fill=TRUE)
  freq <- runif(opt$n_sites)
}else{
  cat("==> Setting to fixed value:",as.numeric(opt$freq),fill=TRUE)
  for (i in 1:opt$n_sites) freq[i] <- as.numeric(opt$freq)
}


### depth
depth = c()
cat("====> Parsing per-individual depth...",fill=TRUE)
if(file.exists(opt$depth)){
  cat("==> Reading values from file:",opt$depth,fill=TRUE)
  depth <- as.numeric(readLines(opt$depth))
  if(opt$n_ind != length(depth)){
    cat("ERROR: number of individuals and depth file do not match", fill=TRUE)
    quit("no",-1)
  }
}else if(opt$depth == "r"){
  cat("==> Using random values",fill=TRUE)
  depth <- runif(opt$n_ind)
}else{
  cat("==> Setting to fixed value:",as.numeric(opt$depth),fill=TRUE)
  for (i in 1:opt$n_ind) depth[i] <- as.numeric(opt$depth)
}



############################ Generate data ############################
set.seed(opt$seed)
### Parameters
n_states=2



############ Populating TRUE transition probabilities
cat("====> Initializing transition probabilities...",fill=TRUE)
a = list()
for (j in 1:opt$n_ind)
  a[[j]] = matrix(c(1-opt$trans*indF[j],opt$trans*indF[j],(1-indF[j])*opt$trans,1-(1-indF[j])*opt$trans),ncol=n_states,byrow=T,dimnames=list(c("nIBD","IBD"),c("nIBD","IBD")))



############ Populating TRUE emission probs
cat("====> Initializing emission probabilities...",fill=TRUE)
e = list()
for (i in 1:opt$n_sites) {
  e[[i]] = matrix(nrow=n_states,ncol=3,dimnames=list(c("nIBD","IBD"),c("AA","Aa","aa")))
  for (k in 1:2){
    f=freq[i]
    e[[i]][k,1] = (1-f)^2+(1-f)*f*(k-1)
    e[[i]][k,2] = 2*(1-f)*f-2*(1-f)*f*(k-1)
    e[[i]][k,3] = f^2+(1-f)*f*(k-1)
  }
}



############ Generating TRUE (hidden) path
cat("====> Generating true path...",fill=TRUE)
path = list();
for (j in 1:opt$n_ind) path[[j]] = numeric(opt$n_sites)

for (j in 1:opt$n_ind){
  path[[j]][1] = sample(0:1,size=1,prob=c(1-indF[[j]],indF[[j]]))
  
  for (i in 2:opt$n_sites)
    path[[j]][i] = sample(0:1,size=1,prob=a[[j]][path[[j]][i-1]+1,])
}

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
geno=list();
for (j in 1:opt$n_ind) geno[[j]] = numeric(opt$n_sites)

for (j in 1:opt$n_ind)
  for (i in 1:opt$n_sites)
    geno[[j]][i]=sample(0:2,size=1,prob=e[[i]][path[[j]][i]+1,])

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
for (j in 1:opt$n_ind) geno_lkl[[j]] = getLikes(geno[[j]],depth[j],opt$error,T,F)

# Print
lkl_out <- paste(opt$out,"glf.gz",sep=".")
if(!is.na(file.info(lkl_out)[,"size"])){
  warning("WARN: LKL file already exists. Skipping...");
} else {
  fh <- gzfile(lkl_out, "w", compression=9)
  for (i in 1:opt$n_sites) {
    for (j in 1:opt$n_ind)
      writeLines(as.character(geno_lkl[[j]][,i]), fh, "\t")
    write("", fh)
  }
  close(fh)
}
