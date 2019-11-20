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
  
  return(list(log(res),depth))
}

#####  Parse command-line arguments
library(optparse)
option_list <- list(make_option(c("-n", "--n_ind"), action="store", type="integer", default=10, help="Number of individuals [%default]"),
                    make_option(c("-s", "--n_sites"), action="store", type="integer", default=1000, help="Number of independent sites [%default]"),
                    make_option(c("-f", "--freq"), action="store", type="character", default="0.1", help="Allele frequencies [%default]"),
                    make_option(c("-x", "--site_pos"), action="store", type="character", default="1", help="Distance between sites [%default]"),
                    make_option(c("-F", "--indF"), action="store", type="character", default="0", help="Per-individual inbreeding coefficients [%default]"),
                    make_option(c("-a", "--alpha"), action="store", type="character", default="0.01", help="Transition parameter [%default]"),
                    make_option(c("-d", "--depth"), action="store", type="character", default="5", help="Sequencing depth [%default]"),
                    make_option(c("-e", "--error"), action="store", type="numeric", default=0.01, help="Error rate [%default]"),
                    make_option(c("--seed"), action="store", type="integer", default=runif(1,1,1e6), help="Seed for random number generator [random]"),
                    make_option(c("-o", "--out"), action="store", type="character", default="sim", help="Output prefix [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

############################ Parsing input arguments ############################
# Print parameters
cat('# Number of individuals:', opt$n_ind, fill=TRUE)
cat('# Number of sites:', opt$n_sites, fill=TRUE)
cat('# Site freqs:', opt$freq, fill=TRUE)
cat('# Site positions:', opt$site_pos, fill=TRUE)
cat('# indF:', opt$indF, fill=TRUE)
cat('# Transition parameter:', opt$alpha, fill=TRUE)
cat('# Depth:', opt$depth, fill=TRUE)
cat('# Seed:', opt$seed, fill=TRUE)
set.seed(opt$seed)
cat('# Output prefix:', opt$out, fill=TRUE)



### indF
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
  indF <- rep(as.numeric(opt$indF),opt$n_ind)
}



### freq
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
  freq <- rep(as.numeric(opt$freq),opt$n_ind)
}



### alpha
cat("====> Parsing per-individual TRANSITION parameter...",fill=TRUE)
if(file.exists(opt$alpha)){
  cat("==> Reading values from file:",opt$alpha,fill=TRUE)
  alpha <- as.numeric(readLines(opt$alpha))
  if(opt$n_ind != length(alpha)){
    cat("ERROR: number of individuals and TRANSITION file do not match", fill=TRUE)
    quit("no",-1)
  }
}else if(opt$alpha == "r"){
  cat("==> Using random values",fill=TRUE)
  alpha <- runif(opt$n_ind)
}else{
  cat("==> Setting to fixed value:",as.numeric(opt$alpha),fill=TRUE)
  alpha <- rep(as.numeric(opt$alpha),opt$n_ind)
}



### depth
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
  depth <- rep(as.numeric(opt$depth),opt$n_ind)
}



### Parse sites' positions
cat("====> Parsing sites positions...",fill=TRUE)
if(file.exists(opt$site_pos)){
  cat("==> Reading values from file:",opt$site_pos,fill=TRUE)
  site_pos <- read.table(opt$site_pos)
  if(opt$n_sites != nrow(site_pos)){
    cat("ERROR: number of sites and DIST file do not match", fill=TRUE)
    quit("no",-1)
  }

  pos_dist <- site_pos[,2]
  for (i in 2:opt$n_sites)
    if(site_pos[i,1] != site_pos[i-1,1]){
      pos_dist[i] <- +Inf
    }else{
      pos_dist[i] <- pos_dist[i] - pos_dist[i-1]
    }
}else if(opt$site_pos == "r"){
  avg_dist <- 1e5 # Avg dist between sampled independent SNPs
  cat("==> Setting to normally distributed with mean:",avg_dist,fill=TRUE)
  pos_dist <- as.integer(rnorm(opt$n_sites,mean=avg_dist,sd=avg_dist/3))
  pos_dist[pos_dist < 1] = 1
}else{
  cat("==> Setting to fixed value:",as.numeric(opt$site_pos),fill=TRUE)
  pos_dist <- rep(as.numeric(opt$site_pos),opt$n_sites)
}



############################ Generate data ############################

############ Generating TRUE (hidden) path
# Distances are converted to Mb (divided by 1e6)
cat("====> Generating true path...",fill=TRUE)
path = list();
for (i in 1:opt$n_ind)
  path[[i]] = get_IBD(pos_dist/1e6, indF[i], alpha[i])


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
true_depth <- geno
for (i in 1:opt$n_ind){
  x <- getLikes(geno[[i]],depth[i],opt$error,T,F)
  geno_lkl[[i]] <- x[[1]]
  true_depth[[i]] <- x[[2]]
 }

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



############ Print sites' positions
cat("====> Printing sites positions...",fill=TRUE)
dist_out <- paste(opt$out,"pos.gz",sep=".")
if(!is.na(file.info(dist_out)[,"size"])){
  warning("WARN: DIST file already exists. Skipping...");
} else {
  x <- as.matrix(as.data.frame(true_depth))
  fh <- gzfile(dist_out, "w", compression=9)
  write.table(paste("chrSIM",cumsum(as.numeric(pos_dist)),rowSums(x),apply(x,1,paste,collapse=","),sep="\t"), fh, quote=FALSE, sep="\n", row.names=FALSE, col.names=FALSE)
  close(fh)
}
