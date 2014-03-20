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

logsum <- function(x) {
  uniq_x = unique(x)
  # Check if all elements -Inf
  if(length(uniq_x)==1 && is.infinite(uniq_x)) return(-Inf)
  
  #function does: log(exp(a)+exp(b)) while protecting for underflow
  return(log(sum(exp(x - max(x)))) + max(x))
}

post_prob <- function(lkl, prior) {
  pp <- lkl#+prior
  pp <- pp-logsum(pp)  
  
  return(pp)
}

callGeno <- function(pp){
  ppN = c(-Inf,-Inf,-Inf)
  ppN[which.max(pp)] <- log(1)

  return(ppN)
}

iter_EM <- function(geno_lkl, a, e, v, callG) {
  n <- length(geno_lkl)
  L <- ncol(geno_lkl[[1]])
  inbreed=0.1
  
  #Define variables
  m <- matrix(nrow=n_states,ncol=L,dimnames=list(c("nIBD","IBD"),c(1:L)))
  V = F = F0 = B = B0 = marginal <- list()
  for (j in 1:n){
    V[[j]] <- m
    F[[j]] <- m
    F0[[j]] <- m[,1]
    F0[[j]][1] = log(1-inbreed)
    F0[[j]][2] = log(inbreed)
    B[[j]] <- m
    B0[[j]] <- m[,1]
    marginal[[j]] <- m
  }
  Px <- numeric(n)
  new_v <- v
  new_a <- a
  new_freq <- freq
  new_e <- e
  
  cat("==> Forward Recursion", fill=TRUE)
  for (j in 1:n)
    for (i in 1:L){
      pp <- post_prob(geno_lkl[[j]][,i],e[[i]][v[[j]][i]+1,])
      if(callG) pp <- callGeno(pp)
      for (l in 1:2){
        if(i == 1)
          F[[j]][l,i] <- logsum(c(F0[[j]][1]+a[[j]][1,l],
                                  F0[[j]][2]+a[[j]][2,l]))+logsum(e[[i]][l,]+pp)
        else
          # logsum(k==1,k==2)
          F[[j]][l,i] <- logsum(c(F[[j]][1,i-1]+a[[j]][1,l],
                                  F[[j]][2,i-1]+a[[j]][2,l]))+logsum(e[[i]][l,]+pp)
      }
    }
  
  cat("==> Backward Recursion", fill=TRUE)
  for (j in 1:n)
    for (i in seq(L,1,by=-1))
      for (k in 1:2){
        if(i == L)
          B[[j]][,i] <- log(1)
        else{
          pp <- post_prob(geno_lkl[[j]][,i+1],e[[i+1]][v[[j]][i+1]+1,])
          if(callG) pp = callGeno(pp)
          # logsum(l==1,l==2)
          B[[j]][k,i] <- logsum(c(a[[j]][k,1]+logsum(e[[i+1]][1,]+pp)+B[[j]][1,i+1],
                                  a[[j]][k,2]+logsum(e[[i+1]][2,]+pp)+B[[j]][2,i+1]))
        }
      }
  cat("> Backward termination", fill=TRUE)
  i=0
  for (j in 1:n){
    pp <- post_prob(geno_lkl[[j]][,i+1],e[[i+1]][v[[j]][i+1]+1,])
    if(callG) pp = callGeno(pp)
    for (k in 1:2)
      B0[[j]][k] <- F0[[j]][k]+logsum(c(a[[j]][k,1]+logsum(e[[i+1]][1,]+pp)+B[[j]][1,i+1],
                                        a[[j]][k,2]+logsum(e[[i+1]][2,]+pp)+B[[j]][2,i+1]))
  }
  
  ### TEST ######
  for (j in 1:n)
    if(abs(logsum(F[[j]][,L])-logsum(B0[[j]])) > 1e-10)
      cat("=======> ",logsum(F[[j]][,L])-logsum(B0[[j]]),fill=T)
  
  cat("==> Marginal probabilities", fill=TRUE)
  for (j in 1:n){
    Px[j] <- logsum(F[[j]][,L])
    marginal[[j]] <- B[[j]]+F[[j]]-Px[j]
  }
  
  if(!opt$path_fixed){
    cat("==> Update most probable states", fill=TRUE)
    if(opt$viterbi){
      cat("> Viterbi algorithm", fill=TRUE)
      for (j in 1:n)
        for (i in 1:L){
          pp <- post_prob(geno_lkl[[j]][,i],e[[i]][v[[j]][i]+1,])
          if(callG) pp <- callGeno(pp)
          for (l in 1:2){
            if(i==1)
              V[[j]][l,i] <- max(c(F0[[j]][1]+a[[j]][1,l],
                                   F0[[j]][2]+a[[j]][2,l]))+logsum(e[[i]][l,]+pp)
            else
              # max(k==1,k==2)
              V[[j]][l,i]=max(V[[j]][1,i-1]+a[[j]][1,l],
                              V[[j]][2,i-1]+a[[j]][2,l])+logsum(e[[i]][l,]+pp)
          }
        }
      cat("> Back-tracking", fill=TRUE)
      for (j in 1:n)
        for (i in 1:L)
          new_v[[j]][i] <- which.max(V[[j]][,i])-1
    }else{
      cat("> Posterior decoding", fill=TRUE)
      for (j in 1:n)
        for (i in 1:L)
          new_v[[j]][i] <- which.max(marginal[[j]][,i])-1
    }
  }else
    cat("==> Most probable states not estimated!", fill=TRUE)
  
  if(!opt$trans_fixed){
    new_a <- list()
    for (j in 1:n)
      new_a[[j]] <- matrix(-Inf,nrow=n_states,ncol=n_states,dimnames=list(c("nIBD","IBD"),c("nIBD","IBD")));
    
    cat("==> Update transition probabilities", fill=TRUE)
    for (j in 1:n)
      for (k in 1:2){
        Pk <- logsum(F[[j]][k,-L]+B[[j]][k,-L])
        for (l in 1:2)
          for (i in 1:(L-1)){
            pp <- post_prob(geno_lkl[[j]][,i+1],e[[i+1]][v[[j]][i+1]+1,])
            if(callG) pp <- callGeno(pp)
            new_a[[j]][k,l] <- logsum(c(new_a[[j]][k,l], F[[j]][k,i]+a[[j]][k,l]+logsum(e[[i+1]][l,]+pp)+B[[j]][l,i+1]-Pk))
          }
      }
  }else
    cat("==> Transition probabilities not estimated!", fill=TRUE)
  
  if(!opt$freq_fixed){
    cat("==> Update allele freqs", fill=TRUE)
    for (i in 1:L) {
      num = den <- 0
      for (j in 1:n) {
        pp <- post_prob(geno_lkl[[j]][,i],e[[i]][v[[j]][i]+1,])
        if(callG) pp <- callGeno(pp)
        num <- num + exp(pp[2]) + exp(pp[3])*(2-exp(marginal[[j]][2,i]))
        den <- den + 2*exp(pp[2]) + exp(logsum(c(pp[1],pp[3])))*(2-exp(marginal[[j]][2,i]))
      }
      new_freq[i] <- num/den
    }
  }else
    cat("==> Allele freqs not estimated!", fill=TRUE)
  
  if(opt$e_BW){
    cat("==> Update emission probabilities with TRUE Baum-Welch", fill=TRUE) #untested
    tmp_e <- matrix(-Inf,nrow=n_states,ncol=3,dimnames=list(c("nIBD","IBD"),c("AA","Aa","aa")))
    
    for (j in 1:n)
      for (i in 1:(L-1)){
        pp <- post_prob(geno_lkl[[j]][,i],e[[1]][v[[j]][i]+1,])        
        if(callG) pp <- callGeno(pp)
        for (k in 1:2){
          Pk <- logsum(F[[j]][k,-L]+B[[j]][k,-L])
          for (x in 0:2)
            tmp_e[k,x+1] <- logsum(c(tmp_e[k,x+1], F[[j]][k,i]+B[[j]][k,i]+pp[x+1]-Pk))
        }
      }
    
    for (i in 1:(L-1))
      new_e[[i]] <- tmp_e/n
  }else{
    cat("==> Update emission probs based on allele freqs", fill=TRUE)
    for (i in 1:L)
      for (k in 1:2){
        f=new_freq[i]
        new_e[[i]][k,1] <- log((1-f)^2+(1-f)*f*(k-1))
        new_e[[i]][k,2] <- log(2*(1-f)*f-2*(1-f)*f*(k-1))
        new_e[[i]][k,3] <- log(f^2+(1-f)*f*(k-1))
      }
  }
  
  cat("==> Update inbreeding coef", fill=TRUE)
  inbreed <- numeric(n)
  for (j in 1:n)
    inbreed[j] = exp(logsum(marginal[[j]][2,]))/L
  
  return(list(Px=Px,a=new_a,e=new_e,marginal=marginal,v=new_v,freq=new_freq,inbreed=inbreed))
}

#####  Parse command-line arguments
library(optparse)
option_list <- list(make_option(c("-g", "--geno"), action="store", type="character", default=NULL, help="Path to input file; it can be either genotypes {0,1,2} or genotype likelihoods [%default]"),
                    make_option(c("--is_lkl"), action="store_true", type="logical", default=FALSE, help="Is input genotype likelihoods? [%default]"),
                    make_option(c("--is_log"), action="store_true", type="logical", default=FALSE, help="Is input in log-scale? [%default]"),
                    make_option(c("--call_geno"), action="store_true", type="logical", default=FALSE, help="Call genotytpes? [%default]"),
                    make_option(c("-f", "--freq"), action="store", type="character", default="r", help="Allele frequencies initial vallues [%default]"),
                    make_option(c("--freq_fixed"), action="store_true", type="logical", default=FALSE, help="Are freqs values to be fixed (not estimated)? [%default]"),
                    make_option(c("-t", "--trans"), action="store", type="character", default="r", help="Transition probabilities initial values [%default]"),
                    make_option(c("--trans_fixed"), action="store_true", type="logical", default=FALSE, help="Are transition values to be fixed (not estimated)? [%default]"),
                    make_option(c("-p", "--path"), action="store", type="character", default="r", help="Hidden path initial values [%default]"),
                    make_option(c("--path_fixed"), action="store_true", type="logical", default=FALSE, help="Are hidden path values to be fixed (not estimated)? [%default]"),
                    make_option(c("--viterbi"), action="store_true", type="logical", default=FALSE, help="Use states from Viterbi? [%default]"),
                    make_option(c("--e_BW"), action="store_true", type="logical", default=FALSE, help="Use emission from normal Baum-Welch? [%default]"),
                    make_option(c("--min_iters"), action="store", type="integer", default=10, help="Minimum number of iterations [%default]"),
                    make_option(c("--max_iters"), action="store", type="integer", default=100, help="Maximum number of iterations [%default]"),
                    make_option(c("-o", "--out"), action="store", type="character", default="output", help="Output prefix [%default]"),
                    make_option(c("-b", "--binary"), action="store_true", type="logical", default=FALSE, help="Is input in binary? [%default]"),
                    make_option(c("-l", "--log"), action="store_true", type="logical", default=FALSE, help="Output log? [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

#opt$geno="sim.geno_lkl"; opt$is_lkl=T; opt$is_log=T; opt$freq="0.1"; opt$trans="0.01-0.01"; opt$path="0"; opt$freq_fixed=F; opt$trans_fixed=F; opt$path_fixed=F; opt$log=T; opt$max_iters=50; opt$viterbi=T; opt$e_BW=F
#opt$geno="sim.geno"; opt$is_lkl=F; opt$is_log=F; opt$freq="0.1"; opt$trans="0.01-0.01"; opt$path="sim.PI"; opt$freq_fixed=F; opt$trans_fixed=F; opt$path_fixed=F; opt$log=T; opt$max_iters=50; opt$viterbi=T; opt$e_BW=F
#opt$geno="sim.geno"; opt$is_lkl=F; opt$is_log=F; opt$freq="0.01"; opt$trans="0.1-0.1"; opt$path="0"; opt$freq_fixed=F; opt$trans_fixed=F; opt$path_fixed=F; opt$log=T; opt$max_iters=50; opt$viterbi=T; opt$e_BW=F

############################ Parsing input arguments ############################
n_states=2

### Read input file
cat("====> Parsing genotype info...",fill=TRUE)
if(!is.null(opt$geno) && file.exists(opt$geno)){
  cat("==> Reading from file:",opt$geno,fill=TRUE)
  geno_lkl <- list();
  
  if(opt$is_lkl){
    cat("==> Genotype likelihoods in", ifelse(opt$is_log,"log-scale","normal-scale"),fill=TRUE)
    
    if(opt$call_geno)
      cat("==> Calling genotype with highest Lkl",fill=TRUE)
    
    fh <- file(opt$geno, "r")
    buf <- strsplit(readLines(fh), "\t")
    close(fh)
    
    n_sites <- length(buf)
    n_ind <- length(buf[[1]])/3
    for (j in 1:n_ind) geno_lkl[[j]] <- matrix(0,3,n_sites)
    
    for (i in 1:n_sites)
      for (j in 1:n_ind){
        gl <- as.numeric(buf[[i]][((j-1)*3+1):((j-1)*3+3)])
        if(!opt$is_log) gl <- log(gl)
        geno_lkl[[j]][,i] <- gl 
      }
    
    rm(buf)
  }else{
    cat("==> Called genotypes",fill=TRUE)
    geno <- as.list(read.table(opt$geno))
    
    n_ind <- length(geno)
    n_sites <- length(geno[[1]])
    for (j in 1:n_ind) geno_lkl[[j]] <- matrix(0,3,n_sites)
    
    for (i in 1:n_sites)
      for (j in 1:n_ind){
        gl <- log(numeric(3))
        gl[geno[[j]][i]+1] <- log(1)
        geno_lkl[[j]][,i] <- gl
      }
    opt$is_log <- TRUE
  }
}else{
  cat("ERROR: cannot read input file!",fill=TRUE)
  quit("no",-1)
}


cat("=> Number individuals:", n_ind,fill=TRUE)
cat("=> Number sites:", n_sites,fill=TRUE)


### freq
freq = c()
cat("====> Parsing initial allele frequencies...",fill=TRUE)
if(file.exists(opt$freq)){
  cat("==> Reading values from file:",opt$freq,fill=TRUE)
  freq <- as.numeric(readLines(opt$freq))
  if(n_sites != length(freq)){
    cat("ERROR: number of sites and freq file do not match", fill=TRUE)
    quit("no",-1)
  }
}else if(opt$freq == "r"){
  cat("==> Using random values",fill=TRUE)
  freq <- runif(n_sites)
}else if(!is.na(as.integer(opt$freq))){
  cat("==> Setting to:",as.numeric(opt$freq),fill=TRUE)
  for (i in 1:n_sites) freq[i] <- as.numeric(opt$freq)
}else{
  cat("ERROR: invalid frequency option", fill=TRUE)
  quit("no",-1)
}
cat("==> Values will",ifelse(opt$freq_fixed,"NOT",""),"be estimated!",fill=TRUE)


### Transicion probs
alpha = beta = c()
cat("====> Parsing initial transition probs...",fill=TRUE)
if(file.exists(opt$trans)){
  cat("==> Reading values from file:",opt$trans,fill=TRUE)
  alpha_beta <- strsplit(readLines(opt$trans),"[-,\t]")
  if(n_ind != length(alpha_beta)){
    cat("ERROR: number of individuals and transitions file do not match", fill=TRUE)
    quit("no",-1)
  }
  
  for (j in 1:n_ind){
    alpha[j] <- as.numeric(alpha_beta[[j]][1])
    beta[j] <- as.numeric(alpha_beta[[j]][2])
  }
}else if(opt$trans == "r"){
  cat("==> Using random values",fill=TRUE)
  alpha <- runif(n_ind)
  beta <- runif(n_ind)
}else {
  for (j in 1:n_ind){
    alpha[j] <- as.numeric(strsplit(opt$trans,"[-, \t]")[[1]][1])
    beta[j] <- as.numeric(strsplit(opt$trans,"[-, \t]")[[1]][2])
  }
  cat("==> Setting to:", alpha[1], beta[1], fill=TRUE)
}
cat("==> Values will",ifelse(opt$trans_fixed,"NOT",""),"be estimated!",fill=TRUE)


### Most probable path
v <- list()
cat("====> Parsing initial path...",fill=TRUE)
if(file.exists(opt$path)){
  cat("==> Reading values from file:",opt$path,fill=TRUE)
  v <- as.list(as.data.frame(t(read.fwf(opt$path, rep.int(1,n_sites)))))
}else if(opt$path == "r"){
  cat("==> Using random states",fill=TRUE)
  for (j in 1:n_ind) v[[j]] <- sample(c(0,1),n_sites,TRUE)
}else if(!is.na(as.integer(opt$path))){
  cat("==> Setting to fixed state:",as.integer(opt$path),fill=TRUE)
  for (j in 1:n_ind) v[[j]] <- rep(as.integer(opt$path), n_sites)
}else{
  cat("ERROR: invalid path option", fill=TRUE)
  quit("no",-1)
}
cat("==> Values will",ifelse(opt$path_fixed,"NOT",""),"be estimated!",fill=TRUE)



############################ Set starting values ############################
### Transition probs
a <- list()
for (j in 1:n_ind) a[[j]] <- matrix(c(1-alpha[[j]],alpha[[j]],beta[[j]],1-beta[[j]]),ncol=n_states,byrow=T,dimnames=list(c("nIBD","IBD"),c("nIBD","IBD")))
a <- lapply(a,log)

### Emission probs
e <- list()
for (s in 1:n_sites) {
  e[[s]] <- matrix(nrow=n_states,ncol=3,dimnames=list(c("nIBD","IBD"),c("AA","Aa","aa")))
  for (k in 1:2){
    f <- freq[s]
    e[[s]][k,1] <- log((1-f)^2+(1-f)*f*(k-1))
    e[[s]][k,2] <- log(2*(1-f)*f-2*(1-f)*f*(k-1))
    e[[s]][k,3] <- log(f^2+(1-f)*f*(k-1))
  }
}

params <- list()
params$a <- a
params$e <- e
params$v <- v
params$Px <- numeric(n_ind)

### main loop
if(opt$log)
  fh <- gzfile(paste(opt$out,"log.gz",sep="."), ifelse(opt$binary,"wb","w"))

cat("====> Inferring parameters...",fill=TRUE)
iter <- 1
old_Px <- rep(-50, n_ind)
while(iter <= opt$max_iters && 
      (iter <= opt$min_iters || max(params$Px - old_Px) > 1)){
  old_Px <- params$Px

  cat("====> Iter:", iter, fill=TRUE)
  params <- iter_EM(geno_lkl,params$a,params$e,params$v,opt$call_geno)

  # Print Lkl to output
  write(params$Px - old_Px, ncolumns=n_ind+1, "", sep="\t")

  if(opt$log){
    # Print Lkl
    if(opt$binary){
      writeBin(params$Px, fh)
    }else
      write(c("//", params$Px), ncolumns=n_ind+1, fh, sep="\t")
    
    # Write most probable path (Viterbi)
    write.table(t(as.data.frame(params$v)), fh, quote=FALSE, sep="", row.names=FALSE, col.names=FALSE)
    if(opt$binary)
	for (j in 1:n_ind)
	    writeBin(as.integer(params$v[[j]]), fh, size=1)

    # Print marginal probs
    for (j in 1:n_ind)
      if(opt$binary){
        writeBin(exp(params$marginal[[j]][2,]), fh)
      }else
        write(exp(params$marginal[[j]][2,]), fh, ncolumns=n_sites, sep="\t")
  }

  iter <- iter + 1
}
if(opt$log)
  close(fh)


# Open filehandle
fh <- file(paste(opt$out,"out",sep="."), "w")
# F
write(params$inbreed, fh, ncolumns=n_ind, sep="\t")
# alpha/beta
for (j in 1:n_ind)
  write(c(exp(params$a[[j]][1,2]), exp(params$a[[j]][2,1])), fh, ncolumns=2, sep="\t")
# Freqs
write(params$freq, fh, ncolumns=1)
# Close filehandle
close(fh)
