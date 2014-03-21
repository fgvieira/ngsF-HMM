

############################ FUNCTIONS ############################
iter_plot <- function(true_path,true_geno,lkl,path,marg_prob){
  n = length(lkl)
  L = length(marg_prob[[1]])
  
  plots=c(1,1)
  if(n>1) plots=c(ceiling(n/2),2)
  par(mfrow=plots, mar=c(2,2,1,1))
  for (i in 1:n){
    ## Prob
    plot(1:L, marg_prob[[i]], ylim=c(0,1), xlab="", ylab="", xaxs="i", yaxs="i", type="n")
    # Plot marginal probs
    lines(1:L, marg_prob[[i]], col="green")
    # Plot most prob state path
    shade_areas(path[[i]], rgb(1,0,0,0.2))
    
    # Add title
    title(main=lkl[[i]], cex.main=0.5)
    
    # TRUE genotypes
    if(length(true_geno) != 0)
      points(true_geno[[i]]/2, pch=".")

    # TRUE state path
    if(length(true_path) != 0)
      shade_areas(true_path[[i]], rgb(0,0,1,0.2))
  }
}

shade_areas <- function(area, color){
  x <- rle(area)
  start <- (cumsum(x$lengths)-x$lengths)[x$values==1] + 1
  end <- cumsum(x$lengths)[x$values==1]
  
  if(length(start) != length(end)){
    cat("ERROR: start and end positions do not match!", fill=TRUE)
    #cat(x$values, fill=TRUE)
    #cat(x$lengths, fill=TRUE)
    cat(area, fill=TRUE)
    cat(start, fill=TRUE)
    cat(end, fill=TRUE)
    quit("no",-1)
  }
  
  # Check if there are any areas to shade
  if(length(start) == 0)
    return(-1);
  
  for (i in 1:length(start))
    polygon(c(start[i],start[i],end[i],end[i]),c(0,1,1,0), border=NA, col=color)
}

#####  Parse command-line arguments
library(optparse)
option_list <- list(make_option(c("-i", "--in_file"), action="store", type="character", default=NULL, help="Path to input file [%default]"),
                    make_option(c("-b", "--binary"), action="store_true", type="logical", default=FALSE, help="Is input in binary? [%default]"),
                    make_option(c("-n", "--n_ind"), action="store", type="integer", default=10, help="Number of individuals [%default]"),
                    make_option(c("-s", "--n_sites"), action="store", type="integer", default=1000, help="Number of sites [%default]"),
                    make_option(c("-g", "--geno"), action="store", type="character", default=NULL, help="Path to file with known/true genotypes (if available) [%default]"),
                    make_option(c("-p", "--path"), action="store", type="character", default=NULL, help="Path to file with known/true paths (if available) [%default]"),
                    make_option(c("--subset"), action="store", type="character", default=NULL, help="Iteration subset to plot (if available) [%default]"),
                    make_option(c("-o", "--out_prefix"), action="store", type="character", default=NULL, help="Output prefix [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

#opt$in_file = "ngsF-HMM/XXX.log.gz"; opt$n_ind=10; opt$n_sites=1000; opt$path="ngsF-HMM/sim.PI"; opt$geno="ngsF-HMM/sim.geno.gz"
#opt$in_file = "indica_allchr.1k.log.gz"; opt$n_ind=13; opt$n_sites=1000; opt$out_prefix="XXX"

############################ Parsing input arguments ############################

cat('### Input arguments', fill=TRUE)
cat('# Input file:', opt$in_file, fill=TRUE)
cat('# Is input file binary?', opt$binary, fill=TRUE)
cat('# Number indiv:', opt$n_ind, fill=TRUE)
cat('# Number sites:', opt$n_sites, fill=TRUE)
cat('# Known genotypes:', opt$geno, fill=TRUE)
cat('# Known path:', opt$path, fill=TRUE)
cat('# Subset:', opt$subset, fill=TRUE)
cat('# Out prefix:', opt$out_prefix, fill=TRUE)


if(!is.null(opt$geno) && file.exists(opt$geno)){
  cat("====> Reading geno file...", fill=TRUE)
  true_geno <- as.list(read.table(opt$geno))
}else
  true_geno <- list()
  
  
if(!is.null(opt$path) && file.exists(opt$path)){
  cat("====> Reading true path file...", fill=TRUE)
  true_path <- as.list(as.data.frame(t(read.fwf(opt$path, rep.int(1, opt$n_sites)))))
}else
  true_path <- list()


if(!is.null(opt$subset)){
  cat("====> Reading subset info...", fill=TRUE)
  subset <- as.numeric(strsplit(opt$subset,"[-:]")[[1]])
}else
  subset <- c()


cat("====> Opening input file stream...", fill=TRUE)
if(!is.null(opt$in_file) && file.exists(opt$in_file)){
  if( opt$binary && file.info(opt$in_file)$size %% 8*opt$n_ind*(1+2*opt$n_sites) != 0 ){
    cat("ERROR: corrupt input file!", fill=TRUE)
    quit("no",-1)
  }
  fh <- file(opt$in_file, ifelse(opt$binary,"rb","r"))
}else{
  cat("ERROR: cannot read input file!", fill=TRUE)
  quit("no",-1)
}



if(is.null(opt$out_prefix)){
  opt$out_prefix <- opt$in_file
  # Remove GZip extension
  opt$out_prefix <- sub(".gz", "", opt$out_prefix, fixed=TRUE)
  # Remove extension
  opt$out_prefix <- sub("\\.[^.]*$", "", opt$out_prefix, perl=TRUE)
}



############################ Plotting data ############################
iter <- 0
pdf(paste(opt$out_prefix,"pdf",sep="."), 2*log10(opt$n_sites))

while(iter <- iter + 1) {
  # Reading Lkl line
  if(opt$binary){
    lkl <- readBin(fh, double(), n=opt$n_ind)
  }else{
    line <- readLines(fh, 1)
    if(length(line) < 1) break
    lkl <- as.numeric(strsplit(line, "\t")[[1]][-1])
  }
  
  # If end of file
  if(length(lkl) < 1) break
  
  # Reading most prob path
  path <- list()
  if(opt$binary){
    for (j in 1:opt$n_ind)
      path[[j]] <- readBin(fh, integer(), n=opt$n_sites)
  }else
    path <- lapply(strsplit(readLines(fh, opt$n_ind), ""), as.numeric)
  
  # Reading post prob
  marg_prob <- list()
  if(opt$binary){
    for (j in 1:opt$n_ind)
      marg_prob[[j]] <- readBin(fh, double(), n=opt$n_sites)
  }else
    marg_prob <- lapply(strsplit(readLines(fh, opt$n_ind), "\t"), as.numeric)
  
  # If multiple of subset
  if(length(subset) == 1)
    if(iter %% subset != 0) next
  
  # If in the interval subset
  if(length(subset) > 1) {
    if(iter < subset[1]) next
    if(iter > subset[2]) break
  }
  
  # Plotting...
  cat("==> Plotting iteration:",iter , fill=TRUE)
  iter_plot(true_path,true_geno,lkl,path,marg_prob)
}

close(fh)
x <- dev.off()
