library(optparse)
library(plyr)

############################ FUNCTIONS ############################
iter_plot <- function(pos, path, marg_prob, plot_sites, true_path, true_geno, lkl, titles){
  n = length(lkl)
  L = length(path[[1]])

  plots=c(1,1)
  if(n>1) plots=c(ceiling(n/2),2)
  
  chrs <- unique(pos[,1])
  for (chr in chrs){
    par(mfrow=plots, mar=c(2,2,1,1))
    chr_pos <- which(pos==chr)

    for(i in 1:n){
      ## Plot
      plot(pos[chr_pos,2], path[[i]][chr_pos], ylim=c(0,1.01), xlab="", ylab="", xaxs="i", yaxs="i", type="n")

      if(plot_sites)
        points(pos[chr_pos,2], rep(0.0, L), pch="|", col='dodgerblue', cex=0.5)
      
      # Add title
      title(main=paste(titles[i],chr,lkl[[i]],sep=" / "), cex.main=0.5)
      
      # Plot marginal probs (GREEN line)
      if(length(marg_prob) != 0)
        lines(pos[chr_pos,2], marg_prob[[i]][chr_pos], col=rgb(0,1,0,0.5), lwd=0.1)
       
      # TRUE genotypes
      if(length(true_geno) != 0)
        points(pos[chr_pos,2], true_geno[[i]][chr_pos]/2, pch=".", col="cyan")
      
      # Shade most prob state path (BLUE)
      shade_areas(path[[i]][chr_pos], pos[chr_pos,2], rgb(0,0,1,0.2))
      
      # TRUE path (RED)
      if(length(true_path) != 0)
        shade_areas(true_path[[i]][chr_pos], pos[chr_pos,2], rgb(1,0,0,0.2), lims=c(0.25,0.75,0.75,0.25))
    }
  }
}

shade_areas <- function(area, pos, color, lims=c(0,1,1,0)){
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
    polygon(c(pos[start[i]],pos[start[i]],pos[end[i]],pos[end[i]]), lims, border=NA, col=color)
}

#####  Parse command-line arguments
option_list <- list(make_option(c("-i", "--in_file"), action="store", type="character", default=NULL, help="Path to input file [%default]"),
                    make_option(c("-b", "--binary"), action="store_true", type="logical", default=FALSE, help="Is input in binary? [%default]"),
                    make_option(c("-n", "--n_ind"), action="store", type="integer", default=10, help="Number of individuals [%default]"),
                    make_option(c("-s", "--n_sites"), action="store", type="integer", default=1000, help="Number of sites [%default]"),
                    make_option(c("-t", "--titles"), action="store", type="character", default=NULL, help="File with per-plot titles."),
                    make_option(c("--pos"), action="store", type="character", default=NULL, help="Path to file with site positions [%default]"),
                    make_option(c("-m", "--marg_prob"), action="store_true", type="logical", default=FALSE, help="Plot marginal probabilities? [%default]"),
                    make_option(c("-g", "--geno"), action="store", type="character", default=NULL, help="Path to file with known/true genotypes (if available) [%default]"),
                    make_option(c("-p", "--path"), action="store", type="character", default=NULL, help="Path to file with known/true paths (if available) [%default]"),
                    make_option(c("--plot_sites"), action="store_true", type="logical", default=FALSE, help="Plot position of SNPs? [%default]"),
                    make_option(c("--subset"), action="store", type="character", default=NULL, help="Iteration subset to plot (if available) [%default]"),
                    make_option(c("-w", "--width"), action="store", type="numeric", default=NULL, help="Each plot width [%default]"),
                    make_option(c("-o", "--out"), action="store", type="character", default=NULL, help="Output file [%default]"),
                    make_option(c("-q", "--quiet"), action="store_true", type="logical", default=FALSE, help="Print info to STDOUT? [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if(is.null(opt$width))
  opt$width <- ceiling(log10(opt$n_sites))

#opt$in_file="ceuMaf05.LD.ngsF-HMM.ibd"; opt$n_ind=7; opt$n_sites=253993; opt$pos="ceuMaf05.LD.pos.gz";opt$marg_prob=TRUE

############################ Parsing input arguments ############################

if(is.null(opt$out)){
  opt$out <- opt$in_file
  # Remove GZip extension
  opt$out <- sub(".gz", "", opt$out, fixed=TRUE)
  # Remove extension
  opt$out <- sub("\\.[^.]*$", "", opt$out, perl=TRUE)
  # Add PDF extension
  opt$out <- paste(opt$out, "pdf", sep=".")
}

if(!opt$quiet){
  cat('### Input arguments', fill=TRUE)
  cat('# Input file:', opt$in_file, fill=TRUE)
  cat('# Is input file binary?', opt$binary, fill=TRUE)
  cat('# Number indiv:', opt$n_ind, fill=TRUE)
  cat('# Number sites:', opt$n_sites, fill=TRUE)
  cat('# Titles:', opt$titles, fill=TRUE)
  cat('# Positions file:', opt$pos, fill=TRUE)
  cat('# Plot marginal probabilities?', opt$marg_prob, fill=TRUE)
  cat('# Known genotypes:', opt$geno, fill=TRUE)
  cat('# Known path:', opt$path, fill=TRUE)
  cat('# Plot site positions?', opt$plot_sites, fill=TRUE)
  cat('# Subset:', opt$subset, fill=TRUE)
  cat('# Width:', opt$width, fill=TRUE)
  cat('# Out file:', opt$out, fill=TRUE)
  cat('# Quiet run?', opt$quiet, fill=TRUE)
}



if(!is.null(opt$geno) && file.exists(opt$geno)){
  if(!opt$quiet)
    cat("====> Reading GENO file...", fill=TRUE)
  true_geno <- as.list(read.table(opt$geno))

  # If there are positions info (at least 3 cols and the second are positions), remove them!
  if( length(true_geno) >= 3 && max(true_geno[[2]]) > 2 ){
    true_geno[[2]] <- NULL
    true_geno[[1]] <- NULL
  }
  
  if(opt$n_ind != length(true_geno) ||
     opt$n_sites != length(true_geno[[1]])){
    cat("ERROR: number of indiv/sites and GENO file do not match!", fill=TRUE)
    quit("no",-1)
  }
}else
  true_geno <- list()


  
if(!is.null(opt$titles)){
  if(!opt$quiet)
    cat("====> Reading TITLES file...", fill=TRUE)
  titles <- apply(as.matrix(read.table(opt$titles, sep="\t")), 1, paste, collapse=" / ")

  if(opt$n_ind != length(titles)){
    cat("ERROR: number of indiv and TITLES file do not match!", fill=TRUE)
    quit("no",-1)
  }
}else
  titles <- matrix()



if(!is.null(opt$path) && file.exists(opt$path)){
  if(!opt$quiet)
    cat("====> Reading PATH file...", fill=TRUE)
  true_path <- lapply(strsplit(readLines(opt$path, opt$n_ind), ""), as.numeric)
  
  if(opt$n_ind != length(true_path) ||
     opt$n_sites != length(true_path[[1]])){
    cat("ERROR: number of indiv/sites and PATH file do not match!", fill=TRUE)
    quit("no",-1)
  }
}else
  true_path <- list()



if(!is.null(opt$subset)){
  if(!opt$quiet)
    cat("====> Reading subset info...", fill=TRUE)
  subset <- as.numeric(strsplit(opt$subset,"[-:/,]")[[1]])
}else
  subset <- c()



if(!is.null(opt$pos)){
  if(file.exists(opt$pos)){
    if(!opt$quiet)
      cat("====> Reading POSITIONS file...", fill=TRUE)
    pos <- read.table(opt$pos, header=FALSE, stringsAsFactors=FALSE, check.names=FALSE)
    pos[,2] <- pos[,2]/1e6
    
    if(opt$n_sites != nrow(pos)){
      cat("ERROR: number of sites and positions file do not match!", fill=TRUE)
      quit("no",-1)
    }
  }else{
    cat("ERROR: cannot open positions file...", fill=TRUE)
    quit("no",-1)
  }
}
 

  
if(!opt$quiet)
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


############################ Plotting data ############################
if(opt$n_ind == 1){
  pdf(opt$out, opt$width, 2)
}else{
  pdf(opt$out, 2*opt$width, 2*opt$n_ind/2)
}

iter <- -1
while(TRUE) {
  iter <- iter + 1
  
  # Reading Lkl line
  if(opt$binary){
    lkl <- readBin(fh, double(), n=opt$n_ind)
  }else{
    line <- readLines(fh, 1)
    if(length(line) < 1) break
    lkl <- strsplit(line, "\t")[[1]][-1]
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
  if(!opt$marg_prob)
    marg_prob <- list();
  
  
  if(length(subset) == 1){
    if(iter < subset[1]) next
    if(iter > subset[1]) break
  }else if(length(subset) == 2){
    if(is.na(subset[1])){
      # If multiple of subset
      if(iter != 1 && iter %% subset[2] != 0) next
    }else{
      # If in the interval subset
      if(iter < subset[1]) next
      if(iter > subset[2]) break
    }
  }
  

  # Plotting...
  if(!opt$quiet)
    cat("> Plotting iter", iter, "...", fill=TRUE)
  iter_plot(pos, path, marg_prob, opt$plot_sites, true_path, true_geno, lkl, titles)
}

close(fh)
x <- dev.off()
