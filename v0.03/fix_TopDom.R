#' Identify Topological Domains from a Hi-C Contact Matrix
#' 
#' @param data A TopDomData object, or the pathname to a normalized
#' Hi-C contact matrix file as read by [readHiC()], that specify N bins.
#' 
#' @param window.size The number of bins to extend (as a non-negative integer).
#' Recommended range is in {5, ..., 20}.
#' 
#' @param outFile (optional) The filename without extension of the three
#' result files optionally produced. See details below.
#' 
#' @param statFilter (logical) Specifies whether non-significant
#' topological-domain boundaries should be dropped or not.
#' 
#' @param ... Additional arguments passed to [readHiC()].
#' 
#' @param debug If `TRUE`, debug output is produced.
#'
#' @return A named list of class `TopDom` with data.frame elements
#' `binSignal`, `domain`, and `bed`.
#' * The `binSignal` data frame (N-by-7) holds mean contact frequency,
#'   local extreme, and p-value for every bin. The first four columns
#'   represent basic bin information given by matrix file, such as
#'   bin id (`id`), chromosome(`chr`), start coordinate (`from.coord`),
#'   and end coordinate (`to.coord`) for each bin.
#'   The last three columns (`local.ext`, `mean.cf`, and `p-value`) represent
#'   computed values by the TopDom algorithm.
#'   The columns are:
#'   - `id`: Bin ID
#'   - `chr`: Chromosome
#'   - `from.coord`: Start coordinate of bin
#'   - `to.coord`: End coordinate of bin
#'   - `local.ext`:
#'      + `-1`: Local minima.
#'      + `-0.5`: Gap region.
#'      + `0`: General bin.
#'      + `1`: Local maxima.
#'   - `mean.cf`: Average of contact frequencies between lower and upper
#'     regions for bin _i = 1,2,...,N_.
#'   - `p-value`: Computed p-value by Wilcox rank sum test.
#'     See Shin et al. (2016) for more details.
#' 
#' * The `domain` data frame (D-by-7):
#'   Every bin is categorized by basic building block, such as gap, domain,
#'   or boundary.
#'   Each row indicates a basic building block.
#'   The first five columns include the basic information about the block,
#'   'tag' column indicates the class of the building block.
#'   - `id`: Identifier of block
#'   - `chr`: Chromosome
#'   - `from.id`: Start bin index of the block
#'   - `from.coord`: Start coordinate of the block
#'   - `to.id`: End bin index of the block
#'   - `to.coord`: End coordinate of the block
#'   - `tag`: Categorized name of the block. Three possible blocks exists:
#'     + `gap`
#'     + `domain`
#'     + `boundary`
#'   - `size`: size of the block
#'
#' * The `bed` data frame (D-by-4) is a representation of the `domain`
#'   data frame in the
#'   [BED file format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1).
#'   It has four columns:
#'   - `chrom`: The name of the chromosome.
#'   - `chromStart`: The starting position of the feature in the chromosome.
#'      The first base in a chromosome is numbered 0.
#'   - `chromEnd`: The ending position of the feature in the chromosome.
#'      The `chromEnd` base is _not_ included in the feature. For example,
#'      the first 100 bases of a chromosome are defined as `chromStart=0`,
#'      `chromEnd=100`, and span the bases numbered 0-99.
#'   - `name`: Defines the name of the BED line. This label is displayed to
#'      the left of the BED line in the
#'      [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgGateway)
#'      window when the track is open to full display mode or directly to
#'      the left of the item in pack mode.
#' 
#' If argument `outFile` is non-`NULL`, then the three elements (`binSignal`,
#' `domain`, and `bed`) returned are also written to tab-delimited files
#' with file names \file{<outFile>.binSignal}, \file{<outFile>.domain}, and
#' \file{<outFile>.bed}, respectively.  None of the files have row names,
#' and all but the BED file have column names.
#'
#' @section Windows size:
#' The `window.size` parameter is by design the only tuning parameter in the
#' TopDom method and affects the amount of smoothing applied when calculating
#' the TopDom bin signals.  The binning window extends symmetrically downstream
#' and upstream from the bin such that the bin signal is the average
#' `window.size^2` contact frequencies.
#' For details, see Equation (1) and Figure 1 in Shin et al. (2016).
#' Typically, the number of identified TDs decreases while their average
#' lengths increase as this window-size parameter increases (Figure 2).
#' The default is `window.size = 5` (bins), which is motivated as:
#' "Considering the previously reported minimum TD size (approx. 200 kb)
#' (Dixon et al., 2012) and our bin size of 40 kb, _w_\[indow.size\] = 5 is a
#' reasonable setting" (Shin et al., 2016).
#'
#' @example incl/TopDom.R
#' 
#' @author Hanjun Shin, Harris Lazaris, and Gangqing Hu.
#' \R package, help, and code refactoring by Henrik Bengtsson.
#' \Edit by Julie Segueni to fix integer overflow bug.
#'
#' @references
#'
#' * Shin et al.,
#'   TopDom: an efficient and deterministic method for identifying
#'   topological domains in genomes,
#'   _Nucleic Acids Research_, 44(7): e70, April 2016.
#'   DOI: 10.1093/nar/gkv1505,
#'   PMCID: [PMC4838359](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4838359/),
#'   PMID: [26704975](https://pubmed.ncbi.nlm.nih.gov/26704975/)
#'
#' * Shin et al., \R script \file{TopDom_v0.0.2.R}, 2017 (originally from
#'   \code{http://zhoulab.usc.edu/TopDom/};
#'   later available on \url{https://github.com/jasminezhoulab/TopDom} via
#'   \url{https://zhoulab.dgsom.ucla.edu/pages/software})
#'
#' * Shin et al., TopDom Manual, 2016-07-08 (original from
#'   \code{http://zhoulab.usc.edu/TopDom/TopDom\%20Manual_v0.0.2.pdf};
#'   later available on \url{https://github.com/jasminezhoulab/TopDom} via
#'   \url{https://zhoulab.dgsom.ucla.edu/pages/software})
#'
#' * Hanjun Shin, Understanding the 3D genome organization in topological
#'   domain level, Doctor of Philosophy Dissertation,
#'   University of Southern California, March 2017,
#'   \url{http://digitallibrary.usc.edu/cdm/ref/collection/p15799coll40/id/347735}
#'
#' * Dixon JR, Selvaraj S, Yue F, Kim A, et al. Topological domains in
#'   mammalian genomes identified by analysis of chromatin interactions.
#'   _Nature_; 485(7398):376-80, April 2012.
#'   DOI: 10.1038/nature11082,
#'   PMCID: [PMC3356448](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3356448/),
#'   PMID: 22495300.
#'
#' @importFrom utils read.table write.table
#' @export
TopDom_0.0.3 <- local({
TopDom <- function( matrix.file, window.size, outFile=NULL, statFilter=T)
{
  if (inherits(matrix.file, "TopDomData")) {
    bins <- matrix.file$bins
    matrix.data <- matrix.file$counts
    n_bins <- as.double(nrow(bins))
    mean.cf <- rep(0, times = n_bins)
    pvalue <- rep(1.0, times = n_bins)
    local.ext <- rep(-0.5, times = n_bins)
  } else {
    print("#########################################################################")
    print("Step 0 : File Read ")
    print("#########################################################################")
    window.size = as.numeric(window.size)
    matdf <- read.table(matrix.file, header=F)
    
    if( ncol(matdf) - nrow(matdf) == 3) {
      colnames(matdf) <- c("chr", "from.coord", "to.coord")
    } else if( ncol(matdf) - nrow(matdf) ==4 ) {
      colnames(matdf) <- c("id", "chr", "from.coord", "to.coord")
    } else {
      print("Unknwon Type of matrix file")
      return(0)
    }
    n_bins = nrow(matdf)
    mean.cf <- rep(0, n_bins)
    pvalue <- rep(1, n_bins)
    
    local.ext = rep(-0.5, n_bins)
    
    bins <- data.frame(id=1:n_bins, 
                       chr=matdf[, "chr"], 
                       from.coord=matdf[, "from.coord"], 
                       to.coord=matdf[, "to.coord"] )
    matrix.data <- as.matrix( matdf[, (ncol(matdf) - nrow(matdf)+1 ):ncol(matdf)] )
    
    print("-- Done!")
    print("Step 0 : Done !!")
  }
  
    
  print("#########################################################################")
  print("Step 1 : Generating binSignals by computing bin-level contact frequencies")
  print("#########################################################################")
  ptm <- proc.time()
  for(i in 1:n_bins)
  {
    diamond = Get.Diamond.Matrix(mat.data=matrix.data, i=i, size=window.size)
    mean.cf[i] = mean(diamond)
  }
  
  eltm = proc.time() - ptm
  print(paste("Step 1 Running Time : ", eltm[3]))
  print("Step 1 : Done !!")
  
  print("#########################################################################")
  print("Step 2 : Detect TD boundaries based on binSignals")
  print("#########################################################################")
  
  ptm = proc.time()
  #gap.idx = Which.Gap.Region(matrix.data=matrix.data)
  #gap.idx = Which.Gap.Region2(mean.cf)
  gap.idx = Which.Gap.Region2(matrix.data=matrix.data, window.size)
  
  proc.regions = Which.process.region(rmv.idx=gap.idx, n_bins=n_bins, min.size=3)
  
  #print(proc.regions)
  
  for( i in 1:nrow(proc.regions))
  {
    start = proc.regions[i, "start"]
    end = proc.regions[i, "end"]
    
    print(paste("Process Regions from ", start, "to", end))
      
    local.ext[start:end] = Detect.Local.Extreme(x=mean.cf[start:end])
  }
  
  eltm = proc.time() - ptm
  print(paste("Step 2 Running Time : ", eltm[3]))
  print("Step 2 : Done !!")
  
  if(statFilter)
  {
    print("#########################################################################")
    print("Step 3 : Statistical Filtering of false positive TD boundaries")
    print("#########################################################################")
  
    ptm = proc.time()  
    print("-- Matrix Scaling....")
    n_bins <- as.double(n_bins)
    print(as.double(n_bins*n_bins))
    scale.matrix.data = matrix.data
    for( i in 1:(2*window.size) )
    {
      #diag(scale.matrix.data[, i:n_bins] ) = scale( diag( matrix.data[, i:n_bins] ) )
      scale.matrix.data[ seq(1+(n_bins*i), as.double(n_bins*n_bins), 1+n_bins) ] = scale( matrix.data[ seq(1+(n_bins*i), as.double(n_bins*n_bins), 1+n_bins) ] )
    }
    
    print("-- Compute p-values by Wilcox Ranksum Test")
    for( i in 1:nrow(proc.regions))
    {
      start = proc.regions[i, "start"]
      end = proc.regions[i, "end"]
      
      print(paste("Process Regions from ", start, "to", end))
      
      pvalue[start:end] <- Get.Pvalue(matrix.data=scale.matrix.data[start:end, start:end], size=window.size, scale=1)
    }
    print("-- Done!")
    
    print("-- Filtering False Positives")
    local.ext[intersect( union(which( local.ext==-1), which(local.ext==-1)), which(pvalue<0.05))] = -2
    local.ext[which(local.ext==-1)] = 0
    local.ext[which(local.ext==-2)] = -1
    print("-- Done!")
    
    eltm = proc.time() - ptm
    print(paste("Step 3 Running Time : ", eltm[3]))
    print("Step 3 : Done!")
  } else pvalue = 0

  domains = Convert.Bin.To.Domain.TMP(bins=bins, 
                                  signal.idx=which(local.ext==-1), 
                                  gap.idx=which(local.ext==-0.5), 
                                  pvalues=pvalue, 
                                  pvalue.cut=0.05)
  
  bins = cbind(bins, 
               local.ext = local.ext,
               mean.cf = mean.cf, 
               pvalue = pvalue)
  
  bedform = domains[, c("chr", "from.coord", "to.coord", "tag")]
  colnames(bedform) = c("chrom", "chromStart", "chromEnd", "name")
  
  if( !is.null(outFile) ) {
    print("#########################################################################")
    print("Writing Files")
    print("#########################################################################")
    
    outBinSignal =  paste(outFile, ".binSignal", sep="")
    print(paste("binSignal File :", outBinSignal) )
    write.table(bins, file=outBinSignal, quote=F, row.names=F, col.names=T, sep="\t")  
  
  
    outDomain = paste(outFile, ".domain", sep="")
    print(paste("Domain File :", outDomain) )
    write.table( domains, file=outDomain, quote=F, row.names=F, col.names=T, sep="\t")
    
    outBed = paste(outFile, ".bed", sep="")
    print(paste("Bed File : ", outBed))
    write.table( bedform, file=outBed, quote=F, row.names=F, col.names=F, sep="\t")
  }
  
  print("Done!!")
  
  print("Job Complete !")
  return(list(binSignal=bins, domain=domains, bed=bedform))
}

# @fn Get.Diamond.Matrix
# @param mat.data : N by N matrix, where each element indicate contact frequency
# @param i :integer, bin index
# @param size : integer, window size to expand from bin
# @retrun : matrix.
Get.Diamond.Matrix <- function(mat.data, i, size)
{
  n_bins = nrow( mat.data )
  if(i==n_bins) return(NA)
  
  lowerbound = max( 1, i-size+1 )
  upperbound = min( i+size, n_bins)
  
  return( mat.data[lowerbound:i, (i+1):upperbound] )
}

# @fn Which.process.region
# @param rmv.idx : vector of idx, remove index vector
# @param n_bins : total number of bins
# @param min.size : minimum size of bins
# @retrun : data.frame of proc.regions
Which.process.region <- function(rmv.idx, n_bins, min.size=3)
{
  gap.idx = rmv.idx
  
  proc.regions = data.frame(start=numeric(0), end=numeric(0))
  proc.set = setdiff(1:n_bins, gap.idx)
  n_proc.set = length(proc.set)
  
  i=1
  while(i < n_proc.set )
  {
    start = proc.set[i]
    j = i+1
    
    while(j <= n_proc.set)
    {
      if( proc.set[j] - proc.set[j-1] <= 1) j = j + 1
      else {
        proc.regions = rbind(proc.regions, c(start=start, end=proc.set[j-1]) )
        i = j
        break
      }
    }
    
    if(j >= n_proc.set ) {
      proc.regions = rbind(proc.regions, c(start=start, end=proc.set[j-1]) )
      break
    }
  }
  
  colnames(proc.regions) = c("start", "end")
  proc.regions <- proc.regions[ which( abs(proc.regions[,"end"] - proc.regions[, "start"]) >= min.size ), ]
  
  return(proc.regions)
}

# @fn Which.Gap.Region
# @breif version 0.0.1 used
# @param matrix.data : n by n matrix
# @return gap index
Which.Gap.Region <- function(matrix.data)
{
  n_bins = nrow(matrix.data)
  gap = rep(0, n_bins)
  
  i=1
  while(i < n_bins)
  {
    j = i + 1
    while( j <= n_bins)
    {
      if( sum( matrix.data[i:j, i:j]) == 0 ) {
        gap[i:j] = -0.5
        j = j+1  
        #if(j-i > 1) gap[i:j]=-0.5
        #j=j+1
      } else break
    }
    
    i = j  
  }
  
  idx = which(gap == -0.5)
  return(idx)
}

# @fn Which.Gap.Region3
# @param matrix.data : n by n matrix
# @return gap index
Which.Gap.Region3 <- function(mean.cf)
{
  n_bins = length(mean.cf)
  gapidx = which(mean.cf==0)
  
  return(gapidx)
}

# @fn Which.Gap.Region2
# @breif version 0.0.2 used
# @param matrix.data : n by n matrix
# @return gap index
Which.Gap.Region2 <- function(matrix.data, w)
{
  n_bins = nrow(matrix.data)
  gap = rep(0, n_bins)
  
  for(i in 1:n_bins)
  {
    if( sum( matrix.data[i, max(1, i-w):min(i+w, n_bins)] ) == 0 ) gap[i]=-0.5
  }
  
  idx = which(gap == -0.5)
  return(idx)
}

# @fn Detect.Local.Extreme
# @param x : original signal to find local minima
# @return vector of local extrme, -1 if the index is local minimum, 1 if the index is local maxima, 0 otherwise.
Detect.Local.Extreme <- function(x)
{
  n_bins = length(x)
  ret = rep(0, n_bins)
  x[is.na(x)]=0
  
  if(n_bins <= 3)
  {
    ret[which.min(x)]=-1
    ret[which.max(x)]=1
    
    return(ret)
  }
  # Norm##################################################3
  new.point = Data.Norm(x=1:n_bins, y=x)
  x=new.point$y
  ##################################################
  cp = Change.Point(x=1:n_bins, y=x)
  
  if( length(cp$cp) <= 2 ) return(ret)
  if( length(cp$cp) == n_bins) return(ret)
  for(i in 2:(length(cp$cp)-1))
  {
    if( x[cp$cp[i]] >= x[cp$cp[i]-1] && x[cp$cp[i]] >= x[cp$cp[i]+1] ) ret[cp$cp[i]] = 1
    else if(x[cp$cp[i]] < x[cp$cp[i]-1] && x[cp$cp[i]] < x[cp$cp[i]+1]) ret[cp$cp[i]] = -1
    
    min.val = min( x[ cp$cp[i-1] ], x[ cp$cp[i] ] )
    max.val = max( x[ cp$cp[i-1] ], x[ cp$cp[i] ] )
    
    if( min( x[cp$cp[i-1]:cp$cp[i]] ) < min.val ) ret[ cp$cp[i-1] - 1 + which.min( x[cp$cp[i-1]:cp$cp[i]] ) ] = -1
    if( max( x[cp$cp[i-1]:cp$cp[i]] ) > max.val ) ret[ cp$cp[i-1] - 1 + which.max( x[cp$cp[i-1]:cp$cp[i]] ) ] = 1
  }

  return(ret)
}

# @fn Data.Norm
# @param x : x axis vector
# @param x : y axis vector
# @return list of normalized x and y
Data.Norm <- function(x, y)
{
  ret.x = rep(0, length(x))
  ret.y = rep(0, length(y))
  
  ret.x[1] = x[1]
  ret.y[1] = y[1]
  
  diff.x = diff(x)
  diff.y = diff(y)
  
  scale.x = 1 / mean( abs(diff(x) ) )
  scale.y = 1 / mean( abs( diff(y) ) )
  
  #print(scale.x)
  #print(scale.y)
  
  for(i in 2:length(x))
  {
    ret.x[i] = ret.x[i-1] + (diff.x[i-1]*scale.x)
    ret.y[i] = ret.y[i-1] + (diff.y[i-1]*scale.y)
  }
  
  return(list(x=ret.x, y=ret.y))
}

# @fn Change.Point
# @param x : x axis vector
# @param x : y axis vector
# @return change point index in x vector,
# Note that the first and the last point will be always change point
Change.Point <- function( x, y )
{
  if( length(x) != length(y)) 
  {
    print("ERROR : The length of x and y should be the same")
    return(0)
  }
  
  n_bins <- length(x)
  Fv <- rep(NA, n_bins)
  Ev <- rep(NA, n_bins)
  cp <- 1
  
  i=1
  Fv[1]=0
  while( i < n_bins )
  {
    j=i+1
    Fv[j] = sqrt( (x[j]-x[i])^2 + (y[j] - y[i] )^2 )
    
    while(j<n_bins)
    {
      j=j+1
      k=(i+1):(j-1)
      Ev[j] = ( sum( abs( (y[j]-y[i] )*x[k] - (x[j] -x[i])*y[k] - (x[i]*y[j]) + (x[j]*y[i]) ) ) / sqrt( (x[j]-x[i])^2 + (y[j] - y[i] )^2 ) )
      Fv[j] = sqrt( (x[j]-x[i])^2 + (y[j] - y[i])^2 ) - ( sum( abs( (y[j]-y[i] )*x[k] - (x[j] -x[i])*y[k] - (x[i]*y[j]) + (x[j]*y[i]) ) ) / sqrt( (x[j]-x[i])^2 + (y[j] - y[i] )^2 ) )
      
      #################################################
      #Not Original Code
      if( is.na(Fv[j]) || is.na(Fv[j-1]) ) {
        j = j-1
        cp <- c(cp, j)
        break
      }
      ####################################################3
      if(Fv[j] < Fv[j-1] ) {
        j = j - 1
        cp <- c(cp, j )
        break
      }
    }
    i=j
  }
  
  cp <- c(cp, n_bins)
  
  return(list(cp=cp, objF=Fv, errF=Ev))
}

# @fn Get.Pvalue
# @param matrix.data : matrix
# @param size : size to extend
# @param scale : scale parameter if necessary. deprecated parameter
# @return computed p-value vector
Get.Pvalue <- function( matrix.data, size, scale=1 )
{
  n_bins = nrow(matrix.data)
  pvalue <- rep(1, n_bins)
  
  for( i in 1:(n_bins-1) )
  {
    dia = as.vector( Get.Diamond.Matrix2(matrix.data, i, size=size) )
    ups = as.vector( Get.Upstream.Triangle(matrix.data, i, size=size) )
    downs = as.vector( Get.Downstream.Triangle(matrix.data, i, size=size) )
    
    wil.test =  wilcox.test(x=dia*scale, y=c(ups, downs), alternative="less", exact=F)
    pvalue[i] = wil.test$p.value  
    
    #print(paste(i, "=", wil.test$p.value) )  
  }
  
  pvalue[ is.na(pvalue) ] = 1
  return(pvalue)
}

# @fn Get.Upstream.Triangle
# @param mat.data : matrix data
# @param i : bin index
# @param size : size of window to extend
# @return upstream triangle matrix
Get.Upstream.Triangle <- function(mat.data, i, size)
{
  n_bins = nrow(mat.data)
  
  lower = max(1, i-size)
  tmp.mat = mat.data[lower:i, lower:i]
  return( tmp.mat[ upper.tri( tmp.mat, diag=F ) ] )
}

# @fn Get.Downstream.Triangle
# @param mat.data : matrix data
# @param i : bin index
# @param size : size of window to extend
# @return downstream triangle matrix
Get.Downstream.Triangle <- function(mat.data, i, size)
{
  n_bins = nrow(mat.data)
  if(i==n_bins) return(NA)
  
  upperbound = min(i+size, n_bins)
  tmp.mat = mat.data[(i+1):upperbound, (i+1):upperbound]
  return( tmp.mat[ upper.tri( tmp.mat, diag=F ) ] )
}

# @fn Get.Diamond.Matrix2
# @param mat.data : matrix data
# @param i : bin index
# @param size : size of window to extend
# @return diamond matrix
Get.Diamond.Matrix2 <- function(mat.data, i, size)
{
  n_bins = nrow(mat.data)
  new.mat = matrix(rep(NA, size*size), nrow=size, ncol=size)
  
  for(k in 1:size)
  {
    if(i-(k-1) >= 1 && i < n_bins)
    {
      lower = min(i+1, n_bins)
      upper = min(i+size, n_bins)
      
      new.mat[size-(k-1), 1:(upper-lower+1)] = mat.data[i-(k-1), lower:upper]
    }
  }
  
  return(new.mat)
}
  
# @fn Convert.Bin.To.Domain
# @param bins : bin information
# @param signal.idx : signal index
# @param signal.idx : gap index
# @param pvalues : pvalue vector
# @param pvalue.cut : pvalue threshold
# @return dataframe storing domain information
Convert.Bin.To.Domain <- function(bins, signal.idx, gap.idx, pvalues=NULL, pvalue.cut=NULL)
{
  n_bins = nrow(bins)
  ret = data.frame(chr=character(0), from.id=numeric(0), from.coord=numeric(0), to.id=numeric(0), to.coord=numeric(0), tag=character(0), size=numeric(0))
  levels( x=ret[, "tag"] ) = c("domain", "gap", "boundary")
  
  rmv.idx = setdiff(1:n_bins, gap.idx)
  proc.region = Which.process.region(rmv.idx, n_bins, min.size=0)
  from.coord = bins[proc.region[, "start"], "from.coord"]
  n_procs = nrow(proc.region)
  gap = data.frame(chr=rep( bins[1, "chr"], n_procs), from.id=rep(0, n_procs), from.coord=from.coord, to.id=rep(0, n_procs), to.coord=rep(0, n_procs), tag=rep("gap", n_procs), size=rep(0, n_procs), stringsAsFactors=F)
  
  rmv.idx = union(signal.idx, gap.idx)
  proc.region = Which.process.region(rmv.idx, n_bins, min.size=0)
  n_procs = nrow(proc.region)
  from.coord = bins[proc.region[, "start"], "from.coord"]
  domain = data.frame(chr=rep( bins[1, "chr"], n_procs), from.id=rep(0, n_procs), from.coord=from.coord, to.id=rep(0, n_procs), to.coord=rep(0, n_procs), tag=rep("domain", n_procs), size=rep(0, n_procs), stringsAsFactors=F)
  
  rmv.idx = setdiff(1:n_bins, signal.idx)
  proc.region = as.data.frame( Which.process.region(rmv.idx, n_bins, min.size=1) )
  n_procs = nrow(proc.region)
  if(n_procs>0)
  {
    from.coord = bins[proc.region[, "start"]+1, "from.coord"]  
    boundary = data.frame(chr=rep( bins[1, "chr"], n_procs), from.id=rep(0, n_procs), from.coord=from.coord, to.id=rep(0, n_procs), to.coord=rep(0, n_procs), tag=rep("boundary", n_procs), size=rep(0, n_procs), stringsAsFactors=F)
    ret = rbind(ret, boundary)
  }
  
  ret = rbind(gap, domain)
  ret = ret[order(ret[,3]), ]
  
  ret[, "to.coord"] = c(ret[2:nrow(ret), "from.coord"], bins[n_bins, "to.coord"])
  ret[, "from.id"] = match( ret[, "from.coord"], bins[, "from.coord"] )
  ret[, "to.id"] = match(ret[, "to.coord"], bins[, "to.coord"])
  ret[, "size"] = ret[,"to.coord"]-ret[,"from.coord"]
  
  if(!is.null(pvalues) && !is.null(pvalue.cut))
  {
    for(i in 1:nrow(ret))
    {
      if(ret[i, "tag"]=="domain")
      {
        domain.bins.idx = ret[i, "from.id"]:ret[i, "to.id"]
        p.value.constr = which( pvalues[ domain.bins.idx ] < pvalue.cut )
        
        if( length(domain.bins.idx) == length(p.value.constr)) ret[i, "tag"] = "boundary"
      }
    }
  }
  
  new.bdr.set = data.frame(chr=character(0), from.id=numeric(0), from.coord=numeric(0), to.id=numeric(0), to.coord=numeric(0), tag=character(0), size=numeric(0))
  stack.bdr = data.frame(chr=character(0), from.id=numeric(0), from.coord=numeric(0), to.id=numeric(0), to.coord=numeric(0), tag=character(0), size=numeric(0))
  
  i=1
  while(i <= nrow(ret))
  {
    if( ret[i, "tag"] == "boundary" )
    {
      stack.bdr = rbind(stack.bdr, ret[i, ])
    } else if(nrow(stack.bdr)>0) {
      new.bdr = data.frame(chr=bins[1, "chr"], 
                           from.id = min( stack.bdr[, "from.id"]), 
                           from.coord=min(stack.bdr[, "from.coord"]),
                           to.id = max( stack.bdr[, "to.id"]),
                           to.coord=max(stack.bdr[, "to.coord"]),
                           tag="boundary",
                           size=max(stack.bdr[, "to.coord"]) - min(stack.bdr[, "from.coord"]))
      new.bdr.set = rbind(new.bdr.set, new.bdr)
      stack.bdr = data.frame(chr=character(0), from.id=numeric(0), from.coord=numeric(0), to.id=numeric(0), to.coord=numeric(0), tag=character(0), size=numeric(0))
    }
    
    i = i + 1
  }
  
  ret = rbind( ret[ ret[, "tag"]!="boundary", ], new.bdr.set )
  ret = ret[order(ret[, "to.coord"]), ]
  
  return(ret)
}


# @fn Convert.Bin.To.Domain
# @param bins : bin information
# @param signal.idx : signal index
# @param signal.idx : gap index
# @param pvalues : pvalue vector
# @param pvalue.cut : pvalue threshold
# @return dataframe storing domain information
Convert.Bin.To.Domain.TMP <- function(bins, signal.idx, gap.idx, pvalues=NULL, pvalue.cut=NULL)
{
  n_bins = nrow(bins)
  ret = data.frame(chr=character(0), from.id=numeric(0), from.coord=numeric(0), to.id=numeric(0), to.coord=numeric(0), tag=character(0), size=numeric(0))
  levels( x=ret[, "tag"] ) = c("domain", "gap", "boundary")
  
  rmv.idx = setdiff(1:n_bins, gap.idx)
  proc.region = Which.process.region(rmv.idx, n_bins, min.size=0)
  from.coord = bins[proc.region[, "start"], "from.coord"]
  n_procs = nrow(proc.region)
  gap = data.frame(chr=rep( bins[1, "chr"], n_procs), from.id=rep(0, n_procs), from.coord=from.coord, to.id=rep(0, n_procs), to.coord=rep(0, n_procs), tag=rep("gap", n_procs), size=rep(0, n_procs), stringsAsFactors=F)
  
  rmv.idx = union(signal.idx, gap.idx)
  proc.region = Which.process.region(rmv.idx, n_bins, min.size=0)
  n_procs = nrow(proc.region)
  from.coord = bins[proc.region[, "start"], "from.coord"]
  domain = data.frame(chr=rep( bins[1, "chr"], n_procs), from.id=rep(0, n_procs), from.coord=from.coord, to.id=rep(0, n_procs), to.coord=rep(0, n_procs), tag=rep("domain", n_procs), size=rep(0, n_procs), stringsAsFactors=F)
  
  rmv.idx = setdiff(1:n_bins, signal.idx)
  proc.region = as.data.frame( Which.process.region(rmv.idx, n_bins, min.size=1) )
  n_procs = nrow(proc.region)
  if(n_procs>0)
  {
    from.coord = bins[proc.region[, "start"]+1, "from.coord"]  
    boundary = data.frame(chr=rep( bins[1, "chr"], n_procs), from.id=rep(0, n_procs), from.coord=from.coord, to.id=rep(0, n_procs), to.coord=rep(0, n_procs), tag=rep("boundary", n_procs), size=rep(0, n_procs), stringsAsFactors=F)
    ret = rbind(ret, boundary)
  }
  
  ret = rbind(gap, domain)
  ret = ret[order(ret[,3]), ]
  
  ret[, "to.coord"] = c(ret[2:nrow(ret), "from.coord"], bins[n_bins, "to.coord"])
  ret[, "from.id"] = match( ret[, "from.coord"], bins[, "from.coord"] )
  ret[, "to.id"] = match(ret[, "to.coord"], bins[, "to.coord"])
  ret[, "size"] = ret[,"to.coord"]-ret[,"from.coord"]
  
  if(!is.null(pvalues) && !is.null(pvalue.cut))
  {
    for(i in 1:nrow(ret))
    {
      if(ret[i, "tag"]=="domain")
      {
        domain.bins.idx = ret[i, "from.id"]:ret[i, "to.id"]
        p.value.constr = which( pvalues[ domain.bins.idx ] < pvalue.cut )
        
        if( length(domain.bins.idx) == length(p.value.constr)) ret[i, "tag"] = "boundary"
      }
    }
  }
  
  return(ret)
}

TopDom
}) ## local()
