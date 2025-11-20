library(GenomicRanges)
library(rtracklayer)

check_coverage_argument <- function(cvg, regions=NULL) {
  if(class(cvg) == "character") {
    stopifnot(file.exists(cvg))
    if(is.null(regions)) {
      cvg <- import(cvg, as="RleList")
    } else {
      cvg <- import(cvg, which=regions, as="RleList")
    }
  }
  stopifnot(is(cvg, "RleList")) 
  cvg
}

grangesApply <- function(regions, cvg, func, ...) {
  stopifnot(class(regions) == "GRanges")
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(viewApply(
           X        = Views(cvg, as(regions, "IntegerRangesList")),
           FUN      = func, ...,
           simplify = FALSE
         ))
  ans[oo] <- ans  # restore original order
  ans
}

grangesBuiltinViewHelper <- function(regions, cvg, func) {
  stopifnot(class(regions) == "GRanges")
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(func(Views(cvg, as(regions, "IntegerRangesList"))))
  ans[oo] <- ans  # restore original order
  ans
}

grangesSums <- function(regions, cvg) {
  grangesBuiltinViewHelper(regions, cvg, viewSums)
}

grangesMeans <- function(regions, cvg) {
  grangesBuiltinViewHelper(regions, cvg, viewMeans)
}

grangesWhichMaxs <- function(regions, cvg) {
  grangesBuiltinViewHelper(regions, cvg, viewWhichMaxs)
}

grangesWhichMins <- function(regions, cvg) {
  grangesBuiltinViewHelper(regions, cvg, viewWhichMins)
}

grangesMaxs <- function(regions, cvg) {
  grangesBuiltinViewHelper(regions, cvg, viewMaxs)
}

grangesMins <- function(regions, cvg) {
  grangesBuiltinViewHelper(regions, cvg, viewMins)
}

grangesMatrix <- function(regions, cvg) {
  stopifnot(class(regions) == "GRanges")
  stopifnot(length(unique(width(regions))) == 1)
  
  reads.list <- grangesApply(regions, cvg, as.numeric)
  reads.m <- matrix(unlist(reads.list), nrow=length(regions), byrow=TRUE)
  
  neg_i <- which(strand(regions) == "-")
  if(length(neg_i) > 0) {
    reads.m[neg_i, ] <- reads.m[neg_i, ncol(reads.m):1]
  }
  reads.m
}

total_signal_rlelist <- function(cvg) {
  stopifnot(is(cvg, "RleList"))
  sum(as.numeric(sum(cvg)))
}

delete_if_older_than <- function(file, timestamp) {
  stopifnot(length(file) == 1)
  if(file.info(file)$mtime < timestamp) {
    if(!file.remove(file)) {
      stop("Could not delete out-of-date cache file: ", file)
    }
  }
}

totalSignal <- function(cvg) {
  if(class(cvg) == "character") {
    stopifnot(length(cvg) == 1 & file.exists(cvg))
    cache.file <- paste0(cvg, ".ts.rds")
    delete_if_older_than(cache.file, file.info(cvg)$mtime)
    if(file.exists(cache.file)) {
      return(readRDS(cache.file))
    } else {
      cvg <- check_coverage_argument(cvg)
      ts <- total_signal_rlelist(cvg)
      saveRDS(ts, file=cache.file)
      return(ts)
    }
  } else {
    return(total_signal_rlelist(cvg))
  }
}

load_object <- function(filename) {
  if(length(grep("rdata$", tolower(filename))) > 0) {
    o <- get(load(filename))    
  } else {
    o <- readRDS(filename)    
  }  
  updateObject(o)
}

filter_chrs <- function(gr, exclude="H|M|U|_") {
  exclude.chrs <- grep(exclude, seqlevels(gr))
  if(length(exclude.chrs) > 0) seqlevels(gr, force=TRUE) <- seqlevels(gr)[-exclude.chrs]
  gr
}


