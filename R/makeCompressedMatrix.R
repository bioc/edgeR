makeCompressedMatrix <- function(x, dims=NULL, byrow=TRUE) 
# Construct a CompressedMatrix from a scalar, vector or matrix.
#
# Created by Aaron Lun, 24 Sep 2016.
# Revised by Lizhong Chen, 9 May 2024.
{
    repeat.row <- repeat.col <- FALSE
	if (is.matrix(x)) {
		if (inherits(x, "CompressedMatrix")) {
			return(x)
		}
        if (!is.null(dims)){
            xdim <- dim(x)
            if(xdim[1]==1L && xdim[2]==1L){
                repeat.row <- repeat.col <- TRUE
		        x <- matrix(x)               
            }
            else if(xdim[1]==1L && xdim[2]>=2L){
                if(xdim[2]!=dims[2]){
                    stop("'dims[2]' should equal length of row vector 'x'")
                }
                if(byrow){
                    repeat.row <- TRUE
                }
                else{
                    dims <- xdim    
                }
            }
            else if(xdim[1]>=2L && xdim[2]==1L){
                if(xdim[1]!=dims[1]){
                    stop("'dims[1]' should equal length of column vector 'x'")
                }
                if(!byrow){
                    repeat.col <- TRUE
                }
                else{
                    dims <- xdim
                }
            }
            else{
                dims <- xdim
            }
        }
        else{
            dims <- dim(x)
        }
	} else if (length(x)==1L) {
        repeat.row <- repeat.col <- TRUE
		x <- matrix(x)
	} else {
		if (!byrow) {
            if (dims[1]!=length(x)) { 
                stop("'dims[1]' should equal length of 'x'")
            }
			x <- cbind(x)
            repeat.col <- TRUE
		} else {
            if (dims[2]!=length(x)) { 
                stop("'dims[2]' should equal length of 'x'")
            }
			x <- rbind(x)
            repeat.row <- TRUE
		}
	}

    dimnames(x) <- NULL
    class(x) <- "CompressedMatrix"
	attr(x, "Dims") <- as.integer(dims)
    attr(x, "repeat.row") <- repeat.row
    attr(x, "repeat.col") <- repeat.col
	return(x)
}

.strip_to_matrix <- function(x) {
	dims <- attr(x, "dim")
	attributes(x) <- NULL
	attr(x, "dim") <- dims
	return(x)
}

dim.CompressedMatrix <- function(x) 
# Getting dimensions.
#
# written by Aaron Lun
# created 21 June 2017
{
	attr(x, "Dims")
}

length.CompressedMatrix <- function(x)
# Getting length. 
#
# written by Aaron Lun
# created 25 January 2018
# last modified 1 March 2018
{
	prod(attr(x,"Dims"))
}

`[.CompressedMatrix` <- function(x, i, j, drop=TRUE)
# Subsetting for CompressedMatrix objects.
#
# written by Aaron Lun
# created 24 September 2016
# last modified 21 June 2017
{
	Nargs <- nargs() - !(missing(drop))
	if (Nargs<3L) {
		return(as.matrix(x)[i])
	}

	raw.mat <- .strip_to_matrix(x)
	row.status <- attr(x, "repeat.row") 
	col.status <- attr(x, "repeat.col")

	if (!row.status && !missing(i)) {
		raw.mat <- raw.mat[i,,drop=FALSE]
	}
	if (!col.status && !missing(j)) {
		raw.mat <- raw.mat[,j,drop=FALSE]
	}

	nr <- nrow(x)
	if (!missing(i)) {
		ref <- seq_len(nr)
		nr <- length(ref[i])
	} 
	nc <- ncol(x)
	if (!missing(j)) {
		ref <- seq_len(nc)
		nc <- length(ref[j])
	}

	class(raw.mat) <- class(x)
	attr(raw.mat, "Dims") <- c(nr, nc)
	attr(raw.mat, "repeat.row") <- row.status
	attr(raw.mat, "repeat.col") <- col.status

	if (drop && 
		((!missing(i) && length(i)==1L) ||
		 (!missing(j) && length(j)==1L))) {
		raw.mat <- as.vector(as.matrix(raw.mat))
	} 
	return(raw.mat)
}

`[<-.CompressedMatrix` <- function(x, i, j, value) 
# Subset assignment for CompressedMatrix objects.
#
# written by Aaron Lun
# created 25 January 2018
{
	ref <- as.matrix(x)
	if (is(value, "CompressedMatrix")) { 
		value <- as.matrix(value)
	}

	if (nargs() < 4L) {
		ref[i] <- value
	} else {
		ref[i,j] <- value
	}
	makeCompressedMatrix(ref, attr(x, "Dims"), TRUE)
}

as.matrix.CompressedMatrix <- function(x, ...) 
# Expanding it to a full matrix.
#
# written by Aaron Lun
# created 26 September 2016
# last modified 21 June 2017
{
	raw.mat <- .strip_to_matrix(x)
	row.status <- attr(x, "repeat.row") 
	col.status <- attr(x, "repeat.col")

	if (row.status) {
		raw.mat <- matrix(raw.mat, nrow(x), ncol(x), byrow=TRUE)				
	} else if (col.status) {
		raw.mat <- matrix(raw.mat, nrow(x), ncol(x))
	} else {
		raw.mat <- as.matrix(raw.mat)
	}
	return(raw.mat)
}

rbind.CompressedMatrix <- function(...) 
# Rbinding things together.
# 
# written by Aaron Lun
# created 21 June 2017	
{
	everything <- list(...)
	nobjects <- length(everything)
	if (nobjects==1) {
		return(everything[[1]])
	}
	all.nr <- sum(unlist(lapply(everything, nrow)))
	
	col.rep <- logical(nobjects)
	row.rep <- logical(nobjects)
	for (i in seq_along(everything)) { 
		x <- everything[[i]]
		col.rep[i] <- attr(x, "repeat.col") 
		row.rep[i] <- attr(x, "repeat.row")
	}

	# If everything is column repeats, we can do a naive concatenation. 
	if (all(col.rep)) { 
		collected.vals <- vector("list", nobjects)
		all.nc <- ncol(everything[[1]])
		for (i in seq_along(everything)) {
			current <- everything[[i]]
			if (!identical(all.nc, ncol(current))) {
				stop("cannot combine CompressedMatrix objects with different number of columns")
			}
			collected.vals[[i]] <- rep(.strip_to_matrix(current), length.out=nrow(current))
		}
		return(makeCompressedMatrix(unlist(collected.vals), dims=c(all.nr, all.nc), byrow=FALSE))
	}

	# If everything is row repeats AND values are all equal, we can just modify the nr.
	if (all(row.rep)) {
		okay <- TRUE
		ref <- .strip_to_matrix(everything[[1]])
		for (i in 2:length(everything)) {
			current <- .strip_to_matrix(everything[[i]])
			if (!isTRUE(all.equal(everything[[i]], ref))) {
				okay <- FALSE
				break
			}
		}
		if (okay) {
			current <- everything[[1]]
			attr(current, "Dims")[1] <- all.nr
			return(current)
		}
	}

	# Otherwise, expanding each element and rbinding them.
	for (i in seq_along(everything)) { 
		everything[[i]] <- as.matrix(everything[[i]])
	} 
	return(makeCompressedMatrix(do.call(rbind, everything)))
}

cbind.CompressedMatrix <- function(...) 
# Cbinding things together.
# 
# written by Aaron Lun
# created 21 June 2017	
{
	everything <- list(...)
	nobjects <- length(everything)
	if (nobjects==1) {
		return(everything[[1]])
	}
	all.nc <- sum(unlist(lapply(everything, ncol)))
	
	col.rep <- logical(nobjects)
	row.rep <- logical(nobjects)
	for (i in seq_along(everything)) { 
		x <- everything[[i]]
		col.rep[i] <- attr(x, "repeat.col") 
		row.rep[i] <- attr(x, "repeat.row")
	}

	# If everything is row repeats, we can do a naive concatenation. 
	if (all(row.rep)) { 
		collected.vals <- vector("list", nobjects)
		all.nr <- nrow(everything[[1]])
		for (i in seq_along(everything)) {
			current <- everything[[i]]
			if (!identical(all.nr, nrow(current))) {
				stop("cannot combine CompressedMatrix objects with different number of rows")
			}
			collected.vals[[i]] <- rep(.strip_to_matrix(current), length.out=ncol(current))
		}
		return(makeCompressedMatrix(unlist(collected.vals), dims=c(all.nr, all.nc), byrow=TRUE))
	}

	# If everything is column repeats AND values are all equal, we can just modify the nc.
	if (all(col.rep)) {
		okay <- TRUE
		ref <- .strip_to_matrix(everything[[1]])
		for (i in 2:length(everything)) {
			current <- .strip_to_matrix(everything[[i]])
			if (!isTRUE(all.equal(everything[[i]], ref))) {
				okay <- FALSE
				break
			}
		}
		if (okay) {
			current <- everything[[1]]
			attr(current, "Dims")[2] <- all.nc
			return(current)
		}
	}

	# Otherwise, expanding each element and rbinding them.
	for (i in seq_along(everything)) { 
		everything[[i]] <- as.matrix(everything[[i]])
	} 
	return(makeCompressedMatrix(do.call(cbind, everything)))
}

Ops.CompressedMatrix <- function(e1, e2)
# A function that performs some binary operation on two CompressedMatrix objects,
# in a manner that best preserves memory usage.
#
# written by Aaron Lun
# created 26 September 2016
# last modified 30 June 2017
{
	if (!inherits(e1, "CompressedMatrix")) {
		e1 <- makeCompressedMatrix(e1, dim(e2), byrow=FALSE) # Promoted to column-major CompressedMatrix 
	}
	if (!inherits(e2, "CompressedMatrix")) {
		e2 <- makeCompressedMatrix(e2, dim(e1), byrow=FALSE)	   
	}
	if (!identical(dim(e1), dim(e2))) {
		stop("CompressedMatrix dimensions should be equal for binary operations")
	}
	
	row.rep <- attr(e1, "repeat.row") && attr(e2, "repeat.row")
	col.rep <- attr(e1, "repeat.col") && attr(e2, "repeat.col")

	if (row.rep || col.rep) { 
		new.dim <- dim(e1)
		e1 <- as.vector(.strip_to_matrix(e1))
		e2 <- as.vector(.strip_to_matrix(e2))
		outcome <- NextMethod(.Generic)
		outcome <- makeCompressedMatrix(outcome, new.dim, byrow=row.rep)
	} else {
		e1 <- as.matrix(e1)
		e2 <- as.matrix(e2)
		outcome <- NextMethod(.Generic)
		outcome <- makeCompressedMatrix(outcome)
	}

	return(outcome)
}

.compressOffsets <- function(y, offset, lib.size=NULL) 
# A convenience function to avoid repeatedly having to write the code below.
# If provided, offsets take precedence over the library size.
# If neither are provided, library sizes are automatically computed
# as the sum of counts in the count matrix 'y'.
# A prefunctory check for finite values is performed in the C++ code.
# If 'offset' is already of the CompressedMatrix class, then 
# we assume it's already gone through this once so we don't do it again.
# Created 3 Oct 2016. Last modified 28 Aug 2024.
{
	if(inherits(offset, "CompressedMatrix")) return(offset)
	if (is.null(offset)) {
		if(is.null(lib.size)) lib.size <- colSums(y)
		offset <- log(lib.size)
	}
	if(!is.double(offset)) storage.mode(offset) <- "double"
	m <- min(offset)
	if(!is.finite(m)) stop("offsets must be finite values")
	makeCompressedMatrix(offset, dim(y), byrow=TRUE)
}

.compressWeights <- function(y, weights=NULL) 
# A convenience function to avoid repeatedly having to write the code below.
# All weights default to 1 if not specified.
# A prefunctory check for finite, positive values is performed in the C code.
# If 'weights' is already a CompressedMatrix, then we assume it's 
# already gone through this and don't do it again.
# Created 3 Oct 2016. Last modified 28 Aug 2024.
{
	if(inherits(weights, "CompressedMatrix")) return(weights)
	if(is.null(weights)) weights <- 1
	if(!is.double(weights)) storage.mode(weights) <- "double"
	m <- min(weights)
	if(is.na(m)) stop("NA weights not allowed")
	if(m <= 0) stop("Weights must be positive")
	makeCompressedMatrix(weights, dim(y), byrow=TRUE)
}

.compressPrior <- function(y, prior.count) 
# Again for the prior counts, checking for non-negative finite values.
# Skipping the check if it's already a CompressedMatrix object.
# Created 3 Oct 2016. Last modified 28 Aug 2024.
{
	if(inherits(prior.count, "CompressedMatrix")) return(prior.count)
	if(!is.double(prior.count)) storage.mode(prior.count) <- "double"
	m <- min(prior.count)
	if(is.na(m)) stop("NA prior counts not allowed")
	if(m < 0) stop("Negative prior counts not allowed")
	makeCompressedMatrix(prior.count, dim(y), byrow=FALSE)
}

.compressDispersions <- function(y, dispersion) 
# Again for the dispersions, checking for non-negative finite values.
# Skipping the check if it's already a CompressedMatrix object.
# Created 3 Oct 2016. Last modified 28 Aug 2024.
{
	if(inherits(dispersion, "CompressedMatrix")) return(dispersion)
	if(!is.double(dispersion)) storage.mode(dispersion) <- "double"
	m <- min(dispersion)
	if(is.na(m)) stop("NA dispersions not allowed")
	if(m < 0) stop("Negative dispersions not allowed")
	makeCompressedMatrix(dispersion, dim(y), byrow=FALSE)
}

