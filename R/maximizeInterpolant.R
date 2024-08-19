maximizeInterpolant <- function( x, y ) 
# Takes an ordered set of spline points and a likelihood matrix where each row 
# corresponds to a tag and each column corresponds to a spline point. It then calculates the 
# position at which the maximum interpolated likelihood occurs for each by solving the derivative
# of the spline function.
# C++ by Aaron Lun.
# C by Lizhong Chen.
# Last modified 13 May 2024.
{
    if (is.vector(y)) {
        y<-rbind(y)
        warning("coverting vector of likelihoods to matrix format for interpolation")
    }
    if (length(x)!=ncol(y)) { 
        stop("number of columns must equal number of spline points")
    } else if (is.unsorted(x) || anyDuplicated(x)) {
        stop("spline points must be unique and sorted")
    }

#	Performing some type checking.
    out<-.Call(.cxx_maximize_interpolant, x, y)
    return(out)
}
