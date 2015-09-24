#include "glm.h"
#include "matvec_check.h"

SEXP R_one_group (SEXP nlib, SEXP ntag, SEXP y, SEXP disp, SEXP offsets, SEXP weights, SEXP max_iterations, SEXP tolerance, SEXP beta) try {
    const int num_tags=asInteger(ntag);
    const int num_libs=asInteger(nlib);
    if (!isNumeric(disp)) { throw std::runtime_error("dispersion matrix must be double precision"); }
    if (num_tags*num_libs !=LENGTH(disp)) { throw std::runtime_error("dimensions of dispersion vector is not equal to number of tags"); }
	if (num_tags*num_libs != LENGTH(y) ) { throw std::runtime_error("dimensions of the count table are not as specified"); }  // Checking that it is an exact division.
  
	if (!isNumeric(beta)) { throw std::runtime_error("beta start vector must be double precision"); }
	const int blen=LENGTH(beta);
	const bool use_old_start=(blen!=0);
	if (use_old_start && blen!=num_tags) { 
		throw std::runtime_error("non-empty start vector must have length equal to the number of tags"); 		
	}
	const double* bsptr=REAL(beta);
   
	const int maxit=asInteger(max_iterations);
	const double tol=asReal(tolerance);
 
    // Setting up some iterators. We provide some flexibility to detecting numeric-ness.
	double *ydptr=NULL;
	int* yiptr=NULL;
	double* yptr=(double*)R_alloc(num_libs, sizeof(double));
	bool is_integer=isInteger(y);
	if (is_integer) { 
		yiptr=INTEGER(y); 
	} else { 
		if (!isNumeric(y)) { throw std::runtime_error("count matrix must be integer or double-precision"); }
		ydptr=REAL(y); 
	}
    matvec_check allo(num_libs, num_tags, offsets, false, "offset");
	const double* const* optr2=allo.access();
	matvec_check allw(num_libs, num_tags, weights, false, "weight", 1);
	const double* const* wptr2=allw.access();
	const double* dptr=REAL(disp);

    // Setting up beta for output. 
	SEXP output=PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(output, 0, allocVector(REALSXP, num_tags));
	SET_VECTOR_ELT(output, 1, allocVector(LGLSXP, num_tags));
    double* bptr=REAL(VECTOR_ELT(output, 0));
	int* cptr=LOGICAL(VECTOR_ELT(output, 1));
	try {
        
    	// Iterating through tags and fitting.
    	int lib=0;
    	for (int tag=0; tag<num_tags; ++tag) {
			if (is_integer) { 
				for (lib=0; lib<num_libs; ++lib) { yptr[lib]=double(yiptr[lib]); }	
				yiptr+=num_libs;
			} else {
				yptr=ydptr;
				ydptr+=num_libs;
			}
			std::pair<double, bool> out=glm_one_group(num_libs, maxit, tol, *optr2,
#ifdef WEIGHTED					
					*wptr2,
#endif					
					yptr, dptr, (use_old_start ? bsptr[tag] : R_NaReal));

			bptr[tag]=out.first;
			cptr[tag]=out.second;
			dptr+=num_libs;
			allo.advance();
			allw.advance();
    	}
	} catch (std::exception& e) { 
		UNPROTECT(1);
		throw; 
	}

	// Returning everything as a list.
    UNPROTECT(1); 
    return output;
} catch (std::exception& e) {
	return mkString(e.what());
}
