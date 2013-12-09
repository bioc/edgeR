#include "utils.h"
#ifndef MATVEC_CHECK_H
#define MATVEC_CHECK_H

class matvec_check {
public:
	matvec_check(const int, const int, SEXP, const bool, const char*, const bool);
	~matvec_check();
	void advance();
	const double* const* access() const;
private:
	const double* mycheck;
	double* temp;
	bool isvec, istran;
	const int nl, nt;
	int tagdex, libdex;
};

#endif
