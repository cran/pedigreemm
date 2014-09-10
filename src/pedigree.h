#ifndef PEDIGREE_H
#define PEDIGREE_H

#include <R.h>
/* Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
 * GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc. */
#include <Rdefines.h>

SEXP pedigree_chol(SEXP x, SEXP ans);
SEXP pedigree_inbreeding(SEXP x);
void get_generation(SEXP sire_in, SEXP dam_in, SEXP id_in, SEXP gene_in, SEXP verbose_in);

#endif /* PEDIGREE_H */
