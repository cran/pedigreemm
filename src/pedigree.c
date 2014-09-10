#include "pedigree.h"

#ifdef ENABLE_NLS		/** Allow for translation of error messages */
#include <libintl.h>
#define _(String) dgettext ("pedigreemm", String)
#else
#define _(String) (String)
#endif

/* When appropriate, alloca is cleaner than malloc/free.  The storage
 * is freed automatically on return from a function. When using gcc the
 * builtin version is much faster. */
#ifdef __GNUC__
# undef alloca
# define alloca(x) __builtin_alloca((x))
#else
/* this is necessary (and sufficient) for Solaris 10: */
# ifdef __sun
#  include <alloca.h>
# endif
#endif

/** alloca n elements of type t */
#define Alloca(n, t)   (t *) alloca( (size_t) ( (n) * sizeof(t) ) )

/* Functions for the pedigree class */

/**
 * Create the left Cholesky factor of the numerator relationship
 * matrix from a pedigree.
 *
 * @param x a pedigree object
 * @param ans T stored as a non-unit, lower dtCMatrix
 *
 * @return ans with elements modified to incorporate D
 */
SEXP pedigree_chol(SEXP x, SEXP ans)
{
    SEXP Sire = GET_SLOT(x, install("sire"));
    int *ai = INTEGER(GET_SLOT(ans, install("i"))),
	*ap = INTEGER(GET_SLOT(ans, install("p"))),
	*dam = INTEGER(GET_SLOT(x, install("dam"))),
	*sire = INTEGER(Sire), j, n = LENGTH(Sire);
    double *ax = REAL(GET_SLOT(ans, install("x"))), *F, Di, tmp;

    setAttrib(ans, install("F"), allocVector(REALSXP, n));
    F = REAL(getAttrib(ans, install("F")));
    for (int i = 0; i < n; i++) {
	int p = sire[i] - 1, q = dam[i] - 1;
	if (sire[i] == NA_INTEGER) {
	    F[i] = 0;
	    Di = (dam[i] == NA_INTEGER) ? 1 : sqrt(0.75 - 0.25 * F[q]);
	} else {
	    if (dam[i] == NA_INTEGER) { /* sire only */
		F[i] = 0;
		Di = sqrt(0.75 - 0.25 * F[p]);
	    } else {		/* both parents in pedigree */
		Di = sqrt(0.5 - 0.25 * (F[p] + F[q]));
		F[i] = NA_REAL;
		if ((ap[i + 1] - ap[i]) > 1) {	  /* skip if no progeny */
		    if (p > q) {j = p; p = q; q = j;} /* ensure p <= q */
		    F[i] = 0;
		    for (int j = 0; j <= p; j++) {
			tmp = 0;
			for (int k = ap[j]; k < ap[j + 1] && ai[k] <= q; k++) {
			    int ii = ai[k];
			    if (ii == p) tmp = ax[k];
			    if (ii == q) F[i] += tmp * ax[k]/2;
			}
		    }
		}
	    }
	}
	for (int j = ap[i]; j < ap[i + 1]; j++) ax[j] *= Di;
    }
    return ans;
}

/**
 * Create the inbreeding coefficients according to the algorithm given
 * in "Comparison of four direct algorithms for computing inbreeding
 * coefficients" by Mehdi Sargolzaei and Hiroaki Iwaisaki, Animal
 * Science Journal (2005) 76, 401--406.  This function is a modified
 * version of the code published in an appendix to that paper.
 *
 * @param x a pedigree object
 *
 * @return a list of the inbreeding coefficients
 */

SEXP pedigree_inbreeding(SEXP x)
{
    SEXP ans, sp = GET_SLOT(x, install("sire"));
    int i, j, t, n = LENGTH(sp), S, D;
    int *SI, *MI,		/* start and minor */
	*sire = INTEGER(sp),
	*dam = INTEGER(GET_SLOT(x, install("dam")));
    double *F = Calloc(n + 1, double), /* inbreeding coefficients */
      *L = Calloc(n + 1, double), *B = Calloc(n + 1, double);
    int *Anc = Calloc(n + 1, int),	/* ancestor */
	*LAP = Calloc(n + 1, int); 	/* longest ancestoral path */
    R_CheckStack();
    
    for (i = 0; i < n; i++) {     /* Replace NA's by zeros */
	if (sire[i] == NA_INTEGER) sire[i] = 0;
	if (dam[i] == NA_INTEGER) dam[i] = 0;
    }
    F[0] =-1; LAP[0] =-1; /* set F and lap for unknown parents */
    for(i = 1, t = -1; i <= n; i++) { 	/* evaluate LAP and its maximum */
	S = sire[i-1]; D = dam[i-1]; /* parents of animal i-1 */
	LAP[i] = ((LAP[S] < LAP[D]) ? LAP[D] : LAP[S]) + 1;
	if (LAP[i] > t) t = LAP[i];
    }
    SI = Calloc(t + 1, int);
    MI = Calloc(t + 1, int);
    for(i = 0; i <= t ; ++i) SI[i] = MI[i] = 0; /* initialize start and minor */
    for(i = 1; i <= n; i++) { 	/* evaluate F */
	S = sire[i-1]; D = dam[i-1]; /* parents of animal i */
	B[i] = 0.5 - 0.25 * (F[S] + F[D]); 
				/* adjust start and minor */
	for (j = 0; j < LAP[i]; j++) {++SI[j]; ++MI[j];} 
	if (S == 0 || D == 0) { /* both parents unknown */
	    F[i] = L[i] = 0; continue;
	}
	if(S == sire[i-2] && D == dam[i-2]) { /* full-sib with last animal */
	    F[i] = F[i-1]; L[i] = L[i-1]; continue;
	}
    
	F[i] = -1; L[i] = 1; 
	t = LAP[i]; /* largest lap group number in the animal's pedigree */
	Anc[MI[t]++] = i; /* initialize Anc and increment MI[t] */
	while(t > -1) { /* from the largest lap group to zero */
	    j = Anc[--MI[t]]; /* next ancestor */
	    S = sire[j-1]; D = dam[j-1]; /* parents of the ancestor */
	    if (S) {
		if (!L[S]) Anc[MI[LAP[S]]++] = S; 
				/* add sire in its lap group in Anc
				 * array if it is not added yet and
				 * increment the minor index for the group */ 
		L[S] += 0.5 * L[j]; /* contribution to sire */
	    }
	    if (D) {
		if (!L[D]) Anc[MI[LAP[D]]++] = D;
		L[D] += 0.5 * L[j]; /* contribution to dam */
	    }
	    F[i] += L[j] * L[j] * B[j];
	    L[j] = 0; /*clear L[j] for the evaluation of the next animal */
	    if (MI[t] == SI[t]) --t; /* move to the next lap group when
				      * all ancestors in group t have been
				      * evaluated */
	} 
    }
    ans = PROTECT(allocVector(REALSXP, n));
    Memcpy(REAL(ans), F + 1, n);

    for (i = 0; i < n; i++) {     /* Restore the NA's */
	if (!sire[i]) sire[i] = NA_INTEGER;
	if (!dam[i]) dam[i] = NA_INTEGER;
    }
    Free(F); Free(L); Free(B); Free(Anc); Free(LAP); Free(SI); Free(MI);
    UNPROTECT(1);
    return ans;
}

/**
 * This is a creplacement for the former R-function "getGenAncestors"
 * It takes integer pointers to the sire, dam, id and generation vectors
 * and it changes the generations directly on the passed R-object, so there is no
 * return value = void.
 */

void calc_generation(int* sire, int* dam, int* id, int* gene, int this_id) {

    int j = this_id;
    int gene_sire, gene_dam;
    
/* If there are no parents assign 0 */
    gene[j] = 0;

/* check first parent */ 
    if(sire[j] != NA_INTEGER) {

      if (gene[sire[j]-1] == NA_INTEGER) {

        calc_generation(sire, dam, id, gene, sire[j]-1);
        gene_sire = 1 + gene[sire[j]-1];

      } else {

          gene_sire = 1 + gene[sire[j]-1];      

        } 
        
      gene[j] = gene_sire;

    }

/* check second parent */ 
    if(dam[j] != NA_INTEGER) {

      if (gene[dam[j]-1] == NA_INTEGER) {

        calc_generation(sire, dam, id, gene, dam[j]-1);
        gene_dam = 1 + gene[dam[j]-1];

      } else {

          gene_dam = 1 + gene[dam[j]-1];      

        }
        
      if (gene_dam > gene[j])  gene[j] = gene_dam;  

    }

}

/**
 * This function iterates over every element of the pedigree and let
 * "calc_generation" calculate the generations. Again without a return value.
 */


void get_generation(SEXP sire_in, SEXP dam_in, SEXP id_in, SEXP gene_in, SEXP verbose_in) {

    int *sire = INTEGER(sire_in),
	*dam = INTEGER(dam_in),
	*id = INTEGER(id_in),
	*gene = INTEGER(gene_in);
 
    int verbose = *INTEGER(verbose_in);
    int n = LENGTH(sire_in);

    for(int i=0; i<n; i++) {
      
      if(verbose) Rprintf("%i \n", i); 

      if(gene[i] == NA_INTEGER) {

        calc_generation(sire, dam, id, gene, i);

      }

    }

}

