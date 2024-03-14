#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void QCD(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void penderiv(void *, void *, void *, void *, void *);

/* .Fortran calls */
extern void F77_NAME(rqbr)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rqfnb)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"QCD",      (DL_FUNC) &QCD,      11},
    {"penderiv", (DL_FUNC) &penderiv,  5},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"rqbr",  (DL_FUNC) &F77_NAME(rqbr),  27},
    {"rqfnb", (DL_FUNC) &F77_NAME(rqfnb), 13},
    {NULL, NULL, 0}
};

void R_init_rqPen(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
