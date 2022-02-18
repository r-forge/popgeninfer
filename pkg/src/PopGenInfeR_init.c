/* Generated using tools::package_native_routine_registration_skeleton(".")
   see: https://stackoverflow.com/questions/42313373/r-cmd-check-note-found-no-calls-to-r-registerroutines-r-usedynamicsymbols
*/
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void SimulateCoalescentTree(void *, void *, void *, void *, void *, void *, void *);
extern void SimulateIslandModel(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"SimulateCoalescentTree", (DL_FUNC) &SimulateCoalescentTree, 7},
    {"SimulateIslandModel",    (DL_FUNC) &SimulateIslandModel,    6},
    {NULL, NULL, 0}
};

void R_init_PopGenInfeR(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
