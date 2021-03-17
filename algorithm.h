#include "bimatrix.h"

//#define eps 1e-5

equilibrium* lemke_howson_gen(double*** tableaus, double** bimatrix, int dim1, int dim2, int pivot, int *npassi, int debug);

eqlist* all_lemke_gen(double*** tableaus, double** bimatrix, int dim1, int dim2, int taboo, eqlist* , int debug);
