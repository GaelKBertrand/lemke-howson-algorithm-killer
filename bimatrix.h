#define _GNU_SOURCE
#include <sys/time.h>

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "equilibria.h"

#define eps 1e-20

typedef struct sh_tab {
  int row;
  int label;
  struct sh_tab* next;
} sh_tableau;


//Imports a bimatrix from a NFG file.
double** gamut_import_bimatrix(FILE *, double *min, int* rdim1, int* rdim2);

//Gets a uniformely random dim1xdim2 bimatrix.
double** get_random_bimatrix_gen(int dim1, int dim2, double *);

//Creates the tableaus starting from the bimatrix
double*** create_systems(double** bimatrix,int dim1, int dim2);

//Adds an offset to all payoffs to have them positive
void positivize_bimatrix(double** bimatrix,int dim1, int dim2, double min);

//Debug output
void view_bimatrix_gen(double**,int dim1, int dim2,FILE*);
void view_tableau_gen(double**,int dim1, int dim2, FILE*);

//Creates a copy of the system
double** system_copy(double**,int);

//Tells if strategy 'strategy' is in the current tableau's base.
int get_pivot_gen(double*** tableaus,int dim1, int dim2, int strategy);

//Returns the tableau in wich the strategy is contained
int get_tableau(int dim1, int dim2, int strategy);

//Returns the column that corresponds to the given strategy
int get_column(int dim1, int dim2, int strategy);

//Memory managment functions
void free_tableaus(double*** tableaus, int dim1, int dim2);
void free_bimatrix(double** bimatrix, int dim1, int dim2);
