#include "bimatrix.h"

/*
  Gets a random dim1 x dim2 bimatrix. The values are drawn from a uniform random 
  distribution with mean zero and extremal values of -1.0 and +1.0. We seed the
  pseudo-random number generator with the system time with microseconds resolution:
  this is necessary because we need to avoid the risk of generating the same
  game in two different executions of the program.
*/

double **get_random_bimatrix_gen(int dim1, int dim2, double *min)
{
  int i;
  double n1, n2;

  double **bimatrix = (double **) malloc(sizeof(double *) * 2 * dim1);
  for (i = 0; i < (2 * dim1); i++)
    bimatrix[i] = (double *) malloc(sizeof(double) * dim2);

  /*
    As we said, we need to use srand48 instead of the usual srand() pseudo-random
    generator.
  */

  struct timeval tim;
  gettimeofday(&tim, NULL);
  srand48((long) (tim.tv_sec * 1000000 + tim.tv_usec));
  *min = 1000000;

  for (i = 0; i < (dim1 * dim2); i++) {
    n1 = 2.0 * drand48() - 1.0;
    n2 = 2.0 * drand48() - 1.0;

    bimatrix[i % dim1][i / dim1] = n1;
    bimatrix[i % dim1 + dim1][i / dim1] = n2;

    *min = *min < (n1 < n2 ? n1 : n2) ? *min : (n1 < n2 ? n1 : n2);
  }

  return bimatrix;
}

/*
  Creates (and allocates necessary memory) the two tableaus needed by the algorithm,
  starting from the bimatrix. 
*/

double*** create_systems(double** bimatrix, int dim1, int dim2) {  
  int i, j;
  
  double*** tableaus = (double***) malloc( 2 * sizeof(double**) );

  //Memory allocation for the two tableaus
  
  tableaus[0] = (double**) malloc( dim1 * sizeof(double*) );
  for(i = 0; i < dim1; i++) {
    tableaus[0][i] = (double*) calloc( (2 + dim1 + dim2), sizeof(double) );
  }
  tableaus[1] = (double**) malloc( dim2 * sizeof(double*) );
  for(i = 0; i < dim2; i++) {
    tableaus[1][i] = (double*) calloc( (2 + dim1 + dim2), sizeof(double) );
  }
  
  /*
    Initialization of the two tableaus. The first column represents the index of the variable,
    with the convention that a negative number represents the slack variable associated with
    the corresponding positive index variable. The second column is the actual first column of
    the tableau, and represent, during the execution of the algorithm, the value of the variable
    in basis for that row.
  */
  
  for (i = 0; i < dim1; i++) {
    tableaus[0][i][0] = - i - 1.0;
    tableaus[0][i][1] = 1.0;
  }
  for (i = 0; i < dim2; i++) {
    tableaus[1][i][0] = - i - dim1 - 1.0;
    tableaus[1][i][1] = 1.0;
  }

  /*
    We now only need to copy the bimatrix in the correct cells in the tableau.
  */
  for (i = 0; i < dim1; i++ ) {
    for (j = (2 + dim1); j<(dim1+dim2+2); j++) {
      tableaus[0][i][j] = - bimatrix[i][j - 2 - dim1];
    }
  }
  for (i = 0; i < dim2; i++) {
    for (j =  (2 + dim2); j<(dim1+dim2+2); j++) {
      tableaus[1][i][j] = - bimatrix[dim1 + ( j - 2 - dim2)][i];
    }
  }

  return tableaus;
}

void view_bimatrix_gen(double** bimatrix, int dim1, int dim2, FILE *f) {
  int i, j;

  fprintf(f,"Bimatrix following:\n\nPlayer A:\n");
  for(i = 0; i < dim1; i++) {
    fprintf(f,"\n");
    for( j = 0; j < dim2; j++) {
      fprintf(f,"%lf ",bimatrix[i][j]);
    }
  }

  fprintf(f,"\n\nPlayer B:\n");
  for(i = dim1; i < 2 * dim1; i++) {
    fprintf(f,"\n");
    for( j = 0; j < dim2; j++) {
      fprintf(f,"%lf ",bimatrix[i][j]);
    }
  }

  fprintf(f,"\n\n");
}

void view_tableau_gen(double** tableau, int dim1, int dim2, FILE *f) {
  int i, j;

  for( i = 0; i < dim1; i++ ) {
    fprintf(f,"\n");
    for( j = 0; j < (2 + dim1 + dim2); j++) {
      fprintf(f,"%lf ",tableau[i][j]);
    }
  }

  fprintf(f,"\n");
}

/*
  Due to some facts we assume in our implementation of the algorithm, we need all of the payoffs to be positive, so
  we simply add an offset to all payoffs, thus having them all > 0.
*/

void positivize_bimatrix(double** bimatrix, int dim1, int dim2, double minimo) {
  int i, j;

  for(i = 0; i < (2 * dim1); i++) {
    for(j = 0; j < dim2; j++) {
      bimatrix[i][j] -= (minimo - 1.0);
    }
  }
}

/*
  Tells if strategy 'strategy' is in base (looking at the current tableau). If not, it returns the same strategy, if it's in base,
  it returns the corresponding slack variable. This is needed by all_lemke, because after each execution of the LH algorithm
  we pivot on every variable from 1 to dim1+dim2, without knowing if that variable is in fact in base or not.
*/

int get_pivot_gen(double*** tableaus, int dim1, int dim2, int strategy) {
  int i;

  for(i = 0; i < dim1; i++) {
    if( tableaus[0][i][0] == strategy )
      return -strategy;
  }

  for(i = 0; i < dim2; i++) {
    if( tableaus[1][i][0] == strategy)
      return -strategy;
  }

  return strategy;
}

//Returns the tableau that contains the given strategy

int get_tableau(int dim1, int dim2, int strategy) {
  if ( strategy > dim1 || (strategy < 0 && strategy >= -dim1) )
    return 0;
  if ( strategy < -dim1 || (strategy > 0 && strategy <= dim1) ) 
    return 1;
  
  return -1;
}

//Returns the column that corresponds to the given strategy

int get_column(int dim1, int dim2, int strategy) {
  
  if( strategy > 0 && strategy <= dim1 ) {
    return (1 + dim2 + strategy );
  }
  if( strategy > 0 && strategy > dim1 ) {
    return (1 + dim1 + strategy - dim1);
  }
  if( strategy < 0 && strategy >= -dim1) {
    return (1 - strategy );
  }
  
  return ( 1 - strategy - dim1 );
}

/*
  This is the main function needed to import a normal form game in Gambit NFG format.
  In fact it accepts only one of the two normal form game formats: only files starting
  with the "NFG 1 D" identification string will be accepted by the program. This, by the
  way, is the format used by GAMUT to export games.
*/

double** gamut_import_bimatrix(FILE *f, double *minimo, int* rdim1, int* rdim2) 
{
  char *buf = (char *) malloc(100 * sizeof(char));
  char c;
  int i, j, tmpn;
  size_t num_bytes = 100;
  int dim1, dim2; 
  double n1, n2;

  /*
    This checks if we are parsing the correct file type.
  */
  getline(&buf,&num_bytes,f);
  if( strncmp(buf,"NFG 1 D",(size_t)7) != 0 ) {
    fprintf(stderr,"NFG file corrupted, aborting\n");
    exit(1);
  }

  /*
    First we need to ignore all comments (and player names) from the .nfg file.
  */
  
  tmpn = 0;
  while(tmpn < 2) { //This is to ignore comments
    c = fgetc(f);
    if( c == '\"' )
      tmpn++;
  }
  tmpn = 0;
  while(tmpn < 2) { //And this ignores player names
    c = fgetc(f);
    if( c == '{' || c == '}')
      tmpn++;
  }
  
  fscanf(f,"%d %d",&dim1,&dim2);
  *rdim1 = dim1; *rdim2 = dim2;
  fgetc(f); fgetc(f);

  /*
    We need now to allocate memory for the bimatrix, as we just determined the dimension of the game.
  */
  double **bimatrix = (double **) malloc(sizeof(double *) * 2 * dim1);
  for (i = 0; i < (2 * dim1); i++)
    bimatrix[i] = (double *) malloc(sizeof(double) * dim2);
  
  *minimo = 1000000;
  
  for(i = 0; i < dim2; i++) {
    for(j = 0; j < dim1; j++) {
      fscanf(f,"%lf %lf ",&n1,&n2);
      *minimo = *minimo < (n1 < n2 ? n1 : n2) ? *minimo : (n1 < n2 ? n1 : n2);
      bimatrix[j][i] = n1; bimatrix[j+dim1][i] = n2;
    }
  }

  free(buf);
  return bimatrix;
}


void free_tableaus(double*** tableaus, int dim1, int dim2) {
  int i;

  for(i=0; i < dim1; i++) {
    free(tableaus[0][i]);
  }
  free(tableaus[0]);

  for(i=0; i < dim2; i++) {
    free(tableaus[1][i]);
  }
  free(tableaus[1]);

  free(tableaus);
}

void free_bimatrix(double** bimatrix, int dim1, int dim2) {
  int i;

  for(i=0; i < (2*dim1); i++) {
    free(bimatrix[i]);
  }
  free(bimatrix);
}
