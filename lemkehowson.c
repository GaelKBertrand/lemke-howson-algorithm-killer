#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <assert.h>
#include <sys/time.h>

#include "algorithm.h"

void single_lemke_exec();
void all_lemke_exec();

int main(int argc, char **argv)
{
  FILE *input;
  int c;
  int sing_l = 0, all_l = 0, readgame = 0, debug_mask = 0, gambit_output = 0, summary = 0;
  double** bimatrix;
  char inputfile[100];
  int startpivot = 1;
  double minimo = 0.0;
  int dim1 = 10, dim2 = 10;

  while ((c = getopt(argc, argv, "p:i:w:l:d:Ghas")) != -1) {
    switch (c) {
    case 'p':
      sing_l = 1;
      startpivot = atoi(optarg);
      break;
    case 'i':
      snprintf(inputfile, 100, "%s", optarg);
      readgame = 1;
      break;
    case 'w':
      dim1 = atoi(optarg);
      break;
    case 'l':
      dim2 = atoi(optarg);
      break;
    case 'a':
      all_l = 1;
      break;
    case 'd':
      debug_mask = atoi(optarg);
      break;
    case 'G':
      gambit_output = 1;
      break;
    case 's':
      summary = 1;
      break;
    case 'h':
      fprintf(stderr, "Usage: ./lemkehowson\n\t\t\t[-i gamefile.NFG (by default generates a random game)]\n\t\t\t[-w DIM1 -l DIM2 (used only to generate a random game of size DIM1xDIM2. Default is 10 x 10)]\n\t\t\t[-p PIVOT (Executes the Lemke-Howson algorithm once, pivoting on strategy PIVOT)]\n\t\t\t[-a (Searches all equilibria reachable by the Lemke-Howson algorithm)]\n\t\t\t[-d DEBUG_LEVEL (Determines the level of debug output)]\n\t\t\t[-G (With this option turned on, the output is similar to that of Gambit, to semplify testing and benchmarking)]\n");
      return 0;
      break;
    default:
      return -1;
    }
  }

/*
  If we don't read the game from a file, by default we generate a uniformely random game
*/

  if (!readgame) {
    bimatrix = (double**) get_random_bimatrix_gen(dim1,dim2, &minimo);
  } 
  else {
    input = fopen(inputfile, "r");
    bimatrix = gamut_import_bimatrix(input, &minimo, &dim1, &dim2);
    fclose(input);
  }

  if( sing_l && all_l ) {
    fprintf(stderr,"You must choose whether to look for a single equilibrium with the Lemke-Howson algorithm with [-p PIVOT] or to have a list of all equilibria reachable by Lemke-Howson (with [-a])\n");
    exit(1);
  }
  
  if( sing_l ) {
    single_lemke_exec(bimatrix,dim1,dim2,startpivot,minimo,gambit_output,summary,debug_mask);
  }
  else if( all_l ) {
    all_lemke_exec(bimatrix,dim1,dim2,minimo,gambit_output,debug_mask);
  }

  return 0;
}

/*
  This way the program executes the Lemke-Howson algorithm one single time, pivoting on the desired strategy,
  on the game specified (it can be a random game or a game imported from a NFG file).
*/

void single_lemke_exec(double** bimatrix, int dim1, int dim2, int pivot, double min, int gambit_output, int summary, int debug_mask) {
  int passi;
  double*** tableaus;

  if( pivot <= 0 || pivot > (dim1+dim2) ) {
    fprintf(stderr,"Starting pivot must be a number between 1 and DIM1 + DIM2\n");
    exit(1);
  }

  positivize_bimatrix(bimatrix,dim1,dim2,min);

  tableaus = create_systems(bimatrix,dim1,dim2);
  
  equilibrium* eq = lemke_howson_gen(tableaus,bimatrix,dim1,dim2,pivot,&passi,debug_mask);

  if(summary) {
    fprintf(stdout,"%d %d\n",passi,eq_size(eq));
  }
  else if(gambit_output) {
    print_equilibrium_gambit(eq,dim1,dim2,stdout);
  }
  else {
    print_equilibrium(eq,stdout);
  }

  if(!summary)
    fprintf(stdout,"Number of complementary pivoting steps performed by the algorithm: %d\n",passi);

  free_equilibrium(eq);
  free_tableaus(tableaus,dim1,dim2);
  free_bimatrix(bimatrix,dim1,dim2);
}

/*
  This way the program enumerates all equilibria reachable by the Lemke-Howson algorithm. This is done
  by recursively executing LH algorithm on each equilibrium found, breaking the recursion when we reach
  an equilibrium we already found before.
*/

void all_lemke_exec(double** bimatrix, int dim1, int dim2, double min, int gambit_output, int debug_mask) {
  double*** tableaus;
  eqlist* found_equilibria;

  positivize_bimatrix(bimatrix,dim1,dim2,min);
  
  tableaus = create_systems(bimatrix,dim1,dim2);

  found_equilibria = all_lemke_gen(tableaus,bimatrix,dim1,dim2,-1,(eqlist*)0,debug_mask);
  
  if(gambit_output) {
    print_eqlist_gambit(found_equilibria,dim1,dim2,stdout);
  }
  else {
    print_eqlist(found_equilibria,stdout);
  }

  free_tableaus(tableaus,dim1,dim2);
  free_bimatrix(bimatrix,dim1,dim2);
  free_eqlist(found_equilibria);
}
