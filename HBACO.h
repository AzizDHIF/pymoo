#ifndef HBACO_h
#define HBACO_h

#include "Common.h"




/*typedef struct {
int *items_nonpris;
  int nombr_nonpris;
  int nombr;float fitness;
	int explored;
	double* f;
	int* d;
}ind;
*/
typedef struct pop_st  /* a population */
{
    int size;
    int maxsize;
    ind **ind_array;
} pop;

typedef struct genetic_op /* a genetic operator (2 point mutations)*/
{
  int p1;
  int p2;
} mut;

double pow(double, double);
double sqrt(double);
ind* ind_copy(ind *i);
void* chk_malloc(size_t size);
pop* create_pop(int maxsize, int dim);
ind* create_ind(int dim);
void complete_free_pop(pop *pp);
int dominates(ind *p_ind_a, ind *p_ind_b);
double random_nb(double,double);
void mutate(ind *x,mut *m);

int max(int a, int b);

#endif
