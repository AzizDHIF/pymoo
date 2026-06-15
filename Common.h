#ifndef Common_h
#define Common_h

typedef struct {

  int *items_nonpris;
 /* int *items_pris;*/
  int nombr_nonpris;
  int nombr;
    float fitness;
	int explored;
	double* v;
	double* f;
	int* d;


}ind;

typedef struct{
  float min;
  float max;
}range;

#endif
