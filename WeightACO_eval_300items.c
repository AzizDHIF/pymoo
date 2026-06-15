#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <conio.h>

#include "Common.h"
#include "HBACO.h"


#define dimension 3
#define NBITEMS 300
/*#define NBANTS 100*/
#define FREQUANCY 80


#ifndef M_PI
 #define M_PI 3.14159265358979323846
#endif

double capacities[dimension];
double weights[dimension][NBITEMS];
double profits[dimension][NBITEMS];

int nf,ni,cardP;
int nombr;
int NBi=250+100;
int nbgeneration=50;
/* MOACO parameters */
/*double meilprofit[dimension];*/
/*double pheromone [NBITEMS];*/
int paretoIni=15000;
int nbants=20;
double rhot=0.90;
int alphat=1;
/*const int nbants=100;*/
int maxcycle=802;
/*double tmin=0.5;*/
double tmax=1.0 , beta=10.0;
int tinit=1;
double referencePoint[dimension];

 ///////////*************


int bruit_rate=10;
double pheromone[NBITEMS];
double eta[dimension];
int inter;
int iseed;

//////////************


/*int alpha=10;*/



pop *archive=NULL;
pop *solutions=NULL;
pop *LSarchive=NULL;
pop *P=NULL;

float rho ;  /* determines the reference point for the hypervolume indicator */
FILE *Wfile;
double vector_weight[dimension];
int dim;
int problem_size;
double OBJ_Weights[dimension][10000];
int nombreLIGNE=0;
int nextLn=0;
int inv=0;

void loadMOKP(char *s){

  FILE* source;

int i,f;

  char cl[20];

  // Opening
   source=fopen(s,"r");

  fscanf(source, " %d %d  \n", &nf,&ni);
  printf( " %d %d  \n ", nf,ni);

  for (f=0;f<nf;f++)
	{	fscanf(source, "%lf  \n ",&capacities[f]);
		/*printf( " %lf  \n \n", capacities[f]);*/

		for (i=0;i<ni;i++)
		{
			fscanf(source, " %s \n", cl); /*printf( " %s  \n \n", cl);*/
			fscanf(source, " %lf  \n", &weights[f][i]); /*printf( " %d \n \n", weights[f][i]);*/
			fscanf(source, "  %lf  \n ", &profits[f][i]); /*printf( " %d \n \n", profits[f][i]);*/

		}
	}
  fclose(source);
}


  /********************************/
 /* Begin memory (des)allocation */
/********************************/

void* chk_malloc(size_t size)
/* Wrapper function for malloc(). Checks for failed allocations. */
{
    void *return_value = malloc(size);
    if(return_value == NULL)
	printf("Selector: Out of memory.");
    return (return_value);
}





pop* create_pop(int maxsize, int nf)
/* Allocates memory for a population. */
{
    int i;
    pop *pp;

    assert(nf >= 0);
    assert(maxsize >= 0);

    pp = (pop*) chk_malloc(sizeof(pop));
    pp->size = 0;
    pp->maxsize = maxsize;
    pp->ind_array = (ind**) chk_malloc(maxsize * sizeof(ind*));

    for (i = 0; i < maxsize; i++)
	pp->ind_array[i] = NULL;

    return (pp);
}


ind* create_ind(int nf)
/* Allocates memory for one individual. */
{
    ind *p_ind;

    assert(nf >= 0);
    assert(ni >=0 );
    assert(NBi >=0 );



    p_ind = (ind*) chk_malloc(sizeof(ind));

    p_ind->items_nonpris=(int*) chk_malloc(ni * sizeof(int));
    /*p_ind->items_pris=(int*) chk_malloc(ni * sizeof(int));*/
    p_ind->nombr_nonpris=0;
    p_ind->nombr=0;
    p_ind->fitness = -1.0;
    p_ind->explored = 0;
    p_ind->f = (double*) chk_malloc(nf * sizeof(double));
    p_ind->v = (double*) chk_malloc(nf * sizeof(double));
    p_ind->d = (int*) chk_malloc(ni * sizeof(int));




    return (p_ind);
}

ind* ind_copy(ind *i){
/* Allocates memory for one individual. */
    ind *p_ind=NULL;
    int k;

    p_ind=create_ind(nf);

    p_ind->nombr_nonpris=i->nombr_nonpris;
    p_ind->nombr=i->nombr;

    /*for (k=0;k<i->nombr;k++)
    {p_ind->items_pris[k]=i->items_pris[k];}*/


    for (k=0;k<i->nombr_nonpris;k++)
    {p_ind->items_nonpris[k]=i->items_nonpris[k];}

    p_ind->fitness = i->fitness;
    p_ind->explored = i->explored;

    for (k=0;k<nf;k++) {p_ind->f[k]=i->f[k];p_ind->v[k]=i->v[k];}

    for (k=0;k<i->nombr;k++) {p_ind->d[k]=i->d[k];}

    return (p_ind);
}

void free_ind(ind *p_ind)
/* Frees memory for given individual. */
{
  assert(p_ind != NULL);
  free(p_ind->items_nonpris);
  free(p_ind->d);
  free(p_ind->v);
  free(p_ind->f);
  free(p_ind);
}

void free_pop(pop *pp)
/* Frees memory for given population. */
{
   if (pp != NULL)
   {
      free(pp->ind_array);
      free(pp);
   }
}

void complete_free_pop(pop *pp)
/* Frees memory for given population and for all individuals in the
   population. */
{
   int i = 0;
   if (pp != NULL)
   {
      if(pp->ind_array != NULL)
      {
         for (i = 0; i < pp->size; i++)
         {
            if (pp->ind_array[i] != NULL)
            {
               free_ind(pp->ind_array[i]);
               pp->ind_array[i] = NULL;
            }
         }
/*printf( " for" );*/
         free(pp->ind_array);
/*printf( " free ind" );*/
      }

      free(pp);
   }
   /*printf( " free " );*/
}



int dominates(ind *p_ind_a, ind *p_ind_b)
/* Determines if one individual dominates another.
   Minimizing fitness values. */
{
    int i;
    int a_is_worse = 0;
    int equal = 1;

     for (i = 0; i < nf && !a_is_worse; i++)
     {
	 a_is_worse = p_ind_a->f[i] < p_ind_b->f[i];
          equal = (p_ind_a->f[i] == p_ind_b->f[i]) && equal;
     }

     return (!equal && !a_is_worse);
}

int non_dominated(ind *p_ind_a, ind *p_ind_b)
/* Determines if one individual is non dominated according to another one.
   Minimizing fitness values.
   -1 dominated, 0 equal, 1 non-dominated */
{
    int i;
    int a_is_good = -1;
    int equal = 1;

     for (i = 0; i < nf; i++)
     {
       if (p_ind_a->f[i] > p_ind_b->f[i]) a_is_good=1;
         equal = (p_ind_a->f[i] == p_ind_b->f[i]) && equal;
     }
     if (equal) return 0;
     return a_is_good;
}


int otherResult(char* file_name,  double mpareto2[][nf]) {
  FILE* fd;

  int f;
  int cardP2 =0;

   if ( (fd=fopen(file_name, "r"))==NULL)
    {	printf("ERREUR: Verifiez le nom de fichier %s", file_name); return 0;	}

  do{

  for (f=0;f<nf;f++)
	{

	  fscanf(fd, "%lf  ",&mpareto2[cardP2][f]);

	}
  cardP2++;
  } while(!feof(fd));

  return cardP2;

}
// fonction de calcul de performance
double Cmesure(double mpareto1[][nf], double mpareto2[][nf], int cardP1, int cardP2)
{
	int i=0,efficace=1,existe=0,f, cdom,cegal,c=0,j;

		printf("carp1 %d cardp2 %d \n ", cardP1, cardP2);
				while (i<cardP2)
				{j=0;
					while (j<cardP1)
					{
					f=0;
					cegal=0;
					cdom=0; //pour detecter si sol cour est dominée
					while (f<nf)
					{
						if (mpareto2[i][f]<mpareto1[j][f])
						{
						cdom++;
						}
						if (mpareto2[i][f]==mpareto1[j][f])
						{
						cegal++;
						}


						f++;
					}
					j++;
					}
					if ((cdom>0)&&((cdom+cegal)==nf))//la sol cou est dominé
						c++;
					i++;
				}

				return ((double) c/ (double) cardP2);

}


 /******************************/
 /* WEIGHTS functions */
/******************************/


void initfile_weights_log()
{
    int F=FREQUANCY/4,i,j,k;
    double lamda1,lamda2,lamda3,lamda4,tmp1,tmp2,tmp3;
    int T1,T2,M1,Zd;
    double C=exp(1);

  Wfile = fopen( "Weights.txt", "a+" );
  fflush(stdout);

if (nf==2)
{
for(i=0;i<=F;i++)
{

lamda1=log((4.0*i/FREQUANCY*C)+cos(2.0*M_PI*i/FREQUANCY));
lamda2=1.0-lamda1;

fprintf(Wfile,"%lf ",lamda1);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda2);
fflush(stdout);
fprintf(Wfile,"\n");
}
/*w1=1.0; w2=0.0;
fprintf(Wfile,"%lf ",w1);
fflush(stdout);
fprintf(Wfile,"%lf ",w2);
fflush(stdout);*/
}

if (nf==3)
{

for(i=0;i<F;i++)
{
lamda1=log((4.0*i/FREQUANCY*C)+(cos(2.0*M_PI*i/FREQUANCY)));

for(j=0;j<F;j++)
{

lamda2=(1.0-lamda1)*log((4.0*j/FREQUANCY*C)+(cos(2.0*M_PI*j/FREQUANCY)));

lamda3=1.0-lamda1-lamda2;

fprintf(Wfile,"%lf ",lamda1);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda2);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda3);
fflush(stdout);
fprintf(Wfile,"\n");
}
}
lamda1=1.0; lamda2=0.0; lamda3=0.0;
fprintf(Wfile,"%lf ",lamda1);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda2);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda3);
fflush(stdout);
}
if (nf==4)
{
for(i=0;i<F;i++)
{
lamda1=log((4.0*i/FREQUANCY*C)+(cos(2.0*M_PI*i/FREQUANCY)));

for(j=0;j<F;j++)
{
lamda2=(1.0-lamda1)*log((4.0*j/FREQUANCY*C)+(cos(2.0*M_PI*j/FREQUANCY)));

for(k=0;k<F;k++)
{
lamda3=(1.0-lamda1-lamda2)*log((4.0*k/FREQUANCY*C)+(cos(2.0*M_PI*k/FREQUANCY)));
lamda4=1.0-lamda1-lamda2-lamda3;

fprintf(Wfile,"%lf ",lamda1);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda2);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda3);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda4);
fflush(stdout);
fprintf(Wfile,"\n");
}
}
}
lamda1=1.0; lamda2=0.0; lamda3=0.0; lamda4=0.0;
fprintf(Wfile,"%lf ",lamda1);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda2);
fflush(stdout);
fprintf(Wfile,"%f ",lamda3);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda4);
fflush(stdout);
}
fclose(Wfile);
}

void initfile_weights()
{
    int F=FREQUANCY/4,t1,t2,t3;
    int FQ=FREQUANCY/2;
    double w1,w2,w3,w4;



  Wfile = fopen( "Weights.txt", "a+" );
  fflush(stdout);

if (nf==2)
{
for(t1=0;t1<F;t1++)
{
w1=fabs(sin(2*M_PI*t1/FREQUANCY));
w2=1.0-w1;

fprintf(Wfile,"%lf ",w1);
fflush(stdout);
fprintf(Wfile,"%lf ",w2);
fflush(stdout);
fprintf(Wfile,"\n");
}
w1=1.0; w2=0.0;
fprintf(Wfile,"%lf ",w1);
fflush(stdout);
fprintf(Wfile,"%lf ",w2);
fflush(stdout);
}

if (nf==3)
{
for(t1=0;t1<F;t1++)
{
w1=fabs(sin(2*M_PI*t1/FREQUANCY));
for(t2=0;t2<F;t2++)
{
w2=(1.0-w1)*fabs(sin(2*M_PI*t2/FREQUANCY));
w3=1.0-w1-w2;
fprintf(Wfile,"%lf ",w1);
		fflush(stdout);
fprintf(Wfile,"%lf ",w2);
fflush(stdout);
fprintf(Wfile,"%lf ",w3);
fflush(stdout);
fprintf(Wfile,"\n");
}
}
w1=1.0; w2=0.0; w3=0.0;
fprintf(Wfile,"%lf ",w1);
fflush(stdout);
fprintf(Wfile,"%lf ",w2);
fflush(stdout);
fprintf(Wfile,"%lf ",w3);
fflush(stdout);
}
if (nf==4)
{
for(t1=0;t1<F;t1++)
{
w1=fabs(sin(2*M_PI*t1/FREQUANCY));
for(t2=0;t2<F;t2++)
{
w2=(1.0-w1)*fabs(sin(2*M_PI*t2/FREQUANCY));
for(t3=0;t3<F;t3++)
{
w3=(1.0-w1-w2)*fabs(sin(2*M_PI*t3/FREQUANCY));
w4=1.0-w1-w2-w3;
fprintf(Wfile,"%lf ",w1);
fflush(stdout);
fprintf(Wfile,"%lf ",w2);
fflush(stdout);
fprintf(Wfile,"%lf ",w3);
fflush(stdout);
fprintf(Wfile,"%lf ",w4);
fflush(stdout);
fprintf(Wfile,"\n");
}
}
}
w1=1.0; w2=0.0; w3=0.0; w4=0.0;
fprintf(Wfile,"%lf ",w1);
fflush(stdout);
fprintf(Wfile,"%lf ",w2);
fflush(stdout);
fprintf(Wfile,"%f ",w3);
fflush(stdout);
fprintf(Wfile,"%lf ",w4);
fflush(stdout);
}
fclose(Wfile);
}

 void read_weights_file (char *s)
 {
     int i,j,caractereLu;
      char cl[20];

      Wfile = fopen( s, "r" );
/*
    do
    {
        caractereLu = fgetc(Wfile);
        printf("%c", caractereLu);
        if (caractereLu == '\n')
            nombreLIGNE++;
    } while(caractereLu != EOF);

       rewind(Wfile);
*/
if (nf==2) {nombreLIGNE=(FREQUANCY/4)+1;}
if (nf==3) {nombreLIGNE=(FREQUANCY/4)*(FREQUANCY/4)+1;}
if (nf==4) {nombreLIGNE=FREQUANCY/4*(FREQUANCY/4)*(FREQUANCY/4)+1;}
printf("nombre ligne %d", nombreLIGNE);
		for (i=0;i<nombreLIGNE;i++)
		{
		    for (j=0;j<nf;j++)
	     {
			fscanf(Wfile, " %lf ", &OBJ_Weights[j][i]); /*printf( " obj_W %lf ", OBJ_Weights[j][i]);*/
		 }
          /*printf("%d",i); printf("\n");*/
	    }

  fclose(Wfile);

 }

double random_nb(double a,double b)
{
   double j= ( rand()/(double)RAND_MAX ) * (b-a) + a;
    return (j);
}





void dynamic_weight_allpop2()
{
    int i;

/*printf("nextln %d",nextLn);*/

if (inv==0){
for(i=0;i<dimension;i++)
{
    vector_weight[i]=OBJ_Weights[i][nextLn];
    /*printf("poids %lf",vector_weight[i]);*/
}
  nextLn++;
}

if (inv==1){
for(i=0;i<dimension;i++)
{
    vector_weight[i]=OBJ_Weights[i][nextLn];
    /*printf("poids %lf",vector_weight[i]);*/
}
nextLn--;
}


if(nextLn==nombreLIGNE)
        {inv=1;}



}


 /******************************/
 /* MOACO functions */
/******************************/
void seed(unsigned int seed)
 {iseed= seed; }


/**************************/
double randFloat(void) {
    /* Uniform random number generator x(n+1)= a*x(n) mod c
    with a = pow(7,5) and c = pow(2,31)-1.
	Copyright (c) Tao Pang 1997. */
   const int ia=16807,ic=2147483647,iq=127773,ir=2836;
    int il,ih,it;
    double rc;
    //extern int iseed;
    ih = iseed/iq;
    il = iseed%iq;
    it = ia*il-ir*ih;
    if (it > 0)	iseed = it;
    else iseed = ic+it;
    rc = ic;
    return(iseed/rc);
}

int choose(double *p, int nbCand){
  /* for i in 0..nbCand-1], p[i] = sum_{j<i} tau[j]^alpha * nu[j]^beta */
  /* returns k with probability (p[k]-p[k-1])/p[nbCand-1] */
  // exemple si on a 5 candidats et candidat1->10, 2->7, 3->3, 4->14, 5->8
  // alors p[0]=10, p[1]=17, p[2]=20, p[3]=34, p[4]=42
  // nbCand = 5

  double f = randFloat();
  int left=0;
  int right=nbCand-1;
  int k;
  double total=p[nbCand-1];
  while (left<right){
    k=(left+right+1)/2;
    if (f<p[k-1]/total) right=k-1;
    else if (f>p[k]/total) left=k+1;
    else return k;
  }
  if (left>=0 && left<nbCand) return left;
  else printf("pb choose\n");
  return k;
}

double monPow(double x, int y){
  double p;
  if (y==0) return 1;
  else if (y==1) return x;
  else if (y==2) return x*x;
  else if (y==3) return x*x*x;
  else if (y==4) {x *= x; return x*x;}
  else if (y==5) {p = x*x; return p*p*x;}
  else if (y==6){x = x*x*x; return x*x;}
  else {
    if (y%2==0){p = monPow(x,y/2); return p*p;}
    else{p = monPow(x,(y-1)/2); return p*p*x;};
  };
}






int irand(int range)
/* Generate a random integer. */
{
    int j;
    j=(int) ((double)range * (double) rand() / (RAND_MAX+1.0));
    return (j);

}
//****fonction pour calcul fitness



int extractPtoArchive(pop *Pareto, pop *Sarchive){
  int i,j,dom;
  int t=Sarchive->size+Pareto->size;
  pop *archiveAndP;
  int convergence_rate=0;


  archiveAndP=create_pop(t,nf);
  for (i=0;i<Sarchive->size;i++){
         /*printf("A %d",Sarchive->size);*/
    archiveAndP->ind_array[i]=Sarchive->ind_array[i]; /* do not copy */
    //archive->ind_array[i]=NULL;  /*!*/
  }

  //archive->size=0;
  for (i=0;i<Pareto->size;i++){
    /*printf("P %d",SP->size);*/

    archiveAndP->ind_array[i+Sarchive->size]=ind_copy(Pareto->ind_array[i]);
  }

  archiveAndP->size=t;
  Sarchive->size=0;

  /* initialize the future size for archive
		     anyway the solutions are lost (-> in archiveAndP) */

  for (i=0;i<t;i++){
    for (j=0;j<t;j++){
      if (i!=j) {
	dom=non_dominated(archiveAndP->ind_array[i],archiveAndP->ind_array[j]);
	if (dom==-1 || (dom==0 && i>j)) j=t+1;
      }
    }
    if(j==t) {
      Sarchive->ind_array[Sarchive->size++]=ind_copy(archiveAndP->ind_array[i]); //save the current
      if (i>=t-Pareto->size) convergence_rate++;
    }
  }
/*printf("%d",Sarchive->size);*/
  complete_free_pop(archiveAndP);
/*printf("encore here");*/

  return convergence_rate;
}


void Initialize_pheromone()
{
    int i;
   for (i=0;i<ni;i++)
	{

		pheromone[i]=tinit;

	}
}





void GwACO(pop *solutions)
{

int  ant, k, i,j,obj, nv, cpt,efficace,cdomine,cdom,cegal,faisable;

int cycle, pos;

int POS,N,nbneg,nbpos, f,t,existe;



/*meilprofit= (double*) chk_malloc(nf * sizeof(double));*/
/*pheromone= (double*) chk_malloc(ni * sizeof(double));*/

int r;
double sum;

cardP=0;
solutions->size=nbants;
			for(ant=0; ant<solutions->size;ant++)
			{

/*printf(" ant %d", ant);*/

/*random_normalisated_weights(); *///random weight
            solutions->ind_array[ant]=create_ind(nf);
			int pris[NBITEMS];
			double capacit[dimension];


			int voisinage [NBITEMS];
			double proba[NBITEMS];



			  for (i=0;i<ni;i++)
			  { pris[i]=0;
			  }

/*for(i=0;i<nf;i++){
    printf("poids %lf",vector_weight[i]);}*/
			  //choisir noeud initial
		pos=(rand()%ni);//items de 0 ŕ ni-1
 /*solutions->ind_array[ant]=create_ind(nf); /*printf("create ind");*/
	solutions->ind_array[ant]->d[0]=pos; /*printf("pos");*/
   solutions->ind_array[ant]->nombr=1; /*printf("nombr");*/
		pris[pos]=1;
		nombr=1;
		cpt=1;
		for(k=0;k<nf;k++)
		{solutions->ind_array[ant]->f[k]=profits[k][pos]; /*printf("profits");*/

		//diminuer poids de pos dans contraintes
	capacit[k]=capacities[k]-weights[k][pos]; /*printf("weights");*/ //ds(f) quantité restante de la ressource f
		}
		do
		{	nv=0;
			//calculer voisinage faisable


			for (obj=0; obj<ni;obj++)
			{
				if (pris[obj]!= 1)
				{
					k=0; faisable=1;
					do{
						if (weights[k][obj]> capacit[k])
							faisable=0;
						k++;
					}while((k<nf)&&(faisable==1));

					if (faisable==1)
					{
					voisinage[nv]=obj;
					nv++;
					}
				}

			}
			//choisir noeud avec probabilité
		if (nv!=0)
	{
	if (nv==1)
		pos=voisinage[0];
	else
	{
		double tot=0,h=0;
		double som[NBITEMS];
        double tmp; double mul;


		for (i=0; i<nv; i++)
			{
				h=0,tmp=1;

			for(j=0;j<nf;j++)//toutes les ressources //hs(oj) ratio: la dureté de l'objet oj par rapport à toutes les contraintes
			{
			h=h+weights[j][voisinage[i]]/capacit[j]; /*printf("voisinage");*/
			}

			for(j=0;j<nf;j++)//toutes les fonctions // ratio profit ressource correspondant à l'inf heuristique
			{
			    eta[j]=0;
			//eta=eta+profits[j][voisinage[i]]/h;/* printf("voisinage 2");*/ //agreration de tt les objectifs
			eta[j]=profits[j][voisinage[i]]/h; //info heuristique pour chaque obj

            mul=beta*vector_weight[j];

            tmp=tmp*pow(eta[j],mul);

		}

	//	eta=profits[f][voisinage[i]]/poids[j][voisinage[i]];//eta dépend que de l'objectif de col //

 	 proba[i]=tmp*monPow(pheromone[voisinage[i]],alphat);
		som[i]= tot + proba[i]; /*printf("som");*/
		tot=som[i]; /*printf("tot %lf",tot);*/
		}/*printf("fin for");*/

/*printf("proba %lf",som);*/
	i=choose(som,nv); /*printf("choose");*/
	pos=voisinage[i]; /*printf("pos 2");*/ //printf("pos %d \n",pos); // pos noeud choisi avec proba
//	printf("cycle %d ant %d card voisinage %d pos cour %d \n",cycle,ant, nv,pos);
	}

			solutions->ind_array[ant]->d[cpt]=pos; /*printf("cpt");*/// ajouter pos à la solution
		for(k=0;k<nf;k++)
			{solutions->ind_array[ant]->f[k]=solutions->ind_array[ant]->f[k]+profits[k][pos]; /*printf("ant");*/ //ajouter profit pos à la sol

			//diminuer poids de pos dans contraintes
			capacit[k]=capacit[k]-weights[k][pos];
			}
/*printf("weights et pos");*/
			nombr++;
			solutions->ind_array[ant]->nombr++; /*printf("nombr solutions");*/ // agmenter nb obj dans sol

			pris[pos]=1;


			cpt++;

		}//if



		}while (nv!=0);////////



} // end for ant

}



void maj_pheromone (pop *Sarchive)
{
    int i,j;

	for (i=0;i<ni;i++)
	{  /* if (pheromone[i]>tmax)
	      /*printf("pheromone%f", pheromone[i]);*/
				/*	pheromone[i]=tmax;

		if (pheromone[i]<tmin)
		pheromone[i]=tmin;*/

		pheromone[i]=rhot*pheromone[i];
    }




	 for(j=0; j<Sarchive->size; j++)
	{


		for(i=0; i<Sarchive->ind_array[j]->nombr; i++){
            pheromone[Sarchive->ind_array[j]->d[i]] =pheromone[Sarchive->ind_array[j]->d[i]]+Sarchive->size;//ajouter qté fixe ŕ chaq sol non dominé
            /*printf(" objet %d",archive->ind_array[j]->d[i]);*/
		}
  /*printf(" pheromone %f"extractPtoArchive(archive,P);*/

	}

}


int main(int argc, char* argv[])
{
	double  duration, durmoy=0;
	clock_t start, finish;
	double	init_time;
    int it=5; int k, i,j;
  
  char *data_path=argv[1];
  loadMOKP(data_path);
  
  /*initfile_weights_log();*/            

  read_weights_file("Weights.txt");

   for(k=1;k<=10;k++)
 {
    init_time=0.0;
    duration=0;
    double mpareto1[15000][dimension];

nombreLIGNE=0; nextLn=0; inv=0;
inter=k;

	FILE *fpareto;
  
  if(strcmp(argv[1],"pymoo\\dataset\\mood_val_dataset\\dataset_0_instance_300_items_3_objectifs.txt")==0){
	fpareto = fopen( "results_val_dataset_0_300_items.txt", "a+" );}

  if(strcmp(argv[1],"pymoo\\dataset\\mood_val_dataset\\dataset_1_instance_300_items_3_objectifs.txt")==0){
	fpareto = fopen( "results_val_dataset_1_300_items.txt", "a+" );}

  if(strcmp(argv[1],"pymoo\\dataset\\mood_val_dataset\\dataset_2_instance_300_items_3_objectifs.txt")==0){
	fpareto = fopen( "results_val_dataset_2_300_items.txt", "a+" );}

  if(strcmp(argv[1],"pymoo\\dataset\\mood_val_dataset\\dataset_3_instance_300_items_3_objectifs.txt")==0){
	fpareto = fopen( "results_val_dataset_3_300_items.txt", "a+" );}

  if(strcmp(argv[1],"pymoo\\dataset\\mood_val_dataset\\dataset_4_instance_300_items_3_objectifs.txt")==0){
	fpareto = fopen( "results_val_dataset_4_300_items.txt", "a+" );}
  

	fprintf(fpareto,"mcycle %d nbants %d alphat %d beta %lf rho %lf tmax %lf\n",maxcycle,nbants,alphat,beta,rhot,tmax);
	fflush(stdout);





seed(inter);

P=create_pop(paretoIni,nf);



  /*  int cardP2=0,j,cardP3=0;
	double pareto2[8000][dimension];//, pareto3[100];

/*for (i=0;i<paretoIni;i++)
{ P->ind_array[i]=create_ind(nf);}*/

/*maxcycle=nombreLIGNE*2;*/

Initialize_pheromone();
int cycle=0;

do{
        printf("cycle %d",cycle);

solutions=create_pop(nbants,nf);
archive=create_pop(paretoIni,nf);


dynamic_weight_allpop2(); //log weight



start=clock();

    GwACO(solutions); /*printf( " MMAS\n" );*/

    extractPtoArchive(solutions,archive);

    maj_pheromone (archive);

    extractPtoArchive(archive,P);

finish = clock(); /*printf( " time" );*/
duration = (double)(finish - start) / CLOCKS_PER_SEC;
init_time+=duration;

complete_free_pop(solutions);
complete_free_pop(archive);

cycle++;


	}while (cycle<maxcycle);

fprintf(fpareto,"cardinalité ensemble Pareto %d \n temps CPU %f \n", P->size, init_time);

	for(i=0;i<P->size;i++)
	{for(j=0;j<nf;j++)
		{
            
        fprintf(fpareto,"%f ",P->ind_array[i]->f[j]);
		fflush(stdout);
		mpareto1[i][j]=P->ind_array[i]->f[j];
		}
		fprintf(fpareto,"\n");

	}




/*fprintf(fpareto,"ACOeps \n");
cardP2= otherResult("reslt ACOEps\\ACOEps750.4\\750.4.30.txt",pareto2);
	fprintf(fpareto,"C Mesure 1 %lf \n", Cmesure(mpareto1, pareto2, cardP, cardP2));
	fprintf(fpareto,"C Mesure 2 %lf \n", Cmesure(pareto2, mpareto1, cardP2, cardP));*/

	/*fprintf(fpareto,"ACOHD \n");
cardP2= otherResult("reslt ACOHD\\ACOHD250.2\\250.2.30.txt",pareto2);

	fprintf(fpareto,"C Mesure 1 %lf \n", Cmesure(mpareto1, pareto2, cardP, cardP2));
	fprintf(fpareto,"C Mesure 2 %lf \n", Cmesure(pareto2, mpareto1, cardP2, cardP));

fprintf(fpareto,"IBEAeps \n");
cardP2= otherResult("C:\\Users\\imen\\Desktop\\projet\\ibea_knapsack\\knapsack\\bin\\Debug\\IBEAeps\\IBEA250.2\\knapsack_output30.txt",pareto2);
	fprintf(fpareto,"C Mesure 1 %lf \n", Cmesure(mpareto1, pareto2, cardP, cardP2));
	fprintf(fpareto,"C Mesure 2 %lf \n", Cmesure(pareto2, mpareto1, cardP2, cardP));

	fprintf(fpareto,"IBEAhd \n");
cardP2= otherResult("C:\\Users\\imen\\Desktop\\projet\\ibea_knapsack\\knapsack\\bin\\Debug\\IBEAHD\\IBEA250.2\\knapsack_output30.txt",pareto2);
	fprintf(fpareto,"C Mesure 1 %lf \n", Cmesure(mpareto1, pareto2, cardP, cardP2));
	fprintf(fpareto,"C Mesure 2 %lf \n", Cmesure(pareto2, mpareto1, cardP2, cardP));

	fprintf(fpareto,"spea2 \n");
cardP2= otherResult("C:\\Users\\imen\\Desktop\\projet\\SPEA2_knapsack\\knapsack\\bin\\Debug\\SPEA2250.2\\knapsack_output30.txt",pareto2);

	fprintf(fpareto,"C Mesure 1 %lf \n", Cmesurerandom_normalisated_weights()(mpareto1, pareto2, cardP, cardP2));
	fprintf(fpareto,"C Mesure 2 %lf \n", Cmesure(pareto2, mpareto1, cardP2, cardP));

	fprintf(fpareto,"HypE \n");
cardP2= otherResult("C:\\Users\\imen\\Desktop\\projet\\HypE_knapsack\\knapsack\\bin\\Debug\\HypE250.2\\knapsack_output30.txt",pareto2);
	fprintf(fpareto,"C Mesure 1 %lf \n", Cmesure(mpareto1, pareto2, cardP, cardP2));
	fprintf(fpareto,"C Mesure 2 %lf \n", Cmesure(pareto2, mpareto1, cardP2, cardP));

	/*fprintf(fpareto,"mACO \n");
cardP2= otherResult("result500.3\\result500.3.1.txt",pareto2);
	fprintf(fpareto,"C Mesure 1 %lf \n", Cmesure(mpareto1, pareto2, cardP, cardP2));
	fprintf(fpareto,"C Mesure 2 %lf \n", Cmesure(pareto2, mpareto1, cardP2, cardP));*/


	complete_free_pop(P);


/*fprintf(fpareto,"Moyenne temps CPU %f \n", durmoy);*/

	fclose(fpareto);
/*printf( " \n before free" );*/

/*free(bounds); printf( " \nfree bounds" );*/





}


return(0);

}




