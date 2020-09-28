#include <stdio.h>
#include <sys/types.h>
#include <pthread.h>
#include <stdlib.h>

#define MAX_THREAD 20

#define NDIM 16

double          a[NDIM][NDIM];
double          b[NDIM][NDIM];
double          c[NDIM][NDIM];

typedef struct
{
	int         rank;
	int         totalProc;
	int         dim;
	double	    (*a)[NDIM][NDIM], 
                (*b)[NDIM][NDIM], 
                (*c)[NDIM][NDIM];
} threadData;



void print_matrix(int dim)
{
	int i,j;

	printf("The %d * %d matrix is\n", dim,dim);
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++)
			printf("%lf ",  c[i][j]);
		printf("\n");
	}
}

void check_matrix(int dim)
{
	int i,j,k;
	int error=0;

	printf("Now checking the results\n");
	for(i=0;i<dim;i++)
		for(j=0;j<dim;j++) {
			double e=0.0;

			for (k=0;k<dim;k++)
				e+=a[i][k]*b[k][j];

			if (e!=c[i][j]) {
				printf("(%d,%d) error\n",i,j);
				error++;
			}
		}
	if (error)
		printf("%d elements error\n",error);
		else
		printf("success\n");
}

void matrixMultiplication(int me_no, int noproc, int n, double a[NDIM][NDIM], double b[NDIM][NDIM], double c[NDIM][NDIM])
{
	int             i,j,k;
	double sum;
	i=me_no;
	while (i<n) {

		for (j = 0; j < n; j++) {
			sum = 0.0;
			for (k = 0; k < n; k++) {
				sum = sum + a[i][k] * b[k][j];
			}
			c[i][j] = sum;

		}
		i+=noproc;
	}
}


//Task for each thread can be initialized here
void * worker(void *arg)
{
	threadData           *p = (threadData *) arg;
	matrixMultiplication(p->rank, p->totalProc, p->dim, *(p->a), *(p->b), *(p->c));
}

int main(int argc, char *argv[])
{
	int             j, k, noproc, me_no;
	double          sum;
	double          t1, t2;

	//Create array of threads "Dynamic"
	pthread_t      *threads;

	//Create struct object for each thread
	threadData           *arg;

	int             n, i;

	//Assign matrix for multiplication
	for (i = 0; i < NDIM; i++)
		for (j = 0; j < NDIM; j++)
		{
			a[i][j] = i + j;
			b[i][j] = i + j;
		}

	if (argc != 2)
	{
		printf("Usage: %s n\n  where n is no. of thread\n", argv[0]);
		exit(1);
	}
	n = atoi(argv[1]);

	if ((n < 1) || (n > MAX_THREAD))
	{
		printf("The no of thread should be between 1 and %d.\n", MAX_THREAD);
		exit(1);
	}
	threads = (pthread_t *) malloc(n * sizeof(pthread_t));

	arg=(threadData *)malloc(sizeof(threadData)*n);
	

	/* Start up thread */

	/* Spawn thread */
	for (i = 0; i < n; i++)
	{
		arg[i].rank      = i;
		arg[i].totalProc = n; 
		arg[i].dim       = NDIM;
		arg[i].a         = &a;
		arg[i].b         = &b;
		arg[i].c         = &c;
		pthread_create(&threads[i], NULL, worker, (void *)(arg+i));
	}

	for (i = 0; i < n; i++)
	{
		pthread_join(threads[i], NULL);

	}
	print_matrix(NDIM);
	check_matrix(NDIM);
	free(arg);
}

