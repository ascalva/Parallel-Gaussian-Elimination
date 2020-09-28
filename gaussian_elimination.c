#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <pthread.h>

#define MAX_THREAD 20

#define NDIM 20

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

double          a[NDIM][NDIM];
//double          b[NDIM];
//double          x[NDIM];

typedef struct
{
	int         totalProc;
	int         dim;
    int         l_indx;
    int         u_indx;
    int         k;
	double	    (*a)[][NDIM], /* *(*a)[], */
                *b; 
} threadData;


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

void initialize_inputs(
        int    N, 
        double (*A)[N], 
        double  *B, 
        double  *X, 
        double (*a_og)[N],
        double  *b_og )
{
    int row, col; 

    printf("\nInitializing...\n");
    for (col = 0; col < NDIM; col++) {
        for (row = 0; row < NDIM; row++) {
            A[row][col] = a[row][col] = a_og[row][col] = (double)rand() / 32768.0;
        }
        B[col] = b_og[col] = (double)rand() / 32768.0;
        X[col] = 0.0;
    }
}

void check_comp(int N, double (*A)[N], double *X, double *result)
{
    // Compute dot products.
    int row, col;
    double res;
    
    printf("\nActual    Comp\n");
    for( row = 0; row < N; row++ )
    {
        // Reset result of dot product.
        res = 0;

        // Compute dot product for specifc row against column vector found through
        // Gaussian Elimination.
        for( col = 0; col < N; col++ )
        {
            res += A[row][col] * X[col];
        }
        printf(" %5.2f   %5.2f   %s\n", result[row], res, 
                ((float)result[row] == (float)res ? "success" : "fail"));
    }
}


void back_sub(int N, double (*U)[N], double *Y, double *X)
{
    int k, i;

    for( k = N - 1; k >= 0; k-- )
    {
        X[k] = Y[k];
        
        for( i = k - 1; i >= 0; i-- )
            Y[i] -= (X[k] * U[i][k]);
    }
}


void comp_ge(int N, int row, int row_lim, int k, double (*A)[N], double *B)
{
    int j, r;
    
    // Only iterate over specified rows (no two threads will access the 
    // rows).
    for( r = row; r < row_lim; r++ )
    {
        for( j = k + 1; j < N; j++ )
        {
            A[r][j] -= (A[r][k] * A[k][j]);
        }

        B[r] -= (A[r][k] * B[k]);
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

//Task for each thread can be initialized here
void * worker(void *arg)
{
	threadData           *p = (threadData *) arg;
    comp_ge(p->dim, p->l_indx, p->u_indx, p->k, *(p->a), (p->b));

    return 0;
}

int main(int argc, char *argv[])
{
	int n;

	//Create array of threads "Dynamic"
	pthread_t      *threads;

	//Create struct object for each thread
	threadData           *arg;

	// Generate matrix A and column vector B such that Ax = B, for some column
    // vector x.
    double (*A)[NDIM]    = malloc(sizeof(double[NDIM][NDIM]));
    double  *B           = malloc(sizeof(double[NDIM]));
    double  *X           = malloc(sizeof(double[NDIM]));

    double (*a_og)[NDIM] = malloc(sizeof(double[NDIM][NDIM]));
    double  *b_og        = malloc(sizeof(double[NDIM]));

	initialize_inputs(NDIM, A, B, X, a_og, b_og);

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

	threads = (pthread_t  *) malloc(n * sizeof(pthread_t ));
	arg     = (threadData *) malloc(n * sizeof(threadData));
	
	/* Start up thread */

    // Define iterators.
    int i, j, k, p;

    // Define intermediate values.
    int n_elements, start_pos, end_pos;

	/* Spawn thread */
    for( k = 0; k < NDIM; k++ )
    {
        for( j = k + 1; j < NDIM; j++ )
        {
            a[k][j] = a[k][j] / a[k][k];
        }

        B[k] = B[k] / a[k][k];

        // Compute number of elements (rows) per processor.
        n_elements = ceil((double)(NDIM - k - 1) / (double)n);
        
        /* Spawn threads for N >= P */
        if( NDIM - k - 1 >= n )
        {
            for( p = 0; p < n; p++ )
            {
                // Compute row limits for threads.
                start_pos        = n_elements * p + k + 1;
                end_pos          = MIN(start_pos + n_elements, NDIM);

                // Setup thread.
                arg[p].totalProc = n; 
                arg[p].dim       = NDIM;
                arg[p].l_indx    = start_pos;
                arg[p].u_indx    = end_pos;
                arg[p].k         = k;
                arg[p].a         = &a;
                arg[p].b         = B;

                pthread_create(&threads[p], NULL, worker, (void *)(arg+p));
            }

            // Join threads.
            for (i = 0; i < n; i++)
            {
            	pthread_join(threads[i], NULL);
            }
        }

        /* Spawn threads for N < P */
        else
        {
            for( p = 0; p < NDIM - k - 1; p++ )
            {
                // Each thread will only compute one row.
                start_pos        = n_elements * p + k + 1;

                // Setup thread.
                arg[p].totalProc = n; 
                arg[p].dim       = NDIM;
                arg[p].l_indx    = start_pos;
                arg[p].u_indx    = start_pos + 1;
                arg[p].k         = k;
                arg[p].a         = &a;
                arg[p].b         = B;
                
                pthread_create(&threads[p], NULL, worker, (void *)(arg+p));
            }

            // Join threads.
            for (i = 0; i < p; i++)
            {
            	pthread_join(threads[i], NULL);
            }
        }
    }

    // Compute back-substitution.
    back_sub(NDIM,a,B,X);
    
    // Check if correct
    check_comp(NDIM, a_og, X, b_og);
    
	free(arg);
}

