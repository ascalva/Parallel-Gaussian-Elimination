/**
 * filename: gaussian_elimination.c
 *
 * @author : Alberto Serrano-Calva (axs4986)
 *
 * purpose : Performs Gaussian Elimination algorithm in parallel using threads.
 *           Solves the equation Ax = y, such that A is an (NxN) matrix, x and y are
 *           both (Nx1) column vectors, where x is unknown.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <pthread.h>
#include <time.h>

#define MAX_THREAD 20

// Define data structure for passing arguments to threads.
typedef struct
{
	int         totalProc;
	int         dim;
    int         l_indx;
    int         u_indx;
    int         k;
	double	    (*a)[],
                 *b; 
} threadData;


// Initialize inputs such that a solution is guaranteed for the equation Ax = B, where
// x is the column vector we're solving for.
void initialize_inputs(
        int    N, 
        double (*A)[N], 
        double  *B, 
        double  *X, 
        double (*a_og)[N],
        double  *b_og )
{
    int row, col; 

    //printf("\nInitializing...\n");
    for (col = 0; col < N; col++) {
        for (row = 0; row < N; row++) {
            A[row][col] = a_og[row][col] = (double)rand() / 32768.0;
        }
        B[col] = b_og[col] = (double)rand() / 32768.0;
        X[col] = 0.0;
    }
}


// Computes B (known) using computed column vector X, if elements match, prints
// 'success', otherwise, 'fail'.
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


// Back substitution function, completely serial-based programming.
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


// Section of Gaussian Elimination that can be computed in parallel, or without
// any locks.
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


//Task for each thread can be initialized here
void * worker(void *arg)
{
	threadData           *p = (threadData *) arg;
    comp_ge(p->dim, p->l_indx, p->u_indx, p->k, (p->a), (p->b));

    return 0;
}


int main(int argc, char *argv[])
{
	int n, ndim;

    // Clock variables.
    clock_t c_begin_init, c_end_init, c_begin_exec, c_end_exec;

	//Create array of threads "Dynamic"
	pthread_t      *threads;

	//Create struct object for each thread
	threadData           *arg;

	if (argc != 3)
	{
		printf("Usage: %s d n\n  where d is number of dimensions and  n is no. of thread\n",
                argv[0]);
		exit(1);
	}

    // Parse CL arguments.
	n    = atoi(argv[2]);
    ndim = atoi(argv[1]);

    // Error checking.
	if ((n < 1) || (n > MAX_THREAD))
	{
		printf("The no of thread should be between 1 and %d.\n", MAX_THREAD);
		exit(1);
	} 
    
    if( ndim < 1 )
    {
        printf("The no of dimensions should be a positive integer.\n");
        exit(1);
    }

    // Print program properties.
    printf("Number of threads:                 %d\n", n);
    printf("Number of dimensions:              %d\n", ndim);

    // Initialize matrices and time.
    c_begin_init         = clock();

	// Generate matrix A and column vector B such that Ax = B, for some column
    // vector x.
    double (*A)[ndim]    = malloc(sizeof(double[ndim][ndim]));
    double  *B           = malloc(sizeof(double[ndim]));
    double  *X           = malloc(sizeof(double[ndim]));

    double (*a_og)[ndim] = malloc(sizeof(double[ndim][ndim]));
    double  *b_og        = malloc(sizeof(double[ndim]));

	initialize_inputs(ndim, A, B, X, a_og, b_og);

    c_end_init           = clock();

    // Allocate thread argument storage.
	threads = (pthread_t  *) malloc(n * sizeof(pthread_t ));
	arg     = (threadData *) malloc(n * sizeof(threadData));
	
	/* Start up thread */

    // Define iterators.
    int i, j, k, p;

    // Define intermediate values.
    int n_elements, start_pos, end_pos;

	/* Spawn thread */

    // Start timer.
    c_begin_exec = clock();
    for( k = 0; k < ndim; k++ )
    {
        for( j = k + 1; j < ndim; j++ )
        {
            A[k][j] = A[k][j] / A[k][k];
        }

        B[k] = B[k] / A[k][k];

        // Compute number of elements (rows) per processor.
        n_elements = ceil((double)(ndim - k - 1) / (double)n);
        
        /* Spawn threads for N >= P */
        if( ndim - k - 1 >= n )
        {
            for( p = 0; p < n; p++ )
            {
                // Compute row limits for threads.
                start_pos        = n_elements * p + k + 1;
                end_pos          = fmin(start_pos + n_elements, ndim);

                // Setup thread.
                arg[p].totalProc = n; 
                arg[p].dim       = ndim;
                arg[p].l_indx    = start_pos;
                arg[p].u_indx    = end_pos;
                arg[p].k         = k;
                arg[p].a         = A;
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
            for( p = 0; p < ndim - k - 1; p++ )
            {
                // Each thread will only compute one row.
                start_pos        = n_elements * p + k + 1;

                // Setup thread.
                arg[p].totalProc = n; 
                arg[p].dim       = ndim;
                arg[p].l_indx    = start_pos;
                arg[p].u_indx    = start_pos + 1;
                arg[p].k         = k;
                arg[p].a         = A;
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
    back_sub(ndim, A, B, X);
    
    // End timer.
    c_end_exec = clock();

    /* Print Analysis */
    // Print times.
    // NOTE: Additional division by factor of 10, result kept being about 10 times 
    //       over of what it should be.
    printf("Time taken to initialize matrices: %3.4f (s)\n", 
            ((double)(c_end_init - c_begin_init)) / (CLOCKS_PER_SEC * 10));

    printf("Time taken to compute vector:      %3.4f (s)\n",
            ((double)(c_end_exec - c_begin_exec)) / (CLOCKS_PER_SEC * 10));

    // Check if correct
    check_comp(ndim, a_og, X, b_og);
    
    // Cleanup
    free(threads);
	free(arg);
    free(A);
    free(B);
    free(X);
    free(a_og);
    free(b_og);
}

