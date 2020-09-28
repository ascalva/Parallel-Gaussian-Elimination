#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int DEF_N=4;

void check_comp(int N, double A[N][N], double X[N], double result[N])
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
void compute_ge_(int N, double A[N][N], double B[N], double y[N])
{
    int k, j, i;

    for( k = 0; k < N; k++ )
    {
        for( j = k + 1; j < N; j++ )
            A[k][j] = A[k][j] / A[k][k];

        B[k] = B[k] / A[k][k];

        for( i = k + 1; i < N; i++ )
        {
            for( j = k + 1; j < N; j++ )
            {
                A[i][j] -= (A[i][k] * A[k][j]);
            }

            B[i]   -= (A[i][k] * B[k]);
        }
    }
}
void compute_ge(int N, double A[N][N], double B[N], double y[N])
{
    int k, j, i;

    for( k = 0; k < N; k++ )
    {
        for( j = k + 1; j < N; j++ )
            A[k][j] = A[k][j] / A[k][k];

        B[k]    = B[k] / A[k][k];
        //A[k][k] = 1;

        for( i = k + 1; i < N; i++ )
        {
            for( j = k + 1; j < N; j++ )
                A[i][j] -= (A[i][k] * A[k][j]);

            B[i]   -= (A[i][k] * B[k]);
            //A[i][k] = 0;
        }
    }
}

void back_sub(int N, double U[N][N], double Y[N], double X[N])
{
    int k, i;

    for( k = N - 1; k >= 0; k-- )
    {
        X[k] = Y[k];
        
        for( i = k - 1; i >= 0; i-- )
            Y[i] -= (X[k] * U[i][k]);
    }
}


void initialize_inputs(int N, double A[N][N], double B[N], double X[N]) 
{
    int row, col;

    printf("\nInitializing...\n");
    for (col = 0; col < N; col++) {
        for (row = 0; row < N; row++)
            A[row][col] = (double)rand() / 32768.0;

        B[col] = (double)rand() / 32768.0;
        X[col] = 0.0;
    }
}

void print_inputs(int N, double A[N][N], double B[N], double X[N])
{
    int row, col;

    printf("\nA =\t");
    for (row = 0; row < N; row++) {
        for (col = 0; col < N; col++) {
            printf("%6.2f%s", A[row][col], (col < N-1) ? " " : "\n\t");
        }
    }
    printf("\nB =\t");
    for (col = 0; col < N; col++) {
        printf("%6.2f%s", B[col], (col < N-1) ? " " : "\n");
    }
}

int main(int argc, char *argv[])
{
    int N = DEF_N;
    
    // Use CL arg if provided for size of matrices (<= 40).
    if( argc == 2 )
    {
        N = atoi( argv[1] );
    }

    // Initialize matrices.
    //double A[N][N], B[N], X[N];
    //double a[N][N], b[N];

    double (*A)[N] = malloc(sizeof(double[N][N]));
    double  *B     = malloc(sizeof(double[N]));
    double  *X     = malloc(sizeof(double[N]));

    double (*a)[N] = malloc(sizeof(double[N][N]));
    double  *b     = malloc(sizeof(double[N]));

    initialize_inputs(N,A,B,X);
   
    // Make deep copy of arrays.
    int row, col;
    for( row = 0; row < N; row++ )
    {
        b[row] = B[row];

        for( col = 0; col < N; col++ )
            a[row][col] = A[row][col];
    }

    // Print original matrices.
    //print_inputs(N,A,B,X);

    // Compute Gaussian Elimination.
    compute_ge(N,A,B,X);

    // Print Upper Triangle Matrix.
    //print_inputs(N,A,B,X);

    // Perform back-substitution.
    back_sub(N,A,B,X);

    check_comp(N, a,X,b);
    
    // Free arrays
    free(A);
    free(B);
    free(X);
    free(a);
    free(b);
}
