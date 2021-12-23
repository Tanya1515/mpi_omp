#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <mpi-ext.h>
#include <signal.h>

#define KILL_PROC 2
int KILL = 0;
MPI_Comm mpi_comm;
size_t startIdx;
size_t step;
int nRows = 8;
int alpha = 3;
int betta = 4;

int* a = NULL;
int* b = NULL;
int* c = NULL;
int* cMpi = NULL;


void printMatrix(int *c, int n)
{
    printf("Begining of printing\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d  ", c[i * n + j]);
        }
        printf("\n");
    }
    printf("End of printing\n");
}


void MultMatrix(int *a, int *b, int *c, int *cMpi, size_t nRows, size_t startIdx, size_t step, int alpha, int betta, int rank)
{
    size_t i, j, k;
    int sum;
    if ((rank == KILL_PROC) && (KILL == 0))
    {
        raise(SIGKILL);
        KILL = 1;
    }
    for (i = 0; i < nRows; i++) {
        for (j = startIdx; j < nRows && j < startIdx + step; j++) {
            sum = 0;
            for (k = 0; k < nRows; k++) {
                sum += a[i * nRows + k] * b[k * nRows + j];
            }
            cMpi[i * nRows + j] = sum * alpha + c[i * nRows + j] * betta;
        }
    }
}


int min (int num1, int num2)
{
    if (num1 > num2)
        return num2;
    else
        return num1;
}


void ElemMatrix(int *c, size_t n, int isRand)
{
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (isRand != 0)
                c[i * n + j] = (i+j+isRand)*(n+i+j);
            else
                c[i * n + j] = 0;
            
        }
    }
}


static void verbose_errhandler(MPI_Comm* comm, int* err, ...) {
    int rank, size, amount_f, len;
    int old_size;
    int old_rank;
    char errstr[MPI_MAX_ERROR_STRING];
    
    MPI_Group group_f;
    MPI_Group group_norm;
    MPI_Comm_rank(mpi_comm, &old_rank);
    MPI_Comm_size(mpi_comm, &old_size);
    int* norm_ranks = malloc(sizeof(int)*old_size);
    int* f_ranks = malloc(sizeof(int)*amount_f);
    //printf ("Amount of processes in communicator with failed processes: %d\n", old_size);
    if (old_rank == 0)
    {
        MPI_Comm_group(mpi_comm, &group_norm);
        MPIX_Comm_failure_ack(mpi_comm);
        MPIX_Comm_failure_get_acked(mpi_comm, &group_f);
        MPI_Group_size(group_f, &amount_f);
        for (int i = 0; i<amount_f; i++)
           f_ranks[i] = i;
        MPI_Group_translate_ranks(group_f, amount_f, f_ranks, group_norm, norm_ranks);
    }
    MPI_Error_string( *err, errstr, &len );
    //printf("Rank %d / %d: Notified of error %s in %d processes\n", rank, size, errstr, amount_f);
    
    MPIX_Comm_shrink(*comm, &mpi_comm);
    MPI_Comm_rank(mpi_comm, &rank);
    MPI_Comm_size(mpi_comm, &size);
    //printf ("Amount of processes in communicator without failed processes: %d\n", size);
    MPI_Barrier(mpi_comm);
    
    if (old_rank == 0)
    {
        startIdx = (nRows / old_size * norm_ranks[0]) + min (norm_ranks[0], (nRows % old_size));
        step = nRows / old_size + (rank < (nRows % old_size));
        MultMatrix(a, b, c, cMpi, nRows, startIdx, step, alpha, betta, rank);
    }
}

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);

            int nProcs, rank;
            MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            mpi_comm = MPI_COMM_WORLD;

            a = malloc(sizeof(int)*nRows*nRows);
            b = malloc(sizeof(int)*nRows*nRows);
            c = malloc(sizeof(int)*nRows*nRows);
            int *cMpiRes = malloc(sizeof(int)*nRows*nRows);
            cMpi = malloc(sizeof(int)*nRows*nRows);

            MPI_Errhandler errh;

            MPI_Comm_create_errhandler(verbose_errhandler, &errh);
            MPI_Comm_set_errhandler(MPI_COMM_WORLD, errh);
            
            if (!rank) {
                ElemMatrix(a, nRows, 1);
                ElemMatrix(b, nRows, 1);
                ElemMatrix(c, nRows, 1);
                ElemMatrix(cMpiRes, nRows, 0);
                ElemMatrix(cMpi, nRows, 0);
            }
            
            MPI_Bcast(a, nRows * nRows, MPI_INT, 0, mpi_comm);
            MPI_Bcast(b, nRows * nRows, MPI_INT, 0, mpi_comm);

            double timerMpi;

            startIdx = (nRows / nProcs * rank) + min (rank, (nRows % nProcs));
            step = nRows / nProcs + (rank < (nRows % nProcs));
            MPI_Barrier(mpi_comm);
            timerMpi = MPI_Wtime();
            MultMatrix(a, b, c, cMpi, nRows, startIdx, step, alpha, betta, rank);
            MPI_Barrier(mpi_comm);
            timerMpi = MPI_Wtime() - timerMpi;
            MPI_Barrier(mpi_comm);

            MPI_Reduce(cMpi, cMpiRes, nRows * nRows, MPI_INT, MPI_SUM, 0, mpi_comm);

            MPI_Comm_size(mpi_comm, &nProcs);
            if (!rank) {
                printf("Size of matrix: %d * %d \n", nRows, nRows);
                printf("Amount of processes, that are alive: %d \n", nProcs);
                printf("Time for matrix computing: %f \n", timerMpi);
            }
            
            if (rank == 0)
            {
                printf("Final matrix: \n");
                printMatrix(cMpiRes, nRows);
            }

            free(a);
            free(b);
            free(c);
            free(cMpi);
            free(cMpiRes);
    MPI_Finalize();
    return 0;
}
