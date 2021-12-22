#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <mpi-ext.h>
#include <signal.h>

void MultMatrix(int *a, int *b, int *c, int *cMpi, size_t nRows, size_t startIdx, size_t step, int alpha, int betta)
{
    size_t i, j, k;
    int sum;
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

void printMatrix(int *c, int n)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d", c[i * n + j]);
        }
    }
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

/*static void verbose_errhandler(MPI_Comm* comm, int* err, ...) {
    MPI_Comm mpi_comm = *comm;
    int rank, size, amount_f, len;
    char errstr[MPI_MAX_ERROR_STRING];
    MPI_Group group_f;

    MPI_Comm_rank(mpi_comm, &rank);
    MPI_Comm_size(mpi_comm, &size);
    printf("Size of communicator: %d\n", size);

    MPIX_Comm_failure_ack(mpi_comm);
    MPIX_Comm_failure_get_acked(mpi_comm, &group_f);
    MPI_Group_size(group_f, &amount_f);

    MPI_Error_string( *err, errstr, &len );
    printf("Rank %d / %d: Notified of error %s in %d processes\n",
           rank, size, errstr, amount_f);
    MPIX_Comm_shrink(mpi_comm, &mpi_comm);

    MPI_Barrier(mpi_comm);

    MPI_Comm_size(mpi_comm, &size);
    printf("Size of communicator without failed process: %d\n", size);

}*/

int min (int num1, int num2)
{
    if (num1 > num2)
        return num2;
    else 
        return num1;
}

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);

        for (int nRows = 64; nRows <= 128; nRows *= 2) {

            int nProcs, rank;
            MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            int *a = malloc(sizeof(int)*nRows*nRows);
            int *b = malloc(sizeof(int)*nRows*nRows);
            int *c = malloc(sizeof(int)*nRows*nRows);
            int *cMpiRes = malloc(sizeof(int)*nRows*nRows);
            int alpha = 3;
            int betta = 4;

            //MPI_Errhandler errh;

            //MPI_Comm_create_errhandler(verbose_errhandler, &errh);
            //MPI_Comm_set_errhandler(MPI_COMM_WORLD, errh);

            //if (rank == 2) {
            //    raise(SIGKILL);
            //}

            if (!rank) {
                ElemMatrix(a, nRows, 1);
                ElemMatrix(b, nRows, 1);
                ElemMatrix(c, nRows, 1);
                ElemMatrix(cMpiRes, nRows, 0);
            }
            MPI_Bcast(a, nRows * nRows, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(b, nRows * nRows, MPI_INT, 0, MPI_COMM_WORLD);

            double timerMpi;

            size_t startIdx = (nRows / nProcs * rank) + min (rank, (nRows % nProcs));
            size_t step = nRows / nProcs + (rank < (nRows % nProcs));
            int *cMpi = malloc(sizeof(int)*nRows*nRows);
            ElemMatrix(cMpi, nRows, 0);
            MPI_Barrier(MPI_COMM_WORLD);
            timerMpi = MPI_Wtime();
            MultMatrix(a, b, c, cMpi, nRows, startIdx, step, alpha, betta);
            MPI_Barrier(MPI_COMM_WORLD);
            timerMpi = MPI_Wtime() - timerMpi;
            MPI_Barrier(MPI_COMM_WORLD);

            MPI_Reduce(cMpi, cMpiRes, nRows * nRows, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

            if (!rank) {
                printf("nRows: %d \n", nRows);
                printf("nProcs: %d \n", nProcs);
                printf("timer: %f \n", timerMpi);
            }

            free(a);
            free(b);
            free(c);
            free(cMpi);
            free(cMpiRes);
        }
    MPI_Finalize();
    return 0;
}



