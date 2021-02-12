#include <iostream>
#include <assert.h>
#include <iomanip>
#include <mpi.h>
#include <cmath>
#include <sstream>
#include <limits.h>
#include "cstdlib"


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
            std::cout << c[i * n + j] << "\t";
        }
        std::cout << std::endl;
    }
}

void ElemMatrix(int *c, size_t n, bool isRand=0)
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

int main(int argc, char *argv[])
{
    char buf[PATH_MAX];
    setbuf(stdout, buf);

    MPI_Init(&argc, &argv);

        for (int nRows = 64; nRows <= 4096; nRows *= 2) {

            int nProcs, rank;
            MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            int *a = new int[nRows * nRows];
            int *b = new int[nRows * nRows];
            int *c = new int[nRows * nRows];
            int *cMpiRes = new int[nRows * nRows];
            int alpha = 3;
            int betta = 4;

            if (!rank) {
                ElemMatrix(a, nRows, 1);
                ElemMatrix(b, nRows, 1);
                ElemMatrix(c, nRows, 1);
                ElemMatrix(cMpiRes, nRows);
            }
            MPI_Bcast(a, nRows * nRows, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(b, nRows * nRows, MPI_INT, 0, MPI_COMM_WORLD);

            double timerMpi;

            size_t startIdx = nRows / nProcs * rank + std::min<size_t>(rank, nRows % nProcs);
            size_t step = nRows / nProcs + (rank < nRows % nProcs);
            int *cMpi = new int[nRows * nRows];
            ElemMatrix(cMpi, nRows);
            MPI_Barrier(MPI_COMM_WORLD);
            timerMpi = MPI_Wtime();
            MultMatrix(a, b, c, cMpi, nRows, startIdx, step, alpha, betta);
            MPI_Barrier(MPI_COMM_WORLD);
            timerMpi = MPI_Wtime() - timerMpi;
            MPI_Barrier(MPI_COMM_WORLD);

            MPI_Reduce(cMpi, cMpiRes, nRows * nRows, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

            if (!rank) {
                std::cout << "nRows: " << nRows << " nProcs: " << nProcs << " timer: " << timerMpi << std::endl;
            }

            delete[] a;
            delete[] b;
            delete[] c;
            delete[] cMpi;
            delete[] cMpiRes;
        }
    MPI_Finalize();
    return 0;
}
