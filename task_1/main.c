#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int wrank, wsize;
    unsigned L = 8000; 
    unsigned K = 100; 
    unsigned Count = L / (2 * K); 
    double* message = malloc(L * sizeof(double));
    int paths[2][7] = {{0, 1, 2, 3, 7, 11, 15}, {0, 4, 8, 12, 13, 14, 15}};
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    if (wrank == 0) {
        for (unsigned i = 0; i < L; i++) {
            message[i] = i;
        }
    }
    if (wrank == 0) {
        for (unsigned i = 0; i < L; i++) {
            printf("First process%lf\n",  message[i]);
        }
    }
    
    for (unsigned i = 0, tag = 1215; i < 2; i++, tag++) {
        for (unsigned j = 1; j < 7; j++) {
            if (paths[i][j] == wrank) {
                unsigned c = 0;
                if (wrank == 15 && i == 1) {
                    c = K * Count;
                }
                for (unsigned k = 0; k < K; k++) {
                    MPI_Recv( message + k * Count + c, Count, MPI_DOUBLE, paths[i][j - 1], tag, MPI_COMM_WORLD, &status);
                }
                break;
            }
        }
    }

    unsigned buff_size = L * sizeof(double) + MPI_BSEND_OVERHEAD;
    double* buff = malloc(buff_size);
    MPI_Buffer_attach(buff, buff_size);
    for (unsigned i = 0, tag = 1215; i < 2; i++, tag++) {
        for (unsigned j = 0; j < 6; j++) {
            if (paths[i][j] == wrank) {
                unsigned c = 0;
                if (wrank == 0 && i == 1) {
                    c = K * Count;
                }
                for (unsigned k = 0; k < K; k++) {
                    MPI_Bsend( message + k * Count + c, Count, MPI_DOUBLE, paths[i][j + 1], tag, MPI_COMM_WORLD);
                }
                break;
            }
        }
    }
    MPI_Buffer_detach(buff, &buff_size);
    free(buff);
    
    if (wrank == 15) {
        for (unsigned i = 0; i < L; i++) {
            printf("%lf\n",  message[i]);
        }
    }
    free( message);

    MPI_Finalize();
    return 0;
}