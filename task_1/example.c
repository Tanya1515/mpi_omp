#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define M 10


int main( int argc, char **argv ) {
    int n;
    int rank, size;
    MPI_Status status;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &size ); 
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if ( rank == 0 ) {
        int blen = M * (sizeof(int) + MPI_BSEND_OVERHEAD);
        int *buf = (int*) malloc(blen);
        MPI_Buffer_attach (buf, blen);
        for(int i = 0; i < M; i ++) {
            n = i;
            MPI_Bsend (&n, 1, MPI_INT, 1, i, MPI_COMM_WORLD );
        } 
        MPI_Buffer_detach(&buf, &blen);
        free(buf);
    }
    else if ( rank == 1) {
        for(int i = 0; i < M; i ++) {
            MPI_Recv (&n, 1, MPI_INT, 0, i, MPI_COMM_WORLD,&status );
        }
    }
    MPI_Finalize();
    return 0;
}
 
