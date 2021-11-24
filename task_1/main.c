#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 5

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, tasks;
    MPI_Comm comm;
    //Функция определения номера процесса
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Функция определения числа процессов в области связи
    MPI_Comm_size(MPI_COMM_WORLD, &tasks);

    int size[2] = {SIZE, SIZE};
    int periodic[2] = {0};
    // создание транспьютерной матрицы
    // Функция создания коммуникатора с декартовой топологией
    // MPI_COMM_WORL - родительский коммуникатор,
    // 2 - число измерений 
    // size - массив размера ndims, в котором задается число процессов вдоль каждого измерения;
    // comm - новый коммуникатор
    MPI_Cart_create(MPI_COMM_WORLD, 2, size, periodic, 0, &comm);
    int coords[2];
    // Функция определения координат процесса по его идентификатору
    // comm - коммуникатор с декартовой топологией 
    // rank - идентификатор процесса 
    // 2 - число измерений
    // coords - координаты процесса в декартовой топологии
    MPI_Cart_coords(comm, rank, 2, coords);
    printf("Coordinates for process %d: (%d, %d)\n", rank, coords[0], coords[1]);
    MPI_Finalize();
    return 0;
}