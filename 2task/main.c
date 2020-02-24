/* /////////////////////////////////////////////////////
Пучкин Данила, 524 группа
Задание 0. Функция-барьер

Закомментированы main и операторы отладочной печати.
Также в комментариях находятся поясняющие слова.

В качестве пересылаемого сообщения выступает
номер процесса, который отправляет сообщение.
///////////////////////////////////////////////////// */


#include <stdio.h>
#include <mpi.h>

void barrier(int rank, int size)
{
  int max2 = 1, minDig = 1, from, next;
  MPI_Status status;
  for (; max2 < size; max2 <<= 1); 

  /*1st step*/
  if ((rank & 1) == 0)
  {
    for (; ((minDig & rank) == 0) && (minDig < size); minDig <<= 1)
    {
      from = rank | minDig;
      if (from < size)
      {
        MPI_Recv(&from, 1, MPI_INT, from, 13, MPI_COMM_WORLD, &status);
        //printf("Recv %d from %d\n", rank, from);
      }
    }
  }
  if (rank > 0)
  { 
    next = rank ^ minDig;
    MPI_Send(&rank, 1, MPI_INT, next, 13, MPI_COMM_WORLD);
    //printf("Send %d to %d\n", rank, next);
  }

  /*2nd step*/
  minDig = 1;
  if (rank > 0)
  {
    for (; (minDig & rank) == 0; minDig <<= 1);
    from = rank-minDig;
    MPI_Recv(&from, 1, MPI_INT, from, 13, MPI_COMM_WORLD, &status);
    //printf("Recv %d from %d\n", rank, from);
  }
  else
    minDig = max2;

  for (minDig >>= 1; minDig > 0; minDig >>= 1)
  {
    next = rank + minDig;
    if (next < size)
    {
      MPI_Send(&rank, 1, MPI_INT, next, 13, MPI_COMM_WORLD);
      //printf("Send %d to %d\n", rank, next);
    }
  }
}

void main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  //printf("Size: %d\n", size);
  //printf("Process: %d\n", rank);
  barrier(rank, size);

  MPI_Finalize();
  return;
}
