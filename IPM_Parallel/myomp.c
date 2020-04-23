#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include "unistd.h"
#include <time.h>

//union_t buf;

int sum = 0, total;
int nt, mt, start, end;
int np = 1, mp = 0;
int N = 32;

int ar[32];

// ============================================== MPI
void mpierr(char *msg, const int n)
{
  puts(msg);
  MPI_Abort(MPI_COMM_WORLD,n);
}

int MyNetInit(int* argc, char*** argv, int* np, int* mp,
              int* nl, char* pname, double* tick)
{
  int i;

  i = MPI_Init(argc,argv);
  if (i != 0){
    fprintf(stderr,"MPI initialization error");
    exit(i);
  }

  MPI_Comm_size(MPI_COMM_WORLD,np);
  MPI_Comm_rank(MPI_COMM_WORLD,mp);
  MPI_Get_processor_name(pname,nl);

  *tick = MPI_Wtick();

  sleep(1);

  return 0;
}

// ============================================== MPI

int MyJob(int start, int end, int nt, int mt){

  int nn, mm;
  nn = np*nt;
  mm = nt*mp + mt;

  // sleep(1); // do some job..
  // printf("nn = %i, mm = %i\n", nn, mm);
  # pragma omp critical
  {
    printf("I am process number = %i\n", mp);
    printf("I do global thread number = %i\n", mm);
  }

  for (int i = start + mt; i < end; i += nt){
    sum += ar[i];
    //sleep(1);
    // printf("thread = %i, i = %i, sum = %i\n",mm, i, sum );
  }
  return 0;
}

int main(int argc, char *argv[]) {
  /* code */
  time_t t1;
  printf("Hello!\n\n");

  for (int j = 0; j < N ; j++){
    ar[j] = 1;
  }

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Comm_rank(MPI_COMM_WORLD,&mp);

  if (mp == 0) {
    start = 0;
    end = N/np;
  }
  else {
    start = N * mp /np;
    end = N * (mp+1) /np;
  }

  // if (mp==0) buf.ddata[0] = sum;
  // MPI_Bcast(buf.ddata,8,MPI_DOUBLE,0,MPI_COMM_WORLD);
  // if (mp>0) sum  = buf.ddata[0];


  printf("Netsize: %d, process: %d\n", np, mp );
  t1 = MPI_Wtime();

  printf("\nstart = %i, end = %i \n", start, end );
  #pragma omp parallel
  {
    int nt = omp_get_max_threads();
    int mt = omp_get_thread_num();
    // printf("PRAGMA: thread number = %i\n", mt);
    MyJob(start, end, nt, mt);
  } // end parallel


  t1 = MPI_Wtime() - t1;

  MPI_Reduce(&sum, &total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  printf("time of execution is : %ld \n", t1);
  MPI_Barrier(MPI_COMM_WORLD);
  printf("for process | %i | final sum is : %d \n", mp, sum);
  MPI_Barrier(MPI_COMM_WORLD);

  sleep(1);
  if (mp == 0) printf("total = %i\n", total );
  MPI_Finalize();


  return 0;
}
