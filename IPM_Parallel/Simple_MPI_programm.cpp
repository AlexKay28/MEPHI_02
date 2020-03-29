#include <stdio.h>
#include <string.h>
#include "mpi.h"

#define BUF_LEN 256

int main(int argc, char *argv[]) {
	int my_rank;
	int p;
	int sourse;
	int dest;
	int tag = 0;
	char message[BUF_LEN];
	MPI_Status status;

	MPI_Init (&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

	// Process 0
	if (my_rank != 0) {
		sprintf (message, "Hello from process %d!", my_rank);
		dest = 0;
		MPI_Send (message, strlen(message) + 1, MPI_CHAR, dest, tag, MPI_COMM_WORD);
	} 
	// All others
	else {
		for (sourse = 1; sourse < p; sourse++) {
			MPI_Recv(message, BUF_LEN, MPI_CHAR, sourse, tag, MPI_COMM_WORD, &status);
			sprintf("%s\n", message);
		}
	}

	MPI_Finalize();
	return 0;
}