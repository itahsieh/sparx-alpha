#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/*----------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
	int rank, size, sts = 0;

	/* Init MPI */
	sts = MPI_Init(&argc, &argv);

	/* Get rank of process */
	if(!sts) sts = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* Get size of process pool */
	if(!sts) sts = MPI_Comm_size(MPI_COMM_WORLD, &size);

	/* Cleanup MPI */
	if(!sts) sts = MPI_Finalize();

	return sts;
}

