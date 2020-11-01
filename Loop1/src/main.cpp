#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

void sum(double* res, uint32_t xSize, uint32_t ySize) {
	for (uint32_t y = 0; y < ySize; y++) {
	  for (uint32_t x = 0; x < xSize; x++) {
	    res[y*xSize + x] = sin(0.00001*res[y*xSize + x]);
	  }
	}
}


void calc(double* arr, uint32_t ySize, uint32_t xSize, int rank, int size)
{
	MPI_Bcast(&ySize, 1, MPI_INT, 0, MPI_COMM_WORLD);	
	MPI_Bcast(&xSize, 1, MPI_INT, 0, MPI_COMM_WORLD);	
	uint32_t loc_ySize = ySize / size;
	if (rank < (int)ySize % size) {
		loc_ySize++;
	}
	double* res = (double*)malloc(loc_ySize * xSize * sizeof(double));

	if (rank == 0) {
		int* displs = (int*)malloc(size * sizeof(int));
		int* sendcounts = (int*)malloc(size * sizeof(int));
		int offset = 0;
		for (int i = 0; i < size; i++) {
			displs[i] = offset;
			sendcounts[i] = ySize / size;
			if (i < (int)ySize % size) {
				sendcounts[i]++;
			}
			sendcounts[i] *= xSize;
			offset += sendcounts[i];
		}
		MPI_Scatterv(arr, sendcounts, displs, MPI_DOUBLE, res, loc_ySize * xSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		sum(res, xSize, loc_ySize);
		MPI_Gatherv(res, loc_ySize * xSize, MPI_DOUBLE, arr, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		free(displs);
		free(sendcounts);
	} else {
		MPI_Scatterv(NULL, NULL, NULL, 0, res, loc_ySize * xSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		sum(res, xSize, loc_ySize);
		MPI_Gatherv(res, loc_ySize * xSize, MPI_DOUBLE, NULL, NULL, NULL, 0, 0, MPI_COMM_WORLD);
	}
	free(res);
}

int main(int argc, char** argv)
{
  int rank = 0, size = 0, buf = 0;
  uint32_t ySize = 0, xSize = 0;
  double* arr = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
  {
    // Check arguments
    if (argc != 3)
    {
      std::cout << "[Error] Usage <inputfile> <output file>\n";
      buf = 1;
      MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Prepare input file
    std::ifstream input(argv[1]);
    if (!input.is_open())
    {
      std::cout << "[Error] Can't open " << argv[1] << " for write\n";
      buf = 1;
      MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Read arguments from input
    input >> ySize >> xSize;
    MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    arr = new double[ySize * xSize];

    for (uint32_t y = 0; y < ySize; y++)
    {
     for (uint32_t x = 0; x < xSize; x++)
      {
        input >> arr[y*xSize + x];
      }
    }
    input.close();
  } else {
    MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (buf != 0)
    {
      return 1;
    }
  }

  calc(arr, ySize, xSize, rank, size);

  if (rank == 0)
  {
    // Prepare output file
    std::ofstream output(argv[2]);
    if (!output.is_open())
    {
      std::cout << "[Error] Can't open " << argv[2] << " for read\n";
      delete arr;
      return 1;
    }
    for (uint32_t y = 0; y < ySize; y++)
    {
      for (uint32_t x = 0; x < xSize; x++)
      {
        output << " " << arr[y*xSize + x];
      }
      output << std::endl;
    }
    output.close();
    delete arr;
  }

  MPI_Finalize();
  return 0;
}
