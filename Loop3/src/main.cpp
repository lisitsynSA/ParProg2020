#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

void sum(double* arr, uint32_t xSize, uint32_t ySize) {
	for (uint32_t x = 0; x < xSize; x++) {
		for (uint32_t y = 4; y < ySize; y++) {
	  		arr[x*ySize + y] = sin(arr[x*ySize + (y-4)]);
		}
	}
}

void calc(double* arr, uint32_t ySize, uint32_t xSize, int rank, int size)
{
	MPI_Bcast(&ySize, 1, MPI_INT, 0, MPI_COMM_WORLD); 
	MPI_Bcast(&xSize, 1, MPI_INT, 0, MPI_COMM_WORLD); 
	uint32_t loc_xSize = xSize / size;
	if (rank < (int)xSize % size) {
		loc_xSize++;
	}
	double* res = (double*)malloc(loc_xSize * ySize * sizeof(double));

	if (rank == 0) {
		double* arr_T = (double *)malloc(ySize * xSize * sizeof(double));
		for(int i = 0; i < (int)ySize; i++) {
			for(int j = 0; j < (int)xSize; j++) {
				arr_T[j*ySize + i] = arr[i*xSize + j];
			}
		}

		int* displs = (int*)malloc(size * sizeof(int));
		int* sendcounts = (int*)malloc(size * sizeof(int));
		int offset = 0;
		for (int i = 0; i < size; i++) {
			displs[i] = offset;
			sendcounts[i] = xSize / size;
			if (i < (int)xSize % size) {
				sendcounts[i]++;
			}
			sendcounts[i] *= ySize;
			offset += sendcounts[i];
		}
		MPI_Scatterv(arr_T, sendcounts, displs, MPI_DOUBLE, res, loc_xSize * ySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		sum(res, loc_xSize, ySize);
		MPI_Gatherv(res, loc_xSize * ySize, MPI_DOUBLE, arr_T, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for(int i = 0; i < (int)xSize; i++) {
			for(int j = 0; j < (int)ySize; j++) {
				arr[j*xSize + i] = arr_T[i*ySize + j];
			}
		}
		free(displs);
		free(sendcounts);
		free(arr_T);
	} else {
		MPI_Scatterv(NULL, NULL, NULL, 0, res, loc_xSize * ySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		sum(res, loc_xSize, ySize);
		MPI_Gatherv(res, loc_xSize * ySize, MPI_DOUBLE, NULL, NULL, NULL, 0, 0, MPI_COMM_WORLD);
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
