#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

void calc(double* frame, uint32_t ySize, uint32_t xSize, double delta, int rank, int size)
{
	MPI_Status status;
	MPI_Bcast (&xSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast (&ySize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast (&delta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	uint32_t mySize = ySize / size;
	uint32_t yFirst = mySize * rank + ySize % size;

	if (rank == 0) {
		yFirst = 1;
		mySize += ySize % size - 1;
	}

	if (rank == size - 1) {
		mySize--;
	}

	uint32_t yLast = yFirst + mySize;

	double diff = 0;

	if (mySize == 0) {
		diff = 2 * delta;
	}

	double* inFrame = (double *)malloc(xSize * ySize * sizeof(double));

	if (rank > 0) {
		frame = (double *)malloc(xSize * ySize * sizeof(double));
	}

	MPI_Bcast (frame, xSize * ySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	memcpy(inFrame, frame, xSize * ySize * sizeof (double));

	for (uint32_t y = yFirst; y < yLast; y++) {
		for (uint32_t x = 1; x < xSize - 1; x++) {
			inFrame[y * xSize + (x)] = (frame[(y + 1) * xSize + x] + frame[(y - 1) * xSize + x] + frame[y * xSize + x + 1] + frame[y * xSize + x - 1]) / 4.0;
			diff += std::abs(inFrame[y * xSize + x] - frame[y * xSize + x]);
		}
	}

	double* p;

	while (diff > delta)
	{
		p = inFrame;
		inFrame = frame;
		frame = p;

		diff = 0;

		if (size > 1) {
			if (rank == 0) {
				MPI_Send (inFrame + (yLast - 1) * xSize, xSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
			} else {
				MPI_Recv (inFrame + (yFirst - 1) * xSize, xSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
				if (rank != size - 1) {
					MPI_Send (inFrame + (yLast - 1) * xSize, xSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
				}
			}
			if (rank == size - 1) {
				MPI_Send (inFrame + yFirst * xSize, xSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
			} else {
				MPI_Recv (inFrame + yLast * xSize, xSize, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
				if (rank != 0) {
					MPI_Send (inFrame + yFirst * xSize, xSize, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
				}
			}
		}



		for (uint32_t y = yFirst; y < yLast; y++) {
			for (uint32_t x = 1; x < xSize - 1; x++) {
				frame[y * xSize + x] = (inFrame[(y + 1) * xSize + x] + inFrame[(y - 1) * xSize + x] + inFrame[y * xSize + x + 1] + inFrame[y * xSize + x - 1]) / 4.0;
				diff += std::abs(frame[y * xSize + x] - inFrame[y * xSize + x]);
			}
		}

		MPI_Allreduce(&diff, &diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	}

	if (rank == 0) {
		yFirst += mySize;
		mySize -= (ySize % size - 1);
		for (int i = 1; i < size; i++) {
			if (i == size - 1) {
				mySize--;
			}
			MPI_Recv (frame + yFirst * xSize, mySize * xSize, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
			yFirst += mySize;
		}
	} else {
		MPI_Send (frame + yFirst * xSize, mySize * xSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	free (inFrame);

	if (rank > 0) {
		free (frame);
	}
}


int main(int argc, char** argv)
{
  int rank = 0, size = 0, status = 0;
  double delta = 0;
  uint32_t ySize = 0, xSize = 0;
  double* frame = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
  {
    // Check arguments
    if (argc != 3)
    {
      std::cout << "[Error] Usage <inputfile> <output file>\n";
      status = 1;
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Prepare input file
    std::ifstream input(argv[1]);
    if (!input.is_open())
    {
      std::cout << "[Error] Can't open " << argv[1] << " for write\n";
      status = 1;
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Read arguments from input
    input >> ySize >> xSize >> delta;
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);

    frame = new double[ySize * xSize];

    for (uint32_t y = 0; y < ySize; y++)
    {
     for (uint32_t x = 0; x < xSize; x++)
      {
        input >> frame[y*xSize + x];
      }
    }
    input.close();
  } else {
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (status != 0)
    {
      return 1;
    }
  }

  calc(frame, ySize, xSize, delta, rank, size);

  if (rank == 0)
  {
    // Prepare output file
    std::ofstream output(argv[2]);
    if (!output.is_open())
    {
      std::cout << "[Error] Can't open " << argv[2] << " for read\n";
      delete frame;
      return 1;
    }
    for (uint32_t y = 0; y < ySize; y++)
    {
      for (uint32_t x = 0; x < xSize; x++)
      {
        output << " " << frame[y*xSize + x];
      }
      output << std::endl;
    }
    output.close();
    delete frame;
  }

  MPI_Finalize();
  return 0;
}
