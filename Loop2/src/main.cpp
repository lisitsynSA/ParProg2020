#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

void sum(double* res, uint32_t xSize, uint32_t ySize) {
  for (uint32_t y = 0; y < ySize - 1; y++) {
    for (uint32_t x = 3; x < xSize; x++) {
      res[y*xSize + x] = sin(0.00001*res[(y+1)*xSize + x - 3]);
    }
  }
}


void calc(double* arr, uint32_t ySize, uint32_t xSize, int rank, int size)
{
  MPI_Status status;
  MPI_Bcast(&ySize, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast(&xSize, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  uint32_t loc_ySize = ySize / size;
  if (rank < (int)(ySize % size)){
    loc_ySize++;
  }
  if (rank == size - 1) {
    loc_ySize--;
  }
  if (rank == 0) {
    uint32_t send_ySize = loc_ySize;
    uint32_t send_yAddr = loc_ySize;
    for (int i = 1; i < size; i++) {
      if (i == (int)(ySize % size)) {
        send_ySize--;
      }
      if (i == size - 1) {
        send_ySize--;
      }
      MPI_Send(&arr[send_yAddr * xSize], (send_ySize + 1) * xSize, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
      send_yAddr += send_ySize;
    }
  } else {
    arr = (double*) malloc((loc_ySize + 1) * xSize * sizeof(double));
    MPI_Recv(arr, (loc_ySize + 1) * xSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
  }
  sum(arr, xSize, loc_ySize + 1);
  if (rank == 0) {
    uint32_t recv_ySize = loc_ySize;
    uint32_t recv_yAddr = loc_ySize;
    for (int i = 1; i < size; i++) {
      if (i == (int)(ySize % size)) {
        recv_ySize--;
      }
      if (i == size - 1) {
        recv_ySize--;
      }
      MPI_Recv(&arr[recv_yAddr * xSize], recv_ySize * xSize, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
      recv_yAddr += recv_ySize;
    }
  } else {
    MPI_Send(arr, loc_ySize * xSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    free(arr);
  }
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
