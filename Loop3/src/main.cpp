#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

void sum(double* arr, uint32_t xSize, uint32_t ySize, uint32_t start, int& iter) {
  iter = 0;
  for (uint32_t y = 4 + start; y < ySize; y += 4) {
    iter++;
    for (uint32_t x = 0; x < xSize; x++) {
      arr[y*xSize + x] = sin(arr[(y - 4)*xSize + x]);
    }
  }
}

void calc(double* arr, uint32_t ySize, uint32_t xSize, int rank, int size)
{
  if (rank < 4) {
    MPI_Status status;
    MPI_Bcast(&ySize, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast(&xSize, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    uint32_t loc_ySize = 0, loc_yAddr = 0;

    if (rank ==0) {
      loc_yAddr = 0;
      loc_ySize = 4 / size;
      if (rank < (int)(4 % size))
        loc_ySize++;
      uint32_t send_ySize = loc_ySize;
      uint32_t send_yAddr = loc_ySize;
      if (size >= 4) {
        size = 4;
      }

      for (int i = 1; i < size; i++) {
        if (i == (int)(4 % size)) {
          send_ySize--;
        } 
        MPI_Send(&send_ySize, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        MPI_Send(&send_yAddr, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        MPI_Send(&arr[send_yAddr * xSize], send_ySize * xSize, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        send_yAddr += send_ySize;
      }
    } else {
      arr = (double*) malloc(ySize * xSize * sizeof(double));
      MPI_Recv(&loc_ySize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&loc_yAddr, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&arr[loc_yAddr * xSize], loc_ySize * xSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    }
    int iter;

    for (uint32_t i = 0; i < loc_ySize; i++) {
      sum(arr, xSize, ySize, loc_yAddr + i, iter);
    }


    if (rank == 0) {
      uint32_t send_ySize = loc_ySize;
      uint32_t send_yAddr = loc_ySize;
      for (int i = 1; i < size; i++) {
        if (i == (int)(4 % size)) {
          send_ySize--;
        } 
        MPI_Recv(&iter, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
        MPI_Datatype floattype;
        MPI_Type_vector(iter, xSize * send_ySize, xSize * 4, MPI_DOUBLE, &floattype);
        MPI_Type_commit(&floattype);
        MPI_Recv(&arr[(send_yAddr + 4) * xSize], 1, floattype, i, 0, MPI_COMM_WORLD, &status);
        send_yAddr += send_ySize;
      }
    } else {
      MPI_Datatype floattype;
      MPI_Type_vector(iter, xSize * loc_ySize, xSize * 4, MPI_DOUBLE, &floattype);
      MPI_Type_commit(&floattype);
      MPI_Send(&iter, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&arr[(loc_yAddr + 4) * xSize], 1, floattype, 0, 0, MPI_COMM_WORLD);
    }

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
