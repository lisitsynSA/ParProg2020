#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

//arr[y*xSize + x] = sin(0.00001*arr[(y + 1)*xSize + x - 3]);
void calc(double* arr, uint32_t ySize, uint32_t xSize, int rank, int size)
{
  int count;
  int mod;
  double* local_array;
  double* send_array;
  if (rank == 0 && size > 0)
  {
    uint64_t xySize;
    xSize -= 3;
    ySize -= 1;
    xySize = xSize * ySize;
    int div = xySize / size;
    mod = xySize % size;
    send_array = new double[xySize];
    
    for(uint32_t i = 0; i < ySize; ++i)
      for(uint32_t j = 0; j < xSize; ++j)
        send_array[i * xSize + j] = arr[(i + 1) * (xSize + 3) + j];
    
    count = div + 1;
    for(int i = 0; i < mod; ++i)
    {
      MPI_Send(&count, 1 , MPI_INT, i + 1, 0, MPI_COMM_WORLD);
      MPI_Send(&send_array[(count - 1) + i * count], count , MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
    }
    
    --count;
    for(int i = mod; i < size - 1; ++i)
    {
      MPI_Send(&count, 1 , MPI_INT, i + 1, 0, MPI_COMM_WORLD);
      MPI_Send(&send_array[count + i * count + mod], count, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
    }

    local_array = send_array;
  }
  if (rank != 0 && size > 0)
  {
    MPI_Recv(&count, 1 , MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    local_array = new double[count];
    MPI_Recv(local_array, count, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  
  for(int i = 0; i < count; ++i)
    local_array[i] = sin(0.00001*(local_array[i]));

  if (rank != 0 && size > 0)
  {
    MPI_Send(local_array, count, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    delete local_array;
  }

  if (rank == 0 && size > 0)
  {
    ++count;
    for(int i = 0; i < mod; ++i)
      MPI_Recv(&send_array[(count - 1) + i * count], count, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    --count;
    for(int i = mod; i < size - 1; ++i)
      MPI_Recv(&send_array[count + i * count + mod], count, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for(uint32_t i = 0; i < ySize; ++i)
      for(uint32_t j = 0; j < xSize; ++j)
        arr[i * (xSize + 3) + (j + 3)] = send_array[i * xSize + j];
    delete send_array;
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
