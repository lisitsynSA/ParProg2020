#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

void calc(double* arr, uint32_t ySize, uint32_t xSize, int rank, int size)
{
  int count;
  int chain_size; //Some threads might needs chain_size-2 calculations. But for easier coding, we will always calculate chain_size-1 calculations
  int mod;
  double* local_array;
  double* send_array;
  
  if(rank == 0 && size > 0)
  {
    chain_size = ySize / 4 + 1;
    if (ySize % 4 == 0)
      --chain_size;
  }
  MPI_Bcast(&chain_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if (rank == 0 && size > 0)
  {
    uint64_t startSize;
    startSize = xSize * 4;
    int div = startSize / size;
    mod = startSize % size;
    send_array = new double[(div + 1) * chain_size];
    count = div + 1;
    for(int i = 0; i < mod; ++i)
    {
      MPI_Send(&count, 1 , MPI_INT, i + 1, 0, MPI_COMM_WORLD);
      for(int j = 0; j < chain_size; ++j)
        for(int k = 0; k < count; ++k)
          send_array[j * count + k] = arr[(count - 1) + i * count + j * 4 * xSize + k];
      MPI_Send(send_array, count * chain_size, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
    }
    
    --count;
    for(int i = mod; i < size - 1; ++i)
    {
      MPI_Send(&count, 1 , MPI_INT, i + 1, 0, MPI_COMM_WORLD);
      for(int j = 0; j < chain_size; ++j)
        for(int k = 0; k < count; ++k)
          send_array[j * count + k] = arr[count + mod + i * count + j * 4 * xSize + k];
      MPI_Send(send_array, count * chain_size, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
    }

    local_array = new double[count * chain_size];
    for(int j = 0; j < chain_size; ++j)
      for(int k = 0; k < count; ++k)
        local_array[j * count + k] = arr[j * 4 * xSize + k];
  }

  if (rank != 0 && size > 0)
  {
    MPI_Recv(&count, 1 , MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    local_array = new double[count * chain_size];
    MPI_Recv(local_array, count * chain_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  for(int i = 0; i < count; ++i)
    for(int j = 0; j < chain_size - 1; ++j)
      local_array[i + (j + 1) * count] = sin(local_array[i + j * count]);

  if (rank != 0 && size > 0)
  {
    MPI_Send(local_array, count * chain_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  if (rank == 0 && size > 0)
  {
    int xySize = xSize * ySize;
    int ind;
    ++count;
    for(int i = 0; i < mod; ++i)
    {
      MPI_Recv(send_array, count * chain_size, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for(int j = 0; j < chain_size; ++j)
        for(int k = 0; k < count; ++k)
          if ((ind = (count - 1) + i * count + j * 4 * xSize + k) < xySize)
	    arr[ind] =  send_array[j * count + k];
    }
    --count;
    for(int i = mod; i < size - 1; ++i)
    {
      MPI_Recv(send_array, count * chain_size, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for(int j = 0; j < chain_size; ++j)
        for(int k = 0; k < count; ++k)
          if ((ind = count + mod + i * count + j * 4 * xSize + k) < xySize)
	    arr[ind] =  send_array[j * count + k];
    }
    for(int j = 0; j < chain_size; ++j)
      for(int k = 0; k < count; ++k)
        if((ind = j * 4 * xSize + k) < xySize)
          arr[j * 4 * xSize + k] =  local_array[j * count + k];
    delete send_array;
  }

  delete local_array;
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
