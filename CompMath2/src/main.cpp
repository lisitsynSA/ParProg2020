#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

void calc(double* frame, uint32_t ySize, uint32_t xSize, double delta, int rank, int size)
{
  int count, mod;
  double* local_array;
  MPI_Bcast(&xSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0 && size > 0)
  {
    count = (ySize - 2) / size;
    mod = (ySize - 2) % size;
    
    ++count;
    for(int i = 0; i < mod; ++i)
    {
      MPI_Send(&count, 1 , MPI_INT, i + 1, 0, MPI_COMM_WORLD);
      MPI_Send(&frame[(1 + (count - 1) + i * count - 1) * xSize], (count + 2) * xSize , MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
    }

    --count;
    for(int i = mod; i < size - 1; ++i)
    {
      MPI_Send(&count, 1 , MPI_INT, i + 1, 0, MPI_COMM_WORLD);
      MPI_Send(&frame[(1 + count + mod + i * count - 1) * xSize], (count + 2) * xSize , MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
    }

    local_array = frame;
  } 
  
  if (rank != 0 && size > 0)
  {
    MPI_Recv(&count, 1 , MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    local_array = new double[(count + 2) * xSize];
    MPI_Recv(local_array, (count + 2) * xSize , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    
  }

  double* up_string = new double[xSize - 2]; //up_string is what was upper
  double* tmp_string = new double[xSize - 1]; //tmp_string[0] is what was lefter
  int is_calculate = 1;
  double diff;
  double diff_sum;
  double* diff_recv;
  if (rank == 0 && size > 0)
    diff_recv = new double[size];
  while(is_calculate)
  {
    diff = 0;
    for(uint32_t x = 1; x < xSize - 1; ++x)
      up_string[x - 1] = local_array[x];
    for(int y = 1; y < count + 1; ++y)
    {
      tmp_string[0] = local_array[y * xSize];
      for(uint32_t x = 1; x < xSize - 1; ++x)
      {
        tmp_string[x] = local_array[y*xSize + x];
        local_array[y*xSize + x] = (up_string[x - 1] + local_array[(y + 1)*xSize + x] +\
                                    local_array[y*xSize + x + 1] + tmp_string[x - 1])/4.0;
        up_string[x - 1] = tmp_string[x];
        diff += std::abs(local_array[y*xSize + x] - tmp_string[x]);
      }
    }

    //0--->ySize
    if (rank != 0 && size > 0)
      MPI_Recv(&local_array[1], xSize - 2 , MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    
    if ((rank != size - 1) && size > 0)
      MPI_Send(&local_array[count * xSize + 1], xSize - 2 , MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    //0<---ySize
    if ((rank != size - 1) && size > 0)
      MPI_Recv(&local_array[(count + 1) * xSize + 1], xSize - 2 , MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    
    if (rank != 0 && size > 0)
      MPI_Send(&local_array[1 * xSize + 1], xSize - 2 , MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
    
    MPI_Gather(&diff, 1, MPI_DOUBLE, diff_recv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0 && size > 0)
    {
      diff_sum = 0;
      for(int i = 0; i < size; ++i)
        diff_sum += diff_recv[i];
      if (diff_sum <= delta)
        is_calculate = 0;
    }
    MPI_Bcast(&is_calculate, 1, MPI_INT, 0, MPI_COMM_WORLD);

  }
  if (rank != 0 && size > 0)
  {
    MPI_Send(&local_array[1 * xSize], xSize * count , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
  if (rank == 0 && size > 0)
  {
    ++count;
    for(int i = 0; i < mod; ++i)
      MPI_Recv(&frame[(1 + (count - 1) + i * count) * xSize], count * xSize , MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    --count;
    for(int i = mod; i < size - 1; ++i)
      MPI_Recv(&frame[(1 + count + mod + i * count) * xSize], count * xSize, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
    
  if (rank == 0 && size > 0)
    delete diff_recv;
  if (rank != 0 && size > 0)
    delete local_array;
  delete up_string;
  delete tmp_string;
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
