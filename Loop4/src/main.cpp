#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>
#include <limits> //Add this for NaN in border_checker

double border_checker(double* array, int xSize, int ySize, int zSize, int x, int y, int z)
{
  if((x < xSize) && (x >= 0) && (y < ySize) && (y >= 0) && (z < zSize) && (z >= 0))
    return array[z * xSize * ySize + y * xSize + x];
  else
    return std::numeric_limits<double>::quiet_NaN();
}

void calc(double* arr, uint32_t zSize, uint32_t ySize, uint32_t xSize, int rank, int size)
{
  int count;
  int chain_size; //Some threads might needs chain_size-2 calculations. But for easier coding, we will always calculate chain_size-1 calculations
  int mod;
  double* local_array;
  double* send_array;
  int* root_x;
  int* root_y;
  int* root_z;
  if(rank == 0 && size > 0)
  {
    chain_size = std::min(std::min(xSize, ySize), zSize); //it's max chain size
  }
  MPI_Bcast(&chain_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if (rank == 0 && size > 0)
  { //roots: z = 0; x = xSize; y = ySize; 3 sides, sideSize - 1
    uint64_t startSize = (xSize - 1) * (ySize - 1) + (xSize - 1) * (zSize - 2) + (ySize - 2) * (zSize - 2);
    root_x = new int[startSize];
    root_y = new int[startSize];
    root_z = new int[startSize];
    int index = 0;
    for(uint32_t i = 1; i < xSize; ++i)
      for(uint32_t j = 1; j < ySize; ++j)
      {
        root_x[index] = i;
	root_y[index] = j;
	root_z[index] = 0;
	++index;
      }
    for(uint32_t i = 1; i < xSize; ++i)
      for(uint32_t j = 1; j < zSize - 1; ++j)
      {
        root_x[index] = i;
	root_y[index] = ySize - 1;
	root_z[index] = j;
	++index;
      }
    for(uint32_t i = 1; i < ySize - 1; ++i)
      for(uint32_t j = 1; j < zSize - 1; ++j)
      {
        root_x[index] = xSize - 1;
	root_y[index] = i;
	root_z[index] = j;
	++index;
      }
    int div = startSize / size;
    mod = startSize % size;
    send_array = new double[(div + 1) * chain_size];
    count = div + 1;
    for(int i = 0; i < mod; ++i)
    {
      MPI_Send(&count, 1 , MPI_INT, i + 1, 0, MPI_COMM_WORLD);
      for(int j = 0; j < chain_size; ++j)
        for(int k = 0; k < count; ++k)
          send_array[j * count + k] = border_checker(arr, xSize, ySize, zSize, 
			  root_x[(count - 1) + i * count + k] - j * 1,
			  root_y[(count - 1) + i * count + k] - j * 1,
			  root_z[(count - 1) + i * count + k] + j * 1);
      MPI_Send(send_array, count * chain_size, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
    }
    
    --count;
    for(int i = mod; i < size - 1; ++i)
    {
      MPI_Send(&count, 1 , MPI_INT, i + 1, 0, MPI_COMM_WORLD);
      for(int j = 0; j < chain_size; ++j)
        for(int k = 0; k < count; ++k)
	  send_array[j * count + k] = border_checker(arr, xSize, ySize, zSize, 
			  root_x[count + mod + i * count + k] - j * 1,
			  root_y[count + mod + i * count + k] - j * 1,
			  root_z[count + mod + i * count + k] + j * 1);

      MPI_Send(send_array, count * chain_size, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
    }

    local_array = new double[count * chain_size];
    for(int j = 0; j < chain_size; ++j)
      for(int k = 0; k < count; ++k)
        local_array[j * count + k] = border_checker(arr, xSize, ySize, zSize, 
			  root_x[k] - j * 1,
			  root_y[k] - j * 1,
			  root_z[k] + j * 1);

  }

  if (rank != 0 && size > 0)
  {
    MPI_Recv(&count, 1 , MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    local_array = new double[count * chain_size];
    MPI_Recv(local_array, count * chain_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  for(int i = 0; i < count; ++i)
    for(int j = 0; j < chain_size - 1; ++j)
    {
      if(local_array[i + (j + 1) * count] == local_array[i + (j + 1) * count]) //Check if not NaN
        local_array[i + (j + 1) * count] = sin(local_array[i + j * count]);
      else
        j = chain_size; //End calculating for this chain
    }

  if (rank != 0 && size > 0)
  {
    MPI_Send(local_array, count * chain_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  if (rank == 0 && size > 0)
  {
    double tmp;
    ++count;
    for(int i = 0; i < mod; ++i)
    {
      MPI_Recv(send_array, count * chain_size, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for(int k = 0; k < count; ++k)
        for(int j = 0; j < chain_size; ++j)  
        {
          if((tmp = send_array[j * count + k]) == (send_array[j * count + k])) //Check if not NaN
            arr[(root_x[(count - 1) + i * count + k] - j * 1) +
                (root_y[(count - 1) + i * count + k] - j * 1) * xSize +
		(root_z[(count - 1) + i * count + k] + j * 1) * xSize * ySize] = tmp;
	  else
            j = chain_size;
        }
    }
    --count;
    for(int i = mod; i < size - 1; ++i)
    {
      MPI_Recv(send_array, count * chain_size, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for(int k = 0; k < count; ++k)
        for(int j = 0; j < chain_size; ++j)  
        {
          if((tmp = send_array[j * count + k]) == (send_array[j * count + k])) //Check if not NaN
            arr[(root_x[count + mod + i * count + k] - j * 1) +
                (root_y[count + mod + i * count + k] - j * 1) * xSize +
		(root_z[count + mod + i * count + k] + j * 1) * xSize * ySize] = tmp;
	  else
            j = chain_size;
        }

    }
    for(int k = 0; k < count; ++k)
      for(int j = 0; j < chain_size; ++j)  
      {
        if((tmp = local_array[j * count + k]) == (local_array[j * count + k])) //Check if not NaN
          arr[(root_x[k] - j * 1) +
              (root_y[k] - j * 1) * xSize +
              (root_z[k] + j * 1) * xSize * ySize] = tmp;
        else
          j = chain_size;
      }

    delete send_array;
    delete root_x;
    delete root_y;
    delete root_z;
  }

  delete local_array;
}

int main(int argc, char** argv)
{
  int rank = 0, size = 0, buf = 0;
  uint32_t zSize = 0, ySize = 0, xSize = 0;
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
    input >> zSize >> ySize >> xSize;
    MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    arr = new double[zSize * ySize * xSize];
    for (uint32_t z = 0; z < zSize; z++) {
      for (uint32_t y = 0; y < ySize; y++) {
        for (uint32_t x = 0; x < xSize; x++) {
          input >> arr[z*ySize*xSize + y*xSize + x];
        }
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

  calc(arr, zSize, ySize, xSize, rank, size);

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

    for (uint32_t z = 0; z < zSize; z++) {
      for (uint32_t y = 0; y < ySize; y++) {
        for (uint32_t x = 0; x < xSize; x++) {
          output << " " << arr[z*ySize*xSize + y*xSize + x];
        }
        output << std::endl;
      }
      output << std::endl;
    }
    output.close();
    delete arr;
  }

  MPI_Finalize();
  return 0;
}
