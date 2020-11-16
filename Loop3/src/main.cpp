#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

#define ROOT 0

void calc(double* arr, uint32_t ySize, uint32_t xSize, int rank, int size)
{
    double* new_arr = NULL;
    if(rank == ROOT) {
        new_arr = (double*)calloc(xSize * ySize, sizeof(double));
        for(size_t y = 0; y < ySize; ++y) {
            for(size_t x = 0; x < xSize; ++x) {
                new_arr[x*ySize + y] = arr[y*xSize + x];
            }
        }
    }

    int32_t step = 0;
    int32_t* range = NULL;
    int32_t* displace = NULL;
    
    MPI_Bcast(&xSize, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&ySize, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD);
    
    if(rank == ROOT) {        
        uint32_t x_step = ceil(xSize / (1.0 * size));
        uint32_t prev_range = 0;
        
        range = (int32_t*)calloc(size, sizeof(int32_t));
        displace = (int32_t*)calloc(size, sizeof(int32_t));

        for(int send_rank = 0; send_rank < size; ++send_rank) {
            // We will send start_y, step, xSize and part of arr
            uint32_t real_step = std::min(x_step, xSize - prev_range);
            range[send_rank] = real_step * ySize;
            displace[send_rank] = prev_range * ySize;
            prev_range += real_step; 
        }
    }
    MPI_Scatter(range, 1, MPI_INT,
            &step, 1, MPI_INT,
            ROOT, MPI_COMM_WORLD);

    double* my_copy = NULL;
    my_copy = (double*) calloc(step, sizeof(double));
    MPI_Scatterv(new_arr, range, displace, MPI_DOUBLE, 
            my_copy, step, MPI_DOUBLE, 
            ROOT, MPI_COMM_WORLD);    
    for(size_t x = 0; x < step / ySize; ++x) {
        for(size_t y = 4; y < ySize; ++y) {
            my_copy[x*ySize + y] = sin(my_copy[x*ySize + y - 4]);
        }
    }

    MPI_Gatherv(my_copy, step, MPI_DOUBLE, 
            new_arr, range, displace, MPI_DOUBLE,
            ROOT, MPI_COMM_WORLD);

    if(rank == ROOT) {
        for(size_t y = 0; y < ySize; ++y) {
            for(size_t x = 0; x < xSize; ++x) {
                arr[y*xSize + x] = new_arr[x*ySize + y];
            }
        }
    }
    
    
    if(rank == ROOT) {
        free(new_arr);
        free(range);
        free(displace);
    }
    free(my_copy); 
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
