#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

#include <cstdlib>

#define ROOT 0
// auntidep
// y = 1; y < ySize;
// x = 0; y < xSize - 3
//arr[(y - 1)*xSize + (x + 3)] = sin(0.00001*arr[y*xSize + x]);

void calc(double* arr, uint32_t ySize, uint32_t xSize, int rank, int size) {
    int32_t step = 0;
    int32_t* range = NULL;
    int32_t* displace = NULL;
    
    MPI_Bcast(&xSize, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&ySize, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD);
    
    if(rank == ROOT) {        
        uint32_t y_step = ceil(ySize / (1.0 * size)) + 1;
        uint32_t prev_range = 0;
        
        range = (int32_t*)calloc(size, sizeof(int32_t));
        displace = (int32_t*)calloc(size, sizeof(int32_t));
        // реализуем перекрытие
        for(int send_rank = 0; send_rank < size; ++send_rank) {
            // We will send start_y, step, xSize and part of arr
            uint32_t real_step = std::min(y_step, ySize - prev_range);
            range[send_rank] = real_step * xSize;
            displace[send_rank] = prev_range * xSize;
            prev_range += y_step - 1; 
        }
    }
    MPI_Scatter(range, 1, MPI_INT,
            &step, 1, MPI_INT,
            ROOT, MPI_COMM_WORLD);
    

    double* my_copy = (double*) calloc(step, sizeof(double));
    MPI_Scatterv(arr, range, displace, MPI_DOUBLE, 
            my_copy, step, MPI_DOUBLE, 
            ROOT, MPI_COMM_WORLD);    
    
    if(rank == ROOT) {
        for(int i = 0; i < size; ++i) {
            range[i] -= (range[i] == 0) ? 0 : xSize;
        }
    }

    if(step != 0) {
        for(uint32_t y = 0; y < step / xSize - 1; ++y) {
            for(uint32_t x = 3; x < xSize; ++x) {
                my_copy[y*xSize + x] = sin(0.00001*my_copy[(y + 1) * xSize + (x - 3)]);
            }
        }
    }
 
    MPI_Gatherv(my_copy, (step == 0)?0:step - xSize, MPI_DOUBLE, 
            arr, range, displace, MPI_DOUBLE,
            ROOT, MPI_COMM_WORLD);

    if(rank == ROOT) {

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
