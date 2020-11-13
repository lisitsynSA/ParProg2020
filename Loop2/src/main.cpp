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
    int* send_dist = NULL;
    int* send_range = NULL;
    
    int* recv_dist = NULL;
    int* recv_range = NULL;

    MPI_Bcast(&ySize, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&xSize, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD);
    int real_rank_size = 0;

    if(rank == ROOT) {
        size_t my_calc_size = ceil((ySize) * 1.0 / size);
        // if  proc_rank > ceil(ySize / my_calc_size) it will do nothing
        real_rank_size = ceil((ySize) * 1.0 / my_calc_size);
        size_t now_size = 0; 
    
        recv_dist = (int*)calloc(size, sizeof(int));
        recv_range = (int*)calloc(size, sizeof(int));

        send_dist = (int*)calloc(size, sizeof(int));
        send_range = (int*)calloc(size, sizeof(int));

        for(int i = 0; i < real_rank_size; ++i) {
            recv_dist[i] = i * xSize * my_calc_size;
            recv_range[i] = xSize * std::min(my_calc_size, ySize - now_size);
            now_size += std::min(my_calc_size, ySize - now_size);
        }
        
        for(int i = 0; i < real_rank_size; ++i) {
            send_range[i] = (i == real_rank_size - 1)? recv_range[i] : recv_range[i] + xSize;
            send_dist[i] = recv_dist[i]; 
        }
    }    
    int send_step = 0;
    int recv_step = 0;


    // Bug with MPI <- i hate it! really!
    MPI_Scatter(
            send_range, 1, MPI_INT,
            &send_step, size, MPI_INT, 
            ROOT, MPI_COMM_WORLD);
    MPI_Scatter(
            recv_range, 1, MPI_INT,
            &recv_step, size, MPI_INT, 
            ROOT, MPI_COMM_WORLD);

    double* my_calc = NULL;
    my_calc = (double*)calloc(send_step, sizeof(double));    

    MPI_Scatterv(
            arr, send_range, send_dist, MPI_DOUBLE,
            my_calc, send_step, MPI_DOUBLE,
            ROOT, MPI_COMM_WORLD);
    
    if((uint32_t)send_step > xSize) {
        for(size_t y = 0; y < send_step / xSize - 1; ++y) {
            for(size_t x = 3; x < xSize; ++x) {
                my_calc[y*xSize + x] = sin(0.00001*my_calc[(y + 1) * xSize + (x - 3)]);
            }
        }
    }

    MPI_Gatherv(
            my_calc, recv_step, MPI_DOUBLE,
            arr, recv_range, recv_dist, MPI_DOUBLE, 
            ROOT, MPI_COMM_WORLD);    


    free(my_calc);
    if(rank == ROOT) {
        free(recv_dist);
        free(recv_range);
        free(send_dist);
        free(send_range);
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
