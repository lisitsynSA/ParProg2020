#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

#define ROOT 0

void form_arr(double* arr, double* copy_arr, uint32_t ySize, uint32_t xSize, int* ran) {
    int pos = 0;
    for(int i = 0; i < 4; ++i) {
        for(uint32_t y = i; y < ySize; y += 4) {
            memcpy(&copy_arr[(pos++) * xSize], &arr[y * xSize], xSize * sizeof(double));
            ran[i] += xSize;
        }
    }
}

void unform_arr(double* arr, double* copy_arr, uint32_t ySize, uint32_t xSize) {
    int pos = 0;
    for(int i = 0; i < 4; ++i) {
        for(uint32_t y = i; y < ySize; y += 4) {
            memcpy(&arr[y * xSize], &copy_arr[(pos++) * xSize], xSize * sizeof(double));
        }
    }
}

void calc(double* arr, uint32_t ySize, uint32_t xSize, int rank, int size)
{
    MPI_Comm comm;
    if(size < 4) {
        if(rank == ROOT) {
            for (uint32_t y = 4; y < ySize; y++) {
                for (uint32_t x = 0; x < xSize; x++) {
                    arr[y*xSize + x] = sin(arr[(y - 4)*xSize + x]);
                }
            }
        }
        return;
    }

    MPI_Comm_split(MPI_COMM_WORLD, (rank < 4) ? 1 : MPI_UNDEFINED, rank, &comm);
    if(rank >= 4) {
        return; 
    }

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    int* ran = NULL;
    int* dis = NULL;
    double* new_array = NULL;
    int step = 0;

    MPI_Bcast(&xSize, 1, MPI_UNSIGNED, ROOT, comm);
    MPI_Bcast(&ySize, 1, MPI_UNSIGNED, ROOT, comm);
    /*
    if(rank == ROOT) {
        printf("%d\n", size);
    }*/
    
    if(rank == ROOT) { 
        ran = (int*)calloc(4, sizeof(int));
        dis = (int*)calloc(4, sizeof(int));
        new_array = (double*)calloc(xSize*ySize, sizeof(double));
        form_arr(arr, new_array, ySize, xSize, ran);
        dis[1] = ran[0];
        dis[2] = ran[1] + dis[1];
        dis[3] = ran[2] + dis[2];
    }

    MPI_Scatter(ran, 1, MPI_INT,
            &step, 1, MPI_INT,
            ROOT, comm);
    //printf("[%d] task = %d, step = %d\n", rank, ySize * xSize, step);

    double* my_copy = (double*)calloc(step, sizeof(double));
    MPI_Scatterv(new_array, ran, dis, MPI_DOUBLE, 
            my_copy, step, MPI_DOUBLE, 
            ROOT, comm);

    for(int x = xSize; x < step; ++x) {
        my_copy[x] = sin(my_copy[x - xSize]);
    }
    
    MPI_Gatherv(my_copy, step, MPI_DOUBLE, 
            new_array, ran, dis, MPI_DOUBLE,
            ROOT, comm);

    if(rank == ROOT) {
        unform_arr(arr, new_array, ySize, xSize);
        free(new_array);
        free(ran);
        free(dis);
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
