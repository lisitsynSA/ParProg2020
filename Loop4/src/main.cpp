#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

#define ROOT 0

// arr[z*ySize*xSize + y*xSize + x] = sin(arr[(z - 1)*ySize*xSize + (y + 1)*xSize + x + 1]);
void calc(double* arr, uint32_t zSize, uint32_t ySize, uint32_t xSize, int rank, int size)
{
    int32_t step = 0;
    int* range = NULL;
    int* displace = NULL;
    int* recv_range = NULL;
    
    MPI_Bcast(&xSize, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&ySize, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&zSize, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD);
     
    if(rank == ROOT) {
        uint32_t y_step = ceil(ySize / (1.0 * size)) + 1;
        uint32_t prev_range = 0;
        range = (int*)calloc(size, sizeof(int));
        displace = (int*)calloc(size, sizeof(int));

        for(int send_rank = 0; send_rank < size; ++send_rank) {
            uint32_t real_step = std::min(y_step, ySize - prev_range);
            range[send_rank] = real_step * xSize;
            displace[send_rank] = prev_range * xSize;
            prev_range += y_step;         
        }

        recv_range = (int*)calloc(size, sizeof(int));
        for(int i = 0; i < size; ++i) {
            recv_range[i] = range[i] - (range[i] == 0) ? 0 : xSize;
        }

    }

    MPI_Scatter(range, 1, MPI_INT,
            &step, 1, MPI_INT,
            ROOT, MPI_COMM_WORLD);

    double* my_copy = (double*)calloc(step, sizeof(double));
    

    for(size_t z = 1; z < zSize; ++z) {
        MPI_Scatterv(arr + (z - 1)*ySize*xSize, range, displace, MPI_DOUBLE, 
                my_copy, step, MPI_DOUBLE, 
                ROOT, MPI_COMM_WORLD);    

        if(step != 0) {
            for(uint32_t y = 0; y < step / xSize - 1; ++y) {
                for(uint32_t x = 0; x < xSize - 1; ++x) {
                    my_copy[y*xSize + x] = sin(my_copy[(y + 1) * xSize + (x + 1)]);
                }
            }
        }

        MPI_Gatherv(my_copy, (step == 0) ? 0 : step - xSize, MPI_DOUBLE, 
                arr + z*ySize*xSize, recv_range, displace, MPI_DOUBLE,
                ROOT, MPI_COMM_WORLD);    
    }

    if(rank == ROOT) {
        free(range);
        free(displace);
        free(recv_range);
    }
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
