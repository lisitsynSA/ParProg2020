#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

#define ROOT 0



#define MPI_recv_start_data(worker, have_task) { \
    MPI_Recv(worker, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, NULL); \
    MPI_Recv(have_task, 1, MPI_INT, *worker, 1, MPI_COMM_WORLD, NULL); \
}

#define MPI_send_start_data(rank, have_task) { \
    MPI_Send(rank, 1, MPI_INT, ROOT, 1, MPI_COMM_WORLD); \
    MPI_Send(have_task, 1, MPI_INT, ROOT, 1, MPI_COMM_WORLD); \
}

#define MPI_recv_res(worker, arr, xSize) { \
    uint32_t ux, uy; \
    MPI_Recv(&ux, 1, MPI_UNSIGNED, worker, 2, MPI_COMM_WORLD, NULL); \
    MPI_Recv(&uy, 1, MPI_UNSIGNED, worker, 2, MPI_COMM_WORLD, NULL); \
    MPI_Recv(&(arr[uy*xSize + ux]), 1, MPI_DOUBLE, worker, 2, MPI_COMM_WORLD, NULL); \
}  

#define MPI_send_res(x, y, arr) { \
    MPI_Send(x, 1, MPI_UNSIGNED, ROOT, 2, MPI_COMM_WORLD); \
    MPI_Send(y, 1, MPI_UNSIGNED, ROOT, 2, MPI_COMM_WORLD); \
    MPI_Send(arr, 1, MPI_DOUBLE, ROOT, 2, MPI_COMM_WORLD); \
}

#define MPI_recv_data(x, y, arr) { \
    MPI_Recv(x, 1, MPI_UNSIGNED, ROOT, 1, MPI_COMM_WORLD, NULL); \
    MPI_Recv(y, 1, MPI_UNSIGNED, ROOT, 1, MPI_COMM_WORLD, NULL); \
    if(*x == uint32_t(-1) && *y == uint32_t(-1)) { \
        break; \
    } \
    MPI_Recv(arr, 1, MPI_DOUBLE, ROOT, 1, MPI_COMM_WORLD, NULL); \
}

#define MPI_send_data(worker, x, y, arr, xSize, need_to_send_arr) { \
    MPI_Send(x, 1, MPI_UNSIGNED, worker, 1, MPI_COMM_WORLD); \
    MPI_Send(y, 1, MPI_UNSIGNED, worker, 1, MPI_COMM_WORLD); \
    if(need_to_send_arr) { \
        MPI_Send(&arr[*y * xSize + *x], 1, MPI_DOUBLE, worker, 1, MPI_COMM_WORLD); \
    } \
}

void calc(double* arr, uint32_t ySize, uint32_t xSize, int rank, int size)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(size == 1) {
        for(uint32_t y = 0; y < ySize; ++y) {
            for(uint32_t x = 0; x < xSize; ++x) {
                arr[y*xSize + x] = sin(0.00001*arr[y*xSize + x]);
            }
        }
        return; 
    }
 
    //we will
    if(rank == ROOT) {
        int worker = 0;
        int have_task = 0;
        for(uint32_t y = 0; y < ySize; ++y) {
            for(uint32_t x = 0; x < xSize; ++x) {
                MPI_recv_start_data(&worker, &have_task);  
                
                if(have_task != 0) { // Recv task
                    MPI_recv_res(worker, arr, xSize);
                }
                
                MPI_send_data(worker, &x, &y, arr, xSize, 1);
            }
        }
        for(int i = 1; i < size; ++i) {
            MPI_recv_start_data(&worker, &have_task);  
            
            if(have_task != 0) { // Recv task
                MPI_recv_res(worker, arr, xSize);
            }
            
            uint32_t x = uint32_t(-1), y = x; 
            MPI_send_data(worker, &x, &y, arr, xSize, 0);
        }
    } else {
        int have_task = 0;
        while(1) {
            MPI_send_start_data(&rank, &have_task);
            have_task = 1;
             
            uint32_t x = 0, y = 0;
            double arr = 0;
            MPI_recv_data(&x, &y, &arr);
            
            arr = sin(0.00001*arr);
            MPI_send_res(&x, &y, &arr);
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
