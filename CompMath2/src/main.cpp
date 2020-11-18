#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

#define ROOT 0

void calc(double* frame, uint32_t ySize, uint32_t xSize, double delta, int rank, int size) { 
    int* send_dist = NULL;
    int* send_range = NULL;
    
    int* recv_dist = NULL;
    int* recv_range = NULL;

    double* new_frame = NULL;

    MPI_Bcast(&ySize, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&xSize, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD);
    int real_rank_size = 0;

    if(rank == ROOT) {
        size_t my_calc_size = ceil((ySize) * 1.0 / size);
        // if  proc_rank > ceil(ySize / my_calc_size) it will do nothing
        real_rank_size = ceil((ySize) * 1.0 / my_calc_size);
        size_t now_size = 0; 
        
        new_frame = (double*)calloc(xSize*ySize, sizeof(double));

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
            if(real_rank_size == 1) {
                send_range[i] = recv_range[i];
            } else if(i == 0) {
                send_range[i] = recv_range[i] + xSize;
            } else if(i == real_rank_size - 1) {
                send_range[i] = recv_range[i] + xSize; 
                send_dist[i] = recv_dist[i] - xSize; 
            } else {
                send_range[i] = recv_range[i] + 2*xSize; 
                send_dist[i] = recv_dist[i] - xSize;
            }
        }
        recv_dist[0] += xSize;
        recv_range[0] -= xSize;
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

    double* my_calc = (double*)calloc(send_step, sizeof(double));
    double* recv_calc = (double*)calloc(send_step, sizeof(double));
    double diff = 0;
    int stop = 0;
    do { 
        MPI_Scatterv(
                frame, send_range, send_dist, MPI_DOUBLE,
                my_calc, send_step, MPI_DOUBLE,
                ROOT, MPI_COMM_WORLD);
        for(size_t y = 1; y < send_step * 1.0 / xSize - 1; ++y) {
            for(size_t x = 1; x < xSize - 1; ++x) {
                recv_calc[y * xSize + x] = (my_calc[(y + 1) * xSize + x] +\
                        my_calc[(y - 1) * xSize + x] +\
                        my_calc[y * xSize + x + 1] +\
                        my_calc[y * xSize + x - 1]) / 4.0; 
            }
        }

        MPI_Gatherv(
                recv_calc + xSize, recv_step, MPI_DOUBLE,
                new_frame, recv_range, recv_dist, MPI_DOUBLE, 
                ROOT, MPI_COMM_WORLD);

        if(rank == ROOT) {
            for (uint32_t y = 1; y < ySize - 1; y++) {
                for (uint32_t x = 1; x < xSize - 1; x++) {
                    diff += std::abs(new_frame[y*xSize + x] - frame[y*xSize + x]);
                }
            }
            if(diff <= delta) {
                stop = 1;
            } else {
                stop = 0;
            }
            diff = 0;
            for (uint32_t y = 1; y < ySize - 1; y++) {
                for (uint32_t x = 1; x < xSize - 1; x++) {
                    frame[y*xSize + x] = new_frame[y*xSize + x];
                }
            } 
        }
        MPI_Bcast(&stop, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
        //printf("[%d] diff = %lf\n", rank, diff);
    } while(!stop);

    free(my_calc);
    free(recv_calc);
    if(rank == ROOT) {
        free(recv_dist);
        free(recv_range);
        free(send_dist);
        free(send_range);
    }
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
