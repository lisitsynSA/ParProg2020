#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>
#include <cassert>

double acceleration(double t)
{
  return sin(t);
}

#ifdef REFERENCE
void calc(double* trace, uint32_t traceSize, double t0, double dt, double y0, double y1, int rank, int size)
#else
void st_calc(double* trace, uint32_t traceSize, double t0, double dt, double y0, double y1, int rank, int size)
#endif
{
  // Sighting shot
  double v0 = 0;
  if (rank == 0 && size > 0)
  {
    trace[0] = y0;
    trace[1] = y0 + dt*v0;
    for (uint32_t i = 2; i < traceSize; i++)
    {
      trace[i] = dt*dt*acceleration(t0 + (i - 1)*dt) + 2*trace[i - 1] - trace[i - 2];
    }
  }

  // The final shot
  if (rank == 0 && size > 0)
  {
    v0 = (y1 - trace[traceSize - 1])/(dt*traceSize);
    trace[0] = y0;
    trace[1] = y0 + dt*v0;
    for (uint32_t i = 2; i < traceSize; i++)
    {
      trace[i] = dt*dt*acceleration(t0 + (i - 1)*dt) + 2*trace[i - 1] - trace[i - 2];
    }
  }
}

#ifndef REFERENCE
void calc(double* trace, uint32_t trace_size, double t0, double dt, double y0, double y1, int rank, int size)
{
  if (size == 0)
    return;
  if (size == 1) {
    st_calc(trace, trace_size, t0, dt, y0, y1, rank, size);
    return;
  }
  // In this implementation cannot parallel further than trace_size
  MPI_Comm my_comm;
  int color = (rank >= (int)trace_size) ? MPI_UNDEFINED : 0;
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &my_comm);
  if (rank >= (int)trace_size)
    return;
  if (size > (int)trace_size)
    size = trace_size;

  MPI_Bcast(&trace_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&t0,         1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dt,         1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  uint32_t trace_per_process = trace_size / size;
  double* my_shot = (double*)calloc(trace_per_process, sizeof(*my_shot));
  assert(my_shot);

  double tau = dt * trace_size / size;

  // Sighting shot

  double _y0 = rank ? 0 : y0;
  double _v0 = 0;
  my_shot[0] = _y0; my_shot[1] = _y0 + dt * _v0;
  double my_t0 = t0 + tau * rank;
  for (uint32_t i = 2; i < trace_per_process; i++)
    my_shot[i] = 2 * my_shot[i - 1] - my_shot[i - 2] + dt * dt * acceleration(my_t0 + (i - 1) * dt);
  double _y1 = my_shot[trace_per_process - 1];
  double _v1 = (my_shot[trace_per_process - 1] - my_shot[trace_per_process - 2]) / dt;

  double py = 0, pv = 0;
  if (rank != 0) {
    MPI_Recv(&py, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, NULL);
    MPI_Recv(&pv, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, NULL);
    _y0 = py; _v0 = pv;
    _y1 += py + pv * tau; _v1 += pv;
  }

  if (rank != size - 1) {
    MPI_Send(&_y1, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    MPI_Send(&_v1, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
  } else
    MPI_Send(&_y1, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

  double v0 = 0;
  if (rank == 0)
  {
    MPI_Recv(&py, 1, MPI_DOUBLE, size - 1, 0, MPI_COMM_WORLD, NULL);
    v0 = (y1 - py) / (dt * trace_size);
  }
  MPI_Bcast(&v0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  _y0 += v0 * tau * rank; _v0 += v0;

  // The final shot

  my_shot[0] = _y0; my_shot[1] = _y0 + dt * _v0;
  for (uint32_t i = 2; i < trace_per_process; i++)
    my_shot[i] = 2 * my_shot[i - 1] - my_shot[i - 2] + dt * dt * acceleration(my_t0 + (i - 1) * dt);

  if (rank)
    MPI_Send(my_shot, trace_per_process, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  else {
    memcpy(trace, my_shot, trace_per_process * sizeof(double));
    for (int i = 1; i < size; i++)
      MPI_Recv(trace + trace_per_process * i, trace_per_process, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, NULL);
  }

  free(my_shot);
}
#endif /* ~REFERENCE */

int main(int argc, char** argv)
{
  int rank = 0, size = 0, status = 0;
  uint32_t traceSize = 0;
  double t0 = 0, t1 = 0, dt = 0, y0 = 0, y1 = 0;
  double* trace = 0;

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
    input >> t0 >> t1 >> dt >> y0 >> y1;
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    traceSize = (t1 - t0)/dt;
    trace = new double[traceSize];

    input.close();
  } else {
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (status != 0)
    {
      return 1;
    }
  }

  calc(trace, traceSize, t0, dt, y0, y1, rank, size);

  if (rank == 0)
  {
    // Prepare output file
    std::ofstream output(argv[2]);
    if (!output.is_open())
    {
      std::cout << "[Error] Can't open " << argv[2] << " for read\n";
      delete trace;
      return 1;
    }

    for (uint32_t i = 0; i < traceSize; i++)
    {
      output << std::fixed << std::setprecision(13) << trace[i] << std::endl;
    }
    output.close();
    delete trace;
  }

  MPI_Finalize();
  return 0;
}
