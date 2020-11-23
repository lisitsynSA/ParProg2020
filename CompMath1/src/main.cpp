#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

double acceleration(double t)
{
  return sin(t);
}

void calc(double* trace, uint32_t traceSize, double t0, double dt, double y0, double y1, int rank, int size)
{
	MPI_Status status;
	MPI_Bcast (&traceSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast (&t0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast (&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	uint32_t mySize = traceSize / size;
	if (rank == 0) {
		mySize += traceSize % size;
	}

	double tau = dt * mySize;

	if (rank != 0) {
		trace = (double *)malloc(mySize * sizeof(double));
	}
	// Sighting shot
	double v0 = 0, v1, myY1;
	t0 += tau * rank;
	if (rank != 0) {
		t0 += dt * (traceSize % size);
	}
	trace[0] = y0;
	trace[1] = y0 + dt*v0;
	for (uint32_t i = 2; i < mySize; i++) {
		trace[i] = dt*dt*acceleration(t0 + (i - 1)*dt) + 2*trace[i - 1] - trace[i - 2];
	}

	v1 = (trace[mySize - 1] - trace[mySize - 2]) / dt;
	myY1 = trace[mySize - 1];

	if (size == 1) {
		v0 = (y1 - myY1) / tau;
	} else {
		double bufY, bufV;
		if (rank == 0) {
			MPI_Send(&myY1, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
			MPI_Send(&v1, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
			MPI_Recv(&bufY, 1, MPI_DOUBLE, size - 1, 0, MPI_COMM_WORLD, &status);	
			v0 = (y1 - bufY) / (traceSize * dt);
			MPI_Bcast(&v0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		} else {
			MPI_Recv(&bufY, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&bufV, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);

			myY1 += bufY + bufV * (mySize - 1) * dt;
			v1 += bufV;

			if (rank != size - 1) {
				MPI_Send(&myY1, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
				MPI_Send(&v1, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);	
			} else {
				MPI_Send(&myY1, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			}
			
			MPI_Bcast(&v0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			y0 = bufY + v0 * tau * rank + v0 * dt * (traceSize % size);
			v0 += bufV;
		}

}		
  // The final shot
	trace[0] = y0;
	trace[1] = y0 + dt*v0;
	for (uint32_t i = 2; i < mySize; i++) {
		 trace[i] = dt*dt*acceleration(t0 + (i - 1)*dt) + 2*trace[i - 1] - trace[i - 2];
	}
	if (rank == 0) {
		mySize -= (traceSize % size);
		for(int i = 1; i < size; i++) {
			MPI_Recv(&trace[i * mySize + (traceSize % size)], mySize, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
		}
	} else {
		MPI_Send(trace, mySize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		free(trace);
	}
}


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
      output << " " << trace[i];
    }
    output << std::endl;
    output.close();
    delete trace;
  }

  MPI_Finalize();
  return 0;
}
