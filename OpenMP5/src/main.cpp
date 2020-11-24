//:set tabstop=2

#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>

void calc(uint32_t xSize, uint32_t ySize, uint32_t iterations, uint32_t num_threads, uint8_t* inputFrame, uint8_t* outputFrame)
{
  uint8_t* neighbours_map[xSize + 2][ySize + 2];  //to easy understand who is neighbour
	#pragma omp parallel for num_threads(num_threads) shared(neighbours_map, inputFrame)
	for(uint64_t i = 0; i < xSize * ySize; ++i)
	{
		uint32_t tmp_x = i % xSize;
		uint32_t tmp_y = i / xSize;
		neighbours_map[tmp_x + 1][tmp_y + 1] = &inputFrame[i];
	}
	#pragma omp parallel for num_threads(num_threads) shared(neighbours_map, inputFrame)
	for(uint32_t i = 0; i < xSize; ++i)
	{
		neighbours_map[i + 1][0] = &inputFrame[i + (ySize-1)*xSize];
		neighbours_map[i + 1][ySize + 1] = &inputFrame[i + 0*xSize];
	}
	#pragma omp parallel for num_threads(num_threads) shared(neighbours_map, inputFrame)
	for(uint32_t i = 0; i < ySize; ++i)
	{
		neighbours_map[0][i + 1] = &inputFrame[(xSize-1) + i*xSize];
		neighbours_map[xSize + 1][i + 1] = &inputFrame[0 + i*xSize];
	}
	neighbours_map[0][0] = &inputFrame[(xSize-1)+(ySize-1)*xSize];
	neighbours_map[0][ySize + 1] = &inputFrame[(xSize-1) + 0*xSize];
	neighbours_map[xSize + 1][0] = &inputFrame[0 + (ySize-1)*xSize];
	neighbours_map[xSize + 1][ySize + 1] = &inputFrame[0 + 0*xSize];
	
	if(iterations == 0) // little kostil
		#pragma omp parallel for num_threads(num_threads) shared(inputFrame, outputFrame)
		for(uint64_t i = 0; i < xSize*ySize; ++i)
			outputFrame[i] = inputFrame[i];

	for(uint32_t iter = 0; iter < iterations; ++iter)
	{
		#pragma omp parallel for num_threads(num_threads) shared(neighbours_map, inputFrame, outputFrame)
		for(uint64_t i = 0; i < xSize*ySize; ++i)
		{
			uint32_t x = i % xSize;
			uint32_t y = i / xSize;
			uint8_t neighbours_sum = *neighbours_map[x][y] + *neighbours_map[x+1][y] +            //123 (1,2)
				*neighbours_map[x+2][y] + *neighbours_map[x][y+1] + *neighbours_map[x+2][y+1] +     //4X5 (3,4,5)
				*neighbours_map[x][y+2] + *neighbours_map[x+1][y+2] + *neighbours_map[x+2][y+2];    //678 (6,7,8)
			if((inputFrame[i] == 0 && neighbours_sum == 3) || \
				(inputFrame[i] == 1 && (neighbours_sum == 2 || neighbours_sum == 3)))
				outputFrame[i] = 1;
			else
				outputFrame[i] = 0;
		}
		#pragma omp parallel for num_threads(num_threads) shared(inputFrame, outputFrame)
		for(uint64_t i = 0; i < xSize*ySize; ++i)
			inputFrame[i] = outputFrame[i];
	}
}

int main(int argc, char** argv)
{
  // Check arguments
  if (argc != 3)
  {
    std::cout << "[Error] Usage <inputfile> <output file>\n";
    return 1;
  }

  // Prepare input file
  std::ifstream input(argv[1]);
  if (!input.is_open())
  {
    std::cout << "[Error] Can't open " << argv[1] << " for write\n";
    return 1;
  }

  // Prepare output file
  std::ofstream output(argv[2]);
  if (!output.is_open())
  {
    std::cout << "[Error] Can't open " << argv[2] << " for read\n";
    input.close();
    return 1;
  }

  // Read arguments from input
  uint32_t xSize = 0, ySize = 0, iterations = 0, num_threads = 0;
  input >> xSize >> ySize >> iterations >> num_threads;
  uint8_t* inputFrame = new uint8_t[xSize*ySize];
  uint8_t* outputFrame = new uint8_t[xSize*ySize];
  for (uint32_t y = 0; y < ySize; y++)
  {
    for (uint32_t x = 0; x < xSize; x++)
    {
      input >> inputFrame[y*xSize + x];
      inputFrame[y*xSize + x] -= '0';
    }
  }

  // Calculation
  calc(xSize, ySize, iterations, num_threads, inputFrame, outputFrame);

  // Write result
  for (uint32_t y = 0; y < ySize; y++)
  {
    for (uint32_t x = 0; x < xSize; x++)
    {
      outputFrame[y*xSize + x] += '0';
      output << " " << outputFrame[y*xSize + x];
    }
    output << "\n";
  }

  // Prepare to exit
  delete outputFrame;
  delete inputFrame;
  output.close();
  input.close();
  return 0;
}
