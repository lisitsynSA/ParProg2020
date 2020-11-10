#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>

void calc(uint32_t xSize, uint32_t ySize, uint32_t iterations, uint32_t num_threads, uint8_t* inputFrame, uint8_t* outputFrame)
{
  int* inFrame = (int*)malloc(xSize*ySize * sizeof(int));
  int* outFrame = (int *)malloc(xSize*ySize * sizeof(int));

  for (uint32_t y = 0; y < ySize; y++) {
    for (uint32_t x = 0; x < xSize; x++) {
      outFrame[y * xSize + x] = (int)inputFrame[y * xSize + x];
    }
  }
  uint32_t block_size = xSize / num_threads;
  uint32_t block_extra = xSize % num_threads;
  #pragma omp parallel num_threads(num_threads)
  {
    uint32_t tid = omp_get_thread_num();

    uint32_t first = tid * block_size;
    uint32_t last = (tid + 1) * block_size;

    if (tid == num_threads - 1) {
      last += block_extra;
    }


    for (uint32_t i = 0; i < iterations; i++) {
      #pragma omp barrier
      #pragma omp single 
      {
        int* ptr = inFrame;
        inFrame = outFrame;
        outFrame = ptr;
      }
      for (uint32_t x = first; x < last; x++) {
          for (uint32_t y = 0; y < ySize; y++) {
            uint32_t count = 0;

            uint32_t x1 = x == 0 ? xSize - 1 : x - 1;
            uint32_t x2 = x == xSize - 1 ? 0 : x + 1;
            uint32_t y1 = y == 0 ? ySize - 1 : y - 1;
            uint32_t y2 = y == ySize - 1 ? 0 : y + 1;


            if (inFrame[y1 * xSize + x] == 1){
              count++;
            }
            if (inFrame[y1 * xSize + x1] == 1) {
              count++;
            }
            if (inFrame[y1 * xSize + x2] == 1) {
              count++;
            }
            if (inFrame[y2 * xSize + x1] == 1) {
              count++;
            }
            if (inFrame[y2 * xSize + x2] == 1) {
              count++;
            }

            if (inFrame[y2 * xSize + x] == 1) {
              count++;
            }
            if (inFrame[y * xSize + x1] == 1) {
              count++;
            }
            if (inFrame[y * xSize + x2] == 1) {
              count++;
            }

            if (inFrame[y * xSize + x] == 0 && count == 3) {
              outFrame[y * xSize + x] = 1;
            } else if (inFrame[y * xSize + x] == 1 && (count < 2 || count > 3)) {
              outFrame[y * xSize + x] = 0;
            } else { 
              outFrame[y * xSize + x] = inFrame[y * xSize + x];
            }
        }
      }
    }
  }

  for (uint32_t y = 0; y < ySize; y++) {
    for (uint32_t x = 0; x < xSize; x++) {
      outputFrame[y * xSize + x] = outFrame[y * xSize + x];
    }
  }
  free(inFrame);
  free(outFrame);
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
  delete[] outputFrame;
  delete[] inputFrame;
  output.close();
  input.close();
  return 0;
}
