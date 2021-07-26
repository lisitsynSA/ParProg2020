#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>

/* Константы для игрового поля */
#define DEAD 0
#define ALIVE 1

/* Флаги для включения и отключения функций */
#define ENABLED 1
#define DISABLED 0

/* Включение режима отладочной печати */
#define DEBUG_PRINT_INIT DISABLED

/* Макрос отладочной печати */
#define DEBUG_PRINT( sequence) do{ if( DEBUG_PRINT_INIT) (sequence);} while(0);

/* Распространненый цикл */
#define FOR_ALL_ARRAY\
 for ( uint32_t i = 0; i < xSize; i++)\
	for( uint32_t j = 0; j < ySize; j++)\

/* Заполнение массива 0 */
#define SET_ARRAY_ZERO( array)\
    FOR_ALL_ARRAY\
   	    ( array)[i][j] = 0;

/* Копирование одного массива в другой */
#define COPY_ARRAY( src, dst)\
    FOR_ALL_ARRAY\
	( dst)[i][j] = ( src)[i][j];

/**
 * Выполняет циклический сдвиг по координатам.
 *
 * Функция void, так как нам необходимо вернуть 2 значения.
 */
void CycleShift( uint32_t xSize_t, uint32_t ySize_t, int* x, int* y)
{
    int xSize = xSize_t;
    int ySize = ySize_t;

    /* Если x за пределами поля, тогда выполняем циклический сдвиг по x*/ 
    if ( !( ( *x <= xSize) && ( *x >= 0)))
    {
	if ( *x < 0)
	    *x = xSize;
	else
	    if ( *x >= xSize)
		*x = 0;
    }
    
    /* Если y за пределами поля, тогда выполняем циклический сдвиг по y */
    if( !(( *y <=  ySize) && ( *y >= 0)))
    {
	if ( *y < 0)
		*y = ySize;
	else
	    if ( *y >= ySize)
		*y = 0;

    }
    
    /* Иначе выполняем циклический сдвиг в зависимости от того, на сколько отличаются координаты от границы */

    }

/**
 * Конвертирует массив в 2D для более удобной работы
 *
 * WARNING Массив удаляется из вызывающей функции
 */
uint8_t** ConvertFrom1Dto2D( uint32_t xSize, uint32_t ySize, uint8_t* input_frame)
{
    uint8_t** temp_frame = (uint8_t**)calloc( xSize, sizeof(uint8_t*));
    for ( uint32_t i = 0; i < xSize; i++)
    {
        temp_frame[i] = (uint8_t*)calloc( ySize, sizeof(uint8_t));
    }
    
    /* x - столбцы, y - строки */
    for ( uint32_t i = 0; i < ySize; i++)
        for ( uint32_t j = 0; j < xSize; j++)
            temp_frame[i][j] = input_frame[i * ySize + j];
    return ( temp_frame);
}

/**
 * Функция для отладки. Рисует промежуточный результат
 */
void DrawField( uint32_t xSize, uint32_t ySize, uint8_t** field)
{
    printf("\n");
     for ( uint32_t i = 0; i < xSize; i++)
	{
	    for ( uint32_t j = 0; j < ySize; j++)
	    {
		printf("%d  ", field[i][j]);
	    }
	    printf("\n");
	}
     printf("\n");
}

/**
 * Переводит матрицу в линейный массив
 */
inline void ConvertFrom2Dto1D( uint32_t xSize, uint32_t ySize, uint8_t** input_frame, uint8_t* output_frame)
{
    
    for ( uint32_t i = 0; i < ySize; i++)
	for( uint32_t j = 0; j < xSize; j++)
	    output_frame[i * xSize + j] = input_frame[i][j];
}

/**
 * Вычисляет количество соседей возле указанной точки (x, y)
 */
int NumOfNeighbours( uint32_t xSize, uint32_t ySize, uint32_t x, uint32_t y, uint8_t** input_frame)
{
    int neighbours_num = 0;

    xSize--;
    ySize--;
    
    /* Необходимо вокруг посмотреть 8 клеток, если там, кто-то есть, то количество соседей увеличиваем на 1 */
    for ( int i = -1; i <= 1; i++)
    {
	for ( int j = -1; j <= 1; j++)
	{
	    int x_shifted = x + i, 
		y_shifted = y + j;

	    CycleShift( xSize, ySize, &x_shifted, &y_shifted);
	    if ( !( ( i == 0) && ( j == 0)))
	    {
		if ( ( input_frame[x_shifted][y_shifted] == ALIVE))
		{
		    neighbours_num++;
		}
	    }
	}
    }

 
    return ( neighbours_num);
}


void calc(uint32_t xSize, uint32_t ySize, uint32_t iterations, uint32_t num_threads, uint8_t* inputFrame, uint8_t* outputFrame)
{
    omp_set_dynamic( 0);
    omp_set_num_threads( num_threads);

    uint8_t** temp_frame = ConvertFrom1Dto2D( xSize, ySize, inputFrame);
   
    uint8_t** result_frame = (uint8_t**)calloc( xSize, sizeof(uint8_t*));
    for ( uint32_t i = 0; i < xSize; i++)
    {
        result_frame[i] = (uint8_t*)calloc( ySize, sizeof(uint8_t));
    }
    
    /* Для каждой ячейки необходимо проверить соседние и посчитать количество соседей */
    for ( uint32_t k = 0; k < iterations; k++)
	{	    
	    #pragma omp parallel for num_threads( num_threads)
	    for ( uint32_t i = 0; i < xSize; i++)
	    {
		#pragma omp parallel for num_threads( num_threads)
		for ( uint32_t j = 0; j < ySize; j++)
		{
		    int neighbour_num = NumOfNeighbours( xSize, ySize, i, j, temp_frame);
		    
    		    /* В пустой клетке, рядом с которой есть ровно 3 соседа, зарождается жизнь */
		    if ( (  neighbour_num == 3) && (temp_frame[i][j] == DEAD))
			result_frame[i][j] = ALIVE;
		    else
		    /* Если у живой клетки 2 или 3 соседа, то она продолжает жить */
		    if ( ( ( neighbour_num == 2) || ( neighbour_num == 3)) && (temp_frame[i][j] == ALIVE))
			result_frame[i][j] = ALIVE;
		    else
			result_frame[i][j] = DEAD;
		}
	    }

	    #pragma omp parallel for num_threads( num_threads)
	    COPY_ARRAY( result_frame, temp_frame);
	}


    ConvertFrom2Dto1D( xSize, ySize, temp_frame, outputFrame);

    #pragma omp parallel for num_threads( num_threads)
    for ( uint32_t i = 0; i < xSize; i++)
    {
	free( temp_frame[i]);
	free( result_frame[i]);
    }

    free( temp_frame);
    free( result_frame);
 
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
