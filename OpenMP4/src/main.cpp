#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <math.h>

#include <vector>
/*
typedef struct {
    double value;
    size_t gr;
} limd_t;

limd_t limd_add(limd_t lv, limd_t rv) {
    size_t new_gr = 0;
    if(lv.gr > rv.gr) {
        new_gr = lv.gr = rv.gr;
        lv.value /= pow(10, lv.gr - rv.gr);
    } else {
        new_gr = rv.gr = lv.gr;
        lv.value /= pow(10, rv.gr - lv.gr);
    }
    limd_t res;
    res.value = lv.value + rv.value;
    res.gr = new_gr;
    return res;
}

limd_t limd_mul(limd_t lv, limd_t rv) {
    limd_t res;
    res.value = lv.value * rv.value;
    res.gr = lv.gr + rv.gr;
    return res;
}

limd_t limd_sub(limd_t lv, limd_t rv) {
    rv.value *= -1;
    return limd_add(lv, rv);
}

limd_t limd_dev(limd_t lv, limd_t rv) {
    rv.value = 1 / rv.value;
    rv.gr *= -1;
    return limd_mul(lv, rv);
}

limd_t limd_init(double value) {
    if(value == 0) {
        return {0, (size_t(-1))};
    }
    limd_t res = {value, 0};
    while(res.value > 10 && res.value < 1) {
        if(res.value > 1.0) {
            res.value /= 10;
            res.gr--;
        } else {
            res.value *= 10;
            res.gr++;
        }
    }
    return res;
}*/

double calc(uint32_t x_last, uint32_t num_threads)
{
    double exp = 0.0;
    std::vector<double> fact(x_last);
    for(auto i = next(fact.begin()); i != fact.end(); i = next(i)) {
        *i = -1;
    }
    fact[0] = 1;
    double c = 0.0;
    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp for
        for(size_t n = 1; n <= x_last; ++n) {
            double my_fact = ( 1.0 / n );
            for(size_t i = n - 1; i >= 0; --i) {
                if(fact[i] == -1) {
                    my_fact *= ( 1.0 / i );
                } else if(fact[i] > 0) {
                    my_fact *= fact[i];
                    fact[n] = my_fact;
                    break;
                } else { // fact[i] < 0
                    break;
                }
            }
            //std::cout << "I HAVE " << my_fact << std::endl;
            /*
            #pragma omp critical 
            {
                double y = my_fact - c;
                double t = exp + y;
                c = (t - exp) - y;
                exp = t;
            }*/
        }
    }
    
    for(auto i = fact.rbegin(); i != fact.rend(); i = next(i)) {
        exp += *i;
    }
    return exp;
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
  uint32_t x_last = 0, num_threads = 0;
  input >> x_last >> num_threads;

  // Calculation
  double res = calc(x_last, num_threads);

  // Write result
  output << std::setprecision(16) << res << std::endl;
  // Prepare to exit
  output.close();
  input.close();
  return 0;
}
