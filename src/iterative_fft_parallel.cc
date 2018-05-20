#include <complex>
#include <valarray>
#include <tbb/tbb.h>
#include "fft.hh"

std::valarray<std::complex<double> > bit_reverse_copy_parallel(std::valarray<std::complex<double> > input, unsigned int length) {
  std::valarray<std::complex<double> > res (length);
  tbb::parallel_for(tbb::blocked_range<unsigned int>(0, length),
                    [&](const tbb::blocked_range<unsigned int>& r) {
                      for (auto i = r.begin(); i < r.end(); ++i)
                        res[reverse(length, i)] = input[i];
                    });

  return res;
}

std::valarray<std::complex<double> > iterative_fft_parallel(std::valarray<std::complex<double> > input) {
  unsigned int length = input.size();
  std::valarray<std::complex<double> > rev_array = bit_reverse_copy_parallel(input, length);

  std::valarray<std::complex<double> > omega (length / 2);

  omega[1] = std::polar(1., -2. * M_PI / length);
  omega[0] = 1;
  tbb::parallel_for(tbb::blocked_range<unsigned int>(2, length / 2),
                    [&](const tbb::blocked_range<unsigned int>& r) {
                      for (unsigned int i = r.begin(); i < r.end(); ++i) {
                        omega[i] = std::pow(omega[1] , i);
                      }
                    });
  int n = 1;
  int a = length / 2;

  for (int j = 0; j < my_log2(length); ++j) {
    for (unsigned int k = 0; k < length; ++k) {
      if (!(k & n)) {
        std::complex<double> firstTmp = rev_array[k];
        std::complex<double> secondTmp = omega[(k * a) % (n * a)] * rev_array[k + n];

        rev_array[k] = firstTmp + secondTmp;
        rev_array[k + n] = firstTmp - secondTmp;
      }
    }
    n *= 2;
    a /= 2;
  }

  return rev_array;
}
