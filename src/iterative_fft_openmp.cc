#include "timer.hh"
#include <complex>
#include <valarray>
#include <iostream>
#include <cmath>
#include <omp.h>
#include "fft.hh"

std::valarray<std::complex<double> > bit_reverse_copy_para(std::valarray<std::complex<double> > input, unsigned int length) {
  std::valarray<std::complex<double> > res (length);
  unsigned int i;
  #pragma omp parallel shared ( input, length) private (i)
  #pragma omp for
  for (i = 0; i < length; ++i) {
    res[reverse(length, i)] = input[i];
  }

  return res;
}

std::valarray<std::complex<double> > iterative_openmp_fft(std::valarray<std::complex<double> > input) {
  unsigned int length = input.size();
  std::valarray<std::complex<double> > rev_array = bit_reverse_copy_para(input, length);

  std::valarray<std::complex<double> > omega (length / 2);

  omega[1] = std::polar(1., -2. * M_PI / length);
  omega[0] = 1;
  unsigned int i;
  //#pragma omp parallel shared (length, omega) private (i)
  //#pragma omp for
  for (i = 2; i < length / 2; ++i) {
    omega[i] = std::pow(omega[1], i);
  }
  int n = 1;
  int a = length / 2;
  unsigned int k;
  for (int j = 0; j < my_log2(length); ++j) {
    #pragma omp parallel shared ( input, length, omega, n, a ) private (k)
    #pragma omp for
    for (k = 0; k < length; ++k) {
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

