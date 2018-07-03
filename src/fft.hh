#include <complex>
#include <valarray>

/* Recursive part */
void difft2_parallel(std::complex<double>* arr, size_t n);
void difft2(std::complex<double>* arr, size_t n);

/* Iterative part */
void stockham_openmp_fft(std::complex<double>* x, size_t n);
std::valarray<std::complex<double>> iterative_fft(std::valarray<std::complex<double>> input);
std::valarray<std::complex<double>> iterative_fft_parallel(std::valarray<std::complex<double>> input);
std::valarray<std::complex<double>> iterative_openmp_fft(std::valarray<std::complex<double>> input);
int my_log2(int N);
unsigned int reverse(int N, int n);
