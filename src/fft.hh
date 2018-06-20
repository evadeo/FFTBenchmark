#include <complex>
#include <valarray>

/* Recursive part */
void difft2_parallel(std::complex<double>* arr, size_t n);
void difft2(std::complex<double>* arr, size_t n);

/* Iterative part */
std::valarray<std::complex<double>> iterative_fft(std::valarray<std::complex<double> > input);
std::valarray<std::complex<double>> iterative_fft_parallel(std::valarray<std::complex<double> > input);

int my_log2(int N);
unsigned int reverse(int N, int n);
void stockham_fft(std::complex<double>* arr, size_t n);
void difft2_avx(std::complex<double>* arr, size_t n);
void fftwBench(std::complex<double>* arr, size_t n);
void fftwBenchParr(std::complex<double>* arr, size_t n);
void my_cufft(std::complex<double>* arr, size_t n);
void do_square(double *a_h, const int N);
