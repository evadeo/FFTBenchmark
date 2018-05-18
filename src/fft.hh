#include <complex>
#include <valarray>

void difft2_parallel(std::complex<double>* arr, size_t n);
void difft2(std::complex<double>* arr, size_t n);
std::valarray<std::complex<double>> iterative_fft(std::valarray<std::complex<double> > input);
