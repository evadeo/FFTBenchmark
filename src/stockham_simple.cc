#include <complex>
#include <cmath>
#include <immintrin.h>
#include <xmmintrin.h>

static void fft0(size_t n, size_t s, bool eo, std::complex<double>* x,
		 std::complex<double>* y)
{
    const size_t m = n/2;
    const double theta0 = 2*M_PI/n;

    if (n == 1)
    {
       if (eo)
	  for (size_t q = 0; q < s; q++)
	     y[q] = x[q];
    }
    else {
       for (size_t p = 0; p < m; p++) {
	  const std::complex<double>& wp = std::complex<double>(cos(p*theta0), -sin(p*theta0));
	  for (size_t q = 0; q < s; q++) {
	     const std::complex<double>& a = x[q + s*(p + 0)];
	     const std::complex<double>& b = x[q + s*(p + m)];
	     y[q + s*(2*p + 0)] =  a + b;
	     y[q + s*(2*p + 1)] = (a - b) * wp;
	  }
       }
       fft0(n/2, 2*s, !eo, y, x);
    }
}

void stockham_fft(std::complex<double>* x, size_t n) // Fourier transform
// n : sequence length
// x : input/output sequence
{
    std::complex<double>* y = new std::complex<double>[n];
    fft0(n, 1, 0, x, y);
    delete[] y;
    for (size_t k = 0; k < n; k++)
       x[k] /= n;
}
