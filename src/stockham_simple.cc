#include <complex>
#include <cmath>
#include <immintrin.h>
#include <xmmintrin.h>

static void fft0(size_t n, size_t s, bool eo, std::complex<double>* x, std::complex<double>* y)
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
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

__m256d mul_cplx(const __m256d vec1, const __m256d vec2)
{
   const __m256d aa = _mm256_unpacklo_pd(vec1, vec1);
   const __m256d bb = _mm256_unpackhi_pd(vec1, vec1);
   const __m256d yx = _mm256_shuffle_pd(vec2, vec2, 5);
   return _mm256_addsub_pd(_mm256_mul_pd(aa, vec2), _mm256_mul_pd(bb, yx));
}

// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)

static void fft0_avx(size_t n, size_t s, bool eo, std::complex<double>* x, std::complex<double>* y)
{
/*   const size_t m = n/2;
   const double theta0 = 2*M_PI/n;

   if (n == 2)
   {
      std::complex<double> *z = eo ? y : x;
      if (s == 1)
      {
	 const __m128d a = _mm_load_pd((double*)x);
      }
      else
      {
	 for (size_t q = 0; q < s; q+=2) {
	    const __m256d a = _mm256_load_pd((double*)x + q*2);
	    const __m256d b = _mm256_load_pd((double*)x + q*2 + 2*s);
	    _mm256_store_pd((double*)z + q*2,  _mm256_add_pd(a + b));
	    _mm256_store_pd((double*)z + q*2 + 2*s,  _mm256_sub_pd(a + b));
	 }
      }
   }
   else  if (n > 2) {
      for (size_t p = 0; p < m; p+= 2) {
	 const std::complex<double>& wp = std::complex<double>(cos(p*theta0), -sin(p*theta0));
	 for (size_t q = 0; q < s; q++) {
	    const std::complex<double>& a = x[q + s*(p + 0)];
	    const std::complex<double>& b = x[q + s*(p + m)];
	    y[q + s*(2*p + 0)] =  a + b;
	    y[q + s*(2*p + 1)] = (a - b) * wp;
	 }
      }
      fft0_avx(n/2, 2*s, !eo, y, x);
      }*/
}

void stockham_fft_avx(std::complex<double>* x, size_t n) // Fourier transform
// n : sequence length
// x : input/output sequence
{
    std::complex<double>* y = new std::complex<double>[n];
    fft0_avx(n, 1, 0, x, y);
    delete[] y;
    for (size_t k = 0; k < n; k++)
       x[k] /= n;
}
