#include <fftw3.h>
#include <complex>

void fftwBench(std::complex<double>* arr, size_t n)
{
   fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof (fftw_complex) * n);
   fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof (fftw_complex) * n);

   for (size_t i = 0; i < n; ++i)
   {
      in[i][0] = arr[i].real();
      in[i][1] = arr[i].imag();
   }

   fftw_plan p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
   fftw_execute(p);

   fftw_destroy_plan(p);
   fftw_free(in);
   fftw_free(out);
}


void fftwBenchParr(std::complex<double>* arr, size_t n)
{
   fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof (fftw_complex) * n);
   fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof (fftw_complex) * n);

   for (size_t i = 0; i < n; ++i)
   {
      in[i][0] = arr[i].real();
      in[i][1] = arr[i].imag();
   }

   fftw_plan p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
   fftw_execute(p);

   fftw_destroy_plan(p);
   fftw_free(in);
   fftw_free(out);
}
