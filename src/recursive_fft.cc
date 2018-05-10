#include <complex>
#include <iostream>
#include "timer.hh"

#define M_PI 3.14159265358979323846


// Separate even & odd elements
void separate (std::complex<double>* a, size_t n) {
   std::complex<double>* b = new std::complex<double>[n/2];
   // FIXME: Merge the two loops
   for (size_t i = 0; i < n/2; ++i) {
      b[i] = a[i*2 + 1];
      a[i] = a[i*2];
   }

   std::move(b, b + n/2, a + n/2); //FIXME: Make it parallel
   delete[] b;
}

void difft2(std::complex<double>* arr, size_t n)
{
   if (n == 1)
      return;
   separate(arr, n);
   difft2(arr, n/2);
   difft2(arr + n/2, n/2);

   for (size_t k = 0; k < n/2; ++k)
   {
      std::complex<double> t = arr[k];
      auto x = std::exp(std::complex<double>(0, -2. * M_PI * k / n));
      arr[k] = t + x * arr[k + n/2];
      arr[k + n/2] = t - x * arr[k + n/2];
   }
}


int main() {
   const unsigned int nSamples = std::pow(2,23);

   std::complex<double> *x = new std::complex<double>[nSamples];
   std::complex<double> *X = new std::complex<double>[nSamples];
   if (!x || !X)
   {
      std::cout << "Malloc error, exiting" << std::endl;
      return 1;
   }
   const int nFreqs = 5;
   double freq[nFreqs] = { 2, 5, 11, 17, 29 }; // known freqs for testing

   // generate samples for testing
   for(size_t i=0; i<nSamples; i++) {
      x[i] = std::complex<double>(0.,0.);
      // sum several known sinusoids into x[]
      for(int j=0; j<nFreqs; j++)
	 x[i] += sin( 2*M_PI*freq[j]*i/nSamples );
      X[i] = x[i];        // copy into X[] for FFT work & result
   }
   // compute fft for this data
   double timer;
   std::cout << "fft" << std::endl;

   {
      auto timerC = scope_timer(timer);
      difft2(X, nSamples);
   }

   std::cout << "----------------------" << std::endl;
   std::cout << "Execution time: " << timer << std::endl;

   /*for (int i = 0; i < 8; ++i)
     std::cout << samples[i] << " ";*/
   std::cout << std::endl << "----------------------" << std::endl;

   delete[] x;
   delete[] X;
   return 0;
}
