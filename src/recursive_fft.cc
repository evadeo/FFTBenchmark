#include <complex>
#include <iostream>
#include "timer.hh"

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
