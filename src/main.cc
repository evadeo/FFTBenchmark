#include "fft.hh"
#include "timer.hh"
#include <iostream>
#include <complex>
#include <valarray>
#include <algorithm>

std::complex<double> *generate_data(unsigned int nSamples)
{
   std::complex<double> *x = new std::complex<double>[nSamples];
   std::complex<double> *X = new std::complex<double>[nSamples];
   if (!x || !X)
   {
      std::cout << "Malloc error, exiting" << std::endl;
      return nullptr;
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
   delete[] x;
   return X;
}

void print_res(std::string name, double time)
{
   std::cout << name << std::endl;
   std::cout << "----------------------" << std::endl;
   std::cout << "Execution time: " << time << std::endl;
   std::cout << std::endl << "----------------------" << std::endl;

}

bool complex_close(std::complex<double> a, std::complex<double> b)
{
   auto c = a - b;
   return ((std::abs(c.real()) < 0.1) && (std::abs(c.imag()) < 0.1));
}

int main () {
   const int power = 19;
   const unsigned int nSamples = std::pow(2, power);
   std::complex<double> *X = generate_data(nSamples);
   if (!X)
      return 1;

   std::complex<double> *Xrec = new std::complex<double>[nSamples];
   std::copy(X, X + nSamples, Xrec);
   std::complex<double> *XparRec = new std::complex<double>[nSamples];
   std::copy(X, X + nSamples, XparRec);


   std::cout << "Size: 2^(" << power << ")" << std::endl;
   std::valarray<std::complex<double>> data(X, nSamples);
   std::valarray<std::complex<double>> output;
   double timer;


   {
      { // Iterative fft
	 auto timerC = scope_timer(timer);
	 output = iterative_fft(data);
      }
      print_res("iterative_fft", timer);
   }

   {
      timer = 0;
      { // Iterative parallel fft
	 auto timerC = scope_timer(timer);
	 auto tmp = output;
	 output = iterative_fft_parallel(data);
      }
      print_res("iterative_fft_tbb", timer);
   }

   {
      timer = 0;
      { // Iterative parallel openmpfft
	 auto timerC = scope_timer(timer);
	 auto tmp = output;
	 output = iterative_openmp_fft(data);
      }
      print_res("iterative_fft_openmp", timer);
   }

   {
      timer = 0;
      { // Recursive
	 auto timerC = scope_timer(timer);
	 difft2(Xrec, nSamples);
      }
      print_res("recursive_fft", timer);
   }

   {
      timer = 0;
      { // Recursive parallel
	 auto timerC = scope_timer(timer);
	 difft2_parallel(XparRec, nSamples);
      }
      print_res("Parallel Recursive FFT TBB", timer);
   }

   {
      timer = 0;
      { // Stockham FFT
	 auto timerC = scope_timer(timer);
	 stockham_fft(XparRec, nSamples);
      }
      print_res("stockham", timer);
   }

   {
      timer = 0;
      { // Stockham openmp FFT
	 auto timerC = scope_timer(timer);
	 stockham_openmp_fft(XparRec, nSamples);
      }
      print_res("stockham openmp", timer);
   }

   {
      timer = 0;
      { // FFTW
	 auto timerC = scope_timer(timer);
	 fftwBench(XparRec, nSamples);
      }
      print_res("FFTW", timer);
   }

   {
      timer = 0;
      {
	 auto timerC = scope_timer(timer);
	 cufft_example(XparRec, nSamples);
      }
      print_res("cuFFT", timer);
   }

   delete[] X;
   delete[] Xrec;
   delete[] XparRec;
   return 0;
}
