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
   const int power = 23;
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
   { // Iterative fft
      auto timerC = scope_timer(timer);
      output = iterative_fft(data);
   }
   print_res("iterative_fft", timer);

   timer = 0;
   { // Recursive
      auto timerC = scope_timer(timer);
      difft2(Xrec, nSamples);
   }
   if (std::equal(std::begin(output),
		  std::end(output),
		  Xrec,
		  complex_close))
      print_res("recursive_fft", timer);
   else {
      std::cout << "Invalid recursive fft" << std::endl;
      for (unsigned int i = 0; i < nSamples; ++i) {
	 auto c = output[i] - Xrec[i];
	 if ((std::abs(c.real()) > 0.001) || (std::abs(c.imag()) > 0.001))
	    std::cout << output[i] << "\t" << Xrec[i] << '\t' << c.real() << '\t' << c.imag() << std::endl;
      }
   }

   timer = 0;
   { // Recursive parallel
      auto timerC = scope_timer(timer);
      difft2_parallel(XparRec, nSamples);
   }
   if (std::equal(std::begin(output),
		  std::end(output),
		  XparRec,
		  complex_close))
      print_res("Parallel Recursive FFT", timer);
   else {
      std::cout << "Invalid parallel recursive fft" << std::endl;
      for (unsigned int i = 0; i < nSamples; ++i) {
	 auto c = output[i] - XparRec[i];
	 if ((std::abs(c.real()) > 0.0001) || (std::abs(c.imag()) > 0.0001))
	    std::cout << output[i] << "\t" << XparRec[i] << '\t' << c.real() << '\t' << c.imag() << std::endl;
      }
   }
/*
   for (unsigned int i = 0; i < 15; ++i)
      std::cout << output[i] << '\n';
*/
   delete[] X;
   delete[] Xrec;
   delete[] XparRec;
   return 0;
}
