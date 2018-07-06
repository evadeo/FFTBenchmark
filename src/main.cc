#include "fft.hh"
#include "timer.hh"
#include <iostream>
#include <complex>
#include <valarray>
#include <algorithm>
#include <stdlib.h>

static std::complex<double>* alloc_vec(unsigned int SZ)
{
   void *v;
   // 32B aligned
   int res = posix_memalign(&v, 32, sizeof (std::complex<double>) * SZ);
   if (res)
      return nullptr;
   return (std::complex<double>*)v;
}

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
   int nRepeat = 1;
   const int power = 23;
   const unsigned int nSamples = std::pow(2, power);
   std::complex<double> *X = generate_data(nSamples);
   if (!X)
      return 1;

   std::complex<double> *Xrec = (std::complex<double>*)alloc_vec(nSamples);
   std::copy(X, X + nSamples, Xrec);
   std::complex<double> *XparRec = new std::complex<double>[nSamples];
   std::copy(X, X + nSamples, XparRec);


   std::cout << "Size: 2^(" << power << ")" << std::endl;
   std::valarray<std::complex<double>> data(X, nSamples);
   std::valarray<std::complex<double>> output;
   double timer;
   bool silent = true;

   double timer_cufft = 0;
   double timer_stockham = 0;
   double timer_iterative = 0;
   double timer_recursive = 0;
   double timer_iterative_parallel = 0;
   double timer_fftw = 0;

   for (int rep = 0; rep < nRepeat; rep++)
   {
      {
	 timer = 0;
	 {
	    auto timerC = scope_timer(timer);
	    cufft_example(XparRec, nSamples);
	 }
	 print_res("cuFFT", timer);
	 timer_cufft += timer;
      }

      /* DOUBLE FREE PROBLEM
	 {
	 timer = 0;
	 { // Recursive AVX
	 auto timerC = scope_timer(timer);
	 difft2_avx(Xrec, nSamples);
	 }
	 if (std::equal(std::begin(output),
	 std::end(output),
	 Xrec,
	 complex_close))
	 print_res("recursive_fft AVX", timer);
	 else {
	 std::cout << "Invalid recursive fft AVX" << std::endl;
	 print_res("recursive_fft AVX", timer);
	 if (!silent)
	 {
	 for (unsigned int i = 0; i < nSamples; ++i) {
	 auto c = output[i] - Xrec[i];
	 if ((std::abs(c.real()) > 0.001) || (std::abs(c.imag()) > 0.001))
	 std::cout << output[i] << "\t" << Xrec[i] << '\t' << c.real() << '\t' << c.imag() << std::endl;
	 }
	 }
	 }
	 }
      */

      {
	 timer = 0;
	 { // Stockham FFT
	    auto timerC = scope_timer(timer);
	    stockham_fft(XparRec, nSamples);
	 }
	 timer_stockham += timer / (double)nRepeat;
	 print_res("stockham", timer);
      }

      {
	 timer = 0;
	 { // Stockham FFT AVX
	    auto timerC = scope_timer(timer);
	    stockham_fft_avx(XparRec, nSamples);
	 }
	 print_res("stockham AVX", timer);
      }


      {
	 { // Iterative fft
	    auto timerC = scope_timer(timer);
	    output = iterative_fft(data);
	 }
	 timer_iterative += timer / (double)nRepeat;
	 print_res("iterative_fft", timer);
      }

      {
	 timer = 0;
	 { // Iterative parallel fft
	    auto timerC = scope_timer(timer);
	    auto tmp = output;
	    output = iterative_fft_parallel(data);

	    if (!std::equal(std::begin(output), std::end(output), std::begin(tmp), std::end(tmp)))
	       std::cout << "invalid results for iterative parallel fft" << std::endl;
	 }
	 timer_iterative_parallel += timer / (double)nRepeat;
	 print_res("iterative_fft_parallel", timer);
      }



      {
	 timer = 0;
	 { // Recursive
	    auto timerC = scope_timer(timer);
	    difft2(Xrec, nSamples);
	 }
	 timer_recursive += timer / (double)nRepeat;
	 if (std::equal(std::begin(output),
			std::end(output),
			Xrec,
			complex_close))
	    print_res("recursive_fft", timer);
	 else {
	    std::cout << "Invalid recursive fft" << std::endl;
	    print_res("recursive_fft", timer);
	    if (!silent)
	       for (unsigned int i = 0; i < nSamples; ++i) {
		  auto c = output[i] - Xrec[i];
		  if ((std::abs(c.real()) > 0.001) || (std::abs(c.imag()) > 0.001))
		     std::cout << output[i] << "\t" << Xrec[i]
			       << '\t' << c.real() << '\t' << c.imag() << std::endl;
	       }
	 }
      }

      {
	 timer = 0;
	 { // Recursive parallel
	    auto timerC = scope_timer(timer);
	    difft2_parallel(XparRec, nSamples);
	 }
	 //timer_iterative_parallel += timer / (double)nRepeat;
	 if (std::equal(std::begin(output),
			std::end(output),
			XparRec,
			complex_close))
	    print_res("Parallel Recursive FFT", timer);
	 else {
	    std::cout << "Invalid parallel recursive fft" << std::endl;
	    if (!silent)
	       for (unsigned int i = 0; i < nSamples; ++i) {
		  auto c = output[i] - XparRec[i];
		  if ((std::abs(c.real()) > 0.0001) || (std::abs(c.imag()) > 0.0001))
		     std::cout << output[i] << "\t" << XparRec[i] << '\t' << c.real()
			       << '\t' << c.imag() << std::endl;
	       }
	 }
      }


      {
	 timer = 0;
	 { // FFTW
	    auto timerC = scope_timer(timer);
	    fftwBench(XparRec, nSamples);
	 }
	 timer_fftw += timer / (double)nRepeat;
	 print_res("FFTW", timer);
      }

      {
	 timer = 0;
	 { // FFTW Parallel
	    auto timerC = scope_timer(timer);
	    fftwBenchParr(XparRec, nSamples);
	 }
	 print_res("FFTW Parallel", timer);
      }
   }

   if (nRepeat > 1)
   {
      print_res("cuFFT - final", timer_cufft / nRepeat);
      print_res("stockham - final", timer_stockham);
      print_res("iterative - final", timer_iterative);
      print_res("iterative_parallel - final", timer_iterative_parallel);
      print_res("recursive - final", timer_recursive);
      print_res("fftw - final", timer_fftw);
   }

   delete[] X;
   delete[] Xrec;
   delete[] XparRec;
   return 0;
}
