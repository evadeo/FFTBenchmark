#include "timer.hh"
#include <complex>
#include <valarray>
#include <iostream>

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


int my_log2(int N) {
  int k = N;
  int i = 0;

  while(k) {
    k >>= 1;
    i++;
  }

  return i - 1;
}

unsigned int reverse(int N, int n) {
  int p = 0;
  for (int j = 0; j <= my_log2(N); ++j) {
    if (n & (1 << (my_log2(N) - j))) {
      p |= 1 << (j - 1);
    }
  }

  return p;
}

std::valarray<std::complex<double> > bit_reverse_copy(std::valarray<std::complex<double> > input, unsigned int length) {
  std::valarray<std::complex<double> > res (length);
  for (unsigned int i = 0; i < length; ++i) {
    res[reverse(length, i)] = input[i];
  }

  return res;
}

std::valarray<std::complex<double> > iterative_fft(std::valarray<std::complex<double> > input) {
  unsigned int length = input.size();
  std::valarray<std::complex<double> > rev_array = bit_reverse_copy(input, length);

  std::valarray<std::complex<double> > omega (length / 2);

  omega[1] = std::polar(1., -2. * M_PI / length);
  omega[0] = 1;
  for (unsigned int i = 2; i < length / 2; ++i) {
    omega[i] = std::pow(omega[1], i);
  }
  int n = 1;
  int a = length / 2;

  for (int j = 0; j < my_log2(length); ++j) {
    for (unsigned int k = 0; k < length; ++k) {
      if (!(k & n)) {
        std::complex<double> firstTmp = rev_array[k];
        std::complex<double> secondTmp = omega[(k * a) % (n * a)] * rev_array[k + n];

        rev_array[k] = firstTmp + secondTmp;
        rev_array[k + n] = firstTmp - secondTmp;
      }
    }
    n *= 2;
    a /= 2;
  }

  return rev_array;
}

int main (int argc, char** argv) {
   const unsigned int nSamples = std::pow(2,13);

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
   std::valarray<std::complex<double> > data(X, nSamples);
   // compute fft for this data
   double timer;
   std::cout << "iterative_fft" << std::endl;

   // forward fft
   {
      auto timerC = scope_timer(timer);
      iterative_fft(data);
   }

   std::cout << "----------------------" << std::endl;
   std::cout << "Execution time: " << timer << std::endl;

   std::cout << std::endl << "----------------------" << std::endl;


   timer = 0;
   std::cout << std::endl << std::endl << "recursive_fft" << std::endl;

   // forward fft
   {
      auto timerC = scope_timer(timer);
      difft2(X, nSamples);
   }

   std::cout << "----------------------" << std::endl;
   std::cout << "Execution time: " << timer << std::endl;
   std::cout << std::endl << "----------------------" << std::endl;


   delete[] x;
   delete[] X;
   return 0;
}
