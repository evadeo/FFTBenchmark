#include <complex>
#include <iostream>

#define M_PI 3.14159265358979323846

// Separate even & odd elements
void separate (std::complex<double>* a, int n) {
   std::complex<double>* b = new std::complex<double>[n/2];
   // FIXME: Merge the two loops
   for (int i = 0; i < n/2; ++i)
      b[i] = a[i*2 + 1];
   for (int i = 0; i < n/2; ++i)
      a[i] = a[i*2];
   std::move(b, b + n/2, a + n/2); //FIXME: Make it parallel
   delete[] b;
}

void difft2(std::complex<double>* arr, int n)
{
   if (n == 1)
      return;
   separate(arr, n);
   difft2(arr, n/2);
   difft2(arr + n/2, n/2);

   for (int k = 0; k < n/2; ++k)
   {
      std::complex<double> t = arr[k];
      auto x = std::exp(std::complex<double>(0, -2. * M_PI * k / n));
      arr[k] = t + x * arr[k + n/2];
      arr[k + n/2] = t - x * arr[k + n/2];
   }
}


int main() {
   const int nSamples = 8;
   std::complex<double> samples[nSamples] = {1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0};
   //std::complex<double> result[nSamples];

   difft2(samples, nSamples);

   std::cout << "fft" << std::endl;
   for (int i = 0; i < 8; ++i)
      std::cout << samples[i] << " ";
   std::cout << std::endl;
   return 0;
}
