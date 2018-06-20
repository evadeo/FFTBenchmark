#include <complex>
#include <unordered_map>
#include <iostream>


// Separate even & odd elements
static void separate (std::complex<double>* a, size_t n) {
   std::complex<double>* b = new std::complex<double>[n/2];
   for (size_t i = 0; i < n/2; ++i) {
      b[i] = a[i*2 + 1];
      a[i] = a[i*2];
   }

   std::move(b, b + n/2, a + n/2);
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
      auto my_exp = std::exp(std::complex<double>(0, -2. * M_PI * k/n));
      arr[k] = t + my_exp * arr[k + n/2];
      arr[k + n/2] = t - my_exp * arr[k + n/2];
   }
}
