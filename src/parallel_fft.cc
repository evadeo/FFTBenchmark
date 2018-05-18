#include <complex>
#include <valarray>
#include <tbb/task_group.h>


// Separate even & odd elements
void separate_parallel(std::complex<double>* a, size_t n) {
   std::complex<double>* b = new std::complex<double>[n/2];
   for (size_t i = 0; i < n/2; ++i) {
      b[i] = a[i*2 + 1];
      a[i] = a[i*2];
   }

   std::move(b, b + n/2, a + n/2); //FIXME: Make it parallel
   delete[] b;
}


void difft2_parallel(std::complex<double>* arr, size_t n)
{
   if (n == 1)
      return;
   tbb::task_group *sp_group = new tbb::task_group;
   separate_parallel(arr, n);
   sp_group->run([arr, n]{ difft2_parallel(arr, n/2); });
   sp_group->run([arr, n]{ difft2_parallel(arr + n/2, n/2); });
   sp_group->wait();

   for (size_t k = 0; k < n/2; ++k)
   {
      std::complex<double> t = arr[k];
      auto x = std::exp(std::complex<double>(0, -2. * M_PI * k / n));
      arr[k] = t + x * arr[k + n/2];
      arr[k + n/2] = t - x * arr[k + n/2];
   }
}
