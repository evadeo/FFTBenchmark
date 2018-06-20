#include <complex>
#include <unordered_map>
#include <iostream>

//#include <omp.h> // OpenMP
#include <stdlib.h>
#include <immintrin.h>
#include <xmmintrin.h>

static std::complex<double>* alloc_vec(unsigned int SZ)
{
   void *v;
   // 32B aligned
   int res = posix_memalign(&v, 32, sizeof (std::complex<double>) * SZ);
   if (res)
      return nullptr;
   return (std::complex<double>*)v;
}

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

/*
void difft2_sse(std::complex<double>* arr, size_t n)
{
   if (n == 1)
      return;
   separate(arr, n);
   difft2(arr, n/2);
   difft2(arr + n/2, n/2);

   for (size_t k = 0; k < n/2; k += 4)
   {
      std::complex<double> t = arr[k];
      auto my_exp = std::exp(std::complex<double>(0, -2. * M_PI * k/n));
      arr[k] = t + my_exp * arr[k + n/2];
      arr[k + n/2] = t - my_exp * arr[k + n/2];
   }
}

*/

static __m256d mm256_mul_cplx(__m256d vec1, __m256d vec2)
{
   __m256d neg = _mm256_setr_pd(1.0,-1.0,1.0,-1.0);
   __m256d vec3 = _mm256_mul_pd(vec1, vec2);
   vec2 = _mm256_permute_pd(vec2, 0x5);
   vec2 = _mm256_mul_pd(vec2, neg);
   return _mm256_hsub_pd(vec3, _mm256_mul_pd(vec1, vec2));
}

void difft2_avx(std::complex<double>* arr, size_t n)
{
   if (n == 1)
      return;
   separate(arr, n);
   difft2_avx(arr, n/2);
   difft2_avx(arr + n/2, n/2);

   double* resTempAlignedArr = (double*)alloc_vec(4);

   for (size_t k = 0; k < n/2; k += 2)
   {
      __m256d tArr;
      __m256d tArrNext;
      if (k % 4 == 0)
      {
	 tArr = _mm256_load_pd((double*)arr+k); // ALIGN to 32
	 tArrNext = _mm256_setr_pd(arr[k+n/2].real(), arr[k+n/2].imag(),
					  arr[k+n/2+1].real(),arr[k+n/2+1].imag()); // ALIGN to 32
      }
      else
      {
	 tArr = _mm256_setr_pd(arr[k].real(), arr[k].imag(),
			       arr[k+1].real(),arr[k+1].imag()); // ALIGN to 32
	 tArrNext = _mm256_setr_pd(arr[k+n/2].real(), arr[k+n/2].imag(),
					  arr[k+n/2+1].real(),arr[k+n/2+1].imag()); // ALIGN to 32
      }

      auto my_exp0 = std::exp(std::complex<double>(0, -2. * M_PI * (k+0)/n));
      auto my_exp1 = std::exp(std::complex<double>(0, -2. * M_PI * (k+1)/n));
      __m256d expArr = _mm256_setr_pd(my_exp0.real(), my_exp0.imag(),
				     my_exp1.real(), my_exp1.imag());

      __m256d res1 = _mm256_add_pd(tArr, mm256_mul_cplx(expArr, tArrNext));
      __m256d res2 = _mm256_sub_pd(tArr, mm256_mul_cplx(expArr, tArrNext));
      //arr[k] = t + my_exp * arr[k + n/2];
      //arr[k + n/2] = t - my_exp * arr[k + n/2];

      if (k % 4 == 0)
      {
	 _mm256_store_pd((double*)arr+k, res1);
	 _mm256_store_pd(resTempAlignedArr, res2);
	 std::copy(resTempAlignedArr, resTempAlignedArr+4, arr+k+n/2);
      }
      else
      {
	 _mm256_store_pd(resTempAlignedArr, res1);
	 std::copy(resTempAlignedArr, resTempAlignedArr+4, arr+k);
	 _mm256_store_pd(resTempAlignedArr, res2);
	 std::copy(resTempAlignedArr, resTempAlignedArr+4, arr+k+n/2);
      }
   }
   free(resTempAlignedArr);
}
