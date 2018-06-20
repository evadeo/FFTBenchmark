#include <complex>
#include <cuComplex.h>

__global__ void square_array(double *a, int N)
{
  int idx = blockIdx.x + threadIdx.x;
  if (idx < N)
    a[idx] = a[idx] * a[idx];
}


void do_square(double *a_h, const int N)
{
  double *a_d;
  size_t size = N * sizeof(double);
  /*
  for (int i = 0; i < N; ++i)
    a_d[i] = a_h[i] * a_h[i];
  */
  cudaMalloc((void **) &a_d, size);
  cudaMemcpy(a_d, a_h, size, cudaMemcpyHostToDevice);

  int blksize = 4;
  int n_blocks = N/blksize + (N % blksize == 0 ? 0 : 1);
  for (int i = 0; i < (N >> 15); ++i)
    square_array<<<n_blocks, blksize>>>(a_d, N);
  cudaMemcpy(a_h, a_d, size, cudaMemcpyDeviceToHost);
  cudaFree(a_d);
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

void difft2_cuda(std::complex<double>* arr, size_t n)
{
   if (n == 1)
      return;
   separate(arr, n);
   difft2_cuda(arr, n/2);
   difft2_cuda(arr + n/2, n/2);

   cuDoubleComplex *devArr;
   size_t size = n * sizeof(cuDoubleComplex);
   cudaMalloc((void **) &devArr, size);
   //cudaMemcpy(devArr, arrHost, size, cudaMemcpyHostToDevice);

   for (size_t k = 0; k < n/2; ++k)
   {
      std::complex<double> t = arr[k];
      auto my_exp = std::exp(std::complex<double>(0, -2. * M_PI * k/n));
      arr[k] = t + my_exp * arr[k + n/2];
      arr[k + n/2] = t - my_exp * arr[k + n/2];
   }
}

void my_cufft(std::complex<double>* arrHost, size_t n)
{
}
