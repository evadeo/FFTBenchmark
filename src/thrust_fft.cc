#include <complex>
#include <thrust/complex.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>


// Separate even & odd elements
static void separate2 (thrust::complex<double> *a, size_t n) {
  //thrust::host_vector<thrust::complex<double>> b(n/2);
  thrust::complex<double>* b = new thrust::complex<double>[n/2];
   for (size_t i = 0; i < n/2; ++i) {
      b[i] = a[i*2 + 1];
      a[i] = a[i*2];
   }

   std::move(b, b + n/2, a + n/2);
   delete[] b;
}

struct doExpFunctor
{
  size_t N;
  __host__ __device__
  doExpFunctor(size_t N_)
  : N(N_)
  {}

  __host__ __device__
  thrust::complex<double> operator()(size_t idx)
  {
    return exp(thrust::complex<double>(0, -2. * M_PI * idx/N));
  }
};

struct dft_cuda_functor
{
  size_t N;

  __host__ __device__
  dft_cuda_functor(size_t N_)
  : N(N_)
  {}

  __host__ __device__
  thrust::complex<double> operator()(thrust)

};

struct applyDFT1
{
   size_t N;
   __host__ __device__
   applyDFT1(size_t N_)
      : N(N_)
   {}

  __host__ __device__
  void dft_cuda(thrust::device_vector<thrust::complex<double>> arr,
				    thrust::device_vector<thrust::complex<double>> exp_arr,
				    size_t N)
  {
  size_t idx = blockIdx.x + threadIdx.x;
  if (idx < N / 2)
  {
    size_t idx_next = blockIdx.x + threadIdx.x + N / 2;
    //auto my_exp = exp(thrust::complex<double>(0, -2. * M_PI * idx/N));
    arr[idx] = arr[idx] + exp_arr[idx] * arr[idx_next];
    arr[idx_next] = arr[idx] - exp_arr[idx] * arr[idx_next];
  }
}


void difft2_cuda(thrust::complex<double> *arr, size_t n)
{
   if (n == 1)
      return;
   separate2(arr, n);
   difft2_cuda(arr, n/2);
   difft2_cuda(arr + n/2, n/2);

   thrust::device_vector<thrust::complex<double>> devArr(arr, arr + n);
   thrust::device_vector<thrust::complex<double>> devExpArr(thrust::counting_iterator<size_t>(0),
							    thrust::counting_iterator<size_t>(n/2));
   thrust::transform(devExpArr.begin(), devExpArr.end(), doExpFunctor(n));

   dft_cuda(devArr, n);
   thrust::copy(devArr.begin(), devArr.end(), arr);

   //dft_cuda(arr, n/2);
   /*
   for (size_t k = 0; k < n/2; ++k)
   {
      std::complex<double> t = arr[k];
      auto my_exp = std::exp(std::complex<double>(0, -2. * M_PI * k/n));
      arr[k] = t + my_exp * arr[k + n/2];
      arr[k + n/2] = t - my_exp * arr[k + n/2];
   }
   */
}


void difft2_cuda_head(std::complex<double>* arr, size_t n)
{
  //thrust::host_vector<thrust::complex<double>> thrArr(n);
  thrust::complex<double> *thrArr = new thrust::complex<double>[n];
  for (size_t i = 0; i < n; ++i)
    thrArr[i] = thrust::complex<double>(arr[i].real(), arr[i].imag());

  difft2_cuda(thrArr, n);

  for (size_t i = 0; i < n; ++i)
    arr[i] = std::complex<double>(thrArr[i].real(), thrArr[i].imag());
}
