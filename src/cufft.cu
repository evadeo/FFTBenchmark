#include <complex>
#include <cuComplex.h>
#include <cufft.h>
#include <iostream>

void cufft_example(std::complex<double>* arrHost, size_t n, int batch)
{
  cufftHandle plan;
  cufftComplex *dataDevice;
  cudaMalloc((void**)&dataDevice, sizeof (cufftComplex) * n * batch);

  cudaMemcpy(dataDevice, (void*)arrHost, n, cudaMemcpyHostToDevice);
  if (cudaGetLastError() != cudaSuccess){
    fprintf(stderr, "Cuda error: %s, Failed to allocate\n", cudaGetErrorString(cudaGetLastError()));
    return;
  }

  if (cufftPlan1d(&plan, n, CUFFT_C2C, batch) != CUFFT_SUCCESS){
    fprintf(stderr, "CUFFT error: Plan creation failed");
    return;
  }

  /* Note:
   *  Identical pointers to input and output arrays implies in-place transformation
   */

  if (cufftExecC2C(plan, dataDevice, dataDevice, CUFFT_FORWARD) != CUFFT_SUCCESS){
    fprintf(stderr, "CUFFT error: ExecC2C Forward failed");
    return;
  }

  if (cudaDeviceSynchronize() != cudaSuccess){
    fprintf(stderr, "Cuda error: Failed to synchronize\n");
    return;
  }

  /*
   *  Divide by number of elements in dataDevice set to get back original data
   */

  cufftDestroy(plan);
  cudaFree(dataDevice);
}
