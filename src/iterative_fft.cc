#include <complex>
#include <valarray>
#include <iostream>

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
  unsigned int n = input.size();
  std::valarray<std::complex<double> > rev_array = bit_reverse_copy(input, n);

  std::valarray<std::complex<double> > omega (n / 2);

  omega[1] = std::polar(1., -2. * M_PI / n);
  omega[0] = 1;
  for (unsigned int i = 2; i < n / 2; ++i) {
    omega[i] = std::pow(omega[1], i);
  }
  int en = 1;
  int a = n / 2;

  for (int j = 0; j < my_log2(n); ++j) {
    for (unsigned int k = 0; k < n; ++k) {
      if (!(k & en)) {
        std::complex<double> temp = rev_array[k];
        std::complex<double> t = omega[(k * a) % (n * a)] * rev_array[k + en];

        rev_array[k] = temp + t;
        rev_array[k + en] = temp - t;
      }
    }
    en *= 2;
    a /= 2;
  }

  return rev_array;
}

int main (int argc, char** argv) {
  const std::complex<double> test[] = { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
  std::valarray<std::complex<double> > data(test, 8);
  std::cout << "go" << std::endl;
  // forward fft
  std::valarray<std::complex<double> > res = iterative_fft(data);

  std::cout << "fft" << std::endl;
  for (int i = 0; i < 8; ++i)
  {
    std::cout << res[i] << std::endl;
  }
  return 0;
}
