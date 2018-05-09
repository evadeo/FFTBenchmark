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
  unsigned int length = input.size();
  std::valarray<std::complex<double> > rev_array = bit_reverse_copy(input, length);

  std::valarray<std::complex<double> > omega (length / 2);

  omega[1] = std::polar(1., -2. * M_PI / length);
  omega[0] = 1;
  for (unsigned int i = 2; i < length / 2; ++i) {
    omega[i] = std::pow(omega[1], i);
  }
  int en = 1;
  int a = length / 2;

  for (int j = 0; j < my_log2(length); ++j) {
    for (unsigned int k = 0; k < length; ++k) {
      if (!(k & n)) {
        std::complex<double> firstTmp = rev_array[k];
        std::complex<double> secondTmp = omega[(k * a) % (n * a)] * rev_array[k + n];

        rev_array[k] = firstTmp + secondTmp;
        rev_array[k + n] = firstTemp - secondTmp;
      }
    }
    n *= 2;
    a /= 2;
  }

  return rev_array;
}

int main (int argc, char** argv) {
  const std::complex<double> test[] = { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
  std::valarray<std::complex<double> > data(test, 8);

  // forward fft
  std::valarray<std::complex<double> > res = iterative_fft(data);

  std::cout << "fft" << std::endl;
  for (int i = 0; i < 8; ++i)
  {
    std::cout << res[i] << std::endl;
  }
  return 0;
}
