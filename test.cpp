#include "tensor.h"
#include <iostream>

int main() {
  Tensor t(3,4,5);

  t[1][1][1] = 3;

  t[0][0][1] = 2;

  std::cout << t[0][0][0];
}