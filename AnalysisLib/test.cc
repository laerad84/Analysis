#include <iostream>
#include "Environment.h"

int
main(){

  std::cout<< "Test" << std::endl;
  int i = GetEnvironment();
  PrintEnvironment();
  return 0;
}
