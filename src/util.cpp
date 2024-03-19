#if __cplusplus < 201402L
#error requires c++14 or higher

#else

#include "../include/util.hpp"

#include <cmath>
#include <iostream>

namespace gquad { namespace util {

static long int count_digit(long int number) {
   return number == 0L ? 1L : static_cast<long int>(std::log10(std::abs(number)) + 1);
}

void print(const std::string& s) {
   std::cout << s << std::endl;
}

void print(const char* s) {
   std::cout << s << std::endl;
}

std::string n_spaces(long int n) {
   return std::string(n, 32);
}

std::string index_of(std::string s, long int n) {
   s += "[";
   s += std::to_string(n);
   s += "]";
   return s;
}

std::string smart_spaces(long int max_ind, long int ind) {
   long int max_digits = count_digit(max_ind);
   long int ind_digits = count_digit(ind);
   return n_spaces(1 + max_digits - ind_digits);
}

}}  // namespace gquad::util

#endif

