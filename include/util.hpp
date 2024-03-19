// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

// Timer and output utilities.

#pragma once

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include <chrono>
#include <cstdio>
#include <iostream>
#include <string>

namespace gquad { namespace util {

////////////////////////////////////////////////////////////////////////////////////////////////////
struct timer {
   explicit timer() : start{std::chrono::high_resolution_clock::now()} {
   }

   void restart() {
      start = std::chrono::high_resolution_clock::now();
   }

   double wall_time() {
      return double(std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::high_resolution_clock::now() - start)
                        .count())
           / 1.e9;
   }

private:
   std::chrono::system_clock::time_point start;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
void print(const T x, const std::string& s) {
   std::cout << s << " = " << x << std::endl;
}

void print(const std::string& s);
void print(const char* s);

////////////////////////////////////////////////////////////////////////////////////////////////////
std::string n_spaces(long int n);
std::string index_of(std::string s, long int n);
std::string smart_spaces(long int max_ind, long int ind);

}}  // namespace gquad::util

#define GEN_QUAD_TIME(a)                                         \
   do {                                                          \
      gquad::util::timer t;                                      \
      a;                                                         \
      std::printf("%s took: %.4e seconds\n", #a, t.wall_time()); \
   } while(0)

#endif

