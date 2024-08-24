// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#pragma once

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <vector>


static inline std::string trace_err(const char* file, const char* func, long int line) {
   return "file: " + std::string(file) + ", function: " + std::string(func)
        + ", line: " + std::to_string(line);
}


#ifdef GEN_QUAD_DEBUG_ON
#define GEN_QUAD_ASSERT_DEBUG(condition) assert(condition)
#else
#define GEN_QUAD_ASSERT_DEBUG(condition) ((void)0)
#define EIGEN_NO_DEBUG
#endif


#define GEN_QUAD_ASSERT_ALWAYS(condition)                           \
   if(!(condition)) {                                               \
      std::cerr << trace_err(__FILE__, __func__, __LINE__) << ": "; \
      std::cerr << "Assertion " << (#condition) << " failed\n";     \
      std::exit(EXIT_FAILURE);                                      \
   }


#define GEN_QUAD_ASSERT_ALWAYS_MSG(condition, msg)                  \
   if(!(condition)) {                                               \
      std::cerr << trace_err(__FILE__, __func__, __LINE__) << ": "; \
      std::cerr << (msg) << "\n";                                   \
      std::exit(EXIT_FAILURE);                                      \
   }


#define GQ_THROW_RUNTIME_ERROR()                                                                 \
   do {                                                                                          \
      RuntimeError exc;                                                                          \
      exc.stack_trace.push_back("Exception thrown: " + trace_err(__FILE__, __func__, __LINE__)); \
      throw exc;                                                                                 \
   } while(0)


#define GQ_THROW_RUNTIME_ERROR_MSG(msg)                                                          \
   do {                                                                                          \
      RuntimeError exc(msg);                                                                     \
      exc.stack_trace.push_back("Exception thrown: " + trace_err(__FILE__, __func__, __LINE__)); \
      throw exc;                                                                                 \
   } while(0)


#define GQ_CATCH_AND_RETHROW()                                                                     \
   catch(ExceptionBase & exc) {                                                                    \
      exc.stack_trace.push_back("Exception rethrown: " + trace_err(__FILE__, __func__, __LINE__)); \
      throw;                                                                                       \
   }


// Note if exception is not ExceptionBase trace_err is not used, since it returns
// std::string, which uses memory allocation
#define GQ_CATCH_LAST_LEVEL()                                                                    \
   catch(ExceptionBase & exc) {                                                                  \
      std::cerr << std::endl;                                                                    \
      exc.stack_trace.push_back("Exception caught: " + trace_err(__FILE__, __func__, __LINE__)); \
      std::cerr << exc.what() << std::endl;                                                      \
      for(auto& tr : exc.stack_trace) std::cerr << tr << std::endl;                              \
   }                                                                                             \
   catch(std::exception & exc) {                                                                 \
      std::cerr << "file: " << __FILE__ << ", function: " << __func__ << ", line: " << __LINE__  \
                << std::endl;                                                                    \
      std::cerr << "std exception caught:" << std::endl;                                         \
      std::cerr << exc.what() << std::endl;                                                      \
   }                                                                                             \
   catch(...) {                                                                                  \
      std::cerr << "file: " << __FILE__ << ", function: " << __func__ << ", line: " << __LINE__  \
                << std::endl;                                                                    \
      std::cerr << "unknown exception caught:" << std::endl;                                     \
   }


namespace gquad {


class ExceptionBase : public std::exception {
public:
   explicit ExceptionBase(const std::string& msg) : msg_(msg) {
   }

   const char* what() const noexcept override {
      return msg_.c_str();
   }

   std::vector<std::string> stack_trace;

private:
   std::string msg_;
};


class RuntimeError : public ExceptionBase {
public:
   RuntimeError() : ExceptionBase{"RUNTIME ERROR"} {
   }

   explicit RuntimeError(const std::string& msg) : ExceptionBase{"RUNTIME ERROR: " + msg} {
   }
};


}  // namespace gquad

#endif

