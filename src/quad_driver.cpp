// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "../include/quad_driver.hpp"

#include "../include/compute_quad.hpp"
#include "../include/function.hpp"
#include "../include/node_elimination.hpp"
#include "../include/nonlinear_solve.hpp"
#include "../include/util.hpp"


#define GEN_QUAD_ASSERT_DEG_DIM_RANGE(degree, dimension)                             \
   GEN_QUAD_ASSERT_ALWAYS_MSG(is_reasonable_degree_and_dimension(degree, dimension), \
                              "Specified degree and/or dimension is "                \
                              "outside of currently acceptable range.\n\n");


#define GEN_QUAD_ASSERT_DEG_DIM_RANGE_VECTOR(degree, dimension) \
   for(auto dg : degree) {                                      \
      GEN_QUAD_ASSERT_DEG_DIM_RANGE(dg, dimension);             \
   }


#define GEN_QUAD_INVOKE_DRIVER_ALL(is_interval, function, ...)                                     \
   try {                                                                                           \
      if(gen_quad_generate_driver_all<is_interval>(&function, __VA_ARGS__) != 0) {                 \
         std::cerr << "file: " << __FILE__ << ", function: " << __func__ << ", line: " << __LINE__ \
                   << ": did not successfully compute all quadratures" << std::endl                \
                   << std::endl;                                                                   \
         return 1;                                                                                 \
      } else {                                                                                     \
         std::cout << std::endl << std::endl;                                                      \
         return 0;                                                                                 \
      }                                                                                            \
   }                                                                                               \
   GQ_CATCH_LAST_LEVEL();                                                                          \
   return 0;


#define GEN_QUAD_INVOKE_DRIVER(function, object, ...)                                              \
   try {                                                                                           \
      auto qh = function(__VA_ARGS__);                                                             \
      if(!qh) {                                                                                    \
         std::cerr << "file: " << __FILE__ << ", function: " << __func__ << ", line: " << __LINE__ \
                   << ": did not successfully compute quadrature" << std::endl                     \
                   << std::endl;                                                                   \
         return nullptr;                                                                           \
      } else {                                                                                     \
         std::cout << std::endl << std::endl;                                                      \
         auto q = std::make_unique<object>(qh->first);                                             \
         return q;                                                                                 \
      }                                                                                            \
   }                                                                                               \
   GQ_CATCH_LAST_LEVEL();                                                                          \
   return nullptr;


namespace gquad {


////////////////////////////////////////////////////////////////////////////////////////////////////
static std::unique_ptr<std::pair<QuadInterval, StdVector<History>>> QuadDriverInterval(gq_int deg);


static std::unique_ptr<std::pair<QuadCube, StdVector<History>>> QuadDriverCube(gq_int deg, gq_int dim,
                                                                               SearchWidth search_width);

static std::unique_ptr<std::pair<QuadSimplex, StdVector<History>>> QuadDriverSimplex(
    gq_int deg, gq_int dim, SearchWidth search_width);


static std::unique_ptr<std::pair<QuadPyramid3D, StdVector<History>>> QuadDriverPyramid3D(
    gq_int deg, SearchWidth search_width);


static std::unique_ptr<std::pair<QuadCubeSimplex, StdVector<History>>> QuadDriverCubeSimplex(
    gq_int deg, gq_int dim1, gq_int dim2, SearchWidth search_width);


static std::unique_ptr<std::pair<QuadSimplexSimplex, StdVector<History>>> QuadDriverSimplexSimplex(
    gq_int deg, gq_int dim1, gq_int dim2, SearchWidth search_width);


static std::unique_ptr<std::pair<QuadOmega2D, StdVector<History>>> QuadDriverPentagon(
    gq_int deg, SearchWidth search_width);


static std::unique_ptr<std::pair<QuadOmega2D, StdVector<History>>> QuadDriverHexagon(
    gq_int deg, SearchWidth search_width);


static std::unique_ptr<std::pair<QuadOmega2D, StdVector<History>>> QuadDriverOmega2D(
    gq_int deg, Omega2D omega, SearchWidth search_width);


////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename D, typename F, typename... Fargs>
static std::unique_ptr<std::pair<D, StdVector<History>>> gen_quad_generate_driver(F f, Fargs&&... fargs);


template <bool is_interval, typename F, typename... DimArgs>
static gq_int gen_quad_generate_driver_all(F func_ptr, const StdVector<gq_int>& deg, DimArgs... dargs);


////////////////////////////////////////////////////////////////////////////////////////////////////
static void OutputAll(double ttotal, QuadDomain& q, StdVector<History>& h);
static void HistToFile(const QuadDomain& q, const StdVector<History>& hist);
static void TimesToFile(double total_time, const QuadDomain& q);
static void TimesToScreen(double total_time);
static void EfficiencyToFile(const std::pair<StdVector<double>, StdVector<std::string>>&, const std::string&);


////////////////////////////////////////////////////////////////////////////////////////////////////
static void PrintDriverInfo(gq_int deg, gq_int dim, SearchWidth search_width, std::string domain_name);
static void PrintDebugAndOmpInfo();
static gq_int hours(double x);
static gq_int minutes(double x);
static double seconds(double x);
static void ResetTimers();


////////////////////////////////////////////////////////////////////////////////////////////////////
gq_int ComputeAndOutputIntervals(const StdVector<gq_int>& deg) {
   GEN_QUAD_ASSERT_DEG_DIM_RANGE_VECTOR(deg, 1);
   GEN_QUAD_INVOKE_DRIVER_ALL(true, QuadDriverInterval, deg);
}


gq_int ComputeAndOutputCubes(const StdVector<gq_int>& deg, gq_int dim, SearchWidth search_width) {
   GEN_QUAD_ASSERT_DEG_DIM_RANGE_VECTOR(deg, dim);
   GEN_QUAD_INVOKE_DRIVER_ALL(false, QuadDriverCube, deg, dim, search_width);
}


gq_int ComputeAndOutputSimplexes(const StdVector<gq_int>& deg, gq_int dim, SearchWidth search_width) {
   GEN_QUAD_ASSERT_DEG_DIM_RANGE_VECTOR(deg, dim);
   GEN_QUAD_INVOKE_DRIVER_ALL(false, QuadDriverSimplex, deg, dim, search_width);
}


gq_int ComputeAndOutputPyramids3D(const StdVector<gq_int>& deg, SearchWidth search_width) {
   GEN_QUAD_ASSERT_DEG_DIM_RANGE_VECTOR(deg, 3);
   GEN_QUAD_INVOKE_DRIVER_ALL(false, QuadDriverPyramid3D, deg, search_width);
}


gq_int ComputeAndOutputCubeSimplexes(const StdVector<gq_int>& deg, gq_int dim1, gq_int dim2,
                                     SearchWidth search_width) {
   GEN_QUAD_ASSERT_DEG_DIM_RANGE_VECTOR(deg, dim1 + dim2);
   GEN_QUAD_INVOKE_DRIVER_ALL(false, QuadDriverCubeSimplex, deg, dim1, dim2, search_width);
}


gq_int ComputeAndOutputSimplexSimplexes(const StdVector<gq_int>& deg, gq_int dim1, gq_int dim2,
                                        SearchWidth search_width) {
   GEN_QUAD_ASSERT_DEG_DIM_RANGE_VECTOR(deg, dim1 + dim2);
   GEN_QUAD_INVOKE_DRIVER_ALL(false, QuadDriverSimplexSimplex, deg, dim1, dim2, search_width);
}


gq_int ComputeAndOutputPentagons(const StdVector<gq_int>& deg, SearchWidth search_width) {
   GEN_QUAD_ASSERT_DEG_DIM_RANGE_VECTOR(deg, 2);
   GEN_QUAD_INVOKE_DRIVER_ALL(false, QuadDriverPentagon, deg, search_width);
}


gq_int ComputeAndOutputHexagons(const StdVector<gq_int>& deg, SearchWidth search_width) {
   GEN_QUAD_ASSERT_DEG_DIM_RANGE_VECTOR(deg, 2);
   GEN_QUAD_INVOKE_DRIVER_ALL(false, QuadDriverHexagon, deg, search_width);
}


gq_int ComputeAndOutputOmega2D(Omega2D omega, const StdVector<gq_int>& deg, SearchWidth search_width) {
   GEN_QUAD_ASSERT_DEG_DIM_RANGE_VECTOR(deg, 2);
   GEN_QUAD_INVOKE_DRIVER_ALL(false, QuadDriverOmega2D, deg, std::move(omega), search_width);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
std::unique_ptr<QuadInterval> Quadrature_Interval(gq_int deg) {
   GEN_QUAD_INVOKE_DRIVER(QuadDriverInterval, QuadInterval, deg);
}


std::unique_ptr<QuadPyramid3D> Quadrature_Pyramid3D(gq_int deg, SearchWidth search_width) {
   GEN_QUAD_INVOKE_DRIVER(QuadDriverPyramid3D, QuadPyramid3D, deg, search_width);
}


std::unique_ptr<QuadCube> Quadrature_Cube(gq_int deg, gq_int dim, SearchWidth search_width) {
   GEN_QUAD_INVOKE_DRIVER(QuadDriverCube, QuadCube, deg, dim, search_width);
}


std::unique_ptr<QuadSimplex> Quadrature_Simplex(gq_int deg, gq_int dim, SearchWidth search_width) {
   GEN_QUAD_INVOKE_DRIVER(QuadDriverSimplex, QuadSimplex, deg, dim, search_width);
}


std::unique_ptr<QuadCubeSimplex> Quadrature_CubeSimplex(gq_int deg, gq_int dim1, gq_int dim2,
                                                        SearchWidth search_width) {
   GEN_QUAD_INVOKE_DRIVER(QuadDriverCubeSimplex, QuadCubeSimplex, deg, dim1, dim2, search_width);
}


std::unique_ptr<QuadSimplexSimplex> Quadrature_SimplexSimplex(gq_int deg, gq_int dim1, gq_int dim2,
                                                              SearchWidth search_width) {
   GEN_QUAD_INVOKE_DRIVER(QuadDriverSimplexSimplex, QuadSimplexSimplex, deg, dim1, dim2, search_width);
}


std::unique_ptr<QuadOmega2D> Quadrature_Pentagon(gq_int deg, SearchWidth search_width) {
   GEN_QUAD_INVOKE_DRIVER(QuadDriverPentagon, QuadOmega2D, deg, search_width);
}


std::unique_ptr<QuadOmega2D> Quadrature_Hexagon(gq_int deg, SearchWidth search_width) {
   GEN_QUAD_INVOKE_DRIVER(QuadDriverHexagon, QuadOmega2D, deg, search_width);
}


std::unique_ptr<QuadOmega2D> Quadrature_Omega2D(Omega2D omega, gq_int deg, SearchWidth search_width) {
   GEN_QUAD_INVOKE_DRIVER(QuadDriverOmega2D, QuadOmega2D, deg, std::move(omega), search_width);
}


//////////////////////////////////////////////////////////////////////////////////////////////////
static std::unique_ptr<std::pair<QuadInterval, StdVector<History>>> QuadDriverInterval(gq_int deg) {
   GEN_QUAD_ASSERT_ALWAYS(deg > 0);
   GEN_QUAD_ASSERT_DEG_DIM_RANGE(deg, 1);
   PrintDriverInfo(deg, 1, SearchWidth{1}, "interval");
   ResetTimers();
   return std::make_unique<std::pair<QuadInterval, StdVector<History>>>(
       std::pair<QuadInterval, StdVector<History>>{QuadInterval{ComputeInterval(deg)}, StdVector<History>{}});
}


static std::unique_ptr<std::pair<QuadCube, StdVector<History>>> QuadDriverCube(gq_int deg, gq_int dim,
                                                                               SearchWidth search_width) {
   GEN_QUAD_ASSERT_ALWAYS(deg > 0);
   GEN_QUAD_ASSERT_ALWAYS(dim > 1);
   GEN_QUAD_ASSERT_DEG_DIM_RANGE(deg, dim);
   PrintDriverInfo(deg, dim, search_width, "cube");
   ResetTimers();
   return gen_quad_generate_driver<QuadCube>(ComputeCube, deg, dim, search_width.width);
}


static std::unique_ptr<std::pair<QuadSimplex, StdVector<History>>> QuadDriverSimplex(
    gq_int deg, gq_int dim, SearchWidth search_width) {
   GEN_QUAD_ASSERT_ALWAYS(deg > 0);
   GEN_QUAD_ASSERT_ALWAYS(dim > 1);
   GEN_QUAD_ASSERT_DEG_DIM_RANGE(deg, dim);
   PrintDriverInfo(deg, dim, search_width, "simplex");
   ResetTimers();
   return gen_quad_generate_driver<QuadSimplex>(ComputeSimplex, deg, dim, search_width.width);
}


static std::unique_ptr<std::pair<QuadPyramid3D, StdVector<History>>> QuadDriverPyramid3D(
    gq_int deg, SearchWidth search_width) {
   GEN_QUAD_ASSERT_ALWAYS(deg > 0);
   GEN_QUAD_ASSERT_DEG_DIM_RANGE(deg, 3);
   PrintDriverInfo(deg, 3, search_width, "pyramid3D");
   ResetTimers();
   return gen_quad_generate_driver<QuadPyramid3D>(ComputePyramid3D, deg, search_width.width);
}


static std::unique_ptr<std::pair<QuadCubeSimplex, StdVector<History>>> QuadDriverCubeSimplex(
    gq_int deg, gq_int dim1, gq_int dim2, SearchWidth search_width) {
   GEN_QUAD_ASSERT_ALWAYS(deg > 0);
   GEN_QUAD_ASSERT_ALWAYS(dim1 > 0);
   GEN_QUAD_ASSERT_ALWAYS(dim2 > 1);
   GEN_QUAD_ASSERT_DEG_DIM_RANGE(deg, dim1 + dim2);
   PrintDriverInfo(deg, dim1 + dim2, search_width, "cubesimplex");
   ResetTimers();
   return gen_quad_generate_driver<QuadCubeSimplex>(ComputeCubeSimplex, deg, dim1, dim2, search_width.width);
}


static std::unique_ptr<std::pair<QuadSimplexSimplex, StdVector<History>>> QuadDriverSimplexSimplex(
    gq_int deg, gq_int dim1, gq_int dim2, SearchWidth search_width) {
   GEN_QUAD_ASSERT_ALWAYS(deg > 0);
   GEN_QUAD_ASSERT_ALWAYS(dim1 > 1);
   GEN_QUAD_ASSERT_ALWAYS(dim2 > 1);
   GEN_QUAD_ASSERT_DEG_DIM_RANGE(deg, dim1 + dim2);
   PrintDriverInfo(deg, dim1 + dim2, search_width, "simplexsimplex");
   ResetTimers();
   return gen_quad_generate_driver<QuadSimplexSimplex>(
       ComputeSimplexSimplex, deg, dim1, dim2, search_width.width);
}


static std::unique_ptr<std::pair<QuadOmega2D, StdVector<History>>> QuadDriverPentagon(
    gq_int deg, SearchWidth search_width) {
   GEN_QUAD_ASSERT_ALWAYS(deg > 0);
   GEN_QUAD_ASSERT_DEG_DIM_RANGE(deg, 2);
   PrintDriverInfo(deg, 2, search_width, "pentagon");
   ResetTimers();
   return gen_quad_generate_driver<QuadOmega2D>(ComputePentagon, deg, search_width.width);
}


static std::unique_ptr<std::pair<QuadOmega2D, StdVector<History>>> QuadDriverHexagon(
    gq_int deg, SearchWidth search_width) {
   GEN_QUAD_ASSERT_ALWAYS(deg > 0);
   GEN_QUAD_ASSERT_DEG_DIM_RANGE(deg, 2);
   PrintDriverInfo(deg, 2, search_width, "hexagon");
   ResetTimers();
   return gen_quad_generate_driver<QuadOmega2D>(ComputeHexagon, deg, search_width.width);
}


static std::unique_ptr<std::pair<QuadOmega2D, StdVector<History>>> QuadDriverOmega2D(
    gq_int deg, Omega2D omega, SearchWidth search_width) {
   GEN_QUAD_ASSERT_ALWAYS(deg > 0);
   GEN_QUAD_ASSERT_DEG_DIM_RANGE(deg, 2);
   PrintDriverInfo(deg, 2, search_width, omega.domain_name());
   ResetTimers();
   return gen_quad_generate_driver<QuadOmega2D>(ComputeOmega2D, std::move(omega), deg, search_width.width);
}


////////////////////////////////////////////////////////////////////////////////////////////////
template <typename D, typename F, typename... Fargs>
static std::unique_ptr<std::pair<D, StdVector<History>>> gen_quad_generate_driver(F f, Fargs&&... fargs) {
   auto quad_hist = f(std::forward<Fargs>(fargs)...);
   return std::make_unique<decltype(quad_hist)>(quad_hist);
}


template <bool is_interval, typename F, typename... DimArgs>
static gq_int gen_quad_generate_driver_all(F func_ptr, const StdVector<gq_int>& deg, DimArgs... dargs) {
   StdVector<double> eff_indexes;
   StdVector<std::string> quad_names;
   std::string dims_file_name;

   for(auto dg : deg) {
      util::timer t;
      auto qh = (*func_ptr)(dg, dargs...);
      double ttotal = t.wall_time();

      if(!qh) {
         return 1;
      }

      if(!is_interval) {
         OutputAll(ttotal, qh->first, qh->second);
      } else {
         QuadToFile(qh->first);
      }

      eff_indexes.push_back(qh->first.efficiency());
      quad_names.push_back(qh->first.quad_file_name());
      dims_file_name = qh->first.quad_dims_file_name();
   }
   EfficiencyToFile({eff_indexes, quad_names}, dims_file_name);
   return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
void QuadToFile(const QuadDomain& q) {
   std::string quad_path = std::string("results/quad_rules/") + q.quad_file_name();
   std::ofstream ofs{quad_path};
   if(!ofs) {
      GQ_THROW_RUNTIME_ERROR_MSG("could not open file for writing quadrature");
   }
   ofs << q << std::endl;
}


static void OutputAll(double ttotal, QuadDomain& q, StdVector<History>& h) {
   QuadToFile(q);
   HistToFile(q, h);
   TimesToFile(ttotal, q);
   TimesToScreen(ttotal);
   util::print(q.relative_exponential_residual(), "relative exponential residual");
}


static void HistToFile(const QuadDomain& q, const StdVector<History>& hist) {
   std::string file_path = std::string("results/history/hist_") + q.quad_file_name();
   const char* c_file_path = file_path.c_str();
   FILE* file = std::fopen(c_file_path, "w");
   if(!file) {
      GQ_THROW_RUNTIME_ERROR_MSG("could not open file for writing history");
   }

   if(hist.empty()) {
      std::fprintf(file, "Empty History\n\n");
      return;
   }

   for(const History& h : hist) {
      std::fprintf(file, "*********************************************************\n\n");
      std::fprintf(file, "dimension                 = %lu\n", h.dim);
      std::fprintf(file, "degree of precision       = %lu\n", h.deg);
      std::fprintf(file, "domain type               = %s\n", h.domain_name.c_str());
      std::fprintf(file, "number of basis functions = %lu\n", h.num_funcs);
      std::fprintf(file, "initial number of nodes   = %lu\n", h.nodes_initial);
      std::fprintf(file, "final number of nodes     = %lu\n", h.nodes_final);
      std::fprintf(file, "optimal number of nodes   = %lu\n", h.nodes_optimal);
      std::fprintf(file, "total eliminations        = %lu\n", h.total_elims);
      std::fprintf(file, "final residual            = %.16e\n", h.res);
      std::fprintf(file, "efficiency                = %lf\n\n", h.efficiency);

      if(h.total_elims <= 0) {
         std::fprintf(file, "no nodes eliminated\n");
      }

      for(gq_int j = 0; j < h.total_elims; ++j) {
         if(h.hlist[j].num_solutions == 1) {
            std::fprintf(file, "%lu solution in wide search\n", h.hlist[j].num_solutions);
         } else {
            std::fprintf(file, "%lu solutions in wide search\n", h.hlist[j].num_solutions);
         }
         if(h.hlist[j].num_fails == 1) {
            std::fprintf(file, "%lu failure in wide search\n", h.hlist[j].num_fails);
         } else {
            std::fprintf(file, "%lu failures in wide search\n", h.hlist[j].num_fails);
         }

         std::fprintf(file, "number of nodes = %lu\n", h.hlist[j].nodes_tot);
         std::fprintf(file, "success_node[%lu] = %lu\n", j, h.hlist[j].success_node);
         std::fprintf(file, "Converged in %lu iterations\n\n", h.hlist[j].success_its);
      }
      std::fprintf(file, "\n");
   }

   std::fclose(file);
}


static void TimesToFile(double total_time, const QuadDomain& q) {
   std::string file_path = std::string("results/times/times_") + q.quad_file_name();
   const char* c_file_path = file_path.c_str();
   FILE* file = std::fopen(c_file_path, "w");
   if(!file) {
      GQ_THROW_RUNTIME_ERROR_MSG("could not open file for writing times");
   }

   double lsq_time_total = LsqSolveTimer::lsq_time_total;
   double predictor_time_total = PredictorTimer::predictor_time_total;
   double jacobian_time_total = EvalJacobian::jacobian_time_total;
   double function_time_total = EvalFunction::function_time_total;
   double main_routines = lsq_time_total + jacobian_time_total + function_time_total + predictor_time_total;

   std::fprintf(file,
                "wall clock time for least squares routine in LeastSquaresNewton: %luh, %lum, %.2fs\n",
                hours(lsq_time_total),
                minutes(lsq_time_total),
                seconds(lsq_time_total));

   std::fprintf(file,
                "wall clock time for GetJacobian routine: %luh, %lum, %.2fs\n",
                hours(jacobian_time_total),
                minutes(jacobian_time_total),
                seconds(jacobian_time_total));

   std::fprintf(file,
                "wall clock time for GetFunction routine: %luh, %lum, %.2fs\n",
                hours(function_time_total),
                minutes(function_time_total),
                seconds(function_time_total));

   std::fprintf(file,
                "wall clock time for predictor routine:   %luh, %lum, %.2fs\n",
                hours(predictor_time_total),
                minutes(predictor_time_total),
                seconds(predictor_time_total));

   std::fprintf(file,
                "These routines combined took %.2f percent of the total computation\n",
                100.0 * main_routines / total_time);

   std::fprintf(file,
                "total wall clock time: %luh, %lum, %.2fs\n\n",
                hours(total_time),
                minutes(total_time),
                seconds(total_time));

   std::fclose(file);
}


static void TimesToScreen(double total_time) {
   double lsq_time_total = LsqSolveTimer::lsq_time_total;
   double jacobian_time_total = EvalJacobian::jacobian_time_total;
   double function_time_total = EvalFunction::function_time_total;
   double predictor_time_total = PredictorTimer::predictor_time_total;
   double main_routines = lsq_time_total + jacobian_time_total + function_time_total + predictor_time_total;

   std::printf("\nwall clock time for least squares routine in LeastSquaresNewton: %luh, %lum, %.2fs\n",
               hours(lsq_time_total),
               minutes(lsq_time_total),
               seconds(lsq_time_total));

   std::printf("wall clock time for GetJacobian routine: %luh, %lum, %.2fs\n",
               hours(jacobian_time_total),
               minutes(jacobian_time_total),
               seconds(jacobian_time_total));

   std::printf("wall clock time for GetFunction routine: %luh, %lum, %.2fs\n",
               hours(function_time_total),
               minutes(function_time_total),
               seconds(function_time_total));

   std::printf("wall clock time for predictor routine:   %luh, %lum, %.2fs\n",
               hours(predictor_time_total),
               minutes(predictor_time_total),
               seconds(predictor_time_total));

   std::printf("These routines combined took %.2f percent of the total computation\n",
               100.0 * main_routines / total_time);

   std::printf("total wall clock time: %luh, %lum, %.2fs\n",
               hours(total_time),
               minutes(total_time),
               seconds(total_time));
}


static void EfficiencyToFile(const std::pair<StdVector<double>, StdVector<std::string>>& data,
                             const std::string& dims_file_name) {
   std::string eff_path = std::string("results/efficiencies/" + dims_file_name);
   std::ofstream ofs{eff_path};
   if(!ofs) {
      GQ_THROW_RUNTIME_ERROR_MSG("could not open file for writing efficiency indexes");
   }

   const auto& eff_indexes = data.first;
   const auto& quad_names = data.second;

   ofs << std::endl;
   for(gq_int i = 0; i < eff_indexes.size(); ++i) {
      ofs << "efficiency for " << quad_names[i] << " = " << eff_indexes[i] << std::endl;
   }

   double ave_eff = std::accumulate(eff_indexes.begin(), eff_indexes.end(), 0.) / eff_indexes.size();
   ofs << std::endl;
   ofs << "average efficiency index = " << ave_eff << std::endl << std::endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
static void PrintDriverInfo(gq_int deg, gq_int dim, SearchWidth search_width, std::string domain_name) {
   std::cout << "\n\n";
   PrintDebugAndOmpInfo();

   std::cout << "computing quadrature rule for " << domain_name << " with deg = " << deg << " dim = " << dim
             << std::endl;
   if(domain_name != "interval") {
      std::cout << "search_width = " << search_width.width << std::endl;
   }
}


void PrintDebugAndOmpInfo() {
   static bool printed{false};
   if(!printed) {
#ifdef GEN_QUAD_DEBUG_ON
      std::printf("DEBUG MODE ON\n\n");
#else
      std::printf("DEBUG MODE OFF\n");
#endif

#ifdef _OPENMP
      std::printf("Running OpenMP with %i threads\n", omp_get_max_threads());
#else
      std::printf("OpenMP off\n");
#endif

#ifdef GEN_QUAD_USE_MKL
      std::printf("MKL enabled\n\n");
#endif

      printed = true;
   }
}


static gq_int hours(double x) {
   return gq_int(x / 3600);
}


static gq_int minutes(double x) {
   return gq_int(x / 60) % 60;
}


static double seconds(double x) {
   return x - gq_int(x / 60) * 60;
}


static void ResetTimers() {
   EvalFunction::function_time_total = 0;
   EvalJacobian::jacobian_time_total = 0;
   PredictorTimer::predictor_time_total = 0;
   LsqSolveTimer::lsq_time_total = 0;
}


}  // namespace gquad

#endif

