#pragma once

#if defined(_WIN32) || defined(__CYGWIN__)
#if defined(SURFACE_STATIC)
#define SURFACE_API
#elif defined(SURFACE_BUILD_DYNAMIC)
#define SURFACE_API __declspec(dllexport)
#else
#define SURFACE_API __declspec(dllimport)
#endif
#else
#if defined(SURFACE_STATIC)
#define SURFACE_API
#elif defined(SURFACE_BUILD_DYNAMIC)
#define SURFACE_API __attribute__((visibility("default")))
#else
#define SURFACE_API
#endif
#endif

#ifdef __cplusplus
extern "C"
{
#endif
    void SURFACE_API minimum_curvature(double *x, double *y, double *z, size_t n,
                                       double xmin, double xmax, double ymin, double ymax,
                                       size_t nx, size_t ny, double *output,
                                       unsigned int downsample_mode,
                                       size_t max_iterations,
                                       double relax,
                                       double alpha,
                                       double b_tension, double i_tension,
                                       double converge_limit,
                                       unsigned int converge_mode,
                                       double *lower,
                                       size_t n_lower,
                                       double *upper,
                                       size_t n_upper,
                                       unsigned char verbosity);

#ifdef __cplusplus
}
#endif