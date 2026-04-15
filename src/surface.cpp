/*
 * This file was modified from the original surface.cpp file of GMT tools by CGEM
 * it was moderatly rewritten to use C++ semantics and stdlib features, to be
 * a standalone function, and callable without using the GMT API.
 * It is also distributed under the GNU Lesser Public License; version 3
 *
 * This cpp version is Copyright (c) 2026 by CGEM
 *
 * The original license statement is below!
 */
/*--------------------------------------------------------------------
 *
 *	Copyright (c) 1991-2026 by the GMT Team (https://www.generic-mapping-tools.org/team.html)
 *	See LICENSE.TXT file for copying and redistribution conditions.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU Lesser General Public License as published by
 *	the Free Software Foundation; version 3 or any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU Lesser General Public License for more details.
 *
 *	Contact info: www.generic-mapping-tools.org
 *--------------------------------------------------------------------*/
/*
 * surface.cpp: a gridding program using splines in tension.
 * accepts xyz Cartesian values and fits a surface to the data.
 * The surface satisfies (1 - T) D4 z - T D2 z = 0,
 * where D4 is the 2-D biharmonic operator, D2 is the
 * 2-D Laplacian, and T is a "tension factor" between 0 and 1.
 * End member T = 0 is the classical minimum curvature
 * surface.  T = 1 gives a harmonic surface.  Use T = 0.25
 * or so for potential data; something more for topography.
 *
 * Interface includes over-relaxation for fast convergence and
 * automatic optimal grid factorization.
 *
 * See reference Smith & Wessel (Geophysics, 3, 293-305, 1990) for details.
 *
 * Authors:	Walter H. F. Smith and Paul Wessel
 * Date:	1-JAN-2010
 * Version:	6 API
 *
 * This version is based on the 6.0 rewrite that:
 * 1. uses scan-line grid structures, so we no longer need to transpose grids
 * 2. keeps node spacing at 1, thus we no longer need complicated strides between
 *    active nodes.  That spacing is now always 1 and we expand the grid as we
 *    go to larger grids (i.e., adding more nodes).
 */

#include <iostream>
#include <algorithm>
#include <vector>
#include <cstddef>
#include <cmath>
#include <limits>
#include <format>
#include <span>
#include <string>
#include <functional>

#include "surface.h"

/* Various constants used in surface */
#define SURFACE_CONV8_LIMIT 1.0E-8
#define SURFACE_OUTSIDE std::numeric_limits<size_t>::max()         /* Index number indicating data is outside usable area */
#define SURFACE_USE_DATA_LIMITS std::numeric_limits<size_t>::max() /* Key to use when specifying n_lower/n_upper to use the data limits */
#define SURFACE_CONV_LIMIT 0.0001                                  /* Default is 100 ppm of data range as convergence criterion */
#define SURFACE_MAX_ITERATIONS 500                                 /* Default iterations at final grid size */
#define SURFACE_OVERRELAXATION 1.4                                 /* Default over-relaxation value */
#define SURFACE_CLOSENESS_FACTOR 0.05                              /* A node is considered known if the nearest data is within 0.05 of a gridspacing of the node */
#define SURFACE_BREAKLINE 1                                        /* Flag for breakline constraints that should overrule data constraints */

namespace surface
{
    enum class NodeStatus : unsigned int
    {
        IS_UNCONSTRAINED = 0, /* Node has no data constraint within its bin box */
        DATA_IS_IN_QUAD1 = 1, /* Nearnest data constraint is in quadrant 1 relative to current node */
        DATA_IS_IN_QUAD2 = 2, /* Nearnest data constraint is in quadrant 2 relative to current node */
        DATA_IS_IN_QUAD3 = 3, /* Nearnest data constraint is in quadrant 3 relative to current node */
        DATA_IS_IN_QUAD4 = 4, /* Nearnest data constraint is in quadrant 4 relative to current node */
        IS_CONSTRAINED = 5,   /* Node has already been set (either data constraint < 5% of grid size or during filling) */
    };

    enum NodeConstraintType : unsigned int
    {
        UNCONSTRAINED = 0,
        CONSTRAINED = 1,
    };

    enum Dimensions
    {
        GMT_X = 0, /* x or lon is in 0th column */
        GMT_Y = 1, /* y or lat is in 1st column */
        GMT_Z = 2  /* z is in 2nd column */
    };

    enum WESNIds
    {
        XLO = 0, /* Index for west or xmin value */
        XHI = 1, /* Index for east or xmax value */
        YLO = 2, /* Index for south or ymin value */
        YHI = 3, /* Index for north or ymax value */
        ZLO = 4, /* Index for zmin value */
        ZHI = 5  /* Index for zmax value */
    };

    struct GridHeader
    {
        /* Variables we document for the API:
         * == Do not change the type of the following three items.
         * == They are copied verbatim to the native grid header and must be 4-byte unsigned ints. */
        size_t n_columns; /* Number of columns */
        size_t n_rows;    /* Number of rows */

        /* -- Here is the possible location for data structure padding:
         *    A double is 8-byte aligned on Windows. */

        /* == The types of the following 12 elements must not be changed.
         * == They are also copied verbatim to the native grid header. */
        double wesn[4]; /* Min/max x and y coordinates */
        double inc[2];  /* x and y increment */

        /* Items not stored in the data file for grids but explicitly used in macros computing node numbers */
        size_t nm;           /* Number of data items in this grid (n_columns * n_rows) [padding is excluded] */
        size_t size;         /* Actual number of items (not bytes) required to hold this grid (= mx * my), per band (for images) */
        unsigned int mx, my; /* Actual dimensions of the grid in memory, allowing for the padding */
        unsigned int pad[4]; /* Padding on west, east, south, north sides [2,2,2,2] */
    };

    /* grd is stored in rows going from west (xmin) to east (xmax)
     * first row in file has yvalue = north (ymax).
     * This is SCANLINE orientation.*/

    /*-----------------------------------------------------------------------------------------
     *	Notes on registration:

        Assume x_min = y_min = 0 and x_max = y_max = 10 and x_inc = y_inc = 1.
        For a normal node grid we have:
            (1) n_columns = (x_max - x_min) / x_inc + 1 = 11
                n_rows = (y_max - y_min) / y_inc + 1 = 11
            (2) node # 0 is at (x,y) = (x_min, y_max) = (0,10) and represents the surface
                value in a box with dimensions (1,1) centered on the node.
        For a pixel grid we have:
            (1) n_columns = (x_max - x_min) / x_inc = 10
                n_rows = (y_max - y_min) / y_inc = 10
            (2) node # 0 is at (x,y) = (x_min + 0.5*x_inc, y_max - 0.5*y_inc) = (0.5, 9.5)
                and represents the surface value in a box with dimensions (1,1)
                centered on the node.
    -------------------------------------------------------------------------------------------*/

    struct Grid
    {                             /* To hold a GMT gmt_grdfloat grid and its header in one container */
        GridHeader header;        /* GMT header for the grid */
        std::vector<double> data; /* grid data */
    };

    template <typename T>
    inline T ijp(const GridHeader &h, const T row, const T col) { return (row + h.pad[YHI]) * h.mx + col + h.pad[XLO]; }

    /* Go from row, col to grid node location, accounting for the 2 boundary rows and columns: */
    template <typename I>
    inline I row_col_to_node(const I row, const I col, const I mx) { return (row + 2) * mx + col + 2; }

    /* Go from row, col to index array position, where no boundary rows and columns involved: */
    template <typename I>
    inline I row_col_to_index(const I row, const I col, const I n_columns) { return row * n_columns + col; }

    /* Go from data x to fractional column x: */
    template <typename T>
    inline T x_to_fcol(const T x, const T x0, const T idx) { return (x - x0) * idx; }

    /* Go from x to grid integer column knowing it is a gridline-registered grid: */
    template <typename T>
    inline size_t x_to_col(const T x, const T x0, const T idx) { return static_cast<size_t>(x_to_fcol(x, x0, idx) + 0.5); }

    /* Go from data y to fractional row y_up measured from south (y_up = 0) towards north (y_up = n_rows-1): */
    template <typename T>
    inline T y_to_frow(const T y, const T y0, const T idy) { return (y - y0) * idy; }

    /* Go from y to row (row = 0 is north) knowing it is a gridline-registered grid: */
    template <typename T>
    inline size_t y_to_row(const T y, const T y0, const T idy, const size_t n_rows) { return n_rows - 1 - x_to_col(y, y0, idy); }

    /* Go from col to x knowing it is a gridline-registered grid: */
    template <typename T, typename I>
    inline T col_to_x(const I col, const T x0, const T x1, const T dx, const I n_columns) { return (col + 1 == n_columns) ? x1 : x0 + col * dx; }

    /* Go from row to y knowing it is a gridline-registered grid: */
    template <typename T, typename I>
    inline T row_to_y(const I row, const T y0, const T y1, const T dy, const I n_rows) { return (row + 1 == n_rows) ? y0 : y1 - row * dy; }

    /* Extract col from index: */
    template <typename I>
    inline I index_to_col(const I index, const I n_columns) { return index % n_columns; }

    /* Extract row from index: */
    template <typename I>
    inline I index_to_row(const I index, const I n_columns) { return index / n_columns; }

    /* Add a ptrdiff_t to a size_t by casting i to ptrdiff_t first*/
    inline size_t offset_ind(const size_t i, const std::ptrdiff_t offset)
    {
        return static_cast<size_t>(static_cast<std::ptrdiff_t>(i) + offset);
    }

    enum OffsetDirection
    { /* Node locations relative to current node, using compass directions */
      N2 = 0,
      NW = 1,
      N1 = 2,
      NE = 3,
      W2 = 4,
      W1 = 5,
      E1 = 6,
      E2 = 7,
      SW = 8,
      S1 = 9,
      SE = 10,
      S2 = 11
    }; /* I.e., 0-11 */

    /* The 4 indices per quadrant refer to points A-D in Figure A-1 in the reference for quadrant 1 */

    static OffsetDirection p[5][4] = {
        /* Indices into C.offset for each of the 4 quadrants, i.e., C.offset[p[quadrant][k]]], k = 0-3 */
        {N2, N2, N2, N2}, /* This row is never used, so allows us to use quadrant 1-4 as first array index directly (N2=0) */
        {NW, W1, S1, SE}, /* Indices for 1st quadrant */
        {SW, S1, E1, NE}, /* Indices for 2nd quadrant */
        {SE, E1, N1, NW}, /* Indices for 3rd quadrant */
        {NE, N1, W1, SW}  /* Indices for 4th quadrant */
    };

    enum class ConvergenceLimitMode : unsigned int
    {
        BY_PERCENT = 0,
        BY_VALUE = 1,
    };

    enum class DownsampleMode : unsigned int
    {
        CLOSEST = 0,
        MEAN = 1,
        MEDIAN = 2,
        // MODE = 3,
    };

    enum class IterMode
    {
        GRID_NODES = 0,
        GRID_DATA = 1
    };

    struct Data
    { /* Data point and index to node it currently constrains  */
        double x, y, z;
        size_t index;
    };

    struct Briggs
    { /* Coefficients in Taylor series for Laplacian(z) a la I. C. Briggs (1974)  */
        double b[6];
    };

    struct SurfaceInfo
    { /* Control structure for surface setup and execution */
        unsigned char verbosity;
        size_t node_sw_corner;              /* Node index of southwest interior grid corner for current stride */
        size_t node_se_corner;              /* Node index of southeast interior grid corner for current stride */
        size_t node_nw_corner;              /* Node index of northwest interior grid corner for current stride */
        size_t node_ne_corner;              /* Node index of northeast interior grid corner for current stride */
        size_t n_empty;                     /* No of unconstrained nodes at initialization  */
        size_t nxny;                        /* Total number of grid nodes without boundaries  */
        size_t mxmy;                        /* Total number of grid nodes with padding */
        size_t total_iterations;            /* Total iterations so far. */
        std::vector<Data> points;           /* All the data constraints */
        std::span<const Data> points_sub;   /* Subset of data for current grid size*/
        DownsampleMode downsample_mode;     /* Method used to downsample data if multiple points fall in the same grid cell. */
        std::vector<Briggs> brgs_coefs;     /* Array with Briggs 6-coefficients per nearest active data constraint */
        std::vector<size_t> factors;        /* Array of these common factors */
        size_t max_iterations;              /* Max iterations per call to iterate */
        ConvergenceLimitMode converge_mode; /* BY_PERCENT if -C set fractional convergence limit [BY_VALUE] */
        size_t current_stride;              /* Current node spacings relative to final spacing  */
        size_t previous_stride;             /* Previous node spacings relative to final spacing  */
        size_t n_columns;                   /* Number of nodes in x-dir. (Final grid) */
        size_t n_rows;                      /* Number of nodes in y-dir. (Final grid) */
        size_t mx;                          /* Width of final grid including padding */
        size_t my;                          /* Height of final grid including padding */
        size_t current_nx;                  /* Number of nodes in x-dir for current stride */
        size_t current_ny;                  /* Number of nodes in y-dir for current stride */
        size_t current_mx;                  /* Number of current nodes in x-dir plus 4 extra columns */
        size_t current_my;                  /* Number of current nodes in x-dir plus 4 extra columns */
        size_t previous_nx;                 /* Number of nodes in x-dir for previous stride */
        size_t previous_ny;                 /* Number of nodes in y-dir for previous stride */
        size_t previous_mx;                 /* Number of current nodes in x-dir plus 4 extra columns */
        size_t current_mxmy;                /* Total number of grid nodes with padding */
        std::ptrdiff_t offset[12];          /* Node-indices shifts of 12 nearby points relative center node */
        std::vector<NodeStatus> status;     /* Array with node status or quadrants */
        std::vector<double> lower_bound;    /* array with upper bounds */
        std::vector<double> upper_bound;    /* array with lower bounds */
        double inc[2];                      /* Size of each grid cell for current grid factor */
        double r_inc[2];                    /* Reciprocal grid spacings  */
        double converge_limit;              /* Convergence limit */
        double tension;                     /* Tension parameter on the surface  */
        double boundary_tension;            /* Tension parameter at the boundary */
        double interior_tension;            /* Tension parameter in the interior */
        double z_min;                       /* min data value (before de-trend)*/
        double z_max;                       /* max data value (before de-trend)*/
        double z_mean;                      /* Mean value of the data constraints z */
        double z_rms;                       /* Root mean square range of z after removing planar trend  */
        double r_z_rms;                     /* Reciprocal of z_rms (to avoid dividing) */
        double plane_icept;                 /* Intercept of best fitting plane to data  */
        double plane_sx;                    /* Slope of best fitting plane to data in x-direction */
        double plane_sy;                    /* Slope of best fitting plane to data in y-direction */
        double coeff[2][12];                /* Coefficients for 12 nearby nodes, for constrained [0] and unconstrained [1] nodes */
        double relax_old, relax_new;        /* Coefficients for relaxation factor to speed up convergence */
        double wesn[4];                     /* Original -R domain as we might have shifted it due to -r */
        double alpha;                       /* Aspect ratio dy/dx (1 for square pixels) */
        double a0_const_1, a0_const_2;      /* Various constants for off gridnode point equations */
        double alpha2, e_m2, one_plus_e2;
        double eps_p2, eps_m2, two_plus_ep2;
        double two_plus_em2;
    };

    /* Evaluate the change in LS plane from (0,0) to (xx,y_up) (so in intercept involved): */
    template <typename T>
    inline double evaluate_trend(const SurfaceInfo &C, const T xx, const T y_up) { return C.plane_sx * xx + C.plane_sy * y_up; }

    /* Evaluate the LS plane at location (xx,y_up) (this includes the intercept): */
    template <typename T>
    inline double evaluate_plane(const SurfaceInfo &C, const T xx, const T y_up) { return C.plane_icept + evaluate_trend(C, xx, y_up); }

    void set_coefficients(SurfaceInfo &C)
    {
        /* These are the coefficients in the finite-difference expressions given
         * by equations (A-4) [UNCONSTRAINED=0] and (A-7) [CONSTRAINED=1] in the reference.
         * Note that the UNCONSTRAINED coefficients are normalized by a0 (20 for no tension/aspects)
         * whereas the CONSTRAINED is used for a partial sum hence the normalization is done when the
         * sum over the Briggs coefficients have been included in iterate. */
        double alpha4, loose, a0;

        if (C.verbosity > 1)
        {
            std::cout << "Set finite-difference coefficients [stride = " << C.current_stride << "]\n"
                      << std::endl;
        }

        loose = 1.0 - C.interior_tension;
        C.alpha2 = C.alpha * C.alpha;
        alpha4 = C.alpha2 * C.alpha2;
        C.eps_p2 = C.alpha2;
        C.eps_m2 = 1.0 / C.alpha2;
        C.one_plus_e2 = 1.0 + C.alpha2;
        C.two_plus_ep2 = 2.0 + 2.0 * C.eps_p2;
        C.two_plus_em2 = 2.0 + 2.0 * C.eps_m2;

        C.e_m2 = 1.0 / C.alpha2;

        a0 = 1.0 / ((6 * alpha4 * loose + 10 * C.alpha2 * loose + 8 * loose - 2 * C.one_plus_e2) + 4 * C.interior_tension * C.one_plus_e2);
        C.a0_const_1 = 2.0 * loose * (1.0 + alpha4);
        C.a0_const_2 = 2.0 - C.interior_tension + 2 * loose * C.alpha2;

        C.coeff[CONSTRAINED][W2] = C.coeff[CONSTRAINED][E2] = -loose;
        C.coeff[CONSTRAINED][N2] = C.coeff[CONSTRAINED][S2] = -loose * alpha4;
        C.coeff[UNCONSTRAINED][W2] = C.coeff[UNCONSTRAINED][E2] = -loose * a0;
        C.coeff[UNCONSTRAINED][N2] = C.coeff[UNCONSTRAINED][S2] = -loose * alpha4 * a0;
        C.coeff[CONSTRAINED][W1] = C.coeff[CONSTRAINED][E1] = 2 * loose * C.one_plus_e2;
        C.coeff[UNCONSTRAINED][W1] = C.coeff[UNCONSTRAINED][E1] = (2 * C.coeff[CONSTRAINED][W1] + C.interior_tension) * a0;
        C.coeff[CONSTRAINED][N1] = C.coeff[CONSTRAINED][S1] = C.coeff[CONSTRAINED][W1] * C.alpha2;
        C.coeff[UNCONSTRAINED][N1] = C.coeff[UNCONSTRAINED][S1] = C.coeff[UNCONSTRAINED][W1] * C.alpha2;
        C.coeff[CONSTRAINED][NW] = C.coeff[CONSTRAINED][NE] = C.coeff[CONSTRAINED][SW] =
            C.coeff[CONSTRAINED][SE] = -2 * loose * C.alpha2;
        C.coeff[UNCONSTRAINED][NW] = C.coeff[UNCONSTRAINED][NE] = C.coeff[UNCONSTRAINED][SW] =
            C.coeff[UNCONSTRAINED][SE] = C.coeff[CONSTRAINED][NW] * a0;

        C.alpha2 *= 2; /* We will need these coefficients times two in the boundary conditions; do the doubling here  */
        C.e_m2 *= 2;
    }

    void set_offset(SurfaceInfo &C)
    {
        /* The offset array holds the offset in 1-D index relative
         * to the current node.  For movement along a row this is
         * always -2, -1, 0, +1, +2 but along a column we move in
         * multiples of current_mx, the extended grid row width,
         * which is current_mx = current_nx + 4.
         */
        C.offset[N2] = -2 * static_cast<std::ptrdiff_t>(C.current_mx); /* N2: 2 rows above */
        C.offset[NW] = -static_cast<std::ptrdiff_t>(C.current_mx) - 1; /* NW: 1 row above and one column left */
        C.offset[N1] = -static_cast<std::ptrdiff_t>(C.current_mx);     /* N1: 1 row above */
        C.offset[NE] = -static_cast<std::ptrdiff_t>(C.current_mx) + 1; /* NE: 1 row above and one column right */
        C.offset[W2] = -2;                                             /* W2: 2 columns left */
        C.offset[W1] = -1;                                             /* W1 : 1 column left */
        C.offset[E1] = +1;                                             /* E1 : 1 column right */
        C.offset[E2] = +2;                                             /* E2 : 2 columns right */
        C.offset[SW] = static_cast<std::ptrdiff_t>(C.current_mx) - 1;  /* SW : 1 row below and one column left */
        C.offset[S1] = static_cast<std::ptrdiff_t>(C.current_mx);      /* S1 : 1 row below */
        C.offset[SE] = static_cast<std::ptrdiff_t>(C.current_mx) + 1;  /* SE : 1 row below and one column right */
        C.offset[S2] = 2 * static_cast<std::ptrdiff_t>(C.current_mx);  /* S2 : 2 rows below */
    }

    void fill_in_forecast(SurfaceInfo &C, Grid &G)
    {

        /* Fills in bilinear estimates into new node locations after grid is expanded.
           These new nodes are marked as unconstrained while the coarser data are considered
           constraints in the next iteration.  We do this in two steps:
             a) We sweep through the grid from last to first node and copy each node to the
                new location due to increased grid dimensions.
             b) Once nodes are in place we sweep through and apply the bilinear interpolation.
         */

        size_t index_00, index_10, index_11, index_01, index_new, current_node, previous_node;
        size_t previous_row, previous_col, i, j, col, row, expand, first;
        double c, sx, sy, sxy, r_prev_size, c_plus_sy_dy, sx_plus_sxy_dy;
        auto status = std::span<NodeStatus>(C.status.begin(), C.current_mxmy);
        auto u = std::span<double>(G.data.begin(), C.current_mxmy);
        std::vector<double> fraction;

        /* First we expand the active grid to allow for more nodes. We do this by
         * looping backwards from last node to first so that the locations we copy
         * the old node values to will always have higher node number than their source.
         * The previous grid solution has dimensions previous_nx x previous_ny while
         * the new grid has dimensions current_nx x current_ny.  We thus loop over
         * the old grid and place these nodes into the new grid.  */

        expand = C.previous_stride / C.current_stride; /* Multiplicity of new nodes in both x and y dimensions */

        if (C.verbosity > 1)
        {
            std::cout << "Expand grid by factor of " << expand << " when going from stride = " << C.previous_stride << " to " << C.current_stride << "\n"
                      << std::endl;
        }

        for (previous_row = C.previous_ny; previous_row-- > 0;)
        {                                /* Loop backward over the previous grid rows */
            row = previous_row * expand; /* Corresponding row in the new extended grid */
            for (previous_col = C.previous_nx; previous_col-- > 0;)
            {                                                                               /* Loop backward over previous grid cols */
                col = previous_col * expand;                                                /* Corresponding col in the new extended grid */
                current_node = row_col_to_node(row, col, C.current_mx);                     /* Current node index */
                previous_node = row_col_to_node(previous_row, previous_col, C.previous_mx); /* Previous node index */
                u[current_node] = u[previous_node];                                         /* Copy the value over */
            }
        }

        /* The active grid has now increased in size and the previous values have been copied to their new nodes.
         * The grid nodes in-between these new "constrained" nodes are partly filled with old values (since
         * we just copied, not moved, the nodes) or zeros (since we expanded the grid into new unused memory).
         * This does not matter since we will now fill in those in-between nodes with a bilinear interpolation
         * based on the coarser (previous) nodes.  At the end all nodes in the active grid are valid, except
         * in the boundary rows/cols.  These are reset by set_BC before the iteration starts.
         */

        /* Precalculate the fractional increments of rows and cols in-between the old constrained rows and cols.
         * These are all fractions between 0 and 1.  E.g., if we quadruple the grid dimensions in x and y then
         * expand == 4 and we need 4 fractions = {0, 0.25, 0.5, 0.75}. */

        fraction.resize(expand);
        r_prev_size = 1.0 / C.previous_stride;
        for (i = 0; i < expand; i++)
            fraction[i] = i * r_prev_size;

        if (C.verbosity > 1)
        {
            std::cout << "Fill in expanded grid by bilinear interpolation [stride = " << C.current_stride << "]\n"
                      << std::endl;
        }

        /* Loop over 4-point "bin squares" from the first northwest bin to the last southeast bin. The bin vertices are the expanded previous nodes */

        for (previous_row = 1; previous_row < C.previous_ny; previous_row++)
        {                                /* Starts at row 1 since it is the baseline for the bin extending up to row 0 (north) */
            row = previous_row * expand; /* Corresponding row in the new extended grid */

            for (previous_col = 0; previous_col < (C.previous_nx - 1); previous_col++)
            {                                /* Stop 1 column short of east since east is the right boundary of last bin */
                col = previous_col * expand; /* Corresponding col in the new extended grid */

                /* Get the indices of the bilinear square defined by nodes {00, 10, 11, 01}, with 00 referring to the current (lower left) node */
                index_00 = row_col_to_node(row, col, C.current_mx); /* Lower left corner of square bin and our origin */
                index_01 = index_00 - expand * C.current_mx;        /* Upper left corner of square bin */
                index_10 = index_00 + expand;                       /* Lower right corner of square bin */
                index_11 = index_01 + expand;                       /* Upper right corner of square bin */

                /* Get bilinear coefficients for interpolation z = c + sx * delta_x + sy * delta_y + sxy * delta_x * delta_y,
                 * which we will use as z = (c + sy * delta_y) + delta_x * (sx + sxy * delta_y).
                 * Below, delta_x and delta_y are obtained via C.fraction[i|j] that we pre-calculated above. */
                c = u[index_00];
                sx = u[index_10] - c;
                sy = u[index_01] - c;
                sxy = u[index_11] - u[index_10] - sy;

                /* Fill in all the denser nodes except the lower-left starting point */

                for (j = 0, first = 1; j < expand; j++)
                {                                        /* Set first = 1 so we skip the first column when j = 0 */
                    c_plus_sy_dy = c + sy * fraction[j]; /* Compute terms that remain constant for this j */
                    sx_plus_sxy_dy = sx + sxy * fraction[j];
                    index_new = index_00 - j * C.current_mx + first; /* Start node on this intermediate row */
                    for (i = first; i < expand; i++, index_new++)
                    { /* Sweep across this row and interpolate */
                        u[index_new] = (c_plus_sy_dy + fraction[i] * sx_plus_sxy_dy);
                        status[index_new] = NodeStatus::IS_UNCONSTRAINED; /* These are considered temporary estimates */
                    }
                    first = 0; /* Reset to 0 for the remainder of the j loop */
                }
                status[index_00] = NodeStatus::IS_CONSTRAINED; /* The previous node values will be kept fixed in the next iterate call */
            }
        }

        /* The loops above exclude the north and east boundaries.  First do linear interpolation along the east edge */
        index_00 = C.node_ne_corner; /* Upper NE node */
        for (previous_row = 1; previous_row < C.previous_ny; previous_row++)
        {                                        /* So first edge is from row = 1 up to row = 0 on eastern edge */
            index_01 = index_00;                 /* Previous lower becomes current upper node */
            index_00 += expand * C.current_mx;   /* Lower node after striding down */
            sy = u[index_01] - u[index_00];      /* Vertical gradient in u toward ymax (for increasing j) */
            index_new = index_00 - C.current_mx; /* Since we start at j = 1 we skip up one row here */
            for (j = 1; j < expand; j++, index_new -= C.current_mx)
            { /* Start at 1 since we skip the constrained index_00 node */
                u[index_new] = u[index_00] + fraction[j] * sy;
                status[index_new] = NodeStatus::IS_UNCONSTRAINED; /* These are considered temporary estimates */
            }
            status[index_00] = NodeStatus::IS_CONSTRAINED; /* The previous node values will be kept fixed in the next iterate call */
        }

        /* Next do linear interpolation along the north edge */
        index_10 = C.node_nw_corner; /* Left NW node */
        for (previous_col = 0; previous_col < (C.previous_nx - 1); previous_col++)
        {                                   /* To ensure last edge ends at col = C.previous_nx-1 */
            index_00 = index_10;            /* Previous right node becomes current left node */
            index_10 = index_00 + expand;   /* Right node after striding to the right */
            sx = u[index_10] - u[index_00]; /* Horizontal gradient in u toward xmax (for increasing i) */
            index_new = index_00 + 1;       /* Start at 1 since we skip the constrained index_00 node */

            for (i = 1; i < expand; i++, index_new++)
            {
                u[index_new] = u[index_00] + fraction[i] * sx;
                status[index_new] = NodeStatus::IS_UNCONSTRAINED; /* These are considered temporary estimates */
            }
            status[index_00] = NodeStatus::IS_CONSTRAINED; /* The previous node values will be kept fixed in the next iterate call */
        }
        /* Finally set the northeast corner to be considered fixed in the next iterate call and our work here is done */
        status[C.node_ne_corner] = NodeStatus::IS_CONSTRAINED;
    }

    /* thunk arg is last argument to compare function */
    struct compare_points
    {
        SurfaceInfo *info;

        bool operator()(const Data &point_1, const Data &point_2) const
        {
            /* Routine for QSORT_R to sort data structure for fast access to data by node location.
                   Sorts on index first, then on radius to node corresponding to index, so that index
                   goes from low to high, and so does radius.  Note: These are simple Cartesian distance
                 * calculations.  The metadata needed to do the calculations are passed via *arg.
                */
            size_t col, row, index_1, index_2;
            double x0, y0, dist_1, dist_2;
            index_1 = point_1.index;
            index_2 = point_2.index;
            if (index_1 != index_2)
                return index_1 < index_2;
            if (index_1 == SURFACE_OUTSIDE)
                return false;
            /* Points are in same grid cell.  First check for breakline points to sort those ahead of data points */
            // if (point_1.kind == SURFACE_BREAKLINE && point_2.kind == 0)
            //     return true;
            // if (point_2.kind == SURFACE_BREAKLINE && point_1.kind == 0)
            //     return false;
            /* Now find the one who is nearest to grid point */
            /* Note: index calculations do not include boundary pad */
            // info = arg; /* Get the needed metadata for distance calculations */
            row = index_to_row(index_1, info->current_nx);
            col = index_to_col(index_1, info->current_nx);
            x0 = col_to_x(col, info->wesn[XLO], info->wesn[XHI], info->inc[GMT_X], info->current_nx);
            y0 = row_to_y(row, info->wesn[YLO], info->wesn[YHI], info->inc[GMT_Y], info->current_ny);
            dist_1 = (point_1.x - x0) * (point_1.x - x0) + (point_1.y - y0) * (point_1.y - y0);
            /* Try to speed things up by first checking if point_2 x-distance from x0 alone exceeds point_1's radial distance */
            dist_2 = (point_2.x - x0) * (point_2.x - x0); /* Just dx^2 */
            if (dist_1 < dist_2)
                return true; /* Don't need to consider the y-distance */
            /* Did not exceed, so now we must finalize the dist_2 calculation by including the y-separation */
            dist_2 += (point_2.y - y0) * (point_2.y - y0);
            return dist_1 < dist_2;
        }
    };

    void smart_divide(SurfaceInfo &C)
    {
        /* Divide grid by its next largest prime factor and shift that setting by one */
        C.current_stride /= C.factors.back();
        C.factors.pop_back();
    }

    void set_index(SurfaceInfo &C, const GridHeader &h)
    {
        /* Recomputes data[k].index for the new value of the stride,
           sorts the data again on index and radii, and throws away
           data which are now outside the usable limits.
           Note: These indices exclude the padding. */
        size_t col, row;
        size_t k_skipped = 0;

        if (C.verbosity > 1)
        {
            std::cout << "Recompute data index for next iteration [stride = " << C.current_stride << "]\n"
                      << std::endl;
        }

        for (auto &pt : C.points)
        {
            col = x_to_col(pt.x, h.wesn[XLO], C.r_inc[GMT_X]);
            row = y_to_row(pt.y, h.wesn[YLO], C.r_inc[GMT_Y], C.current_ny);
            if (col < 0 || col >= C.current_nx || row < 0 || row >= C.current_ny)
            {
                pt.index = SURFACE_OUTSIDE;
                k_skipped++;
            }
            else
                pt.index = row_col_to_index(row, col, C.current_nx);
        }

        // sort the SURFACE_OUTSIDE points to the end
        std::sort(C.points.begin(), C.points.end(), compare_points{&C});
        C.points_sub = std::span<const Data>(C.points.begin(), C.points.size() - k_skipped);
    }

    void solve_Briggs_coefficients(SurfaceInfo &C, double (&b)[6], double xx, double yy, double z)
    {
        /* Given the normalized offset (xx,yy) from current node (value z) we determine the
         * Briggs coefficients b_k, k = 1,5  [Equation (A-6) in the reference]
         * Here, xx, yy are the fractional distances, accounting for any anisotropy.
         * Note b[5] initially contains the sum of the 5 Briggs coefficients but
         * we actually need to divide by it so we do that change here as well.
         * Finally, b[4] will be multiplied with the off-node constraint so we do that here.
         */
        double xx2, yy2, xx_plus_yy, xx_plus_yy_plus_one, inv_xx_plus_yy_plus_one, inv_delta, b_4;

        xx_plus_yy = xx + yy;
        xx_plus_yy_plus_one = 1.0 + xx_plus_yy;
        inv_xx_plus_yy_plus_one = 1.0 / xx_plus_yy_plus_one;
        xx2 = xx * xx;
        yy2 = yy * yy;
        inv_delta = inv_xx_plus_yy_plus_one / xx_plus_yy;
        b[0] = (xx2 + 2.0 * xx * yy + xx - yy2 - yy) * inv_delta;
        b[1] = 2.0 * (yy - xx + 1.0) * inv_xx_plus_yy_plus_one;
        b[2] = 2.0 * (xx - yy + 1.0) * inv_xx_plus_yy_plus_one;
        b[3] = (-xx2 + 2.0 * xx * yy - xx + yy2 + yy) * inv_delta;
        b_4 = 4.0 * inv_delta;
        /* We also need to normalize by the sum of the b[k] values, so sum them here */
        b[5] = b[0] + b[1] + b[2] + b[3] + b_4;
        /* We need to sum k = 0<5 of u[k]*b[k], where u[k] are the nodes of the points A-D,
         * but the k = 4 point (E) is our data constraint.  We multiply that in here, once,
         * add add b[4] to the rest of the sum inside the iteration loop. */
        b[4] = (b_4 * z);

        /* b[5] is part of a denominator so we do the division here instead of inside iterate loop */
        b[5] = (1.0 / (C.a0_const_1 + C.a0_const_2 * b[5]));
    }

    void find_nearest_constraint(SurfaceInfo &C, Grid &G)
    {
        /* Determines the nearest data point per bin and sets the
         * Briggs parameters or, if really close, fixes the node value */
        size_t last_index, node, briggs_index;
        size_t row, col;
        double xx, yy, x0, y0, dx, dy;
        double z_at_node;
        std::vector<double> &u = G.data;
        std::vector<NodeStatus> &status = C.status;
        GridHeader &h = G.header;
        bool is_lower_bounded = C.lower_bound.size() > 0;
        bool is_upper_bounded = C.upper_bound.size() > 0;
        bool limited = is_lower_bounded || is_upper_bounded;

        if (C.verbosity > 1)
        {
            std::cout << "Determine nearest point and set Briggs coefficients [stride = " << C.current_stride << "]"
                      << std::endl;
        }

        // gmt_M_grd_loop(GMT, C.Grid, row, col, node)
        for (row = 0; row < h.n_rows; row++)
            for (col = 0, node = ijp(h, row, (size_t)0); col < h.n_columns; col++, node++)
            { /* Reset status of all interior grid nodes */
                status[node] = NodeStatus::IS_UNCONSTRAINED;
            }

        last_index = std::numeric_limits<size_t>::max(); /* Set to max so that first point is always different and we compute the nearest node for it */
        briggs_index = 0;

        for (auto &pt : C.points_sub)
        { /* Find constraining value  */
            if (pt.index != last_index)
            { /* Moving to the next node to address its nearest data constraint */
                /* Note: Index calculations do not consider the boundary padding */
                row = index_to_row(pt.index, C.current_nx);
                col = index_to_col(pt.index, C.current_nx);
                last_index = pt.index; /* Now this is the last unique index we worked on */
                node = row_col_to_node(row, col, C.current_mx);
                /* Get coordinates of this node */
                x0 = col_to_x(col, h.wesn[XLO], h.wesn[XHI], C.inc[GMT_X], C.current_nx);
                y0 = row_to_y(row, h.wesn[YLO], h.wesn[YHI], C.inc[GMT_Y], C.current_ny);
                /* Get offsets dx,dy of data point location relative to this node (dy is positive up) in fraction of grid increments */
                dx = x_to_fcol(pt.x, x0, C.r_inc[GMT_X]);
                dy = y_to_frow(pt.y, y0, C.r_inc[GMT_Y]);
                /* So dx, dy are here fractions of C.inc[GMT_X] and C-inc[GMT_Y] */
                /* "Really close" will mean within 5% of the current grid spacing from the center node */

                if (std::abs(dx) < SURFACE_CLOSENESS_FACTOR && std::abs(dy) < SURFACE_CLOSENESS_FACTOR)
                { /* Considered close enough to assign fixed value to node */
                    status[node] = NodeStatus::IS_CONSTRAINED;
                    /* Since data constraint is forcibly moved from (dx, dy) to (0,0) we must adjust for
                     * the small change in the planar trend between the two locations, and then
                     * possibly clip the value if constraining surfaces were given.  Note that
                     * dx, dy is in -1/1 range normalized by (current_x|y_inc) so to recover the
                     * corresponding dx,dy in units of current grid fractions we must scale both
                     * dx and dy by current_stride; this is equivalent to scaling the trend.
                     * This trend then is normalized by dividing by the z rms.*/

                    z_at_node = pt.z + C.r_z_rms * C.current_stride * evaluate_trend(C, dx, dy);
                    if (limited)
                    { /* Must use final spacing node index to access the Bound grids */
                        size_t node_final = ijp(G.header, C.current_stride * row, C.current_stride * col);
                        if (is_lower_bounded && !std::isnan(C.lower_bound[node_final]))
                            z_at_node = std::max(C.lower_bound[node_final], z_at_node);
                        if (is_upper_bounded && !std::isnan(C.upper_bound[node_final]))
                            z_at_node = std::min(C.upper_bound[node_final], z_at_node);
                    }
                    u[node] = z_at_node;
                }
                else
                { /* We have a nearby data point in one of the quadrants */
                    /* Note: We must swap dx,dy for 2nd and 4th quadrants and always use absolute values since we are
                       rotating other cases (quadrants 2-4) to look like quadrant 1 */
                    if (dy >= 0.0)
                    { /* Upper two quadrants */
                        if (dx >= 0.0)
                        { /* Both positive, use as is */
                            status[node] = NodeStatus::DATA_IS_IN_QUAD1;
                            xx = dx;
                            yy = dy;
                        }
                        else
                        { /* dx negative, so remove sign, and swap */
                            status[node] = NodeStatus::DATA_IS_IN_QUAD2;
                            yy = -dx;
                            xx = dy;
                        }
                    }
                    else
                    { /* Lower two quadrants where we need to remove sign from dy */
                        if (dx >= 0.0)
                        { /* Also swap x and y */
                            status[node] = NodeStatus::DATA_IS_IN_QUAD4;
                            yy = dx;
                            xx = -dy;
                        }
                        else
                        { /* Just remove both signs */
                            status[node] = NodeStatus::DATA_IS_IN_QUAD3;
                            xx = -dx;
                            yy = -dy;
                        }
                    }
                    auto &b = C.brgs_coefs[briggs_index].b; /* Reference to the Briggs coefficients for this data point */
                    /* Evaluate the Briggs coefficients */
                    solve_Briggs_coefficients(C, b, xx, yy, pt.z);
                    briggs_index++;
                }
            }
        }
    }

    void set_grid_parameters(SurfaceInfo &C, const GridHeader &h)
    {
        /* Set the previous settings to the current settings */
        C.previous_nx = C.current_nx;
        C.previous_mx = C.current_mx;
        C.previous_ny = C.current_ny;
        /* Update the current parameters given the new C.current_stride setting */
        C.current_nx = (C.n_columns - 1) / C.current_stride + 1;
        C.current_ny = (C.n_rows - 1) / C.current_stride + 1;
        C.current_mx = C.current_nx + 4;
        C.current_my = C.current_ny + 4;
        C.current_mxmy = C.current_mx * (C.current_ny + 4); /* Only place where "my" is used */
        C.inc[GMT_X] = C.current_stride * h.inc[GMT_X];
        C.inc[GMT_Y] = C.current_stride * h.inc[GMT_Y];
        C.r_inc[GMT_X] = 1.0 / C.inc[GMT_X];
        C.r_inc[GMT_Y] = 1.0 / C.inc[GMT_Y];
        /* Update the grid node indices of the 4 corners */
        C.node_nw_corner = 2 * C.current_mx + 2;
        C.node_sw_corner = C.node_nw_corner + (C.current_ny - 1) * C.current_mx;
        C.node_se_corner = C.node_sw_corner + C.current_nx - 1;
        C.node_ne_corner = C.node_nw_corner + C.current_nx - 1;
    }

    void setup_data(SurfaceInfo &C, const GridHeader &h, double *xs, double *ys, double *zs, size_t n)
    {
        /* Process input data into data structure */
        size_t col, row;
        double zmin = std::numeric_limits<double>::max();
        double zmax = std::numeric_limits<double>::min();
        double wesn_lim[4];
        // struct GMT_RECORD *In = NULL;

        // GMT_Report(GMT->parent, GMT_MSG_INFORMATION, "Processing input table data\n");
        C.points.reserve(n);
        // C.points = gmt_M_memory(GMT, NULL, C.n_alloc, struct SURFACE_DATA);

        /* Read in xyz data and computes index no and store it in a structure */

        // if ((error = GMT_Set_Columns(GMT->parent, GMT_IN, 3, GMT_COL_FIX_NO_TEXT)) != GMT_NOERROR)
        //     return (error);
        // if (GMT_Init_IO(GMT->parent, GMT_IS_DATASET, GMT_IS_POINT, GMT_IN, GMT_ADD_DEFAULT, 0, options) != GMT_NOERROR) /* Establishes data input */
        //     return (GMT->parent->error);

        C.z_mean = 0.0;
        /* Initially allow points to be within 1 grid spacing of the grid */
        wesn_lim[XLO] = h.wesn[XLO] - C.inc[GMT_X];
        wesn_lim[XHI] = h.wesn[XHI] + C.inc[GMT_X];
        wesn_lim[YLO] = h.wesn[YLO] - C.inc[GMT_Y];
        wesn_lim[YHI] = h.wesn[YHI] + C.inc[GMT_Y];

        // if (GMT_Begin_IO(GMT->parent, GMT_IS_DATASET, GMT_IN, GMT_HEADER_ON) != GMT_NOERROR) /* Enables data input and sets access mode */
        //     return (GMT->parent->error);

        for (size_t i = 0; i < n; ++i)
        {
            double x = xs[i];
            double y = ys[i];
            double z = zs[i];
            if (std::isnan(z))
                continue;
            if (std::isnan(y) || y < wesn_lim[YLO] || y > wesn_lim[YHI])
                continue;
            if (std::isnan(x) || x < wesn_lim[XLO] || x > wesn_lim[XHI])
                continue;
            row = y_to_row(y, h.wesn[YLO], C.r_inc[GMT_Y], C.current_ny);
            if (row < 0 || row >= C.current_ny)
                continue;
            else /* Regular point not at the periodic boundary */
                col = x_to_col(x, h.wesn[XLO], C.r_inc[GMT_X]);
            if (col < 0 || col >= C.current_nx)
                continue;

            Data pt{x, y, z, row_col_to_index(row, col, C.current_nx)};
            C.points.push_back(pt);
            // std::cout << "data point in: " << x << ", " << y << ", " << z << std::endl;
            /* Determine the mean, min and max z-values */
            if (zmin > z)
                zmin = z;
            if (zmax < z)
                zmax = z;
            C.z_mean += z;
        }

        if (C.points.size() == 0)
        {
            std::cerr << "No datapoints inside region, aborting" << std::endl;
            std::exit(1);
        }

        C.z_min = zmin;
        C.z_max = zmax;
        C.z_mean /= C.points.size(); /* Estimate mean data value */
        C.points.shrink_to_fit();
    }

    void surface_set_BCs(const SurfaceInfo &C, std::span<double> u)
    {
        /* Fill in auxiliary boundary rows and columns; see equations (A-8,9,10) in the reference */
        size_t n, n_s, n_n, n_w, n_e; /* Node variables */
        size_t col, row;
        auto &d_n = C.offset; /* Relative changes in node index from present node n */
        double x_0_const = 4.0 * (1.0 - C.boundary_tension) / (2.0 - C.boundary_tension);
        double x_1_const = (3 * C.boundary_tension - 2.0) / (2.0 - C.boundary_tension);
        double y_denom = 2 * C.alpha * (1.0 - C.boundary_tension) + C.boundary_tension;
        double y_0_const = 4 * C.alpha * (1.0 - C.boundary_tension) / y_denom;
        double y_1_const = (C.boundary_tension - 2 * C.alpha * (1.0 - C.boundary_tension)) / y_denom;

        if (C.verbosity > 2)
        {
            std::cout << "Apply boundary conditions to auxiliary nodes [stride = " << C.current_stride << "]" << std::endl;
        }

        /* First set (1-T)d2[]/dn2 + Td[]/dn = 0 along edges */

        for (col = 0, n_s = C.node_sw_corner, n_n = C.node_nw_corner; col < C.current_nx; col++, n_s++, n_n++)
        {
            /* set BC1 along south and north side */
            u[offset_ind(n_s, d_n[S1])] = (y_0_const * u[n_s] + y_1_const * u[offset_ind(n_s, d_n[N1])]); /* South: u_{0,-1} = 2 * u_{0,0} - u_{0,+1} */
            u[offset_ind(n_n, d_n[N1])] = (y_0_const * u[n_n] + y_1_const * u[offset_ind(n_n, d_n[S1])]); /* North: u_{0,+1} = 2 * u_{0,0} - u_{0,-1} */
        }
        for (row = 0, n_w = C.node_nw_corner, n_e = C.node_ne_corner; row < C.current_ny; row++, n_w += C.current_mx, n_e += C.current_mx)
        {
            /* West: u_{-10} = 2 * u_{00} - u_{10}  */
            u[offset_ind(n_w, d_n[W1])] = (x_1_const * u[offset_ind(n_w, d_n[E1])] + x_0_const * u[n_w]);
            /* East: u_{10} = 2 * u_{00} - u_{-10}  */
            u[offset_ind(n_e, d_n[E1])] = (x_1_const * u[offset_ind(n_e, d_n[W1])] + x_0_const * u[n_e]);
        }

        /* Now set d2[]/dxdy = 0 at each of the 4 corners */

        n = C.node_sw_corner; /* Just use shorthand in each expression */
        u[offset_ind(n, d_n[SW])] = u[offset_ind(n, d_n[SE])] + u[offset_ind(n, d_n[NW])] - u[offset_ind(n, d_n[NE])];
        n = C.node_nw_corner;
        u[offset_ind(n, d_n[NW])] = u[offset_ind(n, d_n[NE])] + u[offset_ind(n, d_n[SW])] - u[offset_ind(n, d_n[SE])];
        n = C.node_se_corner;
        u[offset_ind(n, d_n[SE])] = u[offset_ind(n, d_n[SW])] + u[offset_ind(n, d_n[NE])] - u[offset_ind(n, d_n[NW])];
        n = C.node_ne_corner;
        u[offset_ind(n, d_n[NE])] = u[offset_ind(n, d_n[NW])] + u[offset_ind(n, d_n[SE])] - u[offset_ind(n, d_n[SW])];

        /* Now set dC/dn = 0 at each edge */

        for (col = 0, n_s = C.node_sw_corner, n_n = C.node_nw_corner; col < C.current_nx; col++, n_s++, n_n++)
        { /* set BC2 along south and north side */
            /* South side */
            u[offset_ind(n_s, d_n[S2])] = (u[offset_ind(n_s, d_n[N2])] + C.eps_m2 * (u[offset_ind(n_s, d_n[NW])] + u[offset_ind(n_s, d_n[NE])] - u[offset_ind(n_s, d_n[SW])] - u[offset_ind(n_s, d_n[SE])]) + C.two_plus_em2 * (u[offset_ind(n_s, d_n[S1])] - u[offset_ind(n_s, d_n[N1])]));
            /* North side */
            u[offset_ind(n_n, d_n[N2])] = (u[offset_ind(n_n, d_n[S2])] + C.eps_m2 * (u[offset_ind(n_n, d_n[SW])] + u[offset_ind(n_n, d_n[SE])] - u[offset_ind(n_n, d_n[NW])] - u[offset_ind(n_n, d_n[NE])]) + C.two_plus_em2 * (u[offset_ind(n_n, d_n[N1])] - u[offset_ind(n_n, d_n[S1])]));
        }

        for (row = 0, n_w = C.node_nw_corner, n_e = C.node_ne_corner; row < C.current_ny; row++, n_w += C.current_mx, n_e += C.current_mx)
        { /* set BC2 along west and east side */
            /* West side */
            u[offset_ind(n_w, d_n[W2])] = (u[offset_ind(n_w, d_n[E2])] + C.eps_p2 * (u[offset_ind(n_w, d_n[NE])] + u[offset_ind(n_w, d_n[SE])] - u[offset_ind(n_w, d_n[NW])] - u[offset_ind(n_w, d_n[SW])]) + C.two_plus_ep2 * (u[offset_ind(n_w, d_n[W1])] - u[offset_ind(n_w, d_n[E1])]));
            /* East side */
            u[offset_ind(n_e, d_n[E2])] = (u[offset_ind(n_e, d_n[W2])] + C.eps_p2 * (u[offset_ind(n_e, d_n[NW])] + u[offset_ind(n_e, d_n[SW])] - u[offset_ind(n_e, d_n[NE])] - u[offset_ind(n_e, d_n[SE])]) + C.two_plus_ep2 * (u[offset_ind(n_e, d_n[E1])] - u[offset_ind(n_e, d_n[W1])]));
        }
    }

    size_t relaxation_iterate(const SurfaceInfo &C, Grid &G, IterMode mode)
    {
        /* Main finite difference solver */
        size_t node, briggs_index, iteration_count = 0, node_final;
        unsigned int set;
        size_t current_max_iterations = C.max_iterations * C.current_stride;
        size_t col, row, k;
        auto &d_node = C.offset; /* Relative changes in node index from present node */
        // char *mode_name[2] = {"node", "data"};
        bool finished;
        double current_limit = C.converge_limit / C.current_stride;
        double u_change, max_u_change, max_z_change, sum_bk_uk, u_00;
        const GridHeader &h = G.header;

        bool is_lower_bounded = C.lower_bound.size() > 0;
        bool is_upper_bounded = C.upper_bound.size() > 0;
        bool limited = is_lower_bounded || is_upper_bounded;

        // const std::vector<unsigned char> &status = C.status; /* Quadrant or status information for each node */
        auto status = std::span<const NodeStatus>(C.status.begin(), C.current_mxmy); /* Quadrant or status information for each node */
        auto u = std::span<double>(G.data.begin(), C.current_mxmy);

        std::string mode_name = (mode == IterMode::GRID_NODES) ? "Smoothing Update" : "Surface Update";
        if (C.verbosity > 1)
        {
            std::cout << "Start iterations [mode = " << mode_name << " Max iterations = " << current_max_iterations << " stride = " << C.current_stride << "]" << std::endl;
        }

        // sprintf(C.format, "%%4ld\t%%c\t%%8" PRIu64 "\t%s\t%s\t%%10" PRIu64 "\n", GMT->current.setting.format_float_out, GMT->current.setting.format_float_out);
        // if (C.logging)
        //     fprintf(C.fp_log, "%c Grid size = %d Mode = %c Convergence limit = %g -Z%d\n",
        //             GMT->current.setting.io_seg_marker[GMT_OUT], C.current_stride, C.mode_type[mode], current_limit, C.current_stride);

        do
        {
            // copy the current u into u_old...
            // std::copy(u.begin(), u.end(), u_old.begin());

            surface_set_BCs(C, u); /* Set the boundary rows and columns */

            briggs_index = 0;    /* Reset the Briggs constraint table index  */
            max_u_change = -1.0; /* Ensure max_u_change is < 0 for starters */

            /* Now loop over all interior data nodes */

            for (row = 0; row < C.current_ny; row++)
            {                                                 /* Loop over rows */
                node = C.node_nw_corner + row * C.current_mx; /* Node at left side of this row */
                if (limited)
                    node_final = ijp(h, C.current_stride * row, size_t(0));
                for (col = 0; col < C.current_nx; col++, node++)
                { /* Loop over all columns */
                    if (status[node] == NodeStatus::IS_CONSTRAINED)
                    { /* Data constraint fell exactly on the node, keep it as is */
                        continue;
                    }

                    /* Here we must estimate a solution via equations (A-4) [UNCONSTRAINED] or (A-7) [CONSTRAINED] */
                    u_00 = 0.0;                                                                         /* Start with zero, build updated solution for central node */
                    set = (status[node] == NodeStatus::IS_UNCONSTRAINED) ? UNCONSTRAINED : CONSTRAINED; /* Index to C.coeff set to use */
                    for (k = 0; k < 12; k++)
                    { /* This is either equation (A-4) or the corresponding part of (A-7), depending on the value of set */
                        // if (current_max_iterations != C.max_iterations)
                        // {
                        //     printf("u_old: %f, node: %i, d_node:%i, other_node: %i, coef:%f\n", u_old[offset_ind(node, d_node[k])], node, d_node[k], offset_ind(node, d_node[k]), C.coeff[set][k]);
                        // }
                        u_00 += u[offset_ind(node, d_node[k])] * C.coeff[set][k];
                    }
                    if (set == CONSTRAINED)
                    {                                                                    /* Solution is (A-7) and modifications depend on which quadrant the point lies in */
                        auto &b = C.brgs_coefs[briggs_index].b;                          /* Shorthand to this node's Briggs b-array */
                        unsigned int quadrant = static_cast<unsigned int>(status[node]); /* Which quadrant did the point fall in? */
                        for (k = 0, sum_bk_uk = 0.0; k < 4; k++)
                        { /* Sum over b[k]*u[k] for nodes A-D in Fig A-1 */
                            sum_bk_uk += b[k] * u[offset_ind(node, d_node[p[quadrant][k]])];
                        }
                        u_00 = (u_00 + C.a0_const_2 * (sum_bk_uk + b[4])) * b[5]; /* Add point E in Fig A-1 to sum_bk_uk and normalize */
                        briggs_index++;                                           /* Got to next sequential Briggs array index */
                    }
                    /* We now apply the over-relaxation: */
                    u_00 = u[node] * C.relax_old + u_00 * C.relax_new;

                    if (limited)
                    { /* clip to bounds on final spacing.  */
                        /* Must use final spacing node index to access the Bound grids */
                        if (is_lower_bounded && !std::isnan(C.lower_bound[node_final]))
                            u_00 = std::max(u_00, C.lower_bound[node_final]);
                        if (is_upper_bounded && !std::isnan(C.upper_bound[node_final]))
                            u_00 = std::min(u_00, C.upper_bound[node_final]);
                        // this could go in the for loop statement, but is techinally only needed if it is constrained.
                        node_final += C.current_stride;
                    }
                    u_change = std::abs(u_00 - u[node]); /* Change in node value between iterations */
                    u[node] = u_00;                      /* Our updated estimate at this node */
                    if (u_change > max_u_change)
                        max_u_change = u_change; /* Keep track of max u_change across all nodes */
                } /* End of loop over columns */
            } /* End of loop over rows [and possibly threads via OpenMP] */
            iteration_count++;
            // ;                  /* Update iteration counts for this stride and for total */
            max_z_change = max_u_change * C.z_rms; /* Scale max_u_change back into original z units -> max_z_change */
            // GMT_Report(GMT->parent, GMT_MSG_DEBUG, C.format,
            //            C.current_stride, C.mode_type[mode], iteration_count, max_z_change, current_limit, C.total_iterations);
            // if (C.logging)
            //     fprintf(C.fp_log, "%d\t%c\t%" PRIu64 "\t%.8g\t%.8g\t%" PRIu64 "\n", C.current_stride, C.mode_type[mode], iteration_count, max_z_change, current_limit, C.total_iterations);
            finished = (max_z_change <= current_limit || iteration_count >= current_max_iterations);

        } while (!finished);

        surface_set_BCs(C, u); /* Set the boundary rows and columns on the output*/

        if (C.verbosity > 1)
        {
            std::cout << "Iterations " << iteration_count << std::endl;
        }

        return iteration_count;
    }

    void surface_check_errors(const SurfaceInfo &C, const Grid &G)
    {
        /* Compute misfits at original data locations,  This is only done at the
         * final grid resolution, hence current_stride == 1. */

        size_t node;
        size_t row, col;
        auto &d_node = C.offset;
        const std::vector<NodeStatus> &status = C.status;

        double x0, y0, dx, dy, mean_error = 0.0, mean_squared_error = 0.0, z_est, z_err, curvature, c;
        double du_dx, du_dy, d2u_dx2, d2u_dxdy, d2u_dy2, d3u_dx3, d3u_dx2dy, d3u_dxdy2, d3u_dy3;
        double r_inc_x = 1.0 / G.header.inc[GMT_X];
        double r_inc_y = 1.0 / G.header.inc[GMT_Y];

        const std::vector<double> &u = G.data;
        const GridHeader &h = G.header;

        if (C.verbosity > 1)
        {
            std::cout << "Compute rms misfit and curvature." << std::endl;
        }

        // surface_set_BCs(C, u); /* First update the boundary values */

        /* Estimate solution at all data constraints using 3rd order Taylor expansion from nearest node */

        for (auto &pt : C.points)
        {
            row = index_to_row(pt.index, C.n_columns);
            col = index_to_col(pt.index, C.n_columns);
            node = row_col_to_node(row, col, C.mx);
            if (status[node] == NodeStatus::IS_CONSTRAINED)
                continue; /* Since misfit by definition is zero so no point adding it */
            /* Get coordinates of this node */
            x0 = col_to_x(col, h.wesn[XLO], h.wesn[XHI], h.inc[GMT_X], h.n_columns);
            y0 = row_to_y(row, h.wesn[YLO], h.wesn[YHI], h.inc[GMT_Y], h.n_rows);
            /* Get dx,dy of data point away from this node */
            dx = x_to_fcol(pt.x, x0, r_inc_x);
            dy = y_to_frow(pt.y, y0, r_inc_y);

            du_dx = 0.5 * (u[offset_ind(node, d_node[E1])] - u[offset_ind(node, d_node[W1])]);
            du_dy = 0.5 * (u[offset_ind(node, d_node[N1])] - u[offset_ind(node, d_node[S1])]);
            d2u_dx2 = u[offset_ind(node, d_node[E1])] + u[offset_ind(node, d_node[W1])] - 2 * u[node];
            d2u_dy2 = u[offset_ind(node, d_node[N1])] + u[offset_ind(node, d_node[S1])] - 2 * u[node];
            d2u_dxdy = 0.25 * (u[offset_ind(node, d_node[NE])] - u[offset_ind(node, d_node[NW])] - u[offset_ind(node, d_node[SE])] + u[offset_ind(node, d_node[SW])]);
            d3u_dx3 = 0.5 * (u[offset_ind(node, d_node[E2])] - 2 * u[offset_ind(node, d_node[E1])] + 2 * u[offset_ind(node, d_node[W1])] - u[offset_ind(node, d_node[W2])]);
            d3u_dy3 = 0.5 * (u[offset_ind(node, d_node[N2])] - 2 * u[offset_ind(node, d_node[N1])] + 2 * u[offset_ind(node, d_node[S1])] - u[offset_ind(node, d_node[S2])]);
            d3u_dx2dy = 0.5 * ((u[offset_ind(node, d_node[NE])] + u[offset_ind(node, d_node[NW])] - 2 * u[offset_ind(node, d_node[N1])]) - (u[offset_ind(node, d_node[SE])] + u[offset_ind(node, d_node[SW])] - 2 * u[offset_ind(node, d_node[S1])]));
            d3u_dxdy2 = 0.5 * ((u[offset_ind(node, d_node[NE])] + u[offset_ind(node, d_node[SE])] - 2 * u[offset_ind(node, d_node[E1])]) - (u[offset_ind(node, d_node[NW])] + u[offset_ind(node, d_node[SW])] - 2 * u[offset_ind(node, d_node[W1])]));

            /* Compute the 3rd order Taylor approximation from current node */

            z_est = u[node] + dx * (du_dx + dx * ((0.5 * d2u_dx2) + dx * (d3u_dx3 / 6.0))) + dy * (du_dy + dy * ((0.5 * d2u_dy2) + dy * (d3u_dy3 / 6.0))) + dx * dy * (d2u_dxdy) + (0.5 * dx * d3u_dx2dy) + (0.5 * dy * d3u_dxdy2);

            z_err = z_est - pt.z; /* Misfit between surface estimate and observation */
            mean_error += z_err;
            mean_squared_error += (z_err * z_err);
        }
        mean_error /= C.points.size();
        mean_squared_error = sqrt(mean_squared_error / C.points.size());

        /* Compute the total curvature of the grid */

        curvature = 0.0;
        for (row = 0; row < h.n_rows; row++)
            for (col = 0, node = ijp(h, row, (size_t)0); col < h.n_columns; col++, node++)
            {
                c = u[offset_ind(node, d_node[E1])] + u[offset_ind(node, d_node[W1])] + u[offset_ind(node, d_node[N1])] + u[offset_ind(node, d_node[S1])] - 4.0 * u[offset_ind(node, d_node[E1])];
                curvature += (c * c);
            }

        if (C.verbosity > 0)
        {
            std::cout << "Fit info: N data points = " << C.points.size() << " N nodes = " << C.nxny << " mean error = " << mean_error << " rms error = " << mean_squared_error << " curvature = " << curvature << std::endl;
        }
    }

    void surface_remove_planar_trend(SurfaceInfo &C, const GridHeader &h)
    {
        /* Fit LS plane and remove trend from our (x,y,z) input data; we add trend to grid before output.
         * Note: Here, x and y are first converted to fractional grid spacings from 0 to {n_columns,n_rows}-1.
         * Hence the same scheme is used by replace_planar trend (i.e., just use row,col as coordinates).
         * Note: The plane is fit to the original data z-values before normalizing by rms. */
        double a, b, c, d, xx, y_up, zz, sx, sy, sz, sxx, sxy, sxz, syy, syz;
        double r_inc_x = 1.0 / h.inc[GMT_X];
        double r_inc_y = 1.0 / h.inc[GMT_Y];

        sx = sy = sz = sxx = sxy = sxz = syy = syz = 0.0;

        for (auto &pt : C.points)
        {                                                 /* Sum up normal equation terms */
            xx = x_to_fcol(pt.x, h.wesn[XLO], r_inc_x);   /* Distance from west to this point */
            y_up = y_to_frow(pt.y, h.wesn[YLO], r_inc_y); /* Distance from south to this point */
            zz = pt.z;
            sx += xx;
            sy += y_up;
            sz += zz;
            sxx += (xx * xx);
            sxy += (xx * y_up);
            sxz += (xx * zz);
            syy += (y_up * y_up);
            syz += (y_up * zz);
        }

        d = C.points.size() * sxx * syy + 2 * sx * sy * sxy - C.points.size() * sxy * sxy - sx * sx * syy - sy * sy * sxx;

        if (d == 0.0)
        { /* When denominator is zero we have a horizontal plane */
            C.plane_icept = C.plane_sx = C.plane_sy = 0.0;
            return;
        }

        a = sz * sxx * syy + sx * sxy * syz + sy * sxy * sxz - sz * sxy * sxy - sx * sxz * syy - sy * syz * sxx;
        b = C.points.size() * sxz * syy + sz * sy * sxy + sy * sx * syz - C.points.size() * sxy * syz - sz * sx * syy - sy * sy * sxz;
        c = C.points.size() * sxx * syz + sx * sy * sxz + sz * sx * sxy - C.points.size() * sxy * sxz - sx * sx * syz - sz * sy * sxx;

        C.plane_icept = a / d;
        C.plane_sx = b / d;
        C.plane_sy = c / d;

        for (auto &pt : C.points)
        { /* Now remove this plane from the data constraints */
            xx = x_to_fcol(pt.x, h.wesn[XLO], r_inc_x);
            y_up = y_to_frow(pt.y, h.wesn[YLO], r_inc_y);
            pt.z -= evaluate_plane(C, xx, y_up);
        }

        if (C.verbosity > 1)
        {
            std::cout << "Plane fit z = " << C.plane_icept << " + (" << C.plane_sx << " * col) + (" << C.plane_sy << " * row)" << std::endl;
        }
    }

    void surface_restore_planar_trend(const SurfaceInfo &C, Grid &G)
    {
        /* Scale grid back up by the data rms and restore the least-square plane.
         * Note: In determining the plane and in evaluating it, remember that the
         * x and y coordinates needed are not data coordinates but fractional col
         * and row distances from an origin at the lower left (southwest corner).
         * This means the y-values are positive up and increase in the opposite
         * direction than how rows increase. Hence the use of y_up below. */
        size_t row, col;
        size_t node;
        std::vector<double> &u = G.data;
        size_t y_up; /* Measure y up from south in fractional rows */

        for (row = 0; row < G.header.n_rows; row++)
        {
            y_up = G.header.n_rows - row - 1;                      /* # of rows from south (where y_up = 0) to this node */
            node = row_col_to_node(row, (size_t)0, C.current_mx);  /* Node index at left end of interior row */
            for (col = 0; col < G.header.n_columns; col++, node++) /* March across this row */
                u[node] = ((u[node] * C.z_rms) + (evaluate_plane(C, col, y_up)));
        }
    }

    void cull_data(SurfaceInfo &C)
    {
        /* We eliminate data which will become unusable on the final iteration, when current_stride = 1.
           It assumes current_stride = 1 and that set_grid_parameters has been called.
           We sort, mark redundant data as SURFACE_OUTSIDE, and sort again, chopping off the excess.
        */

        size_t last_index = UINTMAX_MAX, n_outside = 0;

        if (C.verbosity > 1)
        {

            std::cout << "Eliminate data points that are not nearest a node." << std::endl;
        }

        /* Sort the data  */

        std::sort(C.points.begin(), C.points.end(), compare_points{&C});

        /* If more than one datum is indexed to the same node, only the first should be kept.
           Mark the additional ones as SURFACE_OUTSIDE.
        */
        std::vector<Data *> index_points;
        for (auto &pt : C.points)
        {
            if (C.downsample_mode == DownsampleMode::CLOSEST)
            {

                if (pt.index == last_index)
                { /* Same node but further away than our guy */
                    pt.index = SURFACE_OUTSIDE;
                    n_outside++;
                    // std ::cout << "Skipping unusable point at (" << C.points[k].x << " " << C.points[k].y << " " << C.points[k].z << ") as ("
                    //            << C.points[last_k].x << " " << C.points[last_k].y << " " << C.points[last_k].z << ") is closer to node " << last_index << std::endl;
                }
                else
                { /* New index, just update last_index */
                    last_index = pt.index;
                }
            }
            else
            {
                // build up a list of all points with the same index
                if (index_points.size() == 0) // Should only happen at the first iteration.
                {
                    last_index = pt.index;
                    index_points.push_back(&pt);
                }
                else if (last_index == pt.index)
                {
                    pt.index = SURFACE_OUTSIDE;
                    index_points.push_back(&pt);
                }
                else
                {
                    // Was a new node! (and not the first iteration)
                    // first do the statistics on the clump (if it's greater than 1)
                    std::size_t n_pts = index_points.size();
                    if (n_pts > 1)
                    {
                        double new_x = 0.0, new_y = 0.0, new_z = 0.0;
                        if (C.downsample_mode == DownsampleMode::MEAN)
                        {
                            for (auto pti : index_points)
                            {
                                new_x += pti->x;
                                new_y += pti->y;
                                new_z += pti->z;
                            }
                            new_x /= n_pts;
                            new_y /= n_pts;
                            new_z /= n_pts;
                        }
                        else
                        {
                            // sort the points
                            std::sort(
                                index_points.begin(),
                                index_points.end(),
                                [](const Data *p1, const Data *p2)
                                {
                                    return p1->z < p2->z;
                                });
                            if (C.downsample_mode == DownsampleMode::MEDIAN)
                            {
                                if (n_pts % 2 == 0)
                                {
                                    auto p1 = index_points[n_pts / 2 - 1];
                                    auto p2 = index_points[n_pts / 2];
                                    new_x = (p1->x + p2->x) / 2;
                                    new_y = (p1->y + p2->y) / 2;
                                    new_z = (p1->z + p2->z) / 2;
                                }
                                else
                                {
                                    auto p1 = index_points[n_pts / 2];
                                    new_x = p1->x;
                                    new_y = p1->y;
                                    new_z = p1->z;
                                }
                            }
                            else // Mode
                            {
                                // auto best = index_points[0];
                                // size_t best_count = 1;

                                // auto current = index_points[0];
                                // size_t count = 1;

                                // for (size_t i = 1; i < n_pts; ++i)
                                // {
                                //     if (index_points[i]->z == current->z)
                                //     {
                                //         count++;
                                //     }
                                //     else
                                //     {
                                //         if (count > best_count)
                                //         {
                                //             best = current;
                                //             best_count = count;
                                //         }
                                //         current = index_points[i];
                                //         count = 1;
                                //     }
                                // }

                                // if (count > best_count)
                                //     best = current;
                                // new_x = best->x;
                                // new_y = best->y;
                                // new_z = best->z;
                            }

                            auto p_first = index_points.front();
                            p_first->x = new_x;
                            p_first->y = new_y;
                            p_first->z = new_z;
                        }
                        // we've updated the first value in the index_points vector
                        // reset it and add this new point
                        last_index = pt.index;
                        index_points.resize(0);
                        index_points.push_back(&pt);
                    }
                }
            }
        }

        if (n_outside)
        { /* Sort again; this time the SURFACE_OUTSIDE points will be sorted to end of the array */

            // QSORT_R(C.points, C.npoints, sizeof(struct SURFACE_DATA), compare_points, &(C.info));
            std::sort(C.points.begin(), C.points.end(), compare_points{&C});
            // C.npoints -= n_outside;                   /* Effectively chopping off the eliminated points */
            C.points.resize(C.points.size() - n_outside); // gmt_M_memory(GMT, C.points, C.npoints, struct SURFACE_DATA); /* Adjust memory accordingly */
            if (C.verbosity > 1)
            {
                std::cout << n_outside << " unusable points were supplied; these will be ignored." << std::endl;
            }
            // GMT_Report(GMT->parent, GMT_MSG_WARNING, "%" PRIu64 " unusable points were supplied; these will be ignored.\n", n_outside);
            // GMT_Report(GMT->parent, GMT_MSG_WARNING, "You should have pre-processed the data with block-mean, -median, or -mode.\n");
            // GMT_Report(GMT->parent, GMT_MSG_WARNING, "Check that previous processing steps write results with enough decimals.\n");
            // GMT_Report(GMT->parent, GMT_MSG_WARNING, "Possibly some data were half-way between nodes and subject to IEEE 754 rounding.\n");
        }
        C.points_sub = std::span<const Data>(C.points);
    }

    bool rescale_z_values(SurfaceInfo &C)
    {
        /* Find and normalize data by their rms value */
        double ssz = 0.0;

        for (auto &pt : C.points)
            ssz += (pt.z * pt.z);
        C.z_rms = sqrt(ssz / C.points.size());
        if (C.verbosity > 1)
        {
            std::cout << "Normalize detrended data constraints by z rms = " << C.z_rms << std::endl;
        }

        if (C.z_rms < SURFACE_CONV8_LIMIT)
        {
            if (C.verbosity > 0)
            {
                std::cout << "Input data lie exactly on a plane." << std::endl;
            }
            C.r_z_rms = C.z_rms = 1.0;
            return true; /* Flag to tell the main to just write out the plane */
        }
        else
            C.r_z_rms = 1.0 / C.z_rms;

        for (auto &pt : C.points)
            pt.z *= C.r_z_rms;

        if (C.converge_limit == 0.0 || C.converge_mode == ConvergenceLimitMode::BY_PERCENT)
        { /* Set default values for convergence criteria */
            double limit = (C.converge_mode == ConvergenceLimitMode::BY_PERCENT) ? C.converge_limit : SURFACE_CONV_LIMIT;
            C.converge_limit = limit * C.z_rms; /* i.e., 100 ppm of L2 scale */
            if (C.verbosity > 1)
            {
                auto ppm = std::llround(limit * 1.0e6);
                std ::cout << "Select default convergence limit of " << C.converge_limit << " (" << ppm << " ppm of L2 scale)" << std::endl;
            }
        }
        return false;
    }

    template <typename T>
    T gcd_euclid(T a, T b)
    {
        /* Returns the greatest common divisor of u and v by Euclid's method.
         * I have experimented also with Stein's method, which involves only
         * subtraction and left/right shifting; Euclid is faster, both for
         * integers of size 0 - 1024 and also for random integers of a size
         * which fits in a long integer.  Stein's algorithm might be better
         * when the integers are HUGE, but for our purposes, Euclid is fine.
         *
         * Walter H. F. Smith, 25 Feb 1992, after D. E. Knuth, vol. II  */

        T u, v, r;

        u = std::max(a, b);
        v = std::min(a, b);

        while (v > 0)
        {
            r = u % v; /* Knuth notes that u < 2v 40% of the time;  */
            u = v;     /* thus we could have tried a subtraction  */
            v = r;     /* followed by an if test to do r = u%v  */
        }
        return (u);
    }

    std::vector<size_t> get_prime_factors(size_t n)
    {
        /* Fills the integer array f with the prime factors of n.
         * Returns the number of locations filled in f, which is
         * one if n is prime.
         *
         * f[] should have been malloc'ed to enough space before
         * calling prime_factors().  We can be certain that f[32]
         * is enough space, for if n fits in a long, then n < 2**32,
         * and so it must have fewer than 32 prime factors.  I think
         * that in general, ceil(log2((double)n)) is enough storage
         * space for f[].
         *
         * Tries 2,3,5 explicitly; then alternately adds 2 or 4
         * to the previously tried factor to obtain the next trial
         * factor.  This is done with the variable two_four_toggle.
         * With this method we try 7,11,13,17,19,23,25,29,31,35,...
         * up to a maximum of sqrt(n).  This shortened list results
         * in 1/3 fewer divisions than if we simply tried all integers
         * between 5 and sqrt(n).  We can reduce the size of the list
         * of trials by an additional 20% by removing the multiples
         * of 5, which are equal to 30m +/- 5, where m >= 1.  Starting
         * from 25, these are found by alternately adding 10 or 20.
         * To do this, we use the variable ten_twenty_toggle.
         *
         * W. H. F. Smith, 26 Feb 1992, after D.E. Knuth, vol. II  */

        size_t current_factor = 0;         /* The factor currently being tried  */
        size_t max_factor;                 /* Don't try any factors bigger than this  */
        size_t two_four_toggle = 0;        /* Used to add 2 or 4 to get next trial factor  */
        size_t ten_twenty_toggle = 0;      /* Used to add 10 or 20 to skip_five  */
        size_t skip_five = 25;             /* Used to skip multiples of 5 in the list  */
        size_t base_factor[3] = {2, 3, 5}; /* Standard factors to try */
        uint64_t m = n;                    /* Used to keep a working copy of n  */

        /* Initialize max_factor  */
        std::vector<size_t> f;

        if (m < 2)
            return f;
        max_factor = static_cast<size_t>(std::sqrt(m));

        /* First find the 2s, 3s, and 5s */
        for (auto k = 0; k < 3; k++)
        {
            current_factor = base_factor[k];
            while (!(m % current_factor))
            {
                m /= current_factor;
                f.push_back(current_factor);
                // f[n_factors++] = current_factor;
            }
            if (m == 1)
                return f;
        }

        /* Unless we have already returned we now try all the rest  */

        while (m > 1 && current_factor <= max_factor)
        {

            /* Current factor is either 2 or 4 more than previous value  */

            if (two_four_toggle)
            {
                current_factor += 4;
                two_four_toggle = 0;
            }
            else
            {
                current_factor += 2;
                two_four_toggle = 1;
            }

            /* If current factor is a multiple of 5, skip it.  But first,
                set next value of skip_five according to 10/20 toggle:  */

            if (current_factor == skip_five)
            {
                if (ten_twenty_toggle)
                {
                    skip_five += 20;
                    ten_twenty_toggle = 0;
                }
                else
                {
                    skip_five += 10;
                    ten_twenty_toggle = 1;
                }
                continue;
            }

            /* Get here when current_factor is not a multiple of 2,3 or 5:  */

            while (!(m % current_factor))
            {
                m /= current_factor;
                f.push_back(current_factor);
                // f[n_factors++] = current_factor;
            }
        }

        /* Get here when all factors up to floor(sqrt(n)) have been tried.  */

        if (m > 1)
            f.push_back(m);
        // f[n_factors++] = (unsigned int)m; /* m is an additional prime factor of n  */

        return f;
    }

    namespace suggestions
    {

        /*! Definition of structure use for finding optimal n_columns/n_rows for surface */
        struct Suggestion
        { /* Used to find top ten list of faster grid dimensions  */
            size_t n_columns;
            size_t n_rows;
            double factor; /* Speed up by a factor of factor  */
        };
        double guess_surface_time(size_t n_columns, size_t n_rows)
        {
            /* Routine to guess a number proportional to the operations
             * required by surface working on a user-desired grid of
             * size n_columns by n_rows, where n_columns = (x_max - x_min)/dx, and same for
             * n_rows.  (That is, one less than actually used in routine.)
             *
             * This is based on the following untested conjecture:
             * 	The operations are proportional to T = nxg*nyg*L,
             *	where L is a measure of the distance that data
             *	constraints must propagate, and nxg, nyg are the
             * 	current size of the grid.
             *	For n_columns,n_rows relatively prime, we will go through only
             * 	one grid cycle, L = max(n_columns,n_rows), and T = n_columns*n_rows*L.
             *	But for n_columns,n_rows whose greatest common divisor is a highly
             * 	composite number, we will have L equal to the division
             * 	step made at each new grid cycle, and nxg,nyg will
             * 	also be smaller than n_columns,n_rows.  Thus we can hope to find
             *	some n_columns,n_rows for which the total value of T is C.small.
             *
             * The above is pure speculation and has not been derived
             * empirically.  In actual practice, the distribution of the
             * data, both spatially and in terms of their values, will
             * have a strong effect on convergence.
             *
             * W. H. F. Smith, 26 Feb 1992.  */

            size_t gcd;      /* Current value of the gcd  */
            size_t nxg, nyg; /* Current value of the grid dimensions  */
            size_t factor;   /* Currently used factor  */
            /* Doubles are used below, even though the values will be integers,
                because the multiplications might reach sizes of O(n**3)  */
            double t_sum;  /* Sum of values of T at each grid cycle  */
            double length; /* Current propagation distance.  */
            std::vector<size_t> factors;

            gcd = gcd_euclid(n_columns, n_rows);
            if (gcd > 1)
            {
                factors = get_prime_factors(gcd);
                nxg = n_columns / gcd;
                nyg = n_rows / gcd;
                if (nxg < 3 || nyg < 3)
                {
                    factor = factors.back();
                    factors.pop_back();
                    gcd /= factor;
                    nxg *= factor;
                    nyg *= factor;
                }
            }
            else
            {
                nxg = n_columns;
                nyg = n_rows;
            }
            length = (double)std::max(nxg, nyg);
            t_sum = nxg * (nyg * length); /* Make it double at each multiply  */

            /* Are there more grid cycles ?  */
            while (gcd > 1)
            {
                factor = factors.back();
                factors.pop_back();
                gcd /= factor;
                nxg *= factor;
                nyg *= factor;
                length = (double)factor;
                t_sum += nxg * (nyg * length);
            }
            return (t_sum);
        }

        std::vector<Suggestion> optimal_dim_for_surface(size_t n_columns, size_t n_rows)
        {
            /* Calls guess_surface_time for a variety of trial grid
             * sizes, where the trials are highly composite numbers
             * with lots of factors of 2, 3, and 5.  The sizes are
             * within the range (n_columns,n_rows) - (2*n_columns, 2*n_rows).  Prints to
             * GMT->session.std[GMT_ERR] the values which are an improvement over the
             * user's original n_columns,n_rows.
             * Should be called with n_columns=(x_max-x_min)/dx, and ditto
             * for n_rows; that is, one smaller than the lattice used
             * in surface.c
             *
             * W. H. F. Smith, 26 Feb 1992.  */

            double users_time;                   /* Time for user's n_columns, n_rows  */
            double current_time;                 /* Time for current nxg, nyg  */
            size_t nxg, nyg;                     /* Guessed by this routine  */
            size_t nx2, ny2, nx3, ny3, nx5, ny5; /* For powers  */
            size_t xstop, ystop;                 /* Set to 2*n_columns, 2*n_rows  */
            std::vector<Suggestion> sug_list;

            users_time = guess_surface_time(n_columns, n_rows);
            xstop = 2 * n_columns;
            ystop = 2 * n_rows;

            for (nx2 = 2; nx2 <= xstop; nx2 *= 2)
            {
                for (nx3 = 1; nx3 <= xstop; nx3 *= 3)
                {
                    for (nx5 = 1; nx5 <= xstop; nx5 *= 5)
                    {
                        nxg = nx2 * nx3 * nx5;
                        if (nxg < n_columns || nxg > xstop)
                            continue;

                        for (ny2 = 2; ny2 <= ystop; ny2 *= 2)
                        {
                            for (ny3 = 1; ny3 <= ystop; ny3 *= 3)
                            {
                                for (ny5 = 1; ny5 <= ystop; ny5 *= 5)
                                {
                                    nyg = ny2 * ny3 * ny5;
                                    if (nyg < n_rows || nyg > ystop)
                                        continue;

                                    current_time = guess_surface_time(nxg, nyg);
                                    if (current_time < users_time)
                                    {
                                        Suggestion sug;
                                        sug.n_columns = nxg;
                                        sug.n_rows = nyg;
                                        sug.factor = users_time / current_time;

                                        sug_list.push_back(sug);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            std::sort(sug_list.begin(), sug_list.end(), [](const Suggestion &a, const Suggestion &b)
                      {
                          return a.factor < b.factor; // Sort in ascending order
                      });
            return sug_list;
        }
    }

    size_t suggest_sizes(SurfaceInfo &C, GridHeader &header)
    {
        /* Calls optimal_dim_for_surface to determine if there are
         * better choices for n_columns, n_rows that might speed up calculations
         * by having many more common factors.
         *
         * W. H. F. Smith, 26 Feb 1992.  */

        size_t n_columns = header.n_columns - 1;
        size_t n_rows = header.n_rows - 1;
        auto sug = suggestions::optimal_dim_for_surface(n_columns, n_rows);
        if (C.verbosity > 1)
        {
            for (auto &s : sug)
            {
                std::cout << "suggestions" << s.n_columns << ", " << s.n_rows << "factor" << std::endl;
            }
        }

        if (sug.size() > 0)
        { /* We did find some suggestions, report them (up to the first 10 suggestions) */
            // char region[GMT_LEN128] = {""}, buffer[GMT_LEN128] = {""};
            // bool lat_bad = false;
            size_t m; //, save_range = GMT->current.io.geo.range;
            // GMT->current.io.geo.range = GMT_IS_GIVEN_RANGE; /* Override this setting explicitly */
            m = sug.back().n_columns - n_columns; /* Additional nodes needed in x to give more factors */
            header.n_columns += m;
            header.wesn[XLO] -= (m / 2) * header.inc[GMT_X];
            header.wesn[XHI] += ((m + 1) / 2) * header.inc[GMT_X];

            m = sug.back().n_rows - n_rows; /* Additional nodes needed in y to give more factors */
            header.n_rows += m;
            header.wesn[YLO] -= (m / 2) * header.inc[GMT_Y];
            header.wesn[YHI] += ((m + 1) / 2) * header.inc[GMT_Y];

            if (C.verbosity > 1)
            {
                std::cout << "Internally padded up to :  " << header.n_columns - 1 << " " << header.n_rows - 1 << std::endl;
            }

            // probably also need to offset wesn by the same amount...
        }
        return sug.size();
    }

    void store_grid(const SurfaceInfo &C, const Grid &G, double *out, size_t nx, size_t ny)
    {
        /* Write output grid to file */
        // char *limit[2] = {"lower", "upper"};
        auto &u = G.data;
        // if ny != C.current_ny - 1, then it was padded to a slightly larger size for convergence.
        // (and same for nx).
        // these should be the shift in row and column coordinates from the start
        size_t nx0 = (C.current_nx - nx) / 2;
        size_t ny0 = (C.current_ny - ny + 1) / 2; // +1 because extra padding is at the end of y
        for (size_t iy = 0, iy_g = ny0; iy < ny; ++iy, ++iy_g)
        {
            size_t ig = C.node_nw_corner + iy_g * C.current_mx + nx0;
            for (size_t ix = 0; ix < nx; ++ix, ++ig) // columns
            {
                out[iy * nx + ix] = u[ig];
            }
        }
    }

    void load_constraints(SurfaceInfo &C,
                          const GridHeader &h,
                          const double *lower, const size_t n_lower,
                          const double *upper, const size_t n_upper,
                          const size_t nx, const size_t ny)
    {
        /* Deal with the constants or grids supplied via -L.  Note: Because we remove a
         * best-fitting plane from the data, even a simple constant constraint will become
         * a plane and thus must be represented on a grid. */
        size_t col, row, y_up;
        size_t node;

        if (C.verbosity > 1)
            std::cout << "Load any data constraint limit grids" << std::endl;

        /* Load lower/upper limits, verify range, deplane, and rescale */
        // initializes lower and upper bound vectors to nan, incase there is any padding from added
        // cells for convergence acceleration.
        size_t nx0 = (C.current_nx - nx) / 2;
        size_t ny0 = (C.current_ny - ny + 1) / 2; // +1 because extra padding is at the end of y
        if (n_lower > 0)
        {
            C.lower_bound.resize(C.mxmy, std::numeric_limits<double>::quiet_NaN());
            if (n_lower == SURFACE_USE_DATA_LIMITS || n_lower == 1)
            {
                // if n_lower == 1 use the same value, else use the data limits
                double lim = n_lower == 1 ? lower[0] : C.z_min;
                if (lim > C.z_min)
                {
                    std::cout << "Your lower bound is larger than the minimum data value." << std::endl;
                }

                for (size_t iy = 0, iy_g = ny0; iy < ny; ++iy, ++iy_g) // rows
                {
                    size_t ig = C.node_nw_corner + iy_g * C.current_mx + nx0;
                    for (size_t ix = 0; ix < nx; ++ix, ++ig) // columns
                    {
                        C.lower_bound[ig] = lim;
                    }
                }
            }
            else
            {
                for (size_t iy = 0, iy_g = ny0; iy < ny; ++iy, ++iy_g) // rows
                {
                    size_t ig = C.node_nw_corner + iy_g * C.current_mx + nx0;
                    for (size_t ix = 0; ix < nx; ++ix, ++ig) // columns
                    {
                        C.lower_bound[ig] = lower[iy * nx + ix];
                    }
                }
            }

            /* Remove best-fitting plane and normalize the bounding values */
            for (row = 0; row < h.n_rows; row++)
            {
                y_up = (h.n_rows - row - 1); /* Require y_up = 0 at south and positive toward north */
                node = row_col_to_node(row, size_t(0), C.mx);
                for (col = 0; col < h.n_columns; col++, node++)
                {
                    if (std::isnan(C.lower_bound[node]))
                        continue;
                    C.lower_bound[node] -= evaluate_plane(C, col, y_up); /* Remove plane */
                    C.lower_bound[node] *= C.r_z_rms;                    /* Normalize residuals */
                }
            }
        }
        if (n_upper > 0)
        {
            C.upper_bound.resize(C.mxmy, std::numeric_limits<double>::quiet_NaN());

            if (n_upper == SURFACE_USE_DATA_LIMITS || n_upper == 1)
            {
                // if n_lower == 1 use the same value, else use the data limits
                double lim = n_upper == 1 ? upper[0] : C.z_max;

                if (lim < C.z_max)
                {
                    std::cout << "Your upper bound is less than the maximum data value." << std::endl;
                }

                for (size_t iy = 0, iy_g = ny0; iy < ny; ++iy, ++iy_g) // rows
                {
                    size_t ig = C.node_nw_corner + iy_g * C.current_mx + nx0;
                    for (size_t ix = 0; ix < nx; ++ix, ++ig) // columns
                    {
                        C.upper_bound[ig] = lim;
                    }
                }
            }
            else
            {
                for (size_t iy = 0, iy_g = ny0; iy < ny; ++iy, ++iy_g) // rows
                {
                    size_t ig = C.node_nw_corner + iy_g * C.current_mx + nx0;
                    for (size_t ix = 0; ix < nx; ++ix, ++ig) // columns
                    {
                        C.upper_bound[ig] = upper[iy * nx + ix];
                    }
                }
            }

            /* Remove best-fitting plane and normalize the bounding values */
            for (row = 0; row < h.n_rows; row++)
            {
                y_up = (h.n_rows - row - 1); /* Require y_up = 0 at south and positive toward north */
                node = row_col_to_node(row, size_t(0), C.mx);
                for (col = 0; col < h.n_columns; col++, node++)
                {
                    if (std::isnan(C.upper_bound[node]))
                        continue;
                    C.upper_bound[node] -= evaluate_plane(C, col, y_up); /* Remove plane */
                    C.upper_bound[node] *= C.r_z_rms;                    /* Normalize residuals */
                }
            }
        }
    }

    void surface(
        double *x, double *y, double *z, size_t n_nodes,
        double xmin, double xmax, double ymin, double ymax, size_t nx, size_t ny,
        double *output,
        DownsampleMode downsample_mode = DownsampleMode::CLOSEST,
        size_t max_iterations = 500, double relax = SURFACE_OVERRELAXATION,
        double alpha = 1.0,
        double b_tension = 0.0, double i_tension = 0.0,
        double converge_limit = 1.0E-4,
        ConvergenceLimitMode converge_mode = ConvergenceLimitMode::BY_PERCENT,
        unsigned char verbosity = 0,
        double *lower = NULL,
        size_t n_lower = 0,
        double *upper = NULL,
        size_t n_upper = 0)
    {

        Grid G;

        struct SurfaceInfo C;
        C.verbosity = verbosity;
        for (size_t i = 0; i < nx * ny; ++i)
        {
            output[i] = 0.0;
        }
        /* Determine if there is a better region that would allow more intermediate resolutions to converge better */
        G.header.n_columns = nx;
        G.header.n_rows = ny;
        double dx = (xmax - xmin) / (nx - 1);
        double dy = (ymax - ymin) / (ny - 1);

        G.header.wesn[XLO] = xmin;
        G.header.wesn[XHI] = xmax;
        G.header.wesn[YLO] = ymin;
        G.header.wesn[YHI] = ymax;
        // padding needed for boundary conditions here
        G.header.pad[XLO] = 2;
        G.header.pad[XHI] = 2;
        G.header.pad[YLO] = 2;
        G.header.pad[YHI] = 2;
        G.header.inc[GMT_X] = dx;
        G.header.inc[GMT_Y] = dy;

        // increase padding a bit to make sure it's all good.
        suggest_sizes(C, G.header);
        G.header.nm = G.header.n_columns * G.header.n_rows;

        G.header.mx = G.header.n_columns + G.header.pad[XLO] + G.header.pad[XHI];
        G.header.my = G.header.n_rows + G.header.pad[YLO] + G.header.pad[YHI];
        G.header.size = G.header.mx * G.header.my;

        C.relax_new = relax;
        C.relax_old = 1.0 - relax;
        C.max_iterations = max_iterations;
        C.boundary_tension = b_tension;
        C.interior_tension = i_tension;
        C.downsample_mode = downsample_mode;
        C.alpha = alpha * G.header.inc[GMT_Y] / G.header.inc[GMT_X]; // alpha is dy/dx?
        C.converge_limit = converge_limit;
        C.converge_mode = converge_mode;
        C.total_iterations = 0;
        C.previous_stride = 0;
        C.tension = 0; // this is unused
        C.z_rms = 1.0;
        C.r_z_rms = 1.0;
        C.n_columns = G.header.n_columns;
        C.n_rows = G.header.n_rows;
        C.nxny = G.header.nm;
        C.mx = G.header.mx;
        C.my = G.header.my;
        C.mxmy = G.header.size;

        std::copy(G.header.wesn, G.header.wesn + 4, C.wesn);

        if (G.header.n_columns < 4 || G.header.n_rows < 4)
        {
            std::cerr << "Grid must have at least 4 nodes in each direction (you have " << G.header.n_columns << " by " << G.header.n_rows << ") - abort." << std::endl;
            std::exit(2);
        }

        /* Determine the initial and intermediate grid dimensions */
        C.current_stride = gcd_euclid(C.n_columns - 1, C.n_rows - 1);

        /* Set current_stride = 1, read data, setting indices.  Then throw
           away data that can't be used in the end game, limiting the
           size of data arrays and Briggs->b[6] structure/array.  */

        C.current_stride = 1;
        set_grid_parameters(C, G.header);

        setup_data(C, G.header, x, y, z, n_nodes);

        cull_data(C); /* Eliminate data points that will not serve as constraints */

        surface_remove_planar_trend(C, G.header); /* Fit best-fitting plane and remove it from the data; plane will be restored at the end */
        bool is_flat = rescale_z_values(C);       /* Divide residual data by their rms value */

        if (is_flat)
        { /* Data lie exactly on a plane; write a grid with the plane and exit */ /* Get a grid of zeros... */
            surface_restore_planar_trend(C, G);
            store_grid(C, G, output, nx, ny);
            return;
        }

        load_constraints(C, G.header, lower, n_lower, upper, n_upper, nx, ny);

        /* Set up factors and reset current_stride to its initial (and largest) value  */
        G.data.resize(C.mxmy, 0.0);
        C.status.resize(C.mxmy, NodeStatus::IS_UNCONSTRAINED);
        C.current_stride = gcd_euclid(C.n_columns - 1, C.n_rows - 1);
        C.factors = get_prime_factors(C.current_stride);

        set_grid_parameters(C, G.header);
        while (C.current_nx < 4 || C.current_ny < 4)
        { /* Must have at least a grid of 4x4 */
            smart_divide(C);
            set_grid_parameters(C, G.header);
        }
        set_offset(C);          /* Initialize the node-jumps across rows for this grid size */
        set_index(C, G.header); /* Determine the nearest data constraint for this grid size */

        /* Now the data are ready to go for the first iteration.  */

        C.brgs_coefs.resize(C.points.size());

        set_coefficients(C); /* Initialize the coefficients needed in the finite-difference expressions */

        /* Here is the main multigrid loop, were we first grid using a coarse grid and the
         * progressively refine the grid until we reach the final configuration. */

        C.previous_stride = C.current_stride;
        find_nearest_constraint(C, G);                                       /* Assign nearest data value to nodes and evaluate Briggs coefficients */
        C.total_iterations += relaxation_iterate(C, G, IterMode::GRID_DATA); /* Grid the data using the data constraints */

        while (C.current_stride > 1)
        {                                                                         /* More intermediate grids remain, go to next */
            smart_divide(C);                                                      /* Set the new current_stride */
            set_grid_parameters(C, G.header);                                     /* Update node book-keeping constants */
            set_offset(C);                                                        /* Reset the node-jumps across rows for this grid size */
            set_index(C, G.header);                                               /* Recompute the index values for the nearest data points */
            fill_in_forecast(C, G);                                               /* Expand the grid and fill it via bilinear interpolation */
            C.total_iterations += relaxation_iterate(C, G, IterMode::GRID_NODES); /* Grid again but only to improve on the bilinear guesses */
            find_nearest_constraint(C, G);                                        /* Assign nearest data value to nodes and evaluate Briggs coefficients */
            C.total_iterations += relaxation_iterate(C, G, IterMode::GRID_DATA);  /* Grid the data but now use the data constraints */
            C.previous_stride = C.current_stride;                                 /* Remember previous stride before we smart-divide again */
        }

        surface_check_errors(C, G); /* Report on mean misfit and curvature */

        surface_restore_planar_trend(C, G); /* Restore the least-square plane we removed earlier */

        store_grid(C, G, output, nx, ny);
    }
}

extern "C"
{
    void minimum_curvature(double *x, double *y, double *z, size_t n,
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
                           unsigned char verbosity)
    {
        size_t max_iters = (max_iterations <= 0) ? 500 : max_iterations;
        double rlx = (relax <= 0.0) ? SURFACE_OVERRELAXATION : std::max(std::min(2.0, relax), 0.0);
        double alph = (alpha <= 0.0) ? 1.0 : alpha;
        double b_tens = std::min(std::max(0.0, b_tension), 1.0);
        double i_tens = std::min(std::max(0.0, i_tension), 1.0);
        double cnvg = (converge_limit <= 0.0) ? SURFACE_CONV_LIMIT : converge_limit;
        auto d_mode = surface::DownsampleMode{std::min(downsample_mode, 2U)};
        auto cnvg_mode = surface::ConvergenceLimitMode{std::min(converge_mode, 1U)};

        try
        {
            surface::surface(
                x, y, z, n, xmin, xmax, ymin, ymax, nx, ny, output,
                d_mode, max_iters, rlx, alph, b_tens, i_tens, cnvg, cnvg_mode, verbosity, lower, n_lower, upper, n_upper);
        }
        catch (std::exception &e)
        {
            std::cerr << e.what() << std::endl;
            std::exit(1);
        }
    }
}