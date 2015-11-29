/*!
************************************************************************
* \file adaptive_filter.c
*
* \brief
*    adaptive filter routines
*
* \author
*  Yuri Vatis               <vatis@tnt.uni-hannover.de>    \n
*
* This software may be used for research purposes only.
* Universitaet Hannover would also appreciate that all those using the 
* software or any extension of it in any way comply with the 
* following condition:
*    Including in any technical report, conference or journal paper
*    produced which uses the software or any extension of it in
*    any way, at least one of bibliographic references:
*    Yuri Vatis, Bernd Edler, Dieu Thanh Nguyen, Jörn Ostermann,
*    "Motion-and Aliasing-compensated Prediction using a two-dimensional non-separable Adaptive Wiener Interpolation Filter",
*    Proc. ICIP 2005, IEEE International Conference on Image Processing, Genova, Italy, September 2005.
*    Yuri Vatis, Bernd Edler, Ingolf Wassermann, Dieu Thanh Nguyen, Jörn Ostermann,
*    "Coding of Coefficients of two-dimensional non-separable Adaptive Wiener Interpolation Filter", 
*    Proc. VCIP 2005, SPIE Visual Communication & Image Processing, Beijing, China, July 2005.
************************************************************************
*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "defines.h"
#include "adaptive_filter.h"

#ifdef ADAPTIVE_FILTER
#include "global.h" 
#ifdef DIRECTIONAL_FILTER
#include <memory.h>
#include "../../lcommon/inc/adaptive_filter_1DAIF.h"
#include "../../lcommon/inc/adaptive_filter_orig.h"
#ifdef EDAIF2
#include "../../lcommon/inc/const_DAIF.h"
#include "../../lcommon/inc/const_RAIF.h"
#endif
#endif
#ifdef EIGHTH_PEL
#include "image.h"
#endif

double STANDARD_2D_FILTER[15][SQR_FILTER] = {
  {   0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  1.0/64.0  ,   -5.0/64.0  ,   52.0/64.0  ,   20.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0
  }, // a_pos
  {   0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  1.0/32.0  ,   -5.0/32.0  ,   20.0/32.0  ,   20.0/32.0  ,   -5.0/32.0  ,   1.0/32.0  ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0
  }, // b_pos
  {   0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  1.0/64.0  ,   -5.0/64.0  ,   20.0/64.0  ,   52.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0
  }, // c_pos
  {   0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   52.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0
  }, // d_pos
  {   0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  1.0/64.0  ,   -5.0/64.0  ,   40.0/64.0  ,   20.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
  0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0
  }, // e_pos
  {   1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0,
  -5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
  52.0/2048.0, -260.0/2048.0, 1040.0/2048.0, 1040.0/2048.0, -260.0/2048.0,  52.0/2048.0,
  20.0/2048.0, -100.0/2048.0,  400.0/2048.0,  400.0/2048.0, -100.0/2048.0,  20.0/2048.0,
  -5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
  1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0
  }, // f_pos
  {   0.0       ,    0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,   0.0       ,
  1.0/64.0  ,   -5.0/64.0  ,   20.0/64.0  ,   40.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
  0.0       ,    0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,   0.0
  }, // g_pos
  {   0.0       ,    0.0       ,    1.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   20.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   20.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    1.0/32.0  ,    0.0       ,    0.0       ,   0.0
  }, // h_pos
  {   1.0/2048.0,   -5.0/2048.0,   52.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0,
  -5.0/2048.0,   25.0/2048.0, -260.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
  20.0/2048.0, -100.0/2048.0, 1040.0/2048.0,  400.0/2048.0, -100.0/2048.0,  20.0/2048.0,
  20.0/2048.0, -100.0/2048.0, 1040.0/2048.0,  400.0/2048.0, -100.0/2048.0,  20.0/2048.0,
  -5.0/2048.0,   25.0/2048.0, -260.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
  1.0/2048.0,   -5.0/2048.0,   52.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0
  }, // i_pos
  {   1.0/1024.0,   -5.0/1024.0,   20.0/1024.0,   20.0/1024.0,   -5.0/1024.0,   1.0/1024.0,
  -5.0/1024.0,   25.0/1024.0, -100.0/1024.0, -100.0/1024.0,   25.0/1024.0,  -5.0/1024.0,
  20.0/1024.0, -100.0/1024.0,  400.0/1024.0,  400.0/1024.0, -100.0/1024.0,  20.0/1024.0,
  20.0/1024.0, -100.0/1024.0,  400.0/1024.0,  400.0/1024.0, -100.0/1024.0,  20.0/1024.0,
  -5.0/1024.0,   25.0/1024.0, -100.0/1024.0, -100.0/1024.0,   25.0/1024.0,  -5.0/1024.0,
  1.0/1024.0,   -5.0/1024.0,   20.0/1024.0,   20.0/1024.0,   -5.0/1024.0,   1.0/1024.0
  }, // j_pos
  {   1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   52.0/2048.0,   -5.0/2048.0,   1.0/2048.0,
  -5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -260.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
  20.0/2048.0, -100.0/2048.0,  400.0/2048.0, 1040.0/2048.0, -100.0/2048.0,  20.0/2048.0,
  20.0/2048.0, -100.0/2048.0,  400.0/2048.0, 1040.0/2048.0, -100.0/2048.0,  20.0/2048.0,
  -5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -260.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
  1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   52.0/2048.0,   -5.0/2048.0,   1.0/2048.0
  }, // i_pos
  {   0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   52.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0
  }, // l_pos
  {   0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  1.0/64.0  ,   -5.0/64.0  ,   40.0/64.0  ,   20.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
  0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0
  }, // m_pos
  {   1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0,
  -5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
  20.0/2048.0, -100.0/2048.0,  400.0/2048.0,  400.0/2048.0, -100.0/2048.0,  20.0/2048.0,
  52.0/2048.0, -260.0/2048.0, 1040.0/2048.0, 1040.0/2048.0, -260.0/2048.0,  52.0/2048.0,
  -5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
  1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0
  }, // n_pos
  {   0.0         ,    0.0         ,  0.0     ,    1.0/64.0  ,    0.0       ,   0.0       ,
  0.0         ,    0.0         ,  0.0     ,   -5.0/64.0  ,    0.0       ,   0.0       ,
  0.0         ,    0.0         ,  0.0     ,   20.0/64.0  ,    0.0       ,   0.0       ,
  1.0/64.0    ,   -5.0/64.0    , 20.0/64.0,   40.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
  0.0         ,    0.0         ,  0.0     ,   -5.0/64.0  ,    0.0       ,   0.0       ,
  0.0         ,    0.0         ,  0.0     ,    1.0/64.0  ,    0.0       ,   0.0
  } // o_pos
};
// separable aif (BEGIN)
static double STANDARD_2D_FILTER_SEP[15][SQR_FILTER] = {
  {1.0 / 64.0, -5.0 / 64.0, 52.0 / 64.0, 20.0 / 64.0, -5.0 / 64.0, 1.0 / 64.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},       // a_pos
  {1.0 / 32.0, -5.0 / 32.0, 20.0 / 32.0, 20.0 / 32.0, -5.0 / 32.0, 1.0 / 32.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},       // b_pos
  {1.0 / 64.0, -5.0 / 64.0, 20.0 / 64.0, 52.0 / 64.0, -5.0 / 64.0, 1.0 / 64.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},       // c_pos
  {1.0 / 64.0, -5.0 / 64.0, 52.0 / 64.0, 20.0 / 64.0, -5.0 / 64.0, 1.0 / 64.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},       // d_pos
  {1.0 / 64.0, -5.0 / 64.0, 52.0 / 64.0, 20.0 / 64.0, -5.0 / 64.0, 1.0 / 64.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},       // e_pos
  {1.0 / 64.0, -5.0 / 64.0, 52.0 / 64.0, 20.0 / 64.0, -5.0 / 64.0, 1.0 / 64.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},       // f_pos
  {1.0 / 64.0, -5.0 / 64.0, 52.0 / 64.0, 20.0 / 64.0, -5.0 / 64.0, 1.0 / 64.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},       // g_pos
  {1.0 / 32.0, -5.0 / 32.0, 20.0 / 32.0, 20.0 / 32.0, -5.0 / 32.0, 1.0 / 32.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},       // h_pos
  {1.0 / 32.0, -5.0 / 32.0, 20.0 / 32.0, 20.0 / 32.0, -5.0 / 32.0, 1.0 / 32.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},       // i_pos
  {1.0 / 32.0, -5.0 / 32.0, 20.0 / 32.0, 20.0 / 32.0, -5.0 / 32.0, 1.0 / 32.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},       // j_pos
  {1.0 / 32.0, -5.0 / 32.0, 20.0 / 32.0, 20.0 / 32.0, -5.0 / 32.0, 1.0 / 32.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},       // k_pos
  {1.0 / 64.0, -5.0 / 64.0, 20.0 / 64.0, 52.0 / 64.0, -5.0 / 64.0, 1.0 / 64.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},       // l_pos
  {1.0 / 64.0, -5.0 / 64.0, 20.0 / 64.0, 52.0 / 64.0, -5.0 / 64.0, 1.0 / 64.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},       // m_pos
  {1.0 / 64.0, -5.0 / 64.0, 20.0 / 64.0, 52.0 / 64.0, -5.0 / 64.0, 1.0 / 64.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},       // n_pos
  {1.0 / 64.0, -5.0 / 64.0, 20.0 / 64.0, 52.0 / 64.0, -5.0 / 64.0, 1.0 / 64.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},       // o_pos
};
// separable aif (END)
int TwoDSymmetricPattern[15][SQR_FILTER]  =   // get one filter from another one, if symmetry properties are used (used e.g. in ExtendFilterCoefficients)
{
  { 0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,
  12, 13, 14, 15, 16, 17,        // a_pos  - not used actually
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0  },
  { 0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,
  12, 13, 14, 14, 13, 12,        // b_pos
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0  },
  { 0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,
  17, 16, 15, 14, 13, 12,        // c_pos
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0  },
  { 0,  0, 12,  0,  0,  0,
  0,  0, 13,  0,  0,  0,
  0,  0, 14,  0,  0,  0,        // d_pos
  0,  0, 15,  0,  0,  0,
  0,  0, 16,  0,  0,  0,
  0,  0, 17,  0,  0,  0  },
  { 0,  1,  2,  3,  4,  5,
  1,  7,  8,  9, 10, 11,
  2,  8, 14, 15, 16, 17,        // e_pos
  3,  9, 15, 21, 22, 23,
  4, 10, 16, 22, 28, 29,
  5, 11, 17, 23, 29, 35  },
  { 0,  1,  2,  2,  1,  0,
  6,  7,  8,  8,  7,  6,
  12, 13, 14, 14, 13, 12,        // f_pos
  18, 19, 20, 20, 19, 18,
  24, 25, 26, 26, 25, 24,
  30, 31, 32, 32, 31, 30  },
  { 5,  4,  3,  2,  1,  0,
  11, 10,  9,  8,  7,  1,
  17, 16, 15, 14,  8,  2,        // g_pos
  23, 22, 21, 15,  9,  3,
  29, 28, 22, 16, 10,  4,
  35, 29, 23, 17, 11,  5  },
  { 0,  0, 12,  0,  0,  0,
  0,  0, 13,  0,  0,  0,
  0,  0, 14,  0,  0,  0,        // h_pos
  0,  0, 14,  0,  0,  0,
  0,  0, 13,  0,  0,  0,
  0,  0, 12,  0,  0,  0  },
  { 0,  6, 12, 18, 24, 30,
  1,  7, 13, 19, 25, 31,
  2,  8, 14, 20, 26, 32,        // i_pos
  2,  8, 14, 20, 26, 32,
  1,  7, 13, 19, 25, 31,
  0,  6, 12, 18, 24, 30  },
  { 0,  1,  2,  2,  1,  0,
  1,  7,  8,  8,  7,  1,
  2,  8, 14, 14,  8,  2,         // j_pos
  2,  8, 14, 14,  8,  2,
  1,  7,  8,  8,  7,  1,
  0,  1,  2,  2,  1,  0  },
  {30, 24, 18, 12,  6,  0,
  31, 25, 19, 13,  7,  1,
  32, 26, 20, 14,  8,  2,         // k_pos
  32, 26, 20, 14,  8,  2,
  31, 25, 19, 13,  7,  1,
  30, 24, 18, 12,  6,  0  },
  { 0,  0, 17,  0,  0,  0,
  0,  0, 16,  0,  0,  0,
  0,  0, 15,  0,  0,  0,         // l_pos
  0,  0, 14,  0,  0,  0,
  0,  0, 13,  0,  0,  0,
  0,  0, 12,  0,  0,  0  },
  { 5, 11, 17, 23, 29, 35,
  4, 10, 16, 22, 28, 29,
  3,  9, 15, 21, 22, 23,        // m_pos
  2,  8, 14, 15, 16, 17,
  1,  7,  8,  9, 10, 11,
  0,  1,  2,  3,  4,  5  },
  {30, 31, 32, 32, 31, 30,
  24, 25, 26, 26, 25, 24,
  18, 19, 20, 20, 19, 18,         // n_pos
  12, 13, 14, 14, 13, 12,
  6,  7,  8,  8,  7,  6,
  0,  1,  2,  2,  1,  0  },
  {35, 29, 23, 17, 11,  5,
  29, 28, 22, 16, 10,  4,
  23, 22, 21, 15,  9,  3,         // o_pos
  17, 16, 15, 14,  8,  2,
  11, 10,  9,  8,  7,  1,
  5,  4,  3,  2,  1,  0  }
};
int FILTER_NEW_SUB_POS[15] = {       a_pos, b_pos, a_pos,
a_pos, e_pos, f_pos, e_pos,
b_pos, f_pos, j_pos, f_pos,
a_pos, e_pos, f_pos, e_pos };
//static int SymmetryPosition[15] = {1,1,0,0,1,1,0,0,0,1,0,0,0,0};

#ifdef E_DAIF
int numQBitsInt[SQR_FILTER_INT-1] = 
{
  12, 12, 12, 12, 12, 
  12, 11, 11, 11, 12, 
  12, 11, 10, 11, 12, 
  12, 11, 11, 11, 12,
  12, 12, 12, 12, 12, 
};

int DiffQFilterOffsetI[15], DiffQFilterOffsetF[15], DiffQFilterOffsetSign[15];
int DiffQFilterOffsetIntI, DiffQFilterOffsetIntF, DiffQFilterOffsetIntSign; 
int DiffQFilterCoeffInt[SQR_FILTER_INT]; 
int FilterFlagInt;
int FILTER_OFFSET_INT = 2; 
double FilterCoefInt[SQR_FILTER_INT];
#endif

#ifdef EAIF
// Number of quantization bits 
int numQBits1DH[SQR_FILTER-1] = 
{
   0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,
  12, 11,  9,  9, 11, 12,
   0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,
};

int numQBits1DV[SQR_FILTER-1] = 
{
   0,  0, 12,  0,  0,  0,
   0,  0, 11,  0,  0,  0,
   0,  0,  9,  0,  0,  0,
   0,  0,  9,  0,  0,  0,
   0,  0, 11,  0,  0,  0,
   0,  0, 12,  0,  0,  0,
};

int numQBits2D[SQR_FILTER-1] = 
{
   0,  0,  0,  0,  0,  0,
   0,  0, 10, 10,  0,  0,
   0, 10,  9,  9, 10,  0,
   0, 10,  9,  9, 10,  0,
   0,  0, 10, 10,  0,  0,
   0,  0,  0,  0,  0,  0,
};
#endif
#ifdef EDAIF2  // Variable Structure AIF
int numQBitsDIF[SQR_FILTER-1] = 
{
   12, 0,  0,  0,  0,  12,
   0, 11,  0, 0,   11,  0,
   0,  0,  9,  9,   0,  0,
   0,  0,  9,  9,   0,  0,
   0, 11,  0,  0,  11,  0,
  12,  0,  0,  0,  0,  12,
};
#endif

#ifndef EDAIF2
int FILTER_SIZE = 6;
int FILTER_OFFSET;
#endif
int UseAllSubpelPositions;
int SubpelPositionsPattern;
int DiffQFilterCoef[15][SQR_FILTER];             // differences to be received
int FilterFlag[15];                              // Flags defining if Filter at the particular position calculated
int POS_EQUATION_NUMBER[15];                     // number of different coefficietns for each sub-pel position
int SymmetryPosition[15];
int NumberOfQBits;                               // number of bits used for quantization of
static double FilterCoef[15][SQR_FILTER];                // SQR(PEL)-1 X SQR(FILTER_SIZE) - Filter Coeffitients
short int FilterCoef_16bits[15][SQR_FILTER];             // AIF in integer representation
short int PredFilterCoef16bits[15][SQR_FILTER];          // AVC in integer representation
static int    QFilterCoef[15][SQR_FILTER];       // SQR(PEL)-1 X SQR(FILTER_SIZE) - quantized filter coefficients
//  static double LastFilterCoef[15][SQR_FILTER];  // filter used for previous frame
static double PredFilterCoef[15][SQR_FILTER];    // predicted filter coefficients
static int    QPredFilterCoef[15][SQR_FILTER];   // quantized predicted filter coefficients

#ifdef DIRECTIONAL_FILTER
int IsDiagonal1D[15];
int nQBits[15] = {7,7,7,7,7,8,7,7,8,8,8,7,7,8,7};

#ifdef EDAIF2
int RealTimeDecomposition;
double CalculatedFilter[MAX_NUM_AIF][SQR_FILTER];
int realSymCommands[4][MAX_NUM_AIF];
int Calc2Filt_Indexes[MAX_NUM_AIF][SQR_FILTER];
int Filt_Indexes_offset[MAX_NUM_AIF];
int nBitsIntRepresent[MAX_NUM_AIF];
int nPairsSets[4] = {4,4,6,6};
int SetHor[4][2] = 
    {
      {a_pos, c_pos}, 
      {e_pos, g_pos}, 
      {i_pos, k_pos}, 
      {m_pos, o_pos}, 
    };
int SetVert[4][2] = 
    {
      {d_pos, l_pos}, 
      {e_pos, m_pos}, 
      {f_pos, n_pos}, 
      {g_pos, o_pos}, 
    };

int SetDiagNW[6][2] = 
    {
      {a_pos, l_pos}, 
      {b_pos, h_pos}, 
      {c_pos, d_pos}, 
      {e_pos, o_pos}, 
      {f_pos, k_pos},
      {i_pos, n_pos}
    };

int SetDiagNE[6][2] = 
    {
      {a_pos, d_pos}, 
      {b_pos, h_pos}, 
      {c_pos, l_pos}, 
      {g_pos, m_pos}, 
      {f_pos, i_pos},
      {k_pos, n_pos},
    };

#endif


void initFilterCustom(int filterID)
{
#ifdef EDAIF2
    UseAllSubpelPositions = 0;
    SubpelPositionsPattern = 0;
    RealTimeDecomposition = 0;
    memset(realSymCommands,0,sizeof(int)*4*MAX_NUM_AIF);
#endif
    
    memset(FilterFlag,0,MAX_NUM_AIF*sizeof(int));
    memset(DiffQFilterCoef, 0, SQR_FILTER*MAX_NUM_AIF*sizeof(int));
    memset(CalculatedFilter, 0, SQR_FILTER*MAX_NUM_AIF*sizeof(double));
    memset(FilterCoef, 0, SQR_FILTER*MAX_NUM_AIF*sizeof(double));


  if (img->ImpType == IMP_INT16)
    FILTCOEF_BITS = 7;
  else 
    FILTCOEF_BITS = DEFAULT_QUANT;
  if (filterID==FILTER_TYPE_2D_NS)
  {
    memcpy(&STANDARD_2D_FILTER[0][0],&STANDARD_2D_FILTER_orig[0][0],sizeof(double)*15*SQR_FILTER);
    memcpy(&TwoDSymmetricPattern[0][0],&TwoDSymmetricPattern_orig[0][0],sizeof(int)*15*SQR_FILTER);

    memcpy(&IsDiagonal1D[0],&IsDiagonal1D_orig[0],sizeof(int)*15);
    memcpy(&POS_EQUATION_NUMBER[0],&POS_EQUATION_NUMBER_orig[0],sizeof(int)*15);
    memcpy(&FILTER_NEW_SUB_POS[0],&FILTER_NEW_SUB_POS_orig[0],sizeof(int)*15);

  }
#ifdef E_DAIF
  else if (filterID==FILTER_TYPE_1D || filterID == FILTER_TYPE_EDAIF)
#else
  else if (filterID==FILTER_TYPE_1D)
#endif
  {
#ifndef EDAIF2
    NumberOfQBits = FILTCOEF_BITS;
    memcpy(TwoDSymmetricPattern,TwoDSymmetricPattern_v4,sizeof(int)*15*SQR_FILTER);
    memcpy(IsDiagonal1D,IsDiagonal1D_v4,sizeof(int)*15);
    memcpy(POS_EQUATION_NUMBER,POS_EQUATION_NUMBER_v4,sizeof(int)*15);
    memcpy(FILTER_NEW_SUB_POS,FILTER_NEW_SUB_POS_v4,sizeof(int)*15);
    memcpy(&SymmetryPosition[0],&SymmetryPosition_v4[0],sizeof(int)*15);
    memcpy(&STANDARD_2D_FILTER[0][0],STANDARD_2D_FILTER_v4,sizeof(double)*15*SQR_FILTER);
#else
    if (img->ImpType == IMP_INT16)
      memcpy(&nBitsIntRepresent[0],&nBitsIntRepresent_DIAF[0],sizeof(int)*MAX_NUM_AIF);
    memcpy(&STANDARD_2D_FILTER[0][0],&STANDARD_2D_FILTER_DAIF[0][0],sizeof(double)*MAX_NUM_AIF*SQR_FILTER);


    memcpy(&IsDiagonal1D[0],&IsDiagonal1D_DAIF[0],sizeof(int)*MAX_NUM_AIF);
    memcpy(&SymmetryPosition[0],&SymmetryPosition_DAIF[0],sizeof(int)*MAX_NUM_AIF);
    memcpy(&POS_EQUATION_NUMBER[0],&POS_EQUATION_NUMBER_DAIF[0],sizeof(int)*MAX_NUM_AIF);
    memcpy(&FILTER_NEW_SUB_POS[0],&FILTER_NEW_SUB_POS_DAIF[0],sizeof(int)*MAX_NUM_AIF);

    memcpy(&Filt_Indexes_offset[0],&Filt_Indexes_offset_DAIF[0],sizeof(int)*MAX_NUM_AIF);
    memcpy(&Calc2Filt_Indexes[0][0],&Calc2Filt_Indexes_DAIF[0][0],sizeof(int)*MAX_NUM_AIF*SQR_FILTER);
    if (img->ImpType != IMP_INT16)
    {
      int k;
      for(k = 0; k < MAX_NUM_AIF; k++)
        nBitsIntRepresent[k]=DEFAULT_QUANT;
    }
#endif
  }
#ifdef EAIF
  else if (filterID == FILTER_TYPE_EAIF)
  {
    NumberOfQBits = FILTCOEF_BITS;
    memcpy(&TwoDSymmetricPattern[0][0],&TwoDSymmetricPattern_v4[0][0],sizeof(int)*15*SQR_FILTER);
    memcpy(&IsDiagonal1D[0],&IsDiagonal1D_v4[0],sizeof(int)*15);
    memcpy(&POS_EQUATION_NUMBER[0],&POS_EQUATION_NUMBER_EAIF[0],sizeof(int)*15);
    memcpy(&FILTER_NEW_SUB_POS[0],&FILTER_NEW_SUB_POS_EAIF[0],sizeof(int)*15);
    memcpy(&SymmetryPosition[0],&SymmetryPosition_EAIF[0],sizeof(int)*15);
    memcpy(&STANDARD_2D_FILTER[0][0],&STANDARD_2D_FILTER_orig[0][0],sizeof(double)*15*SQR_FILTER);
  }
#endif
}
#endif
/*!
************************************************************************
* \brief
*    init adaptive filter structure
************************************************************************
*/
void InitAdaptiveFilter(void)
{
#ifndef EDAIF2
  FILTER_OFFSET = (FILTER_SIZE)/2-1;
#endif

  POS_EQUATION_NUMBER[a_pos] =  6;
  POS_EQUATION_NUMBER[b_pos] =  3;
  POS_EQUATION_NUMBER[c_pos] =  6;

  POS_EQUATION_NUMBER[d_pos] =  6;
  POS_EQUATION_NUMBER[e_pos] = 21;
  POS_EQUATION_NUMBER[f_pos] = 18;
  POS_EQUATION_NUMBER[g_pos] = 21;

  POS_EQUATION_NUMBER[h_pos] =  3;
  POS_EQUATION_NUMBER[i_pos] = 18;
  POS_EQUATION_NUMBER[j_pos] =  6;
  POS_EQUATION_NUMBER[k_pos] = 18;

  POS_EQUATION_NUMBER[l_pos] =  6;
  POS_EQUATION_NUMBER[m_pos] = 21;
  POS_EQUATION_NUMBER[n_pos] = 18;
  POS_EQUATION_NUMBER[o_pos] = 21;
}


/*!
************************************************************************
* \brief
*    free adaptive filter structure
************************************************************************
*/
void FreeAdaptiveFilter(void)
{
  ;
}

/*
************************************************************************
* \brief
*    reset adaptive filter structure
************************************************************************
*/
void ResetAdaptiveFilter(void)
{
  int sub_pel,i;

  for(sub_pel = 0; sub_pel < 15; sub_pel++)
    for(i = 0; i < SQR_FILTER; i++)
    {
      QPredFilterCoef[sub_pel][i] = 0;
      DiffQFilterCoef[sub_pel][i] = 0;
      QFilterCoef    [sub_pel][i] = 0;
      PredFilterCoef [sub_pel][i] = 0.0;
    }
    NumberOfQBits = DEFAULT_QUANT;
#ifdef DIRECTIONAL_FILTER
    if (img->ImpType==IMP_INT16)
      FILTCOEF_BITS = 7;
    else
      FILTCOEF_BITS = DEFAULT_QUANT;
#endif
}

/*
************************************************************************
* \brief
*    return if adaptive filter is used for current frame
************************************************************************
*/
int UseAdaptiveFilterForCurrentFrame(void)
{
  if( (img->type == P_SLICE || img->type == B_SLICE) && img->UseAdaptiveFilter)
    return 1;
  else
    return 0;
}

#ifdef E_DAIF

static void ZeroFilterOffsets()
{
  int sub_pos;

  for(sub_pos = 0; sub_pos < 15; sub_pos++)
  {
    FilterCoef[sub_pos][SQR_FILTER-1] = 0.;
  }
}
#endif

/*!
************************************************************************
* \brief
*    Predict and calculate filter coefficients
************************************************************************
*/
void CalculateFilterCoefficients(void)
{
  int i,j, sub_pel;

  // at first predict b and h position using inter prediction, b and h are already quantized in last frame
  if(img->AdaptiveFilterFlag)
  {
    if(FilterFlag[b_pos])
    {
      for(i=0; i < SQR_FILTER; i++)
      {
        PredFilterCoef[b_pos][i] = STANDARD_2D_FILTER[b_pos][i];
        PredFilterCoef[h_pos][i] = STANDARD_2D_FILTER[h_pos][i];
      }
      DequantizeFilterCoefficients(DiffQFilterCoef[b_pos], FilterCoef[b_pos]);
      // add differences to predicted values
      for(i=FILTER_SIZE*FILTER_OFFSET; i < FILTER_SIZE*(FILTER_OFFSET+1); i++)
      {
        FilterCoef[b_pos][i] += PredFilterCoef[b_pos][i];
      }
    }
    else
    {
      for(i=0; i < SQR_FILTER; i++)
        FilterCoef[b_pos][i] = STANDARD_2D_FILTER[b_pos][i];
    }
    ExtendFilterCoefficientsFloat(h_pos,FilterCoef);

    // predict a_pos:
    if(FilterFlag[a_pos])
    {
      for(i=0; i < SQR_FILTER; i++)
      {
        PredFilterCoef[a_pos][i] = STANDARD_2D_FILTER[a_pos][i];
      }
      DequantizeFilterCoefficients(DiffQFilterCoef[a_pos], FilterCoef[a_pos]);
      // add differences to predicted values
      for(i=FILTER_SIZE*FILTER_OFFSET; i < FILTER_SIZE*(FILTER_OFFSET+1); i++)
      {
        FilterCoef[a_pos][i] += PredFilterCoef[a_pos][i];
      }
    }
    else
    {
      for(i=0; i < SQR_FILTER; i++)
        FilterCoef[a_pos][i] = STANDARD_2D_FILTER[a_pos][i];
    }
    ExtendFilterCoefficientsFloat(c_pos,FilterCoef);
    ExtendFilterCoefficientsFloat(d_pos,FilterCoef);
    ExtendFilterCoefficientsFloat(l_pos,FilterCoef);

    //predict e_pos
    if(FilterFlag[e_pos])
    {
      for(i = 0; i < FILTER_SIZE; i++)
      {
        for(j = 0; j < FILTER_SIZE; j++)
        {
          PredFilterCoef[e_pos][FILTER_SIZE*i+j] =
            FilterCoef[d_pos][FILTER_SIZE*i+FILTER_OFFSET]*FilterCoef[a_pos][FILTER_SIZE*FILTER_OFFSET+j];
        }
      }
      DequantizeFilterCoefficients(DiffQFilterCoef[e_pos], FilterCoef[e_pos]);
      for(i = 0; i < SQR_FILTER; i++)
      {
        FilterCoef[e_pos][i] += PredFilterCoef[e_pos][i];
      }
    }
    else
    {
      for(i=0; i < SQR_FILTER; i++)
        FilterCoef[e_pos][i] = STANDARD_2D_FILTER[e_pos][i];
    }
    ExtendFilterCoefficientsFloat(g_pos,FilterCoef);
    ExtendFilterCoefficientsFloat(m_pos,FilterCoef);
    ExtendFilterCoefficientsFloat(o_pos,FilterCoef);

    //predict f_pos:
    if(FilterFlag[f_pos])
    {
      for(i = 0; i < FILTER_SIZE; i++)
      {
        for(j = 0; j < FILTER_SIZE; j++)
        {
          PredFilterCoef[f_pos][FILTER_SIZE*i+j] =
            FilterCoef[d_pos][FILTER_SIZE*i+FILTER_OFFSET]*FilterCoef[b_pos][FILTER_SIZE*FILTER_OFFSET+j];
        }
      }
      DequantizeFilterCoefficients(DiffQFilterCoef[f_pos], FilterCoef[f_pos]);
      for(i = 0; i < SQR_FILTER; i++)
      {
        FilterCoef[f_pos][i] += PredFilterCoef[f_pos][i];
      }
    }
    else
    {
      for(i=0; i < SQR_FILTER; i++)
        FilterCoef[f_pos][i] = STANDARD_2D_FILTER[f_pos][i];
    }
    ExtendFilterCoefficientsFloat(i_pos,FilterCoef);
    ExtendFilterCoefficientsFloat(k_pos,FilterCoef);
    ExtendFilterCoefficientsFloat(n_pos,FilterCoef);

    //predict j_pos:
    if(FilterFlag[j_pos])
    {
      for(i = 0; i < FILTER_SIZE; i++)
      {
        for(j = 0; j < FILTER_SIZE; j++)
        {
          PredFilterCoef[j_pos][FILTER_SIZE*i+j] =
            FilterCoef[h_pos][FILTER_SIZE*i+FILTER_OFFSET]*FilterCoef[b_pos][FILTER_SIZE*FILTER_OFFSET+j];
        }
      }
      DequantizeFilterCoefficients(DiffQFilterCoef[j_pos], FilterCoef[j_pos]);
      for(i = 0; i < SQR_FILTER; i++)
      {
        FilterCoef[j_pos][i] += PredFilterCoef[j_pos][i];
      }
    }
    else
    {
      for(i=0; i < SQR_FILTER; i++)
        FilterCoef[j_pos][i] = STANDARD_2D_FILTER[j_pos][i];
    }
  }
  else
    for(sub_pel = 0; sub_pel < 15; sub_pel++)
      for(i = 0; i < SQR_FILTER; i++)
        FilterCoef[sub_pel][i] = STANDARD_2D_FILTER[sub_pel][i];
  //  PrintFilterCoefInt(DiffQFilterCoef);
  //  PrintFilterCoefFloat(PredFilterCoef);
  //  PrintFilterCoefFloat(FilterCoef);

#ifdef E_DAIF
  ZeroFilterOffsets();
#endif
}

#ifdef DIRECTIONAL_FILTER
/*!
************************************************************************
* \brief
*    Predict and calculate filter coefficients
************************************************************************
*/
void CalculateFilterCoefficients1DAIF(void)
{
  int i,j;
  int multplier;
  double abs_value; 
  int signFC;

  // at first predict b and h position using inter prediction, b and h are already quantized in last frame
  if(img->AdaptiveFilterFlag)
  {
    if (img->ImpType == IMP_FLOAT32)
    {
      if(FilterFlag[b_pos])
      {
        for(i=0; i < SQR_FILTER; i++)
        {
          PredFilterCoef[b_pos][i] = STANDARD_2D_FILTER[b_pos][i];
          PredFilterCoef[h_pos][i] = STANDARD_2D_FILTER[h_pos][i];
        }
        DequantizeFilterCoefficients(DiffQFilterCoef[b_pos], FilterCoef[b_pos]);
        // add differences to predicted values
        for(i=FILTER_SIZE*FILTER_OFFSET; i < FILTER_SIZE*(FILTER_OFFSET+1); i++)
        {
          FilterCoef[b_pos][i] += PredFilterCoef[b_pos][i];
        }
      }
      else
      {
        for(i=0; i < SQR_FILTER; i++)
          FilterCoef[b_pos][i] = STANDARD_2D_FILTER[b_pos][i];
      }
      ExtendFilterCoefficientsFloat(h_pos,FilterCoef);

      // predict a_pos:
      if(FilterFlag[a_pos])
      {
        for(i=0; i < SQR_FILTER; i++)
        {
          PredFilterCoef[a_pos][i] = STANDARD_2D_FILTER[a_pos][i];
        }
        DequantizeFilterCoefficients(DiffQFilterCoef[a_pos], FilterCoef[a_pos]);
        // add differences to predicted values
        for(i=FILTER_SIZE*FILTER_OFFSET; i < FILTER_SIZE*(FILTER_OFFSET+1); i++)
        {
          FilterCoef[a_pos][i] += PredFilterCoef[a_pos][i];
        }
      }
      else
      {
        for(i=0; i < SQR_FILTER; i++)
          FilterCoef[a_pos][i] = STANDARD_2D_FILTER[a_pos][i];
      }
      ExtendFilterCoefficientsFloat(c_pos,FilterCoef);
      ExtendFilterCoefficientsFloat(d_pos,FilterCoef);
      ExtendFilterCoefficientsFloat(l_pos,FilterCoef);
      //predict e_pos
      for(i=0; i < SQR_FILTER; i++)
      {
        PredFilterCoef[e_pos][i] = STANDARD_2D_FILTER[e_pos][i];
        PredFilterCoef[f_pos][i] = STANDARD_2D_FILTER[f_pos][i];
        PredFilterCoef[j_pos][i] = STANDARD_2D_FILTER[j_pos][i];
      }
      DequantizeFilterCoefficients(DiffQFilterCoef[e_pos], FilterCoef[e_pos]);
      DequantizeFilterCoefficients(DiffQFilterCoef[f_pos], FilterCoef[f_pos]);
      DequantizeFilterCoefficients(DiffQFilterCoef[j_pos], FilterCoef[j_pos]);
      for(i = 0; i < SQR_FILTER; i++)
      {
        FilterCoef[e_pos][i] += PredFilterCoef[e_pos][i];
        FilterCoef[f_pos][i] += PredFilterCoef[f_pos][i];
        FilterCoef[j_pos][i] += PredFilterCoef[j_pos][i];
      }
      ExtendFilterCoefficientsFloat(g_pos,FilterCoef);
      ExtendFilterCoefficientsFloat(i_pos,FilterCoef);
      ExtendFilterCoefficientsFloat(k_pos,FilterCoef);
      ExtendFilterCoefficientsFloat(m_pos,FilterCoef);
      ExtendFilterCoefficientsFloat(n_pos,FilterCoef);
      ExtendFilterCoefficientsFloat(o_pos,FilterCoef);
    }
    else
    {
      multplier = 1<<FILTCOEF_BITS;
      for(j = 0; j < 15; j++)
      {
        multplier = 1<<nQBits[j];
        for(i = 0; i < SQR_FILTER; i++)
        {
          if (STANDARD_2D_FILTER[j][i]>=0) signFC = 1; else signFC = -1;
          abs_value = STANDARD_2D_FILTER[j][i]*multplier;
          abs_value = fabs(abs_value) + 0.5;
          // define a integer representation of the FC.
          PredFilterCoef16bits[j][i] = (short int)abs_value*signFC;
          FilterCoef_16bits[j][i] = PredFilterCoef16bits[j][i] + DiffQFilterCoef[j][i];
          FilterCoef[j][i]= (double)FilterCoef_16bits[j][i]/(double)multplier;
        }
      }
      ExtendFilterCoefficientsShort(h_pos,FilterCoef_16bits);
      ExtendFilterCoefficientsShort(c_pos,FilterCoef_16bits);
      ExtendFilterCoefficientsShort(d_pos,FilterCoef_16bits);
      ExtendFilterCoefficientsShort(l_pos,FilterCoef_16bits);
      ExtendFilterCoefficientsShort(g_pos,FilterCoef_16bits);
      ExtendFilterCoefficientsShort(i_pos,FilterCoef_16bits);
      ExtendFilterCoefficientsShort(k_pos,FilterCoef_16bits);
      ExtendFilterCoefficientsShort(m_pos,FilterCoef_16bits);
      ExtendFilterCoefficientsShort(n_pos,FilterCoef_16bits);
      ExtendFilterCoefficientsShort(o_pos,FilterCoef_16bits);

    }
  }


}
#endif

#ifdef E_DAIF

void CalculateFilterCoefficientsInt(void)
{
  int i, j; 
  int filterHalfLen = (FILTER_SIZE_INT/2)+1; 
  double filterIntFixed[] = 
  {
    41./2048.,  -72./2048.,  61./2048., 0., 0.,   
    -51./2048.,  -68./2048., 246./2048., 0., 0., 
    -20./2048.,  420./2048.,1229./2048., 0., 0.,
    0.,          0.,         0., 0., 0., 
    0.,          0.,         0., 0., 0., 
  };
  double predFilter[SQR_FILTER_INT]; 
  int factor;
  double is; 
  int numQBits;

  int offsetFracCodeLen[] = 
  {
    5, 4, 4, 2, 2, 2, 2, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
  };
  int bitsForOffsetFrac; 

  // predict top-left quadrant from fixed filter and quantize
  for(i = 0; i < filterHalfLen; i++)
  {
    for(j = 0; j < filterHalfLen; j++)
    {
      numQBits = numQBitsInt[i*FILTER_SIZE_INT+j]-1; // 1-bit for sign
      factor = (1<<numQBits); 
      predFilter[i*FILTER_SIZE_INT+j] = filterIntFixed[i*FILTER_SIZE_INT+j];
      is = (double)(DiffQFilterCoeffInt[i*FILTER_SIZE_INT+j])/(double)factor;
      FilterCoefInt[i*FILTER_SIZE_INT+j] = predFilter[i*FILTER_SIZE_INT+j]+is; 
    }
  }

  // predict top-right quadrant from top-left quadrant and quantize
  for(i = 0; i < filterHalfLen; i++)
  {
    for(j = filterHalfLen; j < FILTER_SIZE_INT; j++)
    {
      numQBits = numQBitsInt[i*FILTER_SIZE_INT+j]-1; // 1-bit for sign
      factor = (1<<numQBits); 
      predFilter[i*FILTER_SIZE_INT+j] = FilterCoefInt[i*FILTER_SIZE_INT+(FILTER_SIZE_INT-1-j)];
      is = (double)(DiffQFilterCoeffInt[i*FILTER_SIZE_INT+j])/(double)factor;
      FilterCoefInt[i*FILTER_SIZE_INT+j] = predFilter[i*FILTER_SIZE_INT+j]+is; 
    }
  }

  // predict bottom half from top half and quantize
  for(i = filterHalfLen; i < FILTER_SIZE_INT; i++)
  {
    for(j = 0; j < FILTER_SIZE_INT; j++)
    {
      numQBits = numQBitsInt[i*FILTER_SIZE_INT+j]-1; // 1-bit for sign
      factor = (1<<numQBits); 
      predFilter[i*FILTER_SIZE_INT+j] = FilterCoefInt[(FILTER_SIZE_INT-1-i)*FILTER_SIZE_INT+j];
      is = (double)(DiffQFilterCoeffInt[i*FILTER_SIZE_INT+j])/(double)factor;
      FilterCoefInt[i*FILTER_SIZE_INT+j] = predFilter[i*FILTER_SIZE_INT+j]+is; 
    }
  }

  // full-pel filter offset
  bitsForOffsetFrac = offsetFracCodeLen[DiffQFilterOffsetIntI];
  factor = 1<<bitsForOffsetFrac;
  FilterCoefInt[SQR_FILTER_INT-1] = (double)DiffQFilterOffsetIntI + (double)DiffQFilterOffsetIntF/(double)factor; 
  FilterCoefInt[SQR_FILTER_INT-1] = FilterCoefInt[SQR_FILTER_INT-1] *DiffQFilterOffsetIntSign; 
}

#ifdef EDAIF2
void ExtendDiagonalFilter_DAIF2(int sub_pos, double FCoef[MAX_NUM_AIF][SQR_FILTER])
{
	int j,i,mother_pel;
  int offset,act_i;
	mother_pel = FILTER_NEW_SUB_POS[FILTER_NEW_SUB_POS[sub_pos]];

  for(j=0;j<Filt_Indexes_offset[sub_pos];j++)
  {
      offset = j*POS_EQUATION_NUMBER[sub_pos];
      for(i = 0; i < POS_EQUATION_NUMBER[sub_pos]; i++)
      {
        if(offset%2==0)
          act_i = i;
        else
          act_i = POS_EQUATION_NUMBER[sub_pos]-i-1;
        FCoef[sub_pos][Calc2Filt_Indexes[sub_pos][i+offset]] = CalculatedFilter[mother_pel][act_i];
      }
  }

  FCoef[sub_pos][SQR(FILTER_SIZE)] = CalculatedFilter[mother_pel][POS_EQUATION_NUMBER[sub_pos]];  // copy filter offset
}

void SetCalc2FilterIdx()
{
  int sub_pos;
  for(sub_pos=a_pos; sub_pos <= o_pos; sub_pos++)
  {
    if(FilterFlag[sub_pos]>1)
    {// Init RAIF
      memcpy(Calc2Filt_Indexes[sub_pos],Calc2Filt_Indexes_RAIF[sub_pos],sizeof(int)*SQR_FILTER);
      memcpy(STANDARD_2D_FILTER[sub_pos],STANDARD_2D_FILTER_RAIF[sub_pos],sizeof(double)*SQR_FILTER);
    }else
    {// Init DAIF
      memcpy(Calc2Filt_Indexes[sub_pos],Calc2Filt_Indexes_DAIF[sub_pos],sizeof(int)*SQR_FILTER);
      memcpy(STANDARD_2D_FILTER[sub_pos],STANDARD_2D_FILTER_DAIF[sub_pos],sizeof(double)*SQR_FILTER);
    }
  }
}

void CalculateFilterCoefficients_EDAIF2(void)
{
  int i, j, sub_pos;

  if(img->AdaptiveFilterFlag)
  {
    for(sub_pos=a_pos; sub_pos <= o_pos; sub_pos++)
    {
      if(FilterFlag[sub_pos]>1)
      {// Init RAIF
        memcpy(Calc2Filt_Indexes[sub_pos],Calc2Filt_Indexes_RAIF[sub_pos],sizeof(int)*SQR_FILTER);
        memcpy(STANDARD_2D_FILTER[sub_pos],STANDARD_2D_FILTER_RAIF[sub_pos],sizeof(double)*SQR_FILTER);
      }else
      {// Init DAIF
        memcpy(Calc2Filt_Indexes[sub_pos],Calc2Filt_Indexes_DAIF[sub_pos],sizeof(int)*SQR_FILTER);
        memcpy(STANDARD_2D_FILTER[sub_pos],STANDARD_2D_FILTER_DAIF[sub_pos],sizeof(double)*SQR_FILTER);
      }
    }

    memset(CalculatedFilter, 0, SQR_FILTER*MAX_NUM_AIF*sizeof(double));
    {
      if (img->ImpType == IMP_FLOAT32)
      {
        for(j=0;j<MAX_NUM_AIF;j++)
        {
          NumberOfQBits = nBitsIntRepresent[FILTER_NEW_SUB_POS[FILTER_NEW_SUB_POS[j]]]+1;
          DequantizeFilterCoefficientsEDAIF2(DiffQFilterCoef[j], CalculatedFilter[j],j);
          // Dequantize filter offset
          DequantizeFilterOffsetEDAIF2(DiffQFilterOffsetI[j], DiffQFilterOffsetF[j], DiffQFilterOffsetSign[j], 
            CalculatedFilter[j]);
          CalculatedFilter[j][POS_EQUATION_NUMBER[j]] = CalculatedFilter[j][SQR_FILTER-1];
          CalculatedFilter[j][SQR_FILTER-1] = 0.0;

          ExtendDiagonalFilter_DAIF2(j, FilterCoef);
          for(i=0;i<SQR_FILTER;i++)
              FilterCoef[j][i]+=STANDARD_2D_FILTER[j][i];
        }
      }
      else
      {
        printf("EDAIF2 16-bit implementation is currently not supported!\n");
      }
      for(j=0;j<MAX_NUM_AIF;j++)
      {
        if(FilterFlag[j]%2==0)
        {
          for(i=0;i<SQR_FILTER-1;i++)
          {
            FilterCoef[j][i]=STANDARD_2D_FILTER_orig[j][i];
          }
          FilterCoef[j][SQR_FILTER-1] = 0.0;
        }
      }
    }

    if(FilterFlagInt)
      CalculateFilterCoefficientsInt(); 
    else 
    {
      for(i=0; i < SQR_FILTER_INT; i++)
        FilterCoefInt[i] = 0.;
      FilterCoefInt[12] = 1.0; 
    }
  }
}

#endif
void CalculateFilterCoefficients_EDAIF(void)
{
  int i;
  int sub_pos;
  int FILTCOEF_BITS[SQR_FILTER];

  if(img->AdaptiveFilterFlag)
  {
    for(i = 0; i < SQR(FILTER_SIZE); i++)
      FILTCOEF_BITS[i] = NumberOfQBits; 

    memset(&FilterCoef[0][0], 0, sizeof(double)*15*SQR_FILTER);
    if (img->ImpType == IMP_FLOAT32)
    {
      for(sub_pos=0; sub_pos < 15; sub_pos++)
      {
        if (!FilterFlag[sub_pos])
        {
          memcpy(&FilterCoef[sub_pos][0],&STANDARD_2D_FILTER_orig[sub_pos][0],sizeof(double)*SQR_FILTER);
        }
        else if((SymmetryPosition[sub_pos]))
        {
          DequantizeFilterCoefficientsExt(DiffQFilterCoef[sub_pos], FILTCOEF_BITS, 
            DiffQFilterOffsetI[sub_pos], DiffQFilterOffsetF[sub_pos], DiffQFilterOffsetSign[sub_pos], 
            FilterCoef[sub_pos]);
          for(i=0; i < SQR_FILTER; i++)
          {
            FilterCoef[sub_pos][i] += STANDARD_2D_FILTER[sub_pos][i];
          }
        }
      }
      if(FilterFlag[c_pos])
        ExtendFilterCoefficientsFloat(c_pos,FilterCoef);
      if(FilterFlag[d_pos])
        ExtendFilterCoefficientsFloat(d_pos,FilterCoef);
      if(FilterFlag[g_pos])
        ExtendFilterCoefficientsFloat(g_pos,FilterCoef);
      if(FilterFlag[h_pos])
        ExtendFilterCoefficientsFloat(h_pos,FilterCoef);
      if(FilterFlag[i_pos])
        ExtendFilterCoefficientsFloat(i_pos,FilterCoef);
      if(FilterFlag[k_pos])
        ExtendFilterCoefficientsFloat(k_pos,FilterCoef);
      if(FilterFlag[l_pos])
        ExtendFilterCoefficientsFloat(l_pos,FilterCoef);
      if(FilterFlag[m_pos])
        ExtendFilterCoefficientsFloat(m_pos,FilterCoef);
      if(FilterFlag[n_pos])
        ExtendFilterCoefficientsFloat(n_pos,FilterCoef);
      if(FilterFlag[o_pos])
        ExtendFilterCoefficientsFloat(o_pos,FilterCoef);
      // offset 
      for(sub_pos=0; sub_pos < 15; sub_pos++)
      {
        if(!(SymmetryPosition[sub_pos]))
          FilterCoef[sub_pos][SQR_FILTER-1] = FilterCoef[FILTER_NEW_SUB_POS[sub_pos]][SQR_FILTER-1];
      }
    }
    else
    {
      printf("ERROR: Current implementation of DAIF supports only IMP_FLOAT32\n");
    }

    if(FilterFlagInt)
      CalculateFilterCoefficientsInt(); 
    else 
    {
      for(i=0; i < SQR_FILTER_INT; i++)
        FilterCoefInt[i] = 0.;
      FilterCoefInt[12] = 1.0; 
    }
  }
}
#endif  // E_DAIF

#ifdef EAIF
// flip the coefficients horizontally 
static void flipHoriz(double *coeff, double *coeffT)
{
  int i, j;
  for(i = 0; i < FILTER_SIZE; i++)
    for(j = 0; j < FILTER_SIZE; j++)
      coeffT[i*FILTER_SIZE+j] = coeff[i*FILTER_SIZE+(FILTER_SIZE-1-j)];
}

// flip the coefficients vertically
static void flipVert(double *coeff, double *coeffT)
{
  int i, j;
  for(i = 0; i < FILTER_SIZE; i++)
    for(j = 0; j < FILTER_SIZE; j++)
      coeffT[i*FILTER_SIZE+j] = coeff[(FILTER_SIZE-1-i)*FILTER_SIZE+j];
}

// Propogate the filter coefficients from SymmetryPositions to non-symmetry positions. 
// For example, from a_pos to c_pos, from d_pos to l_pos, etc. 
static void propogateFilterCoef(int SUB_POS)
{
  if(!SymmetryPosition[SUB_POS])
    return; 

  if(SUB_POS == a_pos)
  {
    flipHoriz(FilterCoef[a_pos], FilterCoef[c_pos]);
    FilterCoef[c_pos][SQR_FILTER-1] = FilterCoef[a_pos][SQR_FILTER-1];
    FilterFlag[c_pos] = FilterFlag[a_pos]; 
  }
  else if(SUB_POS == d_pos)
  {
    flipVert (FilterCoef[d_pos], FilterCoef[l_pos]);
    FilterCoef[l_pos][SQR_FILTER-1] = FilterCoef[d_pos][SQR_FILTER-1];
    FilterFlag[l_pos] = FilterFlag[d_pos]; 
  }
  else if(SUB_POS == e_pos)
  {
    flipHoriz(FilterCoef[e_pos], FilterCoef[g_pos]);
    flipVert (FilterCoef[e_pos], FilterCoef[m_pos]);
    flipHoriz(FilterCoef[m_pos], FilterCoef[o_pos]);
    FilterCoef[g_pos][SQR_FILTER-1] = FilterCoef[e_pos][SQR_FILTER-1];
    FilterCoef[m_pos][SQR_FILTER-1] = FilterCoef[e_pos][SQR_FILTER-1];
    FilterCoef[o_pos][SQR_FILTER-1] = FilterCoef[e_pos][SQR_FILTER-1];
    FilterFlag[g_pos] = FilterFlag[e_pos]; 
    FilterFlag[m_pos] = FilterFlag[e_pos]; 
    FilterFlag[o_pos] = FilterFlag[e_pos]; 
  }
  else if(SUB_POS == f_pos)
  {
    flipVert (FilterCoef[f_pos], FilterCoef[n_pos]);
    FilterCoef[n_pos][SQR_FILTER-1] = FilterCoef[f_pos][SQR_FILTER-1];
    FilterFlag[c_pos] = FilterFlag[a_pos]; 
  }
  else if(SUB_POS == i_pos)
  {
    flipHoriz(FilterCoef[i_pos], FilterCoef[k_pos]);
    FilterCoef[k_pos][SQR_FILTER-1] = FilterCoef[i_pos][SQR_FILTER-1];
    FilterFlag[k_pos] = FilterFlag[i_pos]; 
  }

  return; 
}

/*!
************************************************************************
* \brief
*    Calculate the EAIF filter coefficients from the received 
*    DiffQFilterCoef and DiffQFilterOffset. Populate the entire 
*    FilterCoef array
************************************************************************
*/
void CalculateFilterCoefficients_EAIF(void)
{
  int i, sub_pos;
  int *numQBits;

  // at first predict b and h position using inter prediction, b and h are already quantized in last frame
  if(img->AdaptiveFilterFlag)
  {
    for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
    {
      if(FilterFlag[sub_pos] && SymmetryPosition[sub_pos])
      {
        if((sub_pos+1)/4 == 0)  // 1D-filter horizontal
        {
          numQBits = numQBits1DH; 
          for(i = 0; i < SQR_FILTER; i++)
            PredFilterCoef[sub_pos][i] = STANDARD_2D_FILTER_orig[sub_pos][i]; 
        }
        else if ((sub_pos+1)%4 == 0)  // 1D-filter vertical 
        {
          numQBits = numQBits1DV; 
          for(i = 0; i < SQR_FILTER; i++)
            PredFilterCoef[sub_pos][i] = STANDARD_2D_FILTER_orig[sub_pos][i]; 
        }
        else
        {
          numQBits = numQBits2D; 
          for(i = 0; i < SQR_FILTER; i++)
            PredFilterCoef[sub_pos][i] = 0.; 
        }
        DequantizeFilterCoefficientsExt(DiffQFilterCoef[sub_pos], numQBits, 
          DiffQFilterOffsetI[sub_pos], DiffQFilterOffsetF[sub_pos], DiffQFilterOffsetSign[sub_pos], 
          FilterCoef[sub_pos]);
        // add differences to predicted values
        for(i=0; i < SQR_FILTER; i++)
        {
          FilterCoef[sub_pos][i] += PredFilterCoef[sub_pos][i];
        }
      }
      else if(!FilterFlag[sub_pos])
      {
        for(i=0; i < SQR_FILTER; i++)
          FilterCoef[sub_pos][i] = STANDARD_2D_FILTER[sub_pos][i];
        FilterCoef[sub_pos][SQR_FILTER-1] = 0.; 
      }
    }
    propogateFilterCoef(a_pos);
    propogateFilterCoef(d_pos);
    propogateFilterCoef(e_pos);
    propogateFilterCoef(f_pos);
    propogateFilterCoef(i_pos);

    // full-pel filter
    if(FilterFlagInt)
      CalculateFilterCoefficientsInt(); 
    else 
    {
      for(i=0; i < SQR_FILTER_INT; i++)
        FilterCoefInt[i] = 0.;
      FilterCoefInt[12] = 1.0; 
    }
  }
}
#endif

// separable aif (BEGIN)
/*!
************************************************************************
* \brief
*    Predict and calculate filter coefficients for separable filter
************************************************************************
*/
void CalculateFilterCoefficientsSep(void)
{
  int i,sub_pel;
  if(img->AdaptiveFilterFlag)
  {
    for(sub_pel = a_pos; sub_pel <= o_pos; sub_pel++)
    {
      if(FilterFlag[sub_pel])
      {
        DequantizeFilterCoefficients(DiffQFilterCoef[sub_pel], FilterCoef[sub_pel]);
        // add differences to predicted values
        for(i=0; i < FILTER_SIZE; i++)
        {
          FilterCoef[sub_pel][i] += STANDARD_2D_FILTER_SEP[sub_pel][i];
        }
      }
      else
      {
        for(i=0; i < FILTER_SIZE; i++)
          FilterCoef[sub_pel][i] = STANDARD_2D_FILTER_SEP[sub_pel][i];
      }
    }
  }
  else
    for(sub_pel = a_pos; sub_pel <= o_pos; sub_pel++)
      for(i = 0; i < FILTER_SIZE; i++)
        FilterCoef[sub_pel][i] = STANDARD_2D_FILTER_SEP[sub_pel][i];
}
// separable aif (END)
/*!
************************************************************************
* \brief
*   inverse quantization of filter coefficients
*    input:
*        int QFCoef[][] - filter coefficients to be dequantized
*    output:
*        int FCoef[][]   - double filter coefficients
************************************************************************
*/
void DequantizeFilterCoefficients(int QFCoef[SQR_FILTER], double FCoef[SQR_FILTER])
{
  double number_of_steps = pow(2,NumberOfQBits-1);  // 1 bit for sign, NumberOfQBits-1 for amplitude
  int j;
  double is = 0.0;
  for(j=0; j < SQR_FILTER; j++)
  {
    is = (double)(QFCoef[j]);
    is = is/(double)number_of_steps;
    FCoef[j]  = is;  // reconstructed (quantized) filter coefficients
  }
}

#ifdef EDAIF2
void DequantizeFilterCoefficientsEDAIF2(int QFCoef[SQR_FILTER], double FCoef[SQR_FILTER], int sub_pos)
{
  double factor;
  double number_of_steps;
  int j,q_j;
  double is = 0.0;
  for(j=0; j < POS_EQUATION_NUMBER[sub_pos]; j++)
  {
    {
      q_j = Calc2Filt_Indexes[sub_pos][j];

      if((sub_pos == a_pos)||(sub_pos == b_pos)||(sub_pos == c_pos))
      {
        number_of_steps = pow(2,numQBits1DH[q_j]-1);  // 1 bit for sign, NumberOfQBits-1 for amplitude
      }
      else if((sub_pos == d_pos)||(sub_pos == h_pos)||(sub_pos == l_pos))
      {
        number_of_steps = pow(2,numQBits1DV[q_j]-1);  // 1 bit for sign, NumberOfQBits-1 for amplitude
      }
      else if((sub_pos == e_pos)||(sub_pos == g_pos)||(sub_pos == m_pos)||(sub_pos == o_pos))
      {
        number_of_steps = pow(2,numQBitsDIF[q_j]-1);  // 1 bit for sign, NumberOfQBits-1 for amplitude
      }
      else if(FilterFlag[sub_pos]<=1)// DAIF
      {
        number_of_steps = pow(2,numQBitsDIF[q_j]-1);  // 1 bit for sign, NumberOfQBits-1 for amplitude
      }
      else 
      {
        number_of_steps = pow(2,numQBits2D[q_j]-1);  // 1 bit for sign, NumberOfQBits-1 for amplitude
      }
    }
    factor = number_of_steps;
    if(QFCoef[j])
    {
      is = (double)(QFCoef[j])/factor;
    }
    else 
      is = 0.;
    FCoef[j]  = is;		// reconstructed (quantized) filter coefficients
  }
}
#endif

#ifdef E_DAIF

// dequantize coefficients and offset
void DequantizeFilterCoefficientsExt(int QFCoef[SQR_FILTER], int *numQBits, int QFOffsetI, int QFOffsetF, int QFOffsetSign, double FCoef[SQR_FILTER])
{
  double factor; 
  int j;
  double is = 0.0;
  int offsetFracCodeLen[] = 
  {
    5, 4, 4, 2, 2, 2, 2, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
  };
  int bitsForOffsetFrac; 

  for(j=0; j < SQR_FILTER-1; j++)
  {
    if(QFCoef[j])
    {
      factor = pow(2,numQBits[j]-1);
      is = (double)(QFCoef[j])/factor;
    }
    else 
      is = 0.;
    FCoef[j]  = is;		// reconstructed (quantized) filter coefficients
  }

  // offset
  bitsForOffsetFrac = offsetFracCodeLen[QFOffsetI];
  factor = pow(2.,bitsForOffsetFrac);
  FCoef[j] = (double)QFOffsetI + (double)QFOffsetF/factor; 
  FCoef[j] = FCoef[j]*QFOffsetSign; 
}

#endif  // E_DAIF

#ifdef EDAIF2
void DequantizeFilterOffsetEDAIF2(int QFOffsetI, int QFOffsetF, int QFOffsetSign, double FCoef[SQR_FILTER])
{
  double factor; 
  int j;
  int offsetFracCodeLen[] = 
  {
    5, 4, 4, 2, 2, 2, 2, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
  };
  int bitsForOffsetFrac; 

  j = SQR_FILTER-1;

  // offset
  bitsForOffsetFrac = offsetFracCodeLen[QFOffsetI];
  factor = pow(2.,bitsForOffsetFrac);
  FCoef[j] = (double)QFOffsetI + (double)QFOffsetF/factor; 
  FCoef[j] = FCoef[j]*QFOffsetSign; 
}

#endif

/*!
************************************************************************
* \brief
*    quantize filter coefficients
*    input:
*        double FCoef[][] - filter coefficients to be quantized
*    output:
*        int QFCoef[][]   - quantized filter coefficients, 
*                             represented in step numbers
************************************************************************
*/
void QuantizeFilterCoefficients(double FCoef[SQR_FILTER], int QFCoef[SQR_FILTER])
{
  int number_of_steps = 2^(NumberOfQBits-1);  // 1 bit for sign, NumberOfQBits-1 for amplitude
  int sign;
  int j;
  for(j=0; j < SQR_FILTER; j++)
  {
    sign = FCoef[j] >= 0.0 ? 1:-1;
    QFCoef[j] = sign*(int)(fabs(FCoef[j])*number_of_steps + 0.5);  // step number
    FCoef[j] = (double)(QFCoef[j])/(double)number_of_steps; // reconstructed (quantized) filter coefficients
  }
}


/*!
************************************************************************
* \brief
*    extend the filter to 2D (FILTER_SIZE X FILTER_SIZE)
*
*    input:
*      sub_pos - sub-position
*    input output:
*      FCoef[][] - Filter to be extended at symmetric positions,
*                  using FILTER_NEW_SUB_POS
************************************************************************
*/
void ExtendFilterCoefficientsFloat(int sub_pos, double FCoef[15][SQR_FILTER])
{
  int i;
  int filter_type = img->AdaptiveFilterFlag;
#ifndef DIRECTIONAL_FILTER
  for(i = 0; i < SQR_FILTER; i++)
    FCoef[sub_pos][i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i]];
#else
  if(filter_type == FILTER_TYPE_2D_NS)
  {
    for(i = 0; i < SQR_FILTER; i++)
      FCoef[sub_pos][i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i]];
  }

#ifdef E_DAIF
  else if (filter_type == FILTER_TYPE_1D || filter_type == FILTER_TYPE_EDAIF)
#else
  else if (filter_type == FILTER_TYPE_1D)
#endif
  {
    int j;
    if (IsDiagonal1D[sub_pos]==1)
    {
      for(i = 0; i < FILTER_SIZE; i++)
      {
        j = FILTER_SIZE-1-i;
        FCoef[sub_pos][i*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i*FILTER_SIZE + i]];
      }
    }
    else if (IsDiagonal1D[sub_pos]==2)
    {
      for(i = 0; i < FILTER_SIZE; i++)
      {
        j = FILTER_SIZE-1-i;
        //FCoef[sub_pos][j*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][j*FILTER_SIZE + i]];
        FCoef[sub_pos][i*FILTER_SIZE + j] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i*FILTER_SIZE + j]];
      }
    }
    else if(IsDiagonal1D[sub_pos]==3)
    {
      if(sub_pos == j_pos)
      {
        for(i = 0; i < FILTER_SIZE; i++)
        {
          j = FILTER_SIZE-1-i;
          FCoef[sub_pos][j*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][j*FILTER_SIZE + i]];
          FCoef[sub_pos][i*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i*FILTER_SIZE + i]];
        }
      }
      else 
      {
        for(i = 0; i < FILTER_SIZE; i++)
        {
          j = i;
          FCoef[sub_pos][j*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][j*FILTER_SIZE + i]];
          j = FILTER_SIZE-1-i;
          FCoef[sub_pos][i*FILTER_SIZE + j] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i*FILTER_SIZE + j]];
        }
      }
    }
    else 
      for(i = 0; i < SQR_FILTER; i++)
        FCoef[sub_pos][i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i]];

  }
#ifdef EAIF
  else if(filter_type == FILTER_TYPE_EAIF)
  {
		for(i = 0; i < SQR_FILTER-1; i++)
			FCoef[sub_pos][i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i]];
    FCoef[sub_pos][i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][i];
  }
#endif

  return;
#endif
}


/*!
************************************************************************
* \brief
*    extend the filter to 2D (FILTER_SIZE X FILTER_SIZE)
*
*    input:
*      sub_pos - sub-position
*     input output:
*       FCoef[][] - Filter to be extended at symmetric positions,
*                     using FILTER_NEW_SUB_POS
************************************************************************
*/
void ExtendFilterCoefficientsInt(int sub_pos, int FCoef[15][SQR_FILTER])
{
  int i;
  int filter_type = img->AdaptiveFilterFlag;
#ifndef DIRECTIONAL_FILTER
  for(i = 0; i < SQR_FILTER; i++)
    FCoef[sub_pos][i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i]];
#else
  if(filter_type == FILTER_TYPE_2D_NS)
  {
    for(i = 0; i < SQR_FILTER; i++)
      FCoef[sub_pos][i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i]];
  }
#ifdef E_DAIF
  else if (filter_type == FILTER_TYPE_1D || filter_type == FILTER_TYPE_EDAIF )
#else
  else if (filter_type == FILTER_TYPE_1D)
#endif
  {
    int j;
    if (IsDiagonal1D[sub_pos]==1)
    {
      for(i = 0; i < FILTER_SIZE; i++)
      {
        j = FILTER_SIZE-1-i;
        FCoef[sub_pos][i*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i*FILTER_SIZE + i]];
      }
    }
    //else if ((sub_pos == g_pos)||(sub_pos == m_pos))
    else if (IsDiagonal1D[sub_pos]==2)
    {
      for(i = 0; i < FILTER_SIZE; i++)
      {
        j = FILTER_SIZE-1-i;
        //FCoef[sub_pos][j*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][j*FILTER_SIZE + i]];
        FCoef[sub_pos][i*FILTER_SIZE + j] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i*FILTER_SIZE + j]];
      }
    }
    else if(IsDiagonal1D[sub_pos]==3)
    {
      if(sub_pos == j_pos)
      {
        for(i = 0; i < FILTER_SIZE; i++)
        {
          j = FILTER_SIZE-1-i;
          FCoef[sub_pos][j*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][j*FILTER_SIZE + i]];
          FCoef[sub_pos][i*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i*FILTER_SIZE + i]];
        }
      }
      else 
      {
        for(i = 0; i < FILTER_SIZE; i++)
        {
          j = i;
          FCoef[sub_pos][j*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][j*FILTER_SIZE + i]];
          j = FILTER_SIZE-1-i;
          FCoef[sub_pos][i*FILTER_SIZE + j] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i*FILTER_SIZE + j]];
        }
      }
    }
    else 
      for(i = 0; i < SQR_FILTER; i++)
        FCoef[sub_pos][i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i]];
  }
  return;
#endif
}
/*!
************************************************************************
* \brief
*    extend the filter DAIF with respect to symmeties.
*
*    input:
*      sub_pos - sub-position
*     input output:
*       FCoef[][] - Filter to be extended at symmetric positions,
*                     using FILTER_NEW_SUB_POS
************************************************************************
*/
void ExtendFilterCoefficientsShort(int sub_pos, short FCoef[15][SQR_FILTER])
{
  int i;
  int filter_type = img->AdaptiveFilterFlag;
#ifndef DIRECTIONAL_FILTER
  for(i = 0; i < SQR_FILTER; i++)
    FCoef[sub_pos][i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i]];
#else
  if(filter_type == FILTER_TYPE_2D_NS)
  {
    for(i = 0; i < SQR_FILTER; i++)
      FCoef[sub_pos][i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i]];
  }
  else if (filter_type == FILTER_TYPE_1D)
  {
    int j;
    if (IsDiagonal1D[sub_pos]==1)
    {
      for(i = 0; i < FILTER_SIZE; i++)
      {
        j = FILTER_SIZE-1-i;
        FCoef[sub_pos][i*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i*FILTER_SIZE + i]];
      }
    }
    else if (IsDiagonal1D[sub_pos]==2)
    {
      for(i = 0; i < FILTER_SIZE; i++)
      {
        j = FILTER_SIZE-1-i;
        FCoef[sub_pos][i*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i*FILTER_SIZE + i]];
      }
    }
    else if(IsDiagonal1D[sub_pos]==3)
    {
      if(sub_pos == j_pos)
      {
        for(i = 0; i < FILTER_SIZE; i++)
        {
          j = FILTER_SIZE-1-i;
          FCoef[sub_pos][j*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][j*FILTER_SIZE + i]];
          FCoef[sub_pos][i*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i*FILTER_SIZE + i]];
        }
      }
      else 
      {
        for(i = 0; i < FILTER_SIZE; i++)
        {
          j = i;
          FCoef[sub_pos][j*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][j*FILTER_SIZE + i]];
          j = FILTER_SIZE-1-i;
          FCoef[sub_pos][i*FILTER_SIZE + j] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i*FILTER_SIZE + j]];
        }
      }
    }
    else 
      for(i = 0; i < SQR_FILTER; i++)
        FCoef[sub_pos][i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i]];
  }
  return;
#endif
}

#ifdef EIGHTH_PEL
void GetBlockWith2DAIF_eighth_pel( int             ref_idx,
                                  StorablePicture **list,
                                  int             x_pos,
                                  int             y_pos,
                                  int             mvx,
                                  int             mvy,
                                  int             img_width,
                                  int             img_height,
                                  int             block[BLOCK_SIZE][BLOCK_SIZE])
{
  int block2[BLOCK_SIZE][BLOCK_SIZE]; 
  int block3[BLOCK_SIZE][BLOCK_SIZE]; 
  int block4[BLOCK_SIZE][BLOCK_SIZE]; 
  int dmvx = (mvx%8+8)%8;
  int dmvy = (mvy%8+8)%8;


  if(!(dmvx&1) && !(dmvy&1))                                                              //              1/4-pel
    GetBlockWith2DAIF_quarter_pel(ref_idx, list, x_pos, y_pos, mvx/2, mvy/2, img_width,img_height, block);
  else                                                                                                                            //              1/$
  {
    if(dmvx&1 && !(dmvy&1))                                                         //              horizontal 1/8 and vertical 1/4
    {
      GetBlockWith2DAIF_quarter_pel(ref_idx, list, x_pos, y_pos, (mvx+1)/2, mvy/2, img_width,img_height, block);
      GetBlockWith2DAIF_quarter_pel(ref_idx, list, x_pos, y_pos, (mvx-1)/2, mvy/2, img_width,img_height, block2);
      average_block(block, block2);
    }
    else if(!(dmvx&1) && dmvy&1)                                    //              horizontal 1/4 and vertical 1/8
    {
      GetBlockWith2DAIF_quarter_pel(ref_idx, list, x_pos, y_pos, mvx/2, (mvy+1)/2, img_width,img_height, block);
      GetBlockWith2DAIF_quarter_pel(ref_idx, list, x_pos, y_pos, mvx/2, (mvy-1)/2, img_width,img_height, block2);
      average_block(block, block2);
    }
    else if((dmvx+dmvy)%4)// diagonal positions upper left - down right interpolation
    {
      GetBlockWith2DAIF_quarter_pel(ref_idx, list, x_pos, y_pos, (mvx-1)/2, (mvy-1)/2, img_width,img_height, block);
      GetBlockWith2DAIF_quarter_pel(ref_idx, list, x_pos, y_pos, (mvx+1)/2, (mvy+1)/2, img_width,img_height, block2);

      GetBlockWith2DAIF_quarter_pel(ref_idx, list, x_pos, y_pos, (mvx+1)/2, (mvy-1)/2, img_width,img_height, block3);
      GetBlockWith2DAIF_quarter_pel(ref_idx, list, x_pos, y_pos, (mvx-1)/2, (mvy+1)/2, img_width,img_height, block4);
      average_block4(block, block2, block3, block4);
    }
    else // diagonal positions upper right - down left interpolation
    {
      GetBlockWith2DAIF_quarter_pel(ref_idx, list, x_pos, y_pos, (mvx+1)/2, (mvy-1)/2, img_width,img_height, block);
      GetBlockWith2DAIF_quarter_pel(ref_idx, list, x_pos, y_pos, (mvx-1)/2, (mvy+1)/2, img_width,img_height, block2);

      GetBlockWith2DAIF_quarter_pel(ref_idx, list, x_pos, y_pos, (mvx+1)/2, (mvy+1)/2, img_width,img_height, block3);
      GetBlockWith2DAIF_quarter_pel(ref_idx, list, x_pos, y_pos, (mvx-1)/2, (mvy-1)/2, img_width,img_height, block4);
      average_block4(block, block2, block3, block4);
    }
    //    average_block(block, block2);
  }
}

// separable aif (BEGIN)
void GetBlockWithSepAIF_eighth_pel(int             ref_idx, 
                                   StorablePicture **list,
                                   int             x_pos,
                                   int             y_pos,
                                   int             xoff4x4,
                                   int             yoff4x4,
                                   int             mvx,
                                   int             mvy,
                                   int             img_width,
                                   int             img_height,
                                   int             block_size_x,
                                   int             block_size_y,
                                   int             block[MB_BLOCK_SIZE][MB_BLOCK_SIZE])
{
  int block2[MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  int block3[MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  int block4[MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  int dmvx = (mvx%8+8)%8;
  int dmvy = (mvy%8+8)%8;


  if(!(dmvx&1) && !(dmvy&1))                                                              //              1/4-pel
    GetBlockWithSepAIF_quarter_pel(ref_idx, list, x_pos, y_pos, xoff4x4, yoff4x4, mvx/2, mvy/2, img_width,img_height, block_size_x, block_size_y, block); 
  else                                                                                                                            //              1/$
  {
    if(dmvx&1 && !(dmvy&1))                                                         //              horizontal 1/8 and vertical 1/4
    {
      GetBlockWithSepAIF_quarter_pel(ref_idx, list, x_pos, y_pos, xoff4x4, yoff4x4, (mvx+1)/2, mvy/2, img_width,img_height, block_size_x, block_size_y, block); 
      GetBlockWithSepAIF_quarter_pel(ref_idx, list, x_pos, y_pos, xoff4x4, yoff4x4, (mvx-1)/2, mvy/2, img_width,img_height, block_size_x, block_size_y, block2); 
      average_block_2(block, block2);
    }
    else if(!(dmvx&1) && dmvy&1)                                    //              horizontal 1/4 and vertical 1/8
    {
      GetBlockWithSepAIF_quarter_pel(ref_idx, list, x_pos, y_pos, xoff4x4, yoff4x4, mvx/2, (mvy+1)/2, img_width,img_height, block_size_x, block_size_y, block);
      GetBlockWithSepAIF_quarter_pel(ref_idx, list, x_pos, y_pos, xoff4x4, yoff4x4, mvx/2, (mvy-1)/2, img_width,img_height, block_size_x, block_size_y, block2);
      average_block_2(block, block2);
    }
    else if((dmvx+dmvy)%4)// diagonal positions upper left - down right interpolation
    {
      GetBlockWithSepAIF_quarter_pel(ref_idx, list, x_pos, y_pos, xoff4x4, yoff4x4, (mvx-1)/2, (mvy-1)/2, img_width,img_height, block_size_x, block_size_y, block);
      GetBlockWithSepAIF_quarter_pel(ref_idx, list, x_pos, y_pos, xoff4x4, yoff4x4, (mvx+1)/2, (mvy+1)/2, img_width,img_height, block_size_x, block_size_y, block2);

      GetBlockWithSepAIF_quarter_pel(ref_idx, list, x_pos, y_pos, xoff4x4, yoff4x4, (mvx+1)/2, (mvy-1)/2, img_width,img_height, block_size_x, block_size_y, block3);
      GetBlockWithSepAIF_quarter_pel(ref_idx, list, x_pos, y_pos, xoff4x4, yoff4x4, (mvx-1)/2, (mvy+1)/2, img_width,img_height, block_size_x, block_size_y, block4);
      average_block4_2(block, block2, block3, block4);
    }
    else // diagonal positions upper right - down left interpolation
    {
      GetBlockWithSepAIF_quarter_pel(ref_idx, list, x_pos, y_pos, xoff4x4, yoff4x4, (mvx+1)/2, (mvy-1)/2, img_width,img_height, block_size_x, block_size_y, block);
      GetBlockWithSepAIF_quarter_pel(ref_idx, list, x_pos, y_pos, xoff4x4, yoff4x4, (mvx-1)/2, (mvy+1)/2, img_width,img_height, block_size_x, block_size_y, block2);

      GetBlockWithSepAIF_quarter_pel(ref_idx, list, x_pos, y_pos, xoff4x4, yoff4x4, (mvx+1)/2, (mvy+1)/2, img_width,img_height, block_size_x, block_size_y, block3);
      GetBlockWithSepAIF_quarter_pel(ref_idx, list, x_pos, y_pos, xoff4x4, yoff4x4, (mvx-1)/2, (mvy-1)/2, img_width,img_height, block_size_x, block_size_y, block4);
      average_block4_2(block, block2, block3, block4);
    }
    //    average_block(block, block2);
  }
}
// separable aif (END)

void GetBlockWith2DAIF(int             ref_idx,
                       StorablePicture **list,
                       int             x_pos,
                       int             y_pos,
                       int             mvx,
                       int             mvy,
                       int             img_width,
                       int             img_height,
                       int             block[BLOCK_SIZE][BLOCK_SIZE])
{
  if(img->mv_res)
    GetBlockWith2DAIF_eighth_pel (ref_idx, list, x_pos, y_pos, mvx, mvy, img_width,img_height, block);
  else
    GetBlockWith2DAIF_quarter_pel(ref_idx, list, x_pos, y_pos, mvx, mvy, img_width,img_height, block);
}

#ifdef EAIF
/*!
************************************************************************
* \brief
*    Interpolation of 1/4 subpixel with EAIF filter
************************************************************************
*/
void GetBlockWithEAIF(int             ref_idx,
                       StorablePicture **list,
                       int             x_pos,
                       int             y_pos,
                       int             mvx,
                       int             mvy,
											 int             img_width,
                       int             img_height,
                       int             block[BLOCK_SIZE][BLOCK_SIZE])
{
  int i,j, fi, fj;
	int fpos_x, fpos_y;
	double result;
	int is;
  int sub_pos, mvx_sub, mvy_sub;
  int is1D, isFP, is1DH; 
  imgpel sum; 
  int start_x, start_y;
  int end_x, end_y;

  // do not handle 1/8-pixel yet 
  if(img->mv_res)
    return; 
  
  if(mvx >= 0)
    mvx_sub = mvx%4;      						// x-sub-coordinate in a 4x4block
  else
    mvx_sub = (4-abs(mvx)%4)%4;
  if(mvy >= 0)
    mvy_sub = mvy%4;      						// y-sub-coordinate in a 4x4block
  else
    mvy_sub = (4-abs(mvy)%4)%4;
  sub_pos = mvx_sub + 4*mvy_sub-1;    // pos 0..14 in a 4x4 block

  if(sub_pos != -1)
  {
    isFP = 0;
    is1D = ((sub_pos+1)/4==0) || ((sub_pos+1)%4 == 0);
    if(is1D)
      is1DH = (sub_pos+1)/4==0;
    else is1DH = 0;
  }
  else
  {
    isFP = 1; is1D = is1DH = 0; 
  }

  if(sub_pos == a_pos || sub_pos == c_pos)   
  {
  	for(i = 0; i < BLOCK_SIZE; i++)
	  for(j = 0; j < BLOCK_SIZE; j++)
		{
				result = 0.0;
    		fj = FILTER_OFFSET;
        fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
    		for(fi = 0; fi < FILTER_SIZE; fi++)
        {
    			fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
    			result += FilterCoef[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x];
    		}
        result += FilterCoef[sub_pos][SQR_FILTER-1]; 
   			is = (int)(0.5+result);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        block[i][j] = max(0, min(((1<<img->bitdepth_luma)-1),is));
#else
        block[i][j] = max(0, min(255,is));
#endif
    }
  }
  else if (sub_pos == b_pos)
  {
  	for(i = 0; i < BLOCK_SIZE; i++)
	  for(j = 0; j < BLOCK_SIZE; j++)
		{
				result = 0.0;
    		fj = FILTER_OFFSET;
    	  fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
    		for(fi = 0; fi < FILTER_SIZE/2; fi++)
        {
    			fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
          sum = list[ref_idx]->imgY[fpos_y][fpos_x];
    			fpos_x = FindPosition(img_width,  i+x_pos, FILTER_SIZE-1-fi,mvx);
          sum += list[ref_idx]->imgY[fpos_y][fpos_x];
    			result += FilterCoef[sub_pos][FILTER_SIZE*fj+fi]*sum;
    		}
        result += FilterCoef[sub_pos][SQR_FILTER-1]; 
   			is = (int)(0.5+result);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        block[i][j] = max(0, min(((1<<img->bitdepth_luma)-1),is));
#else
        block[i][j] = max(0, min(255,is));
#endif
    }
  }
  else if(sub_pos == d_pos || sub_pos == l_pos) 
  {
  	for(i = 0; i < BLOCK_SIZE; i++)
	  for(j = 0; j < BLOCK_SIZE; j++)
		{
				result = 0.0;
    		fi = FILTER_OFFSET;
        fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
    		for(fj = 0; fj < FILTER_SIZE; fj++)
        {
    			fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
    			result += FilterCoef[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x];
    		}
        result += FilterCoef[sub_pos][SQR_FILTER-1]; 
   			is = (int)(0.5+result);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        block[i][j] = max(0, min(((1<<img->bitdepth_luma)-1),is));
#else
        block[i][j] = max(0, min(255,is));
#endif
    }
  }
  else if(sub_pos == h_pos)
  {
  	for(i = 0; i < BLOCK_SIZE; i++)
	  for(j = 0; j < BLOCK_SIZE; j++)
		{
				result = 0.0;
    		fi = FILTER_OFFSET;
        fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
    		for(fj = 0; fj < FILTER_SIZE/2; fj++)
        {
    			fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
          sum = list[ref_idx]->imgY[fpos_y][fpos_x];
    			fpos_y = FindPosition(img_height, j+y_pos, FILTER_SIZE-1-fj,mvy);
          sum += list[ref_idx]->imgY[fpos_y][fpos_x];
    			result += FilterCoef[sub_pos][FILTER_SIZE*fj+fi]*sum;
    		}
        result += FilterCoef[sub_pos][SQR_FILTER-1]; 
   			is = (int)(0.5+result);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        block[i][j] = max(0, min(((1<<img->bitdepth_luma)-1),is));
#else
        block[i][j] = max(0, min(255,is));
#endif
    }
  }
  else if(sub_pos == e_pos || sub_pos == g_pos || sub_pos == m_pos || sub_pos == o_pos) // 2-D sub-pixel position 
  {
    if(FilterFlag[sub_pos])
    {
      start_x = 1; end_x = FILTER_SIZE-1;
      start_y = 1; end_y = FILTER_SIZE-1;
    }
    else
    {
      start_x = 0; end_x = FILTER_SIZE;
      start_y = 0; end_y = FILTER_SIZE;
    }
	  for(i = 0; i < BLOCK_SIZE; i++)
    for(j = 0; j < BLOCK_SIZE; j++)
	  {
			  result = 0.0;
  		  for(fi = start_x; fi < end_x; fi++)
  		  for(fj = start_y; fj < end_y; fj++)
  		  {
  			  fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
  			  fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
  			  result += FilterCoef[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x];
  		  }
        result += FilterCoef[sub_pos][SQR_FILTER-1]; 
 			  is = (int)(0.5+result);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        block[i][j] = max(0, min(((1<<img->bitdepth_luma)-1),is));
#else
        block[i][j] = max(0, min(255,is));
#endif
    }
  }
  else if(sub_pos == f_pos || sub_pos == n_pos)
  {
    if(FilterFlag[sub_pos])
    {
      start_x = 1; end_x = FILTER_SIZE/2;
      start_y = 1; end_y = FILTER_SIZE-1;
    }
    else
    {
      start_x = 0; end_x = FILTER_SIZE/2;
      start_y = 0; end_y = FILTER_SIZE;
    }
  	for(i = 0; i < BLOCK_SIZE; i++)
	  for(j = 0; j < BLOCK_SIZE; j++)
		{
				result = 0.0;
    		for(fi = start_x; fi < end_x; fi++)
    		for(fj = start_y; fj < end_y; fj++)
    		{
    			fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
    			fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
          sum = list[ref_idx]->imgY[fpos_y][fpos_x];
    			fpos_x = FindPosition(img_width,  i+x_pos, FILTER_SIZE-1-fi,mvx);
    			//fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
          sum += list[ref_idx]->imgY[fpos_y][fpos_x];
    			result += FilterCoef[sub_pos][FILTER_SIZE*fj+fi]*sum;
    		}
        result += FilterCoef[sub_pos][SQR_FILTER-1]; 
   			is = (int)(0.5+result);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        block[i][j] = max(0, min(((1<<img->bitdepth_luma)-1),is));
#else
        block[i][j] = max(0, min(255,is));
#endif
    }
  }
  else if(sub_pos == i_pos || sub_pos == k_pos)
  {
    if(FilterFlag[sub_pos])
    {
      start_x = 1; end_x = FILTER_SIZE-1;
      start_y = 1; end_y = FILTER_SIZE/2;
    }
    else
    {
      start_x = 0; end_x = FILTER_SIZE;
      start_y = 0; end_y = FILTER_SIZE/2;
    }
  	for(i = 0; i < BLOCK_SIZE; i++)
	  for(j = 0; j < BLOCK_SIZE; j++)
		{
				result = 0.0;
    		for(fi = start_x; fi < end_x; fi++)
    		for(fj = start_y; fj < end_y; fj++)
    		{
    			fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
    			fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
          sum = list[ref_idx]->imgY[fpos_y][fpos_x];
    			//fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
    			fpos_y = FindPosition(img_height, j+y_pos, FILTER_SIZE-1-fj,mvy);
          sum += list[ref_idx]->imgY[fpos_y][fpos_x];
    			result += FilterCoef[sub_pos][FILTER_SIZE*fj+fi]*sum;
    		}
        result += FilterCoef[sub_pos][SQR_FILTER-1]; 
   			is = (int)(0.5+result);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        block[i][j] = max(0, min(((1<<img->bitdepth_luma)-1),is));
#else
        block[i][j] = max(0, min(255,is));
#endif
    }
  }
  else if(sub_pos == j_pos)
  {
    if(FilterFlag[sub_pos])
    {
      start_x = 1; end_x = FILTER_SIZE/2;
      start_y = 1; end_y = FILTER_SIZE/2;
    }
    else
    {
      start_x = 0; end_x = FILTER_SIZE/2;
      start_y = 0; end_y = FILTER_SIZE/2;
    }
  	for(i = 0; i < BLOCK_SIZE; i++)
	  for(j = 0; j < BLOCK_SIZE; j++)
		{
				result = 0.0;
    		for(fi = start_x; fi < end_x; fi++)
    		for(fj = start_y; fj < end_y; fj++)
    		{
    			fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
    			fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
          sum = list[ref_idx]->imgY[fpos_y][fpos_x];
    			fpos_x = FindPosition(img_width,  i+x_pos, FILTER_SIZE-1-fi,mvx);
    			//fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
          sum += list[ref_idx]->imgY[fpos_y][fpos_x];
    			//fpos_x = FindPosition(img_width,  i+x_pos, FILTER_SIZE-1-fi,mvx);
    			fpos_y = FindPosition(img_height, j+y_pos, FILTER_SIZE-1-fj,mvy);
          sum += list[ref_idx]->imgY[fpos_y][fpos_x];
    			fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
    			//fpos_y = FindPosition(img_height, j+y_pos, FILTER_SIZE-1-fj,mvy);
          sum += list[ref_idx]->imgY[fpos_y][fpos_x];
    			result += FilterCoef[sub_pos][FILTER_SIZE*fj+fi]*sum;
    		}
        result += FilterCoef[sub_pos][SQR_FILTER-1]; 
   			is = (int)(0.5+result);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        block[i][j] = max(0, min(((1<<img->bitdepth_luma)-1),is));
#else
        block[i][j] = max(0, min(255,is));
#endif
    }
  }
  else if(FilterFlagInt)  // full-pel position w/ filtering
  {
  	for(i = 0; i < BLOCK_SIZE; i++)
	  for(j = 0; j < BLOCK_SIZE; j++)
    {
        result = 0.0;
#ifdef COEFF_SYM_FULLPEL
        for(fi = 0; fi < FILTER_SIZE_INT/2; fi++)
        for(fj = 0; fj < FILTER_SIZE_INT/2; fj++)
        {
          fpos_x = FindPositionInt(img_width,  i+x_pos, fi,mvx);
          fpos_y = FindPositionInt(img_height, j+y_pos, fj,mvy);
          sum = list[ref_idx]->imgY[fpos_y][fpos_x];
          fpos_x = FindPositionInt(img_width,  i+x_pos, FILTER_SIZE_INT-1-fi,mvx);
          //fpos_y = FindPositionInt(img_height, j+y_pos, fj,mvy);
          sum += list[ref_idx]->imgY[fpos_y][fpos_x];
          //fpos_x = FindPositionInt(img_width,  i+x_pos, FILTER_SIZE_INT-1-fi,mvx);
          fpos_y = FindPositionInt(img_height, j+y_pos, FILTER_SIZE_INT-1-fj,mvy);
          sum += list[ref_idx]->imgY[fpos_y][fpos_x];
          fpos_x = FindPositionInt(img_width,  i+x_pos, fi,mvx);
          //fpos_y = FindPositionInt(img_height, j+y_pos, FILTER_SIZE_INT-1-fj,mvy);
          sum += list[ref_idx]->imgY[fpos_y][fpos_x];
          result += FilterCoefInt[FILTER_SIZE_INT*fj+fi]*sum; 
        }

        fi = FILTER_SIZE_INT/2; 
        for(fj = 0; fj < FILTER_SIZE_INT/2; fj++)
        {
          fpos_x = FindPositionInt(img_width,  i+x_pos, fi,mvx);
          fpos_y = FindPositionInt(img_height, j+y_pos, fj,mvy);
          sum = list[ref_idx]->imgY[fpos_y][fpos_x];
          //fpos_x = FindPositionInt(img_width,  i+x_pos, fi,mvx);
          fpos_y = FindPositionInt(img_height, j+y_pos, FILTER_SIZE_INT-1-fj,mvy);
          sum += list[ref_idx]->imgY[fpos_y][fpos_x];
          result += FilterCoefInt[FILTER_SIZE_INT*fj+fi]*sum;
        }

        fj = FILTER_SIZE_INT/2; 
        for(fi = 0; fi < FILTER_SIZE_INT/2; fi++)
        {
          fpos_x = FindPositionInt(img_width,  i+x_pos, fi,mvx);
          fpos_y = FindPositionInt(img_height, j+y_pos, fj,mvy);
          sum = list[ref_idx]->imgY[fpos_y][fpos_x];
          fpos_x = FindPositionInt(img_width,  i+x_pos, FILTER_SIZE_INT-1-fi,mvx);
          //fpos_y = FindPositionInt(img_height, j+y_pos, fj,mvy);
          sum += list[ref_idx]->imgY[fpos_y][fpos_x];
          result += FilterCoefInt[FILTER_SIZE_INT*fj+fi]*sum;
        }

        fi = FILTER_SIZE_INT/2; 
        fj = FILTER_SIZE_INT/2; 
        fpos_x = FindPositionInt(img_width,  i+x_pos, fi,mvx);
        fpos_y = FindPositionInt(img_height, j+y_pos, fj,mvy);
        sum = list[ref_idx]->imgY[fpos_y][fpos_x];
        result += FilterCoefInt[FILTER_SIZE_INT*fj+fi]*sum;
#else
        for(fi = 0; fi < FILTER_SIZE_INT; fi++)
        for(fj = 0; fj < FILTER_SIZE_INT; fj++)
        {
          fpos_x = FindPositionInt(img_width,  i+x_pos, fi,mvx);
          fpos_y = FindPositionInt(img_height, j+y_pos, fj,mvy);
          result += FilterCoefInt[FILTER_SIZE_INT*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]; 
        }
#endif 
        result += FilterCoefInt[SQR_FILTER_INT-1]; 
        is = (int)(0.5+result);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        block[i][j] = max(0, min(((1<<img->bitdepth_luma)-1),is));
#else
        block[i][j] = max(0, min(255,is));
#endif
    }
  }
  else  // full-pel position w/o filtering
  {
  	for(i = 0; i < BLOCK_SIZE; i++)
	  for(j = 0; j < BLOCK_SIZE; j++)
    {
				is = 	list[ref_idx]->imgY[max(0, min(img_height-1,j+y_pos+mvy/4))]
                                 [max(0, min(img_width -1,i+x_pos+mvx/4))];
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        block[i][j] = max(0, min(((1<<img->bitdepth_luma)-1),is));
#else
        block[i][j] = max(0, min(255,is));
#endif
    }
  }
}

#endif 

// separable aif (BEGIN)
void GetBlockWithSepAIF(int             ref_idx,
                        StorablePicture **list,
                        int             x_pos,
                        int             y_pos,
                        int             xoff4x4,
                        int             yoff4x4,
                        int             mvx,
                        int             mvy,
                        int             img_width,
                        int             img_height,
                        int             block_size_x,
                        int             block_size_y,
                        int             block[MB_BLOCK_SIZE][MB_BLOCK_SIZE])
{
  if(img->mv_res)
    GetBlockWithSepAIF_eighth_pel (ref_idx, list, x_pos, y_pos, xoff4x4, yoff4x4, mvx, mvy, img_width,img_height, block_size_x, block_size_y, block); 
  else
    GetBlockWithSepAIF_quarter_pel(ref_idx, list, x_pos, y_pos, xoff4x4, yoff4x4, mvx, mvy, img_width,img_height, block_size_x, block_size_y, block); 
}
// separable aif (END)


/*!
************************************************************************
* \brief
*    Interpolation of 1/4 subpixel with 2D adaptive filter
************************************************************************
*/
void GetBlockWith2DAIF_quarter_pel(int             ref_idx,
                                   StorablePicture **list,
                                   int             x_pos,
                                   int             y_pos,
                                   int             mvx,
                                   int             mvy,
                                   int             img_width,
                                   int             img_height,
                                   int             block[BLOCK_SIZE][BLOCK_SIZE])
#else
void GetBlockWith2DAIF(int             ref_idx,
                       StorablePicture **list,
                       int             x_pos,
                       int             y_pos,
                       int             mvx,
                       int             mvy,
                       int             img_width,
                       int             img_height,
                       int             block[BLOCK_SIZE][BLOCK_SIZE])
#endif
{
  int i,j, fi, fj;
  int fpos_x, fpos_y;
  double result;
  int is;
  int sub_pos, mvx_sub, mvy_sub;
  if(mvx >= 0)
    mvx_sub = mvx%4;                  // x-sub-coordinate in a 4x4block
  else
    mvx_sub = (4-abs(mvx)%4)%4;
  if(mvy >= 0)
    mvy_sub = mvy%4;                  // y-sub-coordinate in a 4x4block
  else
    mvy_sub = (4-abs(mvy)%4)%4;
  sub_pos = mvx_sub + 4*mvy_sub-1;    // pos 0..14 in a 4x4 block

  for(i = 0; i < BLOCK_SIZE; i++)
    for(j = 0; j < BLOCK_SIZE; j++)
    {
      if(sub_pos != -1)
      {
        result = 0.0;
        for(fi = 0; fi < FILTER_SIZE; fi++)
          for(fj = 0; fj < FILTER_SIZE; fj++)
          {
            fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
            fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
            result += FilterCoef[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x];
          }
#ifdef E_DAIF
          result += FilterCoef[sub_pos][SQR_FILTER-1];
#endif  // E_DAIF
          is = (int)(0.5+result);
      }
#ifdef E_DAIF
      else if(FilterFlagInt)  // full-pel position w/ filtering
      {
        result = 0.0;
        for(fi = 0; fi < FILTER_SIZE_INT; fi++)
          for(fj = 0; fj < FILTER_SIZE_INT; fj++)
          {
            fpos_x = FindPositionInt(img_width,  i+x_pos, fi,mvx);
            fpos_y = FindPositionInt(img_height, j+y_pos, fj,mvy);
            result += FilterCoefInt[FILTER_SIZE_INT*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]; 
          }
          result += FilterCoefInt[SQR_FILTER_INT-1]; 
          is = (int)(0.5+result);
      }
#endif
      else
      {
        is =  list[ref_idx]->imgY[max(0, min(img_height-1,j+y_pos+mvy/4))]
        [max(0, min(img_width -1,i+x_pos+mvx/4))];
      }
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
      block[i][j] = max(0, min(((1<<img->bitdepth_luma)-1),is));
#else
      block[i][j] = max(0, min(255,is));
#endif
    }
}
#ifdef DIRECTIONAL_FILTER
void GetBlockWith1DAIF( int             ref_idx,
                       StorablePicture **list,
                       int             x_pos,
                       int             y_pos,
                       int             mvx,
                       int             mvy,
                       int             img_width,
                       int             img_height,
                       int             block[BLOCK_SIZE][BLOCK_SIZE])
{
#ifdef EIGHTH_PEL
  if(img->mv_res)
    //GetBlockWith2DAIF_Sep_eighth_pel (ref_idx, list, x_pos, y_pos, mvx, mvy, img_width,img_height, block);
    printf("There is no 1D-AIF interpolation for decoder yet!\n");
  else
#endif
  {
    if (img->ImpType == IMP_FLOAT32)
      GetBlockWith1DAIF_quarter_pel(ref_idx, list, x_pos, y_pos, mvx, mvy, img_width,img_height, block);
    else if (img->ImpType == IMP_INT16)
      GetBlockWith1DAIF_16bits__quarter_pel(ref_idx, list, x_pos, y_pos, mvx, mvy, img_width,img_height, block);
  }
}
#endif
// separable aif (END)
/*!
************************************************************************
* \brief
*    Interpolation of 1/4 subpixel with 1D adaptive filter
************************************************************************
*/
#ifdef DIRECTIONAL_FILTER
void GetBlockWith1DAIF_quarter_pel(int             ref_idx,
                                   StorablePicture **list,
                                   int             x_pos,
                                   int             y_pos,
                                   int             mvx,
                                   int             mvy,
                                   int             img_width,
                                   int             img_height,
                                   int             block[BLOCK_SIZE][BLOCK_SIZE])
{
  int i,j, fi, fj;
  int fpos_x, fpos_y;
  int sub_pos, mvx_sub, mvy_sub;

  double result;

  if(mvx >= 0)
    mvx_sub = mvx%4;                  // x-sub-coordinate in a 4x4block
  else
    mvx_sub = (4-abs(mvx)%4)%4;
  if(mvy >= 0)
    mvy_sub = mvy%4;                  // y-sub-coordinate in a 4x4block
  else
    mvy_sub = (4-abs(mvy)%4)%4;
  sub_pos = mvx_sub + 4*mvy_sub-1;    // pos 0..14 in a 4x4 block


  for(i = 0; i < BLOCK_SIZE; i++)
    for(j = 0; j < BLOCK_SIZE; j++)
    {
      if(sub_pos != -1)
      {
        result = 0;
        if((sub_pos==a_pos)||(sub_pos==b_pos)||(sub_pos==c_pos))
          // horizontal filtering
        {
          fj = FILTER_OFFSET;
          fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
          for(fi = 0; fi < FILTER_SIZE; fi++)
          {
            fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
            result += (FilterCoef[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]);
          }
#ifdef E_DAIF
          result += FilterCoef[sub_pos][SQR_FILTER-1];
#endif  // E_DAIF
        }
        else if((sub_pos==d_pos)||(sub_pos==h_pos)||(sub_pos==l_pos))
          // vertical filtering
        {
          fi = FILTER_OFFSET;
          fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
          for(fj = 0; fj < FILTER_SIZE; fj++)
          {
            fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
            result += (FilterCoef[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]);
          }
#ifdef E_DAIF
          result += FilterCoef[sub_pos][SQR_FILTER-1];
#endif  // E_DAIF
        }
        else if ((sub_pos==e_pos)||(sub_pos==o_pos))
        {
#ifdef E_DAIF
          if(img->AdaptiveFilterFlag == FILTER_TYPE_EDAIF && !FilterFlag[sub_pos])
          {
            for(fi = 0; fi < FILTER_SIZE; fi++)
              for(fj = 0; fj < FILTER_SIZE; fj++)
              {
                fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
                fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
                result += FilterCoef[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x];
              }
          }
          else
#endif
          for(fj = 0; fj < FILTER_SIZE; fj++) // horizontal.
          {
            fi = fj;
            fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
            fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
            result += (FilterCoef[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]);
          }
#ifdef E_DAIF
          result += FilterCoef[sub_pos][SQR_FILTER-1];
#endif  // E_DAIF
        }
        else if ((sub_pos==g_pos)||(sub_pos==m_pos))
        {
#ifdef E_DAIF
          if(img->AdaptiveFilterFlag == FILTER_TYPE_EDAIF && !FilterFlag[sub_pos])
          {
            for(fi = 0; fi < FILTER_SIZE; fi++)
              for(fj = 0; fj < FILTER_SIZE; fj++)
              {
                fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
                fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
                result += FilterCoef[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x];
              }
          }
          else
#endif
          for(fi = 0; fi < FILTER_SIZE; fi++) // horizontal.
          {
            fj = FILTER_SIZE-1 - fi;
            fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
            fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
            result += (FilterCoef[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]);
          }
#ifdef E_DAIF
          result += FilterCoef[sub_pos][SQR_FILTER-1];
#endif  // E_DAIF
        }
        else if ((sub_pos==f_pos)||(sub_pos==n_pos)||(sub_pos==i_pos)||(sub_pos==k_pos)||(sub_pos==j_pos))
        {
#ifdef E_DAIF
          if(img->AdaptiveFilterFlag == FILTER_TYPE_EDAIF && !FilterFlag[sub_pos])
          {
            for(fi = 0; fi < FILTER_SIZE; fi++)
              for(fj = 0; fj < FILTER_SIZE; fj++)
              {
                fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
                fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
                result += FilterCoef[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x];
              }
          }
          else
#endif
          for(fi = 0; fi < FILTER_SIZE; fi++)
          {
            // We scale sum of every 2 samples beeing accumulated [fi,fi] + [fi+5(4,3),fi]

            fj = fi;
            fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
            fpos_y = FindPosition(img_height, j+y_pos, fi,mvy);
            result += list[ref_idx]->imgY[fpos_y][fpos_x]*FilterCoef[sub_pos][FILTER_SIZE*fi+fi];

            fj = FILTER_SIZE-1 - fi;
            fpos_x = FindPosition(img_width,  i+x_pos, fj,mvx);
            result += list[ref_idx]->imgY[fpos_y][fpos_x]*FilterCoef[sub_pos][FILTER_SIZE*fi+fj];
          }
#ifdef E_DAIF
          result += FilterCoef[sub_pos][SQR_FILTER-1];
#endif  // E_DAIF
        }
        result = (int)(result+0.5);
      }
#ifdef E_DAIF
      else if(FilterFlagInt)  // full-pel position w/ filtering
      {
        result = 0.0;
        for(fi = 0; fi < FILTER_SIZE_INT; fi++)
          for(fj = 0; fj < FILTER_SIZE_INT; fj++)
          {
            fpos_x = FindPositionInt(img_width,  i+x_pos, fi,mvx);
            fpos_y = FindPositionInt(img_height, j+y_pos, fj,mvy);
            result += FilterCoefInt[FILTER_SIZE_INT*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]; 
          }
        result += FilterCoefInt[SQR_FILTER_INT-1]; 
        result = (int)(0.5+result);
      }
#endif
      else
      {
        result = list[ref_idx]->imgY[max(0, min(img_height-1,j+y_pos+mvy/4))]
        [max(0, min(img_width -1,i+x_pos+mvx/4))];
      }

#ifdef  INTERNAL_BIT_DEPTH_INCREASE
      block[i][j] = (imgpel)max(0, min(((1<<img->bitdepth_luma)-1),result));
#else
      block[i][j] = (imgpel)max(0, min(255,result));
#endif

    }
}
#endif


#ifdef EIGHTH_PEL
// separable aif (BEGIN)
/*!
************************************************************************
* \brief
*    Interpolation of 1/4 subpixel with separable adaptive filter
************************************************************************
*/
void GetBlockWithSepAIF_quarter_pel( int             ref_idx,
                                    StorablePicture **list,
                                    int             x_pos,
                                    int             y_pos,
                                    int             xoff4x4,
                                    int             yoff4x4,
                                    int             mvx,
                                    int             mvy,
                                    int             img_width,
                                    int             img_height,
                                    int             block_x,
                                    int             block_y,
                                    int             block[MB_BLOCK_SIZE][MB_BLOCK_SIZE])

#else
void GetBlockWithSepAIF( int             ref_idx,
                        StorablePicture **list,
                        int             x_pos,
                        int             y_pos,
                        int             xoff4x4,
                        int             yoff4x4,
                        int             mvx,
                        int             mvy,
                        int             img_width,
                        int             img_height,
                        int             block_x,
                        int             block_y,
                        int             block[MB_BLOCK_SIZE][MB_BLOCK_SIZE])

#endif
{
  int block_y_ext = block_y + FILTER_SIZE - 1;
  int i,j, fi, fj, is;
  int fpos_x, fpos_y;
  double result=0.0;
  int sub_pos, mvx_sub, mvy_sub;
  if(mvx >= 0)
    mvx_sub = mvx%4;                   // x-sub-coordinate in a 4x4block
  else
    mvx_sub = (4-abs(mvx)%4)%4;
  if(mvy >= 0)
    mvy_sub = mvy%4;                  // y-sub-coordinate in a 4x4block
  else
    mvy_sub = (4-abs(mvy)%4)%4;
  sub_pos = mvx_sub + (mvy_sub<<2)-1;    // pos 0..14 in a 4x4 block

  if(sub_pos == -1) //full-pel pos
  {
    for(i = 0; i < block_x; i++)
      for(j = 0; j < block_y; j++)
        block[i+xoff4x4][j+yoff4x4] = list[ref_idx]->imgY[max(0, min(img_height-1,j+y_pos+mvy/4))][max(0, min(img_width-1,i+x_pos+mvx/4))];
  }
  else 
  {
    if (sub_pos < d_pos)
    {
      for(j = 0; j < block_y; j++)
      {
        fpos_y = FindPosition(img_height, j+y_pos, FILTER_OFFSET, mvy);
        for(i = 0; i < block_x; i++)
        {
          result = 0.0;
          for(fi = 0; fi < FILTER_SIZE; fi++)
          {
            fpos_x = FindPosition(img_width,  i+x_pos, fi, mvx);
            result += FilterCoef[sub_pos][fi]*list[ref_idx]->imgY[fpos_y][fpos_x];
          }
          is = (int)(0.5+result);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          block[i+xoff4x4][j+yoff4x4] = max(0, min(((1<<img->bitdepth_luma)-1),is));
#else
          block[i+xoff4x4][j+yoff4x4] = max(0, min(255,is));
#endif
        }
      }
    }
    else if (sub_pos == d_pos || sub_pos == h_pos || sub_pos == l_pos)
    {
      for(i = 0; i < block_x; i++)
      {
        fpos_x = FindPosition(img_width,  i+x_pos, FILTER_OFFSET,mvx);
        for(j = 0; j < block_y; j++)
        {
          result = 0.0;
          for(fj = 0; fj < FILTER_SIZE; fj++)
          {
            fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
            result += FilterCoef[sub_pos][fj]*list[ref_idx]->imgY[fpos_y][fpos_x];
          }
          is = (int)(0.5+result);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          block[i+xoff4x4][j+yoff4x4] = max(0, min(((1<<img->bitdepth_luma)-1),is));
#else
          block[i+xoff4x4][j+yoff4x4] = max(0, min(255,is));
#endif
        }
      }
    }
    else if (sub_pos == e_pos || sub_pos == i_pos || sub_pos == m_pos)
    {
      for(j = 0; j < block_y_ext; j++)
      {
        fpos_y = FindPosition(img_height, j+y_pos, 0, mvy);
        for(i = 0; i < block_x; i++)
        {
          result = 0.0;
          for(fi = 0; fi < FILTER_SIZE; fi++)
          {
            fpos_x = FindPosition(img_width,  i+x_pos, fi, mvx);
            result += FilterCoef[a_pos][fi]*list[ref_idx]->imgY[fpos_y][fpos_x];
          }
          tmp_coef[j][i] = result;
        } 
      }	
      for(j = 0; j < block_y; j++)
      {
        for(i = 0; i < block_x; i++)
        {
          result = 0.0;
          for(fj = 0; fj < FILTER_SIZE; fj++)
          {
            result += FilterCoef[sub_pos][fj]*tmp_coef[j+fj][i];
          }
          is = (int)(0.5+result);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          block[i+xoff4x4][j+yoff4x4] = max(0, min(((1<<img->bitdepth_luma)-1),is));
#else
          block[i+xoff4x4][j+yoff4x4] = max(0, min(255,is));
#endif
        }
      }
    }
    else if (sub_pos == f_pos || sub_pos == j_pos || sub_pos == n_pos)
    {
      for(j = 0; j < block_y_ext; j++)
      {
        fpos_y = FindPosition(img_height, j+y_pos, 0, mvy);
        for(i = 0; i < block_x; i++)
        {
          result = 0.0;
          for(fi = 0; fi < FILTER_SIZE; fi++)
          {
            fpos_x = FindPosition(img_width,  i+x_pos, fi, mvx);
            result += FilterCoef[b_pos][fi]*list[ref_idx]->imgY[fpos_y][fpos_x];
          }
          tmp_coef[j][i] = result;
        } 
      }	
      for(j = 0; j < block_y; j++)
      {
        for(i = 0; i < block_x; i++)
        {
          result = 0.0;
          for(fj = 0; fj < FILTER_SIZE; fj++)
          {
            result += FilterCoef[sub_pos][fj]*tmp_coef[j+fj][i];
          }
          is = (int)(0.5+result);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          block[i+xoff4x4][j+yoff4x4] = max(0, min(((1<<img->bitdepth_luma)-1),is));
#else
          block[i+xoff4x4][j+yoff4x4] = max(0, min(255,is));
#endif
        }
      }
    }
    else // sup_pos == g_pos || sub_pos == k_pos || sub_pos == o_pos
    {
      for(j = 0; j < block_y_ext; j++)
      {
        fpos_y = FindPosition(img_height, j+y_pos, 0, mvy);
        for(i = 0; i < block_x; i++)
        {
          result = 0.0;
          for(fi = 0; fi < FILTER_SIZE; fi++)
          {
            fpos_x = FindPosition(img_width,  i+x_pos, fi, mvx);
            result += FilterCoef[c_pos][fi]*list[ref_idx]->imgY[fpos_y][fpos_x];
          }
          tmp_coef[j][i] = result;
        } 
      }	
      for(j = 0; j < block_y; j++)
      {
        for(i = 0; i < block_x; i++)
        {
          result = 0.0;
          for(fj = 0; fj < FILTER_SIZE; fj++)
          {
            result += FilterCoef[sub_pos][fj]*tmp_coef[j+fj][i];
          }
          is = (int)(0.5+result);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          block[i+xoff4x4][j+yoff4x4] = max(0, min(((1<<img->bitdepth_luma)-1),is));
#else
          block[i+xoff4x4][j+yoff4x4] = max(0, min(255,is));
#endif
        }
      }
    }
  }
}
// separable aif (END)
#ifdef DIRECTIONAL_FILTER
void GetBlockWith1DAIF_16bits__quarter_pel(int             ref_idx,
                                           StorablePicture **list,
                                           int             x_pos,
                                           int             y_pos,
                                           int             mvx,
                                           int             mvy,
                                           int             img_width,
                                           int             img_height,
                                           int             block[BLOCK_SIZE][BLOCK_SIZE])
{
  int i,j, fi, fj;
  int fpos_x, fpos_y;
  int sub_pos, mvx_sub, mvy_sub;

  short int result_16bits, temp_accum,temp_accum2;
  imgpel result_uns16bits = 0;
  imgpel scale_shift;
  imgpel FQ = FILTCOEF_BITS;
  short int shift_half = (1 << (FQ-1));


  if(mvx >= 0)
    mvx_sub = mvx%4;                  // x-sub-coordinate in a 4x4block
  else
    mvx_sub = (4-abs(mvx)%4)%4;
  if(mvy >= 0)
    mvy_sub = mvy%4;                  // y-sub-coordinate in a 4x4block
  else
    mvy_sub = (4-abs(mvy)%4)%4;
  sub_pos = mvx_sub + 4*mvy_sub-1;    // pos 0..14 in a 4x4 block

  for(i = 0; i < BLOCK_SIZE; i++)
    for(j = 0; j < BLOCK_SIZE; j++)
    {
      if(sub_pos != -1)
      {
        FQ = nQBits[sub_pos];
        shift_half = (1 << (FQ-1));
        result_16bits = 0;
        scale_shift = 0; 
        if((sub_pos==a_pos)||(sub_pos==b_pos)||(sub_pos==c_pos))
          // horizontal filtering
        {
          fj = FILTER_OFFSET;
          fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
          temp_accum = 0;
          for(fi = 0; fi < FILTER_SIZE/2; fi++)
          {
            fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
            temp_accum += (FilterCoef_16bits[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]);
          }
          result_16bits = (short int)max(0,temp_accum);
          temp_accum = 0;
          for(fi = FILTER_SIZE/2; fi < FILTER_SIZE; fi++)
          {
            fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
            temp_accum += (FilterCoef_16bits[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]);
          }
          temp_accum = (short int)max(0,temp_accum);
          result_uns16bits = (imgpel)(result_16bits + temp_accum);

        }
        else if((sub_pos==d_pos)||(sub_pos==h_pos)||(sub_pos==l_pos))
          // vertical filtering
        {
          fi = FILTER_OFFSET;
          fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
          temp_accum = 0;
          for(fj = 0; fj < FILTER_SIZE/2; fj++)
          {
            fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
            temp_accum += (FilterCoef_16bits[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]);
          }
          temp_accum = (short int)max(0,temp_accum);
          result_16bits = temp_accum;
          temp_accum = 0;
          for(fj = FILTER_SIZE/2; fj < FILTER_SIZE; fj++)
          {
            fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
            temp_accum += (FilterCoef_16bits[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]);
          }
          temp_accum = (short int)max(0,temp_accum);
          result_uns16bits = (imgpel)(result_16bits + temp_accum);
        }
        else if ((sub_pos==e_pos)||(sub_pos==o_pos))
        {
          temp_accum = 0;
          for(fj = 0; fj < FILTER_SIZE/2; fj++) // horizontal.
          {
            fi = fj;
            fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
            fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
            temp_accum += (FilterCoef_16bits[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]);
          }
          temp_accum = (short int)max(0,temp_accum);
          result_16bits = temp_accum;

          temp_accum = 0;
          for(fj = FILTER_SIZE/2; fj < FILTER_SIZE; fj++) // horizontal.
          {
            fi = fj;
            fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
            fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
            temp_accum += (FilterCoef_16bits[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]);
          }
          temp_accum = (short int)max(0,temp_accum);
          result_uns16bits = (imgpel)(result_16bits + temp_accum);
        }
        else if ((sub_pos==g_pos)||(sub_pos==m_pos))
        {
          temp_accum = 0;
          for(fi = 0; fi < FILTER_SIZE/2; fi++) // horizontal.
          {
            fj = FILTER_SIZE-1 - fi;
            fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
            fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
            temp_accum += (FilterCoef_16bits[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]);
          }
          temp_accum = (short int)max(0,temp_accum);
          result_16bits = temp_accum;

          temp_accum = 0;
          for(fi = FILTER_SIZE/2; fi < FILTER_SIZE; fi++) // horizontal.
          {
            fj = FILTER_SIZE-1 - fi;
            fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
            fpos_y = FindPosition(img_height, j+y_pos, fj,mvy);
            temp_accum += (FilterCoef_16bits[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]);
          }
          temp_accum = (short int)max(0,temp_accum);
          result_uns16bits = (imgpel)(result_16bits + temp_accum);
        }
        else if ((sub_pos==f_pos)||(sub_pos==n_pos)||(sub_pos==i_pos)||(sub_pos==k_pos)||(sub_pos==j_pos))
        {
          result_16bits= 0;
          temp_accum = temp_accum2 = 0;
          for(fi = 0; fi < FILTER_SIZE/2; fi++) 
          {
            // We scale sum of every 2 samples beeing accumulated [fi,fi] + [fi+5(4,3),fi]
            fj = fi;
            fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
            fpos_y = FindPosition(img_height, j+y_pos, fi,mvy);
            //temp_accum += (FilterCoef_16bits[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]);
            temp_accum += list[ref_idx]->imgY[fpos_y][fpos_x]*FilterCoef_16bits[sub_pos][FILTER_SIZE*fi+fi];

            fj = FILTER_SIZE-1 - fi;
            fpos_x = FindPosition(img_width,  i+x_pos, fj,mvx);
            //fpos_y = FindPosition(img_height, j+y_pos, fi,mvy);
            //temp_accum += (FilterCoef_16bits[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]);
            temp_accum2 += list[ref_idx]->imgY[fpos_y][fpos_x]*FilterCoef_16bits[sub_pos][FILTER_SIZE*fi+fj];
          }
          temp_accum = (short int)max(0,temp_accum);
          temp_accum2 = (short int)max(0,temp_accum2);
          result_uns16bits =  ((imgpel)temp_accum  +  (imgpel)temp_accum2+0)>>1;
          temp_accum = temp_accum2 = 0;
          for(fi = FILTER_SIZE/2; fi < FILTER_SIZE; fi++)
          {
            // We scale sum of every 2 samples beeing accumulated [fi,fi] + [fi+5(4,3),fi]
            fj = fi;
            fpos_x = FindPosition(img_width,  i+x_pos, fi,mvx);
            fpos_y = FindPosition(img_height, j+y_pos, fi,mvy);
            //temp_accum += (FilterCoef_16bits[sub_pos][FILTER_SIZE*fj+fi]*list[ref_idx]->imgY[fpos_y][fpos_x]);
            temp_accum += list[ref_idx]->imgY[fpos_y][fpos_x]*FilterCoef_16bits[sub_pos][FILTER_SIZE*fi+fi];

            fj = FILTER_SIZE-1 - fi;
            fpos_x = FindPosition(img_width,  i+x_pos, fj,mvx);
            temp_accum2 += list[ref_idx]->imgY[fpos_y][fpos_x]*FilterCoef_16bits[sub_pos][FILTER_SIZE*fi+fj];
          }
          temp_accum = (short int)max(0,temp_accum);
          temp_accum2 = (short int)max(0,temp_accum2);
          result_uns16bits = result_uns16bits + (imgpel)(((imgpel)temp_accum  + (imgpel)temp_accum2+0)>>1);
          shift_half = shift_half>>1;
          FQ = FQ-1;
        }
        result_uns16bits = (result_uns16bits+shift_half) >> FQ;
      }
      else
      {
        result_uns16bits =  list[ref_idx]->imgY[max(0, min(img_height-1,j+y_pos+mvy/4))]
        [max(0, min(img_width -1,i+x_pos+mvx/4))];
      }

#ifdef  INTERNAL_BIT_DEPTH_INCREASE
      block[i][j] = (imgpel)max(0, min(((1<<img->bitdepth_luma)-1),result_uns16bits));
#else
      block[i][j] = (imgpel)max(0, min(255,result_uns16bits));
#endif

    }
    //PrintFilterCoefFloat(FilterCoef);
}
#endif

#ifdef E_DAIF
int FindPositionInt(int Dim, int cur_pos, int Offset, int mv)
{
  int position;
  mv =mv/4;
  position = cur_pos-FILTER_OFFSET_INT+Offset+mv;
  return max(0, min(Dim-1, position ) );
}
#endif  // E_DAIF

/*!
************************************************************************
* \brief
*    finds current position in X-direction
*    input:
*      Dim     - X/Y-size of image,
*      cur_pos - current X/Y-coordinate,
*      Offset  - offset in X/Y-direction from left filter position (2D-filter)
*      mv      - mvx/y-motion vector
************************************************************************
*/

int FindPosition(int Dim, int cur_pos, int Offset, int mv)
{
  int position;
  if(mv < 0 && abs(mv)%4!=0)
    mv = mv - 4;
  mv = mv/4;
  position = cur_pos - FILTER_OFFSET + Offset + mv;
  return (max(0, min(Dim-1, position)));
}

/*!
************************************************************************
* \brief
*    decide whether to interpolate with a new or an old filter
*    input:
*     mvx - x-motion vector
*      mvy - y-motion vector
*    output:
*      0 - standard filter
*      1 - adaptive filter
************************************************************************
*/
int CountWithNewFilter(int mvx, int mvy)
{
  int sub_pos, mvx_sub, mvy_sub;
  if(mvx >= 0)
    mvx_sub = mvx%4;                  // x-sub-coordinate in a 4x4block
  else
    mvx_sub = (4-abs(mvx)%4)%4;
  if(mvy >= 0)
    mvy_sub = mvy%4;                  // y-sub-coordinate in a 4x4block
  else
    mvy_sub = (4-abs(mvy)%4)%4;
  sub_pos = mvx_sub + 4*mvy_sub-1;    // pos 0..14 in a 4x4 block
  if(sub_pos == -1)
    return 0;
  else if(FilterFlag[sub_pos] == 0)
    return 0;
  else
    return 1;
}


/*!
************************************************************************
* \brief
*    prints filter coefficient
************************************************************************
*/
void PrintFilterCoefFloat(double FCoef[15][SQR_FILTER])
{
  int i,j;
  for(i = a_pos; i <= o_pos; i++)
  {
    //    printf("Filter_%c=[\n",97+i);
    printf("FilterFlag[%c]=%d\n",97+i,FilterFlag[i]);
    printf("Filter(%d,1:6,1:6)=[\n",i+1);
    for(j = 0; j < SQR(FILTER_SIZE); j++)
    {
      if(!(j%FILTER_SIZE))
        printf("[");
      printf("%12.8f   ",FCoef[i][j]);
      if(j%FILTER_SIZE == FILTER_SIZE-1)
        printf("];\n");
      else
        printf(", ");
    }
    printf("];\n");
  }
}

void PrintFilterCoefFloatFile(double FCoef[15][SQR_FILTER])
{
  int i,j;
  FILE* fid = fopen("dec_filter.txt","at");
  for(i = a_pos; i <= o_pos; i++)
  {
    fprintf(fid,"%c_pos: %d %d\n",97+i,FilterFlag[i], SymmetryPosition[i]);
    for(j = 0; j < SQR_FILTER; j++)
    {
      fprintf(fid, "  %1.6f \t",FCoef[i][j]);
      if(j%FILTER_SIZE == FILTER_SIZE-1)
        fprintf(fid, "\n");
    }
    fprintf(fid, "\n");
  }
  fclose(fid);
}


/*!
************************************************************************
* \brief
*    prints filter coefficient
************************************************************************
*/
void PrintFilterCoefInt(int FCoef[15][SQR_FILTER])
{
  int i,j;
  for(i = a_pos; i <= o_pos; i++)
  {
    printf("%c_pos:\n",97+i);
    for(j = 0; j < SQR_FILTER; j++)
    {
      if(FCoef[i][j] >= 0)
        printf(" ");
      printf("%d   \t",FCoef[i][j]);
      if(j%FILTER_SIZE == 5)
        printf("\n");
    }
    printf("========================================================\n");
  }
}


#ifdef EDAIF2
void ArrangeSymmetryEDAIF2(void)
{
  int * p_SymProf[4] = {&SetHor[0][0],&SetVert[0][0],&SetDiagNW[0][0],&SetDiagNE[0][0]};
  int * SymProf = NULL;
  int h1,h2;

  int j,k,counter;
  counter = 0;
  for(k = 0;k<4;k++)
  {
    SymProf  = p_SymProf[k];
    for (j = 0; j < nPairsSets[k]; j ++)
    {
      h1 = FILTER_NEW_SUB_POS[*(SymProf + j*2 +0)];
      h2 = *(SymProf+ 2*j + 1);
      if(realSymCommands[k][j])
      {
        FILTER_NEW_SUB_POS[h2] = h1;
        counter++;
      }
    }
  }

  for (j=0;j<=o_pos;j++)
  {
    if(counter>16)
      nBitsIntRepresent[j] = DEFAULT_QUANT-1;
    if(counter==20)
      nBitsIntRepresent[j] = DEFAULT_QUANT-2;

    if(FILTER_NEW_SUB_POS[j] != j)
      SymmetryPosition[j] = 0;
    else
      SymmetryPosition[j] = 1;
  }

}
#endif
#endif
