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
*    "Motion-and Aliasing-compensated Prediction using a two-dimensional non-separable 
*    Adaptive Wiener Interpolation Filter",
*    Proc. ICIP 2005, IEEE International Conference on Image Processing, Genova, Italy, September 2005.
*    Yuri Vatis, Bernd Edler, Ingolf Wassermann, Dieu Thanh Nguyen, Jörn Ostermann,
*    "Coding of Coefficients of two-dimensional non-separable Adaptive Wiener Interpolation Filter", 
*    Proc. VCIP 2005, SPIE Visual Communication & Image Processing, Beijing, China, July 2005.
************************************************************************
*/



#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "global.h"
#include "adaptive_filter.h"
#include "memalloc.h"
#include "malloc.h"
#include "defines.h"
#include "mbuffer.h"
#include "image.h"
#include "vlc.h"
#ifdef DIRECTIONAL_FILTER
#include <memory.h>
#include "adaptive_filter_1DAIF.h"
#include "adaptive_filter_orig.h"
#ifdef EDAIF2
#include "const_DAIF.h"
#include "const_RAIF.h"
#endif
#include "header.h"
#endif

#ifdef ADAPTIVE_FILTER
#ifdef EDAIF2
int TwoDEquationPattern[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER]=
#else
int TwoDEquationPattern[15][SQR_FILTER]  =    // get equation number for symmetric filter coefficients
#endif
{
  { 0,  1,  2,  3,  4,  5, 
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,        // a_pos
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0  },
  { 0,  1,  2,  2,  1,  0,
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,        // b_pos
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0  },
  { 5,  4,  3,  2,  1,  0,
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,        // c_pos
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0  },
  { 0,  1,  2,  3,  4,  5,
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,        // d_pos
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0  },
  { 0,  1,  2,  3,  4,  5,
  1,  6,  7,  8,  9, 10,
  2,  7, 11, 12, 13, 14,        // e_pos
  3,  8, 12, 15, 16, 17,
  4,  9, 13, 16, 18, 19,
  5, 10, 14, 17, 19, 20  },
  { 0,  1,  2,  2,  1,  0,
  3,  4,  5,  5,  4,  3,
  6,  7,  8,  8,  7,  6,        // f_pos
  9, 10, 11, 11, 10,  9,
  12, 13, 14, 14, 13, 12,
  15, 16, 17, 17, 16, 15  },
  { 5,  4,  3,  2,  1,  0,
  10,  9,  8,  7,  6,  1,
  14, 13, 12, 11,  7,  2,        // g_pos
  17, 16, 15, 12,  8,  3,
  19, 18, 16, 13,  9,  4,
  20, 19, 17, 14, 10,  5  },
  { 0,  1,  2,  2,  1,  0,
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,        // h_pos
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0  },
  { 0,  3,  6,  9, 12, 15,
  1,  4,  7, 10, 13, 16,
  2,  5,  8, 11, 14, 17,        // i_pos
  2,  5,  8, 11, 14, 17,
  1,  4,  7, 10, 13, 16,
  0,  3,  6,  9, 12, 15  },
  { 0,  1,  2,  2,  1,  0,
  1,  3,  4,  4,  3,  1,
  2,  4,  5,  5,  4,  2,         // j_pos
  2,  4,  5,  5,  4,  2,
  1,  3,  4,  4,  3,  1,
  0,  1,  2,  2,  1,  0  },
  {15, 12,  9,  6,  3,  0,
  16, 13, 10,  7,  4,  1,
  17, 14, 11,  8,  5,  2,         // k_pos
  17, 14, 11,  8,  5,  2,
  16, 13, 10,  7,  4,  1,
  15, 12,  9,  6,  3,  0  },
  { 5,  4,  3,  2,  1,  0,
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,         // l_pos
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0  },
  { 5, 10, 14, 17, 19, 20,
  4,  9, 13, 16, 18, 19,
  3,  8, 12, 15, 16, 17,         // m_pos
  2,  7, 11, 12, 13, 14,
  1,  6,  7,  8,  9, 10,
  0,  1,  2,  3,  4,  5  },
  {15, 16, 17, 17, 16, 15,
  12, 13, 14, 14, 13, 12,
  9, 10, 11, 11, 10,  9,         // n_pos
  6,  7,  8,  8,  7,  6,
  3,  4,  5,  5,  4,  3,
  0,  1,  2,  2,  1,  0  },
  {20, 19, 17, 14, 10,  5,
  19, 18, 16, 13,  9,  4,
  17, 16, 15, 12,  8,  3,         // o_pos
  14, 13, 12, 11,  7,  2,
  10,  9,  8,  7,  6,  1,
  5,  4,  3,  2,  1,  0  }
};
#ifdef EDAIF2
double STANDARD_2D_FILTER[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER]=
#else
double STANDARD_2D_FILTER[15][SQR_FILTER] = 
#endif
{
  { 0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  1.0/64.0  ,   -5.0/64.0  ,   52.0/64.0  ,   20.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0
  }, // a_pos
  { 0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  1.0/32.0  ,   -5.0/32.0  ,   20.0/32.0  ,   20.0/32.0  ,   -5.0/32.0  ,   1.0/32.0  ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0
  }, // b_pos
  { 0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  1.0/64.0  ,   -5.0/64.0  ,   20.0/64.0  ,   52.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0
  }, // c_pos
  { 0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   52.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0
  }, // d_pos
  { 0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  1.0/64.0  ,   -5.0/64.0  ,   40.0/64.0  ,   20.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
  0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0
  }, // e_pos
  { 1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0,
  -5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
  52.0/2048.0, -260.0/2048.0, 1040.0/2048.0, 1040.0/2048.0, -260.0/2048.0,  52.0/2048.0,
  20.0/2048.0, -100.0/2048.0,  400.0/2048.0,  400.0/2048.0, -100.0/2048.0,  20.0/2048.0,
  -5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
  1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0
  }, // f_pos
  { 0.0       ,    0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,   0.0       ,
  1.0/64.0  ,   -5.0/64.0  ,   20.0/64.0  ,   40.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
  0.0       ,    0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,   0.0
  }, // g_pos
  { 0.0       ,    0.0       ,    1.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   20.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   20.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    1.0/32.0  ,    0.0       ,    0.0       ,   0.0
  }, // h_pos
  { 1.0/2048.0,   -5.0/2048.0,   52.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0,
  -5.0/2048.0,   25.0/2048.0, -260.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
  20.0/2048.0, -100.0/2048.0, 1040.0/2048.0,  400.0/2048.0, -100.0/2048.0,  20.0/2048.0,
  20.0/2048.0, -100.0/2048.0, 1040.0/2048.0,  400.0/2048.0, -100.0/2048.0,  20.0/2048.0,
  -5.0/2048.0,   25.0/2048.0, -260.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
  1.0/2048.0,   -5.0/2048.0,   52.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0
  }, // i_pos
  { 1.0/1024.0,   -5.0/1024.0,   20.0/1024.0,   20.0/1024.0,   -5.0/1024.0,   1.0/1024.0,
  -5.0/1024.0,   25.0/1024.0, -100.0/1024.0, -100.0/1024.0,   25.0/1024.0,  -5.0/1024.0,
  20.0/1024.0, -100.0/1024.0,  400.0/1024.0,  400.0/1024.0, -100.0/1024.0,  20.0/1024.0,
  20.0/1024.0, -100.0/1024.0,  400.0/1024.0,  400.0/1024.0, -100.0/1024.0,  20.0/1024.0,
  -5.0/1024.0,   25.0/1024.0, -100.0/1024.0, -100.0/1024.0,   25.0/1024.0,  -5.0/1024.0,
  1.0/1024.0,   -5.0/1024.0,   20.0/1024.0,   20.0/1024.0,   -5.0/1024.0,   1.0/1024.0
  }, // j_pos
  { 1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   52.0/2048.0,   -5.0/2048.0,   1.0/2048.0,
  -5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -260.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
  20.0/2048.0, -100.0/2048.0,  400.0/2048.0, 1040.0/2048.0, -100.0/2048.0,  20.0/2048.0,
  20.0/2048.0, -100.0/2048.0,  400.0/2048.0, 1040.0/2048.0, -100.0/2048.0,  20.0/2048.0,
  -5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -260.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
  1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   52.0/2048.0,   -5.0/2048.0,   1.0/2048.0
  }, // i_pos
  { 0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   52.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0
  }, // l_pos
  { 0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  1.0/64.0  ,   -5.0/64.0  ,   40.0/64.0  ,   20.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
  0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0
  }, // m_pos
  { 1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0,
  -5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
  20.0/2048.0, -100.0/2048.0,  400.0/2048.0,  400.0/2048.0, -100.0/2048.0,  20.0/2048.0,
  52.0/2048.0, -260.0/2048.0, 1040.0/2048.0, 1040.0/2048.0, -260.0/2048.0,  52.0/2048.0,
  -5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
  1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0
  }, // n_pos
  { 0.0       ,    0.0       ,    0.0     ,      1.0/64.0  ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0     ,     -5.0/64.0  ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0     ,     20.0/64.0  ,    0.0       ,   0.0       ,
  1.0/64.0  ,   -5.0/64.0  ,   20.0/64.0,     40.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
  0.0       ,    0.0       ,    0.0     ,     -5.0/64.0  ,    0.0       ,   0.0       ,
  0.0       ,    0.0       ,    0.0     ,      1.0/64.0  ,    0.0       ,   0.0
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
#ifdef EDAIF2
int FILTER_NEW_SUB_POS[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA]=
#else
int FILTER_NEW_SUB_POS[SQR_FILTER] = 
#endif
{        a_pos, b_pos, a_pos,
a_pos, e_pos, f_pos, e_pos,
b_pos, f_pos, j_pos, f_pos,
a_pos, e_pos, f_pos, e_pos };

int TwoDSymmetricPattern[15][SQR_FILTER]  =    // get one filter from another one, if symmetry properties are used (used e.g. in ExtendFilterCoefficients)
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

#ifdef E_DAIF
int numQBitsInt[SQR_FILTER_INT-1] = 
{
  12, 12, 12, 12, 12, 
  12, 11, 11, 11, 12, 
  12, 11, 10, 11, 12, 
  12, 11, 11, 11, 12,
  12, 12, 12, 12, 12, 
};   

int UseEquation_DAIF[3][SQR_FILTER-1] = 
{
  {
    1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1,
  },
  {
    0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 0,
      0, 0, 0, 1, 0, 0,
      0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0,
    },
    {
      1, 0, 0, 0, 0, 1,
        0, 1, 0, 0, 1, 0,
        0, 0, 1, 1, 0, 0,
        0, 0, 1, 1, 0, 0,
        0, 1, 0, 0, 1, 0,
        1, 0, 0, 0, 0, 1,
    },
};
#endif

#ifdef EAIF
int UseEquation_EAIF[3][SQR_FILTER-1]  =
{
  {
    0, 0, 0, 0, 0, 0, 
      0, 0, 1, 1, 0, 0,
      0, 1, 1, 1, 1, 0,
      0, 1, 1, 1, 1, 0,
      0, 0, 1, 1, 0, 0,
      0, 0, 0, 0, 0, 0,
  },
  {
    0, 0, 0, 0, 0, 0, 
      0, 0, 1, 1, 0, 0,
      0, 1, 1, 1, 1, 0,
      0, 1, 1, 1, 1, 0,
      0, 0, 1, 1, 0, 0,
      0, 0, 0, 0, 0, 0,
    },
    {
      0, 0, 0, 0, 0, 0, 
        0, 0, 1, 1, 0, 0,
        0, 1, 1, 1, 1, 0,
        0, 1, 1, 1, 1, 0,
        0, 0, 1, 1, 0, 0,
        0, 0, 0, 0, 0, 0,
    },
};
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

FilterEntity SetAIF[3];
double costMCP[4][MAX_NUM_AIF+1]; 
int g_nBitsThresh;

int nPairsSets[4] = {4,4,6,6};
int SetHor[4][2] =     {
  {a_pos, c_pos}, 
  {e_pos, g_pos}, 
  {i_pos, k_pos}, 
  {m_pos, o_pos}, 
};
int SetVert[4][2] =     {
  {d_pos, l_pos}, 
  {e_pos, m_pos}, 
  {f_pos, n_pos}, 
  {g_pos, o_pos}, 
};

int SetDiagNW[6][2] =     {
  {a_pos, l_pos}, 
  {b_pos, h_pos}, 
  {c_pos, d_pos}, 
  {e_pos, o_pos}, 
  {f_pos, k_pos},
  {i_pos, n_pos}
};

int SetDiagNE[6][2] =     {
  {a_pos, d_pos}, 
  {b_pos, h_pos}, 
  {c_pos, l_pos}, 
  {g_pos, m_pos}, 
  {f_pos, i_pos},
  {k_pos, n_pos},
};
int realSymCommands[4][MAX_NUM_AIF];

int RealTimeDecomposition;
double g_AIFDecompThresh;
int  UseAllSubpelPositions;                  // 1 if FilterFlag for all independent positions is 1
int SubpelPositionsPattern;


int Calc2Filt_Indexes[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];
int Filt_Indexes_offset[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];

//		System of Linear equations coordinates
//int TwoDEquationPattern[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];
int FILTER_NEW_SUB_POS[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];

//		Storage for correlation matrixes
double NxN_Matrix[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER][SQR_FILTER];			// SQR(PEL)-1 X SQR(FILTER_SIZE) X SQR(FILTER_SIZE) - matrix
double N_Vector[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];										// SQR(PEL)-1 X SQR(FILTER_SIZE) - vector

//	Filter coefficient storages
double STANDARD_2D_FILTER[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];
double CalculatedFilter[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];							// temporal storage
double FilterCoef[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];									// estimated filter coefficients
short int FilterCoef16bits[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];	

double DiffFilterCoef[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];									
int DiffQFilterCoef[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];				// differences to be transmitted
int POS_EQUATION_NUMBER[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];								// number of different coefficietns for each sub-pel position
int FilterFlag[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];													// Flags defining if Filter at the particular position calculated

int SymmetryPosition[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];	// if 0, the position is copied from another one
int Is1DPosition[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];// if 0, the position is a 2D one
int IsDiagonal1D[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA]; // if 1, the filter alighned NW-SE, 2 - NE-SW
int nBitsIntRepresent[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];

#endif

int NumberOfQBits;                              // number of bits used for quantization of the filter coefficients

#ifndef EDAIF2
int DiffQFilterCoef[15][SQR_FILTER];        // differences to be transmitted
int POS_EQUATION_NUMBER[15];                // number of different coefficietns for each sub-pel position
int FilterFlag[15];                          // Flags defining if Filter at the particular position calculated
int SymmetryPosition[15] = {1,1,0,0,1,1,0,0,0,1,0,0,0,0,0};  // if 0, the position is copied from another one
int Is1DPosition[15] = {1,1,1,1,0,0,0,1,0,0,0,1,0,0,0}; // if 0, the position is a 2D one

#ifdef DIRECTIONAL_FILTER
int IsDiagonal1D[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // if 1 - the filter alighned NW-SE, 2 - NE-SW
short int FilterCoef16bits[15][SQR_FILTER];                  // estimated filter coefficients
int FILTER_SIZE = 6;    // is not checked for other values
int FILTER_OFFSET;
#endif

double NxN_Matrix[15][SQR_FILTER][SQR_FILTER];      // SQR(PEL)-1 X SQR(FILTER_SIZE) X SQR(FILTER_SIZE) - matrix
double N_Vector[15][SQR_FILTER];                    // SQR(PEL)-1 X SQR(FILTER_SIZE) - vector
double FilterCoef[15][SQR_FILTER];                  // estimated filter coefficients
double  DiffFilterCoef[15][SQR_FILTER];             // differences
double CalculatedFilter[15][SQR_FILTER];            // Filter for different positions
#endif
int nQBits[15] = {7,7,7,7,7,8,7,7,8,8,8,7,7,8,7};    // Number of bits for each coefficient representation, 16-bit DAIF
int    QFilterCoef[15][SQR_FILTER];                 // SQR(PEL)-1 X SQR(FILTER_SIZE) - quantized filter coefficients
double PredFilterCoef[15][SQR_FILTER];              // predicted filter coefficients

int **MBRateStorage;                                // storage for number of bits per macroblock for part. mode coded with standard filter
int *BitsPerFrameCounter;                           // storage for number of bits per macroblock for part. mode coded with standard filter
int number_of_macroblocks;
int current_poc;

double cost_frame[2];

#ifdef E_DAIF
double NxN_Matrix_Int[SQR_FILTER_INT][SQR_FILTER_INT];
double N_Vector_Int  [SQR_FILTER_INT];
double FilterCoefInt [SQR_FILTER_INT];
double CalculatedFilterInt [SQR_FILTER_INT];
int FilterFlagInt; 
int FILTER_OFFSET_INT;

int DiffQFilterOffsetI[15], DiffQFilterOffsetF[15];
int DiffQFilterOffsetIntI, DiffQFilterOffsetIntF; 
int DiffQFilterCoeffInt[SQR_FILTER_INT]; 
int UseEquation[3][SQR_FILTER-1];
#endif

#ifdef DIRECTIONAL_FILTER
/*! 
*************************************************************************************
* \brief
*   Initialize the predefined AIF constants. 
*  
* \para <title>
*    <paragraph>
*
* \para
*    <another paragraph>
*
* \author
*    - Dmytro Rusanovskyy                   <dmytro.rusanovskyy@tut.fi>
*************************************************************************************
*/
void SwapAIF(int filterID, int dir)
{
  FilterEntity *p_Filter = &SetAIF[filterID];
  if(dir)// copy from global
  {
    p_Filter->RealTimeDecomposition = RealTimeDecomposition;
    p_Filter->UseAllSubpelPositions = UseAllSubpelPositions;
    p_Filter->SubpelPositionsPattern = SubpelPositionsPattern;
    memcpy(p_Filter->realSymCommands,realSymCommands,sizeof(int)*4*MAX_NUM_AIF);

    memcpy(p_Filter->STANDARD_2D_FILTER,STANDARD_2D_FILTER,sizeof(double)*MAX_NUM_AIF*SQR_FILTER);
    memcpy(p_Filter->DiffQFilterCoef,DiffQFilterCoef,sizeof(int)*MAX_NUM_AIF*SQR_FILTER);
    memcpy(p_Filter->TwoDEquationPattern,TwoDEquationPattern,sizeof(int)*MAX_NUM_AIF*SQR_FILTER);
    memcpy(p_Filter->Calc2Filt_Indexes,Calc2Filt_Indexes,sizeof(int)*MAX_NUM_AIF*SQR_FILTER);
    memcpy(p_Filter->IsDiagonal1D,IsDiagonal1D,sizeof(int)*MAX_NUM_AIF);
    memcpy(p_Filter->Is1DPosition,Is1DPosition,sizeof(int)*MAX_NUM_AIF);
    memcpy(p_Filter->SymmetryPosition,SymmetryPosition,sizeof(int)*MAX_NUM_AIF);
    memcpy(p_Filter->POS_EQUATION_NUMBER,POS_EQUATION_NUMBER,sizeof(int)*MAX_NUM_AIF);
    memcpy(p_Filter->FILTER_NEW_SUB_POS,FILTER_NEW_SUB_POS,sizeof(int)*MAX_NUM_AIF);
    memcpy(p_Filter->Filt_Indexes_offset,Filt_Indexes_offset,sizeof(int)*MAX_NUM_AIF);
    memcpy(p_Filter->N_Vector,N_Vector, SQR_FILTER*(MAX_NUM_AIF)*sizeof(double));
    memcpy(p_Filter->NxN_Matrix,NxN_Matrix, SQR_FILTER*SQR_FILTER*(MAX_NUM_AIF)*sizeof(double));
    memcpy(p_Filter->FilterFlag,FilterFlag,sizeof(int)*MAX_NUM_AIF);
    memcpy(p_Filter->FilterCoef,FilterCoef, SQR_FILTER*(MAX_NUM_AIF)*sizeof(double));
    memcpy(p_Filter->FilterCoef16bits,FilterCoef16bits, SQR_FILTER*(MAX_NUM_AIF)*sizeof(short int));
    memcpy(p_Filter->nBitsIntRepresent,nBitsIntRepresent,sizeof(int)*MAX_NUM_AIF);
  }
  else // copy to global
  {
    RealTimeDecomposition = p_Filter->RealTimeDecomposition;
    UseAllSubpelPositions = p_Filter->UseAllSubpelPositions;
    SubpelPositionsPattern = p_Filter->SubpelPositionsPattern;
    memcpy(realSymCommands,p_Filter->realSymCommands,sizeof(int)*4*MAX_NUM_AIF);

    memcpy(STANDARD_2D_FILTER,p_Filter->STANDARD_2D_FILTER,sizeof(double)*MAX_NUM_AIF*SQR_FILTER);
    memcpy(DiffQFilterCoef,p_Filter->DiffQFilterCoef,sizeof(int)*MAX_NUM_AIF*SQR_FILTER);
    memcpy(TwoDEquationPattern,p_Filter->TwoDEquationPattern,sizeof(int)*MAX_NUM_AIF*SQR_FILTER);
    memcpy(Calc2Filt_Indexes,p_Filter->Calc2Filt_Indexes,sizeof(int)*MAX_NUM_AIF*SQR_FILTER);
    memcpy(IsDiagonal1D,p_Filter->IsDiagonal1D,sizeof(int)*MAX_NUM_AIF);
    memcpy(Is1DPosition,p_Filter->Is1DPosition,sizeof(int)*MAX_NUM_AIF);
    memcpy(SymmetryPosition,p_Filter->SymmetryPosition,sizeof(int)*MAX_NUM_AIF);
    memcpy(POS_EQUATION_NUMBER,p_Filter->POS_EQUATION_NUMBER,sizeof(int)*MAX_NUM_AIF);
    memcpy(FILTER_NEW_SUB_POS,p_Filter->FILTER_NEW_SUB_POS,sizeof(int)*MAX_NUM_AIF);
    memcpy(Filt_Indexes_offset,p_Filter->Filt_Indexes_offset,sizeof(int)*MAX_NUM_AIF);
    memcpy(N_Vector, p_Filter->N_Vector,SQR_FILTER*(MAX_NUM_AIF)*sizeof(double));
    memcpy(NxN_Matrix, p_Filter->NxN_Matrix,SQR_FILTER*SQR_FILTER*(MAX_NUM_AIF)*sizeof(double));
    memcpy(FilterFlag,p_Filter->FilterFlag,sizeof(int)*MAX_NUM_AIF);
    memcpy(FilterCoef, p_Filter->FilterCoef,SQR_FILTER*(MAX_NUM_AIF)*sizeof(double));
    memcpy(FilterCoef16bits, p_Filter->FilterCoef16bits,SQR_FILTER*(MAX_NUM_AIF)*sizeof(short int));
    memcpy(nBitsIntRepresent,p_Filter->nBitsIntRepresent,sizeof(int)*MAX_NUM_AIF);
  }
}
void initEDAIF2(void)
{
  int i,k,j;
  RealTimeDecomposition = 0;
  memset(realSymCommands,0,sizeof(int)*4*MAX_NUM_AIF);

  for(k = 0; k < MAX_NUM_AIF+MAX_NUM_AIF_EXTRA; k++)
  {
    for(i = 0; i < SQR_FILTER; i++)
    {
      for(j = 0; j < SQR_FILTER; j++)
      {
        NxN_Matrix[k][i][j] = 0.0;
      }
      N_Vector[k][i]       = 0.0;
    }
    FILTER_NEW_SUB_POS[k] = k;

    if (input->ImpType==0)
      nBitsIntRepresent[k]=DEFAULT_QUANT;
    else
      nBitsIntRepresent[k]=nBitsIntRepresent_DIAF[k];
  }


  {
    memcpy(TwoDEquationPattern,TwoDEquationPattern_RAIF,sizeof(int)*MAX_NUM_AIF*SQR_FILTER);
    memcpy(Calc2Filt_Indexes,Calc2Filt_Indexes_RAIF,sizeof(int)*MAX_NUM_AIF*SQR_FILTER);
    memcpy(STANDARD_2D_FILTER,STANDARD_2D_FILTER_RAIF,sizeof(double)*MAX_NUM_AIF*SQR_FILTER);

    memcpy(IsDiagonal1D,IsDiagonal1D_RAIF,sizeof(int)*MAX_NUM_AIF);
    memcpy(Is1DPosition,Is1DPosition_RAIF,sizeof(int)*MAX_NUM_AIF);
    memcpy(SymmetryPosition,SymmetryPosition_RAIF,sizeof(int)*MAX_NUM_AIF);
    memcpy(POS_EQUATION_NUMBER,POS_EQUATION_NUMBER_RAIF,sizeof(int)*MAX_NUM_AIF);
    memcpy(FILTER_NEW_SUB_POS,FILTER_NEW_SUB_POS_RAIF,sizeof(int)*MAX_NUM_AIF);
    memcpy(Filt_Indexes_offset,Filt_Indexes_offset_RAIF,sizeof(int)*MAX_NUM_AIF);
    memcpy(&UseEquation[0][0],&UseEquation_EAIF[0][0],sizeof(int)*3*(SQR_FILTER-1));
    SwapAIF(1, 1);// copy global settings to the storage
  }// Init RAIF filters set

  {
    memcpy(TwoDEquationPattern,TwoDEquationPattern_DAIF,sizeof(int)*MAX_NUM_AIF*SQR_FILTER);
    memcpy(Calc2Filt_Indexes,Calc2Filt_Indexes_DAIF,sizeof(int)*MAX_NUM_AIF*SQR_FILTER);
    memcpy(STANDARD_2D_FILTER,STANDARD_2D_FILTER_DAIF,sizeof(double)*MAX_NUM_AIF*SQR_FILTER);

    memcpy(IsDiagonal1D,IsDiagonal1D_DAIF,sizeof(int)*MAX_NUM_AIF);
    memcpy(Is1DPosition,Is1DPosition_DAIF,sizeof(int)*MAX_NUM_AIF);
    memcpy(SymmetryPosition,SymmetryPosition_DAIF,sizeof(int)*MAX_NUM_AIF);
    memcpy(POS_EQUATION_NUMBER,POS_EQUATION_NUMBER_DAIF,sizeof(int)*MAX_NUM_AIF);
    memcpy(FILTER_NEW_SUB_POS,FILTER_NEW_SUB_POS_DAIF,sizeof(int)*MAX_NUM_AIF);
    memcpy(Filt_Indexes_offset,Filt_Indexes_offset_DAIF,sizeof(int)*MAX_NUM_AIF);
    memcpy(&UseEquation[0][0],&UseEquation_DAIF[0][0],sizeof(int)*3*(SQR_FILTER-1));
    SwapAIF(0, 1);// copy global settings to the storage
  }// Init DAIF filters set

  // copy DAIF into global vars
  SwapAIF(0, 0);
  //  add-on RAIF settings
  {
    memcpy(TwoDEquationPattern[f2_pos], TwoDEquationPattern_RAIF[f_pos],SQR_FILTER*sizeof(int));
    memcpy(TwoDEquationPattern[i2_pos], TwoDEquationPattern_RAIF[i_pos],3*SQR_FILTER*sizeof(int));
    memcpy(TwoDEquationPattern[n2_pos], TwoDEquationPattern_RAIF[n_pos],SQR_FILTER*sizeof(int));

    memcpy(Calc2Filt_Indexes[f2_pos], Calc2Filt_Indexes_RAIF[f_pos],SQR_FILTER*sizeof(int));
    memcpy(Calc2Filt_Indexes[i2_pos], Calc2Filt_Indexes_RAIF[i_pos],3*SQR_FILTER*sizeof(int));
    memcpy(Calc2Filt_Indexes[n2_pos], Calc2Filt_Indexes_RAIF[n_pos],SQR_FILTER*sizeof(int));

    memcpy(STANDARD_2D_FILTER[f2_pos], STANDARD_2D_FILTER_RAIF[f_pos],SQR_FILTER*sizeof(double));
    memcpy(STANDARD_2D_FILTER[i2_pos], STANDARD_2D_FILTER_RAIF[i_pos],3*SQR_FILTER*sizeof(double));
    memcpy(STANDARD_2D_FILTER[n2_pos], STANDARD_2D_FILTER_RAIF[n_pos],SQR_FILTER*sizeof(double));

    for(k = MAX_NUM_AIF; k < MAX_NUM_AIF+MAX_NUM_AIF_EXTRA; k++)
    {
      Filt_Indexes_offset[k] = 1;
      SymmetryPosition[k] = 1;
      Is1DPosition[k] = 0;
      IsDiagonal1D[k] = 0;
      POS_EQUATION_NUMBER[k] = 12;
      if (input->ImpType==0)
        nBitsIntRepresent[k]=DEFAULT_QUANT;
      else
        nBitsIntRepresent[k]=nBitsIntRepresent_RAIF[k];
    }
  }
}



void initFilterCustom(int filterID)
{
  UseAllSubpelPositions = 0;
  SubpelPositionsPattern = 0;

  memset(FilterFlag,0,MAX_NUM_SUBPELS_AIF*sizeof(int));
  memset(DiffQFilterCoef, 0, MAX_AIF_SUPPORT*MAX_NUM_SUBPELS_AIF*sizeof(int));
  memset(FilterCoef, 0, MAX_AIF_SUPPORT*MAX_NUM_SUBPELS_AIF*sizeof(double));

  if  ((filterID == FILTER_TYPE_2D_NS) || (filterID == FILTER_TYPE_2D_S))
  {
    memcpy(STANDARD_2D_FILTER,STANDARD_2D_FILTER_orig,sizeof(double)*MAX_NUM_SUBPELS_AIF*SQR_FILTER);
    memcpy(TwoDEquationPattern,TwoDEquationPattern_orig,sizeof(int)*MAX_NUM_SUBPELS_AIF*SQR_FILTER);
    memcpy(TwoDSymmetricPattern,TwoDSymmetricPattern_orig,sizeof(int)*MAX_NUM_SUBPELS_AIF*SQR_FILTER);

    memcpy(IsDiagonal1D,IsDiagonal1D_orig,sizeof(int)*15);
    memcpy(Is1DPosition,Is1DPosition_orig,sizeof(int)*15);
    memcpy(SymmetryPosition,SymmetryPosition_orig,sizeof(int)*15);
    memcpy(POS_EQUATION_NUMBER,POS_EQUATION_NUMBER_orig,sizeof(int)*15);
    memcpy(FILTER_NEW_SUB_POS,FILTER_NEW_SUB_POS_orig,sizeof(int)*MAX_NUM_SUBPELS_AIF);
  }
#ifdef E_DAIF 
  else if (filterID == FILTER_TYPE_1D || filterID == FILTER_TYPE_EDAIF)
#else
  else if (filterID == FILTER_TYPE_1D)
#endif
  {
    memcpy(&STANDARD_2D_FILTER[0][0],&STANDARD_2D_FILTER_v4[0][0],sizeof(double)*MAX_NUM_SUBPELS_AIF*SQR_FILTER);
    memcpy(&TwoDEquationPattern[0][0],&TwoDEquationPattern_v4[0][0],sizeof(int)*MAX_NUM_SUBPELS_AIF*SQR_FILTER);
    memcpy(&TwoDSymmetricPattern[0][0],&TwoDSymmetricPattern_v4[0][0],sizeof(int)*MAX_NUM_SUBPELS_AIF*SQR_FILTER);
    memcpy(&IsDiagonal1D[0],&IsDiagonal1D_v4[0],sizeof(int)*MAX_NUM_SUBPELS_AIF);
    memcpy(&Is1DPosition[0],&Is1DPosition_v4[0],sizeof(int)*MAX_NUM_SUBPELS_AIF);
    memcpy(&SymmetryPosition[0],&SymmetryPosition_v4[0],sizeof(int)*MAX_NUM_SUBPELS_AIF);
    memcpy(&POS_EQUATION_NUMBER[0],&POS_EQUATION_NUMBER_v4[0],sizeof(int)*MAX_NUM_SUBPELS_AIF);
    memcpy(&FILTER_NEW_SUB_POS[0],&FILTER_NEW_SUB_POS_v4[0],sizeof(int)*MAX_NUM_SUBPELS_AIF);
    memcpy(&UseEquation[0][0],&UseEquation_DAIF[0][0],sizeof(int)*3*(SQR_FILTER-1));
  }
#ifdef EAIF
  else if (filterID == FILTER_TYPE_EAIF)
  {
    memcpy(&STANDARD_2D_FILTER[0][0],&STANDARD_2D_FILTER_orig[0][0],sizeof(double)*MAX_NUM_SUBPELS_AIF*SQR_FILTER);
    memcpy(&TwoDEquationPattern[0][0],&TwoDEquationPattern_EAIF[0][0],sizeof(int)*MAX_NUM_SUBPELS_AIF*SQR_FILTER);
    memcpy(&TwoDSymmetricPattern[0][0],&TwoDSymmetricPattern_v4[0][0],sizeof(int)*MAX_NUM_SUBPELS_AIF*SQR_FILTER);
    memcpy(&IsDiagonal1D[0],&IsDiagonal1D_v4[0],sizeof(int)*MAX_NUM_SUBPELS_AIF);
    memcpy(&Is1DPosition[0],&Is1DPosition_orig[0],sizeof(int)*MAX_NUM_SUBPELS_AIF);
    memcpy(&SymmetryPosition[0],&SymmetryPosition_EAIF[0],sizeof(int)*MAX_NUM_SUBPELS_AIF);
    memcpy(&POS_EQUATION_NUMBER[0],&POS_EQUATION_NUMBER_EAIF[0],sizeof(int)*MAX_NUM_SUBPELS_AIF);
    memcpy(&FILTER_NEW_SUB_POS[0],&FILTER_NEW_SUB_POS_EAIF[0],sizeof(int)*MAX_NUM_SUBPELS_AIF);
    memcpy(&UseEquation[0][0],&UseEquation_EAIF[0][0],sizeof(int)*3*(SQR_FILTER-1));
  }
#endif
}

#endif
/*!
************************************************************************
* \brief
*    init adaptive filter structures
************************************************************************
*/
void InitAdaptiveFilter(InputParameters *input)
{
  int block_x_nr;
  int block_y_nr;
  if(!input->UseAdaptiveFilter)
  {
    return;
  }
#ifdef DIRECTIONAL_FILTER
  if (input->UseAdaptiveFilter!= FILTER_TYPE_1D)
    input->ImpType = IMP_FLOAT32;
  if (input->ImpType == IMP_INT16)
    FILTCOEF_BITS = 7;
  else 
    FILTCOEF_BITS = DEFAULT_QUANT;
#endif
  NumberOfQBits = DEFAULT_QUANT;
  if(NumberOfQBits < 5 || NumberOfQBits > 16)
  {
    NumberOfQBits=DEFAULT_QUANT;
    printf("setting NumberOfQBits to default value...%d\n", NumberOfQBits);
  }
#ifndef EDAIF2
  FILTER_OFFSET = (FILTER_SIZE)/2-1;
#endif
#ifdef E_DAIF
  FILTER_OFFSET_INT = (FILTER_SIZE_INT)/2;
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

  number_of_macroblocks = img->width*img->height/256;
  BitsPerFrameCounter = calloc(number_of_macroblocks, sizeof(int));
  get_mem2Dint(&MBRateStorage, number_of_macroblocks, MAXMODE);
  block_y_nr = img->height/4;
  block_x_nr = img->width/4;
  BitsPerFrameCounter = calloc(number_of_macroblocks, sizeof(int));
  get_mem2Dint(&MBRateStorage, number_of_macroblocks, MAXMODE);
  current_poc=-1;

}

/*!
************************************************************************
* \brief
*    free adaptive filter structure
************************************************************************
*/
void FreeAdaptiveFilter()
{
  int block_x_nr;
  int block_y_nr;
  if(!input->UseAdaptiveFilter)
    return;
  if(MBRateStorage)
    free_mem2Dint(MBRateStorage);
  if(BitsPerFrameCounter)
    free(BitsPerFrameCounter);
  block_y_nr = img->height/4;
  block_x_nr = img->width/4;
}


/*!
************************************************************************
* \brief
*    resets all filter data
************************************************************************
*/
void ResetAdaptiveFilter()
{
  int i, j, k;
  if(!UseAdaptiveFilterForCurrentFrame())
  {
    return;
  }
  else
  {
    for(k = 0; k < 15; k++)
    {
      for(i = 0; i < SQR_FILTER; i++)
      {
        for(j = 0; j < SQR_FILTER; j++)
        {
          NxN_Matrix[k][i][j] = 0.0;
        }
        N_Vector[k][i]       = 0.0;
        FilterCoef[k][i]     = 0.0;
        QFilterCoef[k][i]    = 0;
        PredFilterCoef[k][i] = 0.0;
      }
      FilterFlag[k] = 0;
    }
#ifdef E_DAIF
    for(i = 0; i < SQR_FILTER_INT; ++i)
    {
      for(j = 0; j < SQR_FILTER_INT; ++j)
      {
        NxN_Matrix_Int[i][j] = 0.0;
      }
      N_Vector_Int[i] = 0.0;
      FilterCoefInt[i] = 0.0;
    }
    FilterFlagInt = 0;
#endif 
  }
  NumberOfQBits = DEFAULT_QUANT;
}


/*!
************************************************************************
* \brief
*    return if adaptive filter is used for current frame
************************************************************************
*/
int UseAdaptiveFilterForCurrentFrame(void)
{
  if((img->type == P_SLICE || img->type == B_SLICE) && input->UseAdaptiveFilter)
    return 1;
  else
    return 0;
}


/*!
************************************************************************
* \brief
*    sets necessary filter data for current macroblock
************************************************************************
*/
void SetAdaptiveFilter(void)
{
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int i, j, x, y, xj, yi;
  float b_scaling=1.0;
  int run=0, rerun=0;
  int equation=0;
  int mvx_subF=0, mvy_subF=0;
  int mvx_subB=0, mvy_subB=0;
  int mvx=0,  mvy=0;
  int mvxF=0, mvyF=0;
  int mvxB=0, mvyB=0;
  int sub_pos=0;
  int sub_posF=0;
  int sub_posB=0;
  int ref_frame=0;
  int ref_frameF=-1;
  int ref_frameB=-1;
  int subblock=0;
  int filter_pos_x1, filter_pos_y1, filter_pos_x2, filter_pos_y2;
  int x_orig, y_orig;       // absolute distances from current macroblock
  int h1, h2, h3;
  int *new_filter_pos, *new_eq, *new_sub_pos;
  int number_equation;
  imgpel **RecY=NULL;
  imgpel **RecYF=NULL;
  imgpel **RecYB=NULL;
  new_filter_pos = &h1;
  new_eq = &h2;
  new_sub_pos = &h3;
#ifdef EIGHTH_PEL
  if(img->mv_res)
  {
    SetAdaptiveFilterEighthpel();
    return;
  }
#endif
  if(IS_INTRA(currMB)) //intra macroblocks are not used for calculation of the filter coeffs.
    return;

  x_orig = MB_BLOCK_SIZE*(img->block_x/4);
  y_orig = MB_BLOCK_SIZE*(img->block_y/4);
  for(subblock = 0; subblock < 16; subblock++)
  {
    x = x_orig+4*(subblock%4);
    y = y_orig+4*(subblock/4);
    mvxF = enc_picture->mv[LIST_0][y/4][x/4][0];
    mvyF = enc_picture->mv[LIST_0][y/4][x/4][1];
    if(mvxF >= 0)
      mvx_subF = mvxF%4;      // x-sub-coordinate in a 4x4block
    else
      mvx_subF = (4-abs(mvxF)%4)%4;
    if(mvyF >= 0)
      mvy_subF = mvyF%4;      // y-sub-coordinate in a 4x4block
    else
      mvy_subF = (4-abs(mvyF)%4)%4;

    sub_posF = mvx_subF + 4*mvy_subF-1;    // pos 0..14 in a 4x4 block
    ref_frameF = enc_picture->ref_idx[LIST_0][y/4][x/4];
    if(ref_frameF != -1)
      RecYF  = listX[LIST_0][ref_frameF]->imgY;
    if(img->type == B_SLICE)
    {
      mvxB= enc_picture->mv[LIST_1][y/4][x/4][0];
      mvyB= enc_picture->mv[LIST_1][y/4][x/4][1];
      if(mvxB >= 0)
        mvx_subB = mvxB%4;      // x-sub-coordinate in a 4x4block
      else
        mvx_subB = (4-abs(mvxB)%4)%4;
      if(mvyB >= 0)
        mvy_subB = mvyB%4;      // y-sub-coordinate in a 4x4block
      else
        mvy_subB = (4-abs(mvyB)%4)%4;

      sub_posB = mvx_subB + 4*mvy_subB-1;    // pos 0..14 in a 4x4 block
      ref_frameB = enc_picture->ref_idx[LIST_1][y/4][x/4];
      if(ref_frameB != -1)
        RecYB  = listX[LIST_1][ref_frameB]->imgY;
      rerun = 1;
    }
    else
      rerun = 0;
    if(ref_frameB != -1 && ref_frameF != -1)
      b_scaling = 0.5;
    else
      b_scaling = 1.0;
    for(run = 0; run <= rerun; run++)
    {
      if(run == 0)
      {
        sub_pos = sub_posF;
        mvx = mvxF;
        mvy = mvyF;
        RecY = RecYF;
        ref_frame = ref_frameF;
      }
      else
      {
        sub_pos = sub_posB;
        mvx = mvxB;
        mvy = mvyB;
        RecY = RecYB;
        ref_frame = ref_frameB;
      }

      if(sub_pos != -1 && ref_frame != -1 ) // if no full-pel motion
      {
#ifndef DIRECTIONAL_FILTER
        if((sub_pos+1)%4 == 0 || (sub_pos+1)/4 == 0)      // if line or column, then only 1D-Filter,
          number_equation = FILTER_SIZE;                  // hor or vert. resp.
        else
          number_equation = SQR(FILTER_SIZE);
#else
        if (Is1DPosition_orig[sub_pos])// separate branch for Z18 filter
          number_equation = FILTER_SIZE;          // 1D case
        else
          number_equation = SQR(FILTER_SIZE);
#endif

        for(yi = 0; yi < 4; yi++)    //y
        {
          for(xj = 0; xj < 4; xj++)  //x
          {
            for(equation = 0; equation < number_equation; equation++)
            {
              if(number_equation == SQR(FILTER_SIZE))
              {
                filter_pos_x1 = FindPosition(img->width, x+xj,equation%FILTER_SIZE,mvx);          //current derivation (h_ij)
                filter_pos_y1 = FindPosition(img->height,y+yi,equation/FILTER_SIZE,mvy);          //current derivation (h_ij)
#ifdef DIRECTIONAL_FILTER
                if (IsDiagonal1D[sub_pos]==0)
#endif
                {
                  for(i = 0; i < FILTER_SIZE; i++)
                  {
                    for(j = 0; j < FILTER_SIZE; j++)
                    {
                      filter_pos_x2 = FindPosition(img->width, x+xj,j,mvx);
                      filter_pos_y2 = FindPosition(img->height,y+yi,i,mvy);
                      ReorderPosition(new_filter_pos,  new_eq,   new_sub_pos,
                        FILTER_SIZE*i+j, equation, sub_pos);
                      NxN_Matrix[*new_sub_pos][*new_eq][*new_filter_pos] += b_scaling*
                        RecY[filter_pos_y1][filter_pos_x1]*
                        RecY[filter_pos_y2][filter_pos_x2];
                    }
                  }
                }
#ifdef EAIF
                else if(input->UseAdaptiveFilter == FILTER_TYPE_EAIF)
                {
                  if(!UseEquation[0][equation]) 
                    continue; 
                  for(i = 0; i < FILTER_SIZE; i++)
                  {
                    for(j = 0; j < FILTER_SIZE; j++)
                    {
                      if(!UseEquation[0][FILTER_SIZE*i+j])
                        continue; 
                      filter_pos_x2 = FindPosition(img->width, x+xj,j,mvx);
                      filter_pos_y2 = FindPosition(img->height,y+yi,i,mvy);
                      ReorderPosition(new_filter_pos,  new_eq,   new_sub_pos,
                        FILTER_SIZE*i+j, equation, sub_pos);
                      NxN_Matrix[*new_sub_pos][*new_eq][*new_filter_pos] += b_scaling *
                        RecY[filter_pos_y1][filter_pos_x1] *
                        RecY[filter_pos_y2][filter_pos_x2];
                    }
                  }
                }
#endif
#ifdef DIRECTIONAL_FILTER
                else if (IsDiagonal1D[sub_pos]==1)//NW-SE
                {
                  // Count correlation functions only for diagonal NW-SE
                  if (!((equation==0)||(equation==7)||(equation==14)||
                    (equation==21)||(equation==28)||(equation==35)))
                    continue;
                  for(i = 0; i < FILTER_SIZE; i++)
                  {
                    j = i;
                    filter_pos_x2 = FindPosition(img->width, x+xj,j,mvx);
                    filter_pos_y2 = FindPosition(img->height,y+yi,i,mvy);
                    ReorderPosition(new_filter_pos,  new_eq,   new_sub_pos,
                      FILTER_SIZE*i+j, equation, sub_pos);
                    NxN_Matrix[*new_sub_pos][*new_eq][*new_filter_pos] += b_scaling*
                      RecY[filter_pos_y1][filter_pos_x1]*
                      RecY[filter_pos_y2][filter_pos_x2];
                  }
                }else if (IsDiagonal1D[sub_pos]==2)//NE-SW
                {
                  // Count correlation functions only for diagonal NE-SW
                  if (!((equation==5)||(equation==10)||(equation==15)||
                    (equation==20)||(equation==25)||(equation==30)))

                    continue;
                  for(i = FILTER_SIZE-1; i >=0; i--)
                  {
                    j = FILTER_SIZE-1 - i;
                    filter_pos_x2 = FindPosition(img->width, x+xj,j,mvx);
                    filter_pos_y2 = FindPosition(img->height,y+yi,i,mvy);
                    ReorderPosition(new_filter_pos,  new_eq,   new_sub_pos,
                      FILTER_SIZE*i+j, equation, sub_pos);
                    NxN_Matrix[*new_sub_pos][*new_eq][*new_filter_pos] += b_scaling*
                      RecY[filter_pos_y1][filter_pos_x1]*
                      RecY[filter_pos_y2][filter_pos_x2];
                  }
                }
                else if (IsDiagonal1D[sub_pos]==3)//NW-SE & NE-SW
                {
                  // Count correlation functions only for diagonal NW-SE
                  if (!((equation==0)||(equation==7)||(equation==14)||
                    (equation== 5)||(equation==10)||(equation==15)||
                    (equation==20)||(equation==25)||(equation==30)||
                    (equation==21)||(equation==28)||(equation==35)))
                    continue;
                  for(i = 0; i < FILTER_SIZE; i++)
                  {
                    j = i;
                    filter_pos_x2 = FindPosition(img->width, x+xj,j,mvx);
                    filter_pos_y2 = FindPosition(img->height,y+yi,i,mvy);
                    ReorderPosition(new_filter_pos,  new_eq,   new_sub_pos,
                      FILTER_SIZE*i+j, equation, sub_pos);
                    NxN_Matrix[*new_sub_pos][*new_eq][*new_filter_pos] += b_scaling*
                      RecY[filter_pos_y1][filter_pos_x1]*
                      RecY[filter_pos_y2][filter_pos_x2];
                    {//  Compute XCors for diagonal filter
                      j = FILTER_SIZE-1 - i;
                      filter_pos_x2 = FindPosition(img->width, x+xj,j,mvx);
                      filter_pos_y2 = FindPosition(img->height,y+yi,i,mvy);
                      ReorderPosition(new_filter_pos,  new_eq,   new_sub_pos,
                        FILTER_SIZE*i+j, equation, sub_pos);
                      NxN_Matrix[*new_sub_pos][*new_eq][*new_filter_pos] += b_scaling*
                        RecY[filter_pos_y1][filter_pos_x1]*
                        RecY[filter_pos_y2][filter_pos_x2];
                    }
                  }
                }
#endif
              }
              else   // 1D-Filter case
              {
                if((sub_pos+1) / 4 == 0)  // horizontal
                {
                  filter_pos_x1 = FindPosition(img->width, x+xj,equation%FILTER_SIZE,mvx);                                               //current derivation (h_ij)
                  filter_pos_y1 = FindPosition(img->height,y+yi,FILTER_OFFSET,mvy);                                               //current derivation (h_ij)
                  for(i = 0; i < FILTER_SIZE; i++)
                  {
                    filter_pos_x2 = FindPosition(img->width, x+xj, i,mvx);
                    filter_pos_y2 = FindPosition(img->height,y+yi, FILTER_OFFSET,mvy);
                    ReorderPosition(new_filter_pos,   new_eq, new_sub_pos,
                      i, equation,   sub_pos);

                    NxN_Matrix[*new_sub_pos][*new_eq][*new_filter_pos] += b_scaling*
                      RecY[filter_pos_y1][filter_pos_x1]*
                      RecY[filter_pos_y2][filter_pos_x2];
                  }
                }
                else // vertical
                {
                  filter_pos_x1 = FindPosition(img->width, x+xj,FILTER_OFFSET,mvx);                                                  //current derivation (h_ij)
                  filter_pos_y1 = FindPosition(img->height,y+yi,equation%FILTER_SIZE,mvy);                                              //current derivation (h_ij)
                  for(i = 0; i < FILTER_SIZE; i++)
                  {
                    filter_pos_x2 = FindPosition(img->width, x+xj, FILTER_OFFSET,mvx);
                    filter_pos_y2 = FindPosition(img->height,y+yi, i,mvy);
                    ReorderPosition(new_filter_pos,   new_eq, new_sub_pos,
                      i, equation,   sub_pos);
                    NxN_Matrix[*new_sub_pos][*new_eq][*new_filter_pos] += b_scaling*
                      RecY[filter_pos_y1][filter_pos_x1]*
                      RecY[filter_pos_y2][filter_pos_x2];
                  }
                }
              }

#ifdef EAIF
              if(input->UseAdaptiveFilter == FILTER_TYPE_EAIF && !Is1DPosition[sub_pos])
              {
                if(!UseEquation[0][equation])
                  continue; 
              }
              else
#endif
#ifdef DIRECTIONAL_FILTER
                if (IsDiagonal1D[sub_pos]==1)//NW-SE
                {
                  // Count correlation functions only for diagonal NW-SE
                  if (!((equation==0)||(equation==7)||(equation==14)||
                    (equation==21)||(equation==28)||(equation==35)))
                    continue;
                }
                else if (IsDiagonal1D[sub_pos]==2)//NE-SW
                {
                  // Count correlation functions only for diagonal NW-SE
                  if (!((equation==5)||(equation==10)||(equation==15)||
                    (equation==20)||(equation==25)||(equation==30)))
                    continue;
                }
                else if (IsDiagonal1D[sub_pos]==3)//NW-SE & NE-SW
                {
                  if (!((equation==0)||(equation==7)||(equation==14)||
                    (equation== 5)||(equation==10)||(equation==15)||
                    (equation==20)||(equation==25)||(equation==30)||
                    (equation==21)||(equation==28)||(equation==35)))
                    continue;
                }
#endif
                ReorderPosition(new_filter_pos,   new_eq, new_sub_pos,
                  1, equation,   sub_pos);       // 1 is not important
                N_Vector[*new_sub_pos][*new_eq] += b_scaling*imgY_org[y+yi][x+xj]*
                  RecY[filter_pos_y1][filter_pos_x1];
            } // equation loop
#ifdef E_DAIF
#ifdef EAIF
            if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF || input->UseAdaptiveFilter == FILTER_TYPE_EAIF)
#else
            if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)
#endif
            {
              // one more coefficient to solve to get filter offset
              if(number_equation == SQR(FILTER_SIZE))   // 2-D case
              {
                int N = POS_EQUATION_NUMBER[sub_pos];

                {
                  for(i = 0; i < FILTER_SIZE; i++)
                    for(j = 0; j < FILTER_SIZE; j++)
                    {
                      if(!UseEquation[IsDiagonal1D[sub_pos]-1][FILTER_SIZE*i+j])
                        continue; 
                      filter_pos_x2 = FindPosition(img->width, x+xj,j,mvx);
                      filter_pos_y2 = FindPosition(img->height,y+yi,i,mvy);
                      *new_sub_pos = FILTER_NEW_SUB_POS[sub_pos];
                      *new_eq = N;
                      *new_filter_pos = TwoDEquationPattern[sub_pos][FILTER_SIZE*i+j];
                      NxN_Matrix[*new_sub_pos][*new_eq][*new_filter_pos] += b_scaling*
                        RecY[filter_pos_y2][filter_pos_x2];
                      NxN_Matrix[*new_sub_pos][*new_filter_pos][*new_eq] += b_scaling*
                        RecY[filter_pos_y2][filter_pos_x2];
                    }
                }
                *new_sub_pos = FILTER_NEW_SUB_POS[sub_pos];
                NxN_Matrix[*new_sub_pos][N][N] += 1.; 
                N_Vector[*new_sub_pos][N] += b_scaling*imgY_org[y+yi][x+xj];
              }
              else    // 1-D case
              {
                int N = POS_EQUATION_NUMBER[sub_pos];
                if((sub_pos+1) / 4 == 0)	// horizontal
                {
                  for(i = 0; i < FILTER_SIZE; i++)
                  {
                    filter_pos_x2 = FindPosition(img->width, x+xj, i,mvx);
                    filter_pos_y2 = FindPosition(img->height,y+yi, FILTER_OFFSET,mvy);
                    *new_sub_pos = FILTER_NEW_SUB_POS[sub_pos];
                    *new_eq = N;
                    *new_filter_pos = TwoDEquationPattern[sub_pos][i];
                    NxN_Matrix[*new_sub_pos][*new_eq][*new_filter_pos] += b_scaling*
                      RecY[filter_pos_y2][filter_pos_x2];
                    NxN_Matrix[*new_sub_pos][*new_filter_pos][*new_eq] += b_scaling*
                      RecY[filter_pos_y2][filter_pos_x2];
                  }
                }
                else // vertical
                {
                  for(i = 0; i < FILTER_SIZE; i++)
                  {
                    filter_pos_x2 = FindPosition(img->width, x+xj, FILTER_OFFSET,mvx);
                    filter_pos_y2 = FindPosition(img->height,y+yi, i,mvy);
                    *new_sub_pos = FILTER_NEW_SUB_POS[sub_pos];
                    *new_eq = N;
                    *new_filter_pos = TwoDEquationPattern[sub_pos][i];
                    NxN_Matrix[*new_sub_pos][*new_eq][*new_filter_pos] += b_scaling*
                      RecY[filter_pos_y2][filter_pos_x2];
                    NxN_Matrix[*new_sub_pos][*new_filter_pos][*new_eq] += b_scaling*
                      RecY[filter_pos_y2][filter_pos_x2];
                  }
                }
                *new_sub_pos = FILTER_NEW_SUB_POS[sub_pos];
                NxN_Matrix[*new_sub_pos][N][N] += 1.; 
                N_Vector[*new_sub_pos][N] += b_scaling*imgY_org[y+yi][x+xj];
              }
            }
#endif
          }
        }
      }
#ifdef E_DAIF
#ifdef EAIF
      else if((input->UseAdaptiveFilter == FILTER_TYPE_EDAIF || input->UseAdaptiveFilter == FILTER_TYPE_EAIF)&&(ref_frame != -1))
#else
      else if((input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)&&(ref_frame != -1))
#endif
        // full-pel position 
      {
        number_equation = SQR(FILTER_SIZE_INT);
        for(yi = 0; yi < 4; yi++)    //y
        {
          for(xj = 0; xj < 4; xj++)  //x
          {
            for(equation = 0; equation < number_equation; equation++)
            {
              filter_pos_x1 = FindPositionInt(img->width, x+xj,equation%FILTER_SIZE_INT,mvx);					//current derivation (h_ij)
              filter_pos_y1 = FindPositionInt(img->height,y+yi,equation/FILTER_SIZE_INT,mvy);					//current derivation (h_ij)
              for(i = 0; i < FILTER_SIZE_INT; i++)
              {
                for(j = 0; j < FILTER_SIZE_INT; j++)
                {
                  filter_pos_x2 = FindPositionInt(img->width, x+xj,j,mvx);
                  filter_pos_y2 = FindPositionInt(img->height,y+yi,i,mvy);
                  *new_eq         = equation; 
                  *new_filter_pos = FILTER_SIZE_INT*i+j;
                  NxN_Matrix_Int[*new_eq][*new_filter_pos] += b_scaling *
                    RecY[filter_pos_y1][filter_pos_x1] *
                    RecY[filter_pos_y2][filter_pos_x2];
                }
              }
              *new_eq         = equation; 
              N_Vector_Int[*new_eq] += b_scaling * imgY_org[y+yi][x+xj] *
                RecY[filter_pos_y1][filter_pos_x1];
            }   // equation loop

            // one more coefficient for offset
            // equation = SQR(FILTER_SIZE_INT)
            {
              int N; 
              for(i = 0; i < FILTER_SIZE_INT; i++)
              {
                for(j = 0; j < FILTER_SIZE_INT; j++)
                {
                  filter_pos_x2 = FindPositionInt(img->width, x+xj,j,mvx);
                  filter_pos_y2 = FindPositionInt(img->height,y+yi,i,mvy);

                  *new_eq         = equation; 
                  *new_filter_pos = FILTER_SIZE_INT*i+j;
                  NxN_Matrix_Int[*new_eq][*new_filter_pos] += b_scaling *
                    RecY[filter_pos_y2][filter_pos_x2];
                  NxN_Matrix_Int[*new_filter_pos][*new_eq] += b_scaling *
                    RecY[filter_pos_y2][filter_pos_x2];
                }
              }

              N = SQR(FILTER_SIZE_INT);
              NxN_Matrix_Int[N][N] += 1.; 
              N_Vector_Int[N] += b_scaling * imgY_org[y+yi][x+xj];
            }
          }
        }
      }
#endif  // E_DAIF 
#ifdef EDAIF2
      if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF){
        switch (sub_pos)
        {
        case f_pos:
          sub_pos = f2_pos;
          break;
        case i_pos:
          sub_pos = i2_pos;
          break;
        case j_pos:
          sub_pos = j2_pos;
          break;
        case k_pos:
          sub_pos = k2_pos;
          break;
        case n_pos:
          sub_pos = n2_pos;
          break;
        default:
          continue;
          break;
        }
        {
          number_equation = SQR(FILTER_SIZE);
          for(yi = 0; yi < 4; yi++)    //y
          {
            for(xj = 0; xj < 4; xj++)  //x
            {
              for(equation = 0; equation < number_equation; equation++)
              {
                if (TwoDEquationPattern[sub_pos][equation]==-1)
                  continue;
                filter_pos_x1 = FindPosition(img->width, x+xj,equation%FILTER_SIZE,mvx);          //current derivation (h_ij)
                filter_pos_y1 = FindPosition(img->height,y+yi,equation/FILTER_SIZE,mvy);          //current derivation (h_ij)

                for(i = 0; i < FILTER_SIZE; i++)
                {
                  for(j = 0; j < FILTER_SIZE; j++)
                  {
                    if (TwoDEquationPattern[sub_pos][i*FILTER_SIZE+j]==-1)
                      continue;
                    filter_pos_x2 = FindPosition(img->width, x+xj,j,mvx);
                    filter_pos_y2 = FindPosition(img->height,y+yi,i,mvy);
                    ReorderPosition(new_filter_pos,  new_eq,   new_sub_pos,
                      FILTER_SIZE*i+j, equation, sub_pos);
                    NxN_Matrix[*new_sub_pos][*new_eq][*new_filter_pos] += b_scaling*
                      RecY[filter_pos_y1][filter_pos_x1]*
                      RecY[filter_pos_y2][filter_pos_x2];
                  }
                }
                ReorderPosition(new_filter_pos,   new_eq, new_sub_pos,
                  1, equation,   sub_pos);       // 1 is not important
                N_Vector[*new_sub_pos][*new_eq] += b_scaling*imgY_org[y+yi][x+xj]*
                  RecY[filter_pos_y1][filter_pos_x1];
              }
#ifdef E_DAIF
#ifdef EAIF
              if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF || input->UseAdaptiveFilter == FILTER_TYPE_EAIF)
#else
              if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)
#endif
              {
                // one more coefficient to solve to get filter offset
                if(number_equation == SQR(FILTER_SIZE))   // 2-D case
                {
                  int N = POS_EQUATION_NUMBER[sub_pos];

                  for(i = 0; i < FILTER_SIZE; i++)
                  {
                    for(j = 0; j < FILTER_SIZE; j++)
                    {
                      if (TwoDEquationPattern[sub_pos][equation]==-1)
                        continue;
                      filter_pos_x2 = FindPosition(img->width, x+xj,j,mvx);
                      filter_pos_y2 = FindPosition(img->height,y+yi,i,mvy);
                      *new_sub_pos = FILTER_NEW_SUB_POS[sub_pos];
                      *new_eq = N;
                      *new_filter_pos = TwoDEquationPattern[sub_pos][FILTER_SIZE*i+j];
                      NxN_Matrix[*new_sub_pos][*new_eq][*new_filter_pos] += b_scaling*
                        RecY[filter_pos_y2][filter_pos_x2];
                      NxN_Matrix[*new_sub_pos][*new_filter_pos][*new_eq] += b_scaling*
                        RecY[filter_pos_y2][filter_pos_x2];
                    }
                  }
                  *new_sub_pos = FILTER_NEW_SUB_POS[sub_pos];
                  NxN_Matrix[*new_sub_pos][N][N] += 1.; 
                  N_Vector[*new_sub_pos][N] += b_scaling*imgY_org[y+yi][x+xj];
                }
              }
#endif
            }
          }
        }
      }
#endif
    }
  }
}

// separable aif (BEGIN)
/*!
************************************************************************
* \brief
*    sets necessary horizontal filter data for current macroblock
************************************************************************
*/
void SetAdaptiveFilterHor(void)
{
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int i, x, y, xj, yi;
  float b_scaling=1.0;
  int run=0, rerun=0;
  int equation=0, eq_new = 0, i_new = 0;
  int mvx_subF=0, mvy_subF=0;
  int mvx_subB=0, mvy_subB=0;
  int mvx=0,  mvy=0;
  int mvxF=0, mvyF=0;
  int mvxB=0, mvyB=0;
  int sub_pos=0;
  int sub_posF=0;
  int sub_posB=0;
  int ref_frame=0;
  int ref_frameF=-1;
  int ref_frameB=-1;
  int subblock=0;
  int filter_pos_x1, filter_pos_y1, filter_pos_x2;
  int x_orig, y_orig;       // absolute distances from current macroblock
  int number_equation;
  imgpel **RecY=NULL;
  imgpel **RecYF=NULL;
  imgpel **RecYB=NULL;
#ifdef EIGHTH_PEL
  if(img->mv_res)
  {
    SetAdaptiveFilterEighthpelHor();
    return;
  }
#endif
  if(IS_INTRA(currMB)) //intra macroblocks are not used for calculation of the filter coeffs.
    return;

  x_orig = MB_BLOCK_SIZE*(img->block_x/4);
  y_orig = MB_BLOCK_SIZE*(img->block_y/4);
  for(subblock = 0; subblock < 16; subblock++)
  {
    x = x_orig+4*(subblock%4);
    y = y_orig+4*(subblock/4);
    mvxF = enc_picture->mv[LIST_0][y/4][x/4][0];
    mvyF = enc_picture->mv[LIST_0][y/4][x/4][1];
    if(mvxF >= 0)
      mvx_subF = mvxF%4;      // x-sub-coordinate in a 4x4block
    else
      mvx_subF = (4-abs(mvxF)%4)%4;
    if(mvyF >= 0)
      mvy_subF = mvyF%4;      // y-sub-coordinate in a 4x4block
    else
      mvy_subF = (4-abs(mvyF)%4)%4;

    sub_posF = mvx_subF + 4*mvy_subF-1;    // pos 0..14 in a 4x4 block
    ref_frameF = enc_picture->ref_idx[LIST_0][y/4][x/4];
    if(ref_frameF != -1)
      RecYF  = listX[LIST_0][ref_frameF]->imgY;
    if(img->type == B_SLICE)
    {
      mvxB= enc_picture->mv[LIST_1][y/4][x/4][0];
      mvyB= enc_picture->mv[LIST_1][y/4][x/4][1];
      if(mvxB >= 0)
        mvx_subB = mvxB%4;      // x-sub-coordinate in a 4x4block
      else
        mvx_subB = (4-abs(mvxB)%4)%4;
      if(mvyB >= 0)
        mvy_subB = mvyB%4;      // y-sub-coordinate in a 4x4block
      else
        mvy_subB = (4-abs(mvyB)%4)%4;

      sub_posB = mvx_subB + 4*mvy_subB-1;    // pos 0..14 in a 4x4 block
      ref_frameB = enc_picture->ref_idx[LIST_1][y/4][x/4];
      if(ref_frameB != -1)
        RecYB  = listX[LIST_1][ref_frameB]->imgY;
      rerun = 1;
    }
    else
      rerun = 0;
    if(ref_frameB != -1 && ref_frameF != -1)
      b_scaling = 0.5;
    else
      b_scaling = 1.0;
    for(run = 0; run <= rerun; run++)
    {
      if(run == 0)
      {
        sub_pos = sub_posF;
        mvx = mvxF;
        mvy = mvyF;
        RecY = RecYF;
        ref_frame = ref_frameF;
      }
      else
      {
        sub_pos = sub_posB;
        mvx = mvxB;
        mvy = mvyB;
        RecY = RecYB;
        ref_frame = ref_frameB;
      }

      if(sub_pos != -1 && ref_frame != -1 && sub_pos < d_pos) 
      {
        number_equation = FILTER_SIZE;    
        for(yi = 0; yi < 4; yi++)    //y
        {
          filter_pos_y1 = FindPosition(img->height,y+yi,FILTER_OFFSET,mvy); 
          for(xj = 0; xj < 4; xj++)  //x
          {
            for(equation = 0; equation < number_equation; equation++)
            {
              if(sub_pos == a_pos || sub_pos == c_pos) 
              {
                filter_pos_x1 = FindPosition(img->width, x+xj,equation,mvx); 
                for(i = 0; i < FILTER_SIZE; i++)
                {
                  filter_pos_x2 = FindPosition(img->width, x+xj, i,mvx);
                  NxN_Matrix[sub_pos][equation][i] += b_scaling*
                    RecY[filter_pos_y1][filter_pos_x1]*
                    RecY[filter_pos_y1][filter_pos_x2];
                }
                N_Vector[sub_pos][equation] += b_scaling*imgY_org[y+yi][x+xj]*
                  RecY[filter_pos_y1][filter_pos_x1];
              }
              else // b_pos 
              {
                filter_pos_x1 = FindPosition(img->width, x+xj,equation,mvx);           
                eq_new = TwoDEquationPattern[b_pos][equation];
                for(i = 0; i < FILTER_SIZE; i++)
                {
                  filter_pos_x2 = FindPosition(img->width, x+xj, i,mvx);
                  i_new = TwoDEquationPattern[b_pos][i];
                  NxN_Matrix[sub_pos][eq_new][i_new] += b_scaling*
                    RecY[filter_pos_y1][filter_pos_x1]*
                    RecY[filter_pos_y1][filter_pos_x2];
                }
                N_Vector[sub_pos][eq_new] += b_scaling*imgY_org[y+yi][x+xj]*
                  RecY[filter_pos_y1][filter_pos_x1];
              }
            }
          }
        }
      }
    }
  }
}

/*!
************************************************************************
* \brief
*    sets necessary vertical filter data for current macroblock
************************************************************************
*/
void SetAdaptiveFilterVer(void)
{
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int i, x, y, xj, yi;
  float b_scaling=1.0;
  int run=0, rerun=0;
  int equation=0, eq_new = 0, i_new = 0;
  int mvx_subF=0, mvy_subF=0;
  int mvx_subB=0, mvy_subB=0;
  int mvx=0,  mvy=0;
  int mvxF=0, mvyF=0;
  int mvxB=0, mvyB=0;
  int sub_pos=0;
  int sub_posF=0;
  int sub_posB=0;
  int ref_frame=0;
  int ref_frameF=-1;
  int ref_frameB=-1;
  int subblock=0;
  int filter_pos_x1, filter_pos_y1, filter_pos_y2;
  int x_orig, y_orig;       // absolute distances from current macroblock
  int number_equation;
  double **RecY=NULL; 
  double **RecYF=NULL; 
  double **RecYB=NULL; 
  imgpel **RecY_FullPel = NULL;
  imgpel **RecYF_FullPel = NULL;
  imgpel **RecYB_FullPel = NULL;

#ifdef EIGHTH_PEL
  if(img->mv_res)
  {
    SetAdaptiveFilterEighthpelVer();
    return;
  }
#endif
  if(IS_INTRA(currMB)) //intra macroblocks are not used for calculation of the filter coeffs.
    return;

  x_orig = MB_BLOCK_SIZE*(img->block_x/4);
  y_orig = MB_BLOCK_SIZE*(img->block_y/4);
  for(subblock = 0; subblock < 16; subblock++)
  {
    x = x_orig+4*(subblock%4);
    y = y_orig+4*(subblock/4);
    mvxF = enc_picture->mv[LIST_0][y/4][x/4][0];
    mvyF = enc_picture->mv[LIST_0][y/4][x/4][1];
    if(mvxF >= 0)
      mvx_subF = mvxF%4;      // x-sub-coordinate in a 4x4block
    else
      mvx_subF = (4-abs(mvxF)%4)%4;
    if(mvyF >= 0)
      mvy_subF = mvyF%4;      // y-sub-coordinate in a 4x4block
    else
      mvy_subF = (4-abs(mvyF)%4)%4;

    sub_posF = mvx_subF + 4*mvy_subF-1;    // pos 0..14 in a 4x4 block
    ref_frameF = enc_picture->ref_idx[LIST_0][y/4][x/4];
    if(ref_frameF != -1)
    {
      RecYF = listX[LIST_0][ref_frameF]->imgY_ups_aif_hor;
      RecYF_FullPel = listX[LIST_0][ref_frameF]->imgY;
    }
    if(img->type == B_SLICE)
    {
      mvxB= enc_picture->mv[LIST_1][y/4][x/4][0];
      mvyB= enc_picture->mv[LIST_1][y/4][x/4][1];
      if(mvxB >= 0)
        mvx_subB = mvxB%4;      // x-sub-coordinate in a 4x4block
      else
        mvx_subB = (4-abs(mvxB)%4)%4;
      if(mvyB >= 0)
        mvy_subB = mvyB%4;      // y-sub-coordinate in a 4x4block
      else
        mvy_subB = (4-abs(mvyB)%4)%4;

      sub_posB = mvx_subB + 4*mvy_subB-1;    // pos 0..14 in a 4x4 block
      ref_frameB = enc_picture->ref_idx[LIST_1][y/4][x/4];
      if(ref_frameB != -1)
      {
        RecYB = listX[LIST_1][ref_frameB]->imgY_ups_aif_hor;
        RecYB_FullPel = listX[LIST_1][ref_frameB]->imgY;
      }
      rerun = 1;
    }
    else
      rerun = 0;
    if(ref_frameB != -1 && ref_frameF != -1)
      b_scaling = 0.5;
    else
      b_scaling = 1.0;
    for(run = 0; run <= rerun; run++)
    {
      if(run == 0)
      {
        sub_pos = sub_posF;
        mvx = mvxF;
        mvy = mvyF;
        RecY = RecYF;
        RecY_FullPel = RecYF_FullPel;
        ref_frame = ref_frameF;
      }
      else
      {
        sub_pos = sub_posB;
        mvx = mvxB;
        mvy = mvyB;
        RecY = RecYB;
        RecY_FullPel = RecYB_FullPel;
        ref_frame = ref_frameB;
      }

      if(sub_pos != -1 && ref_frame != -1 ) // if no full-pel motion
      {
        number_equation = FILTER_SIZE;      
        if (sub_pos == d_pos)
        {
          for (xj = 0; xj < 4; xj++)     
          {
            filter_pos_x1 = FindPosition (img->width, x + xj, FILTER_OFFSET, mvx); 
            for (yi = 0; yi < 4; yi++)      
            {
              for (equation = 0; equation < number_equation; equation++)
              {
                filter_pos_y1 = FindPosition (img->height, y + yi, equation, mvy);       
                for (i = 0; i < FILTER_SIZE; i++)
                {
                  filter_pos_y2 = FindPosition (img->height, y + yi, i, mvy);
                  NxN_Matrix[sub_pos][equation][i] +=
                    b_scaling * RecY_FullPel[filter_pos_y1][filter_pos_x1] * RecY_FullPel[filter_pos_y2][filter_pos_x1];
                }
                N_Vector[sub_pos][equation] += b_scaling * imgY_org[y + yi][x + xj] * RecY_FullPel[filter_pos_y1][filter_pos_x1];
              }
            }
          }
        }
        else if (sub_pos == h_pos)
        {
          for (xj = 0; xj < 4; xj++)   
          {
            filter_pos_x1 = FindPosition (img->width, x + xj, FILTER_OFFSET, mvx); 
            for (yi = 0; yi < 4; yi++) 
            {
              for (equation = 0; equation < number_equation; equation++)
              {
                filter_pos_y1 = FindPosition (img->height, y + yi, equation, mvy);        
                eq_new = TwoDEquationPattern[b_pos][equation]; 
                for (i = 0; i < FILTER_SIZE; i++)
                {
                  filter_pos_y2 = FindPosition (img->height, y + yi, i, mvy);
                  i_new = TwoDEquationPattern[b_pos][i]; 
                  NxN_Matrix[sub_pos][eq_new][i_new] +=
                    b_scaling * RecY_FullPel[filter_pos_y1][filter_pos_x1] * RecY_FullPel[filter_pos_y2][filter_pos_x1];
                }
                N_Vector[sub_pos][eq_new] += b_scaling * imgY_org[y + yi][x + xj] * RecY_FullPel[filter_pos_y1][filter_pos_x1];
              }
            }
          }
        }
        else if (sub_pos == l_pos)
        {
          for (xj = 0; xj < 4; xj++)  
          {
            filter_pos_x1 = FindPosition (img->width, x + xj, FILTER_OFFSET, mvx); 
            for (yi = 0; yi < 4; yi++)     
            {
              for (equation = 0; equation < number_equation; equation++)
              {
                filter_pos_y1 = FindPosition (img->height, y + yi, equation, mvy);        
                for (i = 0; i < FILTER_SIZE; i++)
                {
                  filter_pos_y2 = FindPosition (img->height, y + yi, i, mvy);
                  NxN_Matrix[d_pos][number_equation-equation-1][FILTER_SIZE-i-1] +=
                    b_scaling * RecY_FullPel[filter_pos_y1][filter_pos_x1] * RecY_FullPel[filter_pos_y2][filter_pos_x1];
                }
                N_Vector[d_pos][number_equation-equation-1] += b_scaling * imgY_org[y + yi][x + xj] * RecY_FullPel[filter_pos_y1][filter_pos_x1];
              }
            }
          }
        }
        else if (sub_pos > l_pos)   
        {
          for (xj = 0; xj < 4; xj++)     
          {
            filter_pos_x1 = FindPosition (img->width, x + xj, FILTER_OFFSET, mvx);    
            filter_pos_x1 = (filter_pos_x1 << 2) + ((sub_pos + 1) % 4); 
            for (yi = 0; yi < 4; yi++)      
            {
              for (equation = 0; equation < number_equation; equation++)
              {
                filter_pos_y1 = FindPosition (img->height, y + yi, equation, mvy);        
                for (i = 0; i < FILTER_SIZE; i++)
                {
                  filter_pos_y2 = FindPosition (img->height, y + yi, i, mvy);
                  NxN_Matrix[sub_pos-8][number_equation-equation-1][FILTER_SIZE-i-1] += 
                    b_scaling * RecY[filter_pos_y1 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4] * 
                    RecY[filter_pos_y2 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4];
                }
                N_Vector[sub_pos-8][number_equation-equation-1] +=
                  b_scaling * imgY_org[y + yi][x + xj] * RecY[filter_pos_y1 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4];
              }
            }
          }
        }
        else if (sub_pos == i_pos || sub_pos == j_pos || sub_pos == k_pos)
        {
          for (xj = 0; xj < 4; xj++)     
          {
            filter_pos_x1 = FindPosition (img->width, x + xj, FILTER_OFFSET, mvx);   
            filter_pos_x1 = (filter_pos_x1 << 2) + ((sub_pos + 1) % 4); 
            for (yi = 0; yi < 4; yi++)       
            {
              for (equation = 0; equation < number_equation; equation++)
              {
                filter_pos_y1 = FindPosition (img->height, y + yi, equation, mvy);        
                eq_new = TwoDEquationPattern[b_pos][equation];  
                for (i = 0; i < FILTER_SIZE; i++)
                {
                  filter_pos_y2 = FindPosition (img->height, y + yi, i, mvy);
                  i_new = TwoDEquationPattern[b_pos][i]; 
                  NxN_Matrix[sub_pos][eq_new][i_new] +=
                    b_scaling * RecY[filter_pos_y1 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4] * 
                    RecY[filter_pos_y2 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4];
                }
                N_Vector[sub_pos][eq_new] +=
                  b_scaling * imgY_org[y + yi][x + xj] * RecY[filter_pos_y1 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4];
              }
            }
          }
        }
        else 
        {
          for (xj = 0; xj < 4; xj++)   
          {
            filter_pos_x1 = FindPosition (img->width, x + xj, FILTER_OFFSET, mvx);  
            filter_pos_x1 = (filter_pos_x1 << 2) + ((sub_pos + 1) % 4);  
            for (yi = 0; yi < 4; yi++)      
            {
              for (equation = 0; equation < number_equation; equation++)
              {
                filter_pos_y1 = FindPosition (img->height, y + yi, equation, mvy);      
                for (i = 0; i < FILTER_SIZE; i++)
                {
                  filter_pos_y2 = FindPosition (img->height, y + yi, i, mvy);
                  NxN_Matrix[sub_pos][equation][i] +=
                    b_scaling * RecY[filter_pos_y1 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4] * 
                    RecY[filter_pos_y2 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4];
                }
                N_Vector[sub_pos][equation] +=
                  b_scaling * imgY_org[y + yi][x + xj] * RecY[filter_pos_y1 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4];
              }
            }
          }
        }
      }
    }
  }
}
// separable aif (END)
/*!
************************************************************************
* \brief
*    sets necessary filter data for current macroblock in 1/8-pel case
************************************************************************
*/

void SetAdaptiveFilterEighthpel(void)
{
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int i, j, x, y, xj, yi;
  float b_scaling=1.0;
  int run=0, rerun=0;
  int run8=0, rerun8 = 1;
  int equation=0;
  int mvx_subF=0, mvy_subF=0;
  int mvx_subB=0, mvy_subB=0;
  int mvx=0,  mvy=0;
  int mvxF=0, mvyF=0;
  int mvxB=0, mvyB=0;
  int mvxF8, mvyF8, mvxB8, mvyB8;
  int tmp_x_f, tmp_y_f, tmp_x_b, tmp_y_b;
  int sub_pos=0;
  int sub_posF=0;
  int sub_posB=0;
  int ref_frame=0;
  int ref_frameF=-1;
  int ref_frameB=-1;
  int subblock=0;
  int filter_pos_x1, filter_pos_y1, filter_pos_x2, filter_pos_y2;
  int x_orig, y_orig;       // absolute distances from current macroblock
  int h1, h2, h3;
  int *new_filter_pos, *new_eq, *new_sub_pos;
  int number_equation;
  imgpel **RecY=NULL;
  imgpel **RecYF=NULL;
  imgpel **RecYB=NULL;
  new_filter_pos = &h1;
  new_eq = &h2;
  new_sub_pos = &h3;

  if(IS_INTRA(currMB)) //intra macroblocks are not used for calculation of the filter coeffs.
    return;

  x_orig = MB_BLOCK_SIZE*(img->block_x/4);
  y_orig = MB_BLOCK_SIZE*(img->block_y/4);
  for(run8=0;run8<rerun8;run8++)
  {
    for(subblock = 0; subblock < 16; subblock++)
    {
      x = x_orig+4*(subblock%4);
      y = y_orig+4*(subblock/4);
      mvxF8 = enc_picture->mv[LIST_0][y/4][x/4][0];
      mvyF8 = enc_picture->mv[LIST_0][y/4][x/4][1];
      if(mvxF8>=0)
        tmp_x_f = mvxF8%8;
      else
        tmp_x_f = 8-(abs(mvxF8)%8);
      if(mvyF8>=0)
        tmp_y_f = mvyF8%8;
      else
        tmp_y_f = 8-(abs(mvyF8)%8);

      if(run8 == 0)
      {
        if(  abs(mvxF8)%2 == 0  && abs(mvyF8)%2 == 0 )    // quarter-pel
        {
          mvxF = mvxF8/2;
          mvyF = mvyF8/2;
        }
        else if(abs(mvxF8)%2 == 0)    // quarter-pel in x
        {
          mvxF = mvxF8/2;
          mvyF = (mvyF8-1)/2;
        }
        else if(abs(mvyF8)%2 == 0)    // quarter-pel in y
        {
          mvxF = (mvxF8-1)/2;
          mvyF = mvyF8/2;
        }
        else    // upper right - down left
          // or upper left - down right interpolation
        {
          if((tmp_x_f+tmp_y_f)%4 == 0)    // upper right-down left
          {
            mvxF = (mvxF8+1)/2;
            mvyF = (mvyF8-1)/2;
          }
          else    // upper left - down right
          {
            mvxF = (mvxF8-1)/2;
            mvyF = (mvyF8-1)/2;
          }
        }
      }
      else // if(run8 == 0)
      {
        if(  abs(mvxF8)%2 == 0  && abs(mvyF8)%2 == 0 )    // quarter-pel
        {
          mvxF = mvxF8/2;
          mvyF = mvyF8/2;
        }
        else if(abs(mvxF8)%2 == 0)    // quarter-pel in x
        {
          mvxF = mvxF8/2;
          mvyF = (mvyF8+1)/2;
        }
        else if(abs(mvyF8)%2 == 0)    // quarter-pel in y
        {
          mvxF = (mvxF8+1)/2;
          mvyF = mvyF8/2;
        }
        else    // upper right - down left
          // or upper left - down right interpolation
        {
          if((tmp_x_f+tmp_y_f)%4 == 0)        // upper right - down left interpolation
          {
            mvxF = (mvxF8-1)/2;
            mvyF = (mvyF8+1)/2;
          }
          else            // upper left - down right interpolation
          {
            mvxF = (mvxF8+1)/2;
            mvyF = (mvyF8+1)/2;
          }
        }
      }

      if(mvxF >= 0)
        mvx_subF = mvxF%4;      // x-sub-coordinate in a 4x4block
      else
        mvx_subF = (4-abs(mvxF)%4)%4;
      if(mvyF >= 0)
        mvy_subF = mvyF%4;      // y-sub-coordinate in a 4x4block
      else
        mvy_subF = (4-abs(mvyF)%4)%4;

      sub_posF = mvx_subF + 4*mvy_subF-1;    // pos 0..14 in a 4x4 block
      ref_frameF = enc_picture->ref_idx[LIST_0][y/4][x/4];
      if(ref_frameF != -1)
        RecYF  = listX[LIST_0][ref_frameF]->imgY;
      if(img->type == B_SLICE)
      {
        mvxB8= enc_picture->mv[LIST_1][y/4][x/4][0];
        mvyB8= enc_picture->mv[LIST_1][y/4][x/4][1];
        if(mvxB8>=0)
          tmp_x_b = mvxB8%8;
        else
          tmp_x_b = 8-(abs(mvxB8)%8);
        if(mvyB8>=0)
          tmp_y_b = mvyB8%8;
        else
          tmp_y_b = 8-(abs(mvyB8)%8);
        if(run8 == 0)
        {
          if(  abs(mvxF8)%2 == 0  && abs(mvyF8)%2 == 0 )    // quarter-pel
          {
            mvxB = mvxB8/2;
            mvyB = mvyB8/2;
          }
          else if(abs(mvxB8)%2 == 0)    // quarter-pel in x
          {
            mvxB = mvxB8/2;
            mvyB = (mvyB8-1)/2;
          }
          else if(abs(mvyB8)%2 == 0)    // quarter-pel in y
          {
            mvxB = (mvxB8-1)/2;
            mvyB = mvyB8/2;
          }
          else    // upper right - down left
            // or upper left - down right interpolation
          {
            if((tmp_x_b+tmp_y_b)%4 == 0)    // upper right-down left
            {
              mvxB = (mvxB8+1)/2;
              mvyB = (mvyB8-1)/2;
            }
            else    // upper left - down right
            {
              mvxB = (mvxB8-1)/2;
              mvyB = (mvyB8-1)/2;
            }
          }
        }
        else
        {
          if(  abs(mvxB8)%2 == 0  && abs(mvyB8)%2 == 0 )    // quarter-pel
          {
            mvxB = mvxB8/2;
            mvyB = mvyB8/2;
          }
          else if(abs(mvxB8)%2 == 0)    // quarter-pel in x
          {
            mvxB = mvxB8/2;
            mvyB = (mvyB8+1)/2;
          }
          else if(abs(mvyB8)%2 == 0)    // quarter-pel in y
          {
            mvxB = (mvxB8+1)/2;
            mvyB = mvyB8/2;
          }
          else    // upper right - down left
            // or upper left - down right interpolation
          {
            if((tmp_x_b+tmp_y_b)%4 == 0)        // upper right - down left interpolation
            {
              mvxB = (mvxB8-1)/2;
              mvyB = (mvyB8+1)/2;
            }
            else            // upper left - down right interpolation
            {
              mvxB = (mvxB8+1)/2;
              mvyB = (mvyB8+1)/2;
            }
          }
        }

        if(mvxB >= 0)
          mvx_subB = mvxB%4;      // x-sub-coordinate in a 4x4block
        else
          mvx_subB = (4-abs(mvxB)%4)%4;
        if(mvyB >= 0)
          mvy_subB = mvyB%4;      // y-sub-coordinate in a 4x4block
        else
          mvy_subB = (4-abs(mvyB)%4)%4;

        sub_posB = mvx_subB + 4*mvy_subB-1;    // pos 0..14 in a 4x4 block
        ref_frameB = enc_picture->ref_idx[LIST_1][y/4][x/4];
        if(ref_frameB != -1)
          RecYB  = listX[LIST_1][ref_frameB]->imgY;
        rerun = 1;
      }
      else
        rerun = 0;
      if(ref_frameB != -1 && ref_frameF != -1)
        b_scaling = 0.5;
      else
        b_scaling = 1.0;
      for(run = 0; run <= rerun; run++)
      {
        if(run == 0)
        {
          sub_pos = sub_posF;
          mvx = mvxF;
          mvy = mvyF;
          RecY = RecYF;
          ref_frame = ref_frameF;
        }
        else
        {
          sub_pos = sub_posB;
          mvx = mvxB;
          mvy = mvyB;
          RecY = RecYB;
          ref_frame = ref_frameB;
        }

        if(sub_pos != -1 && ref_frame != -1 ) // if no full-pel motion
        {
          if((sub_pos+1)%4 == 0 || (sub_pos+1)/4 == 0)      // if line or column, then only 1D-Filter,
            number_equation = FILTER_SIZE;          // hor or vert. resp.
          else
            number_equation = SQR(FILTER_SIZE);

          for(yi = 0; yi < 4; yi++)    //y
          {
            for(xj = 0; xj < 4; xj++)  //x
            {
              for(equation = 0; equation < number_equation; equation++)
              {
                if(number_equation == SQR(FILTER_SIZE))
                {
                  filter_pos_x1 = FindPosition(img->width, x+xj,equation%FILTER_SIZE,mvx);          //current derivation (h_ij)
                  filter_pos_y1 = FindPosition(img->height,y+yi,equation/FILTER_SIZE,mvy);          //current derivation (h_ij)
                  for(i = 0; i < FILTER_SIZE; i++)
                  {
                    for(j = 0; j < FILTER_SIZE; j++)
                    {
                      filter_pos_x2 = FindPosition(img->width, x+xj,j,mvx);
                      filter_pos_y2 = FindPosition(img->height,y+yi,i,mvy);
                      ReorderPosition(new_filter_pos,  new_eq,   new_sub_pos,
                        FILTER_SIZE*i+j, equation, sub_pos);
                      NxN_Matrix[*new_sub_pos][*new_eq][*new_filter_pos] += 0.5*b_scaling*
                        RecY[filter_pos_y1][filter_pos_x1]*
                        RecY[filter_pos_y2][filter_pos_x2];
                    }
                  }
                }
                else   // 1D-Filter case
                {
                  if((sub_pos+1) / 4 == 0)  // horizontal
                  {
                    filter_pos_x1 = FindPosition(img->width, x+xj,equation%FILTER_SIZE,mvx);                                               //current derivation (h_ij)
                    filter_pos_y1 = FindPosition(img->height,y+yi,FILTER_OFFSET,mvy);                                               //current derivation (h_ij)
                    for(i = 0; i < FILTER_SIZE; i++)
                    {
                      filter_pos_x2 = FindPosition(img->width, x+xj, i,mvx);
                      filter_pos_y2 = FindPosition(img->height,y+yi, FILTER_OFFSET,mvy);
                      ReorderPosition(new_filter_pos,   new_eq, new_sub_pos,
                        i, equation,   sub_pos);

                      NxN_Matrix[*new_sub_pos][*new_eq][*new_filter_pos]+= 0.5*b_scaling*RecY[filter_pos_y1][filter_pos_x1]*RecY[filter_pos_y2][filter_pos_x2];
                    }
                  }
                  else // vertical
                  {
                    filter_pos_x1 = FindPosition(img->width, x+xj,FILTER_OFFSET,mvx);                                                  //current derivation (h_ij)
                    filter_pos_y1 = FindPosition(img->height,y+yi,equation%FILTER_SIZE,mvy);                                              //current derivation (h_ij)
                    for(i = 0; i < FILTER_SIZE; i++)
                    {
                      filter_pos_x2 = FindPosition(img->width, x+xj, FILTER_OFFSET,mvx);
                      filter_pos_y2 = FindPosition(img->height,y+yi, i,mvy);
                      ReorderPosition(new_filter_pos,   new_eq, new_sub_pos,
                        i, equation,   sub_pos);
                      NxN_Matrix[*new_sub_pos][*new_eq][*new_filter_pos] += 0.5*b_scaling*RecY[filter_pos_y1][filter_pos_x1]*RecY[filter_pos_y2][filter_pos_x2];
                    }
                  }
                }
                ReorderPosition(new_filter_pos,   new_eq, new_sub_pos,
                  1, equation,   sub_pos);       // 1 is not important
                N_Vector[*new_sub_pos][*new_eq] += 0.5*b_scaling*imgY_org[y+yi][x+xj]*
                  RecY[filter_pos_y1][filter_pos_x1];
              }
            }
          }
        }
      }
    }
  }
}

// separable aif (BEGIN)
#ifdef EIGHTH_PEL 
/*!
*****************************************************************************************
* \brief
*    sets necessary filter data for current macroblock in 1/8-pel case (separable filter)
*****************************************************************************************
*/

void SetAdaptiveFilterEighthpelHor(void)
{
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int i, x, y, xj, yi;
  float b_scaling=1.0;
  int run=0, rerun=0;
  int run8=0, rerun8 = 1;
  int equation=0, eq_new = 0, i_new = 0;
  int mvx_subF=0, mvy_subF=0;
  int mvx_subB=0, mvy_subB=0;
  int mvx=0,  mvy=0;
  int mvxF=0, mvyF=0;
  int mvxB=0, mvyB=0;
  int mvxF8, mvyF8, mvxB8, mvyB8;
  int tmp_x_f, tmp_y_f, tmp_x_b, tmp_y_b;
  int sub_pos=0;
  int sub_posF=0;
  int sub_posB=0;
  int ref_frame=0;
  int ref_frameF=-1;
  int ref_frameB=-1;
  int subblock=0;
  int filter_pos_x1, filter_pos_y1, filter_pos_x2;
  int x_orig, y_orig;       // absolute distances from current macroblock
  int h1, h2, h3;
  int *new_filter_pos, *new_eq, *new_sub_pos;
  int number_equation;
  imgpel **RecY=NULL;
  imgpel **RecYF=NULL;
  imgpel **RecYB=NULL;
  new_filter_pos = &h1;
  new_eq = &h2;
  new_sub_pos = &h3;
  if(IS_INTRA(currMB)) //intra macroblocks are not used for calculation of the filter coeffs.
    return;

  x_orig = MB_BLOCK_SIZE*(img->block_x/4);
  y_orig = MB_BLOCK_SIZE*(img->block_y/4);
  for(run8=0;run8<rerun8;run8++)
  {
    for(subblock = 0; subblock < 16; subblock++)
    {
      x = x_orig+4*(subblock%4);
      y = y_orig+4*(subblock/4);
      mvxF8 = enc_picture->mv[LIST_0][y/4][x/4][0];
      mvyF8 = enc_picture->mv[LIST_0][y/4][x/4][1];
      if(mvxF8>=0)
        tmp_x_f = mvxF8%8;
      else
        tmp_x_f = 8-(abs(mvxF8)%8);
      if(mvyF8>=0)
        tmp_y_f = mvyF8%8;
      else
        tmp_y_f = 8-(abs(mvyF8)%8);

      if(run8 == 0)
      {
        if(  abs(mvxF8)%2 == 0  && abs(mvyF8)%2 == 0 )    // quarter-pel
        {
          mvxF = mvxF8/2;
          mvyF = mvyF8/2;
        }
        else if(abs(mvxF8)%2 == 0)    // quarter-pel in x
        {
          mvxF = mvxF8/2;
          mvyF = (mvyF8-1)/2;
        }
        else if(abs(mvyF8)%2 == 0)    // quarter-pel in y
        {
          mvxF = (mvxF8-1)/2;
          mvyF = mvyF8/2;
        }
        else    // upper right - down left
          // or upper left - down right interpolation
        {
          if((tmp_x_f+tmp_y_f)%4 == 0)    // upper right-down left
          {
            mvxF = (mvxF8+1)/2;
            mvyF = (mvyF8-1)/2;
          }
          else    // upper left - down right
          {
            mvxF = (mvxF8-1)/2;
            mvyF = (mvyF8-1)/2;
          }
        }
      }
      else // if(run8 == 0)
      {
        if(  abs(mvxF8)%2 == 0  && abs(mvyF8)%2 == 0 )    // quarter-pel
        {
          mvxF = mvxF8/2;
          mvyF = mvyF8/2;
        }
        else if(abs(mvxF8)%2 == 0)    // quarter-pel in x
        {
          mvxF = mvxF8/2;
          mvyF = (mvyF8+1)/2;
        }
        else if(abs(mvyF8)%2 == 0)    // quarter-pel in y
        {
          mvxF = (mvxF8+1)/2;
          mvyF = mvyF8/2;
        }
        else    // upper right - down left
          // or upper left - down right interpolation
        {
          if((tmp_x_f+tmp_y_f)%4 == 0)        // upper right - down left interpolation
          {
            mvxF = (mvxF8-1)/2;
            mvyF = (mvyF8+1)/2;
          }
          else            // upper left - down right interpolation
          {
            mvxF = (mvxF8+1)/2;
            mvyF = (mvyF8+1)/2;
          }
        }
      }      
      if(mvxF >= 0)
        mvx_subF = mvxF%4;      // x-sub-coordinate in a 4x4block
      else
        mvx_subF = (4-abs(mvxF)%4)%4;
      if(mvyF >= 0)
        mvy_subF = mvyF%4;      // y-sub-coordinate in a 4x4block
      else
        mvy_subF = (4-abs(mvyF)%4)%4;

      sub_posF = mvx_subF + 4*mvy_subF-1;    // pos 0..14 in a 4x4 block
      ref_frameF = enc_picture->ref_idx[LIST_0][y/4][x/4];
      if(ref_frameF != -1)
        RecYF  = listX[LIST_0][ref_frameF]->imgY;
      if(img->type == B_SLICE)
      {
        mvxB8= enc_picture->mv[LIST_1][y/4][x/4][0];
        mvyB8= enc_picture->mv[LIST_1][y/4][x/4][1];
        if(mvxB8>=0)
          tmp_x_b = mvxB8%8;
        else
          tmp_x_b = 8-(abs(mvxB8)%8);
        if(mvyB8>=0)
          tmp_y_b = mvyB8%8;
        else
          tmp_y_b = 8-(abs(mvyB8)%8);
        if(run8 == 0)
        {
          if(  abs(mvxF8)%2 == 0  && abs(mvyF8)%2 == 0 )    // quarter-pel
          {
            mvxB = mvxB8/2;
            mvyB = mvyB8/2;
          }
          else if(abs(mvxB8)%2 == 0)    // quarter-pel in x
          {
            mvxB = mvxB8/2;
            mvyB = (mvyB8-1)/2;
          }
          else if(abs(mvyB8)%2 == 0)    // quarter-pel in y
          {
            mvxB = (mvxB8-1)/2;
            mvyB = mvyB8/2;
          }
          else    // upper right - down left
            // or upper left - down right interpolation
          {
            if((tmp_x_b+tmp_y_b)%4 == 0)    // upper right-down left
            {
              mvxB = (mvxB8+1)/2;
              mvyB = (mvyB8-1)/2;
            }
            else    // upper left - down right
            {
              mvxB = (mvxB8-1)/2;
              mvyB = (mvyB8-1)/2;
            }
          }
        }
        else
        {
          if(  abs(mvxB8)%2 == 0  && abs(mvyB8)%2 == 0 )    // quarter-pel
          {
            mvxB = mvxB8/2;
            mvyB = mvyB8/2;
          }
          else if(abs(mvxB8)%2 == 0)    // quarter-pel in x
          {
            mvxB = mvxB8/2;
            mvyB = (mvyB8+1)/2;
          }
          else if(abs(mvyB8)%2 == 0)    // quarter-pel in y
          {
            mvxB = (mvxB8+1)/2;
            mvyB = mvyB8/2;
          }
          else    // upper right - down left
            // or upper left - down right interpolation
          {
            if((tmp_x_b+tmp_y_b)%4 == 0)        // upper right - down left interpolation
            {
              mvxB = (mvxB8-1)/2;
              mvyB = (mvyB8+1)/2;
            }
            else            // upper left - down right interpolation
            {
              mvxB = (mvxB8+1)/2;
              mvyB = (mvyB8+1)/2;
            }
          }
        }
        if(mvxB >= 0)
          mvx_subB = mvxB%4;      // x-sub-coordinate in a 4x4block
        else
          mvx_subB = (4-abs(mvxB)%4)%4;
        if(mvyB >= 0)
          mvy_subB = mvyB%4;      // y-sub-coordinate in a 4x4block
        else
          mvy_subB = (4-abs(mvyB)%4)%4;

        sub_posB = mvx_subB + 4*mvy_subB-1;    // pos 0..14 in a 4x4 block
        ref_frameB = enc_picture->ref_idx[LIST_1][y/4][x/4];
        if(ref_frameB != -1)
          RecYB  = listX[LIST_1][ref_frameB]->imgY;
        rerun = 1;
      }
      else
        rerun = 0;
      if(ref_frameB != -1 && ref_frameF != -1)
        b_scaling = 0.5;
      else
        b_scaling = 1.0;
      for(run = 0; run <= rerun; run++)
      {
        if(run == 0)
        {
          sub_pos = sub_posF;
          mvx = mvxF;
          mvy = mvyF;
          RecY = RecYF;
          ref_frame = ref_frameF;
        }
        else
        {
          sub_pos = sub_posB;
          mvx = mvxB;
          mvy = mvyB;
          RecY = RecYB;
          ref_frame = ref_frameB;
        }

        if(sub_pos != -1 && ref_frame != -1 && sub_pos < d_pos)
        {
          number_equation = FILTER_SIZE;            
          for(yi = 0; yi < 4; yi++)    //y
          {
            for(xj = 0; xj < 4; xj++)  //x
            {
              for(equation = 0; equation < number_equation; equation++)
              {
                if(sub_pos == a_pos || sub_pos == c_pos)
                {
                  filter_pos_x1 = FindPosition(img->width, x+xj,equation,mvx);                                             
                  filter_pos_y1 = FindPosition(img->height,y+yi,FILTER_OFFSET,mvy);                                           
                  for(i = 0; i < FILTER_SIZE; i++)
                  {
                    filter_pos_x2 = FindPosition(img->width, x+xj, i,mvx);
                    NxN_Matrix[sub_pos][equation][i]+= 0.5*b_scaling*RecY[filter_pos_y1][filter_pos_x1]*RecY[filter_pos_y1][filter_pos_x2];
                  }
                  N_Vector[sub_pos][equation] += 0.5*b_scaling*imgY_org[y+yi][x+xj]*
                    RecY[filter_pos_y1][filter_pos_x1];
                }
                else // b_pos
                {
                  filter_pos_x1 = FindPosition(img->width, x+xj,equation,mvx);                                   
                  filter_pos_y1 = FindPosition(img->height,y+yi,FILTER_OFFSET,mvy);
                  eq_new = TwoDEquationPattern[b_pos][equation];
                  for(i = 0; i < FILTER_SIZE; i++)
                  {
                    filter_pos_x2 = FindPosition(img->width, x+xj, i,mvx);
                    i_new = TwoDEquationPattern[b_pos][i];  
                    NxN_Matrix[sub_pos][eq_new][i_new]+= 0.5*b_scaling*RecY[filter_pos_y1][filter_pos_x1]*RecY[filter_pos_y1][filter_pos_x2];
                  }
                  N_Vector[sub_pos][eq_new] += 0.5*b_scaling*imgY_org[y+yi][x+xj]*
                    RecY[filter_pos_y1][filter_pos_x1];
                }
              }
            }
          }
        }
      }
    }
  }
}

/*!
****************************************************************************************
* \brief
*    sets necessary filter data for current macroblock in 1/8-pel case (separable filter
****************************************************************************************
*/

void SetAdaptiveFilterEighthpelVer(void)
{
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int i, x, y, xj, yi;
  float b_scaling=1.0;
  int run=0, rerun=0;
  int run8=0, rerun8 = 1;
  int equation=0, eq_new = 0, i_new = 0;
  int mvx_subF=0, mvy_subF=0;
  int mvx_subB=0, mvy_subB=0;
  int mvx=0,  mvy=0;
  int mvxF=0, mvyF=0;
  int mvxB=0, mvyB=0;
  int mvxF8, mvyF8, mvxB8, mvyB8;
  int tmp_x_f, tmp_y_f, tmp_x_b, tmp_y_b;
  int sub_pos=0;
  int sub_posF=0;
  int sub_posB=0;
  int ref_frame=0;
  int ref_frameF=-1;
  int ref_frameB=-1;
  int subblock=0;
  int filter_pos_x1, filter_pos_y1, filter_pos_y2;
  int x_orig, y_orig;       // absolute distances from current macroblock
  int h1, h2, h3;
  int *new_filter_pos, *new_eq, *new_sub_pos;
  int number_equation;
  double **RecY=NULL; 
  double **RecYF=NULL; 
  double **RecYB=NULL; 
  imgpel **RecY_FullPel = NULL;
  imgpel **RecYF_FullPel = NULL;
  imgpel **RecYB_FullPel = NULL;
  new_filter_pos = &h1;
  new_eq = &h2;
  new_sub_pos = &h3;
  if(IS_INTRA(currMB)) //intra macroblocks are not used for calculation of the filter coeffs.
    return;

  x_orig = MB_BLOCK_SIZE*(img->block_x/4);
  y_orig = MB_BLOCK_SIZE*(img->block_y/4);
  for(run8=0;run8<rerun8;run8++)
  {
    for(subblock = 0; subblock < 16; subblock++)
    {
      x = x_orig+4*(subblock%4);
      y = y_orig+4*(subblock/4);
      mvxF8 = enc_picture->mv[LIST_0][y/4][x/4][0];
      mvyF8 = enc_picture->mv[LIST_0][y/4][x/4][1];
      if(mvxF8>=0)
        tmp_x_f = mvxF8%8;
      else
        tmp_x_f = 8-(abs(mvxF8)%8);
      if(mvyF8>=0)
        tmp_y_f = mvyF8%8;
      else
        tmp_y_f = 8-(abs(mvyF8)%8);

      if(run8 == 0)
      {
        if(  abs(mvxF8)%2 == 0  && abs(mvyF8)%2 == 0 )    // quarter-pel
        {
          mvxF = mvxF8/2;
          mvyF = mvyF8/2;
        }
        else if(abs(mvxF8)%2 == 0)    // quarter-pel in x
        {
          mvxF = mvxF8/2;
          mvyF = (mvyF8-1)/2;
        }
        else if(abs(mvyF8)%2 == 0)    // quarter-pel in y
        {
          mvxF = (mvxF8-1)/2;
          mvyF = mvyF8/2;
        }
        else    // upper right - down left
          // or upper left - down right interpolation
        {
          if((tmp_x_f+tmp_y_f)%4 == 0)    // upper right-down left
          {
            mvxF = (mvxF8+1)/2;
            mvyF = (mvyF8-1)/2;
          }
          else    // upper left - down right
          {
            mvxF = (mvxF8-1)/2;
            mvyF = (mvyF8-1)/2;
          }
        }
      }
      else // if(run8 == 0)
      {
        if(  abs(mvxF8)%2 == 0  && abs(mvyF8)%2 == 0 )    // quarter-pel
        {
          mvxF = mvxF8/2;
          mvyF = mvyF8/2;
        }
        else if(abs(mvxF8)%2 == 0)    // quarter-pel in x
        {
          mvxF = mvxF8/2;
          mvyF = (mvyF8+1)/2;
        }
        else if(abs(mvyF8)%2 == 0)    // quarter-pel in y
        {
          mvxF = (mvxF8+1)/2;
          mvyF = mvyF8/2;
        }
        else    // upper right - down left
          // or upper left - down right interpolation
        {
          if((tmp_x_f+tmp_y_f)%4 == 0)        // upper right - down left interpolation
          {
            mvxF = (mvxF8-1)/2;
            mvyF = (mvyF8+1)/2;
          }
          else            // upper left - down right interpolation
          {
            mvxF = (mvxF8+1)/2;
            mvyF = (mvyF8+1)/2;
          }
        }
      }      
      if(mvxF >= 0)
        mvx_subF = mvxF%4;      // x-sub-coordinate in a 4x4block
      else
        mvx_subF = (4-abs(mvxF)%4)%4;
      if(mvyF >= 0)
        mvy_subF = mvyF%4;      // y-sub-coordinate in a 4x4block
      else
        mvy_subF = (4-abs(mvyF)%4)%4;

      sub_posF = mvx_subF + 4*mvy_subF-1;    // pos 0..14 in a 4x4 block
      ref_frameF = enc_picture->ref_idx[LIST_0][y/4][x/4];
      if(ref_frameF != -1)
      {
        RecYF  = listX[LIST_0][ref_frameF]->imgY_ups_aif_hor;
        RecYF_FullPel = listX[LIST_0][ref_frameF]->imgY;
      }
      if(img->type == B_SLICE)
      {
        mvxB8= enc_picture->mv[LIST_1][y/4][x/4][0];
        mvyB8= enc_picture->mv[LIST_1][y/4][x/4][1];
        if(mvxB8>=0)
          tmp_x_b = mvxB8%8;
        else
          tmp_x_b = 8-(abs(mvxB8)%8);
        if(mvyB8>=0)
          tmp_y_b = mvyB8%8;
        else
          tmp_y_b = 8-(abs(mvyB8)%8);
        if(run8 == 0)
        {
          if(  abs(mvxF8)%2 == 0  && abs(mvyF8)%2 == 0 )    // quarter-pel
          {
            mvxB = mvxB8/2;
            mvyB = mvyB8/2;
          }
          else if(abs(mvxB8)%2 == 0)    // quarter-pel in x
          {
            mvxB = mvxB8/2;
            mvyB = (mvyB8-1)/2;
          }
          else if(abs(mvyB8)%2 == 0)    // quarter-pel in y
          {
            mvxB = (mvxB8-1)/2;
            mvyB = mvyB8/2;
          }
          else    // upper right - down left
            // or upper left - down right interpolation
          {
            if((tmp_x_b+tmp_y_b)%4 == 0)    // upper right-down left
            {
              mvxB = (mvxB8+1)/2;
              mvyB = (mvyB8-1)/2;
            }
            else    // upper left - down right
            {
              mvxB = (mvxB8-1)/2;
              mvyB = (mvyB8-1)/2;
            }
          }
        }
        else
        {
          if(  abs(mvxB8)%2 == 0  && abs(mvyB8)%2 == 0 )    // quarter-pel
          {
            mvxB = mvxB8/2;
            mvyB = mvyB8/2;
          }
          else if(abs(mvxB8)%2 == 0)    // quarter-pel in x
          {
            mvxB = mvxB8/2;
            mvyB = (mvyB8+1)/2;
          }
          else if(abs(mvyB8)%2 == 0)    // quarter-pel in y
          {
            mvxB = (mvxB8+1)/2;
            mvyB = mvyB8/2;
          }
          else    // upper right - down left
            // or upper left - down right interpolation
          {
            if((tmp_x_b+tmp_y_b)%4 == 0)        // upper right - down left interpolation
            {
              mvxB = (mvxB8-1)/2;
              mvyB = (mvyB8+1)/2;
            }
            else            // upper left - down right interpolation
            {
              mvxB = (mvxB8+1)/2;
              mvyB = (mvyB8+1)/2;
            }
          }
        }
        if(mvxB >= 0)
          mvx_subB = mvxB%4;      // x-sub-coordinate in a 4x4block
        else
          mvx_subB = (4-abs(mvxB)%4)%4;
        if(mvyB >= 0)
          mvy_subB = mvyB%4;      // y-sub-coordinate in a 4x4block
        else
          mvy_subB = (4-abs(mvyB)%4)%4;

        sub_posB = mvx_subB + 4*mvy_subB-1;    // pos 0..14 in a 4x4 block
        ref_frameB = enc_picture->ref_idx[LIST_1][y/4][x/4];
        if(ref_frameB != -1)
        {
          RecYB  = listX[LIST_1][ref_frameB]->imgY_ups_aif_hor;
          RecYB_FullPel = listX[LIST_1][ref_frameB]->imgY;
        }
        rerun = 1;
      }
      else
        rerun = 0;
      if(ref_frameB != -1 && ref_frameF != -1)
        b_scaling = 0.5;
      else
        b_scaling = 1.0;
      for(run = 0; run <= rerun; run++)
      {
        if(run == 0)
        {
          sub_pos = sub_posF;
          mvx = mvxF;
          mvy = mvyF;
          RecY = RecYF;
          RecY_FullPel = RecYF_FullPel;
          ref_frame = ref_frameF;
        }
        else
        {
          sub_pos = sub_posB;
          mvx = mvxB;
          mvy = mvyB;
          RecY = RecYB;
          RecY_FullPel = RecYB_FullPel;
          ref_frame = ref_frameB;
        }

        if(sub_pos != -1 && ref_frame != -1)
        {
          number_equation = FILTER_SIZE;      
          if (sub_pos == d_pos)
          {
            for (yi = 0; yi < 4; yi++)     
            {
              for (xj = 0; xj < 4; xj++)    
              {
                for (equation = 0; equation < number_equation; equation++)
                {
                  filter_pos_x1 = FindPosition (img->width, x + xj, FILTER_OFFSET, mvx);   
                  filter_pos_y1 = FindPosition (img->height, y + yi, equation, mvy);       
                  for (i = 0; i < FILTER_SIZE; i++)
                  {
                    filter_pos_y2 = FindPosition (img->height, y + yi, i, mvy);
                    NxN_Matrix[sub_pos][equation][i] +=
                      0.5*b_scaling * RecY_FullPel[filter_pos_y1][filter_pos_x1] * RecY_FullPel[filter_pos_y2][filter_pos_x1];
                  }
                  N_Vector[sub_pos][equation] += 0.5*b_scaling * imgY_org[y + yi][x + xj] * RecY_FullPel[filter_pos_y1][filter_pos_x1];
                }
              }
            }
          }
          else if (sub_pos == h_pos)
          {
            for (yi = 0; yi < 4; yi++)     
            {
              for (xj = 0; xj < 4; xj++)   
              {
                for (equation = 0; equation < number_equation; equation++)
                {
                  filter_pos_x1 = FindPosition (img->width, x + xj, FILTER_OFFSET, mvx);    
                  filter_pos_y1 = FindPosition (img->height, y + yi, equation, mvy);        
                  eq_new = TwoDEquationPattern[b_pos][equation]; 
                  for (i = 0; i < FILTER_SIZE; i++)
                  {
                    filter_pos_y2 = FindPosition (img->height, y + yi, i, mvy);
                    i_new = TwoDEquationPattern[b_pos][i]; 
                    NxN_Matrix[sub_pos][eq_new][i_new] +=
                      0.5*b_scaling * RecY_FullPel[filter_pos_y1][filter_pos_x1] * RecY_FullPel[filter_pos_y2][filter_pos_x1];
                  }
                  N_Vector[sub_pos][eq_new] += 0.5*b_scaling * imgY_org[y + yi][x + xj] * RecY_FullPel[filter_pos_y1][filter_pos_x1];
                }
              }
            }
          }
          else if (sub_pos == l_pos)
          {
            for (yi = 0; yi < 4; yi++)     
            {
              for (xj = 0; xj < 4; xj++)   
              {
                for (equation = 0; equation < number_equation; equation++)
                {
                  filter_pos_x1 = FindPosition (img->width, x + xj, FILTER_OFFSET, mvx);    
                  filter_pos_y1 = FindPosition (img->height, y + yi, equation, mvy);        
                  for (i = 0; i < FILTER_SIZE; i++)
                  {
                    filter_pos_y2 = FindPosition (img->height, y + yi, i, mvy);
                    NxN_Matrix[d_pos][number_equation-equation-1][FILTER_SIZE-i-1] +=
                      0.5*b_scaling * RecY_FullPel[filter_pos_y1][filter_pos_x1] * RecY_FullPel[filter_pos_y2][filter_pos_x1];
                  }
                  N_Vector[d_pos][number_equation-equation-1] += 0.5*b_scaling * imgY_org[y + yi][x + xj] * RecY_FullPel[filter_pos_y1][filter_pos_x1];
                }
              }
            }
          }
          else if (sub_pos > l_pos)   
          {
            for (yi = 0; yi < 4; yi++)      
            {
              for (xj = 0; xj < 4; xj++)    
              {
                for (equation = 0; equation < number_equation; equation++)
                {
                  filter_pos_x1 = FindPosition (img->width, x + xj, FILTER_OFFSET, mvx);    
                  filter_pos_x1 = (filter_pos_x1 << 2) + ((sub_pos + 1) % 4);
                  filter_pos_y1 = FindPosition (img->height, y + yi, equation, mvy);        
                  for (i = 0; i < FILTER_SIZE; i++)
                  {
                    filter_pos_y2 = FindPosition (img->height, y + yi, i, mvy);
                    NxN_Matrix[sub_pos-8][number_equation-equation-1][FILTER_SIZE-i-1] += 
                      0.5*b_scaling * RecY[filter_pos_y1 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4] * 
                      RecY[filter_pos_y2 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4];
                  }
                  N_Vector[sub_pos-8][number_equation-equation-1] +=
                    0.5*b_scaling * imgY_org[y + yi][x + xj] * RecY[filter_pos_y1 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4];
                }
              }
            }
          }
          else if (sub_pos == i_pos || sub_pos == j_pos || sub_pos == k_pos)
          {
            for (yi = 0; yi < 4; yi++)      
            {
              for (xj = 0; xj < 4; xj++)    
              {
                for (equation = 0; equation < number_equation; equation++)
                {
                  filter_pos_x1 = FindPosition (img->width, x + xj, FILTER_OFFSET, mvx);    
                  filter_pos_x1 = (filter_pos_x1 << 2) + ((sub_pos + 1) % 4);
                  filter_pos_y1 = FindPosition (img->height, y + yi, equation, mvy);        
                  eq_new = TwoDEquationPattern[b_pos][equation];  
                  for (i = 0; i < FILTER_SIZE; i++)
                  {
                    filter_pos_y2 = FindPosition (img->height, y + yi, i, mvy);
                    i_new = TwoDEquationPattern[b_pos][i]; 
                    NxN_Matrix[sub_pos][eq_new][i_new] +=
                      0.5*b_scaling * RecY[filter_pos_y1 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4] * 
                      RecY[filter_pos_y2 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4];
                  }
                  N_Vector[sub_pos][eq_new] +=
                    0.5*b_scaling * imgY_org[y + yi][x + xj] * RecY[filter_pos_y1 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4];
                }
              } 
            } 
          }
          else 
          {
            for (yi = 0; yi < 4; yi++)      
            {
              for (xj = 0; xj < 4; xj++)   
              {
                for (equation = 0; equation < number_equation; equation++)
                {
                  filter_pos_x1 = FindPosition (img->width, x + xj, FILTER_OFFSET, mvx);
                  filter_pos_x1 = (filter_pos_x1 << 2) + ((sub_pos + 1) % 4);
                  filter_pos_y1 = FindPosition (img->height, y + yi, equation, mvy);      
                  for (i = 0; i < FILTER_SIZE; i++)
                  { 
                    filter_pos_y2 = FindPosition (img->height, y + yi, i, mvy);
                    NxN_Matrix[sub_pos][equation][i] +=
                      0.5*b_scaling * RecY[filter_pos_y1 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4] * 
                      RecY[filter_pos_y2 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4];
                  }
                  N_Vector[sub_pos][equation] +=
                    0.5*b_scaling * imgY_org[y + yi][x + xj] * RecY[filter_pos_y1 + IMG_PAD_SIZE][filter_pos_x1 + IMG_PAD_SIZE * 4];
                }
              }
            }
          }
        }
      }
    }
  }
}
#endif 
// separable aif (END)

/*!
************************************************************************
* \brief
*    finds current position in X-direction
*    input:
*      DimX   - X-size of image,
*      curX   - current X-coordinate,
*      Offset - offset in X-direction from left filter position (2D-filter)
*      mvx    - x-motion vector
************************************************************************
*/
int FindPosition(int Dim, int cur_pos, int Offset, int mv)
{
  int position;
  if(mv < 0 && abs(mv)%4!=0)
    mv = mv-4;
  mv =mv/4;
  position = cur_pos-FILTER_OFFSET+Offset+mv;
  return max(0, min(Dim-1, position ) );
}

#ifdef E_DAIF
int FindPositionInt(int Dim, int cur_pos, int Offset, int mv)
{
  int position;
  mv = mv / 4;
  position = cur_pos - FILTER_OFFSET_INT + Offset + mv;
  return max(0, min(Dim - 1, position));
}
#endif  // E_DAIF 

/*!
************************************************************************
* \brief
*    reorder position, depending on current sub-position
*    in order to guarantee symmetric filter
*    input:
*      cur_xy   - current unsymmetric filter position
*      equation - current filter position
*      sub_pos  - sub-position, depending on mvx and mvy
*    output:
*      *new_filter_pos   - new symmetric filter position
*      *new_eq   - new equation
*      *new_sub_pos - new sub_position in a symmetric filter
************************************************************************
*/

void ReorderPosition(int *new_filter_pos, int *new_eq,  int *new_sub_pos,
                     int cur_xy,          int equation, int sub_pos)
{
  *new_sub_pos    = FILTER_NEW_SUB_POS[sub_pos];
  *new_eq         = TwoDEquationPattern[sub_pos][equation];
  *new_filter_pos = TwoDEquationPattern[sub_pos][cur_xy];
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
*                   using FILTER_NEW_SUB_POS
************************************************************************
*/
void ExtendFilterCoefficientsDouble(int sub_pos, double FCoef[15][SQR_FILTER])
{
  int i;
  if(SymmetryPosition[sub_pos] == 0)    //check if it is allowed to change original filter coefficients
  {
    for(i = 0; i < SQR_FILTER; i++)
      FCoef[sub_pos][i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i]];
#ifdef E_DAIF
    FCoef[sub_pos][SQR_FILTER-1] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][SQR_FILTER-1];
#endif
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
*       FCoef[][] - Filter to be extended at symmetric positions,
*                   using FILTER_NEW_SUB_POS
************************************************************************
*/
void ExtendFilterCoefficientsInt(int sub_pos, int ICoef[15][SQR_FILTER])
{
  int i;
  if(SymmetryPosition[sub_pos] == 0)    //check if it is allowed to change original filter coefficients
  {
    for(i = 0; i < SQR_FILTER; i++)
      ICoef[sub_pos][i] = ICoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i]];
  }
}


/*!
************************************************************************
* \brief

*      adapt 1D filter to 2D filter by inserting zeros in 1st,2nd,4th,5th and 6th raw or column
*      and copying the 3d one
*
*    input:
*      sub_pos - sub-position
*     input output:
*       FCoef[][] - Filter to be extended at symmetric positions,
*                     using FILTER_NEW_SUB_POS
************************************************************************
*/
void Position1DFilter(int sub_pos, double FCoef[SQR_FILTER])
{
  int i;
  if((sub_pos+1) % 4 == 0)     // vertical
  {
    for(i = SQR(FILTER_SIZE)-1; i>=0; i--)
      if(i%FILTER_SIZE == FILTER_OFFSET)
        FCoef[i]=FCoef[i/FILTER_SIZE];
      else
        FCoef[i]=0.0;
  }
  else if((sub_pos+1) / 4 == 0)     // horizontal
  {
    for(i = SQR(FILTER_SIZE)-1; i>=0; i--)
      if(i/FILTER_SIZE == FILTER_OFFSET)
        FCoef[i]=FCoef[i%FILTER_SIZE];
      else
        FCoef[i]=0.0;
  }
#ifdef DIRECTIONAL_FILTER
  else if (IsDiagonal1D[sub_pos]==1)
  {
    int j;
    for(i=0;i< FILTER_SIZE;i++)
      for(j=0;j<FILTER_SIZE;j++)
        if (i!=j)
          FCoef[i*FILTER_SIZE+j]=0.0;
  }
  else if (IsDiagonal1D[sub_pos]==2)
  {
    int j;
    for(i=0;i<FILTER_SIZE;i++)
      for(j=FILTER_SIZE-1;j>=0;j--)
        if (i!=(FILTER_SIZE-1 - j))
          FCoef[i*FILTER_SIZE+j]=0.0;
  }  
  else if (IsDiagonal1D[sub_pos]==3)
  {
    int j;
    for(i=0;i<FILTER_SIZE;i++)
      for(j=FILTER_SIZE-1;j>=0;j--)
        if (!((i==j)||(i==FILTER_SIZE-1 - j)))
          FCoef[i*FILTER_SIZE+j]=0.0;
  }

#endif
  else
    return;
}

#ifdef E_DAIF
void FindFilterCoefInt(void)
{
  int i, j, k;
  double **alpha, **beta;
  double *y;
  int N = SQR_FILTER_INT;

  get_mem1Ddouble(&y,N);
  get_mem2Ddouble(&alpha,N,N);
  get_mem2Ddouble(&beta,N,N);

  {
    FilterFlagInt=1;
    //1.Step : LU Decomposition

    //initialize the ii-components
    for(j = 0; j < N; j++)
    {
      for(i = 0; i < N; i++)
      {
        beta[i][j] = 0.0;
        alpha[i][j] = i==j ? 1.0:0.0;
      }
      y[j] = 0.0;
      CalculatedFilterInt[j] = 0.0;
    }
    for(j = 0; j < N; j++)
    {
      for(i = 0; i <=j; i++)
      {
        beta[i][j] = NxN_Matrix_Int[i][j];
        for(k = 0; k < i; k++)
        {
          beta[i][j] -= alpha[i][k]*beta[k][j];
        }
      }
      if(fabs(beta[j][j])<0.000001)
        FilterFlagInt=0;
      else
      {
        for(i = j+1; i < N; i++)
        {
          alpha[i][j] = NxN_Matrix_Int[i][j];
          for(k = 0; k < j; k++)
          {
            alpha[i][j] -= alpha[i][k]*beta[k][j];
          }
          alpha[i][j] = alpha[i][j]/beta[j][j];
        }
      }
    }

    if(FilterFlagInt == 0)
    {
      // default to no filtering on integer position
      for(i = 0; i < SQR_FILTER_INT; i++)
        FilterCoefInt[i]= 0;
      FilterCoefInt[(FILTER_SIZE_INT/2)*FILTER_SIZE_INT+(FILTER_SIZE_INT/2)]= 1.; 
    }
    else 
    {
      //2. Step: Solving L*y = b, forward substitution
      y[0]=N_Vector_Int[0];
      for(i = 1; i < N; i++)
      {
        y[i] = N_Vector_Int[i];
        for(j = 0; j < i; j++)
        {
          y[i] -= alpha[i][j]*y[j];
        }
      }
      //3. Step: Solving U*x = y, backsubstitution
      CalculatedFilterInt[N-1] = y[N-1]/beta[N-1][N-1];
      for(i = N-2; i >= 0; i--)
      {
        CalculatedFilterInt[i] = y[i];
        for(j = i+1; j < N; j++)
        {
          CalculatedFilterInt[i] -= beta[i][j]*CalculatedFilterInt[j];
        }
        CalculatedFilterInt[i] = CalculatedFilterInt[i]/beta[i][i];
      }
      for(i = 0; i < SQR_FILTER_INT; i++)
        FilterCoefInt[i] = CalculatedFilterInt[i];
    }
  }

  free_mem1Ddouble(y);
  free_mem2Ddouble(alpha);
  free_mem2Ddouble(beta);
}

static void PropogateFilterFlag()
{
  FilterFlag[c_pos] = FilterFlag[d_pos] = FilterFlag[l_pos] = FilterFlag[a_pos];
  FilterFlag[h_pos] = FilterFlag[b_pos];
  FilterFlag[g_pos] = FilterFlag[m_pos] = FilterFlag[o_pos] = FilterFlag[e_pos];
  FilterFlag[i_pos] = FilterFlag[k_pos] = FilterFlag[n_pos] = FilterFlag[f_pos];
}

static void ZeroFilterOffsets()
{
  int sub_pos;

  for(sub_pos = 0; sub_pos < 15; sub_pos++)
  {
    FilterCoef[sub_pos][SQR_FILTER-1] = 0.;
  }
}
#endif

#ifdef EAIF
// flip the coefficients horizontally
static void flipHoriz(double *coeff, double *coeffT)
{
  int i, j;
  for(i = 0; i < FILTER_SIZE; i++)
    for(j = 0; j < FILTER_SIZE; j++)
      coeffT[i*FILTER_SIZE+j] = coeff[i*FILTER_SIZE+(FILTER_SIZE-1-j)];
  coeffT[SQR_FILTER-1] = coeff[SQR_FILTER-1];
}

// flip the coefficients vertically
static void flipVert(double *coeff, double *coeffT)
{
  int i, j;
  for(i = 0; i < FILTER_SIZE; i++)
    for(j = 0; j < FILTER_SIZE; j++)
      coeffT[i*FILTER_SIZE+j] = coeff[(FILTER_SIZE-1-i)*FILTER_SIZE+j];
  coeffT[SQR_FILTER-1] = coeff[SQR_FILTER-1];
}

// propogate the FilterFlag from SymmetryPositions to Non-symmetry positions
static void PropogateFilterFlagEAIF()
{
  FilterFlag[c_pos] = FilterFlag[a_pos];
  FilterFlag[g_pos] = FilterFlag[e_pos];
  FilterFlag[k_pos] = FilterFlag[i_pos];
  FilterFlag[l_pos] = FilterFlag[d_pos];
  FilterFlag[m_pos] = FilterFlag[e_pos];
  FilterFlag[n_pos] = FilterFlag[f_pos];
  FilterFlag[o_pos] = FilterFlag[e_pos];
}
#endif

/*!
************************************************************************
* \brief
*    finds filter coefficients, solving the equation
*                      NxN_Matrix*FilterCoef=N_Vector
************************************************************************
*/
void FindFilterCoef(void)
{
  int i, j, k;
  double **alpha, **beta;
  double *y;
  int N = SQR(FILTER_SIZE);
  int sub_pos;
  if(!UseAdaptiveFilterForCurrentFrame())
    return;

#ifdef E_DAIF
  N = SQR_FILTER; // Offset 
#endif

  get_mem1Ddouble(&y,N);
  get_mem2Ddouble(&alpha,N,N);
  get_mem2Ddouble(&beta,N,N);
  for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
  {
    if(!SymmetryPosition[sub_pos])
    {
      if((FilterFlag[sub_pos] = FilterFlag[FILTER_NEW_SUB_POS[sub_pos]]))
      {
        if(Is1DPosition[sub_pos])
        {
          for(i = 0; i < FILTER_SIZE; i++)
          {
            FilterCoef[sub_pos][i]=
              CalculatedFilter[FILTER_NEW_SUB_POS[sub_pos]][TwoDEquationPattern[sub_pos][i]];
          }
#ifdef E_DAIF
          FilterCoef[sub_pos][SQR_FILTER-1]=
            CalculatedFilter[FILTER_NEW_SUB_POS[sub_pos]][POS_EQUATION_NUMBER[sub_pos]];
#endif
        }
        else
        {
          for(i = 0; i < SQR(FILTER_SIZE); i++)
          {
#ifdef EAIF
            if(input->UseAdaptiveFilter != FILTER_TYPE_EAIF || UseEquation[0][i])
#endif
              FilterCoef[sub_pos][i]=
              CalculatedFilter[FILTER_NEW_SUB_POS[sub_pos]][TwoDEquationPattern[sub_pos][i]];
#ifdef EAIF
            else if(input->UseAdaptiveFilter == FILTER_TYPE_EAIF)
              FilterCoef[sub_pos][i]=0.;
#endif
          }
#ifdef E_DAIF
          FilterCoef[sub_pos][SQR_FILTER-1]=
            CalculatedFilter[FILTER_NEW_SUB_POS[sub_pos]][POS_EQUATION_NUMBER[sub_pos]];
#endif
        }
      }
#ifdef E_DAIF
      else if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)
      {
        for(i = 0; i < SQR_FILTER; i++)
          FilterCoef[sub_pos][i]= STANDARD_2D_FILTER_orig[sub_pos][i];
      }
#endif
      else
        for(i = 0; i < SQR_FILTER; i++)
          FilterCoef[sub_pos][i]= STANDARD_2D_FILTER[sub_pos][i];

      continue;
    }

    N = POS_EQUATION_NUMBER[sub_pos];
#ifdef E_DAIF
#ifdef EAIF
    if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF || input->UseAdaptiveFilter == FILTER_TYPE_EAIF)
#else
    if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)
#endif
      N++;  // Offset
#endif 

    FilterFlag[sub_pos]=1;
    //1.Step : LU Decomposition

    //initialize the ii-components
    for(j = 0; j < N; j++)
    {
      for(i = 0; i < N; i++)
      {
        beta[i][j] = 0.0;
        alpha[i][j] = i==j ? 1.0:0.0;
      }
      y[j] = 0.0;
      CalculatedFilter[sub_pos][j] = 0.0;
    }
    for(j = 0; j < N; j++)
    {
      for(i = 0; i <=j; i++)
      {
        beta[i][j] = NxN_Matrix[sub_pos][i][j];
        for(k = 0; k < i; k++)
        {
          beta[i][j] -= alpha[i][k]*beta[k][j];
        }
      }
      if(fabs(beta[j][j])<0.000001)
        FilterFlag[sub_pos]=0;
      else
      {
        for(i = j+1; i < N; i++)
        {
          alpha[i][j] = NxN_Matrix[sub_pos][i][j];
          for(k = 0; k < j; k++)
          {
            alpha[i][j] -= alpha[i][k]*beta[k][j];
          }
          alpha[i][j] = alpha[i][j]/beta[j][j];
        }
      }
    }
    /*
    for(i = 0; i < N; i++)
    for(j = 0; j < N; j++)
    A_decompose[i][j] = (i<=j ? beta[i][j]:alpha[i][j]);
    */
    if(FilterFlag[sub_pos]==0)
    {
#ifdef E_DAIF
      for(i = 0; i < SQR(FILTER_SIZE); i++)
        FilterCoef[sub_pos][i]= 
        input->UseAdaptiveFilter == FILTER_TYPE_EDAIF ? STANDARD_2D_FILTER_orig[sub_pos][i]:STANDARD_2D_FILTER[sub_pos][i];  
      FilterCoef[sub_pos][i] = 0.0;
#else
      for(i = 0; i < SQR_FILTER; i++)
        FilterCoef[sub_pos][i]= STANDARD_2D_FILTER[sub_pos][i];
#endif  // E_DAIF
      continue;
    }
    //2. Step: Solving L*y = b, forward substitution
    y[0]=N_Vector[sub_pos][0];
    for(i = 1; i < N; i++)
    {
      y[i] = N_Vector[sub_pos][i];
      for(j = 0; j < i; j++)
      {
        y[i] -= alpha[i][j]*y[j];
      }
    }
    //3. Step: Solving U*x = y, backsubstitution
    CalculatedFilter[sub_pos][N-1] = y[N-1]/beta[N-1][N-1];
    for(i = N-2; i >= 0; i--)
    {
      CalculatedFilter[sub_pos][i] = y[i];
      for(j = i+1; j < N; j++)
      {
        CalculatedFilter[sub_pos][i] -= beta[i][j]*CalculatedFilter[sub_pos][j];
      }
      CalculatedFilter[sub_pos][i] = CalculatedFilter[sub_pos][i]/beta[i][i];
    }
#ifdef  DIRECTIONAL_FILTER
#ifdef E_DAIF
    if (input->UseAdaptiveFilter == FILTER_TYPE_1D)
#endif
      ApplySumRestriction(input->ImpType, sub_pos);
    // FilterFlag is changed in ApplySumRestriction
#ifdef EAIF
    if(input->UseAdaptiveFilter == FILTER_TYPE_EAIF)
    {
      for(i = 0; i < SQR(FILTER_SIZE); i++)
      {
        FilterCoef[sub_pos][i]=CalculatedFilter[sub_pos][TwoDEquationPattern[sub_pos][i]];
      }
      FilterCoef[sub_pos][i]=CalculatedFilter[sub_pos][N-1];
    }
    else
#endif
      if(FilterFlag[sub_pos]!=0)
      {
        for(i = 0; i < SQR_FILTER; i++)
        {
          FilterCoef[sub_pos][i]=
            CalculatedFilter[FILTER_NEW_SUB_POS[sub_pos]][TwoDEquationPattern[sub_pos][i]];
        }
      }
      else
      {
        for(i = 0; i < SQR_FILTER; i++)
          FilterCoef[sub_pos][i]= STANDARD_2D_FILTER[sub_pos][i];
      }
#else
    for(i = 0; i < SQR_FILTER; i++)
      FilterCoef[sub_pos][i]=
      CalculatedFilter[FILTER_NEW_SUB_POS[sub_pos]][TwoDEquationPattern[sub_pos][i]];
#endif
  }

  // positioning 1D filter
  for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
  {
#ifndef DIRECTIONAL_FILTER
    if(FilterFlag[sub_pos])
      Position1DFilter(sub_pos, FilterCoef[sub_pos]);
#else
    if((FilterFlag[sub_pos])&&(Is1DPosition[sub_pos]))
      Position1DFilter(sub_pos, FilterCoef[sub_pos]);
#endif
  }
  // calculate UseAllPositions
  SubpelPositionsPattern =   (FilterFlag[j_pos]<<4) + (FilterFlag[f_pos]<<3) + (FilterFlag[e_pos]<<2) + 
    (FilterFlag[b_pos]<<1) +  FilterFlag[a_pos];
  if( SubpelPositionsPattern == 31)
    UseAllSubpelPositions = 1;          // for indicating that every sub-pel position is used
  else
    UseAllSubpelPositions  = 0;

#ifndef DIRECTIONAL_FILTER
  PredictFilterCoefficients();
#else
  if(input->UseAdaptiveFilter==FILTER_TYPE_2D_NS)
    PredictFilterCoefficients();
  else if(input->UseAdaptiveFilter==FILTER_TYPE_1D)
  {
    PredictFilterCoefficients1DAIF();
    if (input->ImpType==IMP_INT16)
      RepresentCoefsIn16bits(input->ImpType);
  }
#endif
#ifdef E_DAIF
  else if(input->UseAdaptiveFilter==FILTER_TYPE_EDAIF)
  {
    PredictFilterCoefficients1DAIF();
  }
#ifdef EAIF
  else if(input->UseAdaptiveFilter == FILTER_TYPE_EAIF)
  {
    PredictFilterCoefficientsEAIF();
  }
#endif

#ifdef EAIF
  if(input->UseAdaptiveFilter != FILTER_TYPE_EDAIF && input->UseAdaptiveFilter != FILTER_TYPE_EAIF)
#else
  if(input->UseAdaptiveFilter != FILTER_TYPE_EDAIF)
#endif
  {
    ZeroFilterOffsets();
  }
#endif
  //  PrintFilterCoefInt(DiffQFilterCoef);
  //  PrintFilterCoefDouble(PredFilterCoef);
  //  PrintFilterCoefDouble(FilterCoef);

  free_mem1Ddouble(y);
  free_mem2Ddouble(alpha);
  free_mem2Ddouble(beta);

#ifdef E_DAIF
#ifdef EAIF
  if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF || input->UseAdaptiveFilter == FILTER_TYPE_EAIF)
#else
  if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)
#endif
  {
    FindFilterCoefInt(); 
    QuantizeFilterOffsets(); 

    if(FilterFlagInt)
      QuantizeIntegerFilter(); 
  }

  if(input->UseAdaptiveFilter!=FILTER_TYPE_2D_NS && input->UseAdaptiveFilter != FILTER_TYPE_2D_S)
  {
    // Filter decision    
    AIFDecision(&FilterFlagInt, FilterFlag, SymmetryPosition);

#ifdef EAIF
    if(input->UseAdaptiveFilter == FILTER_TYPE_EAIF)
      PropogateFilterFlagEAIF();
    else
#endif
      PropogateFilterFlag();

    for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
    {
      if(FilterFlag[sub_pos] == 0)
      {
        for(i = 0; i < SQR(FILTER_SIZE); i++)
          FilterCoef[sub_pos][i] = 
          input->UseAdaptiveFilter == FILTER_TYPE_EDAIF ? STANDARD_2D_FILTER_orig[sub_pos][i] : STANDARD_2D_FILTER[sub_pos][i];
        FilterCoef[sub_pos][SQR_FILTER-1] = 0.;
      }
    }
  }

#ifdef EAIF
  if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF || input->UseAdaptiveFilter == FILTER_TYPE_EAIF)
#else
  if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)
#endif
  {
    if(!FilterFlagInt)
    {
      for(i = 0; i < SQR_FILTER_INT; ++i)
        FilterCoefInt[i] = 0.0; 
      FilterCoefInt[12] = 1.0; 
    }
  }

  {
    SubpelPositionsPattern =   (FilterFlag[j_pos]<<4) + (FilterFlag[f_pos]<<3) + (FilterFlag[e_pos]<<2) + 
      (FilterFlag[b_pos]<<1) +  FilterFlag[a_pos];
    if( SubpelPositionsPattern == 31)
      UseAllSubpelPositions = 1;          // for indicating that every sub-pel position is used
    else
      UseAllSubpelPositions  = 0;
  }
#endif
}

#ifdef EDAIF2
void FindFilterCoef2(void)
{
  int i, j, k;
  double **alpha, **beta;
  double *y;
  int N = SQR(FILTER_SIZE);
  int sub_pos;
  if(!UseAdaptiveFilterForCurrentFrame())
    return;

#ifdef E_DAIF
  N = SQR_FILTER; // Offset 
#endif

  get_mem1Ddouble(&y,N);
  get_mem2Ddouble(&alpha,N,N);
  get_mem2Ddouble(&beta,N,N);

  memset(CalculatedFilter, 0, SQR_FILTER*MAX_NUM_AIF*sizeof(double));
  for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
  {
    if((!SymmetryPosition[sub_pos])||(FilterFlag[sub_pos]%2==0))
    {// if filter is static, or it is to be copied from another location - skip it!
      continue;
    }

    N = POS_EQUATION_NUMBER[sub_pos]+1;

    //FilterFlag[sub_pos]=1;
    //1.Step : LU Decomposition

    //initialize the ii-components
    for(j = 0; j < N; j++)
    {
      for(i = 0; i < N; i++)
      {
        beta[i][j] = 0.0;
        alpha[i][j] = i==j ? 1.0:0.0;
      }
      y[j] = 0.0;
      CalculatedFilter[sub_pos][j] = 0.0;
    }
    for(j = 0; j < N; j++)
    {
      for(i = 0; i <=j; i++)
      {
        beta[i][j] = NxN_Matrix[sub_pos][i][j];
        for(k = 0; k < i; k++)
        {
          beta[i][j] -= alpha[i][k]*beta[k][j];
        }
      }
      if(fabs(beta[j][j])<0.000001)
        FilterFlag[sub_pos]=0;// Filter flag initially either 3 or 1, so use static filter of the same structure
      else
      {
        for(i = j+1; i < N; i++)
        {
          alpha[i][j] = NxN_Matrix[sub_pos][i][j];
          for(k = 0; k < j; k++)
          {
            alpha[i][j] -= alpha[i][k]*beta[k][j];
          }
          alpha[i][j] = alpha[i][j]/beta[j][j];
        }
      }
    }
    if(FilterFlag[sub_pos]%2==0)
    {
      continue;
    }
    //2. Step: Solving L*y = b, forward substitution
    y[0]=N_Vector[sub_pos][0];
    for(i = 1; i < N; i++)
    {
      y[i] = N_Vector[sub_pos][i];
      for(j = 0; j < i; j++)
      {
        y[i] -= alpha[i][j]*y[j];
      }
    }
    //3. Step: Solving U*x = y, backsubstitution
    CalculatedFilter[sub_pos][N-1] = y[N-1]/beta[N-1][N-1];
    for(i = N-2; i >= 0; i--)
    {
      CalculatedFilter[sub_pos][i] = y[i];
      for(j = i+1; j < N; j++)
      {
        CalculatedFilter[sub_pos][i] -= beta[i][j]*CalculatedFilter[sub_pos][j];
      }
      CalculatedFilter[sub_pos][i] = CalculatedFilter[sub_pos][i]/beta[i][i];
    }
  }

  FindFilterCoefInt(); 
  if(FilterFlagInt)
    QuantizeIntegerFilter(); 

  PredictFilterCoefficientsDAIF2(1);
  QuantizeFilterOffsets(); 


  AIFDecision2(&FilterFlagInt, FilterFlag, SymmetryPosition,FILTER_NEW_SUB_POS);  // Shall be rewritten to support variable symmetrices!

  PredictFilterCoefficientsDAIF2(2);// do the final filter propagation
  QuantizeFilterOffsets(); 


  if(!FilterFlagInt)
  {
    for(i = 0; i < SQR_FILTER_INT; ++i)
      FilterCoefInt[i] = 0.0; 
    FilterCoefInt[12] = 1.0; 
  }

  free_mem1Ddouble(y);
  free_mem2Ddouble(alpha);
  free_mem2Ddouble(beta);
}
void FilterDecomposition(int filterType)
{    
  int i,j,k,counter;
  double thresh1;
  double diff1, abs_sum, f1,f2; 
  int h1, h2;

  int * SymProf = NULL;
  int SymProfCount[4] = {0,0,0,0};
  double SymProfDist[4] = {0.0,0.0,0.0,0.0};
  int * p_SymProf[4] = {&SetHor[0][0],&SetVert[0][0],&SetDiagNW[0][0],&SetDiagNE[0][0]};

  FilterEntity *p_Filter = &SetAIF[filterType];

  p_Filter->RealTimeDecomposition = 0;

  thresh1 = (g_nBitsThresh/20000)*(g_nBitsThresh/20000)+1;
  if (thresh1 <  1) 
    thresh1 = 1;
  g_AIFDecompThresh = thresh1;
  thresh1 = 1/g_AIFDecompThresh;
  //    printf("%f\n",thresh1);

  memset(p_Filter->realSymCommands,0,4*MAX_NUM_AIF*sizeof(int));

  for(k = 0;k<4;k++)
  {
    SymProf  = p_SymProf[k];
    SymProfCount[k] = 0;

    for (j = 0; j < nPairsSets[k]; j ++)
    {
      diff1 = 0.0;
      abs_sum  = 0.0;
      h1 = *(SymProf + j*2 +0);
      h2 = *(SymProf+ 2*j + 1);
      if ((p_Filter->FilterFlag[h1]%2==0)||(p_Filter->FilterFlag[h2]%2==0))
      {
        p_Filter->realSymCommands[k][j] = 1;
        continue;
      }
      if (p_Filter->FilterFlag[h1]!=p_Filter->FilterFlag[h2])
      {
        continue;
      }

      for(i = 0; i < p_Filter->POS_EQUATION_NUMBER[h1]; i++)
      {
        f1 = p_Filter->CalculatedFilter[h1][i];
        f2 = p_Filter->CalculatedFilter[h2][i];
        diff1 += (f1 - f2)*(f1- f2);
      }

      diff1  = sqrt(diff1);
      if (diff1 < thresh1)
      {
        p_Filter->realSymCommands[k][j] = 1;
      }else
        p_Filter->realSymCommands[k][j] = 0;

      SymProfDist[k]+=diff1;
      SymProfCount[k]++;
    }
    if(SymProfCount[k]!=0)
      SymProfDist[k] = SymProfDist[k]/SymProfCount[k];
  }

  counter = 0;
  for(k = 0;k<4;k++)
  {
    SymProf  = p_SymProf[k];
    for (j = 0; j < nPairsSets[k]; j ++)
    {
      h1 = p_Filter->FILTER_NEW_SUB_POS[*(SymProf + j*2 +0)];
      h2 = *(SymProf+ 2*j + 1);
      if(p_Filter->realSymCommands[k][j])
      {
        p_Filter->FILTER_NEW_SUB_POS[h2] = h1;
        counter++;
      }
    }
  }

  for (j=0;j<=o_pos;j++)
  {
    if(input->ImpType== 0)
    {
      if(counter>16)
        p_Filter->nBitsIntRepresent[j] = DEFAULT_QUANT-1;
      if(counter==20)
        p_Filter->nBitsIntRepresent[j] = DEFAULT_QUANT-2;
    }

    if(p_Filter->FILTER_NEW_SUB_POS[j] != j)
    {
      p_Filter->SymmetryPosition[j] = 0;
    }
  }


  for (j = 0; j <=o_pos; j ++)
  {
    h1 = p_Filter->FILTER_NEW_SUB_POS[p_Filter->FILTER_NEW_SUB_POS[j]];
    h2 = j;
    if (!p_Filter->SymmetryPosition[j])
    {
      for (k = 0; k<SQR_FILTER; k++)
      {
        p_Filter->N_Vector[h1][k] += p_Filter->N_Vector[h2][k];
        p_Filter->N_Vector[h2][k] = 0.0;
        for (i = 0; i<SQR_FILTER; i++)
        {
          p_Filter->NxN_Matrix[h1][k][i] += p_Filter->NxN_Matrix[h2][k][i];
          p_Filter->NxN_Matrix[h2][k][i] = 0.0;
        }
      }
    }
  }

  p_Filter->RealTimeDecomposition = 1;

}

void FilterAssignment(int filterType)
{
  FilterEntity *p_Filter = &SetAIF[filterType];

  //    Copy Global Filter Configuration into the SetAIF[ID]
  p_Filter->RealTimeDecomposition = 0;
  memset(p_Filter->realSymCommands,0,4*MAX_NUM_AIF*sizeof(int));
  memcpy(p_Filter->N_Vector, N_Vector, SQR_FILTER*(MAX_NUM_AIF)*sizeof(double));
  memcpy(p_Filter->NxN_Matrix, NxN_Matrix, SQR_FILTER*SQR_FILTER*(MAX_NUM_AIF)*sizeof(double));

  memcpy(p_Filter->FilterCoef, FilterCoef,SQR_FILTER*MAX_NUM_AIF*sizeof(double));
  memcpy(p_Filter->CalculatedFilter, CalculatedFilter,SQR_FILTER*MAX_NUM_AIF*sizeof(double));
  memcpy(p_Filter->FILTER_NEW_SUB_POS, FILTER_NEW_SUB_POS_DAIF,sizeof(int)*MAX_NUM_AIF);
  memcpy(p_Filter->SymmetryPosition, SymmetryPosition,sizeof(int)*MAX_NUM_AIF);
  memcpy(p_Filter->FilterFlag, FilterFlag,sizeof(int)*MAX_NUM_AIF);
  memcpy(p_Filter->POS_EQUATION_NUMBER,POS_EQUATION_NUMBER,sizeof(int)*MAX_NUM_AIF);

  if (filterType == 1)
  {
    // we have SetAIF[1] already setted to RAIF, up_date the statistics and filter coefs.
    memcpy(p_Filter->FilterCoef[f_pos],FilterCoef[f2_pos], SQR_FILTER*sizeof(double));
    memcpy(p_Filter->FilterCoef[i_pos],FilterCoef[i2_pos], 3*SQR_FILTER*sizeof(double));
    memcpy(p_Filter->FilterCoef[n_pos],FilterCoef[n2_pos], SQR_FILTER*sizeof(double));

    memcpy(p_Filter->CalculatedFilter[f_pos],CalculatedFilter[f2_pos], SQR_FILTER*sizeof(double));
    memcpy(p_Filter->CalculatedFilter[i_pos],CalculatedFilter[i2_pos], 3*SQR_FILTER*sizeof(double));
    memcpy(p_Filter->CalculatedFilter[n_pos],CalculatedFilter[n2_pos], SQR_FILTER*sizeof(double));

    memcpy(p_Filter->N_Vector[f_pos],N_Vector[f2_pos], SQR_FILTER*sizeof(double));
    memcpy(p_Filter->N_Vector[i_pos],N_Vector[i2_pos], 3*SQR_FILTER*sizeof(double));
    memcpy(p_Filter->N_Vector[n_pos],N_Vector[n2_pos], SQR_FILTER*sizeof(double));

    memcpy(p_Filter->NxN_Matrix[f_pos],NxN_Matrix[f2_pos], SQR_FILTER*SQR_FILTER*sizeof(double));
    memcpy(p_Filter->NxN_Matrix[i_pos],NxN_Matrix[i2_pos], 3*SQR_FILTER*SQR_FILTER*sizeof(double));
    memcpy(p_Filter->NxN_Matrix[n_pos],NxN_Matrix[n2_pos], SQR_FILTER*SQR_FILTER* sizeof(double));

    p_Filter->FilterFlag[f_pos] = 3*FilterFlag[f2_pos];
    p_Filter->FilterFlag[i_pos] = 3*FilterFlag[i2_pos];
    p_Filter->FilterFlag[j_pos] = 3*FilterFlag[j2_pos];
    p_Filter->FilterFlag[k_pos] = 3*FilterFlag[k2_pos];
    p_Filter->FilterFlag[n_pos] = 3*FilterFlag[n2_pos];
  }
}
void selectAIF_MCP()
{
  int sub_pos, new_sup_pos;
  FilterEntity *p_Filter = &SetAIF[1];


  SwapAIF(0,0);// copy DAIF into global storage
  // Compare DAIF vs. RAIF for 5,8,9,19,13 sub-pels
  // and up-date global if RAIF is better.
  for(sub_pos = 0; sub_pos<=o_pos;sub_pos++)
  {
    if(((sub_pos==f_pos)||(sub_pos==i_pos)||(sub_pos==j_pos)||(sub_pos==k_pos)||(sub_pos==n_pos)))
    {
      if (MIN(costMCP[2][sub_pos],costMCP[3][sub_pos])< MIN(costMCP[0][sub_pos],costMCP[1][sub_pos]))
      {// Use RAIF
        p_Filter = &SetAIF[1];
        FilterFlag[sub_pos] = FilterFlag[sub_pos]*3;
      }
      else
      {
        p_Filter = &SetAIF[0];
        FilterFlag[sub_pos] = FilterFlag[sub_pos]*1;
      }
      //Copy all from from NEW_POS_SYM location!
      new_sup_pos = p_Filter->FILTER_NEW_SUB_POS[FILTER_NEW_SUB_POS[sub_pos]];

      //  Cancel found symmetries
      SymmetryPosition[sub_pos]=1;// we will re-do symmetry check
      memcpy(CalculatedFilter[sub_pos], p_Filter->CalculatedFilter[new_sup_pos],SQR_FILTER*sizeof(double));
      memcpy(N_Vector[sub_pos], p_Filter->N_Vector[new_sup_pos],SQR_FILTER*sizeof(double));
      memcpy(NxN_Matrix[sub_pos], p_Filter->NxN_Matrix[new_sup_pos],SQR_FILTER*SQR_FILTER*sizeof(double));

      memcpy(FilterCoef[sub_pos], p_Filter->FilterCoef[sub_pos],SQR_FILTER*sizeof(double));        
      memcpy(TwoDEquationPattern[sub_pos],p_Filter->TwoDEquationPattern[sub_pos],sizeof(int)*SQR_FILTER);
      memcpy(Calc2Filt_Indexes[sub_pos],p_Filter->Calc2Filt_Indexes[sub_pos],sizeof(int)*SQR_FILTER);
      memcpy(STANDARD_2D_FILTER[sub_pos],p_Filter->STANDARD_2D_FILTER[sub_pos],sizeof(double)*SQR_FILTER);
      IsDiagonal1D[sub_pos] = p_Filter->IsDiagonal1D[sub_pos];
      Is1DPosition[sub_pos] = p_Filter->Is1DPosition[sub_pos];
      Filt_Indexes_offset[sub_pos] = p_Filter->Filt_Indexes_offset[sub_pos];
    }
  }
  SwapAIF(2,1);
}

void PreProcessFilterCoef(void)
{
  int i, j, k;
  double **alpha, **beta;
  double *y;
  int N = SQR_FILTER;
  int sub_pos;
  int counter = 0;

  if(!UseAdaptiveFilterForCurrentFrame())
    return;

  get_mem1Ddouble(&y,N);
  get_mem2Ddouble(&alpha,N,N);
  get_mem2Ddouble(&beta,N,N);
  memset(CalculatedFilter, 0, SQR_FILTER*(MAX_NUM_AIF+MAX_NUM_AIF_EXTRA)*sizeof(double));
  for(sub_pos = a_pos; sub_pos <= n2_pos; sub_pos++)
  {
    N = POS_EQUATION_NUMBER[sub_pos]+1;// compute required number of taps + filter DC offset
    FilterFlag[sub_pos]=1;
    //1.Step : LU Decomposition

    //initialize the ii-components
    for(j = 0; j < N; j++)
    {
      for(i = 0; i < N; i++)
      {
        beta[i][j] = 0.0;
        alpha[i][j] = i==j ? 1.0:0.0;
      }
      y[j] = 0.0;
      CalculatedFilter[sub_pos][j] = 0.0;
    }
    for(j = 0; j < N; j++)
    {
      for(i = 0; i <=j; i++)
      {
        beta[i][j] = NxN_Matrix[sub_pos][i][j];
        for(k = 0; k < i; k++)
        {
          beta[i][j] -= alpha[i][k]*beta[k][j];
        }
      }
      if(fabs(beta[j][j])<0.000001)
        FilterFlag[sub_pos]=0;
      else
      {
        for(i = j+1; i < N; i++)
        {
          alpha[i][j] = NxN_Matrix[sub_pos][i][j];
          for(k = 0; k < j; k++)
          {
            alpha[i][j] -= alpha[i][k]*beta[k][j];
          }
          alpha[i][j] = alpha[i][j]/beta[j][j];
        }
      }
    }

    if(FilterFlag[sub_pos]==0)
    {
      continue;
    }
    //2. Step: Solving L*y = b, forward substitution
    y[0]=N_Vector[sub_pos][0];
    for(i = 1; i < N; i++)
    {
      y[i] = N_Vector[sub_pos][i];
      for(j = 0; j < i; j++)
      {
        y[i] -= alpha[i][j]*y[j];
      }
    }
    //3. Step: Solving U*x = y, backsubstitution
    CalculatedFilter[sub_pos][N-1] = y[N-1]/beta[N-1][N-1];
    for(i = N-2; i >= 0; i--)
    {
      CalculatedFilter[sub_pos][i] = y[i];
      for(j = i+1; j < N; j++)
      {
        CalculatedFilter[sub_pos][i] -= beta[i][j]*CalculatedFilter[sub_pos][j];
      }
      CalculatedFilter[sub_pos][i] = CalculatedFilter[sub_pos][i]/beta[i][i];
    }

    if (FilterFlag[sub_pos])
      counter ++;
  }
  PredictFilterCoefficientsDAIF2(0);

  if (counter!=0)
  {
    FilterAssignment(0);
    GetErrorMCP(&FilterFlagInt,0);
    FilterAssignment(1);
    GetErrorMCP(&FilterFlagInt,1);
    FilterAssignment(2);
    selectAIF_MCP();
    FilterDecomposition(2);
  }

  free_mem1Ddouble(y);
  free_mem2Ddouble(alpha);
  free_mem2Ddouble(beta);
}

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

void PredictFilterCoefficientsDAIF2(int mode)
{
  int i,j,counter;
  double PredFilterCoef[MAX_NUM_AIF + MAX_NUM_AIF_EXTRA][SQR_FILTER];	
  double DiffFilterCoef[MAX_NUM_AIF + MAX_NUM_AIF_EXTRA][SQR_FILTER];
  //int TempArray[SQR_FILTER];
  int NUM_SUB_POS = (mode==0) ? MAX_NUM_AIF + MAX_NUM_AIF_EXTRA : MAX_NUM_AIF;

  counter = 0;
  memset(FilterCoef, 0, SQR_FILTER*(MAX_NUM_AIF + MAX_NUM_AIF_EXTRA)*sizeof(double));
  memset(DiffQFilterCoef, 0, SQR_FILTER*(MAX_NUM_AIF + MAX_NUM_AIF_EXTRA)*sizeof(int));
  memset(FilterCoef16bits, 0, SQR_FILTER*(MAX_NUM_AIF + MAX_NUM_AIF_EXTRA)*sizeof(short int));


  // Checkig then coefs on the range conditions
  if(mode!=0)
  {
    for(j = 0; j < NUM_SUB_POS; j++)
    {
      // Check 1: Is the mother sub-pos is active?
      if (FilterFlag[FILTER_NEW_SUB_POS[j]]%2==0)// Standard filters are used for FilterFlag = 0/2 
      {
        FilterFlag[j]= FilterFlag[FILTER_NEW_SUB_POS[j]];
        continue;
      }
      for(i = 0; i < POS_EQUATION_NUMBER[j]; i++)
      {
        if (fabs(CalculatedFilter[j][i])>1.0)
        {
          FilterFlag[j]= 0;  // if adaptive filter was used earlier, set the standard filter
        }
      }

      // Check 3: Are triplets of filter coefs out of range?
      ApplySumRestriction(input->ImpType, j);
    }
  }

  // Extend coefs from the symmetry
  for(j= 0; j< NUM_SUB_POS; j++)
  {
    if(FilterFlag[j]%2==0)
    {
      for(i = 0; i < POS_EQUATION_NUMBER[j]; i++)
        CalculatedFilter[j][i] = STANDARD_2D_FILTER[j][Calc2Filt_Indexes[j][i]];
      CalculatedFilter[j][POS_EQUATION_NUMBER[j]] = 0.0;
    }
    else
      counter++;
    ExtendDiagonalFilter_DAIF2(j, FilterCoef);
  }


  for(j=0; j < NUM_SUB_POS; j++)
  {
    NumberOfQBits = nBitsIntRepresent[FILTER_NEW_SUB_POS[FILTER_NEW_SUB_POS[j]]] + 1;
    //printf("%d ",NumberOfQBits);
    if(input->ImpType== 0)
    {
      for(i=0; i < SQR_FILTER; i++)
      {
        PredFilterCoef[j][i] = STANDARD_2D_FILTER[j][i];
        DiffFilterCoef[j][i] = FilterCoef[j][i] - PredFilterCoef[j][i];
      }

      // NOTE: we don't quantize here filter offset
      QuantizeFilterCoefficientsEDAIF2(DiffFilterCoef[j], DiffQFilterCoef[j],j);
      for(i=0; i < SQR_FILTER; i++)
        FilterCoef[j][i] = DiffFilterCoef[j][i] + PredFilterCoef[j][i];
    }
    else
    {
      printf("16-bit interpolation is currently not supported!\n");
    }
  }

  if(mode == 2)
  {
    for(j= 0; j< NUM_SUB_POS; j++)
    {
      if(FilterFlag[j]%2==0)
      {
        for(i = 0; i < SQR_FILTER-1; i++)
          FilterCoef[j][i] = STANDARD_2D_FILTER_orig[j][i];
        FilterCoef[j][SQR_FILTER-1] = 0.0;
      }
    }
  }

  if (counter==NUM_SUB_POS)
    UseAllSubpelPositions = 1;
  else 
    UseAllSubpelPositions = 0;

  if (counter ==0)
  {
    memcpy(FilterCoef,STANDARD_2D_FILTER_orig,sizeof(double)*NUM_SUB_POS*SQR_FILTER);
    return;
  }
}

#endif
// separable aif (BEGIN)
/*!
************************************************************************
* \brief
*    calculation of filter coefficients for horizontal filter 
************************************************************************
*/
void FindFilterCoefHor (void)
{
  int i, j, k;
  double **alpha, **beta;
  double *y;
  int N = FILTER_SIZE;
  int sub_pos;
  if (!UseAdaptiveFilterForCurrentFrame ())
    return;
  get_mem1Ddouble (&y, N);
  get_mem2Ddouble (&alpha, N, N);
  get_mem2Ddouble (&beta, N, N);
  for (sub_pos = a_pos; sub_pos <= c_pos; sub_pos++)
  {
    N = (sub_pos == b_pos) ? FILTER_SIZE>>1 : FILTER_SIZE;
    FilterFlag[sub_pos] = 1;
    //1.Step : LU Decomposition
    //initialize the ii-components
    for (j = 0; j < N; j++)
    {
      for (i = 0; i < N; i++)
      {
        beta[i][j] = 0.0;
        alpha[i][j] = i == j ? 1.0 : 0.0;
      }
      y[j] = 0.0;
      CalculatedFilter[sub_pos][j] = 0.0;
    }
    for (j = 0; j < N; j++)
    {
      for (i = 0; i <= j; i++)
      {
        beta[i][j] = NxN_Matrix[sub_pos][i][j];
        for (k = 0; k < i; k++)
        {
          beta[i][j] -= alpha[i][k] * beta[k][j];
        }
      }
      if (fabs (beta[j][j]) < 0.000001)
      {
        FilterFlag[sub_pos] = 0;
        for (i = 0; i < SQR_FILTER; i++) 
          FilterCoef[sub_pos][i] = STANDARD_2D_FILTER_SEP[sub_pos][i]; 
        continue; 
      }
      else
      {
        for (i = j + 1; i < N; i++)
        {
          alpha[i][j] = NxN_Matrix[sub_pos][i][j];
          for (k = 0; k < j; k++)
          {
            alpha[i][j] -= alpha[i][k] * beta[k][j];
          }
          alpha[i][j] = alpha[i][j] / beta[j][j];
        }
      }
    }
    //2. Step: Solving L*y = b, forward substitution
    y[0] = N_Vector[sub_pos][0];
    for (i = 1; i < N; i++)
    {
      y[i] = N_Vector[sub_pos][i];
      for (j = 0; j < i; j++)
      {
        y[i] -= alpha[i][j] * y[j];
      }
    }
    //3. Step: Solving U*x = y, backsubstitution
    CalculatedFilter[sub_pos][N - 1] = y[N - 1] / beta[N - 1][N - 1];
    for (i = N - 2; i >= 0; i--)
    {
      CalculatedFilter[sub_pos][i] = y[i];
      for (j = i + 1; j < N; j++)
      {
        CalculatedFilter[sub_pos][i] -= beta[i][j] * CalculatedFilter[sub_pos][j];
      }
      CalculatedFilter[sub_pos][i] = CalculatedFilter[sub_pos][i] / beta[i][i];
    }
    if (sub_pos == b_pos)
    {
      FilterCoef[b_pos][0] = CalculatedFilter[b_pos][0];
      FilterCoef[b_pos][1] = CalculatedFilter[b_pos][1];
      FilterCoef[b_pos][2] = CalculatedFilter[b_pos][2];
      FilterCoef[b_pos][3] = CalculatedFilter[b_pos][2];
      FilterCoef[b_pos][4] = CalculatedFilter[b_pos][1];
      FilterCoef[b_pos][5] = CalculatedFilter[b_pos][0];
    }
    else
    {
      for (i = 0; i < SQR_FILTER; i++)
      {
        FilterCoef[sub_pos][i] = CalculatedFilter[sub_pos][i];
      }
    }
  }
  PredictFilterCoefficientsHor ();        
  free_mem1Ddouble (y);
  free_mem2Ddouble (alpha);
  free_mem2Ddouble (beta);
}

/*!
************************************************************************
* \brief
*    calculation of filter coefficients for vertical filter 
************************************************************************
*/
void FindFilterCoefVer (void)
{
  int i, j, k;
  double **alpha, **beta;
  double *y;
  int N = FILTER_SIZE;
  int sub_pos;
  if (!UseAdaptiveFilterForCurrentFrame ())
    return;
  get_mem1Ddouble (&y, N);
  get_mem2Ddouble (&alpha, N, N);
  get_mem2Ddouble (&beta, N, N);
  for (sub_pos = d_pos; sub_pos <= k_pos; sub_pos++)
  {
    if (sub_pos == h_pos || sub_pos == i_pos || sub_pos == j_pos || sub_pos == k_pos)
      N = FILTER_SIZE>>1;          
    else
      N = FILTER_SIZE; 
    FilterFlag[sub_pos] = 1;
    //1.Step : LU Decomposition
    //initialize the ii-components
    for (j = 0; j < N; j++)
    {
      for (i = 0; i < N; i++)
      {
        beta[i][j] = 0.0;
        alpha[i][j] = i == j ? 1.0 : 0.0;
      }
      y[j] = 0.0;
      CalculatedFilter[sub_pos][j] = 0.0;
    }
    for (j = 0; j < N; j++)
    {
      for (i = 0; i <= j; i++)
      {
        beta[i][j] = NxN_Matrix[sub_pos][i][j];
        for (k = 0; k < i; k++)
        {
          beta[i][j] -= alpha[i][k] * beta[k][j];
        }
      }
      if (fabs (beta[j][j]) < 0.000001)
      {
        FilterFlag[sub_pos] = 0;
        for (i = 0; i < SQR_FILTER; i++)  
          FilterCoef[sub_pos][i] = STANDARD_2D_FILTER_SEP[sub_pos][i];
        continue;
      }
      else
      {
        for (i = j + 1; i < N; i++)
        {
          alpha[i][j] = NxN_Matrix[sub_pos][i][j];
          for (k = 0; k < j; k++)
          {
            alpha[i][j] -= alpha[i][k] * beta[k][j];
          }
          alpha[i][j] = alpha[i][j] / beta[j][j];
        }
      }
    }
    //2. Step: Solving L*y = b, forward substitution
    y[0] = N_Vector[sub_pos][0];
    for (i = 1; i < N; i++)
    {
      y[i] = N_Vector[sub_pos][i];
      for (j = 0; j < i; j++)
      {
        y[i] -= alpha[i][j] * y[j];
      }
    }
    //3. Step: Solving U*x = y, backsubstitution
    CalculatedFilter[sub_pos][N - 1] = y[N - 1] / beta[N - 1][N - 1];
    for (i = N - 2; i >= 0; i--)
    {
      CalculatedFilter[sub_pos][i] = y[i];
      for (j = i + 1; j < N; j++)
      {
        CalculatedFilter[sub_pos][i] -= beta[i][j] * CalculatedFilter[sub_pos][j];
      }
      CalculatedFilter[sub_pos][i] = CalculatedFilter[sub_pos][i] / beta[i][i];
    }
    if (sub_pos == h_pos || sub_pos == i_pos || sub_pos == j_pos || sub_pos == k_pos)
    {
      FilterCoef[sub_pos][0] = CalculatedFilter[sub_pos][0];
      FilterCoef[sub_pos][1] = CalculatedFilter[sub_pos][1];
      FilterCoef[sub_pos][2] = CalculatedFilter[sub_pos][2];
      FilterCoef[sub_pos][3] = CalculatedFilter[sub_pos][2];
      FilterCoef[sub_pos][4] = CalculatedFilter[sub_pos][1];
      FilterCoef[sub_pos][5] = CalculatedFilter[sub_pos][0];
    }
    else
    {
      for (i = 0; i < SQR_FILTER; i++)
        FilterCoef[sub_pos][i] = CalculatedFilter[sub_pos][i];  
    }
  }
  for (sub_pos = l_pos; sub_pos <= o_pos; sub_pos++)
  {
    for (i = 0; i < FILTER_SIZE; i++) 
    {
      FilterCoef[sub_pos][FILTER_SIZE-i-1] = FilterCoef[sub_pos-8][i];
    }
    FilterFlag[sub_pos] = FilterFlag[sub_pos-8]; 
  }

  // calculate UseAllPositions
  SubpelPositionsPattern =  (FilterFlag[k_pos] << 10) + (FilterFlag[j_pos] << 9) + (FilterFlag[i_pos] << 8) + 
    (FilterFlag[h_pos] << 7) + (FilterFlag[g_pos] << 6) + (FilterFlag[f_pos] << 5) + 
    (FilterFlag[e_pos] << 4) + (FilterFlag[d_pos] << 3) + (FilterFlag[c_pos] << 2) + 
    (FilterFlag[b_pos] << 1) + FilterFlag[a_pos];

  if (SubpelPositionsPattern == 2047)
    UseAllSubpelPositions = 1;        // for indicating that every sub-pel position is used
  else
    UseAllSubpelPositions = 0;

  PredictFilterCoefficientsVer ();  
  //  PrintFilterCoefDoubleSep(FilterCoef); 

  free_mem1Ddouble (y);
  free_mem2Ddouble (alpha);
  free_mem2Ddouble (beta);
}
// separable aif (END)

#ifdef E_DAIF
void quantize1Offset(double *offsetD, int *offsetQI, int *offsetQF)
{
  double offset; 
  int sign; 
  int factor; 
  int offsetInt; 
  int bitsForOffsetInt, bitsForOffsetFrac; 
  double offsetFrac; 
  int offsetFracQ; 
  int offsetIntCodeLen[] = 
  {
    1,  3,  3,  5,  5,  5,  5, 
    7,  7,  7,  7,  7,  7,  7,  7,
    9,  9,  9,  9,  9,  9,  9,  9,
    9,  9,  9,  9,  9,  9,  9,  9,
    11, 11, 11, 11, 11, 11, 11, 11,
    11, 11, 11, 11, 11, 11, 11, 11,
    11, 11, 11, 11, 11, 11, 11, 11,
    11, 11, 11, 11, 11, 11, 11, 11,
  };
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

  offset = *offsetD; 
  sign = offset >= 0. ? 1:-1; 
  offset = offset*sign; 
  // integer part
  offsetInt = (int)offset; 
  bitsForOffsetInt  = offsetIntCodeLen [offsetInt];
  // cap maximum integer offset
  if(offsetInt > 62) 
  {
    offsetInt = 62; 
  }

  bitsForOffsetFrac = offsetFracCodeLen[offsetInt]; 
  // send fraction if any bits left 
  if(bitsForOffsetFrac)
  {
    offsetFrac = offset-offsetInt; 
    factor = (1<<bitsForOffsetFrac);
    offsetFracQ = (int)(offsetFrac*factor+0.5); 
    if(offsetFracQ == (1<<bitsForOffsetFrac))
    {
      offsetInt++;
      offsetFracQ = 0;
    }
    *offsetQI = offsetInt; 
    *offsetQF = offsetFracQ;
    offset = (double)(offsetInt)+(double)offsetFracQ/(double)factor;
    *offsetD = sign*offset; 
  }
  else
  {
    *offsetQI = offsetInt;
    *offsetQF = -1;
    *offsetD = sign*offsetInt; 
  }
}

void QuantizeFilterOffsets()
{
  int pos; 
  for(pos = a_pos; pos <= o_pos; pos++)
  {
    quantize1Offset(&FilterCoef[pos][SQR(FILTER_SIZE)], &DiffQFilterOffsetI[pos], &DiffQFilterOffsetF[pos]); 
  }
  quantize1Offset(&FilterCoefInt[SQR(FILTER_SIZE_INT)], &DiffQFilterOffsetIntI, &DiffQFilterOffsetIntF); 
}

void QuantizeIntegerFilter()
{
  int i, j;
  double diff, is; 
  int diffInt, sign, numBits, factor; 

  int filterHalfLen = (FILTER_SIZE_INT/2)+1; 
  double filterIntFixed[] = 
  {
    41./2048.,  -72./2048.,  61./2048., 0., 0.,   
    -51./2048.,  -68./2048., 246./2048., 0., 0., 
    -20./2048.,  420./2048.,1229./2048., 0., 0.,
    0.,          0.,         0., 0., 0., 
    0.,          0.,         0., 0., 0., 
  };

  // predict top-left quadrant from fixed filter and quantize
  numBits = 8; factor = (1<<numBits);
  for(i = 0; i < filterHalfLen; i++)
  {
    for(j = 0; j < filterHalfLen; j++)
    {
      numBits = numQBitsInt[i*FILTER_SIZE_INT+j]-1; // 1-bit for sign
      factor = (1<<numBits); 
      diff = FilterCoefInt[i*FILTER_SIZE_INT+j]-filterIntFixed[i*FILTER_SIZE_INT+j];
      sign = (diff>0.) ? 1:-1; 
      diff = diff*sign; 
      diffInt = (int)(diff*(double)factor+0.5); 
      if(diffInt >= factor) 
        diffInt = factor-1;
      diffInt = diffInt*sign;
      DiffQFilterCoeffInt[i*FILTER_SIZE_INT+j] = diffInt;
      is = (double)diffInt/(double)factor; 
      FilterCoefInt[i*FILTER_SIZE_INT+j] = filterIntFixed[i*FILTER_SIZE_INT+j]+is;
    }
  }

  // predict top-right quadrant from top-left quadrant and quantize
  for(i = 0; i < filterHalfLen; i++)
  {
    for(j = filterHalfLen; j < FILTER_SIZE_INT; j++)
    {
      numBits = numQBitsInt[i*FILTER_SIZE_INT+j]-1; // 1-bit for sign
      factor = (1<<numBits); 
      diff = FilterCoefInt[i*FILTER_SIZE_INT+j]-FilterCoefInt[i*FILTER_SIZE_INT+(FILTER_SIZE_INT-1-j)];
      sign = (diff>0.) ? 1:-1; 
      diff = diff*sign; 
      diffInt = (int)(diff*(double)factor+0.5); 
      if(diffInt >= factor) 
        diffInt = factor-1;
      diffInt = diffInt*sign;
      DiffQFilterCoeffInt[i*FILTER_SIZE_INT+j] = diffInt; 
      is  = (double)diffInt/(double)factor; 
      FilterCoefInt[i*FILTER_SIZE_INT+j] = FilterCoefInt[i*FILTER_SIZE_INT+(FILTER_SIZE_INT-1-j)]+is;
    }
  }

  // predict bottom half from top half and quantize
  for(i = filterHalfLen; i < FILTER_SIZE_INT; i++)
  {
    for(j = 0; j < FILTER_SIZE_INT; j++)
    {
      numBits = numQBitsInt[i*FILTER_SIZE_INT+j]-1; // 1-bit for sign
      factor = (1<<numBits); 
      diff = FilterCoefInt[i*FILTER_SIZE_INT+j]-FilterCoefInt[(FILTER_SIZE_INT-1-i)*FILTER_SIZE_INT+j];
      sign = (diff>0.) ? 1:-1; 
      diff = diff*sign; 
      diffInt = (int)(diff*(double)factor+0.5); 
      if(diffInt >= factor) 
        diffInt = factor-1;
      diffInt = diffInt*sign;
      DiffQFilterCoeffInt[i*FILTER_SIZE_INT+j] = diffInt; 
      is = (double)diffInt/(double)factor; 
      FilterCoefInt[i*FILTER_SIZE_INT+j] = FilterCoefInt[(FILTER_SIZE_INT-1-i)*FILTER_SIZE_INT+j]+is;
    }
  }
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
#ifdef EDAIF2
void QuantizeFilterCoefficientsEDAIF2(double FCoef[SQR_FILTER], int QFCoef[SQR_FILTER], int sub_pos)
{
  //double number_of_steps = pow(2,NumberOfQBits-1);  // 1 bit for sign, NumberOfQBits-1 for amplitude
  double number_of_steps;
  int sign,maxMag;
  int j;
#ifdef E_DAIF
  for(j=0; j < SQR(FILTER_SIZE); j++)
#else
  for(j=0; j < SQR_FILTER; j++)
#endif
  {
    if((sub_pos == a_pos)||(sub_pos == b_pos)||(sub_pos == c_pos))
    {
      number_of_steps = pow(2,numQBits1DH[j]-1);  // 1 bit for sign, NumberOfQBits-1 for amplitude
    }
    else if((sub_pos == d_pos)||(sub_pos == h_pos)||(sub_pos == l_pos))
    {
      number_of_steps = pow(2,numQBits1DV[j]-1);  // 1 bit for sign, NumberOfQBits-1 for amplitude
    }
    else if((sub_pos == e_pos)||(sub_pos == g_pos)||(sub_pos == m_pos)||(sub_pos == o_pos))
    {
      number_of_steps = pow(2,numQBitsDIF[j]-1);  // 1 bit for sign, NumberOfQBits-1 for amplitude
    }
    else if((sub_pos == f2_pos)||(sub_pos == i2_pos)||(sub_pos == j2_pos)||(sub_pos == k2_pos) || (sub_pos == n2_pos))
    {
      number_of_steps = pow(2,numQBits2D[j]-1);  // 1 bit for sign, NumberOfQBits-1 for amplitude
    }
    else if(FilterFlag[sub_pos]<=1)// DAIF
    {
      number_of_steps = pow(2,numQBitsDIF[j]-1);  // 1 bit for sign, NumberOfQBits-1 for amplitude
    }
    else 
    {
      number_of_steps = pow(2,numQBits2D[j]-1);  // 1 bit for sign, NumberOfQBits-1 for amplitude
    }

    maxMag = (int)number_of_steps-1;
    sign = FCoef[j] >= 0.0 ? 1:-1;
    QFCoef[j] = (int)(fabs(FCoef[j])*number_of_steps + 0.5);		// step number
    QFCoef[j] = max(0, min(maxMag, QFCoef[j]));
    QFCoef[j] = sign*QFCoef[j];
    FCoef[j] = (double)(QFCoef[j])/number_of_steps; // reconstructed (quantized) filter coefficients
  }
}

#endif
void QuantizeFilterCoefficients(double FCoef[SQR_FILTER], int QFCoef[SQR_FILTER])
{
  double number_of_steps = pow(2,NumberOfQBits-1);  // 1 bit for sign, NumberOfQBits-1 for amplitude
  int sign;
  int j;
#ifdef E_DAIF
  for(j=0; j < SQR(FILTER_SIZE); j++)
#else
  for(j=0; j < SQR_FILTER; j++)
#endif
  {
    sign = FCoef[j] >= 0.0 ? 1:-1;
    QFCoef[j] = sign*(int)(fabs(FCoef[j])*number_of_steps + 0.5);    // step number
#ifdef DIRECTIONAL_FILTER
    if (input->ImpType == IMP_FLOAT32)
#endif
      FCoef[j] = (double)(QFCoef[j])/number_of_steps; // reconstructed (quantized) filter coefficients
  }
}

/*!
************************************************************************
* \brief
*    predict filter coefficients
*
************************************************************************
*/
void PredictFilterCoefficients(void)
{
  int i,j;
  // at first predict b and h position using inter prediction, b and h predicted values are already quantized in last frame
  for(i=0; i < SQR_FILTER; i++)
    PredFilterCoef[b_pos][i] = STANDARD_2D_FILTER[b_pos][i];

  for(i = 0; i < SQR_FILTER; i++)
    DiffFilterCoef[b_pos][i] = FilterCoef[b_pos][i] - PredFilterCoef[b_pos][i];
  QuantizeFilterCoefficients(DiffFilterCoef[b_pos], DiffQFilterCoef[b_pos]);

  for(i = 0; i < SQR_FILTER; i++)
    FilterCoef[b_pos][i] = DiffFilterCoef[b_pos][i] + PredFilterCoef[b_pos][i];

  ExtendFilterCoefficientsDouble(h_pos, FilterCoef);

  // predict a using inter-prediction

  for(i=0; i < SQR_FILTER; i++)
    PredFilterCoef[a_pos][i] = STANDARD_2D_FILTER[a_pos][i];

  for(i = 0; i < SQR_FILTER; i++)
    DiffFilterCoef[a_pos][i] = FilterCoef[a_pos][i] - PredFilterCoef[a_pos][i];
  QuantizeFilterCoefficients(DiffFilterCoef[a_pos], DiffQFilterCoef[a_pos]);
  for(i = 0; i < SQR_FILTER; i++)
    FilterCoef[a_pos][i] = DiffFilterCoef[a_pos][i] + PredFilterCoef[a_pos][i];

  ExtendFilterCoefficientsDouble(c_pos, FilterCoef);
  ExtendFilterCoefficientsDouble(d_pos, FilterCoef);
  ExtendFilterCoefficientsDouble(l_pos, FilterCoef);

  // predict e_pos from a_pos and d_pos by multiplication
  // predict f_pos from b_pos and d_pos by multiplication
  // predict j_pos from b_pos and h_pos by multiplication
  for(i = 0; i < FILTER_SIZE; i++)
  {
    for(j = 0; j < FILTER_SIZE; j++)
    {
      PredFilterCoef[e_pos][FILTER_SIZE*i+j] =
        FilterCoef[d_pos][FILTER_SIZE*i+FILTER_OFFSET]*FilterCoef[a_pos][FILTER_SIZE*FILTER_OFFSET+j];
      PredFilterCoef[f_pos][FILTER_SIZE*i+j] =
        FilterCoef[d_pos][FILTER_SIZE*i+FILTER_OFFSET]*FilterCoef[b_pos][FILTER_SIZE*FILTER_OFFSET+j];
      PredFilterCoef[j_pos][FILTER_SIZE*i+j] =
        FilterCoef[h_pos][FILTER_SIZE*i+FILTER_OFFSET]*FilterCoef[b_pos][FILTER_SIZE*FILTER_OFFSET+j];
    }
  }
  for(i = 0; i < SQR_FILTER; i++)
  {
    if(FilterFlag[e_pos])
      DiffFilterCoef[e_pos][i] = FilterCoef[e_pos][i] - PredFilterCoef[e_pos][i];
    if(FilterFlag[f_pos])
      DiffFilterCoef[f_pos][i] = FilterCoef[f_pos][i] - PredFilterCoef[f_pos][i];
    if(FilterFlag[j_pos])
      DiffFilterCoef[j_pos][i] = FilterCoef[j_pos][i] - PredFilterCoef[j_pos][i];
  }
  if(FilterFlag[e_pos])
    QuantizeFilterCoefficients(DiffFilterCoef[e_pos], DiffQFilterCoef[e_pos]);
  if(FilterFlag[f_pos])
    QuantizeFilterCoefficients(DiffFilterCoef[f_pos], DiffQFilterCoef[f_pos]);
  if(FilterFlag[j_pos])
    QuantizeFilterCoefficients(DiffFilterCoef[j_pos], DiffQFilterCoef[j_pos]);
  for(i = 0; i < SQR_FILTER; i++)
  {
    if(FilterFlag[e_pos])
      FilterCoef[e_pos][i] = DiffFilterCoef[e_pos][i] + PredFilterCoef[e_pos][i];
    if(FilterFlag[f_pos])
      FilterCoef[f_pos][i] = DiffFilterCoef[f_pos][i] + PredFilterCoef[f_pos][i];
    if(FilterFlag[j_pos])
      FilterCoef[j_pos][i] = DiffFilterCoef[j_pos][i] + PredFilterCoef[j_pos][i];
  }
  if(FilterFlag[e_pos])
  {
    ExtendFilterCoefficientsDouble(g_pos, FilterCoef);
    ExtendFilterCoefficientsDouble(m_pos, FilterCoef);
    ExtendFilterCoefficientsDouble(o_pos, FilterCoef);
  }
  if(FilterFlag[f_pos])
  {
    ExtendFilterCoefficientsDouble(i_pos, FilterCoef);
    ExtendFilterCoefficientsDouble(k_pos, FilterCoef);
    ExtendFilterCoefficientsDouble(n_pos, FilterCoef);
  }
}

#ifdef DIRECTIONAL_FILTER
/*! 
*************************************************************************************
* \brief
*   Extend directional 1D-AIF filter coefficients 
*  
* \para <title>
*    <paragraph>
*
* \para
*    <another paragraph>
*
* \author
*    - Dmytro Rusanovskyy                   <dmytro.rusanovskyy@tut.fi>
*************************************************************************************
*/
void ExtendDiagonalFilter(int sub_pos, double FCoef[15][SQR_FILTER])
{
  int i,j;

  if(SymmetryPosition[sub_pos] == 0)    //check if it is allowed to change original filter coefficients
  {
    if(IsDiagonal1D[sub_pos]==2)// NE-SW
    {
      for(i = 0; i < FILTER_SIZE; i++)
      {
        j = FILTER_SIZE-1-i;
        FCoef[sub_pos][j*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][j*FILTER_SIZE + i]];
      }
    }
    else if (IsDiagonal1D[sub_pos]==1)
    {
      for(i = 0; i < FILTER_SIZE; i++)
      {
        j = FILTER_SIZE-1-i;
        FCoef[sub_pos][i*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i*FILTER_SIZE + i]];
      }
    }
    else if (IsDiagonal1D[sub_pos]==3)
    {
      for(i = 0; i < FILTER_SIZE; i++)
      {
        FCoef[sub_pos][i*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][i*FILTER_SIZE + i]];
        j = FILTER_SIZE-1-i;
        FCoef[sub_pos][j*FILTER_SIZE + i] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][TwoDSymmetricPattern[sub_pos][j*FILTER_SIZE + i]];
      }
    }

#ifdef E_DAIF
    FCoef[sub_pos][SQR_FILTER-1] = FCoef[FILTER_NEW_SUB_POS[sub_pos]][SQR_FILTER-1];
#endif
  }
}


#ifdef EAIF

// quantize the filter coefficients for EAIF
void QuantizeFilterCoefficientsEAIF(double FCoef[SQR_FILTER], int QFCoef[SQR_FILTER], int numQBits[])
{
  double number_of_steps;  
  int sign;
  int j;
  int maxMag; 
  for(j=0; j < SQR(FILTER_SIZE); j++)
  {
    number_of_steps = pow(2,numQBits[j]-1); // 1 bit for sign, NumberOfQBits-1 for amplitude
    maxMag = (int)number_of_steps-1;
    sign = FCoef[j] >= 0.0 ? 1:-1;
    QFCoef[j] = (int)(fabs(FCoef[j])*number_of_steps + 0.5);		// step number
    QFCoef[j] = max(0, min(maxMag, QFCoef[j]));
    QFCoef[j] = sign*QFCoef[j];
    FCoef[j] = (double)(QFCoef[j])/number_of_steps; // reconstructed (quantized) filter coefficients
  }
}

// predict the filter coeffiicents for EAIF
void PredictFilterCoefficientsEAIF(void)
{
  int i;

  int *numQBits; 

  int sub_pos; 

  for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
  {
    if(!SymmetryPosition[sub_pos])
      continue; 

    if(!FilterFlag[sub_pos])
      continue; 

    if(Is1DPosition_orig[sub_pos])
    {
      for(i=0; i < SQR(FILTER_SIZE); i++)
        PredFilterCoef[sub_pos][i] = STANDARD_2D_FILTER_orig[sub_pos][i]; 
      // offset
      PredFilterCoef[sub_pos][i] = 0.; 

      for(i = 0; i < SQR_FILTER; i++)
        DiffFilterCoef[sub_pos][i] = FilterCoef[sub_pos][i] - PredFilterCoef[sub_pos][i];

      if((sub_pos+1)/4==0) // horizontal 
        numQBits = numQBits1DH;
      else // if((sub_pos+1)%4==0)
        numQBits = numQBits1DV; // vertical 
      QuantizeFilterCoefficientsEAIF(DiffFilterCoef[sub_pos], DiffQFilterCoef[sub_pos], numQBits);

      for(i = 0; i < SQR_FILTER; i++)
        FilterCoef[sub_pos][i] = DiffFilterCoef[sub_pos][i] + PredFilterCoef[sub_pos][i];
    }
    else
    {
      for(i=0; i < SQR(FILTER_SIZE); i++)
        PredFilterCoef[sub_pos][i] = 0.; 
      // offset
      PredFilterCoef[sub_pos][i] = 0.; 

      for(i = 0; i < SQR_FILTER; i++)
        DiffFilterCoef[sub_pos][i] = FilterCoef[sub_pos][i] - PredFilterCoef[sub_pos][i];

      numQBits = numQBits2D;
      QuantizeFilterCoefficientsEAIF(DiffFilterCoef[sub_pos], DiffQFilterCoef[sub_pos], numQBits);

      for(i = 0; i < SQR_FILTER; i++)
        FilterCoef[sub_pos][i] = DiffFilterCoef[sub_pos][i] + PredFilterCoef[sub_pos][i];
    }

  }

  flipHoriz(FilterCoef[a_pos], FilterCoef[c_pos]);
  flipVert (FilterCoef[d_pos], FilterCoef[l_pos]);
  flipHoriz(FilterCoef[e_pos], FilterCoef[g_pos]);
  flipHoriz(FilterCoef[i_pos], FilterCoef[k_pos]);
  flipVert (FilterCoef[e_pos], FilterCoef[m_pos]);
  flipVert (FilterCoef[f_pos], FilterCoef[n_pos]);
  flipHoriz(FilterCoef[m_pos], FilterCoef[o_pos]);
}
#endif

/*! 
*************************************************************************************
* \brief
*   Predict 1D-AIF from static diagonal filter
*  
* \para <title>
*    <paragraph>
*
* \para
*    <another paragraph>
*
* \author
*    - Dmytro Rusanovskyy                   <dmytro.rusanovskyy@tut.fi>
*************************************************************************************
*/

void PredictFilterCoefficients1DAIF(void)
{
  int i;

  // at first predict b and h position using inter prediction, b and h predicted values are already quantized in last frame
  for(i=0; i < SQR_FILTER; i++)
    PredFilterCoef[b_pos][i] = STANDARD_2D_FILTER[b_pos][i];

  for(i = 0; i < SQR_FILTER; i++)
    DiffFilterCoef[b_pos][i] = FilterCoef[b_pos][i] - PredFilterCoef[b_pos][i];

  QuantizeFilterCoefficients(DiffFilterCoef[b_pos], DiffQFilterCoef[b_pos]);
  if (input->ImpType==IMP_FLOAT32)
    for(i = 0; i < SQR_FILTER; i++)
      FilterCoef[b_pos][i] = DiffFilterCoef[b_pos][i] + PredFilterCoef[b_pos][i];

  ExtendFilterCoefficientsDouble(h_pos, FilterCoef);

  // predict 1D filter "a"  using inter-prediction

  for(i=0; i < SQR_FILTER; i++)
    PredFilterCoef[a_pos][i] = STANDARD_2D_FILTER[a_pos][i];

  for(i = 0; i < SQR_FILTER; i++)
    DiffFilterCoef[a_pos][i] = FilterCoef[a_pos][i] - PredFilterCoef[a_pos][i];

  QuantizeFilterCoefficients(DiffFilterCoef[a_pos], DiffQFilterCoef[a_pos]);
  if (input->ImpType==IMP_FLOAT32)
    for(i = 0; i < SQR_FILTER; i++)
      FilterCoef[a_pos][i] = DiffFilterCoef[a_pos][i] + PredFilterCoef[a_pos][i];

  ExtendFilterCoefficientsDouble(c_pos, FilterCoef);
  ExtendFilterCoefficientsDouble(d_pos, FilterCoef);
  ExtendFilterCoefficientsDouble(l_pos, FilterCoef);

  /*-----------------  Diagonal filters  -----------------*/
  for(i=0; i < SQR_FILTER; i++)
  {
    PredFilterCoef[e_pos][i] = STANDARD_2D_FILTER[e_pos][i];
    PredFilterCoef[f_pos][i] = STANDARD_2D_FILTER[f_pos][i];
    PredFilterCoef[j_pos][i] = STANDARD_2D_FILTER[j_pos][i];

    DiffFilterCoef[e_pos][i] = FilterCoef[e_pos][i] - PredFilterCoef[e_pos][i];
    DiffFilterCoef[f_pos][i] = FilterCoef[f_pos][i] - PredFilterCoef[f_pos][i];
    DiffFilterCoef[j_pos][i] = FilterCoef[j_pos][i] - PredFilterCoef[j_pos][i];
  }

  QuantizeFilterCoefficients(DiffFilterCoef[e_pos], DiffQFilterCoef[e_pos]);
  QuantizeFilterCoefficients(DiffFilterCoef[f_pos], DiffQFilterCoef[f_pos]);
  QuantizeFilterCoefficients(DiffFilterCoef[j_pos], DiffQFilterCoef[j_pos]);

  if (input->ImpType==IMP_FLOAT32)  
    for(i = 0; i < SQR_FILTER; i++)
    {
      FilterCoef[e_pos][i] = DiffFilterCoef[e_pos][i] + PredFilterCoef[e_pos][i];
      FilterCoef[f_pos][i] = DiffFilterCoef[f_pos][i] + PredFilterCoef[f_pos][i];
      FilterCoef[j_pos][i] = DiffFilterCoef[j_pos][i] + PredFilterCoef[j_pos][i];

      FilterCoef[g_pos][i] = 0.0;
      FilterCoef[i_pos][i] = 0.0;
      FilterCoef[k_pos][i] = 0.0;
      FilterCoef[m_pos][i] = 0.0;
      FilterCoef[n_pos][i] = 0.0;
      FilterCoef[o_pos][i] = 0.0;
    }
    ExtendDiagonalFilter(g_pos, FilterCoef);
    ExtendDiagonalFilter(i_pos, FilterCoef);
    ExtendDiagonalFilter(k_pos, FilterCoef);
    ExtendDiagonalFilter(m_pos, FilterCoef);
    ExtendDiagonalFilter(n_pos, FilterCoef);
    ExtendDiagonalFilter(o_pos, FilterCoef);
}
#endif
// separable aif (BEGIN)
/*!
************************************************************************
* \brief
*    predict filter coefficients for horizontal filter
*
************************************************************************
*/
void PredictFilterCoefficientsHor (void)      
{
  int i, pos;

  for (pos = a_pos; pos <= c_pos; pos++)
  {
    for (i = 0; i < FILTER_SIZE; i++)
      DiffFilterCoef[pos][i] = FilterCoef[pos][i] - STANDARD_2D_FILTER_SEP[pos][i];
    QuantizeFilterCoefficients (DiffFilterCoef[pos], DiffQFilterCoef[pos]);
    for (i = 0; i < FILTER_SIZE; i++)
      FilterCoef[pos][i] = DiffFilterCoef[pos][i] + STANDARD_2D_FILTER_SEP[pos][i];
  }
}

/*!
************************************************************************
* \brief
*    predict filter coefficients for vertical filter
*
************************************************************************
*/
void PredictFilterCoefficientsVer (void)     
{
  int i, pos;

  for (pos = d_pos; pos <= o_pos; pos++)
  {
    for (i = 0; i < FILTER_SIZE; i++)
      DiffFilterCoef[pos][i] = FilterCoef[pos][i] - STANDARD_2D_FILTER_SEP[pos][i];
    QuantizeFilterCoefficients (DiffFilterCoef[pos], DiffQFilterCoef[pos]);
    for (i = 0; i < FILTER_SIZE; i++)
      FilterCoef[pos][i] = DiffFilterCoef[pos][i] + STANDARD_2D_FILTER_SEP[pos][i];
  }
}
// separable aif (END)

/*!
************************************************************************
* \brief
*    Upsample 4 times with new calculated filter coefficients.
*    Store the old upsampled frames.
************************************************************************
*/
void UnifiedOneForthPixWithNewFilter(void)
{

  int i, j,ii, jj;
  unsigned frame;
  imgpel **out4Y_aif, **out4Y;
  imgpel **refY;
  int img_width, img_height;
  int sub_pos, pos_x, pos_y, x_sub, y_sub;
  double is;
#ifdef E_DAIF
  imgpel *ref11_aif;
#endif

  if(current_poc == img->ThisPOC)
    return;// avoids multiple interolating of the same frames
  else
    current_poc = img->ThisPOC;
  for(frame = 0; frame < dpb.ref_frames_in_buffer; frame++)
  {
    img_width  = dpb.fs_ref[frame]->frame->size_x;
    img_height = dpb.fs_ref[frame]->frame->size_y;
    out4Y_aif  = dpb.fs_ref[frame]->frame->imgY_ups_aif;
    out4Y      = dpb.fs_ref[frame]->frame->imgY_ups;
    refY       = dpb.fs_ref[frame]->frame->imgY;
#ifdef E_DAIF
    ref11_aif = dpb.fs_ref[frame]->frame->imgY_11_aif;
#endif

    for(i=-4*IMG_PAD_SIZE; i < (IMG_PAD_SIZE+img_height)*4; i++)
      for(j=-4*IMG_PAD_SIZE; j < (IMG_PAD_SIZE+img_width)*4; j++)
      {
        x_sub = (j+4*IMG_PAD_SIZE)%4;       // x-sub-coordinate in a 4x4block
        y_sub = (i+4*IMG_PAD_SIZE)%4;       // y-sub-coordinate in a 4x4block
        sub_pos = x_sub + 4*y_sub;          // pos 1..15 in a 4x4 block
        is = 0.0;
        if(sub_pos)
        {
          for(ii = 0; ii < FILTER_SIZE; ii++)
            for(jj = 0; jj < FILTER_SIZE; jj++)
            {
              pos_y = max(0, min(img_height-1, i/4-FILTER_OFFSET+ii ) );
              if(i<0 && pos_y > 0)
                pos_y -=1;
              pos_x = max(0, min(img_width -1, j/4-FILTER_OFFSET+jj ) );
              if(j<0 && pos_x > 0)
                pos_x -=1;
              is += (FilterCoef[sub_pos-1][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
              //is += (STANDARD_2D_FILTER[sub_pos-1][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
            }
#ifdef E_DAIF
            is += FilterCoef[sub_pos-1][SQR(FILTER_SIZE)];
#endif
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
            out4Y_aif[i+4*IMG_PAD_SIZE][j+4*IMG_PAD_SIZE] = max(0,min(((1<<img->bitdepth_luma)-1),(int)(is+0.5)));
#else
            out4Y_aif[i+4*IMG_PAD_SIZE][j+4*IMG_PAD_SIZE] = max(0,min(255,(int)(is+0.5)));
#endif
        }
        else
        {
#ifdef E_DAIF
          if(FilterFlagInt)
          {
            for(ii = 0; ii < FILTER_SIZE_INT; ii++)
              for(jj = 0; jj < FILTER_SIZE_INT; jj++)
              {
                pos_y = max(0, min(img_height-1, i/4-FILTER_OFFSET_INT+ii ) );
                pos_x = max(0, min(img_width -1, j/4-FILTER_OFFSET_INT+jj ) );
                is += (FilterCoefInt[FILTER_SIZE_INT*ii+jj] * refY[pos_y][pos_x]);
              }

              is += FilterCoefInt[SQR(FILTER_SIZE_INT)];

#ifdef INTERNAL_BIT_DEPTH_INCREASE
              out4Y_aif[i+4*IMG_PAD_SIZE][j+4*IMG_PAD_SIZE] = max(0,min(((1<<img->bitdepth_luma)-1),(int)(is+0.5)));
#else
              out4Y_aif[i+4*IMG_PAD_SIZE][j+4*IMG_PAD_SIZE] = max(0,min(255,(int)(is+0.5)));
#endif  // INTERNAL_BIT_DEPTH_INCREASE
          }
          else
          {
            out4Y_aif[i+4*IMG_PAD_SIZE][j+4*IMG_PAD_SIZE] =
              refY[max(0,min(img_height-1,i/4))][max(0,min(img_width-1,j/4))];
          }
#else
          {
            out4Y_aif[i+4*IMG_PAD_SIZE][j+4*IMG_PAD_SIZE] =
              refY[max(0,min(img_height-1,i/4))][max(0,min(img_width-1,j/4))];
          }
#endif  // E_DAIF
        }
      }

#ifdef E_DAIF
      // Generate 1/1th pel representation (used for integer pel MV search)
      {
        int x, y, yy , y_pos;

        for (y = 0; y < img_height; y++)
        {
          yy = (y + IMG_PAD_SIZE)<<2;
          y_pos = y * img_width;
          for (x = 0; x < img_width; x++)
            ref11_aif[y_pos + x] = out4Y_aif[yy][(x + IMG_PAD_SIZE)<<2];
        }
      }
#endif  // E_DAIF
  }
  //  FindSAD();
  //  PrintFilterCoefFloatFile(FilterCoef);
  //save2file_images_QPEL(0);
}
#ifdef DIRECTIONAL_FILTER
/*! 
*************************************************************************************
* \brief
*   Up-sample reference frame(s) with directional (diagonal) adaptive interpolation filter 1D-AIF
*  
*
* \param <parameter>
*    No input parameters
*
* \return
*    No output return
*
* \note
*   1. Uses adaptive filters stored in FilterCoef[15][SQR_FILTER]
*   2. UnifiedOneForthPixWith_1DAIF_float is a faster alternative to UnifiedOneForthPixWithNewFilter function if 1D-AIF is used for itnerpolation. 
*   3. UnifiedOneForthPixWith_1DAIF_float and UnifiedOneForthPixWithNewFilter provide identical interpolated frame(s). 
*
* \author
*    - Dmytro Rusanovskyy                   <dmytro.rusanovskyy@tut.fi>
*************************************************************************************
*/

void UnifiedOneForthPixWith_1DAIF_float(void)
{

  int i, j,ii, jj;
  unsigned frame;
  imgpel **out4Y_aif, **out4Y;
  imgpel **refY;
  int img_width, img_height;
  int ref_pos_x,ref_pos_y;
  int ref_pos_x1,ref_pos_y1;
  int ref_pos_x2,ref_pos_y2;
  double is, is0, is1,is2,is3;
  imgpel or_pel;
  if(current_poc == img->ThisPOC)
    return;// avoids multiple interolating of the same frames
  else
    current_poc = img->ThisPOC;

  for(frame = 0; frame < dpb.ref_frames_in_buffer; frame++)
  {
    img_width  = dpb.fs_ref[frame]->frame->size_x;
    img_height = dpb.fs_ref[frame]->frame->size_y;
    out4Y_aif  = dpb.fs_ref[frame]->frame->imgY_ups_aif;
    out4Y      = dpb.fs_ref[frame]->frame->imgY_ups;
    refY       = dpb.fs_ref[frame]->frame->imgY;
    for(i=0; i < (2*IMG_PAD_SIZE+img_height)*4; i+=4)
      for(j=0; j < (2*IMG_PAD_SIZE+img_width)*4; j+=4)
      {
        //  Copy Integer-pel sample
        ref_pos_x = j/4- IMG_PAD_SIZE;
        ref_pos_y = i/4- IMG_PAD_SIZE;
        ref_pos_x = max(0, min(img_width -1, ref_pos_x));
        ref_pos_y = max(0, min(img_height -1, ref_pos_y));
        out4Y_aif[i][j] = refY[ref_pos_y][ref_pos_x];

        //  a,b,c
        is = is0 = is1 = is2 = 0;
        ii = FILTER_SIZE*FILTER_OFFSET;
        ref_pos_y = i/4 - IMG_PAD_SIZE;
        ref_pos_y = max(0, min(img_height -1, ref_pos_y));
        for(jj = 0; jj < FILTER_SIZE; jj++)  // horizontal.
        {
          ref_pos_x = j/4-FILTER_OFFSET + jj - IMG_PAD_SIZE;
          ref_pos_x = max(0, min(img_width -1, ref_pos_x));
          or_pel = (imgpel)refY[ref_pos_y][ref_pos_x];
          is0 += (FilterCoef[a_pos][ii+jj]*or_pel);
          is1 += (FilterCoef[b_pos][ii+jj]*or_pel);
          is2 += (FilterCoef[c_pos][ii+jj]*or_pel);
        }
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y_aif[i][j+1] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is0+0.5)));
        out4Y_aif[i][j+2] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is1+0.5)));
        out4Y_aif[i][j+3] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is2+0.5)));
#else
        out4Y_aif[i][j+1] = (imgpel)max(0,min(255,(int)(is0+0.5)));
        out4Y_aif[i][j+2] = (imgpel)max(0,min(255,(int)(is1+0.5)));
        out4Y_aif[i][j+3] = (imgpel)max(0,min(255,(int)(is2+0.5)));
#endif

        //  d,h,l
        is = is0 = is1 = is2 = 0;
        jj = FILTER_OFFSET;
        ref_pos_x = j/4  - IMG_PAD_SIZE;
        ref_pos_x = max(0, min(img_width -1, ref_pos_x));

        for(ii = 0; ii < FILTER_SIZE; ii++)  // horizontal.
        {
          ref_pos_y = i/4-FILTER_OFFSET + ii - IMG_PAD_SIZE;
          ref_pos_y = max(0, min(img_height -1, ref_pos_y));
          or_pel = (imgpel)refY[ref_pos_y][ref_pos_x];
          is0 += (FilterCoef[d_pos][ii*FILTER_SIZE+jj]*or_pel);
          is1 += (FilterCoef[h_pos][ii*FILTER_SIZE+jj]*or_pel);
          is2 += (FilterCoef[l_pos][ii*FILTER_SIZE+jj]*or_pel);

        }
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y_aif[i+1][j] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is0+0.5)));
        out4Y_aif[i+2][j] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is1+0.5)));
        out4Y_aif[i+3][j] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is2+0.5)));
#else
        out4Y_aif[i+1][j] = (imgpel)max(0,min(255,(int)(is0+0.5)));
        out4Y_aif[i+2][j] = (imgpel)max(0,min(255,(int)(is1+0.5)));
        out4Y_aif[i+3][j] = (imgpel)max(0,min(255,(int)(is2+0.5)));
#endif

        //j 
        is = 0;
        ref_pos_x = j/4-FILTER_OFFSET- IMG_PAD_SIZE;
        ref_pos_y = i/4-FILTER_OFFSET- IMG_PAD_SIZE;

        ref_pos_y1 = max(0, min(img_height -1, ref_pos_y));
        ref_pos_x1 = max(0, min(img_width -1, ref_pos_x));
        ref_pos_y2 = max(0, min(img_height -1, ref_pos_y+5));
        ref_pos_x2 = max(0, min(img_width -1, ref_pos_x+5 ));

        is += (refY[ref_pos_y1][ref_pos_x1]+ refY[ref_pos_y2][ref_pos_x1] 
        + refY[ref_pos_y1][ref_pos_x2]+ refY[ref_pos_y2][ref_pos_x2]) 
          * FilterCoef[j_pos][0];

        ref_pos_x = j/4-FILTER_OFFSET + 1- IMG_PAD_SIZE;
        ref_pos_y = i/4-FILTER_OFFSET + 1- IMG_PAD_SIZE;
        ref_pos_y1 = max(0, min(img_height -1, ref_pos_y));
        ref_pos_x1 = max(0, min(img_width -1, ref_pos_x));
        ref_pos_y2 = max(0, min(img_height -1, ref_pos_y+3));
        ref_pos_x2 = max(0, min(img_width -1, ref_pos_x+3 ));

        is += (refY[ref_pos_y1][ref_pos_x1]+ refY[ref_pos_y2][ref_pos_x1] 
        + refY[ref_pos_y1][ref_pos_x2]+ refY[ref_pos_y2][ref_pos_x2]) 
          * FilterCoef[j_pos][FILTER_SIZE+1];

        ref_pos_x = j/4-FILTER_OFFSET + 2- IMG_PAD_SIZE;
        ref_pos_y = i/4-FILTER_OFFSET + 2- IMG_PAD_SIZE;
        ref_pos_y1 = max(0, min(img_height -1, ref_pos_y));
        ref_pos_x1 = max(0, min(img_width -1, ref_pos_x));
        ref_pos_y2 = max(0, min(img_height -1, ref_pos_y+1));
        ref_pos_x2 = max(0, min(img_width -1, ref_pos_x+1 ));

        is += (refY[ref_pos_y1][ref_pos_x1]+ refY[ref_pos_y2][ref_pos_x1] 
        + refY[ref_pos_y1][ref_pos_x2]+ refY[ref_pos_y2][ref_pos_x2]) 
          * FilterCoef[j_pos][FILTER_SIZE*2+2];
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y_aif[i+2][j+2] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is+0.5)));
#else
        out4Y_aif[i+2][j+2] = (imgpel)max(0,min(255,(int)(is+0.5)));
#endif


        //e,g,m,o
        is0 = is1 = is2 = is3 = 0;
        ref_pos_x = j/4-FILTER_OFFSET- IMG_PAD_SIZE;
        ref_pos_y = i/4-FILTER_OFFSET- IMG_PAD_SIZE;
        for(ii = 0; ii < FILTER_SIZE; ii++)  // horizontal.
        {
          jj = ii;// filter NW-SE direction
          ref_pos_y1 = max(0, min(img_height -1, ref_pos_y+ii));
          ref_pos_x1 = max(0, min(img_width -1, ref_pos_x + ii));
          is0 += FilterCoef[e_pos][ii*FILTER_SIZE+jj]*refY[ref_pos_y1][ref_pos_x1]; //e
          is1 += FilterCoef[o_pos][ii*FILTER_SIZE+jj]*refY[ref_pos_y1][ref_pos_x1]; //o

          jj = FILTER_SIZE-1 - ii;// filter NE-SW direction
          ref_pos_x2 = max(0, min(img_width -1, ref_pos_x + jj ));
          is2 += FilterCoef[g_pos][ii*FILTER_SIZE+jj]*refY[ref_pos_y1][ref_pos_x2]; 
          is3 += FilterCoef[m_pos][ii*FILTER_SIZE+jj]*refY[ref_pos_y1][ref_pos_x2]; 
        }
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y_aif[i+1][j+1] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is0+0.5)));
        out4Y_aif[i+3][j+3] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is1+0.5)));
        out4Y_aif[i+1][j+3] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is2+0.5)));
        out4Y_aif[i+3][j+1] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is3+0.5)));
#else
        out4Y_aif[i+1][j+1] = (imgpel)max(0,min(255,(int)(is0+0.5)));
        out4Y_aif[i+3][j+3] = (imgpel)max(0,min(255,(int)(is1+0.5)));
        out4Y_aif[i+1][j+3] = (imgpel)max(0,min(255,(int)(is2+0.5)));
        out4Y_aif[i+3][j+1] = (imgpel)max(0,min(255,(int)(is3+0.5)));
#endif

        //  f,n,i,k
        is0 = is1 = is2 = is3 = 0;
        //  set NW 0 positions.
        ref_pos_x = j/4-FILTER_OFFSET- IMG_PAD_SIZE;
        ref_pos_y = i/4-FILTER_OFFSET- IMG_PAD_SIZE;

        for(ii = 0; ii < FILTER_SIZE; ii++)  // horizontal.
        {
          jj = FILTER_SIZE-1 - ii;

          ref_pos_y1 = max(0, min(img_height -1, ref_pos_y+ii));
          ref_pos_y2 = max(0, min(img_height -1, ref_pos_y+jj));
          ref_pos_x1 = max(0, min(img_width -1, ref_pos_x + ii));
          ref_pos_x2 = max(0, min(img_width -1, ref_pos_x + jj ));

          is0 += FilterCoef[f_pos][ii*FILTER_SIZE+ii]*(refY[ref_pos_y1][ref_pos_x1] + refY[ref_pos_y1][ref_pos_x2]);//f
          is1 += FilterCoef[n_pos][ii*FILTER_SIZE+ii]*(refY[ref_pos_y1][ref_pos_x1] + refY[ref_pos_y1][ref_pos_x2]);//n
          is2 += FilterCoef[i_pos][ii*FILTER_SIZE+ii]*(refY[ref_pos_y1][ref_pos_x1] + refY[ref_pos_y2][ref_pos_x1]);//f
          is3 += FilterCoef[k_pos][ii*FILTER_SIZE+ii]*(refY[ref_pos_y1][ref_pos_x1] + refY[ref_pos_y2][ref_pos_x1]);//k
        }
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y_aif[i+1][j+2] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is0+0.5)));
        out4Y_aif[i+3][j+2] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is1+0.5)));
        out4Y_aif[i+2][j+1] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is2+0.5)));
        out4Y_aif[i+2][j+3] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is3+0.5)));
#else
        out4Y_aif[i+1][j+2] = (imgpel)max(0,min(255,(int)(is0+0.5)));
        out4Y_aif[i+3][j+2] = (imgpel)max(0,min(255,(int)(is1+0.5)));
        out4Y_aif[i+2][j+1] = (imgpel)max(0,min(255,(int)(is2+0.5)));
        out4Y_aif[i+2][j+3] = (imgpel)max(0,min(255,(int)(is3+0.5)));
#endif
      }
  }
}
#endif

#ifdef E_DAIF
/*! 
*************************************************************************************
* \brief
*   Up-sample reference frame(s) with ENHANCED directional adaptive interpolation filter (EDAIF)
*  
*
* \param <parameter>
*    No input parameters
*
* \return
*    No output return
*
* \note
*   1. Uses adaptive filters stored in FilterCoef[15][SQR_FILTER]
*   2. UnifiedOneForthPixWith_EDAIF is a faster alternative to UnifiedOneForthPixWithNewFilter function if EDAIF is used for interpolation. 
*   3. UnifiedOneForthPixWith_EDAIF and UnifiedOneForthPixWithNewFilter provide identical interpolated frame(s). 
*
* \author
*    - Yan Ye                   <yye@qualcomm.com>
*************************************************************************************
*/

void UnifiedOneForthPixWith_EDAIF(void)
{

  int i, j,ii, jj;
  unsigned frame;
  imgpel **out4Y_aif, **out4Y;
  imgpel **refY;
  int img_width, img_height;
  int ref_pos_x,ref_pos_y;
  int ref_pos_x1,ref_pos_y1;
  int ref_pos_x2,ref_pos_y2;
  double is, is0, is1,is2,is3;
  imgpel or_pel;
  imgpel *ref11_aif;

  if(current_poc == img->ThisPOC)
    return;// avoids multiple interolating of the same frames
  else
    current_poc = img->ThisPOC;

  for(frame = 0; frame < dpb.ref_frames_in_buffer; frame++)
  {
    img_width  = dpb.fs_ref[frame]->frame->size_x;
    img_height = dpb.fs_ref[frame]->frame->size_y;
    out4Y_aif  = dpb.fs_ref[frame]->frame->imgY_ups_aif;
    out4Y      = dpb.fs_ref[frame]->frame->imgY_ups;
    refY       = dpb.fs_ref[frame]->frame->imgY;
    ref11_aif = dpb.fs_ref[frame]->frame->imgY_11_aif;

    for(i=0; i < (2*IMG_PAD_SIZE+img_height)*4; i+=4)
      for(j=0; j < (2*IMG_PAD_SIZE+img_width)*4; j+=4)
      {
        if(!FilterFlagInt)
        {
          ref_pos_x = j/4- IMG_PAD_SIZE;
          ref_pos_y = i/4- IMG_PAD_SIZE;
          ref_pos_x = max(0, min(img_width -1, ref_pos_x));
          ref_pos_y = max(0, min(img_height -1, ref_pos_y));
          out4Y_aif[i][j] = refY[ref_pos_y][ref_pos_x];
        }
        else
        {
          is = 0.;
          for(ii = 0; ii < FILTER_SIZE_INT; ii++)
            for(jj = 0; jj < FILTER_SIZE_INT; jj++)
            {
              ref_pos_x = max(0, min(img_width -1, j/4-IMG_PAD_SIZE-FILTER_OFFSET_INT+jj ) );
              ref_pos_y = max(0, min(img_height-1, i/4-IMG_PAD_SIZE-FILTER_OFFSET_INT+ii ) );
              is += (FilterCoefInt[FILTER_SIZE_INT*ii+jj] * refY[ref_pos_y][ref_pos_x]);
            }

            is += FilterCoefInt[SQR(FILTER_SIZE_INT)];

#ifdef INTERNAL_BIT_DEPTH_INCREASE
            out4Y_aif[i][j] = max(0,min(((1<<img->bitdepth_luma)-1),(int)(is+0.5)));
#else
            out4Y_aif[i][j] = max(0,min(255,(int)(is+0.5)));
#endif  // INTERNAL_BIT_DEPTH_INCREASE
        }

        //  a,b,c
        is = is0 = is1 = is2 = 0;
        ii = FILTER_SIZE*FILTER_OFFSET;
        ref_pos_y = i/4 - IMG_PAD_SIZE;
        ref_pos_y = max(0, min(img_height -1, ref_pos_y));
        for(jj = 0; jj < FILTER_SIZE; jj++)  // horizontal.
        {
          ref_pos_x = j/4-FILTER_OFFSET + jj - IMG_PAD_SIZE;
          ref_pos_x = max(0, min(img_width -1, ref_pos_x));
          or_pel = (imgpel)refY[ref_pos_y][ref_pos_x];
          is0 += (FilterCoef[a_pos][ii+jj]*or_pel);
          is1 += (FilterCoef[b_pos][ii+jj]*or_pel);
          is2 += (FilterCoef[c_pos][ii+jj]*or_pel);
        }
        is0 += FilterCoef[a_pos][SQR(FILTER_SIZE)];
        is1 += FilterCoef[b_pos][SQR(FILTER_SIZE)];
        is2 += FilterCoef[c_pos][SQR(FILTER_SIZE)];
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y_aif[i][j+1] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is0+0.5)));
        out4Y_aif[i][j+2] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is1+0.5)));
        out4Y_aif[i][j+3] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is2+0.5)));
#else
        out4Y_aif[i][j+1] = (imgpel)max(0,min(255,(int)(is0+0.5)));
        out4Y_aif[i][j+2] = (imgpel)max(0,min(255,(int)(is1+0.5)));
        out4Y_aif[i][j+3] = (imgpel)max(0,min(255,(int)(is2+0.5)));
#endif

        //  d,h,l
        is = is0 = is1 = is2 = 0;
        jj = FILTER_OFFSET;
        ref_pos_x = j/4  - IMG_PAD_SIZE;
        ref_pos_x = max(0, min(img_width -1, ref_pos_x));

        for(ii = 0; ii < FILTER_SIZE; ii++)  // horizontal.
        {
          ref_pos_y = i/4-FILTER_OFFSET + ii - IMG_PAD_SIZE;
          ref_pos_y = max(0, min(img_height -1, ref_pos_y));
          or_pel = (imgpel)refY[ref_pos_y][ref_pos_x];
          is0 += (FilterCoef[d_pos][ii*FILTER_SIZE+jj]*or_pel);
          is1 += (FilterCoef[h_pos][ii*FILTER_SIZE+jj]*or_pel);
          is2 += (FilterCoef[l_pos][ii*FILTER_SIZE+jj]*or_pel);
        }
        is0 += FilterCoef[d_pos][SQR(FILTER_SIZE)];
        is1 += FilterCoef[h_pos][SQR(FILTER_SIZE)];
        is2 += FilterCoef[l_pos][SQR(FILTER_SIZE)];
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y_aif[i+1][j] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is0+0.5)));
        out4Y_aif[i+2][j] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is1+0.5)));
        out4Y_aif[i+3][j] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is2+0.5)));
#else
        out4Y_aif[i+1][j] = (imgpel)max(0,min(255,(int)(is0+0.5)));
        out4Y_aif[i+2][j] = (imgpel)max(0,min(255,(int)(is1+0.5)));
        out4Y_aif[i+3][j] = (imgpel)max(0,min(255,(int)(is2+0.5)));
#endif

        //j 
        if(FilterFlag[j_pos])
        {
          is = 0;
          ref_pos_x = j/4-FILTER_OFFSET- IMG_PAD_SIZE;
          ref_pos_y = i/4-FILTER_OFFSET- IMG_PAD_SIZE;

          ref_pos_y1 = max(0, min(img_height -1, ref_pos_y));
          ref_pos_x1 = max(0, min(img_width -1, ref_pos_x));
          ref_pos_y2 = max(0, min(img_height -1, ref_pos_y+5));
          ref_pos_x2 = max(0, min(img_width -1, ref_pos_x+5 ));

          is += (refY[ref_pos_y1][ref_pos_x1]+ refY[ref_pos_y2][ref_pos_x1] 
          + refY[ref_pos_y1][ref_pos_x2]+ refY[ref_pos_y2][ref_pos_x2]) 
            * FilterCoef[j_pos][0];

          ref_pos_x = j/4-FILTER_OFFSET + 1- IMG_PAD_SIZE;
          ref_pos_y = i/4-FILTER_OFFSET + 1- IMG_PAD_SIZE;
          ref_pos_y1 = max(0, min(img_height -1, ref_pos_y));
          ref_pos_x1 = max(0, min(img_width -1, ref_pos_x));
          ref_pos_y2 = max(0, min(img_height -1, ref_pos_y+3));
          ref_pos_x2 = max(0, min(img_width -1, ref_pos_x+3 ));

          is += (refY[ref_pos_y1][ref_pos_x1]+ refY[ref_pos_y2][ref_pos_x1] 
          + refY[ref_pos_y1][ref_pos_x2]+ refY[ref_pos_y2][ref_pos_x2]) 
            * FilterCoef[j_pos][FILTER_SIZE+1];

          ref_pos_x = j/4-FILTER_OFFSET + 2- IMG_PAD_SIZE;
          ref_pos_y = i/4-FILTER_OFFSET + 2- IMG_PAD_SIZE;
          ref_pos_y1 = max(0, min(img_height -1, ref_pos_y));
          ref_pos_x1 = max(0, min(img_width -1, ref_pos_x));
          ref_pos_y2 = max(0, min(img_height -1, ref_pos_y+1));
          ref_pos_x2 = max(0, min(img_width -1, ref_pos_x+1 ));

          is += (refY[ref_pos_y1][ref_pos_x1]+ refY[ref_pos_y2][ref_pos_x1] 
          + refY[ref_pos_y1][ref_pos_x2]+ refY[ref_pos_y2][ref_pos_x2]) 
            * FilterCoef[j_pos][FILTER_SIZE*2+2];
          is += FilterCoef[j_pos][SQR(FILTER_SIZE)];
        }
        else
        {
          is = 0.; 
          for(ii = 0; ii < FILTER_SIZE; ii++)
            for(jj = 0; jj < FILTER_SIZE; jj++)
            {
              ref_pos_y = max(0, min(img_height-1, i/4-IMG_PAD_SIZE-FILTER_OFFSET+ii ) );
              if(i<IMG_PAD_SIZE && ref_pos_y > 0)
                ref_pos_y -=1;
              ref_pos_x = max(0, min(img_width -1, j/4-IMG_PAD_SIZE-FILTER_OFFSET+jj ) );
              if(j<IMG_PAD_SIZE && ref_pos_x > 0)
                ref_pos_x -=1;
              is += (FilterCoef[j_pos][FILTER_SIZE*ii+jj]*refY[ref_pos_y][ref_pos_x]);
            }
        }
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y_aif[i+2][j+2] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is+0.5)));
#else
        out4Y_aif[i+2][j+2] = (imgpel)max(0,min(255,(int)(is+0.5)));
#endif


        //e,g,m,o
        if(FilterFlag[e_pos])
        {
          is0 = is1 = is2 = is3 = 0;
          ref_pos_x = j/4-FILTER_OFFSET- IMG_PAD_SIZE;
          ref_pos_y = i/4-FILTER_OFFSET- IMG_PAD_SIZE;
          for(ii = 0; ii < FILTER_SIZE; ii++)  // horizontal.
          {
            jj = ii;// filter NW-SE direction
            ref_pos_y1 = max(0, min(img_height -1, ref_pos_y+ii));
            ref_pos_x1 = max(0, min(img_width -1, ref_pos_x + ii));
            is0 += FilterCoef[e_pos][ii*FILTER_SIZE+jj]*refY[ref_pos_y1][ref_pos_x1]; //e
            is1 += FilterCoef[o_pos][ii*FILTER_SIZE+jj]*refY[ref_pos_y1][ref_pos_x1]; //o

            jj = FILTER_SIZE-1 - ii;// filter NE-SW direction
            ref_pos_x2 = max(0, min(img_width -1, ref_pos_x + jj ));
            is2 += FilterCoef[g_pos][ii*FILTER_SIZE+jj]*refY[ref_pos_y1][ref_pos_x2]; 
            is3 += FilterCoef[m_pos][ii*FILTER_SIZE+jj]*refY[ref_pos_y1][ref_pos_x2]; 
          }
          is0 += FilterCoef[e_pos][SQR(FILTER_SIZE)];
          is1 += FilterCoef[o_pos][SQR(FILTER_SIZE)];
          is2 += FilterCoef[g_pos][SQR(FILTER_SIZE)];
          is3 += FilterCoef[m_pos][SQR(FILTER_SIZE)];
        }
        else
        {
          is0 = 0.; is1 = 0.; is2 = 0.; is3 = 0.;
          for(ii = 0; ii < FILTER_SIZE; ii++)
            for(jj = 0; jj < FILTER_SIZE; jj++)
            {
              ref_pos_y = max(0, min(img_height-1, i/4-IMG_PAD_SIZE-FILTER_OFFSET+ii ) );
              if(i<IMG_PAD_SIZE && ref_pos_y > 0)
                ref_pos_y -=1;
              ref_pos_x = max(0, min(img_width -1, j/4-IMG_PAD_SIZE-FILTER_OFFSET+jj ) );
              if(j<IMG_PAD_SIZE && ref_pos_x > 0)
                ref_pos_x -=1;
              is0 += (FilterCoef[e_pos][FILTER_SIZE*ii+jj]*refY[ref_pos_y][ref_pos_x]);
              is1 += (FilterCoef[o_pos][FILTER_SIZE*ii+jj]*refY[ref_pos_y][ref_pos_x]);
              is2 += (FilterCoef[g_pos][FILTER_SIZE*ii+jj]*refY[ref_pos_y][ref_pos_x]);
              is3 += (FilterCoef[m_pos][FILTER_SIZE*ii+jj]*refY[ref_pos_y][ref_pos_x]);
            }
        }
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y_aif[i+1][j+1] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is0+0.5)));
        out4Y_aif[i+3][j+3] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is1+0.5)));
        out4Y_aif[i+1][j+3] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is2+0.5)));
        out4Y_aif[i+3][j+1] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is3+0.5)));
#else
        out4Y_aif[i+1][j+1] = (imgpel)max(0,min(255,(int)(is0+0.5)));
        out4Y_aif[i+3][j+3] = (imgpel)max(0,min(255,(int)(is1+0.5)));
        out4Y_aif[i+1][j+3] = (imgpel)max(0,min(255,(int)(is2+0.5)));
        out4Y_aif[i+3][j+1] = (imgpel)max(0,min(255,(int)(is3+0.5)));
#endif

        //  f,n,i,k
        if(FilterFlag[f_pos])
        {
          is0 = is1 = is2 = is3 = 0;
          //  set NW 0 positions.
          ref_pos_x = j/4-FILTER_OFFSET- IMG_PAD_SIZE;
          ref_pos_y = i/4-FILTER_OFFSET- IMG_PAD_SIZE;

          for(ii = 0; ii < FILTER_SIZE; ii++)  // horizontal.
          {
            jj = FILTER_SIZE-1 - ii;

            ref_pos_y1 = max(0, min(img_height -1, ref_pos_y+ii));
            ref_pos_y2 = max(0, min(img_height -1, ref_pos_y+jj));
            ref_pos_x1 = max(0, min(img_width -1, ref_pos_x + ii));
            ref_pos_x2 = max(0, min(img_width -1, ref_pos_x + jj ));

            is0 += FilterCoef[f_pos][ii*FILTER_SIZE+ii]*(refY[ref_pos_y1][ref_pos_x1] + refY[ref_pos_y1][ref_pos_x2]);//f
            is1 += FilterCoef[n_pos][ii*FILTER_SIZE+ii]*(refY[ref_pos_y1][ref_pos_x1] + refY[ref_pos_y1][ref_pos_x2]);//n
            is2 += FilterCoef[i_pos][ii*FILTER_SIZE+ii]*(refY[ref_pos_y1][ref_pos_x1] + refY[ref_pos_y2][ref_pos_x1]);//f
            is3 += FilterCoef[k_pos][ii*FILTER_SIZE+ii]*(refY[ref_pos_y1][ref_pos_x1] + refY[ref_pos_y2][ref_pos_x1]);//k
          }
          is0 += FilterCoef[f_pos][SQR(FILTER_SIZE)];
          is1 += FilterCoef[n_pos][SQR(FILTER_SIZE)];
          is2 += FilterCoef[i_pos][SQR(FILTER_SIZE)];
          is3 += FilterCoef[k_pos][SQR(FILTER_SIZE)];
        }
        else
        {
          is0 = 0.; is1 = 0.; is2 = 0.; is3 = 0.;
          for(ii = 0; ii < FILTER_SIZE; ii++)
            for(jj = 0; jj < FILTER_SIZE; jj++)
            {
              ref_pos_y = max(0, min(img_height-1, i/4-IMG_PAD_SIZE-FILTER_OFFSET+ii ) );
              if(i<IMG_PAD_SIZE && ref_pos_y > 0)
                ref_pos_y -=1;
              ref_pos_x = max(0, min(img_width -1, j/4-IMG_PAD_SIZE-FILTER_OFFSET+jj ) );
              if(j<IMG_PAD_SIZE && ref_pos_x > 0)
                ref_pos_x -=1;
              is0 += (FilterCoef[f_pos][FILTER_SIZE*ii+jj]*refY[ref_pos_y][ref_pos_x]);
              is1 += (FilterCoef[n_pos][FILTER_SIZE*ii+jj]*refY[ref_pos_y][ref_pos_x]);
              is2 += (FilterCoef[i_pos][FILTER_SIZE*ii+jj]*refY[ref_pos_y][ref_pos_x]);
              is3 += (FilterCoef[k_pos][FILTER_SIZE*ii+jj]*refY[ref_pos_y][ref_pos_x]);
            }
        }

#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y_aif[i+1][j+2] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is0+0.5)));
        out4Y_aif[i+3][j+2] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is1+0.5)));
        out4Y_aif[i+2][j+1] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is2+0.5)));
        out4Y_aif[i+2][j+3] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is3+0.5)));
#else
        out4Y_aif[i+1][j+2] = (imgpel)max(0,min(255,(int)(is0+0.5)));
        out4Y_aif[i+3][j+2] = (imgpel)max(0,min(255,(int)(is1+0.5)));
        out4Y_aif[i+2][j+1] = (imgpel)max(0,min(255,(int)(is2+0.5)));
        out4Y_aif[i+2][j+3] = (imgpel)max(0,min(255,(int)(is3+0.5)));
#endif
      }

      // Generate 1/1th pel representation (used for integer pel MV search)
      {
        int x, y, yy , y_pos;

        for (y = 0; y < img_height; y++)
        {
          yy = (y + IMG_PAD_SIZE)<<2;
          y_pos = y * img_width;
          for (x = 0; x < img_width; x++)
            ref11_aif[y_pos + x] = out4Y_aif[yy][(x + IMG_PAD_SIZE)<<2];
        }
      }
  }
}
#endif

#ifdef EAIF
/*! 
*************************************************************************************
* \brief
*   Up-sample reference frame(s) with ENHANCED directional adaptive interpolation filter (EDAIF)
*  
*
* \param <parameter>
*    No input parameters
*
* \return
*    No output return
*
* \note
*   1. Uses adaptive filters stored in FilterCoef[15][SQR_FILTER]
*   2. UnifiedOneForthPixWith_EAIF is a faster alternative to UnifiedOneForthPixWithNewFilter function if EAIF is used for interpolation. 
*   3. UnifiedOneForthPixWith_EAIF and UnifiedOneForthPixWithNewFilter provide identical interpolated frame(s). 
*
* \author
*    - Yan Ye                   <yye@qualcomm.com>
*************************************************************************************
*/

void UnifiedOneForthPixWith_EAIF(void)
{

  int i, j,ii, jj;
  unsigned frame;
  imgpel **out4Y_aif, **out4Y;
  imgpel **refY;
  int img_width, img_height;
  double is, is0, is1,is2,is3;
  imgpel or_pel;
  imgpel *ref11_aif;
  int pos_x, pos_y; 

  if(current_poc == img->ThisPOC)
    return;// avoids multiple interolating of the same frames
  else
    current_poc = img->ThisPOC;

  for(frame = 0; frame < dpb.ref_frames_in_buffer; frame++)
  {
    img_width  = dpb.fs_ref[frame]->frame->size_x;
    img_height = dpb.fs_ref[frame]->frame->size_y;
    out4Y_aif  = dpb.fs_ref[frame]->frame->imgY_ups_aif;
    out4Y      = dpb.fs_ref[frame]->frame->imgY_ups;
    refY       = dpb.fs_ref[frame]->frame->imgY;
    ref11_aif = dpb.fs_ref[frame]->frame->imgY_11_aif;

    for(i=0; i < (2*IMG_PAD_SIZE+img_height)*4; i+=4)
      for(j=0; j < (2*IMG_PAD_SIZE+img_width)*4; j+=4)
      {
        if(!FilterFlagInt)
        {
          pos_x = j/4- IMG_PAD_SIZE;
          pos_y = i/4- IMG_PAD_SIZE;
          pos_x = max(0, min(img_width -1, pos_x));
          pos_y = max(0, min(img_height -1, pos_y));
          out4Y_aif[i][j] = refY[pos_y][pos_x];
        }
        else
        {
          is = 0.;
          for(ii = 0; ii < FILTER_SIZE_INT; ii++)
            for(jj = 0; jj < FILTER_SIZE_INT; jj++)
            {
              pos_x = max(0, min(img_width -1, j/4-IMG_PAD_SIZE-FILTER_OFFSET_INT+jj ) );
              pos_y = max(0, min(img_height-1, i/4-IMG_PAD_SIZE-FILTER_OFFSET_INT+ii ) );
              is += (FilterCoefInt[FILTER_SIZE_INT*ii+jj] * refY[pos_y][pos_x]);
            }

            is += FilterCoefInt[SQR(FILTER_SIZE_INT)];

#ifdef INTERNAL_BIT_DEPTH_INCREASE
            out4Y_aif[i][j] = max(0,min(((1<<img->bitdepth_luma)-1),(int)(is+0.5)));
#else
            out4Y_aif[i][j] = max(0,min(255,(int)(is+0.5)));
#endif  // INTERNAL_BIT_DEPTH_INCREASE
        }

        //  a,b,c
        is = is0 = is1 = is2 = 0;
        ii = FILTER_SIZE*FILTER_OFFSET;
        pos_y = i/4 - IMG_PAD_SIZE;
        pos_y = max(0, min(img_height -1, pos_y));
        for(jj = 0; jj < FILTER_SIZE; jj++)  // horizontal.
        {
          pos_x = j/4-FILTER_OFFSET + jj - IMG_PAD_SIZE;
          pos_x = max(0, min(img_width -1, pos_x));
          or_pel = (imgpel)refY[pos_y][pos_x];
          is0 += (FilterCoef[a_pos][ii+jj]*or_pel);
          is1 += (FilterCoef[b_pos][ii+jj]*or_pel);
          is2 += (FilterCoef[c_pos][ii+jj]*or_pel);
        }
        is0 += FilterCoef[a_pos][SQR(FILTER_SIZE)];
        is1 += FilterCoef[b_pos][SQR(FILTER_SIZE)];
        is2 += FilterCoef[c_pos][SQR(FILTER_SIZE)];
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y_aif[i][j+1] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is0+0.5)));
        out4Y_aif[i][j+2] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is1+0.5)));
        out4Y_aif[i][j+3] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is2+0.5)));
#else
        out4Y_aif[i][j+1] = (imgpel)max(0,min(255,(int)(is0+0.5)));
        out4Y_aif[i][j+2] = (imgpel)max(0,min(255,(int)(is1+0.5)));
        out4Y_aif[i][j+3] = (imgpel)max(0,min(255,(int)(is2+0.5)));
#endif

        //  d,h,l
        is = is0 = is1 = is2 = 0;
        jj = FILTER_OFFSET;
        pos_x = j/4  - IMG_PAD_SIZE;
        pos_x = max(0, min(img_width -1, pos_x));

        for(ii = 0; ii < FILTER_SIZE; ii++)  // horizontal.
        {
          pos_y = i/4-FILTER_OFFSET + ii - IMG_PAD_SIZE;
          pos_y = max(0, min(img_height -1, pos_y));
          or_pel = (imgpel)refY[pos_y][pos_x];
          is0 += (FilterCoef[d_pos][ii*FILTER_SIZE+jj]*or_pel);
          is1 += (FilterCoef[h_pos][ii*FILTER_SIZE+jj]*or_pel);
          is2 += (FilterCoef[l_pos][ii*FILTER_SIZE+jj]*or_pel);
        }
        is0 += FilterCoef[d_pos][SQR(FILTER_SIZE)];
        is1 += FilterCoef[h_pos][SQR(FILTER_SIZE)];
        is2 += FilterCoef[l_pos][SQR(FILTER_SIZE)];
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y_aif[i+1][j] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is0+0.5)));
        out4Y_aif[i+2][j] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is1+0.5)));
        out4Y_aif[i+3][j] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is2+0.5)));
#else
        out4Y_aif[i+1][j] = (imgpel)max(0,min(255,(int)(is0+0.5)));
        out4Y_aif[i+2][j] = (imgpel)max(0,min(255,(int)(is1+0.5)));
        out4Y_aif[i+3][j] = (imgpel)max(0,min(255,(int)(is2+0.5)));
#endif

        //j 
        if(FilterFlag[j_pos])
        {
          is = 0.; 
          for(ii = 1; ii < FILTER_SIZE-1; ii++)
            for(jj = 1; jj < FILTER_SIZE-1; jj++)
            {
              pos_y = max(0, min(img_height-1, i/4-IMG_PAD_SIZE-FILTER_OFFSET+ii ) );
              if(i<IMG_PAD_SIZE && pos_y > 0)
                pos_y -=1;
              pos_x = max(0, min(img_width -1, j/4-IMG_PAD_SIZE-FILTER_OFFSET+jj ) );
              if(j<IMG_PAD_SIZE && pos_x > 0)
                pos_x -=1;
              is += (FilterCoef[j_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
            }
            is += FilterCoef[j_pos][SQR(FILTER_SIZE)];
        }
        else
        {
          is = 0.; 
          for(ii = 0; ii < FILTER_SIZE; ii++)
            for(jj = 0; jj < FILTER_SIZE; jj++)
            {
              pos_y = max(0, min(img_height-1, i/4-IMG_PAD_SIZE-FILTER_OFFSET+ii ) );
              if(i<IMG_PAD_SIZE && pos_y > 0)
                pos_y -=1;
              pos_x = max(0, min(img_width -1, j/4-IMG_PAD_SIZE-FILTER_OFFSET+jj ) );
              if(j<IMG_PAD_SIZE && pos_x > 0)
                pos_x -=1;
              is += (FilterCoef[j_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
            }
        }
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y_aif[i+2][j+2] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is+0.5)));
#else
        out4Y_aif[i+2][j+2] = (imgpel)max(0,min(255,(int)(is+0.5)));
#endif

        //e,g,m,o
        if(FilterFlag[e_pos])
        {
          is0 = 0.; is1 = 0.; is2 = 0.; is3 = 0.;
          for(ii = 1; ii < FILTER_SIZE-1; ii++)
            for(jj = 1; jj < FILTER_SIZE-1; jj++)
            {
              pos_y = max(0, min(img_height-1, i/4-IMG_PAD_SIZE-FILTER_OFFSET+ii ) );
              if(i<IMG_PAD_SIZE && pos_y > 0)
                pos_y -=1;
              pos_x = max(0, min(img_width -1, j/4-IMG_PAD_SIZE-FILTER_OFFSET+jj ) );
              if(j<IMG_PAD_SIZE && pos_x > 0)
                pos_x -=1;
              is0 += (FilterCoef[e_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
              is1 += (FilterCoef[o_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
              is2 += (FilterCoef[g_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
              is3 += (FilterCoef[m_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
            }
            is0 += FilterCoef[e_pos][SQR(FILTER_SIZE)];
            is1 += FilterCoef[o_pos][SQR(FILTER_SIZE)];
            is2 += FilterCoef[g_pos][SQR(FILTER_SIZE)];
            is3 += FilterCoef[m_pos][SQR(FILTER_SIZE)];
        }
        else
        {
          is0 = 0.; is1 = 0.; is2 = 0.; is3 = 0.;
          for(ii = 0; ii < FILTER_SIZE; ii++)
            for(jj = 0; jj < FILTER_SIZE; jj++)
            {
              pos_y = max(0, min(img_height-1, i/4-IMG_PAD_SIZE-FILTER_OFFSET+ii ) );
              if(i<IMG_PAD_SIZE && pos_y > 0)
                pos_y -=1;
              pos_x = max(0, min(img_width -1, j/4-IMG_PAD_SIZE-FILTER_OFFSET+jj ) );
              if(j<IMG_PAD_SIZE && pos_x > 0)
                pos_x -=1;
              is0 += (FilterCoef[e_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
              is1 += (FilterCoef[o_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
              is2 += (FilterCoef[g_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
              is3 += (FilterCoef[m_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
            }
        }
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y_aif[i+1][j+1] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is0+0.5)));
        out4Y_aif[i+3][j+3] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is1+0.5)));
        out4Y_aif[i+1][j+3] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is2+0.5)));
        out4Y_aif[i+3][j+1] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is3+0.5)));
#else
        out4Y_aif[i+1][j+1] = (imgpel)max(0,min(255,(int)(is0+0.5)));
        out4Y_aif[i+3][j+3] = (imgpel)max(0,min(255,(int)(is1+0.5)));
        out4Y_aif[i+1][j+3] = (imgpel)max(0,min(255,(int)(is2+0.5)));
        out4Y_aif[i+3][j+1] = (imgpel)max(0,min(255,(int)(is3+0.5)));
#endif

        //  f,n
        if(FilterFlag[f_pos])
        {
          is0 = 0.; is1 = 0.; 
          for(ii = 1; ii < FILTER_SIZE-1; ii++)
            for(jj = 1; jj < FILTER_SIZE-1; jj++)
            {
              pos_y = max(0, min(img_height-1, i/4-IMG_PAD_SIZE-FILTER_OFFSET+ii ) );
              if(i<IMG_PAD_SIZE && pos_y > 0)
                pos_y -=1;
              pos_x = max(0, min(img_width -1, j/4-IMG_PAD_SIZE-FILTER_OFFSET+jj ) );
              if(j<IMG_PAD_SIZE && pos_x > 0)
                pos_x -=1;
              is0 += (FilterCoef[f_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
              is1 += (FilterCoef[n_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
            }
            is0 += FilterCoef[f_pos][SQR(FILTER_SIZE)];
            is1 += FilterCoef[n_pos][SQR(FILTER_SIZE)];
        }
        else
        {
          is0 = 0.; is1 = 0.; 
          for(ii = 0; ii < FILTER_SIZE; ii++)
            for(jj = 0; jj < FILTER_SIZE; jj++)
            {
              pos_y = max(0, min(img_height-1, i/4-IMG_PAD_SIZE-FILTER_OFFSET+ii ) );
              if(i<IMG_PAD_SIZE && pos_y > 0)
                pos_y -=1;
              pos_x = max(0, min(img_width -1, j/4-IMG_PAD_SIZE-FILTER_OFFSET+jj ) );
              if(j<IMG_PAD_SIZE && pos_x > 0)
                pos_x -=1;
              is0 += (FilterCoef[f_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
              is1 += (FilterCoef[n_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
            }
        }
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y_aif[i+1][j+2] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is0+0.5)));
        out4Y_aif[i+3][j+2] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is1+0.5)));
#else
        out4Y_aif[i+1][j+2] = (imgpel)max(0,min(255,(int)(is0+0.5)));
        out4Y_aif[i+3][j+2] = (imgpel)max(0,min(255,(int)(is1+0.5)));
#endif

        //  i,k
        if(FilterFlag[i_pos])
        {
          is2 = 0.; is3 = 0.;
          for(ii = 1; ii < FILTER_SIZE-1; ii++)
            for(jj = 1; jj < FILTER_SIZE-1; jj++)
            {
              pos_y = max(0, min(img_height-1, i/4-IMG_PAD_SIZE-FILTER_OFFSET+ii ) );
              if(i<IMG_PAD_SIZE && pos_y > 0)
                pos_y -=1;
              pos_x = max(0, min(img_width -1, j/4-IMG_PAD_SIZE-FILTER_OFFSET+jj ) );
              if(j<IMG_PAD_SIZE && pos_x > 0)
                pos_x -=1;
              is2 += (FilterCoef[i_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
              is3 += (FilterCoef[k_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
            }
            is2 += FilterCoef[i_pos][SQR(FILTER_SIZE)];
            is3 += FilterCoef[k_pos][SQR(FILTER_SIZE)];
        }
        else
        {
          is2 = 0.; is3 = 0.;
          for(ii = 0; ii < FILTER_SIZE; ii++)
            for(jj = 0; jj < FILTER_SIZE; jj++)
            {
              pos_y = max(0, min(img_height-1, i/4-IMG_PAD_SIZE-FILTER_OFFSET+ii ) );
              if(i<IMG_PAD_SIZE && pos_y > 0)
                pos_y -=1;
              pos_x = max(0, min(img_width -1, j/4-IMG_PAD_SIZE-FILTER_OFFSET+jj ) );
              if(j<IMG_PAD_SIZE && pos_x > 0)
                pos_x -=1;
              is2 += (FilterCoef[i_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
              is3 += (FilterCoef[k_pos][FILTER_SIZE*ii+jj]*refY[pos_y][pos_x]);
            }
        }
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y_aif[i+2][j+1] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is2+0.5)));
        out4Y_aif[i+2][j+3] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),(int)(is3+0.5)));
#else
        out4Y_aif[i+2][j+1] = (imgpel)max(0,min(255,(int)(is2+0.5)));
        out4Y_aif[i+2][j+3] = (imgpel)max(0,min(255,(int)(is3+0.5)));
#endif
      }

      // Generate 1/1th pel representation (used for integer pel MV search)
      {
        int x, y, yy , y_pos;

        for (y = 0; y < img_height; y++)
        {
          yy = (y + IMG_PAD_SIZE)<<2;
          y_pos = y * img_width;
          for (x = 0; x < img_width; x++)
            ref11_aif[y_pos + x] = out4Y_aif[yy][(x + IMG_PAD_SIZE)<<2];
        }
      }
  }
}
#endif

// separable aif (BEGIN)
/*!
************************************************************************
* \brief
*    Horizontal upsampling with new calculated filter coefficients.
************************************************************************
*/
void UnifiedOneForthPixWithNewFilterHor (void)
{
  int i, j, jj;
  unsigned frame;
  double **out4Y_aif_hor; 
  imgpel **refY;
  int img_width, img_height;
  int sub_pos, pos_x;
  int pad_size_times_4 = 4 * IMG_PAD_SIZE; 
  double is;
  if(current_poc == img->ThisPOC)
    return;// avoids multiple interolating of the same frames
  else
    current_poc = img->ThisPOC;

  for (frame = 0; frame < dpb.ref_frames_in_buffer; frame++)
  {
    img_width     = dpb.fs_ref[frame]->frame->size_x;
    img_height    = dpb.fs_ref[frame]->frame->size_y;
    out4Y_aif_hor = dpb.fs_ref[frame]->frame->imgY_ups_aif_hor; 
    refY = dpb.fs_ref[frame]->frame->imgY;
    for (i = -IMG_PAD_SIZE; i < IMG_PAD_SIZE + img_height; i++)
    {
      for (j = -pad_size_times_4; j < (IMG_PAD_SIZE + img_width) * 4; j++) 
      {
        sub_pos = (j + pad_size_times_4) % 4; // x-sub-coordinate in a 4x4block 
        is = 0.0;
        if (sub_pos)
        {
          for (jj = 0; jj < FILTER_SIZE; jj++)
          {
            pos_x = max (0, min (img_width - 1, j / 4 - FILTER_OFFSET + jj));
            if (j < 0 && pos_x > 0)
              pos_x -= 1;

            is += (FilterCoef[sub_pos - 1][jj] * refY[max (0, min (img_height - 1, i))][pos_x]);
          }
          out4Y_aif_hor[i + IMG_PAD_SIZE][j + pad_size_times_4] = is; 
        }
        else
          out4Y_aif_hor[i + IMG_PAD_SIZE][j + pad_size_times_4] = 
          refY[max (0, min (img_height - 1, i))][max (0, min (img_width - 1, j / 4))];
      }
    }
  }
}

/*!
************************************************************************
* \brief
*    Vertical upsampling with new calculated filter coefficients.
************************************************************************
*/
void UnifiedOneForthPixWithNewFilterVer (void)
{
  int i, j,ii;
  unsigned frame;
  imgpel **out4Y_aif, **out4Y;
  imgpel **refY;
  double **out4Y_aif_hor; 
  int pad_size_times_4 = 4 * IMG_PAD_SIZE; 
  int img_width, img_height;
  int sub_pos, pos_y, x_sub, y_sub;
  double is;
  for (frame = 0; frame < dpb.ref_frames_in_buffer; frame++)
  {
    img_width = dpb.fs_ref[frame]->frame->size_x;
    img_height = dpb.fs_ref[frame]->frame->size_y;
    out4Y_aif_hor = dpb.fs_ref[frame]->frame->imgY_ups_aif_hor;     
    out4Y_aif = dpb.fs_ref[frame]->frame->imgY_ups_aif;
    out4Y = dpb.fs_ref[frame]->frame->imgY_ups;
    refY = dpb.fs_ref[frame]->frame->imgY;
    for(i=-pad_size_times_4; i < (IMG_PAD_SIZE+img_height)*4; i++) 
    {
      for(j=-pad_size_times_4; j < (IMG_PAD_SIZE+img_width)*4; j++) 
      {
        x_sub = (j+pad_size_times_4)%4;       // x-sub-coordinate in a 4x4block 
        y_sub = (i+pad_size_times_4)%4;       // y-sub-coordinate in a 4x4block 
        sub_pos = x_sub + 4*y_sub;          // pos 1..15 in a 4x4 block
        is = 0.0;
        if (sub_pos > 3)
        {
          for (ii = 0; ii < FILTER_SIZE; ii++)
          {
            pos_y = max (-IMG_PAD_SIZE, min (img_height + IMG_PAD_SIZE - 1, i / 4 - FILTER_OFFSET + ii));
            if (i < 0 && pos_y > 0)
              pos_y -= 1;

            is += (FilterCoef[sub_pos - 1][ii] * out4Y_aif_hor[pos_y + IMG_PAD_SIZE][j + pad_size_times_4]);    
          }
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          out4Y_aif[i + pad_size_times_4][j + pad_size_times_4] = max(0,min(((1<<img->bitdepth_luma)-1),(int)(is+0.5))); 
#else
          out4Y_aif[i + pad_size_times_4][j + pad_size_times_4] = max(0, min (255, (int) (is + 0.5))); 
#endif
        }
        else if (!sub_pos)
        {
          out4Y_aif[i + pad_size_times_4][j + pad_size_times_4] = 
            refY[max (0, min (img_height - 1, i / 4))][max (0, min (img_width - 1, j / 4))];
        }
        else
        {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          out4Y_aif[i + pad_size_times_4][j + pad_size_times_4] = max(0,min(((1<<img->bitdepth_luma)-1),(int)(out4Y_aif_hor[i / 4 + IMG_PAD_SIZE][j + pad_size_times_4]+0.5))); 
#else
          out4Y_aif[i + pad_size_times_4][j + pad_size_times_4] = max(0, min (255, (int) (out4Y_aif_hor[i / 4 + IMG_PAD_SIZE][j + pad_size_times_4] + 0.5))); 
#endif
        }
      }
    }
  }
}
// separable aif (END)

/*!
************************************************************************
* \brief
*    swap standard and new calculated
*    interpolated upsamled frames
************************************************************************
*/
void SwapUpsampledFrames(void)
{
  unsigned frame;

  for(frame = 0; frame < dpb.ref_frames_in_buffer; frame++)
  {
    dpb.fs_ref[frame]->frame->p_imgY_ups_aif = dpb.fs_ref[frame]->frame->imgY_ups_aif;
    dpb.fs_ref[frame]->frame->imgY_ups_aif = dpb.fs_ref[frame]->frame->imgY_ups;
    dpb.fs_ref[frame]->frame->imgY_ups = dpb.fs_ref[frame]->frame->p_imgY_ups_aif;
#ifdef E_DAIF
#ifdef EAIF
    if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF || input->UseAdaptiveFilter == FILTER_TYPE_EAIF)
#else
    if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)
#endif
    {
      dpb.fs_ref[frame]->frame->p_imgY_11_aif  = dpb.fs_ref[frame]->frame->imgY_11_aif;
      dpb.fs_ref[frame]->frame->imgY_11_aif    = dpb.fs_ref[frame]->frame->imgY_11;
      dpb.fs_ref[frame]->frame->imgY_11        = dpb.fs_ref[frame]->frame->p_imgY_11_aif;
    }
#endif  // E_DAIF 
  }
}

/*!
************************************************************************
* \brief
*    prints filter coefficient
************************************************************************
*/
void PrintFilterCoefDouble(double FCoef[15][SQR_FILTER])
{
  int i,j;
  for(i = a_pos; i <= o_pos; i++)
  {
    //    printf("Filter_%c=[\n",97+i);
    printf("FilterFlag[%c]=%d\n",97+i,FilterFlag[i]);
    printf("Filter(%d,1:6,1:6)=[\n",i+1);
    for(j = 0; j < SQR(FILTER_SIZE)+1; j++)
    {
      if(!(j%FILTER_SIZE))
        printf("[");
      printf("%12.6f   ",FCoef[i][j]);
      if(j%FILTER_SIZE == FILTER_SIZE-1)
        printf("];\n");
      else
        printf(", ");
    }
    printf("];\n");
  }
}

/*!
************************************************************************
* \brief
*    prints filter coefficients
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

/*!
************************************************************************
* \brief
*    prints filter coefficient
************************************************************************
*/
void PrintFilterCoefDoubleSep (double FCoef[15][SQR_FILTER])
{
  int i,j;
  printf("\n");
  for(i = a_pos; i <= o_pos; i++)
  {
    for(j = 0; j < FILTER_SIZE; j++)
    {
      if(!(j%FILTER_SIZE))
        printf("[");
      printf("%12.8f,",FCoef[i][j]);
    }
    printf("]");
  }
}

/*! 
*************************************************************************************
* \brief
*   Evaluate the Motion Compensated Prediction (SAD and related MCP_Cost) and select filter with minimal cost
*  
*
* \param <parameter>
*    No input parameters
*
* \return
*    Decission on the cost funciton
*
* \note
*    Supports: 
*    The standard AVC interpolation filter and AIF
*    P-Frames and B-frames
*
* \para <title>
*    <paragraph>
*
* \para
*    <another paragraph>
* \author
*    - Dmytro Rusanovskyy                   <dmytro.rusanovskyy@tut.fi>
*************************************************************************************
*/

int getDecisionOnAIF_MCP(void)
{
  int result;
  int nBits=0;
  int x, y;
  int mvx,mvy,subpel_indx;
  int mvxF=0,mvyF=0;
  int mvxB=0,mvyB=0;
  int i,j;
  int img_width = img->width, img_height=img->height;

  int ref_frameF, ref_frameB;

  imgpel **out4Y_aifF=NULL, **out4Y_aifB=NULL, **out4YF=NULL, **out4YB=NULL;
  int curr_mb_nr;
  int max_mb_nr = (img->width*img->height)/256;
  int n_interMB,interMB;
  double standardFilterError, adaptiveFilterError;

  int orig=0, st_rec=0, ad_rec=0;
  int out4Y_width  = img->width*4;
  int out4Y_height = img->height*4;
  int *st_err;
  int *ad_err;
#ifdef E_DAIF
  double lambda; 
#ifdef RDO_Q
  int qp = input->UseRDO_Q ? img->masterQP:img->qp;
#else
  int qp = img->qp;
#endif

  if ((img->type==B_SLICE) && img->nal_reference_idc)
  {
    lambda = img->lambda_me[5][qp];  
  }
  else
  {
    lambda = img->lambda_me[img->type][qp];
  }
#endif

  n_interMB = 0;
  standardFilterError=0.0;
  adaptiveFilterError=0.0;

  st_err = (int*)calloc(max_mb_nr, sizeof(int));
  ad_err = (int*)calloc(max_mb_nr, sizeof(int));
  for(i = 0; i < max_mb_nr; i++)
    ad_err[i]    = st_err[i] = 0;

  img_width  = dpb.fs_ref[0]->frame->size_x;
  img_height = dpb.fs_ref[0]->frame->size_y;
  for(curr_mb_nr = 0; curr_mb_nr < max_mb_nr; curr_mb_nr++)
  {
    x = curr_mb_nr*16%img->width;
    y = (curr_mb_nr*16/img_width)*16;
    interMB = 0;

    for(i = 0; i < 16; i++)
    {
      for(j = 0; j < 16;  j++)
      {
        ref_frameF = enc_picture->ref_idx[LIST_0][(y+j)/4][(x+i)/4];
        ref_frameB = enc_picture->ref_idx[LIST_1][(y+j)/4][(x+i)/4];
        subpel_indx = -1;
        if(ref_frameF != -1)
        {
          interMB = 1;
          out4Y_aifF   = listX[LIST_0][ref_frameF]->imgY_ups_aif;
          out4YF       = listX[LIST_0][ref_frameF]->imgY_ups;
          mvxF = enc_picture->mv[LIST_0][(y+j)/4][(x+i)/4][0];
          mvyF = enc_picture->mv[LIST_0][(y+j)/4][(x+i)/4][1];

          mvx = mvxF;
          mvy = mvyF;
          mvx = (mvx >=0) ? (mvx%4) : ((4 + mvx%4)%4);
          mvy = (mvy >=0) ? (mvy%4) : ((4 + mvy%4)%4);
          subpel_indx = mvx + 4*mvy - 1;
        }
        if(ref_frameB != -1)
        {
          interMB = 1;
          out4Y_aifB   = listX[LIST_1][ref_frameB]->imgY_ups_aif;
          out4YB        = listX[LIST_1][ref_frameB]->imgY_ups;
          mvxB         = enc_picture->mv[LIST_1][(y+j)/4][(x+i)/4][0];
          mvyB         = enc_picture->mv[LIST_1][(y+j)/4][(x+i)/4][1];
          mvx = mvxB;
          mvy = mvyB;
          mvx = (mvx >=0) ? (mvx%4) : ((4 + mvx%4)%4);
          mvy = (mvy >=0) ? (mvy%4) : ((4 + mvy%4)%4);
          subpel_indx = mvx + 4*mvy - 1;
        }

        orig   = (int)(imgY_org[y+j][x+i]);
        if(ref_frameB != -1 && ref_frameF != -1)
        {
          st_rec = (int)(out4YF    [max(0,min(out4Y_height,4*(y+j)+4*IMG_PAD_SIZE+mvyF))]
          [max(0,min(out4Y_width, 4*(x+i)+4*IMG_PAD_SIZE+mvxF))]+
            out4YB    [max(0,min(out4Y_height,4*(y+j)+4*IMG_PAD_SIZE+mvyB))]
          [max(0,min(out4Y_width, 4*(x+i)+4*IMG_PAD_SIZE+mvxB))])/2;
          ad_rec = (int)(out4Y_aifF[max(0,min(out4Y_height,4*(y+j)+4*IMG_PAD_SIZE+mvyF))]
          [max(0,min(out4Y_width, 4*(x+i)+4*IMG_PAD_SIZE+mvxF))]+
            out4Y_aifB[max(0,min(out4Y_height,4*(y+j)+4*IMG_PAD_SIZE+mvyB))]
          [max(0,min(out4Y_width, 4*(x+i)+4*IMG_PAD_SIZE+mvxB))])/2;
        }
        else if(ref_frameB != -1)
        {
          st_rec = (int)(out4YB    [max(0,min(out4Y_height,4*(y+j)+4*IMG_PAD_SIZE+mvyB))]
          [max(0,min(out4Y_width, 4*(x+i)+4*IMG_PAD_SIZE+mvxB))]);
          ad_rec = (int)(out4Y_aifB[max(0,min(out4Y_height,4*(y+j)+4*IMG_PAD_SIZE+mvyB))]
          [max(0,min(out4Y_width, 4*(x+i)+4*IMG_PAD_SIZE+mvxB))]);
        }
        else if(ref_frameF != -1)
        {
          st_rec = (int)(out4YF    [max(0,min(out4Y_height,4*(y+j)+4*IMG_PAD_SIZE+mvyF))]
          [max(0,min(out4Y_width, 4*(x+i)+4*IMG_PAD_SIZE+mvxF))]);
          ad_rec = (int)(out4Y_aifF[max(0,min(out4Y_height,4*(y+j)+4*IMG_PAD_SIZE+mvyF))]

          [max(0,min(out4Y_width, 4*(x+i)+4*IMG_PAD_SIZE+mvxF))]);

        }
        st_err[curr_mb_nr] += abs(orig-st_rec);
        ad_err[curr_mb_nr] += abs(orig-ad_rec);
      }
    }
    if (interMB)
      n_interMB++;

    standardFilterError += (double)st_err[curr_mb_nr];
    adaptiveFilterError += (double)ad_err[curr_mb_nr];

  }

  if (input->UseAdaptiveFilter == 2) // separable aif
    nBits = nBitsAIF_Sep();
  else if (input->UseAdaptiveFilter == 1)
    nBits = nBitsAIF();
#ifdef DIRECTIONAL_FILTER
  else if(input->UseAdaptiveFilter == 3)
    nBits = sendCoefsAIF(NULL);
#endif
#ifdef E_DAIF
  else if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)
    nBits = sendCoefs_EDAIF(NULL);
#endif

  cost_frame[0] = (double)standardFilterError/n_interMB;
#ifdef E_DAIF // lambda fix 
  cost_frame[1] = (double)adaptiveFilterError/n_interMB + (double)nBits*lambda/n_interMB;
#else
  cost_frame[1] = (double)adaptiveFilterError/n_interMB + nBits*(double)img->lambda_me[img->type][img->qp])/n_interMB;
#endif

  //  printf("AVC_cost: %.2f\n",cost_frame[0]);
  //  printf("AIF_cost: %.2f\n",cost_frame[1]);

  if (cost_frame[0]<=cost_frame[1])
    result = 0;
  else 
    result = 1;

  free(st_err);
  free(ad_err);

  return result;
}
#ifdef DIRECTIONAL_FILTER
/*! 
*************************************************************************************
* \brief
*   Up-sample reference frame(s) with directional (diagonal) adaptive interpolation filter 1D-AIF, 16-bit implementation
*  
*
* \param <parameter>
*    No input parameters
*
* \return
*    No output return
*
* \note
*   1. Activated with user input paramter ImpType=1
*   2. Uses adaptive filters stored in FilterCoef16bits[15][SQR_FILTER]
*   3. FilterCoef16bits filters is designed with RepresentCoefsIn16bits and ApplySumRestriction functions.
*
* \author
*    - Dmytro Rusanovskyy                   <dmytro.rusanovskyy@tut.fi>
*************************************************************************************
*/
void UnifiedOneForthPixWith_1DAIF_int16(void)
{

  int i, j,ii, jj,filt_coord;
  unsigned frame;
  imgpel **out4Y_aif, **out4Y;
  imgpel **refY;
  int img_width, img_height;
  int ref_pos_x,ref_pos_y;
  int ref_pos_x1,ref_pos_y1;
  int ref_pos_x2;

  unsigned short int is, is0, is1,is2,is3,is4;
  short int reg0,reg1,reg2,reg3,reg4,reg5;
  short int reg0_L, reg1_L, reg2_L, reg3_L, reg4_L;
  short int reg0_R, reg1_R, reg2_R, reg3_R, reg4_R;

  imgpel or_pel;
  short int FILTCOEF_BITS = 7;
  short int shift_half = (1 << (FILTCOEF_BITS-1));

  if (input->ImpType != IMP_INT16)
  {
    printf("Incorrect use of DAIF-16bit interpolation!\n");
    return;  // no restriction with floatin point implementation
  }
  current_poc = img->ThisPOC;

  for(frame = 0; frame < dpb.ref_frames_in_buffer; frame++)
  {
    img_width  = dpb.fs_ref[frame]->frame->size_x;
    img_height = dpb.fs_ref[frame]->frame->size_y;
    out4Y_aif  = dpb.fs_ref[frame]->frame->imgY_ups_aif;
    out4Y      = dpb.fs_ref[frame]->frame->imgY_ups;
    refY       = dpb.fs_ref[frame]->frame->imgY;
    for(i=0; i < (2*IMG_PAD_SIZE+img_height)*4; i+=4)
      for(j=0; j < (2*IMG_PAD_SIZE+img_width)*4; j+=4)
      {
        //  Copy Integer-pel sample
        ref_pos_x = j/4- IMG_PAD_SIZE;
        ref_pos_y = i/4- IMG_PAD_SIZE;
        ref_pos_x = max(0, min(img_width -1, ref_pos_x));
        ref_pos_y = max(0, min(img_height -1, ref_pos_y));
        out4Y_aif[i][j] = refY[ref_pos_y][ref_pos_x];

        //  a,b,c
        {
          is = is0 = is1 = is2 = 0;
          reg0_L = reg1_L = reg2_L = 0;
          reg0_R = reg1_R = reg2_R = 0;

          ii = FILTER_SIZE*FILTER_OFFSET;
          ref_pos_y = i/4 - IMG_PAD_SIZE;
          ref_pos_y = max(0, min(img_height -1, ref_pos_y));
          reg0 = reg1 = reg2 = 0;
          reg3 = reg4 = reg5 = 0;
          for(jj = 0; jj < FILTER_SIZE/2; jj++)  // left_half
          {
            filt_coord = ii+jj;
            ref_pos_x = j/4-FILTER_OFFSET + jj - IMG_PAD_SIZE;
            ref_pos_x = max(0, min(img_width -1, ref_pos_x));
            or_pel = (imgpel)refY[ref_pos_y][ref_pos_x];

            reg0_L += (FilterCoef16bits[a_pos][filt_coord]*or_pel);//a
            reg1_L += (FilterCoef16bits[b_pos][filt_coord]*or_pel);//b
            reg2_L += (FilterCoef16bits[c_pos][filt_coord]*or_pel);//c
          }
          reg0_L = (short int)max(0,(reg0_L));
          reg1_L = (short int)max(0,(reg1_L));
          reg2_L = (short int)max(0,(reg2_L));

          for(jj = FILTER_SIZE/2;jj < FILTER_SIZE; jj++)  // right_half
          {
            filt_coord = ii+jj;
            ref_pos_x = j/4-FILTER_OFFSET + jj - IMG_PAD_SIZE;
            ref_pos_x = max(0, min(img_width -1, ref_pos_x));
            or_pel = (imgpel)refY[ref_pos_y][ref_pos_x];

            reg0_R += (FilterCoef16bits[a_pos][filt_coord]*or_pel);//a
            reg1_R += (FilterCoef16bits[b_pos][filt_coord]*or_pel);//b
            reg2_R += (FilterCoef16bits[c_pos][filt_coord]*or_pel);//c
          }
          reg0_R = (short int)max(0,(reg0_R));
          reg1_R = (short int)max(0,(reg1_R));
          reg2_R = (short int)max(0,(reg2_R));

          is0 = (unsigned short)reg0_L+(unsigned short)reg0_R;
          is1 = (unsigned short)reg1_L+(unsigned short)reg1_R;
          is2 = (unsigned short)reg2_L+(unsigned short)reg2_R;


          is0 = (is0+shift_half) >> FILTCOEF_BITS;
          is1 = (is1+shift_half) >> FILTCOEF_BITS;
          is2 = (is2+shift_half) >> FILTCOEF_BITS;

#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          out4Y_aif[i][j+1] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),is0));//a
          out4Y_aif[i][j+2] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),is1));//b
          out4Y_aif[i][j+3] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),is2));//c
#else
          out4Y_aif[i][j+1] = (imgpel)max(0,min(255,is0));//a
          out4Y_aif[i][j+2] = (imgpel)max(0,min(255,is1));//b
          out4Y_aif[i][j+3] = (imgpel)max(0,min(255,is2));//c
#endif
        }
        //  d,h,l
        {
          is = is0 = is1 = is2 = 0;
          reg0_L = reg1_L = reg2_L = 0;
          reg0_R = reg1_R = reg2_R = 0;
          jj = FILTER_OFFSET;
          ref_pos_x = j/4  - IMG_PAD_SIZE;
          ref_pos_x = max(0, min(img_width -1, ref_pos_x));

          for(ii = 0; ii < FILTER_SIZE/2; ii++)  // left_half
          {
            filt_coord = ii*FILTER_SIZE+jj;
            ref_pos_y = i/4-FILTER_OFFSET + ii - IMG_PAD_SIZE;
            ref_pos_y = max(0, min(img_height -1, ref_pos_y));
            or_pel = (imgpel)refY[ref_pos_y][ref_pos_x];

            reg0_L += (FilterCoef16bits[d_pos][filt_coord]*or_pel);//d
            reg1_L += (FilterCoef16bits[h_pos][filt_coord]*or_pel);//h
            reg2_L += (FilterCoef16bits[l_pos][filt_coord]*or_pel);//l
          }
          reg0_L = (short int)max(0,(reg0_L));
          reg1_L = (short int)max(0,(reg1_L));
          reg2_L = (short int)max(0,(reg2_L));

          for(ii = FILTER_SIZE/2; ii < FILTER_SIZE; ii++)  // horizontal.
          {
            filt_coord = ii*FILTER_SIZE+jj;
            ref_pos_y = i/4-FILTER_OFFSET + ii - IMG_PAD_SIZE;
            ref_pos_y = max(0, min(img_height -1, ref_pos_y));
            or_pel = (imgpel)refY[ref_pos_y][ref_pos_x];

            reg0_R += (FilterCoef16bits[d_pos][filt_coord]*or_pel);//d
            reg1_R += (FilterCoef16bits[h_pos][filt_coord]*or_pel);//h
            reg2_R += (FilterCoef16bits[l_pos][filt_coord]*or_pel);//l
          }
          reg0_R = (short int)max(0,(reg0_R));
          reg1_R = (short int)max(0,(reg1_R));
          reg2_R = (short int)max(0,(reg2_R));

          is0 = (unsigned short)reg0_L+(unsigned short)reg0_R;
          is1 = (unsigned short)reg1_L+(unsigned short)reg1_R;
          is2 = (unsigned short)reg2_L+(unsigned short)reg2_R;


          is0 = (is0+shift_half) >> FILTCOEF_BITS;
          is1 = (is1+shift_half) >> FILTCOEF_BITS;
          is2 = (is2+shift_half) >> FILTCOEF_BITS;

#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          out4Y_aif[i+1][j] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),is0));//d
          out4Y_aif[i+2][j] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),is1));//h
          out4Y_aif[i+3][j] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),is2));//l
#else
          out4Y_aif[i+1][j] = (imgpel)max(0,min(255,is0));//d
          out4Y_aif[i+2][j] = (imgpel)max(0,min(255,is1));//h
          out4Y_aif[i+3][j] = (imgpel)max(0,min(255,is2));//l
#endif
        }

        //e,g,m,o
        {
          is0 = is1 = is2 = is3 = 0;
          reg0_L = reg1_L = reg2_L = reg3_L = 0;
          reg0_R = reg1_R = reg2_R = reg3_R = 0;

          ref_pos_x = j/4-FILTER_OFFSET- IMG_PAD_SIZE;
          ref_pos_y = i/4-FILTER_OFFSET- IMG_PAD_SIZE;
          for(ii = 0; ii < FILTER_SIZE/2; ii++)  // horizontal.
          {
            jj = ii;// filter NW-SE direction
            filt_coord = ii*FILTER_SIZE+jj;
            ref_pos_y1 = max(0, min(img_height -1, ref_pos_y+ii));
            ref_pos_x1 = max(0, min(img_width -1, ref_pos_x + ii));
            or_pel = (imgpel)refY[ref_pos_y1][ref_pos_x1];

            reg0_L += FilterCoef16bits[e_pos][filt_coord]*or_pel; //e
            reg1_L += FilterCoef16bits[o_pos][filt_coord]*or_pel; //o

            jj = FILTER_SIZE-1 - ii;// filter NE-SW direction
            filt_coord = ii*FILTER_SIZE+jj;
            ref_pos_x2 = max(0, min(img_width -1, ref_pos_x + jj ));
            or_pel = (imgpel)refY[ref_pos_y1][ref_pos_x2];

            reg2_L += FilterCoef16bits[g_pos][filt_coord]*or_pel; 
            reg3_L += FilterCoef16bits[m_pos][filt_coord]*or_pel; 
          }
          reg0_L = (short int)max(0,reg0_L);
          reg1_L = (short int)max(0,reg1_L);
          reg2_L = (short int)max(0,reg2_L);
          reg3_L = (short int)max(0,reg3_L);

          for(ii = FILTER_SIZE/2; ii < FILTER_SIZE; ii++)  // horizontal.
          {
            jj = ii;// filter NW-SE direction
            filt_coord = ii*FILTER_SIZE+jj;
            ref_pos_y1 = max(0, min(img_height -1, ref_pos_y+ii));
            ref_pos_x1 = max(0, min(img_width -1, ref_pos_x + ii));
            or_pel = (imgpel)refY[ref_pos_y1][ref_pos_x1];

            reg0_R += FilterCoef16bits[e_pos][filt_coord]*or_pel; //e
            reg1_R += FilterCoef16bits[o_pos][filt_coord]*or_pel; //o

            jj = FILTER_SIZE-1 - ii;// filter NE-SW direction
            filt_coord = ii*FILTER_SIZE+jj;
            ref_pos_x2 = max(0, min(img_width -1, ref_pos_x + jj ));
            or_pel = (imgpel)refY[ref_pos_y1][ref_pos_x2];

            reg2_R += FilterCoef16bits[g_pos][filt_coord]*or_pel; 
            reg3_R += FilterCoef16bits[m_pos][filt_coord]*or_pel; 
          }
          reg0_R = (short int)max(0,reg0_R);
          reg1_R = (short int)max(0,reg1_R);
          reg2_R = (short int)max(0,reg2_R);
          reg3_R = (short int)max(0,reg3_R);

          is0 = (unsigned short)reg0_L+(unsigned short)reg0_R;
          is1 = (unsigned short)reg1_L+(unsigned short)reg1_R;
          is2 = (unsigned short)reg2_L+(unsigned short)reg2_R;
          is3 = (unsigned short)reg3_L+(unsigned short)reg3_R;

          is0 = (is0+shift_half) >> FILTCOEF_BITS; is1 = (is1+shift_half) >> FILTCOEF_BITS;
          is2 = (is2+shift_half) >> FILTCOEF_BITS; is3 = (is3+shift_half) >> FILTCOEF_BITS;

#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          out4Y_aif[i+1][j+1] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),is0));//e
          out4Y_aif[i+3][j+3] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),is1));//o
          out4Y_aif[i+1][j+3] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),is2));//g
          out4Y_aif[i+3][j+1] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),is3));//m
#else
          out4Y_aif[i+1][j+1] = (imgpel)max(0,min(255,is0));//e
          out4Y_aif[i+3][j+3] = (imgpel)max(0,min(255,is1));//o
          out4Y_aif[i+1][j+3] = (imgpel)max(0,min(255,is2));//g
          out4Y_aif[i+3][j+1] = (imgpel)max(0,min(255,is3));//m
#endif
        }

        //j,f,k,n
        {
          is0 = is1 = is2 = is3 = is4 =0;
          reg0_L = reg1_L = reg2_L = reg3_L = reg4_L = 0;
          reg0_R = reg1_R = reg2_R = reg3_R = reg4_R = 0;

          ref_pos_x = j/4-FILTER_OFFSET- IMG_PAD_SIZE;
          ref_pos_y = i/4-FILTER_OFFSET- IMG_PAD_SIZE;
          for(ii = 0; ii < FILTER_SIZE/2; ii++)  // horizontal.
          {
            jj = ii;// filter NW-SE direction
            filt_coord = ii*FILTER_SIZE+jj;
            ref_pos_y1 = max(0, min(img_height -1, ref_pos_y+ii));
            ref_pos_x1 = max(0, min(img_width -1, ref_pos_x + ii));
            or_pel = (imgpel)refY[ref_pos_y1][ref_pos_x1];

            reg0_L += FilterCoef16bits[f_pos][filt_coord]*or_pel; //f
            reg1_L += FilterCoef16bits[i_pos][filt_coord]*or_pel; //i
            reg2_L += FilterCoef16bits[j_pos][filt_coord]*or_pel; //j
            reg3_L += FilterCoef16bits[k_pos][filt_coord]*or_pel; //k
            reg4_L += FilterCoef16bits[n_pos][filt_coord]*or_pel; //n

            jj = FILTER_SIZE-1 - ii;// filter NE-SW direction
            filt_coord = ii*FILTER_SIZE+jj;
            ref_pos_x2 = max(0, min(img_width -1, ref_pos_x + jj ));
            or_pel = (imgpel)refY[ref_pos_y1][ref_pos_x2];

            reg0_R += FilterCoef16bits[f_pos][filt_coord]*or_pel; //f
            reg1_R += FilterCoef16bits[i_pos][filt_coord]*or_pel; //i
            reg2_R += FilterCoef16bits[j_pos][filt_coord]*or_pel; //g
            reg3_R += FilterCoef16bits[k_pos][filt_coord]*or_pel; //k
            reg4_R += FilterCoef16bits[n_pos][filt_coord]*or_pel; //n
          }

          reg0_L = (short int)max(0,reg0_L);       reg1_L = (short int)max(0,reg1_L);
          reg2_L = (short int)max(0,reg2_L);      reg3_L = (short int)max(0,reg3_L);
          reg4_L = (short int)max(0,reg4_L);

          reg0_R = (short int)max(0,reg0_R);      reg1_R = (short int)max(0,reg1_R);
          reg2_R = (short int)max(0,reg2_R);      reg3_R = (short int)max(0,reg3_R);
          reg4_R = (short int)max(0,reg4_R);

          //  first 6-tap are computed and stored into unsifned short registers      
          is0 = (unsigned short)reg0_L+(unsigned short)reg0_R;//f
          is1 = (unsigned short)reg1_L+(unsigned short)reg1_R;//i
          is2 = (unsigned short)reg2_L+(unsigned short)reg2_R;//g
          is3 = (unsigned short)reg3_L+(unsigned short)reg3_R;//k
          is4 = (unsigned short)reg4_L+(unsigned short)reg4_R;//n
          is0 = (is0)>>1; is1 = (is1)>>1;  is2 = (is2)>>1; is3 = (is3)>>1; is4 = (is4)>>1;

          reg0_L = reg1_L = reg2_L = reg3_L = reg4_L = 0;
          reg0_R = reg1_R = reg2_R = reg3_R = reg4_R = 0;
          for(ii = FILTER_SIZE/2; ii < FILTER_SIZE; ii++)  
          {
            jj = ii;// filter NW-SE direction
            filt_coord = ii*FILTER_SIZE+jj;
            ref_pos_y1 = max(0, min(img_height -1, ref_pos_y+ii));
            ref_pos_x1 = max(0, min(img_width -1, ref_pos_x + ii));
            or_pel = (imgpel)refY[ref_pos_y1][ref_pos_x1];

            reg0_L += FilterCoef16bits[f_pos][filt_coord]*or_pel; //f
            reg1_L += FilterCoef16bits[i_pos][filt_coord]*or_pel; //i
            reg2_L += FilterCoef16bits[j_pos][filt_coord]*or_pel; //j
            reg3_L += FilterCoef16bits[k_pos][filt_coord]*or_pel; //k
            reg4_L += FilterCoef16bits[n_pos][filt_coord]*or_pel; //n

            jj = FILTER_SIZE-1 - ii;// filter NE-SW direction
            filt_coord = ii*FILTER_SIZE+jj;
            ref_pos_x2 = max(0, min(img_width -1, ref_pos_x + jj ));
            or_pel = (imgpel)refY[ref_pos_y1][ref_pos_x2];

            reg0_R += FilterCoef16bits[f_pos][filt_coord]*or_pel; //f
            reg1_R += FilterCoef16bits[i_pos][filt_coord]*or_pel; //i
            reg2_R += FilterCoef16bits[j_pos][filt_coord]*or_pel; //j
            reg3_R += FilterCoef16bits[k_pos][filt_coord]*or_pel; //k
            reg4_R += FilterCoef16bits[n_pos][filt_coord]*or_pel; //n
          }
          reg0_L = (short int)max(0,reg0_L);       reg1_L = (short int)max(0,reg1_L);
          reg2_L = (short int)max(0,reg2_L);      reg3_L = (short int)max(0,reg3_L);
          reg4_L = (short int)max(0,reg4_L);

          reg0_R = (short int)max(0,reg0_R);      reg1_R = (short int)max(0,reg1_R);
          reg2_R = (short int)max(0,reg2_R);      reg3_R = (short int)max(0,reg3_R);
          reg4_R = (short int)max(0,reg4_R);

          is0 += (((unsigned short)reg0_L + (unsigned short)reg0_R))>>1;//f
          is1 += (((unsigned short)reg1_L + (unsigned short)reg1_R))>>1;//i
          is2 += (((unsigned short)reg2_L + (unsigned short)reg2_R))>>1;//j
          is3 += (((unsigned short)reg3_L + (unsigned short)reg3_R))>>1;//k
          is4 += (((unsigned short)reg4_L + (unsigned short)reg4_R))>>1;//n

          is0 = (is0+shift_half) >> (FILTCOEF_BITS); is1 = (is1+shift_half) >> (FILTCOEF_BITS);
          is2 = (is2+shift_half) >> (FILTCOEF_BITS); is3 = (is3+shift_half) >> (FILTCOEF_BITS);
          is4 = (is4+shift_half) >> (FILTCOEF_BITS);

#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          out4Y_aif[i+1][j+2] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),is0));//f
          out4Y_aif[i+2][j+1] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),is1));//i
          out4Y_aif[i+2][j+2] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),is2));//j
          out4Y_aif[i+2][j+3] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),is3));//k
          out4Y_aif[i+3][j+2] = (imgpel)max(0,min(((1<<img->bitdepth_luma)-1),is4));//n
#else
          out4Y_aif[i+1][j+2] = (imgpel)max(0,min(255,is0));//f
          out4Y_aif[i+2][j+1] = (imgpel)max(0,min(255,is1));//i
          out4Y_aif[i+2][j+2] = (imgpel)max(0,min(255,is2));//j
          out4Y_aif[i+2][j+3] = (imgpel)max(0,min(255,is3));//k
          out4Y_aif[i+3][j+2] = (imgpel)max(0,min(255,is4));//n
#endif

        }
      }
  }
}
/*! 
*************************************************************************************
* \brief
*   Represent filter coefficients of DAIF with 16-bit integer values
*  
*
* \param <parameter>
*    No input parameters
*
* \return
*    No output return
*
* \note
*   1. Activated with user input paramter ImpType=1
*   2. Uses adaptive filters stored in FilterCoef16bits[15][SQR_FILTER]
*  
* \author
*    - Dmytro Rusanovskyy                   <dmytro.rusanovskyy@tut.fi>
*************************************************************************************
*/
void RepresentCoefsIn16bits(int ImpType)
{
  double temp;
  int signFC;
  int i,j,multplier;
  short int PredFilterCoef16bits[15][SQR_FILTER];  
  if (ImpType != IMP_INT16)
    return;  // no restriction with floatin point implementation

  for(j = 0; j < 15; j++)
  {
    multplier = 1<<nQBits[j];
    for(i = 0; i < SQR_FILTER; i++)
    {
      if (FilterCoef[j][i]>=0) signFC = 1; else signFC = -1;

      temp = FilterCoef[j][i]*multplier;
      temp = fabs(temp) + 0.5;
      //  define a integer representation of the FC.
      FilterCoef16bits[j][i] = (short int)temp*signFC;
      FilterCoef[j][i]= (double)FilterCoef16bits[j][i]/(double)multplier;

      if (STANDARD_2D_FILTER[j][i]>=0) signFC = 1; else signFC = -1;
      temp = STANDARD_2D_FILTER[j][i]*multplier;
      temp = fabs(temp) + 0.5;
      //  define a integer representation of the FC.
      PredFilterCoef16bits[j][i] = (short int)temp*signFC;
      DiffQFilterCoef[j][i] = FilterCoef16bits[j][i] - PredFilterCoef16bits[j][i];
    }
  }
}
/*! 
*************************************************************************************
* \brief
*   Apply the "SUM restriction" on the DAIF filter coefficients. 
*  
*
* \param <parameter>
*    No input parameters
*
* \return
*    No output return
*
* \note
*   1. Activated with user input paramter ImpType=1
*   2. Uses adaptive filters stored in CalculatedFilter[15][SQR_FILTER]
*   3. If sum of equi-sign filter coefficients are >=128, and 16-bits implemention is active - substitute with static directional filter
*  
* \author
*    - Dmytro Rusanovskyy                   <dmytro.rusanovskyy@tut.fi>
*************************************************************************************
*/
void ApplySumRestriction(int ImpType, int sub_pos)
{
  int temp_pos,temp_neg;
  int i, j, i2,i3,num_coefs, multplier;

  if (ImpType != IMP_INT16)
    return;

  j = sub_pos;
  multplier = 1<<FILTCOEF_BITS;
  num_coefs = min(3,POS_EQUATION_NUMBER[sub_pos]);

  for (i2=0;i2< (int)POS_EQUATION_NUMBER[sub_pos]/3;i2++)
  {
    temp_pos = temp_neg = 0;
    for(i3 = 0; i3 < num_coefs; i3++)
    {
      i = 3*i2 + i3;
      if (CalculatedFilter[j][i]>=0) 
        temp_pos += (int)(fabs(CalculatedFilter[j][i]*multplier)+0.5);
      else 
        temp_neg += (int)(fabs(CalculatedFilter[j][i]*multplier)+0.5);

    }
    if ((temp_pos >= 128)||(temp_neg >= 128))
    {
#ifndef EDAIF2
      FilterFlag[sub_pos] = 0;
#else
      // If coefs of adaptive filter out of dynamical range, assign static filter of the same fitler structure (DAIF or RAIF)
      FilterFlag[sub_pos]= (FilterFlag[sub_pos]%2!=0)?(FilterFlag[sub_pos]-1):FilterFlag[sub_pos]; 
#endif
#ifdef _DEBUG
      printf("1D-AIF-16bits: %d sub_pel sum of f_coefs >128! AIF filter is substituted with static Directional Interpolation Filter!\n",sub_pos);
#endif
      return;
    }
  }
}

#endif
void PrintFilterCoefFloatFile(double FCoef[15][SQR_FILTER])
{
  int i,j;
  FILE* fid = fopen("enc_filter.txt","at");
  for(i = a_pos; i <= o_pos; i++)
  {
    fprintf(fid,"%c_pos: %d %d\n",97+i,FilterFlag[i], SymmetryPosition[i]);
    //fprintf(fid, "FilterFlag[%c]=%d\n",97+i,FilterFlag[i]);
    //fprintf(fid, "Filter(%d,1:6,1:6)=[\n",i+1);
    for(j = 0; j < SQR_FILTER; j++)
    {
      //if(!(j%FILTER_SIZE))
      //fprintf(fid, "[");
      fprintf(fid, "  %1.6f \t",FCoef[i][j]);
      if(j%FILTER_SIZE == FILTER_SIZE-1)
        fprintf(fid, "\n");
      //else
      //        fprintf(fid, ", ");
    }
    fprintf(fid, "\n");
  }
  fclose(fid);
}

void save2file_images_QPEL(int refFrID){
  FILE* fid;
  imgpel ** yBuf_aif;
  int i,j;
  char fileName[13];

  int extBufWidth  = (img->width + 2*IMG_PAD_SIZE)*4; 
  int extBufHeight  = (img->height + 2*IMG_PAD_SIZE)*4;


  yBuf_aif = dpb.fs_ref[0]->frame->imgY_ups_aif;
  sprintf(fileName,"%d_QPEL.txt",refFrID);

  fid = fopen(fileName,"wb");
  for(j = 0; j  <  extBufHeight ;  j ++) 
  {
    for(i = 0; i  <  extBufWidth;  i ++){ 
      fwrite(&yBuf_aif[j][i],sizeof(imgpel),1,fid);
    }
  }
  fclose(fid);

}
#endif

