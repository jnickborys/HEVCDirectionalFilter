/*!
 ***************************************************************************
 * \file transform8x8.c
 *
 * \brief
 *    8x8 transform functions
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details) 
 *    - Yuri Vatis                      <vatis@hhi.de>
 *    - Jan Muenster                    <muenster@hhi.de>
 *    - Lowell Winger                   <lwinger@lsil.com>
 * \date
 *    12. October 2003
 **************************************************************************
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "global.h"

#include "image.h"
#include "mb_access.h"
#include "elements.h"
#include "cabac.h"
#include "vlc.h"

#include "transform8x8.h"

int   cofAC8x8_chroma[2][4][2][18];

#ifdef USE_INTRA_MDDT
long         quant_stat_rd8x8[64];
#endif

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))


const int quant_coef8[6][8][8] = 
{
  { 
    {13107, 12222,  16777,  12222,  13107,  12222,  16777,  12222},
    {12222, 11428,  15481,  11428,  12222,  11428,  15481,  11428},
    {16777, 15481,  20972,  15481,  16777,  15481,  20972,  15481},
    {12222, 11428,  15481,  11428,  12222,  11428,  15481,  11428},
    {13107, 12222,  16777,  12222,  13107,  12222,  16777,  12222},
    {12222, 11428,  15481,  11428,  12222,  11428,  15481,  11428},
    {16777, 15481,  20972,  15481,  16777,  15481,  20972,  15481},
    {12222, 11428,  15481,  11428,  12222,  11428,  15481,  11428}
  },
  {
    {11916, 11058,  14980,  11058,  11916,  11058,  14980,  11058},
    {11058, 10826,  14290,  10826,  11058,  10826,  14290,  10826},
    {14980, 14290,  19174,  14290,  14980,  14290,  19174,  14290},
    {11058, 10826,  14290,  10826,  11058,  10826,  14290,  10826},
    {11916, 11058,  14980,  11058,  11916,  11058,  14980,  11058},
    {11058, 10826,  14290,  10826,  11058,  10826,  14290,  10826},
    {14980, 14290,  19174,  14290,  14980,  14290,  19174,  14290},
    {11058, 10826,  14290,  10826,  11058,  10826,  14290,  10826}
  },
  {
    {10082, 9675,   12710,  9675,   10082,  9675, 12710,  9675},
    {9675,  8943,   11985,  8943,   9675,   8943, 11985,  8943},
    {12710, 11985,  15978,  11985,  12710,  11985,  15978,  11985},
    {9675,  8943,   11985,  8943,   9675,   8943, 11985,  8943},
    {10082, 9675,   12710,  9675,   10082,  9675, 12710,  9675},
    {9675,  8943,   11985,  8943,   9675, 8943, 11985,  8943},
    {12710, 11985,  15978,  11985,  12710,  11985,  15978,  11985},
    {9675,  8943,   11985,  8943,   9675, 8943, 11985,  8943}
  },
  {
    {9362,  8931, 11984,  8931, 9362, 8931, 11984,  8931},
    {8931,  8228, 11259,  8228, 8931, 8228, 11259,  8228},
    {11984, 11259,  14913,  11259,  11984,  11259,  14913,  11259},
    {8931,  8228, 11259,  8228, 8931, 8228, 11259,  8228},
    {9362,  8931, 11984,  8931, 9362, 8931, 11984,  8931},
    {8931,  8228, 11259,  8228, 8931, 8228, 11259,  8228},
    {11984, 11259,  14913,  11259,  11984,  11259,  14913,  11259},
    {8931,  8228, 11259,  8228, 8931, 8228, 11259,  8228}
  },
  {
    {8192,  7740, 10486,  7740, 8192, 7740, 10486,  7740},
    {7740,  7346, 9777, 7346, 7740, 7346, 9777, 7346},
    {10486, 9777, 13159,  9777, 10486,  9777, 13159,  9777},
    {7740,  7346, 9777, 7346, 7740, 7346, 9777, 7346},
    {8192,  7740, 10486,  7740, 8192, 7740, 10486,  7740},
    {7740,  7346, 9777, 7346, 7740, 7346, 9777, 7346},
    {10486, 9777, 13159,  9777, 10486,  9777, 13159,  9777},
    {7740,  7346, 9777, 7346, 7740, 7346, 9777, 7346}
  },
  {
    {7282,  6830, 9118, 6830, 7282, 6830, 9118, 6830},
    {6830,  6428, 8640, 6428, 6830, 6428, 8640, 6428},
    {9118,  8640, 11570,  8640, 9118, 8640, 11570,  8640},
    {6830,  6428, 8640, 6428, 6830, 6428, 8640, 6428},
    {7282,  6830, 9118, 6830, 7282, 6830, 9118, 6830},
    {6830,  6428, 8640, 6428, 6830, 6428, 8640, 6428},
    {9118,  8640, 11570,  8640, 9118, 8640, 11570,  8640},
    {6830,  6428, 8640, 6428, 6830, 6428, 8640, 6428}
  }
};


const int dequant_coef8[6][8][8] = 
{
  {
    {20,  19, 25, 19, 20, 19, 25, 19},
    {19,  18, 24, 18, 19, 18, 24, 18},
    {25,  24, 32, 24, 25, 24, 32, 24},
    {19,  18, 24, 18, 19, 18, 24, 18},
    {20,  19, 25, 19, 20, 19, 25, 19},
    {19,  18, 24, 18, 19, 18, 24, 18},
    {25,  24, 32, 24, 25, 24, 32, 24},
    {19,  18, 24, 18, 19, 18, 24, 18}
  },
  {
    {22,  21, 28, 21, 22, 21, 28, 21},
    {21,  19, 26, 19, 21, 19, 26, 19},
    {28,  26, 35, 26, 28, 26, 35, 26},
    {21,  19, 26, 19, 21, 19, 26, 19},
    {22,  21, 28, 21, 22, 21, 28, 21},
    {21,  19, 26, 19, 21, 19, 26, 19},
    {28,  26, 35, 26, 28, 26, 35, 26},
    {21,  19, 26, 19, 21, 19, 26, 19}
  },
  {
    {26,  24, 33, 24, 26, 24, 33, 24},
    {24,  23, 31, 23, 24, 23, 31, 23},
    {33,  31, 42, 31, 33, 31, 42, 31},
    {24,  23, 31, 23, 24, 23, 31, 23},
    {26,  24, 33, 24, 26, 24, 33, 24},
    {24,  23, 31, 23, 24, 23, 31, 23},
    {33,  31, 42, 31, 33, 31, 42, 31},
    {24,  23, 31, 23, 24, 23, 31, 23}
  },
  {
    {28,  26, 35, 26, 28, 26, 35, 26},
    {26,  25, 33, 25, 26, 25, 33, 25},
    {35,  33, 45, 33, 35, 33, 45, 33},
    {26,  25, 33, 25, 26, 25, 33, 25},
    {28,  26, 35, 26, 28, 26, 35, 26},
    {26,  25, 33, 25, 26, 25, 33, 25},
    {35,  33, 45, 33, 35, 33, 45, 33},
    {26,  25, 33, 25, 26, 25, 33, 25}
  },
  {
    {32,  30, 40, 30, 32, 30, 40, 30},
    {30,  28, 38, 28, 30, 28, 38, 28},
    {40,  38, 51, 38, 40, 38, 51, 38},
    {30,  28, 38, 28, 30, 28, 38, 28},
    {32,  30, 40, 30, 32, 30, 40, 30},
    {30,  28, 38, 28, 30, 28, 38, 28},
    {40,  38, 51, 38, 40, 38, 51, 38},
    {30,  28, 38, 28, 30, 28, 38, 28}
  },
  {
    {36,  34, 46, 34, 36, 34, 46, 34},
    {34,  32, 43, 32, 34, 32, 43, 32},
    {46,  43, 58, 43, 46, 43, 58, 43},
    {34,  32, 43, 32, 34, 32, 43, 32},
    {36,  34, 46, 34, 36, 34, 46, 34},
    {34,  32, 43, 32, 34, 32, 43, 32},
    {46,  43, 58, 43, 46, 43, 58, 43},
    {34,  32, 43, 32, 34, 32, 43, 32}
  }

};


//! single scan pattern
const byte SNGL_SCAN8x8[64][2] = {
  {0,0}, {1,0}, {0,1}, {0,2}, {1,1}, {2,0}, {3,0}, {2,1}, 
  {1,2}, {0,3}, {0,4}, {1,3}, {2,2}, {3,1}, {4,0}, {5,0},
  {4,1}, {3,2}, {2,3}, {1,4}, {0,5}, {0,6}, {1,5}, {2,4},
  {3,3}, {4,2}, {5,1}, {6,0}, {7,0}, {6,1}, {5,2}, {4,3},
  {3,4}, {2,5}, {1,6}, {0,7}, {1,7}, {2,6}, {3,5}, {4,4},
  {5,3}, {6,2}, {7,1}, {7,2}, {6,3}, {5,4}, {4,5}, {3,6},
  {2,7}, {3,7}, {4,6}, {5,5}, {6,4}, {7,3}, {7,4}, {6,5},
  {5,6}, {4,7}, {5,7}, {6,6}, {7,5}, {7,6}, {6,7}, {7,7}
};


//! field scan pattern
const byte FIELD_SCAN8x8[64][2] = {   // 8x8
  {0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {0,3}, {0,4}, {1,2}, 
  {2,0}, {1,3}, {0,5}, {0,6}, {0,7}, {1,4}, {2,1}, {3,0}, 
  {2,2}, {1,5}, {1,6}, {1,7}, {2,3}, {3,1}, {4,0}, {3,2}, 
  {2,4}, {2,5}, {2,6}, {2,7}, {3,3}, {4,1}, {5,0}, {4,2}, 
  {3,4}, {3,5}, {3,6}, {3,7}, {4,3}, {5,1}, {6,0}, {5,2}, 
  {4,4}, {4,5}, {4,6}, {4,7}, {5,3}, {6,1}, {6,2}, {5,4}, 
  {5,5}, {5,6}, {5,7}, {6,3}, {7,0}, {7,1}, {6,4}, {6,5}, 
  {6,6}, {6,7}, {7,2}, {7,3}, {7,4}, {7,5}, {7,6}, {7,7}
};


//! array used to find expensive coefficients
#ifdef ADAPTIVE_QUANTIZATION
const byte COEFF_COST8x8[3][64] =
{
  {3,3,3,3,2,2,2,2,2,2,2,2,1,1,1,1,
  1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
  {9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
   9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
   9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
   9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9},
  {1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
};
#else
const byte COEFF_COST8x8[2][64] =
{
  {3,3,3,3,2,2,2,2,2,2,2,2,1,1,1,1,
  1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
  {9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
   9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
   9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
   9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9}
};
#endif

#ifdef RDO_Q
const int estErr8x8[6][8][8]={
  {
    {6553600, 6677056, 6400000, 6677056, 6553600, 6677056, 6400000, 6677056}, 
    {6677056, 6765201, 6658560, 6765201, 6677056, 6765201, 6658560, 6765201}, 
    {6400000, 6658560, 6553600, 6658560, 6400000, 6658560, 6553600, 6658560}, 
    {6677056, 6765201, 6658560, 6765201, 6677056, 6765201, 6658560, 6765201}, 
    {6553600, 6677056, 6400000, 6677056, 6553600, 6677056, 6400000, 6677056}, 
    {6677056, 6765201, 6658560, 6765201, 6677056, 6765201, 6658560, 6765201}, 
    {6400000, 6658560, 6553600, 6658560, 6400000, 6658560, 6553600, 6658560}, 
    {6677056, 6765201, 6658560, 6765201, 6677056, 6765201, 6658560, 6765201} 
  },
  {
    {7929856, 8156736, 8028160, 8156736, 7929856, 8156736, 8028160, 8156736}, 
    {8156736, 7537770, 7814560, 7537770, 8156736, 7537770, 7814560, 7537770}, 
    {8028160, 7814560, 7840000, 7814560, 8028160, 7814560, 7840000, 7814560}, 
    {8156736, 7537770, 7814560, 7537770, 8156736, 7537770, 7814560, 7537770}, 
    {7929856, 8156736, 8028160, 8156736, 7929856, 8156736, 8028160, 8156736}, 
    {8156736, 7537770, 7814560, 7537770, 8156736, 7537770, 7814560, 7537770}, 
    {8028160, 7814560, 7840000, 7814560, 8028160, 7814560, 7840000, 7814560}, 
    {8156736, 7537770, 7814560, 7537770, 8156736, 7537770, 7814560, 7537770} 
  },
  {
    {11075584, 10653696, 11151360, 10653696, 11075584, 10653696, 11151360, 10653696}, 
    {10653696, 11045652, 11109160, 11045652, 10653696, 11045652, 11109160, 11045652}, 
    {11151360, 11109160, 11289600, 11109160, 11151360, 11109160, 11289600, 11109160}, 
    {10653696, 11045652, 11109160, 11045652, 10653696, 11045652, 11109160, 11045652}, 
    {11075584, 10653696, 11151360, 10653696, 11075584, 10653696, 11151360, 10653696}, 
    {10653696, 11045652, 11109160, 11045652, 10653696, 11045652, 11109160, 11045652}, 
    {11151360, 11109160, 11289600, 11109160, 11151360, 11109160, 11289600, 11109160}, 
    {10653696, 11045652, 11109160, 11045652, 10653696, 11045652, 11109160, 11045652} 
  },
  {
    {12845056, 12503296, 12544000, 12503296, 12845056, 12503296, 12544000, 12503296}, 
    {12503296, 13050156, 12588840, 13050156, 12503296, 13050156, 12588840, 13050156}, 
    {12544000, 12588840, 12960000, 12588840, 12544000, 12588840, 12960000, 12588840}, 
    {12503296, 13050156, 12588840, 13050156, 12503296, 13050156, 12588840, 13050156}, 
    {12845056, 12503296, 12544000, 12503296, 12845056, 12503296, 12544000, 12503296}, 
    {12503296, 13050156, 12588840, 13050156, 12503296, 13050156, 12588840, 13050156}, 
    {12544000, 12588840, 12960000, 12588840, 12544000, 12588840, 12960000, 12588840}, 
    {12503296, 13050156, 12588840, 13050156, 12503296, 13050156, 12588840, 13050156} 
  },
  {
    {16777216, 16646400, 16384000, 16646400, 16777216, 16646400, 16384000, 16646400}, 
    {16646400, 16370116, 16692640, 16370116, 16646400, 16370116, 16692640, 16370116}, 
    {16384000, 16692640, 16646400, 16692640, 16384000, 16692640, 16646400, 16692640}, 
    {16646400, 16370116, 16692640, 16370116, 16646400, 16370116, 16692640, 16370116}, 
    {16777216, 16646400, 16384000, 16646400, 16777216, 16646400, 16384000, 16646400}, 
    {16646400, 16370116, 16692640, 16370116, 16646400, 16370116, 16692640, 16370116}, 
    {16384000, 16692640, 16646400, 16692640, 16384000, 16692640, 16646400, 16692640}, 
    {16646400, 16370116, 16692640, 16370116, 16646400, 16370116, 16692640, 16370116} 
  },
  {
    {21233664, 21381376, 21667840, 21381376, 21233664, 21381376, 21667840, 21381376}, 
    {21381376, 21381376, 21374440, 21381376, 21381376, 21381376, 21374440, 21381376}, 
    {21667840, 21374440, 21529600, 21374440, 21667840, 21374440, 21529600, 21374440}, 
    {21381376, 21381376, 21374440, 21381376, 21381376, 21381376, 21374440, 21381376}, 
    {21233664, 21381376, 21667840, 21381376, 21233664, 21381376, 21667840, 21381376}, 
    {21381376, 21381376, 21374440, 21381376, 21381376, 21381376, 21374440, 21381376}, 
    {21667840, 21374440, 21529600, 21374440, 21667840, 21374440, 21529600, 21374440}, 
    {21381376, 21381376, 21374440, 21381376, 21381376, 21381376, 21374440, 21381376} 
  }
};
#endif

/*! 
 *************************************************************************************
 * \brief
 *    8x8 Intra mode decision for a macroblock
 *************************************************************************************
 */

int Mode_Decision_for_new_Intra8x8Macroblock (double lambda, int *min_cost)
{
  int  cbp=0, b8, cost8x8;

  *min_cost = (int)floor(6.0 * lambda + 0.4999);

  for (b8=0; b8<4; b8++)
  {
    if (Mode_Decision_for_new_8x8IntraBlocks (b8, lambda, &cost8x8))
    {
      cbp |= (1<<b8);
    }
    *min_cost += cost8x8;
  }

  return cbp;
}

/*! 
 *************************************************************************************
 * \brief
 *    8x8 Intra mode decision for a macroblock
 *************************************************************************************
 */

int Mode_Decision_for_new_8x8IntraBlocks (int b8, double lambda, int *min_cost)
{
  int     ipmode, best_ipmode = 0, i, j, k, x, y, cost, dummy;
  int     c_nz, nonzero = 0, diff[64];
  imgpel  rec8x8[8][8];
  double  rdcost = 0.0;
  int     block4x4_x, block4x4_y;
  int     block_x     = 8*(b8 & 0x01);
  int     block_y     = 8*(b8 >> 1);
  int     pic_pix_x   = img->pix_x+block_x;
  int     pic_pix_y   = img->pix_y+block_y;
  int     pic_opix_x   = img->opix_x+block_x;
  int     pic_opix_y   = img->opix_y+block_y;
  int     pic_block_x = pic_pix_x/4;
  int     pic_block_y = pic_pix_y/4;
  double  min_rdcost  = 1e30;
  imgpel    **imgY_orig  = imgY_org;
  extern  int ****cofAC8x8; 
  int fadjust8x8[2][16][16];
  int left_available, up_available, all_available;

  char   upMode;
  char   leftMode;
  int     mostProbableMode;  

  PixelPos left_block;
  PixelPos top_block;

  // Residue Color Transform
  int residue_R, residue_G, residue_B;
  int rate, temp, b4;
  int64 distortion;
  Macroblock     *currMB       = &img->mb_data[img->current_mb_nr];
  int c_ipmode = currMB->c_ipred_mode;
  int rec8x8_c[2][4][4][4];
#ifdef USE_INTRA_MDDT
  long quant_stat_best8x8[64];
#endif

  getLuma4x4Neighbour(img->current_mb_nr, block_x/4, block_y/4, -1,  0, &left_block);
  getLuma4x4Neighbour(img->current_mb_nr, block_x/4, block_y/4,  0, -1, &top_block);

  if (input->UseConstrainedIntraPred)
  {
    top_block.available  = top_block.available ? img->intra_block [top_block.mb_addr] : 0;
    left_block.available = left_block.available ? img->intra_block [left_block.mb_addr] : 0;
  }

  if(b8 >> 1)
    upMode    =  top_block.available ? img->ipredmode8x8[top_block.pos_y ][top_block.pos_x ] : -1; 
  else
    upMode    =  top_block.available ? img->ipredmode   [top_block.pos_y ][top_block.pos_x ] : -1;
  if(b8 & 0x01)
    leftMode  = left_block.available ? img->ipredmode8x8[left_block.pos_y][left_block.pos_x] : -1;
  else
    leftMode  = left_block.available ? img->ipredmode[left_block.pos_y][left_block.pos_x] : -1;

  mostProbableMode  = (upMode < 0 || leftMode < 0) ? DC_PRED : upMode < leftMode ? upMode : leftMode;

  *min_cost = INT_MAX;

  //===== INTRA PREDICTION FOR 8x8 BLOCK =====
  intrapred_luma8x8 (pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);

  //===== LOOP OVER ALL 8x8 INTRA PREDICTION MODES =====
  for (ipmode=0; ipmode<NO_INTRA_PMODE; ipmode++)
  {

    if( (ipmode==DC_PRED) ||
        ((ipmode==VERT_PRED||ipmode==VERT_LEFT_PRED||ipmode==DIAG_DOWN_LEFT_PRED) && up_available ) ||
        ((ipmode==HOR_PRED||ipmode==HOR_UP_PRED) && left_available ) ||
        (all_available) )
    {
      if (!input->rdopt)
      {
        for (k=j=0; j<8; j++)
          for (i=0; i<8; i++, k++)
          {
            diff[k] = imgY_orig[pic_opix_y+j][pic_opix_x+i] - img->mprr_3[ipmode][j][i];
          }
        cost  = (ipmode == mostProbableMode) ? 0 : (int)floor(4 * lambda );
        cost += SATD8X8 (diff, input->hadamard);
        if (cost < *min_cost)
        {
          best_ipmode = ipmode;
          *min_cost   = cost;
        }
      }
      else
      {
        // Residue Color Transform
        if(!img->residue_transform_flag)
        {
          // get prediction and prediction error
          for (j=0; j<8; j++)
          {
            memcpy(&img->mpr[block_y+j][block_x],img->mprr_3[ipmode][j], 8 * sizeof(imgpel));
            for (i=0; i<8; i++)
            {
              img->m7[j][i] = imgY_orig[pic_opix_y+j][pic_opix_x+i] - img->mprr_3[ipmode][j][i];
            }
          }
          //===== store the coding state =====
          //store_coding_state_cs_cm();
          // get and check rate-distortion cost
          
          if ((rdcost = RDCost_for_8x8IntraBlocks (&c_nz, b8, ipmode, lambda, min_rdcost, mostProbableMode)) < min_rdcost)
          {
#ifdef USE_INTRA_MDDT
            if(input->UseIntraMDDT)
              memcpy(quant_stat_best8x8, quant_stat_rd8x8, 64*sizeof(long));
#endif
            //--- set coefficients ---
            for(k=0; k<4; k++) // do 4x now
            {
              for (j=0; j<2; j++)
                memcpy(cofAC8x8[b8][k][j],img->cofAC[b8][k][j], 65 * sizeof(int));
            }
            
            //--- set reconstruction ---
            for (y=0; y<8; y++)
            {
              memcpy(rec8x8[y],&enc_picture->imgY[pic_pix_y+y][pic_pix_x], 8 * sizeof(imgpel));
            }

            if (img->AdaptiveRounding)
            {
              for (j=block_y; j<block_y + 8; j++)
                memcpy(&fadjust8x8[1][j][block_x],&img->fadjust8x8[1][j][block_x], 8 * sizeof(int));
            }
            
            //--- flag if dct-coefficients must be coded ---
            nonzero = c_nz;
            
            //--- set best mode update minimum cost ---
            min_rdcost  = rdcost;
            best_ipmode = ipmode;
          }
          reset_coding_state_cs_cm();  
        }
        else
        {

          for (j=0; j<8; j++)
          for (i=0; i<8; i++)
          {
            residue_B = imgUV_org[0][pic_opix_y+j][pic_opix_x+i] - img->mprr_c[0][c_ipmode][block_y+j][block_x+i];
            residue_G = imgY_org[pic_opix_y+j][pic_opix_x+i] - img->mprr_3[ipmode][j][i];
            residue_R = imgUV_org[1][pic_opix_y+j][pic_opix_x+i] - img->mprr_c[1][c_ipmode][block_y+j][block_x+i];

            /* Forward Residue Transform */
            resTrans_R[j][i] = residue_R-residue_B;
            temp = residue_B+(resTrans_R[j][i]>>1);
            resTrans_B[j][i] = residue_G-temp;
            resTrans_G[j][i] = temp+(resTrans_B[j][i]>>1);
          }

          for (j=0; j<8; j++)
          for (i=0; i<8; i++)
          {
            img->m7[j][i]  = resTrans_G[j][i];
          }

          //store_coding_state_cs_cm();
          rate = (int) RDCost_for_8x8IntraBlocks (&c_nz, b8, ipmode, lambda, min_rdcost, mostProbableMode);
          reset_coding_state_cs_cm();

          for (j=0; j<8; j++)
            for (i=0; i<8; i++)
            {
              rec_resG[j][i] = img->m7[j][i];            
            }


          for(b4=0;b4<4;b4++)
          {
            
            block4x4_x = 4*(b4 & 0x01);
            block4x4_y = 4*(b4 >> 1);
            
            for (j=0; j<4; j++)
              for (i=0; i<4; i++)
              {
                img->m7[j][i]  = resTrans_B[j+block4x4_y][i+block4x4_x];
              }
            rate += RDCost_for_4x4Blocks_Chroma (b8+4, b4, 0);
            for (j=0; j<4; j++)
              for (i=0; i<4; i++)
              {
                rec_resB[j+block4x4_y][i+block4x4_x] = img->m7[j][i];
                img->m7[j][i]  = resTrans_R[j+block4x4_y][i+block4x4_x];
              }
            rate += RDCost_for_4x4Blocks_Chroma (b8+8, b4, 1);
            for (j=0; j<4; j++)
              for (i=0; i<4; i++)
              {
                rec_resR[j+block4x4_y][i+block4x4_x] = img->m7[j][i];
              }
          }
          reset_coding_state_cs_cm();

          for (j=0; j<8; j++)
          {
            for (i=0; i<8; i++)
            {
              /* Inverse Residue Transform */
              temp      = rec_resG[j][i]-(rec_resB[j][i]>>1);
              residue_G = rec_resB[j][i]+temp;
              residue_B = temp - (rec_resR[j][i]>>1);
              residue_R = residue_B+rec_resR[j][i];
              enc_picture->imgUV[0][pic_pix_y+j][pic_pix_x+i] = min(img->max_imgpel_value_uv,max(0,residue_B+(int)img->mprr_c[0][c_ipmode][block_y+j][block_x+i]));
              enc_picture->imgY[pic_pix_y+j][pic_pix_x+i]     = min(img->max_imgpel_value,max(0,residue_G+(int)img->mprr_3[ipmode][j][i]));
              enc_picture->imgUV[1][pic_pix_y+j][pic_pix_x+i] = min(img->max_imgpel_value_uv,max(0,residue_R+(int)img->mprr_c[1][c_ipmode][block_y+j][block_x+i]));
            }
          }
          //===== get distortion (SSD) of 8x8 block =====
          distortion = 0;
          for (y=0; y<8; y++)
            for (x=pic_pix_x; x<pic_pix_x+8; x++)
            {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
              distortion += SQR_DEPTH(imgY_org    [pic_pix_y+y][x], enc_picture->imgY    [pic_pix_y+y][x], input->InputBitDepth, img->BitDepthIncrease);
              distortion += SQR_DEPTH(imgUV_org[0][pic_pix_y+y][x], enc_picture->imgUV[0][pic_pix_y+y][x], input->InputBitDepth, img->BitDepthIncreaseChroma);
              distortion += SQR_DEPTH(imgUV_org[1][pic_pix_y+y][x], enc_picture->imgUV[1][pic_pix_y+y][x], input->InputBitDepth, img->BitDepthIncreaseChroma);
#else
              distortion += img->quad[imgY_org    [pic_pix_y+y][x] - enc_picture->imgY    [pic_pix_y+y][x]];
              distortion += img->quad[imgUV_org[0][pic_pix_y+y][x] - enc_picture->imgUV[0][pic_pix_y+y][x]];
              distortion += img->quad[imgUV_org[1][pic_pix_y+y][x] - enc_picture->imgUV[1][pic_pix_y+y][x]];
#endif
            }
          rdcost = (double)distortion + lambda*(double)rate;

          if (rdcost < min_rdcost)
          {
            //--- set coefficients ---
            for (j=0; j<2; j++)
              for (i=0; i<65;i++)  
                for(k=0; k<4; k++) //do 4x now
                  cofAC8x8[b8][k][j][i]=img->cofAC[b8][k][j][i]; //k vs 0

            for(b4=0; b4<4; b4++)
            {
              block4x4_x = 4*(b4 & 0x01);
              block4x4_y = 4*(b4 >> 1);

              for (j=0; j<2; j++)
                for (i=0; i<18;i++)  cofAC8x8_chroma[0][b4][j][i]=img->cofAC[b8+4][b4][j][i];
              for (j=0; j<2; j++)
                for (i=0; i<18;i++)  cofAC8x8_chroma[1][b4][j][i]=img->cofAC[b8+8][b4][j][i];

              for (i=0; i<2; i++)
              { //uv
                dc_level[i][2*(b8 & 0x01)+(b4 & 0x01)][2*(b8 >> 1)+(b4 >> 1)] = dc_level_temp[i][2*(b8 & 0x01)+(b4 & 0x01)][2*(b8 >> 1)+(b4 >> 1)];
                cbp_chroma_block[i][2*(b8 & 0x01)+(b4 & 0x01)][2*(b8 >> 1)+(b4 >> 1)] = cbp_chroma_block_temp[i][2*(b8 & 0x01)+(b4 & 0x01)][2*(b8 >> 1)+(b4 >> 1)];
                //--- set reconstruction ---
                for (y=0; y<4; y++)
                  for (x=0; x<4; x++)  rec8x8_c[i][b4][y][x] = enc_picture->imgUV[i][pic_pix_y+y+block4x4_y][pic_pix_x+x+block4x4_x];
              }
            }

            //--- set reconstruction ---
            for (y=0; y<8; y++)
              for (x=0; x<8; x++)  
                rec8x8[y][x] = enc_picture->imgY[pic_pix_y+y][pic_pix_x+x];

            //--- flag if dct-coefficients must be coded ---
            nonzero = c_nz;

            //--- set best mode update minimum cost ---
            min_rdcost  = rdcost;
            best_ipmode = ipmode;
          }
         }
      }
    }
  }

#ifdef USE_INTRA_MDDT
  if(input->UseIntraMDDT)
    for (j=0; j<64; j++)
    {
      img->quant_stat8x8[j]+=quant_stat_best8x8[j];
    }
#endif

  //===== set intra mode prediction =====
  img->ipredmode8x8[pic_block_y][pic_block_x] = best_ipmode;
  currMB->intra_pred_modes8x8[4*b8] = (mostProbableMode == best_ipmode) 
    ? -1 
    : (best_ipmode < mostProbableMode ? best_ipmode : best_ipmode-1);

  for(j = img->mb_y*4+(b8 >> 1)*2; j < img->mb_y*4+(b8 >> 1)*2 + 2; j++)   //loop 4x4s in the subblock for 8x8 prediction setting
   memset(&img->ipredmode8x8[j][img->mb_x*4+(b8 & 0x01)*2], best_ipmode, 2 * sizeof(char));


  if (!input->rdopt)
  {
    // Residue Color Transform
    if(!img->residue_transform_flag)
    {
      // get prediction and prediction error
      for (j=0; j<8; j++)
      {
        memcpy(&img->mpr[block_y+j][block_x],img->mprr_3[best_ipmode][j], 8 * sizeof(imgpel));
        for (i=0; i<8; i++)
        {
          img->m7[j][i] = imgY_orig[pic_opix_y+j][pic_opix_x+i] - img->mprr_3[best_ipmode][j][i];
        }
      } 
#ifdef USE_INTRA_MDDT
      if(input->UseIntraMDDT)
      {
        nonzero = (best_ipmode == 2) ?
          dct_luma8x8     (b8, &dummy, 1, best_ipmode) : 
          klt_luma8x8_sep_fast (b8, &dummy, best_ipmode);
      }
      else
        nonzero = dct_luma8x8 (b8, &dummy, 1, best_ipmode); 
#else
      nonzero = dct_luma8x8 (b8, &dummy, 1);
#endif
    }
    else 
    {
      for (j=0; j<8; j++)
      {
        for (i=0; i<8; i++)
        {
          img->mpr[block_y+j][block_x+i]  = img->mprr_3[best_ipmode][j][i];
          residue_B = imgUV_org[0][pic_opix_y+j][pic_opix_x+i] - img->mprr_c[0][c_ipmode][block_y+j][block_x+i];
          residue_G = imgY_org[pic_opix_y+j][pic_opix_x+i] - img->mprr_3[best_ipmode][j][i];
          residue_R = imgUV_org[1][pic_opix_y+j][pic_opix_x+i] - img->mprr_c[1][c_ipmode][block_y+j][block_x+i];
          
          /* Forward Residue Transform */
          resTrans_R[j][i] = residue_R-residue_B;
          temp = residue_B+(resTrans_R[j][i]>>1);
          resTrans_B[j][i] = residue_G-temp;
          resTrans_G[j][i] = temp+(resTrans_B[j][i]>>1);
        }
      } 
      for (j=0; j<8; j++)
      {
        for (i=0; i<8; i++)
          
          img->m7[j][i]  = resTrans_G[j][i];
      }
      
#ifdef USE_INTRA_MDDT
      nonzero = dct_luma8x8 (b8, &dummy, 1, c_ipmode);
#else
      nonzero = dct_luma8x8 (b8, &dummy, 1);
#endif
      
      for (j=0; j<8; j++)
      {
        for (i=0; i<8; i++)
          rec_resG[j][i] = img->m7[j][i];
      }
      
      for(b4=0;b4<4;b4++)
      {

        block4x4_x = 4*(b4 & 0x01);
        block4x4_y = 4*(b4 >> 1);
        
        for (j=0; j<4; j++)
        {
          for (i=0; i<4; i++)                
            img->m7[j][i]  = resTrans_B[j+block4x4_y][i+block4x4_x];
        }
        cbp_chroma_block[0][2*(b8 & 0x01)+(b4 & 0x01)][2*(b8 >> 1)+(b4 >> 1)] = dct_chroma4x4 (0, b8+4, b4);
        dc_level        [0][2*(b8 & 0x01)+(b4 & 0x01)][2*(b8 >> 1)+(b4 >> 1)] = dc_level_temp[0][2*(b8 & 0x01)+(b4 & 0x01)][2*(b8 >> 1)+(b4 >> 1)];
        for (j=0; j<4; j++)
        {
          for (i=0; i<4; i++)
          {
            rec_resB[j+block4x4_y][i+block4x4_x] = img->m7[j][i];
            img->m7[j][i]  = resTrans_R[j+block4x4_y][i+block4x4_x];
          }
        }
        cbp_chroma_block[1][2*(b8 & 0x01)+(b4 & 0x01)][2*(b8 >> 1)+(b4 >> 1)] = dct_chroma4x4 (1, b8+8, b4);
        dc_level        [1][2*(b8 & 0x01)+(b4 & 0x01)][2*(b8 >> 1)+(b4 >> 1)] = dc_level_temp[1][2*(b8 & 0x01)+(b4 & 0x01)][2*(b8 >> 1)+(b4 >> 1)];
        for (j=0; j<4; j++)
        {
          for (i=0; i<4; i++)
            rec_resR[j+block4x4_y][i+block4x4_x] = img->m7[j][i];
        }
      }
      
      for (j=0; j<8; j++)
      {
        for (i=0; i<8; i++)
        {
          /* Inverse Residue Transform */
          temp      = rec_resG[j][i]-(rec_resB[j][i]>>1);
          residue_G = rec_resB[j][i]+temp;
          residue_B = temp - (rec_resR[j][i]>>1);
          residue_R = residue_B+rec_resR[j][i];
          enc_picture->imgUV[0][pic_pix_y+j][pic_pix_x+i] = min(img->max_imgpel_value_uv,max(0,residue_B+(int)img->mprr_c[0][c_ipmode][block_y+j][block_x+i]));
          enc_picture->imgY[pic_pix_y+j][pic_pix_x+i]     = min(img->max_imgpel_value,max(0,residue_G+(int)img->mprr_3[best_ipmode][j][i]));
          enc_picture->imgUV[1][pic_pix_y+j][pic_pix_x+i] = min(img->max_imgpel_value_uv,max(0,residue_R+(int)img->mprr_c[1][c_ipmode][block_y+j][block_x+i]));
        }
      }
    }
  }
  else
  {
    //===== restore coefficients =====
    for(k=0; k<4; k++) // do 4x now    
    {
      for (j=0; j<2; j++)
        memcpy(img->cofAC[b8][k][j],cofAC8x8[b8][k][j], 65 * sizeof(int));
    }

    if (img->AdaptiveRounding)
    {
      for (j=0; j<8; j++)
        memcpy(&img->fadjust8x8[1][block_y+j][block_x], &fadjust8x8[1][block_y+j][block_x], 8 * sizeof(int));
    }
    
    // Residue Color Transform
    if(img->residue_transform_flag)
    for(b4=0; b4<4; b4++)
    {
      for (j=0; j<2; j++)
      for (i=0; i<18;i++)  
        img->cofAC[b8+4][b4][j][i]=cofAC8x8_chroma[0][b4][j][i];
      for (j=0; j<2; j++)
      for (i=0; i<18;i++)  
        img->cofAC[b8+8][b4][j][i]=cofAC8x8_chroma[1][b4][j][i];
    }

    //===== restore reconstruction and prediction (needed if single coeffs are removed) =====
    for (y=0; y<8; y++) 
    {
      memcpy(&enc_picture->imgY[pic_pix_y+y][pic_pix_x], rec8x8[y], 8 * sizeof(imgpel));
      memcpy(&  img->mpr[block_y+y][block_x], img->mprr_3[best_ipmode][y], 8 * sizeof(imgpel));
    }

    // Residue Color Transform
      if(img->residue_transform_flag)
      {
        for(b4=0; b4<4; b4++)
        {
          block4x4_x = 4*(b4 & 0x01);
          block4x4_y = 4*(b4>>1);
          for (i=0; i<2; i++)
          { //uv
            //--- set reconstruction ---
            for (y=0; y<4; y++)
              memcpy(&enc_picture->imgUV[i][pic_pix_y + block4x4_y + y][pic_pix_x + block4x4_x], rec8x8_c[i][b4][y], 4 * sizeof(imgpel));
          }
        }
      }
  }

  return nonzero;
}



// Notation for comments regarding prediction and predictors.
// The pels of the 4x4 block are labelled a..p. The predictor pels above
// are labelled A..H, from the left I..P, and from above left X, as follows:
//
//  Z  A  B  C  D  E  F  G  H  I  J  K  L  M   N  O  P  
//  Q  a1 b1 c1 d1 e1 f1 g1 h1
//  R  a2 b2 c2 d2 e2 f2 g2 h2
//  S  a3 b3 c3 d3 e3 f3 g3 h3
//  T  a4 b4 c4 d4 e4 f4 g4 h4
//  U  a5 b5 c5 d5 e5 f5 g5 h5
//  V  a6 b6 c6 d6 e6 f6 g6 h6
//  W  a7 b7 c7 d7 e7 f7 g7 h7
//  X  a8 b8 c8 d8 e8 f8 g8 h8


// Predictor array index definitions
#define P_Z (PredPel[0])
#define P_A (PredPel[1])
#define P_B (PredPel[2])
#define P_C (PredPel[3])
#define P_D (PredPel[4])
#define P_E (PredPel[5])
#define P_F (PredPel[6])
#define P_G (PredPel[7])
#define P_H (PredPel[8])
#define P_I (PredPel[9])
#define P_J (PredPel[10])
#define P_K (PredPel[11])
#define P_L (PredPel[12])
#define P_M (PredPel[13])
#define P_N (PredPel[14])
#define P_O (PredPel[15])
#define P_P (PredPel[16])
#define P_Q (PredPel[17])
#define P_R (PredPel[18])
#define P_S (PredPel[19])
#define P_T (PredPel[20])
#define P_U (PredPel[21])
#define P_V (PredPel[22])
#define P_W (PredPel[23])
#define P_X (PredPel[24])

/*!
 ************************************************************************
 * \brief
 *    Make intra 8x8 prediction according to all 9 prediction modes.
 *    The routine uses left and upper neighbouring points from
 *    previous coded blocks to do this (if available). Notice that
 *    inaccessible neighbouring points are signalled with a negative
 *    value in the predmode array .
 *
 *  \par Input:
 *     Starting point of current 8x8 block image posision
 *
 *  \par Output:
 *      none
 ************************************************************************
 */
void intrapred_luma8x8(int img_x,int img_y, int *left_available, int *up_available, int *all_available)
{
  int i,j;
  int s0;
  int PredPel[25];  // array of predictor pels
  imgpel **imgY = enc_picture->imgY;  // For MB level frame/field coding tools -- set default to imgY

  int ioff = (img_x & 15);
  int joff = (img_y & 15);
  int mb_nr=img->current_mb_nr;

  PixelPos pix_a[8];
  PixelPos pix_b, pix_c, pix_d;

  int block_available_up;
  int block_available_left;
  int block_available_up_left;
  int block_available_up_right;

  for (i=0;i<8;i++)
  {
    getNeighbour(mb_nr, ioff -1 , joff +i , 1, &pix_a[i]);
  }

  getNeighbour(mb_nr, ioff    , joff -1 , 1, &pix_b);
  getNeighbour(mb_nr, ioff +8 , joff -1 , 1, &pix_c);
  getNeighbour(mb_nr, ioff -1 , joff -1 , 1, &pix_d);
  
  pix_c.available = pix_c.available &&!(ioff == 8 && joff == 8);

  if (input->UseConstrainedIntraPred)
  {
    for (i=0, block_available_left=1; i<8;i++)
      block_available_left  &= pix_a[i].available ? img->intra_block[pix_a[i].mb_addr]: 0;
    block_available_up       = pix_b.available ? img->intra_block [pix_b.mb_addr] : 0;
    block_available_up_right = pix_c.available ? img->intra_block [pix_c.mb_addr] : 0;
    block_available_up_left  = pix_d.available ? img->intra_block [pix_d.mb_addr] : 0;
  }
  else
  {
    block_available_left     = pix_a[0].available;
    block_available_up       = pix_b.available;
    block_available_up_right = pix_c.available;
    block_available_up_left  = pix_d.available;
  }

  *left_available = block_available_left;
  *up_available   = block_available_up;
  *all_available  = block_available_up && block_available_left && block_available_up_left;

  i = (img_x & 15);
  j = (img_y & 15);

  // form predictor pels
  // form predictor pels
  if (block_available_up)
  {
    P_A = imgY[pix_b.pos_y][pix_b.pos_x+0];
    P_B = imgY[pix_b.pos_y][pix_b.pos_x+1];
    P_C = imgY[pix_b.pos_y][pix_b.pos_x+2];
    P_D = imgY[pix_b.pos_y][pix_b.pos_x+3];
    P_E = imgY[pix_b.pos_y][pix_b.pos_x+4];
    P_F = imgY[pix_b.pos_y][pix_b.pos_x+5];
    P_G = imgY[pix_b.pos_y][pix_b.pos_x+6];
    P_H = imgY[pix_b.pos_y][pix_b.pos_x+7];
  }
  else
  {
    P_A = P_B = P_C = P_D = P_E = P_F = P_G = P_H = img->dc_pred_value_luma;
  }

  if (block_available_up_right)
  {
    P_I = imgY[pix_c.pos_y][pix_c.pos_x+0];
    P_J = imgY[pix_c.pos_y][pix_c.pos_x+1];
    P_K = imgY[pix_c.pos_y][pix_c.pos_x+2];
    P_L = imgY[pix_c.pos_y][pix_c.pos_x+3];
    P_M = imgY[pix_c.pos_y][pix_c.pos_x+4];
    P_N = imgY[pix_c.pos_y][pix_c.pos_x+5];
    P_O = imgY[pix_c.pos_y][pix_c.pos_x+6];
    P_P = imgY[pix_c.pos_y][pix_c.pos_x+7];

  }
  else
  {
    P_I = P_J = P_K = P_L = P_M = P_N = P_O = P_P = P_H;
  }

  if (block_available_left)
  {
    P_Q = imgY[pix_a[0].pos_y][pix_a[0].pos_x];
    P_R = imgY[pix_a[1].pos_y][pix_a[1].pos_x];
    P_S = imgY[pix_a[2].pos_y][pix_a[2].pos_x];
    P_T = imgY[pix_a[3].pos_y][pix_a[3].pos_x];
    P_U = imgY[pix_a[4].pos_y][pix_a[4].pos_x];
    P_V = imgY[pix_a[5].pos_y][pix_a[5].pos_x];
    P_W = imgY[pix_a[6].pos_y][pix_a[6].pos_x];
    P_X = imgY[pix_a[7].pos_y][pix_a[7].pos_x];
  }
  else
  {
    P_Q = P_R = P_S = P_T = P_U = P_V = P_W = P_X = img->dc_pred_value_luma;
  }

  if (block_available_up_left)
  {
    P_Z = imgY[pix_d.pos_y][pix_d.pos_x];
  }
  else
  {
    P_Z = img->dc_pred_value_luma;
  }

  for(i=0;i<9;i++)
    img->mprr_3[i][0][0]=-1;

  LowPassForIntra8x8Pred(&(P_Z), block_available_up_left, block_available_up, block_available_left);

  ///////////////////////////////
  // make DC prediction
  ///////////////////////////////
  s0 = 0;
  if (block_available_up && block_available_left)
  {   
    // no edge
    s0 = (P_A + P_B + P_C + P_D + P_E + P_F + P_G + P_H + P_Q + P_R + P_S + P_T + P_U + P_V + P_W + P_X + 8) >> 4;
  }
  else if (!block_available_up && block_available_left)
  {
    // upper edge
    s0 = (P_Q + P_R + P_S + P_T + P_U + P_V + P_W + P_X + 4) >> 3;             
  }
  else if (block_available_up && !block_available_left)
  {
    // left edge
    s0 = (P_A + P_B + P_C + P_D + P_E + P_F + P_G + P_H + 4) >> 3;             
  }
  else //if (!block_available_up && !block_available_left)
  {
    // top left corner, nothing to predict from
    s0 = img->dc_pred_value_luma;                           
  }
  
  // store DC prediction
  for (j=0; j < 2*BLOCK_SIZE; j++)
  {
    for (i=0; i < 2*BLOCK_SIZE; i++)
    {
      img->mprr_3[DC_PRED][i][j] = s0;
    }
  }

  
  ///////////////////////////////
  // make horiz and vert prediction
  ///////////////////////////////

  for (i=0; i < 2*BLOCK_SIZE; i++)
  {
    img->mprr_3[VERT_PRED][0][i] = 
    img->mprr_3[VERT_PRED][1][i] = 
    img->mprr_3[VERT_PRED][2][i] = 
    img->mprr_3[VERT_PRED][3][i] = 
    img->mprr_3[VERT_PRED][4][i] = 
    img->mprr_3[VERT_PRED][5][i] = 
    img->mprr_3[VERT_PRED][6][i] = 
    img->mprr_3[VERT_PRED][7][i] = (&P_A)[i];
    img->mprr_3[HOR_PRED][i][0]  = 
    img->mprr_3[HOR_PRED][i][1]  = 
    img->mprr_3[HOR_PRED][i][2]  = 
    img->mprr_3[HOR_PRED][i][3]  = 
    img->mprr_3[HOR_PRED][i][4]  = 
    img->mprr_3[HOR_PRED][i][5]  = 
    img->mprr_3[HOR_PRED][i][6]  = 
    img->mprr_3[HOR_PRED][i][7]  = (&P_Q)[i];
  }

  if(!block_available_up)img->mprr_3[VERT_PRED][0][0]=-1;
  if(!block_available_left)img->mprr_3[HOR_PRED][0][0]=-1;

  ///////////////////////////////////
  // make diagonal down left prediction
  ///////////////////////////////////
  if (block_available_up) 
  {
    // Mode DIAG_DOWN_LEFT_PRED
    img->mprr_3[DIAG_DOWN_LEFT_PRED][0][0] = (P_A + P_C + 2*(P_B) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_LEFT_PRED][0][1] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][1][0] = (P_B + P_D + 2*(P_C) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_LEFT_PRED][0][2] =
    img->mprr_3[DIAG_DOWN_LEFT_PRED][1][1] =
    img->mprr_3[DIAG_DOWN_LEFT_PRED][2][0] = (P_C + P_E + 2*(P_D) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_LEFT_PRED][0][3] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][1][2] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][2][1] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][3][0] = (P_D + P_F + 2*(P_E) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_LEFT_PRED][0][4] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][1][3] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][2][2] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][3][1] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][4][0] = (P_E + P_G + 2*(P_F) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_LEFT_PRED][0][5] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][1][4] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][2][3] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][3][2] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][4][1] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][5][0] = (P_F + P_H + 2*(P_G) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_LEFT_PRED][0][6] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][1][5] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][2][4] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][3][3] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][4][2] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][5][1] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][6][0] = (P_G + P_I + 2*(P_H) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_LEFT_PRED][0][7] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][1][6] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][2][5] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][3][4] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][4][3] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][5][2] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][6][1] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][7][0] = (P_H + P_J + 2*(P_I) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_LEFT_PRED][1][7] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][2][6] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][3][5] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][4][4] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][5][3] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][6][2] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][7][1] = (P_I + P_K + 2*(P_J) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_LEFT_PRED][2][7] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][3][6] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][4][5] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][5][4] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][6][3] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][7][2] = (P_J + P_L + 2*(P_K) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_LEFT_PRED][3][7] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][4][6] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][5][5] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][6][4] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][7][3] = (P_K + P_M + 2*(P_L) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_LEFT_PRED][4][7] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][5][6] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][6][5] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][7][4] = (P_L + P_N + 2*(P_M) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_LEFT_PRED][5][7] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][6][6] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][7][5] = (P_M + P_O + 2*(P_N) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_LEFT_PRED][6][7] = 
    img->mprr_3[DIAG_DOWN_LEFT_PRED][7][6] = (P_N + P_P + 2*(P_O) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_LEFT_PRED][7][7] = (P_O + 3*(P_P) + 2) >> 2;

    ///////////////////////////////////
    // make vertical left prediction
    ///////////////////////////////////
    img->mprr_3[VERT_LEFT_PRED][0][0] = (P_A + P_B + 1) >> 1;
    img->mprr_3[VERT_LEFT_PRED][0][1] = 
    img->mprr_3[VERT_LEFT_PRED][2][0] = (P_B + P_C + 1) >> 1;
    img->mprr_3[VERT_LEFT_PRED][0][2] = 
    img->mprr_3[VERT_LEFT_PRED][2][1] = 
    img->mprr_3[VERT_LEFT_PRED][4][0] = (P_C + P_D + 1) >> 1;
    img->mprr_3[VERT_LEFT_PRED][0][3] = 
    img->mprr_3[VERT_LEFT_PRED][2][2] = 
    img->mprr_3[VERT_LEFT_PRED][4][1] = 
    img->mprr_3[VERT_LEFT_PRED][6][0] = (P_D + P_E + 1) >> 1;
    img->mprr_3[VERT_LEFT_PRED][0][4] = 
    img->mprr_3[VERT_LEFT_PRED][2][3] = 
    img->mprr_3[VERT_LEFT_PRED][4][2] = 
    img->mprr_3[VERT_LEFT_PRED][6][1] = (P_E + P_F + 1) >> 1;
    img->mprr_3[VERT_LEFT_PRED][0][5] = 
    img->mprr_3[VERT_LEFT_PRED][2][4] = 
    img->mprr_3[VERT_LEFT_PRED][4][3] = 
    img->mprr_3[VERT_LEFT_PRED][6][2] = (P_F + P_G + 1) >> 1;
    img->mprr_3[VERT_LEFT_PRED][0][6] = 
    img->mprr_3[VERT_LEFT_PRED][2][5] = 
    img->mprr_3[VERT_LEFT_PRED][4][4] = 
    img->mprr_3[VERT_LEFT_PRED][6][3] = (P_G + P_H + 1) >> 1;
    img->mprr_3[VERT_LEFT_PRED][0][7] = 
    img->mprr_3[VERT_LEFT_PRED][2][6] = 
    img->mprr_3[VERT_LEFT_PRED][4][5] = 
    img->mprr_3[VERT_LEFT_PRED][6][4] = (P_H + P_I + 1) >> 1;
    img->mprr_3[VERT_LEFT_PRED][2][7] = 
    img->mprr_3[VERT_LEFT_PRED][4][6] = 
    img->mprr_3[VERT_LEFT_PRED][6][5] = (P_I + P_J + 1) >> 1;
    img->mprr_3[VERT_LEFT_PRED][4][7] = 
    img->mprr_3[VERT_LEFT_PRED][6][6] = (P_J + P_K + 1) >> 1;
    img->mprr_3[VERT_LEFT_PRED][6][7] = (P_K + P_L + 1) >> 1;
    img->mprr_3[VERT_LEFT_PRED][1][0] = (P_A + P_C + 2*P_B + 2) >> 2;
    img->mprr_3[VERT_LEFT_PRED][1][1] = 
    img->mprr_3[VERT_LEFT_PRED][3][0] = (P_B + P_D + 2*P_C + 2) >> 2;
    img->mprr_3[VERT_LEFT_PRED][1][2] = 
    img->mprr_3[VERT_LEFT_PRED][3][1] = 
    img->mprr_3[VERT_LEFT_PRED][5][0] = (P_C + P_E + 2*P_D + 2) >> 2;
    img->mprr_3[VERT_LEFT_PRED][1][3] = 
    img->mprr_3[VERT_LEFT_PRED][3][2] = 
    img->mprr_3[VERT_LEFT_PRED][5][1] = 
    img->mprr_3[VERT_LEFT_PRED][7][0] = (P_D + P_F + 2*P_E + 2) >> 2;
    img->mprr_3[VERT_LEFT_PRED][1][4] = 
    img->mprr_3[VERT_LEFT_PRED][3][3] = 
    img->mprr_3[VERT_LEFT_PRED][5][2] = 
    img->mprr_3[VERT_LEFT_PRED][7][1] = (P_E + P_G + 2*P_F + 2) >> 2;
    img->mprr_3[VERT_LEFT_PRED][1][5] = 
    img->mprr_3[VERT_LEFT_PRED][3][4] = 
    img->mprr_3[VERT_LEFT_PRED][5][3] = 
    img->mprr_3[VERT_LEFT_PRED][7][2] = (P_F + P_H + 2*P_G + 2) >> 2;
    img->mprr_3[VERT_LEFT_PRED][1][6] = 
    img->mprr_3[VERT_LEFT_PRED][3][5] = 
    img->mprr_3[VERT_LEFT_PRED][5][4] = 
    img->mprr_3[VERT_LEFT_PRED][7][3] = (P_G + P_I + 2*P_H + 2) >> 2;
    img->mprr_3[VERT_LEFT_PRED][1][7] = 
    img->mprr_3[VERT_LEFT_PRED][3][6] = 
    img->mprr_3[VERT_LEFT_PRED][5][5] = 
    img->mprr_3[VERT_LEFT_PRED][7][4] = (P_H + P_J + 2*P_I + 2) >> 2;
    img->mprr_3[VERT_LEFT_PRED][3][7] = 
    img->mprr_3[VERT_LEFT_PRED][5][6] = 
    img->mprr_3[VERT_LEFT_PRED][7][5] = (P_I + P_K + 2*P_J + 2) >> 2;
    img->mprr_3[VERT_LEFT_PRED][5][7] = 
    img->mprr_3[VERT_LEFT_PRED][7][6] = (P_J + P_L + 2*P_K + 2) >> 2;
    img->mprr_3[VERT_LEFT_PRED][7][7] = (P_K + P_M + 2*P_L + 2) >> 2;
  }

  ///////////////////////////////////
  // make diagonal down right prediction
  ///////////////////////////////////
  if (block_available_up && block_available_left && block_available_up_left) 
  {
    // Mode DIAG_DOWN_RIGHT_PRED
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][7][0] = (P_X + P_V + 2*(P_W) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][6][0] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][7][1] = (P_W + P_U + 2*(P_V) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][5][0] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][6][1] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][7][2] = (P_V + P_T + 2*(P_U) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][4][0] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][5][1] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][6][2] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][7][3] = (P_U + P_S + 2*(P_T) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][3][0] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][4][1] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][5][2] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][6][3] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][7][4] = (P_T + P_R + 2*(P_S) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][2][0] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][3][1] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][4][2] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][5][3] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][6][4] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][7][5] = (P_S + P_Q + 2*(P_R) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][1][0] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][2][1] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][3][2] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][4][3] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][5][4] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][6][5] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][7][6] = (P_R + P_Z + 2*(P_Q) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][0][0] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][1][1] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][2][2] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][3][3] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][4][4] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][5][5] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][6][6] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][7][7] = (P_Q + P_A + 2*(P_Z) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][0][1] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][1][2] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][2][3] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][3][4] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][4][5] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][5][6] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][6][7] = (P_Z + P_B + 2*(P_A) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][0][2] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][1][3] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][2][4] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][3][5] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][4][6] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][5][7] = (P_A + P_C + 2*(P_B) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][0][3] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][1][4] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][2][5] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][3][6] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][4][7] = (P_B + P_D + 2*(P_C) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][0][4] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][1][5] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][2][6] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][3][7] = (P_C + P_E + 2*(P_D) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][0][5] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][1][6] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][2][7] = (P_D + P_F + 2*(P_E) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][0][6] = 
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][1][7] = (P_E + P_G + 2*(P_F) + 2) >> 2;
    img->mprr_3[DIAG_DOWN_RIGHT_PRED][0][7] = (P_F + P_H + 2*(P_G) + 2) >> 2;

  ///////////////////////////////////
  // make vertical right prediction
  ///////////////////////////////////
    img->mprr_3[VERT_RIGHT_PRED][0][0] = 
    img->mprr_3[VERT_RIGHT_PRED][2][1] = 
    img->mprr_3[VERT_RIGHT_PRED][4][2] = 
    img->mprr_3[VERT_RIGHT_PRED][6][3] = (P_Z + P_A + 1) >> 1;
    img->mprr_3[VERT_RIGHT_PRED][0][1] = 
    img->mprr_3[VERT_RIGHT_PRED][2][2] = 
    img->mprr_3[VERT_RIGHT_PRED][4][3] = 
    img->mprr_3[VERT_RIGHT_PRED][6][4] = (P_A + P_B + 1) >> 1;
    img->mprr_3[VERT_RIGHT_PRED][0][2] = 
    img->mprr_3[VERT_RIGHT_PRED][2][3] = 
    img->mprr_3[VERT_RIGHT_PRED][4][4] = 
    img->mprr_3[VERT_RIGHT_PRED][6][5] = (P_B + P_C + 1) >> 1;
    img->mprr_3[VERT_RIGHT_PRED][0][3] = 
    img->mprr_3[VERT_RIGHT_PRED][2][4] = 
    img->mprr_3[VERT_RIGHT_PRED][4][5] = 
    img->mprr_3[VERT_RIGHT_PRED][6][6] = (P_C + P_D + 1) >> 1;
    img->mprr_3[VERT_RIGHT_PRED][0][4] = 
    img->mprr_3[VERT_RIGHT_PRED][2][5] = 
    img->mprr_3[VERT_RIGHT_PRED][4][6] = 
    img->mprr_3[VERT_RIGHT_PRED][6][7] = (P_D + P_E + 1) >> 1;
    img->mprr_3[VERT_RIGHT_PRED][0][5] = 
    img->mprr_3[VERT_RIGHT_PRED][2][6] = 
    img->mprr_3[VERT_RIGHT_PRED][4][7] = (P_E + P_F + 1) >> 1;
    img->mprr_3[VERT_RIGHT_PRED][0][6] = 
    img->mprr_3[VERT_RIGHT_PRED][2][7] = (P_F + P_G + 1) >> 1;
    img->mprr_3[VERT_RIGHT_PRED][0][7] = (P_G + P_H + 1) >> 1;
    img->mprr_3[VERT_RIGHT_PRED][1][0] = 
    img->mprr_3[VERT_RIGHT_PRED][3][1] = 
    img->mprr_3[VERT_RIGHT_PRED][5][2] = 
    img->mprr_3[VERT_RIGHT_PRED][7][3] = (P_Q + P_A + 2*P_Z + 2) >> 2;
    img->mprr_3[VERT_RIGHT_PRED][1][1] = 
    img->mprr_3[VERT_RIGHT_PRED][3][2] = 
    img->mprr_3[VERT_RIGHT_PRED][5][3] = 
    img->mprr_3[VERT_RIGHT_PRED][7][4] = (P_Z + P_B + 2*P_A + 2) >> 2;
    img->mprr_3[VERT_RIGHT_PRED][1][2] = 
    img->mprr_3[VERT_RIGHT_PRED][3][3] = 
    img->mprr_3[VERT_RIGHT_PRED][5][4] = 
    img->mprr_3[VERT_RIGHT_PRED][7][5] = (P_A + P_C + 2*P_B + 2) >> 2;
    img->mprr_3[VERT_RIGHT_PRED][1][3] = 
    img->mprr_3[VERT_RIGHT_PRED][3][4] = 
    img->mprr_3[VERT_RIGHT_PRED][5][5] = 
    img->mprr_3[VERT_RIGHT_PRED][7][6] = (P_B + P_D + 2*P_C + 2) >> 2;
    img->mprr_3[VERT_RIGHT_PRED][1][4] = 
    img->mprr_3[VERT_RIGHT_PRED][3][5] = 
    img->mprr_3[VERT_RIGHT_PRED][5][6] = 
    img->mprr_3[VERT_RIGHT_PRED][7][7] = (P_C + P_E + 2*P_D + 2) >> 2;
    img->mprr_3[VERT_RIGHT_PRED][1][5] = 
    img->mprr_3[VERT_RIGHT_PRED][3][6] = 
    img->mprr_3[VERT_RIGHT_PRED][5][7] = (P_D + P_F + 2*P_E + 2) >> 2;
    img->mprr_3[VERT_RIGHT_PRED][1][6] = 
    img->mprr_3[VERT_RIGHT_PRED][3][7] = (P_E + P_G + 2*P_F + 2) >> 2;
    img->mprr_3[VERT_RIGHT_PRED][1][7] = (P_F + P_H + 2*P_G + 2) >> 2;
    img->mprr_3[VERT_RIGHT_PRED][2][0] =
    img->mprr_3[VERT_RIGHT_PRED][4][1] =
    img->mprr_3[VERT_RIGHT_PRED][6][2] = (P_R + P_Z + 2*P_Q + 2) >> 2;
    img->mprr_3[VERT_RIGHT_PRED][3][0] =
    img->mprr_3[VERT_RIGHT_PRED][5][1] =
    img->mprr_3[VERT_RIGHT_PRED][7][2] = (P_S + P_Q + 2*P_R + 2) >> 2;
    img->mprr_3[VERT_RIGHT_PRED][4][0] =
    img->mprr_3[VERT_RIGHT_PRED][6][1] = (P_T + P_R + 2*P_S + 2) >> 2;
    img->mprr_3[VERT_RIGHT_PRED][5][0] =
    img->mprr_3[VERT_RIGHT_PRED][7][1] = (P_U + P_S + 2*P_T + 2) >> 2;
    img->mprr_3[VERT_RIGHT_PRED][6][0] = (P_V + P_T + 2*P_U + 2) >> 2;
    img->mprr_3[VERT_RIGHT_PRED][7][0] = (P_W + P_U + 2*P_V + 2) >> 2;

  ///////////////////////////////////
  // make horizontal down prediction
  ///////////////////////////////////
    
    img->mprr_3[HOR_DOWN_PRED][0][0] = 
    img->mprr_3[HOR_DOWN_PRED][1][2] = 
    img->mprr_3[HOR_DOWN_PRED][2][4] = 
    img->mprr_3[HOR_DOWN_PRED][3][6] = (P_Q + P_Z + 1) >> 1;
    img->mprr_3[HOR_DOWN_PRED][1][0] = 
    img->mprr_3[HOR_DOWN_PRED][2][2] = 
    img->mprr_3[HOR_DOWN_PRED][3][4] = 
    img->mprr_3[HOR_DOWN_PRED][4][6] = (P_R + P_Q + 1) >> 1;
    img->mprr_3[HOR_DOWN_PRED][2][0] = 
    img->mprr_3[HOR_DOWN_PRED][3][2] = 
    img->mprr_3[HOR_DOWN_PRED][4][4] = 
    img->mprr_3[HOR_DOWN_PRED][5][6] = (P_S + P_R + 1) >> 1;
    img->mprr_3[HOR_DOWN_PRED][3][0] = 
    img->mprr_3[HOR_DOWN_PRED][4][2] = 
    img->mprr_3[HOR_DOWN_PRED][5][4] = 
    img->mprr_3[HOR_DOWN_PRED][6][6] = (P_T + P_S + 1) >> 1;
    img->mprr_3[HOR_DOWN_PRED][4][0] = 
    img->mprr_3[HOR_DOWN_PRED][5][2] = 
    img->mprr_3[HOR_DOWN_PRED][6][4] = 
    img->mprr_3[HOR_DOWN_PRED][7][6] = (P_U + P_T + 1) >> 1;
    img->mprr_3[HOR_DOWN_PRED][5][0] = 
    img->mprr_3[HOR_DOWN_PRED][6][2] = 
    img->mprr_3[HOR_DOWN_PRED][7][4] = (P_V + P_U + 1) >> 1;
    img->mprr_3[HOR_DOWN_PRED][6][0] = 
    img->mprr_3[HOR_DOWN_PRED][7][2] = (P_W + P_V + 1) >> 1;
    img->mprr_3[HOR_DOWN_PRED][7][0] = (P_X + P_W + 1) >> 1;
    img->mprr_3[HOR_DOWN_PRED][0][1] =
    img->mprr_3[HOR_DOWN_PRED][1][3] =
    img->mprr_3[HOR_DOWN_PRED][2][5] =
    img->mprr_3[HOR_DOWN_PRED][3][7] = (P_Q + P_A + 2*P_Z + 2) >> 2;
    img->mprr_3[HOR_DOWN_PRED][1][1] =
    img->mprr_3[HOR_DOWN_PRED][2][3] =
    img->mprr_3[HOR_DOWN_PRED][3][5] =
    img->mprr_3[HOR_DOWN_PRED][4][7] = (P_Z + P_R + 2*P_Q + 2) >> 2;
    img->mprr_3[HOR_DOWN_PRED][2][1] =
    img->mprr_3[HOR_DOWN_PRED][3][3] =
    img->mprr_3[HOR_DOWN_PRED][4][5] =
    img->mprr_3[HOR_DOWN_PRED][5][7] = (P_Q + P_S + 2*P_R + 2) >> 2;
    img->mprr_3[HOR_DOWN_PRED][3][1] =
    img->mprr_3[HOR_DOWN_PRED][4][3] =
    img->mprr_3[HOR_DOWN_PRED][5][5] =
    img->mprr_3[HOR_DOWN_PRED][6][7] = (P_R + P_T + 2*P_S + 2) >> 2;
    img->mprr_3[HOR_DOWN_PRED][4][1] =
    img->mprr_3[HOR_DOWN_PRED][5][3] =
    img->mprr_3[HOR_DOWN_PRED][6][5] =
    img->mprr_3[HOR_DOWN_PRED][7][7] = (P_S + P_U + 2*P_T + 2) >> 2;
    img->mprr_3[HOR_DOWN_PRED][5][1] =
    img->mprr_3[HOR_DOWN_PRED][6][3] =
    img->mprr_3[HOR_DOWN_PRED][7][5] = (P_T + P_V + 2*P_U + 2) >> 2;
    img->mprr_3[HOR_DOWN_PRED][6][1] =
    img->mprr_3[HOR_DOWN_PRED][7][3] = (P_U + P_W + 2*P_V + 2) >> 2;
    img->mprr_3[HOR_DOWN_PRED][7][1] = (P_V + P_X + 2*P_W + 2) >> 2;
    img->mprr_3[HOR_DOWN_PRED][0][2] = 
    img->mprr_3[HOR_DOWN_PRED][1][4] = 
    img->mprr_3[HOR_DOWN_PRED][2][6] = (P_Z + P_B + 2*P_A + 2) >> 2;
    img->mprr_3[HOR_DOWN_PRED][0][3] = 
    img->mprr_3[HOR_DOWN_PRED][1][5] = 
    img->mprr_3[HOR_DOWN_PRED][2][7] = (P_A + P_C + 2*P_B + 2) >> 2;
    img->mprr_3[HOR_DOWN_PRED][0][4] = 
    img->mprr_3[HOR_DOWN_PRED][1][6] = (P_B + P_D + 2*P_C + 2) >> 2;
    img->mprr_3[HOR_DOWN_PRED][0][5] = 
    img->mprr_3[HOR_DOWN_PRED][1][7] = (P_C + P_E + 2*P_D + 2) >> 2;
    img->mprr_3[HOR_DOWN_PRED][0][6] = (P_D + P_F + 2*P_E + 2) >> 2;
    img->mprr_3[HOR_DOWN_PRED][0][7] = (P_E + P_G + 2*P_F + 2) >> 2;
  }

  ///////////////////////////////////
  // make horizontal up prediction
  ///////////////////////////////////
  if (block_available_left)
  {
    img->mprr_3[HOR_UP_PRED][0][0] = (P_Q + P_R + 1) >> 1;
    img->mprr_3[HOR_UP_PRED][1][0] =
    img->mprr_3[HOR_UP_PRED][0][2] = (P_R + P_S + 1) >> 1;
    img->mprr_3[HOR_UP_PRED][2][0] =
    img->mprr_3[HOR_UP_PRED][1][2] =
    img->mprr_3[HOR_UP_PRED][0][4] = (P_S + P_T + 1) >> 1;
    img->mprr_3[HOR_UP_PRED][3][0] =
    img->mprr_3[HOR_UP_PRED][2][2] =
    img->mprr_3[HOR_UP_PRED][1][4] =
    img->mprr_3[HOR_UP_PRED][0][6] = (P_T + P_U + 1) >> 1;
    img->mprr_3[HOR_UP_PRED][4][0] =
    img->mprr_3[HOR_UP_PRED][3][2] =
    img->mprr_3[HOR_UP_PRED][2][4] =
    img->mprr_3[HOR_UP_PRED][1][6] = (P_U + P_V + 1) >> 1;
    img->mprr_3[HOR_UP_PRED][5][0] =
    img->mprr_3[HOR_UP_PRED][4][2] =
    img->mprr_3[HOR_UP_PRED][3][4] =
    img->mprr_3[HOR_UP_PRED][2][6] = (P_V + P_W + 1) >> 1;
    img->mprr_3[HOR_UP_PRED][6][0] =
    img->mprr_3[HOR_UP_PRED][5][2] =
    img->mprr_3[HOR_UP_PRED][4][4] =
    img->mprr_3[HOR_UP_PRED][3][6] = (P_W + P_X + 1) >> 1;
    img->mprr_3[HOR_UP_PRED][4][6] =
    img->mprr_3[HOR_UP_PRED][4][7] =
    img->mprr_3[HOR_UP_PRED][5][4] =
    img->mprr_3[HOR_UP_PRED][5][5] =
    img->mprr_3[HOR_UP_PRED][5][6] =
    img->mprr_3[HOR_UP_PRED][5][7] =
    img->mprr_3[HOR_UP_PRED][6][2] =
    img->mprr_3[HOR_UP_PRED][6][3] =
    img->mprr_3[HOR_UP_PRED][6][4] =
    img->mprr_3[HOR_UP_PRED][6][5] =
    img->mprr_3[HOR_UP_PRED][6][6] =
    img->mprr_3[HOR_UP_PRED][6][7] =
    img->mprr_3[HOR_UP_PRED][7][0] =
    img->mprr_3[HOR_UP_PRED][7][1] =
    img->mprr_3[HOR_UP_PRED][7][2] =
    img->mprr_3[HOR_UP_PRED][7][3] =
    img->mprr_3[HOR_UP_PRED][7][4] =
    img->mprr_3[HOR_UP_PRED][7][5] =
    img->mprr_3[HOR_UP_PRED][7][6] =
    img->mprr_3[HOR_UP_PRED][7][7] = P_X;
    img->mprr_3[HOR_UP_PRED][6][1] =
    img->mprr_3[HOR_UP_PRED][5][3] =
    img->mprr_3[HOR_UP_PRED][4][5] =
    img->mprr_3[HOR_UP_PRED][3][7] = (P_W + 3*P_X + 2) >> 2;
    img->mprr_3[HOR_UP_PRED][5][1] =
    img->mprr_3[HOR_UP_PRED][4][3] =
    img->mprr_3[HOR_UP_PRED][3][5] =
    img->mprr_3[HOR_UP_PRED][2][7] = (P_X + P_V + 2*P_W + 2) >> 2;
    img->mprr_3[HOR_UP_PRED][4][1] =
    img->mprr_3[HOR_UP_PRED][3][3] =
    img->mprr_3[HOR_UP_PRED][2][5] =
    img->mprr_3[HOR_UP_PRED][1][7] = (P_W + P_U + 2*P_V + 2) >> 2;
    img->mprr_3[HOR_UP_PRED][3][1] =
    img->mprr_3[HOR_UP_PRED][2][3] =
    img->mprr_3[HOR_UP_PRED][1][5] =
    img->mprr_3[HOR_UP_PRED][0][7] = (P_V + P_T + 2*P_U + 2) >> 2;
    img->mprr_3[HOR_UP_PRED][2][1] =
    img->mprr_3[HOR_UP_PRED][1][3] =
    img->mprr_3[HOR_UP_PRED][0][5] = (P_U + P_S + 2*P_T + 2) >> 2;
    img->mprr_3[HOR_UP_PRED][1][1] =
    img->mprr_3[HOR_UP_PRED][0][3] = (P_T + P_R + 2*P_S + 2) >> 2;
    img->mprr_3[HOR_UP_PRED][0][1] = (P_S + P_Q + 2*P_R + 2) >> 2;
  }
}

/*! 
 *************************************************************************************
 * \brief
 *    Prefiltering for Intra8x8 prediction
 *************************************************************************************
 */
void LowPassForIntra8x8Pred(int *PredPel, int block_up_left, int block_up, int block_left)
{
  int i;
  int LoopArray[25];
 

  for(i = 0; i < 25; i++)
     LoopArray[i] = PredPel[i] ;

   if(block_up)
  {
    if(block_up_left) 
    {
      LoopArray[1] = ((&P_Z)[0] + ((&P_Z)[1]<<1) + (&P_Z)[2] + 2)>>2;
    }
    else
      LoopArray[1] = ((&P_Z)[1] + ((&P_Z)[1]<<1) + (&P_Z)[2] + 2)>>2; 


    for(i = 2; i <16; i++)
    {
      LoopArray[i] = ((&P_Z)[i-1] + ((&P_Z)[i]<<1) + (&P_Z)[i+1] + 2)>>2;
    }
    LoopArray[16] = (P_P + (P_P<<1) + P_O + 2)>>2;
  }

  if(block_up_left) 
  {
    
    if(block_up && block_left)
    {
        LoopArray[0] = (P_Q + (P_Z<<1) + P_A +2)>>2;
    }
    else
    {
      if(block_up)
        LoopArray[0] = (P_Z + (P_Z<<1) + P_A +2)>>2;
      else
        if(block_left)
          LoopArray[0] = (P_Z + (P_Z<<1) + P_Q +2)>>2;
    }

  }

  if(block_left)
  {
    if(block_up_left)
      LoopArray[17] = (P_Z + (P_Q<<1) + P_R + 2)>>2; 
    else
      LoopArray[17] = (P_Q + (P_Q<<1) + P_R + 2)>>2;

    for(i = 18; i <24; i++)
    {
      LoopArray[i] = ((&P_Z)[i-1] + ((&P_Z)[i]<<1) + (&P_Z)[i+1] + 2)>>2;
    }
    LoopArray[24] = (P_W + (P_X<<1) + P_X + 2)>>2;
  }

  for(i = 0; i < 25; i++)
    PredPel[i] = LoopArray[i];
}





/*! 
 *************************************************************************************
 * \brief
 *    R-D Cost for an 8x8 Intra block
 *************************************************************************************
 */

double RDCost_for_8x8IntraBlocks(int *nonzero, int b8, int ipmode, double lambda, double min_rdcost, int mostProbableMode)
{
  double  rdcost = 0.0;
  int     dummy, x, y, rate;
  int64   distortion  = 0;
  int     block_x     = 8*(b8 & 0x01);
  int     block_y     = 8*(b8 >> 1);
  int     pic_pix_x   = img->pix_x+block_x;
  int     pic_pix_y   = img->pix_y+block_y;
  int     pic_opix_y  = img->opix_y+block_y;
  imgpel    **imgY_orig  = imgY_org;
  imgpel    **imgY       = enc_picture->imgY;

  Slice          *currSlice    =  img->currentSlice;
  Macroblock     *currMB       = &img->mb_data[img->current_mb_nr];
  SyntaxElement  *currSE       = &img->MB_SyntaxElements[currMB->currSEnr];
  const int      *partMap      = assignSE2partition[input->partition_mode];
  DataPartition  *dataPart;

  //===== perform DCT, Q, IQ, IDCT, Reconstruction =====
  dummy = 0;

#ifdef USE_INTRA_MDDT
  if(input->UseIntraMDDT)
  {
    *nonzero = (ipmode == 2) ? 
      dct_luma8x8 (b8, &dummy, 1, ipmode) : 
      klt_luma8x8_sep_fast (b8, &dummy, ipmode);
  }
  else
    *nonzero = dct_luma8x8 (b8, &dummy, 1, ipmode) ; 
#else
  *nonzero = dct_luma8x8 (b8, &dummy, 1);
#endif

  //===== get distortion (SSD) of 8x8 block =====
  for (y=0; y<8; y++)
    for (x=pic_pix_x; x<pic_pix_x+8; x++)  
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
      distortion += SQR_DEPTH(imgY_orig[pic_opix_y+y][x], imgY[pic_pix_y+y][x], input->BitDepthLuma, img->BitDepthIncrease);
#else
      distortion += img->quad [imgY_orig[pic_opix_y+y][x] - imgY[pic_pix_y+y][x]];
#endif

  //===== RATE for INTRA PREDICTION MODE  (SYMBOL MODE MUST BE SET TO UVLC) =====
  currSE->value1 = (mostProbableMode == ipmode) ? -1 : ipmode < mostProbableMode ? ipmode : ipmode-1;

  //--- set position and type ---
  currSE->context = b8;
  currSE->type    = SE_INTRAPREDMODE;

  //--- choose data partition ---
  if (img->type!=B_SLICE)
    dataPart = &(currSlice->partArr[partMap[SE_INTRAPREDMODE]]);
  else
    dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);

  //--- encode and update rate ---
  if (input->symbol_mode == UVLC)
    writeSyntaxElement_Intra4x4PredictionMode(currSE, dataPart);
  else
  {
    currSE->writing = writeIntraPredMode_CABAC;
    dataPart->writeSyntaxElement (currSE, dataPart);
  }

  rate = currSE->len;
  currSE++;
  currMB->currSEnr++;

  //===== RATE for LUMINANCE COEFFICIENTS =====

  if (input->symbol_mode == UVLC)
  {
    int b4;
    for(b4=0; b4<4; b4++)
      rate  += writeCoeff4x4_CAVLC (LUMA, b8, b4, 0);
  }
  else
  {
    rate  += writeLumaCoeff8x8_CABAC (b8, 1);
  }

  rdcost = (double)distortion + lambda*(double)rate;

  if(img->residue_transform_flag)
    return (double)rate;
  else
    return rdcost;
}


/*!
 ************************************************************************
 * \brief
 *    The routine performs transform,quantization,inverse transform, adds the diff.
 *    to the prediction and writes the result to the decoded luma frame. Includes the
 *    RD constrained quantization also.
 *
 * \par Input:
 *    b8: Block position inside a macro block (0,1,2,3).
 *
 * \par Output:
 *    nonzero: 0 if no levels are nonzero.  1 if there are nonzero levels.  
 *    coeff_cost: Counter for nonzero coefficients, used to discard expensive levels.
 ************************************************************************
 */

#define MC(coeff) ((coeff)&3)


#ifdef RDO_Q //TREL_CAVLC
extern const byte SNGL_SCAN[16][2];
extern const byte FIELD_SCAN[16][2];
void TrellisCAVLC_Q(int m4[4][4], levelDataStruct *levelData, int *levelTrellis, int block_type, int b8, int b4, int coeff_num, double lambda);

void TrellisCAVLC8x8(int img_m7[16][16], int q_bits, int qp_rem, int **levelscale, int **leveloffset, int *levelTrellis, int b8, double lambda)
{
  int i, j, i4, j4, coeff_ctr, MCcoeff, scan_poss4x4;
  int m4[4][4][4];
  int level;
  int levelTrellis4x4[4][16];
  levelDataStruct levelData[64], levelData4x4[4][16]; 
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
  double  normFact=pow(2,(2*Q_BITS_8+9))*(1<<(2*img->BitDepthIncrease));
#else
  double  normFact=pow(2,(2*Q_BITS_8+9));
#endif
  double err;
  short is_field_mode = (img->field_picture || ( img->MbaffFrameFlag && img->mb_data[img->current_mb_nr].mb_field));
  const byte (*pos_scan4x4)[2] = is_field_mode ? FIELD_SCAN : SNGL_SCAN;
  int lowerInt, ii;

  //if(input->UseRDO_Q)
  {
    for (coeff_ctr=0;coeff_ctr < 64;coeff_ctr++)
    {
      if (is_field_mode) 
      {  // Alternate scan for field coding
        i=FIELD_SCAN8x8[coeff_ctr][0];
        j=FIELD_SCAN8x8[coeff_ctr][1];
      }
      else 
      {
        i=SNGL_SCAN8x8[coeff_ctr][0];
        j=SNGL_SCAN8x8[coeff_ctr][1];
      }

      levelData[coeff_ctr].levelDouble=absm(img_m7[j][i]*levelscale[i][j]);
      level=(int) (levelData[coeff_ctr].levelDouble>>q_bits);

      lowerInt=(((int)levelData[coeff_ctr].levelDouble-(level<<q_bits))<(1<<(q_bits-1)))? 1 : 0;

      levelData[coeff_ctr].level[0]=0;
      if (level==0 && lowerInt==1)
      {
        levelData[coeff_ctr].noLevels=1;
      }
      else if (level==0 && lowerInt==0)
      {
        levelData[coeff_ctr].level[1] = level+1;
        levelData[coeff_ctr].noLevels=2;
      }
      else if (level>0 && lowerInt==1)
      {
        if(level > 1)
        {
          levelData[coeff_ctr].level[1] = level-1;
          levelData[coeff_ctr].level[2] = level;
          levelData[coeff_ctr].noLevels=3;
        }
        else
        {
          levelData[coeff_ctr].level[1] = level;
          levelData[coeff_ctr].noLevels=2;
        }
      }
      else
      {
        levelData[coeff_ctr].level[1] = level;
        levelData[coeff_ctr].level[2] = level+1;
        levelData[coeff_ctr].noLevels=3;
      }

      for (ii=0; ii<levelData[coeff_ctr].noLevels; ii++)
      {
        err=(double)(levelData[coeff_ctr].level[ii]<<q_bits)-(double)levelData[coeff_ctr].levelDouble;
        levelData[coeff_ctr].errLevel[ii]=(err*err*(double)estErr8x8[qp_rem][i][j])/normFact; 
      }

      if(levelData[coeff_ctr].noLevels == 1)
        levelData[coeff_ctr].pre_level = 0;
      else
        levelData[coeff_ctr].pre_level = (absm (img_m7[j][i]) * levelscale[i][j] + leveloffset[i][j]) >> q_bits;
    }
  }

  for (coeff_ctr=0;coeff_ctr < 64;coeff_ctr++)
  {
    if (is_field_mode) 
    {  // Alternate scan for field coding
      i=FIELD_SCAN8x8[coeff_ctr][0];
      j=FIELD_SCAN8x8[coeff_ctr][1];
    }
    else 
    {
      i=SNGL_SCAN8x8[coeff_ctr][0];
      j=SNGL_SCAN8x8[coeff_ctr][1];
    }
	
    MCcoeff = MC(coeff_ctr);    //4x4 block index

    scan_poss4x4 = coeff_ctr>>2; // coeff index in a 4x4 block

    //split an 8x8 block into 4 4x4 blocks
    i4=pos_scan4x4[scan_poss4x4][0];
    j4=pos_scan4x4[scan_poss4x4][1];

    m4[MCcoeff][j4][i4] = img_m7[j][i];
    levelData4x4[MCcoeff][scan_poss4x4] = levelData[coeff_ctr];
  }


  //RDO for each 4x4 block
  for(MCcoeff=0; MCcoeff<4; MCcoeff++)
  {
    TrellisCAVLC_Q(m4[MCcoeff], levelData4x4[MCcoeff], levelTrellis4x4[MCcoeff], LUMA, b8, MCcoeff, 16, lambda);
  }

  //put quantized values into an 8x8 block
  for (coeff_ctr=0;coeff_ctr < 64;coeff_ctr++)
  {
    MCcoeff = MC(coeff_ctr);    
    scan_poss4x4 = coeff_ctr>>2;
    levelTrellis[coeff_ctr] = levelTrellis4x4[MCcoeff][scan_poss4x4];
  }
}
#endif


#ifdef USE_INTRA_MDDT
int dct_luma8x8(int b8,int *coeff_cost, int intra, int ipmode)
#else
int dct_luma8x8(int b8,int *coeff_cost, int intra)
#endif
{
  int sign(int a,int b);

  int i,j,ilev,coeff_ctr;
  int level,scan_pos,run;
  int nonzero;
  int qp_per,qp_rem,q_bits;
  int dq_lshift = 0, dq_rshift = 0, dq_round = 0;

  int block_x = 8*(b8 & 0x01);
  int block_y = 8*(b8 >> 1);
  int*  ACLevel = img->cofAC[b8][0][0];
  int*  ACRun   = img->cofAC[b8][0][1];
  int m6[8][8];
  int a[8], b[8];
  int scan_poss[4],runs[4];
  int pix_x, pix_y, ipix_y;
  int **levelscale,**leveloffset;
  int **invlevelscale;
  int MCcoeff;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  short is_field_mode = (img->field_picture || ( img->MbaffFrameFlag && currMB->mb_field));
#ifdef RDO_Q
  levelDataStruct levelData[64];
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
  double  lambda_md=0, normFact=pow(2,(2*Q_BITS_8+9))*(1<<(2*img->BitDepthIncrease));
#else
  double  lambda_md=0, normFact=pow(2,(2*Q_BITS_8+9));
#endif
  double err;
  int lowerInt, levelTrellis[64], ii, kStart=0, kStop=0, noCoeff;
#endif

  Boolean lossless_qpprime = (Boolean)((img->qp + img->bitdepth_luma_qp_scale)==0 && img->lossless_qpprime_flag==1);
  
  qp_per    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)/6;
  qp_rem    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)%6;
  q_bits    = Q_BITS_8+qp_per;

#ifdef RDO_Q
  if(input->UseRDO_Q)
  {
    if ((img->type==B_SLICE) && img->nal_reference_idc)
    {
      lambda_md = img->lambda_md[5][img->masterQP];  
    }
    else
    {
      lambda_md = img->lambda_md[img->type][img->masterQP];
    }
  }
#endif

#ifdef ADAPTIVE_QUANTIZATION
  if(img->slice_fractional_quant_flag)
  {
    levelscale    = LevelScale8x8Luma_IAQMS   [img->mb_iaqms_idx][intra][qp_rem];
    invlevelscale = InvLevelScale8x8Luma_IAQMS[img->mb_iaqms_idx][intra][qp_rem];
  }
  else
  {
    levelscale    = LevelScale8x8Luma[intra][qp_rem];
    invlevelscale = InvLevelScale8x8Luma[intra][qp_rem];
  }
  leveloffset   = LevelOffset8x8Luma[intra][qp_per];
#else
  levelscale    = LevelScale8x8Luma[intra][qp_rem];
  leveloffset   = LevelOffset8x8Luma[intra][qp_per];
  invlevelscale = InvLevelScale8x8Luma[intra][qp_rem];
#endif

  if (qp_per < 6)
  {
    dq_rshift = 6 - qp_per;
    dq_round  = 1<<(5-qp_per);
  }
  else
    dq_lshift = qp_per - 6;
    

  // horizontal transform
  if (!lossless_qpprime) 
  {
    for( i=0; i<8; i++)
    {
      a[0] = img->m7[i][0] + img->m7[i][7];
      a[1] = img->m7[i][1] + img->m7[i][6];
      a[2] = img->m7[i][2] + img->m7[i][5];
      a[3] = img->m7[i][3] + img->m7[i][4];
      
      b[0] = a[0] + a[3];
      b[1] = a[1] + a[2];
      b[2] = a[0] - a[3];
      b[3] = a[1] - a[2];
      
      a[4] = img->m7[i][0] - img->m7[i][7];
      a[5] = img->m7[i][1] - img->m7[i][6];
      a[6] = img->m7[i][2] - img->m7[i][5];
      a[7] = img->m7[i][3] - img->m7[i][4];
      
      b[4]= a[5] + a[6] + ((a[4]>>1) + a[4]);
      b[5]= a[4] - a[7] - ((a[6]>>1) + a[6]);
      b[6]= a[4] + a[7] - ((a[5]>>1) + a[5]);
      b[7]= a[5] - a[6] + ((a[7]>>1) + a[7]);
      
      m6[0][i] = b[0] + b[1];
      m6[2][i] = b[2] + (b[3]>>1);
      m6[4][i] = b[0] - b[1];
      m6[6][i] = (b[2]>>1) - b[3];
      m6[1][i] =   b[4] + (b[7]>>2);
      m6[3][i] =   b[5] + (b[6]>>2);
      m6[5][i] =   b[6] - (b[5]>>2);
      m6[7][i] = - b[7] + (b[4]>>2);
      
    }
    // vertical transform
    for( i=0; i<8; i++)
    {
      a[0] = m6[i][0] + m6[i][7];
      a[1] = m6[i][1] + m6[i][6];
      a[2] = m6[i][2] + m6[i][5];
      a[3] = m6[i][3] + m6[i][4];
      
      b[0] = a[0] + a[3];
      b[1] = a[1] + a[2];
      b[2] = a[0] - a[3];
      b[3] = a[1] - a[2];
      
      a[4] = m6[i][0] - m6[i][7];
      a[5] = m6[i][1] - m6[i][6];
      a[6] = m6[i][2] - m6[i][5];
      a[7] = m6[i][3] - m6[i][4];
      
      b[4]= a[5] + a[6] + ((a[4]>>1) + a[4]);
      b[5]= a[4] - a[7] - ((a[6]>>1) + a[6]);
      b[6]= a[4] + a[7] - ((a[5]>>1) + a[5]);
      b[7]= a[5] - a[6] + ((a[7]>>1) + a[7]);
      
      img->m7[0][i] = b[0] + b[1];
      img->m7[2][i] = b[2] + (b[3]>>1);
      img->m7[4][i] = b[0] - b[1];
      img->m7[6][i] = (b[2]>>1) - b[3];
      img->m7[1][i] =   b[4] + (b[7]>>2);
      img->m7[3][i] =   b[5] + (b[6]>>2);
      img->m7[5][i] =   b[6] - (b[5]>>2);
      img->m7[7][i] = - b[7] + (b[4]>>2);
    }
  }


#ifdef RDO_Q
//#ifdef TREL_CAVLC
  if(input->UseRDO_Q && active_pps->entropy_coding_mode_flag == UVLC)
    TrellisCAVLC8x8(img->m7, q_bits, qp_rem, levelscale, leveloffset, levelTrellis, b8, lambda_md);
//#endif
  if(input->UseRDO_Q && active_pps->entropy_coding_mode_flag == CABAC)
  {
    noCoeff=0;
    for (coeff_ctr=0;coeff_ctr < 64;coeff_ctr++)
    {
      if (is_field_mode) 
      {  // Alternate scan for field coding
        i=FIELD_SCAN8x8[coeff_ctr][0];
        j=FIELD_SCAN8x8[coeff_ctr][1];
      }
#ifdef USE_INTRA_MDDT
      else if(input->UseIntraMDDT && intra && !img->residue_transform_flag)
      {
        i=img->scanOrder8x8[ipmode][coeff_ctr][0];
        j=img->scanOrder8x8[ipmode][coeff_ctr][1];
      }
#endif       
      else 
      {
        i=SNGL_SCAN8x8[coeff_ctr][0];
        j=SNGL_SCAN8x8[coeff_ctr][1];
      }

      levelData[coeff_ctr].levelDouble = absm(img->m7[j][i] * levelscale[i][j]);
      level = (int)(levelData[coeff_ctr].levelDouble >> q_bits);
      lowerInt = (((int)levelData[coeff_ctr].levelDouble - (level << q_bits)) < (1 << (q_bits - 1)))? 1: 0;

      levelData[coeff_ctr].level[0]=0;
      if (level==0 && lowerInt==1)
      {
        levelData[coeff_ctr].noLevels=1;
      }
      else if (level==0 && lowerInt==0)
      {
        levelData[coeff_ctr].level[1] = level+1;
        levelData[coeff_ctr].noLevels=2;
        kStop=coeff_ctr;
        noCoeff++;
      }
      else if (level>0 && lowerInt==1)
      {
        levelData[coeff_ctr].level[1] = level;
        levelData[coeff_ctr].noLevels=2;
        kStop=coeff_ctr;
        noCoeff++;
      }
      else
      {
        levelData[coeff_ctr].level[1] = level;
        levelData[coeff_ctr].level[2] = level+1;
        levelData[coeff_ctr].noLevels=3;
        kStop=coeff_ctr;
        kStart=coeff_ctr;
        noCoeff++;
      }

      for (ii=0; ii<levelData[coeff_ctr].noLevels; ii++)
      {
        err=(double)(levelData[coeff_ctr].level[ii]<<q_bits)-(double)levelData[coeff_ctr].levelDouble;
        levelData[coeff_ctr].errLevel[ii]=(err*err*(double)estErr8x8[qp_rem][i][j])/normFact; 
      }
    }
    est_writeRunLevel_CABAC(levelData, levelTrellis, LUMA_8x8, lambda_md, kStart, kStop, noCoeff, 0);
  }
#endif

  // Quant  
  nonzero=FALSE;
  
  run=-1;
  scan_pos=0;
  
  runs[0]=runs[1]=runs[2]=runs[3]=-1;
  scan_poss[0]=scan_poss[1]=scan_poss[2]=scan_poss[3]=0;
  
  for (coeff_ctr=0;coeff_ctr < 64;coeff_ctr++)
  {
    
    if (is_field_mode) 
    {  // Alternate scan for field coding
      i=FIELD_SCAN8x8[coeff_ctr][0];
      j=FIELD_SCAN8x8[coeff_ctr][1];
    }
#ifdef USE_INTRA_MDDT
    else if(input->UseIntraMDDT && intra && !img->residue_transform_flag)
    {
      i=img->scanOrder8x8[ipmode][coeff_ctr][0];
      j=img->scanOrder8x8[ipmode][coeff_ctr][1];
    }
#endif 
    else 
    {
      i=SNGL_SCAN8x8[coeff_ctr][0];
      j=SNGL_SCAN8x8[coeff_ctr][1];
    }
    MCcoeff = MC(coeff_ctr);
    run++;
    ilev=0;
    
    runs[MCcoeff]++;
    
#ifdef RDO_Q
    if(input->UseRDO_Q)
      level = levelTrellis[coeff_ctr];
    else
    {
      if(lossless_qpprime)
        level = absm (img->m7[j][i]);
      else 
        level = (absm (img->m7[j][i]) * levelscale[i][j] + leveloffset[i][j]) >> q_bits;
    }
#else
    if(lossless_qpprime)
      level = absm (img->m7[j][i]);
    else 
      level = (absm (img->m7[j][i]) * levelscale[i][j] + leveloffset[i][j]) >> q_bits;
#endif
  
    if (img->AdaptiveRounding)
    {
      if (lossless_qpprime || level == 0 )
      {
        img->fadjust8x8[intra][block_y+j][block_x+i] = 0;
      }
      else 
      {
        img->fadjust8x8[intra][block_y + j][block_x + i] = 
          (AdaptRndWeight * (absm (img->m7[j][i]) * levelscale[i][j] - (level << q_bits)) + (1<< (q_bits))) >> (q_bits + 1);       
      }
    }

    if (level != 0)
    {
      nonzero=TRUE;
      
      if (currMB->luma_transform_size_8x8_flag && input->symbol_mode == UVLC)
      {
        *coeff_cost += (level > 1 || lossless_qpprime) ? MAX_VALUE : COEFF_COST8x8[input->disthres][runs[MCcoeff]];

        img->cofAC[b8][MCcoeff][0][scan_poss[MCcoeff]] = sign(level,img->m7[j][i]);
        img->cofAC[b8][MCcoeff][1][scan_poss[MCcoeff]] = runs[MCcoeff];
        ++scan_poss[MCcoeff];
        runs[MCcoeff]=-1;
      }
      else
      {
        *coeff_cost += (level > 1 || lossless_qpprime) ? MAX_VALUE :COEFF_COST8x8[input->disthres][run];
        ACLevel[scan_pos] = sign(level,img->m7[j][i]);
        ACRun  [scan_pos] = run;
        ++scan_pos;
        run=-1;                     // reset zero level counter
      }      
      
      level = sign(level, img->m7[j][i]);
      if(lossless_qpprime)
      {
        ilev = level;
      }
      else 
      {
        if (qp_per>=6)
          ilev = level*invlevelscale[i][j]<<dq_lshift; // dequantization
        else
          ilev = (level*invlevelscale[i][j] + dq_round)>>dq_rshift; // dequantization
      }
    }
    if(!lossless_qpprime)
      img->m7[j][i] = ilev;
  }
  if (!currMB->luma_transform_size_8x8_flag || input->symbol_mode != UVLC)
    ACLevel[scan_pos] = 0;
  else
  {
    for(i=0; i<4; i++)
      img->cofAC[b8][i][0][scan_poss[i]] = 0;
  }
 
  
  //    Inverse Transform
  // horizontal inverse transform
  if (!lossless_qpprime)
  {
    for( i=0; i<8; i++)
    {
      a[0] = img->m7[i][0] + img->m7[i][4];
      a[4] = img->m7[i][0] - img->m7[i][4];
      a[2] = (img->m7[i][2]>>1) - img->m7[i][6];
      a[6] = img->m7[i][2] + (img->m7[i][6]>>1);
      
      b[0] = a[0] + a[6];
      b[2] = a[4] + a[2];
      b[4] = a[4] - a[2];
      b[6] = a[0] - a[6];
      
      a[1] = -img->m7[i][3] + img->m7[i][5] - img->m7[i][7] - (img->m7[i][7]>>1);
      a[3] =  img->m7[i][1] + img->m7[i][7] - img->m7[i][3] - (img->m7[i][3]>>1);
      a[5] = -img->m7[i][1] + img->m7[i][7] + img->m7[i][5] + (img->m7[i][5]>>1);
      a[7] =  img->m7[i][3] + img->m7[i][5] + img->m7[i][1] + (img->m7[i][1]>>1);
      
      b[1] = a[1] + (a[7]>>2);
      b[7] = -(a[1]>>2) + a[7];
      b[3] = a[3] + (a[5]>>2);
      b[5] = (a[3]>>2) - a[5];
      
      m6[0][i] = b[0] + b[7];
      m6[1][i] = b[2] + b[5];
      m6[2][i] = b[4] + b[3];
      m6[3][i] = b[6] + b[1];
      m6[4][i] = b[6] - b[1];
      m6[5][i] = b[4] - b[3];
      m6[6][i] = b[2] - b[5];
      m6[7][i] = b[0] - b[7];
    }
    
    // vertical inverse transform
    for( i=0; i<8; i++)
    {
      a[0] =  m6[i][0] + m6[i][4];
      a[4] =  m6[i][0] - m6[i][4];
      a[2] = (m6[i][2]>>1) - m6[i][6];
      a[6] =  m6[i][2] + (m6[i][6]>>1);
      
      b[0] = a[0] + a[6];
      b[2] = a[4] + a[2];
      b[4] = a[4] - a[2];
      b[6] = a[0] - a[6];
      
      a[1] = -m6[i][3] + m6[i][5] - m6[i][7] - (m6[i][7]>>1);
      a[3] =  m6[i][1] + m6[i][7] - m6[i][3] - (m6[i][3]>>1);
      a[5] = -m6[i][1] + m6[i][7] + m6[i][5] + (m6[i][5]>>1);
      a[7] =  m6[i][3] + m6[i][5] + m6[i][1] + (m6[i][1]>>1);
      
      b[1] =   a[1] + (a[7]>>2);
      b[7] = -(a[1]>>2) + a[7];
      b[3] =   a[3] + (a[5]>>2);
      b[5] =  (a[3]>>2) - a[5];
      
      img->m7[0][i] = b[0] + b[7];
      img->m7[1][i] = b[2] + b[5];
      img->m7[2][i] = b[4] + b[3];
      img->m7[3][i] = b[6] + b[1];
      img->m7[4][i] = b[6] - b[1];
      img->m7[5][i] = b[4] - b[3];
      img->m7[6][i] = b[2] - b[5];
      img->m7[7][i] = b[0] - b[7];
    }
  }
  
  if (!img->residue_transform_flag)
  {
    for( j=0; j<2*BLOCK_SIZE; j++)
    {
      pix_y = block_y+j;    
      ipix_y = img->pix_y + pix_y;
      for( i=0; i<2*BLOCK_SIZE; i++)
      {
        pix_x = block_x+i;
        if(lossless_qpprime)
          img->m7[j][i] = img->m7[j][i]+img->mpr[pix_y][block_x+i];
        else
          img->m7[j][i] = clip1a((img->m7[j][i]+((long)img->mpr[pix_y][pix_x] << DQ_BITS_8)+DQ_ROUND_8)>>DQ_BITS_8);
        enc_picture->imgY[ipix_y][img->pix_x + pix_x]=img->m7[j][i];
      }
    }
  }
  else if(!lossless_qpprime)
  {
    for( j=0; j<2*BLOCK_SIZE; j++)
      for( i=0; i<2*BLOCK_SIZE; i++)
        img->m7[j][i] =(img->m7[j][i]+DQ_ROUND_8)>>DQ_BITS_8;
  }
  
  //  Decoded block moved to frame memory
  
  return nonzero;
}

#ifdef USE_INTRA_MDDT

const int KLTRow8[9][8][8] = 
{ // 0 
{
{ -22,   57,   63,  -59,   53,  -41,   23,   13},
{ -30,   65,   41,    7,  -46,   59,  -55,  -29},
{ -38,   51,  -16,   62,  -37,  -17,   66,   47},
{ -45,   30,  -58,   31,   46,  -46,  -36,  -60},
{ -51,    2,  -57,  -44,   31,   55,  -18,   64},
{ -55,  -27,  -16,  -59,  -51,    1,   54,  -58},
{ -56,  -48,   32,    4,  -37,  -60,  -57,   41},
{ -51,  -48,   51,   51,   55,   45,   27,  -17},
},
{ // 1
{ -18,  -38,   50,  -57,   57,  -56,   45,  -24},
{ -30,  -58,   57,  -27,  -13,   45,  -63,   44},
{ -38,  -57,   14,   49,  -58,   17,   45,  -59},
{ -44,  -42,  -42,   55,   22,  -58,   -3,   64},
{ -50,  -15,  -65,  -20,   53,   35,  -40,  -59},
{ -54,   19,  -37,  -63,  -37,   27,   60,   45},
{ -57,   49,   17,  -15,  -50,  -66,  -52,  -28},
{ -56,   60,   53,   49,   49,   37,   21,   10},
},
{ // 2
{ -24,  -55,   65,  -51,   47,  -49,   31,  -20},
{ -34,  -62,   45,    1,  -33,   60,  -55,   41},
{ -41,  -50,  -19,   62,  -48,   -6,   46,  -60},
{ -47,  -30,  -60,   30,   45,  -47,  -16,   65},
{ -50,    0,  -53,  -50,   40,   51,  -28,  -59},
{ -53,   29,  -14,  -58,  -54,    7,   63,   47},
{ -54,   50,   26,    1,  -37,  -62,  -66,  -26},
{ -50,   53,   50,   55,   54,   40,   31,    8},
},
{ // 3
{ -31,   61,   64,  -55,   48,  -37,   25,   12},
{ -39,   63,   30,   13,  -49,   60,  -54,  -30},
{ -45,   44,  -29,   63,  -35,  -24,   61,   47},
{ -49,   17,  -62,   19,   56,  -39,  -37,  -59},
{ -51,  -12,  -50,  -52,   25,   60,  -10,   64},
{ -51,  -37,   -7,  -55,  -58,   -9,   51,  -58},
{ -49,  -50,   36,   10,  -27,  -58,  -64,   43},
{ -44,  -49,   52,   55,   51,   48,   33,  -18},
},
{ // 4
{ -11,  -25,   42,  -56,   61,  -61,   49,  -30},
{ -22,  -49,   61,  -46,    4,   40,  -61,   48},
{ -34,  -61,   34,   30,  -62,   25,   40,  -58},
{ -44,  -53,  -24,   59,   10,  -60,    1,   61},
{ -52,  -25,  -63,   -4,   57,   33,  -40,  -56},
{ -58,   15,  -45,  -62,  -30,   26,   59,   44},
{ -60,   50,   16,  -17,  -48,  -61,  -54,  -29},
{ -55,   59,   55,   51,   47,   37,   25,   12},
},
{ // 5
{ -12,   28,   44,  -60,   63,  -61,   42,   20},
{ -24,   53,   60,  -40,   -7,   51,  -62,  -36},
{ -36,   62,   27,   37,  -59,   11,   54,   51},
{ -46,   49,  -30,   56,   20,  -57,  -18,  -62},
{ -52,   19,  -64,  -11,   53,   41,  -27,   63},
{ -57,  -19,  -42,  -62,  -35,   17,   55,  -53},
{ -59,  -50,   20,  -14,  -46,  -58,  -55,   36},
{ -54,  -57,   56,   50,   48,   38,   27,  -15},
},
{ // 6
{ -10,  -27,   45,  -61,   54,  -63,   49,  -25},
{ -22,  -50,   60,  -43,    0,   42,  -65,   45},
{ -34,  -61,   31,   34,  -59,   28,   44,  -57},
{ -44,  -52,  -24,   57,   10,  -62,   -2,   63},
{ -52,  -25,  -61,   -7,   60,   26,  -39,  -59},
{ -57,   15,  -45,  -59,  -27,   32,   59,   47},
{ -60,   49,   13,  -18,  -56,  -58,  -50,  -29},
{ -55,   59,   57,   50,   50,   32,   21,   11},
},
{ // 7
{ -30,   60,   64,  -56,   50,  -38,   25,   12},
{ -38,   63,   32,   11,  -48,   60,  -55,  -30},
{ -44,   46,  -26,   63,  -35,  -22,   61,   48},
{ -48,   19,  -62,   22,   53,  -41,  -35,  -60},
{ -51,  -10,  -52,  -51,   27,   58,  -13,   64},
{ -52,  -35,   -9,  -56,  -57,   -7,   52,  -57},
{ -51,  -50,   36,    9,  -29,  -59,  -62,   42},
{ -45,  -49,   51,   54,   52,   48,   32,  -18},
},
{ // 8
{ -13,  -33,   47,  -61,   58,  -57,   45,  -25},
{ -25,  -55,   60,  -34,   -8,   47,  -61,   45},
{ -36,  -61,   22,   44,  -60,   13,   42,  -58},
{ -46,  -47,  -36,   54,   24,  -58,   -7,   63},
{ -52,  -17,  -64,  -18,   53,   40,  -34,  -58},
{ -57,   20,  -36,  -60,  -41,   21,   60,   46},
{ -58,   49,   20,  -14,  -44,  -62,  -58,  -29},
{ -54,   57,   55,   51,   49,   37,   26,   10},
},
};

const int KLTCol8[9][8][8] = 
{
{  // 0
{ -22,  -32,  -40,  -44,  -49,  -53,  -55,  -55},
{ -37,  -52,  -56,  -43,  -15,   19,   51,   64},
{  43,   54,   24,  -34,  -70,  -46,   14,   53},
{ -60,  -31,   39,   61,   -5,  -59,  -23,   50},
{  64,   -9,  -58,   10,   53,  -31,  -52,   48},
{ -58,   54,   16,  -62,   39,   22,  -56,   30},
{  35,  -61,   50,  -12,  -36,   62,  -55,   23},
{ -21,   45,  -59,   61,  -57,   49,  -32,   13},
},
{  // 1
{ -26,  -37,  -42,  -43,  -47,  -54,  -55,  -49},
{  66,   66,   37,   16,    2,  -19,  -49,  -57},
{ -53,  -38,   32,   79,   48,  -20,  -42,  -18},
{ -30,   28,   63,    9,  -67,  -61,   11,   48},
{  68,  -26,  -43,   29,   16,  -56,  -23,   67},
{ -43,   57,   -9,  -51,   53,   -1,  -58,   48},
{  29,  -57,   60,  -29,  -23,   55,  -59,   30},
{ -17,   34,  -51,   60,  -62,   55,  -40,   18},
},
{ // 2
{ -29,  -38,  -43,  -46,  -49,  -51,  -51,  -49},
{  52,   61,   50,   26,   -4,  -32,  -53,  -55},
{  55,   43,   -9,  -58,  -64,  -22,   33,   50},
{ -57,   -6,   59,   40,  -37,  -61,   -3,   53},
{  58,  -34,  -50,   40,   37,  -50,  -36,   51},
{ -42,   55,   -5,  -54,   52,    5,  -62,   45},
{  33,  -62,   57,  -24,  -25,   55,  -56,   27},
{  17,  -36,   52,  -59,   62,  -55,   39,  -16},
},
{ // 3
{ -13,  -27,  -37,  -44,  -52,  -57,  -58,  -53},
{ -27,  -52,  -62,  -50,  -17,   21,   51,   56},
{  43,   61,   27,  -35,  -67,  -37,   21,   50},
{ -56,  -39,   36,   56,  -15,  -62,  -12,   55},
{  67,   -3,  -59,   18,   50,  -37,  -44,   49},
{ -60,   46,   13,  -58,   40,   18,  -61,   39},
{  44,  -61,   47,  -10,  -34,   59,  -56,   26},
{ -24,   43,  -57,   62,  -59,   49,  -32,   12},
},
{ // 4
{ -11,  -22,  -34,  -44,  -53,  -59,  -60,  -53},
{  22,   47,   62,   54,   25,  -15,  -51,  -59},
{  39,   62,   38,  -23,  -64,  -43,   16,   54},
{ -58,  -47,   29,   59,   -5,  -58,  -16,   52},
{  68,    1,  -59,   14,   52,  -33,  -46,   47},
{ -62,   46,   14,  -57,   37,   22,  -60,   38},
{  45,  -63,   49,  -12,  -32,   58,  -54,   25},
{  21,  -40,   55,  -63,   61,  -50,   32,  -12},
},
{ // 5
{ -11,  -22,  -32,  -42,  -51,  -58,  -61,  -57},
{ -24,  -45,  -60,  -55,  -31,    9,   49,   62},
{  42,   62,   40,  -17,  -62,  -47,   11,   52},
{ -59,  -45,   28,   61,   -1,  -58,  -21,   50},
{  63,    2,  -59,    7,   57,  -28,  -51,   48},
{ -62,   43,   24,  -61,   31,   27,  -59,   34},
{  46,  -64,   45,   -5,  -38,   60,  -51,   23},
{ -24,   45,  -58,   61,  -57,   48,  -31,   12},
},
{ // 6
{ -13,  -26,  -37,  -45,  -52,  -56,  -59,  -54},
{  27,   53,   62,   48,   18,  -19,  -51,  -57},
{  40,   58,   29,  -30,  -65,  -43,   17,   57},
{ -57,  -45,   34,   60,   -9,  -60,  -15,   50},
{  65,   -3,  -62,   19,   53,  -36,  -43,   45},
{ -66,   51,    7,  -54,   43,   12,  -55,   37},
{  42,  -61,   54,  -25,  -19,   54,  -59,   29},
{  17,  -34,   50,  -60,   63,  -56,   39,  -15},
},
{ // 7
{ -14,  -27,  -37,  -44,  -52,  -56,  -58,  -55},
{ -29,  -52,  -61,  -48,  -20,   18,   51,   60},
{  43,   59,   29,  -30,  -67,  -42,   17,   53},
{ -55,  -39,   32,   60,   -6,  -62,  -20,   54},
{  66,   -1,  -61,   10,   54,  -30,  -48,   46},
{ -62,   47,   17,  -60,   37,   22,  -59,   35},
{  43,  -61,   47,  -10,  -35,   60,  -55,   25},
{ -23,   44,  -57,   61,  -58,   49,  -31,   12},
},
{ // 8
{ -29,  -37,  -42,  -46,  -49,  -51,  -52,  -50},
{  67,   68,   44,   13,  -15,  -33,  -44,  -45},
{  66,   21,  -43,  -67,  -43,    1,   37,   47},
{ -54,   23,   61,    1,  -62,  -46,   15,   54},
{  45,  -50,  -22,   59,   12,  -64,  -25,   55},
{ -34,   58,  -28,  -35,   61,  -15,  -59,   49},
{  24,  -53,   63,  -41,   -5,   49,  -63,   31},
{  10,  -27,   44,  -57,   64,  -62,   46,  -19},
},
};

static int KLTCol8inProd[9][8], KLTRow8inProd[9][8];

/*! 
*************************************************************************************
* \brief
*   Compute inner product for the input 8x8 matrix 
*
* \para compute_inner_product8x8()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
void compute_inner_product8x8(int matrix[8][8], int isLeft, int *inner_prod)
{
  int i, j;

  if(isLeft)
    for(i = 0; i < 8; i++)
    {
      inner_prod[i] = 0;
      for(j = 0; j < 4; j++)
      {
        inner_prod[i] += matrix[i][2*j+1]*matrix[i][2*j];
      }
    }
  else
    for(i = 0; i < 8; i++)
    {
      inner_prod[i] = 0;
      for(j = 0; j < 4; j++)
      {
        inner_prod[i] += matrix[2*j+1][i]*matrix[2*j][i];
      }
    }
}

/*! 
*************************************************************************************
* \brief
*   Precompute inner products for all 8x8 transform matrices
*
* \para precompute_all_inner_product8x8()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
void precompute_all_inner_product8x8()
{
  int i, j, k;
  int temp2D[8][8];


  for(i = 0; i < 9; i++)
  {
    if(i == 2) // handled by DCT 
      continue;

    for(j = 0; j < 8; j++)
      for(k = 0; k < 8; k++)
        temp2D[j][k] = KLTCol8[i][j][k];
    compute_inner_product8x8(temp2D, 1, KLTCol8inProd[i]);

    for(j = 0; j < 8; j++)
      for(k = 0; k < 8; k++)
        temp2D[j][k] = KLTRow8[i][j][k];
    compute_inner_product8x8(temp2D, 0, KLTRow8inProd[i]);
  }
}

/*! 
*************************************************************************************
* \brief
*   Perform fast mode-dependent 8x8 directional transform 
*
* \para klt_luma8x8_sep_fast()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/

int klt_luma8x8_sep_fast(int b8,int *coeff_cost, int ipmode)
{
  int sign(int a,int b);

  int i,j,k,coeff_ctr;
  int scan_pos,run;
  int nonzero;
  int qp_per,qp_rem, q_bits;

  int block_x = 8*(b8 & 0x01);
  int block_y = 8*(b8 >> 1);
  int*  ACLevel = img->cofAC[b8][0][0];
  int*  ACRun   = img->cofAC[b8][0][1];

  int coeff [2*BLOCK_SIZE][2*BLOCK_SIZE], ilev [2*BLOCK_SIZE][2*BLOCK_SIZE], 
      temp2D[2*BLOCK_SIZE][2*BLOCK_SIZE], y2D  [2*BLOCK_SIZE][2*BLOCK_SIZE]; 
  int inprod[BLOCK_SIZE*2];
  int level;

  int m=4, n=10, KLTprec=14;
  int qpOffset;
  int roundOffset;
  int R[6] = {40, 45, 50, 57, 63, 71};
  int Q[6] = {410,   364,   328,   287,   260,   231};

  long shift, shiftBy2;

  static int first=1;

  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int MCcoeff;
  int runs[4], scan_poss[4]; 

#ifdef RDO_Q
  levelDataStruct levelData[64];
  double lambda_md = 0.0;
  double err;
  int lowerInt;
  int levelRDOQ[64];
  int kStart = 0;
  int kStop = 0;
  int noCoeff;
  double normFact = pow(2, 2 * (KLTprec + m + n) - 15);
#endif


  qp_per    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)/6; //CHA-VG01
  qp_rem    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)%6; //CHA-VG01
  q_bits    = Q_BITS + qp_per;

  qpOffset=n+qp_per+KLTprec-2;
  roundOffset=(1<<(qpOffset))/3;

#ifdef RDO_Q
  if(input->UseRDO_Q)
  {
    if ((img->type == B_SLICE) && img->nal_reference_idc)
    {
      lambda_md = img->lambda_md[5][img->masterQP];  
    }
    else
    {
      lambda_md = img->lambda_md[img->type][img->masterQP];
    }
  }
#endif

  if (first==1)
  {
    first=0;
    img->delta_shift8x8=1<<(n+qp_per+KLTprec-11); 
    img->QPVal8x8=1<<(n+qp_per+KLTprec-3);  
    for (i=0; i<64; i++)
    {
      img->qp_shift8x8[i]=(1<<(n+qp_per+KLTprec-2))/3;
    }
  }

  for (j=0;j<(2*BLOCK_SIZE);j++)
    for(k = 0; k < 2*BLOCK_SIZE; k++)
      temp2D[j][k] = img->m7[j][k];

  compute_inner_product8x8(temp2D, 0, inprod);

  // column transform first
  for (i=0;i<(2*BLOCK_SIZE);i++)
  {
    for (j=0;j<(2*BLOCK_SIZE);j++)
    {
      temp2D[i][j] = 0;
      for(k = 0; k < BLOCK_SIZE; k++)
        temp2D[i][j]+=
        (KLTCol8[ipmode][i][2*k]+img->m7[2*k+1][j])*(KLTCol8[ipmode][i][2*k+1]+img->m7[2*k][j]);
      temp2D[i][j] = temp2D[i][j] - KLTCol8inProd[ipmode][i] - inprod[j];
    }
  }

  compute_inner_product8x8(temp2D, 1, inprod);

  // row transform second
  for (i=0;i<(2*BLOCK_SIZE);i++)
  {
    for (j=0;j<(2*BLOCK_SIZE);j++)
    {
      coeff[i][j] = 0;
      for(k = 0; k < BLOCK_SIZE; k++)
        coeff[i][j]+=
        (temp2D[i][2*k+1]+KLTRow8[ipmode][2*k][j])*(temp2D[i][2*k]+KLTRow8[ipmode][2*k+1][j]);
      coeff[i][j] = coeff[i][j] - inprod[i] - KLTRow8inProd[ipmode][j];
    }
  }

#ifdef RDO_Q
  if(input->UseRDO_Q)
  {
    noCoeff = 0;
    for (coeff_ctr = 0; coeff_ctr < 64; coeff_ctr++)
    {
      i=img->scanOrder8x8[ipmode][coeff_ctr][0];
      j=img->scanOrder8x8[ipmode][coeff_ctr][1];

      levelData[coeff_ctr].levelDouble = (int64)absm(coeff[j][i]) * Q[qp_rem];
      level = (int)(levelData[coeff_ctr].levelDouble >> qpOffset);
      lowerInt = ((levelData[coeff_ctr].levelDouble - (level << qpOffset)) < ((int64) 1 << (qpOffset - 1)))? 1: 0;
      
      levelData[coeff_ctr].level[0] = 0;

      if (level == 0 && lowerInt == 1)
      {
        levelData[coeff_ctr].noLevels = 1;
      }
      else if (level == 0 && lowerInt == 0)
      {
        levelData[coeff_ctr].level[1] = level + 1;
        levelData[coeff_ctr].noLevels = 2;
        kStop = coeff_ctr;
        noCoeff++;
      }
      else if (level > 0 && lowerInt == 1)
      {
        levelData[coeff_ctr].level[1] = level;
        levelData[coeff_ctr].noLevels = 2;
        kStop = coeff_ctr;
        noCoeff++;
      }
      else
      {
        levelData[coeff_ctr].level[1] = level;
        levelData[coeff_ctr].level[2] = level + 1;
        levelData[coeff_ctr].noLevels = 3;
        kStop = coeff_ctr;
        kStart = coeff_ctr;
        noCoeff++;
      }

      for (k = 0; k < levelData[coeff_ctr].noLevels; k++)
      {      
        err = (double)(levelData[coeff_ctr].level[k] << qpOffset) - (double) levelData[coeff_ctr].levelDouble;
        err = err * R[qp_rem];
        levelData[coeff_ctr].errLevel[k] = (err * err) / normFact; 
      }
    }

    est_writeRunLevel_CABAC(levelData, levelRDOQ, LUMA_8x8, lambda_md, kStart, kStop, noCoeff, 0);
  }
#endif

  nonzero=FALSE;

  run=-1;
  scan_pos=0;

  runs[0]=runs[1]=runs[2]=runs[3]=-1;
  scan_poss[0]=scan_poss[1]=scan_poss[2]=scan_poss[3]=0;

  for (coeff_ctr=0;coeff_ctr < 64;coeff_ctr++)
  {   
    i=img->scanOrder8x8[ipmode][coeff_ctr][0];
    j=img->scanOrder8x8[ipmode][coeff_ctr][1];

    run++;
    ilev[j][i]=0;

    MCcoeff = (coeff_ctr&3);
    runs[MCcoeff]++;

#ifdef RDO_Q
      if(input->UseRDO_Q)
      {
        level = levelRDOQ[coeff_ctr];
      }
      else
#endif
    if (img->AdaptiveRounding)
    {
      level = (int)((absm((int64)coeff[j][i]) * Q[qp_rem] + img->qp_shift8x8[coeff_ctr]) >> qpOffset);
      if (level != 0)
      {
        quant_stat_rd8x8[coeff_ctr] = (long)(absm((int64)coeff[j][i]) * Q[qp_rem] - ((int64)level << qpOffset));
      }
      else if (2 * absm((int64)coeff[j][i]) * Q[qp_rem] > ((int64)1 << qpOffset))
      {
        quant_stat_rd8x8[coeff_ctr] = (long)(2 * absm((int64)coeff[j][i]) * Q[qp_rem] - ((int64)1 << qpOffset));
      }
    }
    else
      {
        level = (int)((absm((int64)coeff[j][i]) * Q[qp_rem] + roundOffset) >> qpOffset);
      }

    if (level != 0)
    {
      nonzero=TRUE;
      if (currMB->luma_transform_size_8x8_flag && input->symbol_mode == UVLC)
      {
        *coeff_cost += (level > 1) ? MAX_VALUE : COEFF_COST8x8[input->disthres][runs[MCcoeff]];
        img->cofAC[b8][MCcoeff][0][scan_poss[MCcoeff]] = sign(level, coeff[j][i]);
        img->cofAC[b8][MCcoeff][1][scan_poss[MCcoeff]] = runs[MCcoeff];
        ++scan_poss[MCcoeff];
        runs[MCcoeff]=-1;
      }
      else
      {
        if (level > 1)
          *coeff_cost += MAX_VALUE;                // set high cost, shall not be discarded
        else
          *coeff_cost += COEFF_COST8x8[input->disthres][run];

        ACLevel[scan_pos] = sign(level, coeff[j][i]);
        ACRun  [scan_pos] = run;
        ++scan_pos;
        run=-1;                     // reset zero level counter
      }

      ilev[j][i] = sign(level, coeff[j][i]) * (R[qp_rem] << qp_per); 
    }
    else
      ilev[j][i] = 0;
  }
  if (currMB->luma_transform_size_8x8_flag && input->symbol_mode == UVLC)
  {
    for(i=0; i<4; i++)
      img->cofAC[b8][i][0][scan_poss[i]] = 0;
  }
  else
    ACLevel[scan_pos] = 0;

  // inverse transform, column first
  for (i=0;i<(2*BLOCK_SIZE);i++)
  {
    for (j=0;j<(2*BLOCK_SIZE);j++)
    {
      temp2D[i][j] = 0;
      for(k = 0; k < 2*BLOCK_SIZE; k++)
        temp2D[i][j]+=KLTCol8[ipmode][k][i]*(ilev[k][j]);
    }
  }

  // inverse transform, row second 
  for (i=0;i<(2*BLOCK_SIZE);i++)
  {
    for (j=0;j<(2*BLOCK_SIZE);j++)
    {
      y2D[i][j] = 0;
      for(k = 0; k < 2*BLOCK_SIZE; k++)
        y2D[i][j]+=temp2D[i][k]*KLTRow8[ipmode][j][k];
    }
  }

  shift=KLTprec+m+2; shiftBy2=1<<(shift-1);
  for (j=0;j<2*BLOCK_SIZE;j++)
  {
    for (i=0;i<2*BLOCK_SIZE;i++)
    {
      int temp;

      temp=(y2D[j][i]+(img->mpr[j+block_y][i+block_x]<<shift)+shiftBy2)>>shift;
      temp=min(img->max_imgpel_value,max(0,temp));
      enc_picture->imgY[img->pix_y+block_y+j][img->pix_x+block_x+i]= (imgpel) temp;
    }
  }

  return nonzero;
}

#endif
