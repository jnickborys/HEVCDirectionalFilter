
/*!
 ***************************************************************************
 * \file transform8x8.c
 *
 * \brief
 *    8x8 transform functions
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details) 
 *    - Yuri Vatis          <vatis@hhi.de>
 *    - Jan Muenster        <muenster@hhi.de>
 *
 * \date
 *    12. October 2003
 **************************************************************************
 */

#include <stdlib.h>
#include "transform8x8.h"



#define Q_BITS_8        16
#define DQ_BITS_8       6 
#define DQ_ROUND_8      (1<<(DQ_BITS_8-1))

#ifdef _NEW_8x8_ARRAYS_INCLUDED_
//! single scan pattern
const byte SNGL_SCAN8x8[64][2] = {
  {0,0}, {1,0}, {0,1}, {0,2}, {1,1}, {2,0}, {3,0}, {2,1}, {1,2}, {0,3}, {0,4}, {1,3}, {2,2}, {3,1}, {4,0}, {5,0},
  {4,1}, {3,2}, {2,3}, {1,4}, {0,5}, {0,6}, {1,5}, {2,4}, {3,3}, {4,2}, {5,1}, {6,0}, {7,0}, {6,1}, {5,2}, {4,3},
  {3,4}, {2,5}, {1,6}, {0,7}, {1,7}, {2,6}, {3,5}, {4,4}, {5,3}, {6,2}, {7,1}, {7,2}, {6,3}, {5,4}, {4,5}, {3,6},
  {2,7}, {3,7}, {4,6}, {5,5}, {6,4}, {7,3}, {7,4}, {6,5}, {5,6}, {4,7}, {5,7}, {6,6}, {7,5}, {7,6}, {6,7}, {7,7}
};

//! field scan pattern
const byte FIELD_SCAN8x8[64][2] = {
 {0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {0,3}, {0,4}, {1,2}, {2,0}, {2,1}, {1,3}, {0,5}, {0,6}, {1,4}, {2,2}, {3,0},
 {3,1}, {2,3}, {1,5}, {0,7}, {1,6}, {2,4}, {3,2}, {4,0}, {4,1}, {3,3}, {2,5}, {1,7}, {2,6}, {3,4}, {4,2}, {5,0},
 {5,1}, {4,3}, {3,5}, {2,7}, {3,6}, {4,4}, {5,2}, {6,0}, {6,1}, {5,3}, {4,5}, {3,7}, {4,6}, {5,4}, {6,2}, {7,0},
 {7,1}, {6,3}, {5,5}, {4,7}, {5,6}, {6,4}, {7,2}, {7,3}, {6,5}, {5,7}, {6,6}, {7,4}, {7,5}, {6,7}, {7,6}, {7,7}
};


//! array used to find expencive coefficients
const byte COEFF_COST8x8[64] =
{
  3,3,3,3,2,2,2,2,2,2,2,2,1,1,1,1,
  1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

#endif

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
 ************************************************************************
 */
int intrapred8x8( struct img_par *img,  //!< image parameters
                  int b8)

{
  int i,j;
  int s0;
  int PredPel[25];  // array of predictor pels
  imgpel **imgY = dec_picture->imgY;  // For MB level frame/field coding tools -- set default to imgY

  int mb_nr=img->current_mb_nr;

  PixelPos pix_a[8];
  PixelPos pix_b, pix_c, pix_d;

  int block_available_up;
  int block_available_left;
  int block_available_up_left;
  int block_available_up_right;
  int img_block_x = (img->mb_x)*4 + 2*(b8%2);
  int img_block_y = (img->mb_y)*4 + 2*(b8/2);
  int ioff = (b8%2)*8;
  int joff = (b8/2)*8;

  byte predmode = img->ipredmode[img_block_x][img_block_y];

  for (i=0;i<8;i++)
  {
    getNeighbour(mb_nr, ioff -1 , joff +i , 1, &pix_a[i]);
  }

  getNeighbour(mb_nr, ioff    , joff -1 , 1, &pix_b);
  getNeighbour(mb_nr, ioff +8 , joff -1 , 1, &pix_c);
  getNeighbour(mb_nr, ioff -1 , joff -1 , 1, &pix_d);
  
  pix_c.available = pix_c.available &&!(ioff == 8 && joff == 8);

  if (active_pps->constrained_intra_pred_flag)
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

//  *left_available = block_available_left;
//  *up_available   = block_available_up;
//  *all_available  = block_available_up && block_available_left && block_available_up_left;

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
  
  LowPassForIntra8x8Pred(&(P_Z), block_available_up_left, block_available_up, block_available_left);

//img->mpr[x][y]
  switch(predmode)
  {
  case DC_PRED:
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
    for(i = 0; i < 2*BLOCK_SIZE; i++)
      for(j = 0; j < 2*BLOCK_SIZE; j++)
        img->mpr[i+ioff][j+joff] = s0;
    break;

  case VERT_PRED:
    if (!block_available_up)
      printf ("warning: Intra_8x8_Vertical prediction mode not allowed at mb %d\n",img->current_mb_nr);

    for (i=0; i < 2*BLOCK_SIZE; i++)
    {
      img->mpr[i+ioff][0+joff] = 
      img->mpr[i+ioff][1+joff] = 
      img->mpr[i+ioff][2+joff] = 
      img->mpr[i+ioff][3+joff] = 
      img->mpr[i+ioff][4+joff] = 
      img->mpr[i+ioff][5+joff] = 
      img->mpr[i+ioff][6+joff] = 
      img->mpr[i+ioff][7+joff] = (&P_A)[i];
    }
    break;
  case HOR_PRED:
    if (!block_available_left)
      printf ("warning: Intra_8x8_Horizontal prediction mode not allowed at mb %d\n",img->current_mb_nr);

    for (j=0; j < 2*BLOCK_SIZE; j++)
    {
      img->mpr[0+ioff][j+joff]  = 
      img->mpr[1+ioff][j+joff]  = 
      img->mpr[2+ioff][j+joff]  = 
      img->mpr[3+ioff][j+joff]  = 
      img->mpr[4+ioff][j+joff]  = 
      img->mpr[5+ioff][j+joff]  = 
      img->mpr[6+ioff][j+joff]  = 
      img->mpr[7+ioff][j+joff]  = (&P_Q)[j];
    }
    break;

  case DIAG_DOWN_LEFT_PRED:
    if (!block_available_up)
      printf ("warning: Intra_8x8_Diagonal_Down_Left prediction mode not allowed at mb %d\n",img->current_mb_nr);
    // Mode DIAG_DOWN_LEFT_PRED
    img->mpr[0+ioff][0+joff] = (P_A + P_C + 2*(P_B) + 2) >> 2;
    img->mpr[0+ioff][1+joff] = 
    img->mpr[1+ioff][0+joff] = (P_B + P_D + 2*(P_C) + 2) >> 2;
    img->mpr[0+ioff][2+joff] =
    img->mpr[1+ioff][1+joff] =
    img->mpr[2+ioff][0+joff] = (P_C + P_E + 2*(P_D) + 2) >> 2;
    img->mpr[0+ioff][3+joff] = 
    img->mpr[1+ioff][2+joff] = 
    img->mpr[2+ioff][1+joff] = 
    img->mpr[3+ioff][0+joff] = (P_D + P_F + 2*(P_E) + 2) >> 2;
    img->mpr[0+ioff][4+joff] = 
    img->mpr[1+ioff][3+joff] = 
    img->mpr[2+ioff][2+joff] = 
    img->mpr[3+ioff][1+joff] = 
    img->mpr[4+ioff][0+joff] = (P_E + P_G + 2*(P_F) + 2) >> 2;
    img->mpr[0+ioff][5+joff] = 
    img->mpr[1+ioff][4+joff] = 
    img->mpr[2+ioff][3+joff] = 
    img->mpr[3+ioff][2+joff] = 
    img->mpr[4+ioff][1+joff] = 
    img->mpr[5+ioff][0+joff] = (P_F + P_H + 2*(P_G) + 2) >> 2;
    img->mpr[0+ioff][6+joff] = 
    img->mpr[1+ioff][5+joff] = 
    img->mpr[2+ioff][4+joff] = 
    img->mpr[3+ioff][3+joff] = 
    img->mpr[4+ioff][2+joff] = 
    img->mpr[5+ioff][1+joff] = 
    img->mpr[6+ioff][0+joff] = (P_G + P_I + 2*(P_H) + 2) >> 2;
    img->mpr[0+ioff][7+joff] = 
    img->mpr[1+ioff][6+joff] = 
    img->mpr[2+ioff][5+joff] = 
    img->mpr[3+ioff][4+joff] = 
    img->mpr[4+ioff][3+joff] = 
    img->mpr[5+ioff][2+joff] = 
    img->mpr[6+ioff][1+joff] = 
    img->mpr[7+ioff][0+joff] = (P_H + P_J + 2*(P_I) + 2) >> 2;
    img->mpr[1+ioff][7+joff] = 
    img->mpr[2+ioff][6+joff] = 
    img->mpr[3+ioff][5+joff] = 
    img->mpr[4+ioff][4+joff] = 
    img->mpr[5+ioff][3+joff] = 
    img->mpr[6+ioff][2+joff] = 
    img->mpr[7+ioff][1+joff] = (P_I + P_K + 2*(P_J) + 2) >> 2;
    img->mpr[2+ioff][7+joff] = 
    img->mpr[3+ioff][6+joff] = 
    img->mpr[4+ioff][5+joff] = 
    img->mpr[5+ioff][4+joff] = 
    img->mpr[6+ioff][3+joff] = 
    img->mpr[7+ioff][2+joff] = (P_J + P_L + 2*(P_K) + 2) >> 2;
    img->mpr[3+ioff][7+joff] = 
    img->mpr[4+ioff][6+joff] = 
    img->mpr[5+ioff][5+joff] = 
    img->mpr[6+ioff][4+joff] = 
    img->mpr[7+ioff][3+joff] = (P_K + P_M + 2*(P_L) + 2) >> 2;
    img->mpr[4+ioff][7+joff] = 
    img->mpr[5+ioff][6+joff] = 
    img->mpr[6+ioff][5+joff] = 
    img->mpr[7+ioff][4+joff] = (P_L + P_N + 2*(P_M) + 2) >> 2;
    img->mpr[5+ioff][7+joff] = 
    img->mpr[6+ioff][6+joff] = 
    img->mpr[7+ioff][5+joff] = (P_M + P_O + 2*(P_N) + 2) >> 2;
    img->mpr[6+ioff][7+joff] = 
    img->mpr[7+ioff][6+joff] = (P_N + P_P + 2*(P_O) + 2) >> 2;
    img->mpr[7+ioff][7+joff] = (P_O + 3*(P_P) + 2) >> 2;
    break;

  case VERT_LEFT_PRED:
    if (!block_available_up)
      printf ("warning: Intra_4x4_Vertical_Left prediction mode not allowed at mb %d\n",img->current_mb_nr);

    img->mpr[0+ioff][0+joff] = (P_A + P_B + 1) >> 1;
    img->mpr[1+ioff][0+joff] = 
    img->mpr[0+ioff][2+joff] = (P_B + P_C + 1) >> 1;
    img->mpr[2+ioff][0+joff] = 
    img->mpr[1+ioff][2+joff] = 
    img->mpr[0+ioff][4+joff] = (P_C + P_D + 1) >> 1;
    img->mpr[3+ioff][0+joff] = 
    img->mpr[2+ioff][2+joff] = 
    img->mpr[1+ioff][4+joff] = 
    img->mpr[0+ioff][6+joff] = (P_D + P_E + 1) >> 1;
    img->mpr[4+ioff][0+joff] = 
    img->mpr[3+ioff][2+joff] = 
    img->mpr[2+ioff][4+joff] = 
    img->mpr[1+ioff][6+joff] = (P_E + P_F + 1) >> 1;
    img->mpr[5+ioff][0+joff] = 
    img->mpr[4+ioff][2+joff] = 
    img->mpr[3+ioff][4+joff] = 
    img->mpr[2+ioff][6+joff] = (P_F + P_G + 1) >> 1;
    img->mpr[6+ioff][0+joff] = 
    img->mpr[5+ioff][2+joff] = 
    img->mpr[4+ioff][4+joff] = 
    img->mpr[3+ioff][6+joff] = (P_G + P_H + 1) >> 1;
    img->mpr[7+ioff][0+joff] = 
    img->mpr[6+ioff][2+joff] = 
    img->mpr[5+ioff][4+joff] = 
    img->mpr[4+ioff][6+joff] = (P_H + P_I + 1) >> 1;
    img->mpr[7+ioff][2+joff] = 
    img->mpr[6+ioff][4+joff] = 
    img->mpr[5+ioff][6+joff] = (P_I + P_J + 1) >> 1;
    img->mpr[7+ioff][4+joff] = 
    img->mpr[6+ioff][6+joff] = (P_J + P_K + 1) >> 1;
    img->mpr[7+ioff][6+joff] = (P_K + P_L + 1) >> 1;
    img->mpr[0+ioff][1+joff] = (P_A + P_C + 2*P_B + 2) >> 2;
    img->mpr[1+ioff][1+joff] = 
    img->mpr[0+ioff][3+joff] = (P_B + P_D + 2*P_C + 2) >> 2;
    img->mpr[2+ioff][1+joff] = 
    img->mpr[1+ioff][3+joff] = 
    img->mpr[0+ioff][5+joff] = (P_C + P_E + 2*P_D + 2) >> 2;
    img->mpr[3+ioff][1+joff] = 
    img->mpr[2+ioff][3+joff] = 
    img->mpr[1+ioff][5+joff] = 
    img->mpr[0+ioff][7+joff] = (P_D + P_F + 2*P_E + 2) >> 2;
    img->mpr[4+ioff][1+joff] = 
    img->mpr[3+ioff][3+joff] = 
    img->mpr[2+ioff][5+joff] = 
    img->mpr[1+ioff][7+joff] = (P_E + P_G + 2*P_F + 2) >> 2;
    img->mpr[5+ioff][1+joff] = 
    img->mpr[4+ioff][3+joff] = 
    img->mpr[3+ioff][5+joff] = 
    img->mpr[2+ioff][7+joff] = (P_F + P_H + 2*P_G + 2) >> 2;
    img->mpr[6+ioff][1+joff] = 
    img->mpr[5+ioff][3+joff] = 
    img->mpr[4+ioff][5+joff] = 
    img->mpr[3+ioff][7+joff] = (P_G + P_I + 2*P_H + 2) >> 2;
    img->mpr[7+ioff][1+joff] = 
    img->mpr[6+ioff][3+joff] = 
    img->mpr[5+ioff][5+joff] = 
    img->mpr[4+ioff][7+joff] = (P_H + P_J + 2*P_I + 2) >> 2;
    img->mpr[7+ioff][3+joff] = 
    img->mpr[6+ioff][5+joff] = 
    img->mpr[5+ioff][7+joff] = (P_I + P_K + 2*P_J + 2) >> 2;
    img->mpr[7+ioff][5+joff] = 
    img->mpr[6+ioff][7+joff] = (P_J + P_L + 2*P_K + 2) >> 2;
    img->mpr[7+ioff][7+joff] = (P_K + P_M + 2*P_L + 2) >> 2;
    break;


  case DIAG_DOWN_RIGHT_PRED:
    if ((!block_available_up)||(!block_available_left)||(!block_available_up_left))
      printf ("warning: Intra_8x8_Diagonal_Down_Right prediction mode not allowed at mb %d\n",img->current_mb_nr);

    // Mode DIAG_DOWN_RIGHT_PRED
    img->mpr[0+ioff][7+joff] = (P_X + P_V + 2*(P_W) + 2) >> 2;
    img->mpr[0+ioff][6+joff] = 
    img->mpr[1+ioff][7+joff] = (P_W + P_U + 2*(P_V) + 2) >> 2;
    img->mpr[0+ioff][5+joff] = 
    img->mpr[1+ioff][6+joff] = 
    img->mpr[2+ioff][7+joff] = (P_V + P_T + 2*(P_U) + 2) >> 2;
    img->mpr[0+ioff][4+joff] = 
    img->mpr[1+ioff][5+joff] = 
    img->mpr[2+ioff][6+joff] = 
    img->mpr[3+ioff][7+joff] = (P_U + P_S + 2*(P_T) + 2) >> 2;
    img->mpr[0+ioff][3+joff] = 
    img->mpr[1+ioff][4+joff] = 
    img->mpr[2+ioff][5+joff] = 
    img->mpr[3+ioff][6+joff] = 
    img->mpr[4+ioff][7+joff] = (P_T + P_R + 2*(P_S) + 2) >> 2;
    img->mpr[0+ioff][2+joff] = 
    img->mpr[1+ioff][3+joff] = 
    img->mpr[2+ioff][4+joff] = 
    img->mpr[3+ioff][5+joff] = 
    img->mpr[4+ioff][6+joff] = 
    img->mpr[5+ioff][7+joff] = (P_S + P_Q + 2*(P_R) + 2) >> 2;
    img->mpr[0+ioff][1+joff] = 
    img->mpr[1+ioff][2+joff] = 
    img->mpr[2+ioff][3+joff] = 
    img->mpr[3+ioff][4+joff] = 
    img->mpr[4+ioff][5+joff] = 
    img->mpr[5+ioff][6+joff] = 
    img->mpr[6+ioff][7+joff] = (P_R + P_Z + 2*(P_Q) + 2) >> 2;
    img->mpr[0+ioff][0+joff] = 
    img->mpr[1+ioff][1+joff] = 
    img->mpr[2+ioff][2+joff] = 
    img->mpr[3+ioff][3+joff] = 
    img->mpr[4+ioff][4+joff] = 
    img->mpr[5+ioff][5+joff] = 
    img->mpr[6+ioff][6+joff] = 
    img->mpr[7+ioff][7+joff] = (P_Q + P_A + 2*(P_Z) + 2) >> 2;
    img->mpr[1+ioff][0+joff] = 
    img->mpr[2+ioff][1+joff] = 
    img->mpr[3+ioff][2+joff] = 
    img->mpr[4+ioff][3+joff] = 
    img->mpr[5+ioff][4+joff] = 
    img->mpr[6+ioff][5+joff] = 
    img->mpr[7+ioff][6+joff] = (P_Z + P_B + 2*(P_A) + 2) >> 2;
    img->mpr[2+ioff][0+joff] = 
    img->mpr[3+ioff][1+joff] = 
    img->mpr[4+ioff][2+joff] = 
    img->mpr[5+ioff][3+joff] = 
    img->mpr[6+ioff][4+joff] = 
    img->mpr[7+ioff][5+joff] = (P_A + P_C + 2*(P_B) + 2) >> 2;
    img->mpr[3+ioff][0+joff] = 
    img->mpr[4+ioff][1+joff] = 
    img->mpr[5+ioff][2+joff] = 
    img->mpr[6+ioff][3+joff] = 
    img->mpr[7+ioff][4+joff] = (P_B + P_D + 2*(P_C) + 2) >> 2;
    img->mpr[4+ioff][0+joff] = 
    img->mpr[5+ioff][1+joff] = 
    img->mpr[6+ioff][2+joff] = 
    img->mpr[7+ioff][3+joff] = (P_C + P_E + 2*(P_D) + 2) >> 2;
    img->mpr[5+ioff][0+joff] = 
    img->mpr[6+ioff][1+joff] = 
    img->mpr[7+ioff][2+joff] = (P_D + P_F + 2*(P_E) + 2) >> 2;
    img->mpr[6+ioff][0+joff] = 
    img->mpr[7+ioff][1+joff] = (P_E + P_G + 2*(P_F) + 2) >> 2;
    img->mpr[7+ioff][0+joff] = (P_F + P_H + 2*(P_G) + 2) >> 2;
    break;

  case  VERT_RIGHT_PRED:/* diagonal prediction -22.5 deg to horizontal plane */
    if ((!block_available_up)||(!block_available_left)||(!block_available_up_left))
      printf ("warning: Intra_8x8_Vertical_Right prediction mode not allowed at mb %d\n",img->current_mb_nr);

    img->mpr[0+ioff][0+joff] = 
    img->mpr[1+ioff][2+joff] = 
    img->mpr[2+ioff][4+joff] = 
    img->mpr[3+ioff][6+joff] = (P_Z + P_A + 1) >> 1;
    img->mpr[1+ioff][0+joff] = 
    img->mpr[2+ioff][2+joff] = 
    img->mpr[3+ioff][4+joff] = 
    img->mpr[4+ioff][6+joff] = (P_A + P_B + 1) >> 1;
    img->mpr[2+ioff][0+joff] = 
    img->mpr[3+ioff][2+joff] = 
    img->mpr[4+ioff][4+joff] = 
    img->mpr[5+ioff][6+joff] = (P_B + P_C + 1) >> 1;
    img->mpr[3+ioff][0+joff] = 
    img->mpr[4+ioff][2+joff] = 
    img->mpr[5+ioff][4+joff] = 
    img->mpr[6+ioff][6+joff] = (P_C + P_D + 1) >> 1;
    img->mpr[4+ioff][0+joff] = 
    img->mpr[5+ioff][2+joff] = 
    img->mpr[6+ioff][4+joff] = 
    img->mpr[7+ioff][6+joff] = (P_D + P_E + 1) >> 1;
    img->mpr[5+ioff][0+joff] = 
    img->mpr[6+ioff][2+joff] = 
    img->mpr[7+ioff][4+joff] = (P_E + P_F + 1) >> 1;
    img->mpr[6+ioff][0+joff] = 
    img->mpr[7+ioff][2+joff] = (P_F + P_G + 1) >> 1;
    img->mpr[7+ioff][0+joff] = (P_G + P_H + 1) >> 1;
    img->mpr[0+ioff][1+joff] = 
    img->mpr[1+ioff][3+joff] = 
    img->mpr[2+ioff][5+joff] = 
    img->mpr[3+ioff][7+joff] = (P_Q + P_A + 2*P_Z + 2) >> 2;
    img->mpr[1+ioff][1+joff] = 
    img->mpr[2+ioff][3+joff] = 
    img->mpr[3+ioff][5+joff] = 
    img->mpr[4+ioff][7+joff] = (P_Z + P_B + 2*P_A + 2) >> 2;
    img->mpr[2+ioff][1+joff] = 
    img->mpr[3+ioff][3+joff] = 
    img->mpr[4+ioff][5+joff] = 
    img->mpr[5+ioff][7+joff] = (P_A + P_C + 2*P_B + 2) >> 2;
    img->mpr[3+ioff][1+joff] = 
    img->mpr[4+ioff][3+joff] = 
    img->mpr[5+ioff][5+joff] = 
    img->mpr[6+ioff][7+joff] = (P_B + P_D + 2*P_C + 2) >> 2;
    img->mpr[4+ioff][1+joff] = 
    img->mpr[5+ioff][3+joff] = 
    img->mpr[6+ioff][5+joff] = 
    img->mpr[7+ioff][7+joff] = (P_C + P_E + 2*P_D + 2) >> 2;
    img->mpr[5+ioff][1+joff] = 
    img->mpr[6+ioff][3+joff] = 
    img->mpr[7+ioff][5+joff] = (P_D + P_F + 2*P_E + 2) >> 2;
    img->mpr[6+ioff][1+joff] = 
    img->mpr[7+ioff][3+joff] = (P_E + P_G + 2*P_F + 2) >> 2;
    img->mpr[7+ioff][1+joff] = (P_F + P_H + 2*P_G + 2) >> 2;
    img->mpr[0+ioff][2+joff] =
    img->mpr[1+ioff][4+joff] =
    img->mpr[2+ioff][6+joff] = (P_R + P_Z + 2*P_Q + 2) >> 2;
    img->mpr[0+ioff][3+joff] =
    img->mpr[1+ioff][5+joff] =
    img->mpr[2+ioff][7+joff] = (P_S + P_Q + 2*P_R + 2) >> 2;
    img->mpr[0+ioff][4+joff] =
    img->mpr[1+ioff][6+joff] = (P_T + P_R + 2*P_S + 2) >> 2;
    img->mpr[0+ioff][5+joff] =
    img->mpr[1+ioff][7+joff] = (P_U + P_S + 2*P_T + 2) >> 2;
    img->mpr[0+ioff][6+joff] = (P_V + P_T + 2*P_U + 2) >> 2;
    img->mpr[0+ioff][7+joff] = (P_W + P_U + 2*P_V + 2) >> 2;
    break;

  case  HOR_DOWN_PRED:/* diagonal prediction -22.5 deg to horizontal plane */
    if ((!block_available_up)||(!block_available_left)||(!block_available_up_left))
      printf ("warning: Intra_8x8_Horizontal_Down prediction mode not allowed at mb %d\n",img->current_mb_nr);
    
    img->mpr[0+ioff][0+joff] = 
    img->mpr[2+ioff][1+joff] = 
    img->mpr[4+ioff][2+joff] = 
    img->mpr[6+ioff][3+joff] = (P_Q + P_Z + 1) >> 1;
    img->mpr[0+ioff][1+joff] = 
    img->mpr[2+ioff][2+joff] = 
    img->mpr[4+ioff][3+joff] = 
    img->mpr[6+ioff][4+joff] = (P_R + P_Q + 1) >> 1;
    img->mpr[0+ioff][2+joff] = 
    img->mpr[2+ioff][3+joff] = 
    img->mpr[4+ioff][4+joff] = 
    img->mpr[6+ioff][5+joff] = (P_S + P_R + 1) >> 1;
    img->mpr[0+ioff][3+joff] = 
    img->mpr[2+ioff][4+joff] = 
    img->mpr[4+ioff][5+joff] = 
    img->mpr[6+ioff][6+joff] = (P_T + P_S + 1) >> 1;
    img->mpr[0+ioff][4+joff] = 
    img->mpr[2+ioff][5+joff] = 
    img->mpr[4+ioff][6+joff] = 
    img->mpr[6+ioff][7+joff] = (P_U + P_T + 1) >> 1;
    img->mpr[0+ioff][5+joff] = 
    img->mpr[2+ioff][6+joff] = 
    img->mpr[4+ioff][7+joff] = (P_V + P_U + 1) >> 1;
    img->mpr[0+ioff][6+joff] = 
    img->mpr[2+ioff][7+joff] = (P_W + P_V + 1) >> 1;
    img->mpr[0+ioff][7+joff] = (P_X + P_W + 1) >> 1;
    img->mpr[1+ioff][0+joff] =
    img->mpr[3+ioff][1+joff] =
    img->mpr[5+ioff][2+joff] =
    img->mpr[7+ioff][3+joff] = (P_Q + P_A + 2*P_Z + 2) >> 2;
    img->mpr[1+ioff][1+joff] =
    img->mpr[3+ioff][2+joff] =
    img->mpr[5+ioff][3+joff] =
    img->mpr[7+ioff][4+joff] = (P_Z + P_R + 2*P_Q + 2) >> 2;
    img->mpr[1+ioff][2+joff] =
    img->mpr[3+ioff][3+joff] =
    img->mpr[5+ioff][4+joff] =
    img->mpr[7+ioff][5+joff] = (P_Q + P_S + 2*P_R + 2) >> 2;
    img->mpr[1+ioff][3+joff] =
    img->mpr[3+ioff][4+joff] =
    img->mpr[5+ioff][5+joff] =
    img->mpr[7+ioff][6+joff] = (P_R + P_T + 2*P_S + 2) >> 2;
    img->mpr[1+ioff][4+joff] =
    img->mpr[3+ioff][5+joff] =
    img->mpr[5+ioff][6+joff] =
    img->mpr[7+ioff][7+joff] = (P_S + P_U + 2*P_T + 2) >> 2;
    img->mpr[1+ioff][5+joff] =
    img->mpr[3+ioff][6+joff] =
    img->mpr[5+ioff][7+joff] = (P_T + P_V + 2*P_U + 2) >> 2;
    img->mpr[1+ioff][6+joff] =
    img->mpr[3+ioff][7+joff] = (P_U + P_W + 2*P_V + 2) >> 2;
    img->mpr[1+ioff][7+joff] = (P_V + P_X + 2*P_W + 2) >> 2;
    img->mpr[2+ioff][0+joff] = 
    img->mpr[4+ioff][1+joff] = 
    img->mpr[6+ioff][2+joff] = (P_Z + P_B + 2*P_A + 2) >> 2;
    img->mpr[3+ioff][0+joff] = 
    img->mpr[5+ioff][1+joff] = 
    img->mpr[7+ioff][2+joff] = (P_A + P_C + 2*P_B + 2) >> 2;
    img->mpr[4+ioff][0+joff] = 
    img->mpr[6+ioff][1+joff] = (P_B + P_D + 2*P_C + 2) >> 2;
    img->mpr[5+ioff][0+joff] = 
    img->mpr[7+ioff][1+joff] = (P_C + P_E + 2*P_D + 2) >> 2;
    img->mpr[6+ioff][0+joff] = (P_D + P_F + 2*P_E + 2) >> 2;
    img->mpr[7+ioff][0+joff] = (P_E + P_G + 2*P_F + 2) >> 2;
    break;

  case  HOR_UP_PRED:/* diagonal prediction -22.5 deg to horizontal plane */
    if (!block_available_left)
      printf ("warning: Intra_8x8_Horizontal_Up prediction mode not allowed at mb %d\n",img->current_mb_nr);

    img->mpr[0+ioff][0+joff] = (P_Q + P_R + 1) >> 1;
    img->mpr[0+ioff][1+joff] =
    img->mpr[2+ioff][0+joff] = (P_R + P_S + 1) >> 1;
    img->mpr[0+ioff][2+joff] =
    img->mpr[2+ioff][1+joff] =
    img->mpr[4+ioff][0+joff] = (P_S + P_T + 1) >> 1;
    img->mpr[0+ioff][3+joff] =
    img->mpr[2+ioff][2+joff] =
    img->mpr[4+ioff][1+joff] =
    img->mpr[6+ioff][0+joff] = (P_T + P_U + 1) >> 1;
    img->mpr[0+ioff][4+joff] =
    img->mpr[2+ioff][3+joff] =
    img->mpr[4+ioff][2+joff] =
    img->mpr[6+ioff][1+joff] = (P_U + P_V + 1) >> 1;
    img->mpr[0+ioff][5+joff] =
    img->mpr[2+ioff][4+joff] =
    img->mpr[4+ioff][3+joff] =
    img->mpr[6+ioff][2+joff] = (P_V + P_W + 1) >> 1;
    img->mpr[0+ioff][6+joff] =
    img->mpr[2+ioff][5+joff] =
    img->mpr[4+ioff][4+joff] =
    img->mpr[6+ioff][3+joff] = (P_W + P_X + 1) >> 1;
    img->mpr[6+ioff][4+joff] =
      img->mpr[7+ioff][4+joff] =
      img->mpr[4+ioff][5+joff] =
      img->mpr[5+ioff][5+joff] =
      img->mpr[6+ioff][5+joff] =
      img->mpr[7+ioff][5+joff] =
      img->mpr[2+ioff][6+joff] =
      img->mpr[3+ioff][6+joff] =
      img->mpr[4+ioff][6+joff] =
      img->mpr[5+ioff][6+joff] =
      img->mpr[6+ioff][6+joff] =
      img->mpr[7+ioff][6+joff] =
      img->mpr[0+ioff][7+joff] =
      img->mpr[1+ioff][7+joff] =
      img->mpr[2+ioff][7+joff] =
      img->mpr[3+ioff][7+joff] =
      img->mpr[4+ioff][7+joff] =
      img->mpr[5+ioff][7+joff] =
      img->mpr[6+ioff][7+joff] =
      img->mpr[7+ioff][7+joff] = P_X;
    img->mpr[1+ioff][6+joff] =
      img->mpr[3+ioff][5+joff] =
      img->mpr[5+ioff][4+joff] =
      img->mpr[7+ioff][3+joff] = (P_W + 3*P_X + 2) >> 2;
    img->mpr[1+ioff][5+joff] =
      img->mpr[3+ioff][4+joff] =
      img->mpr[5+ioff][3+joff] =
      img->mpr[7+ioff][2+joff] = (P_X + P_V + 2*P_W + 2) >> 2;
    img->mpr[1+ioff][4+joff] =
      img->mpr[3+ioff][3+joff] =
      img->mpr[5+ioff][2+joff] =
      img->mpr[7+ioff][1+joff] = (P_W + P_U + 2*P_V + 2) >> 2;
    img->mpr[1+ioff][3+joff] =
      img->mpr[3+ioff][2+joff] =
      img->mpr[5+ioff][1+joff] =
      img->mpr[7+ioff][0+joff] = (P_V + P_T + 2*P_U + 2) >> 2;
    img->mpr[1+ioff][2+joff] =
      img->mpr[3+ioff][1+joff] =
      img->mpr[5+ioff][0+joff] = (P_U + P_S + 2*P_T + 2) >> 2;
    img->mpr[1+ioff][1+joff] =
      img->mpr[3+ioff][0+joff] = (P_T + P_R + 2*P_S + 2) >> 2;
    img->mpr[1+ioff][0+joff] = (P_S + P_Q + 2*P_R + 2) >> 2;
    break;
    
  default:
    printf("Error: illegal intra_4x4 prediction mode: %d\n",predmode);
    return SEARCH_SYNC;
    break;
  }
  return DECODING_OK;
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
 ***********************************************************************
 * \brief
 *    Inverse 8x8 transformation
 ***********************************************************************
 */
void itrans8x8(struct img_par *img, //!< image parameters
              int ioff,            //!< index to 4x4 block
              int joff)            //!<
{
  int i,j;
  int m6[8][8];
  Boolean lossless_qpprime = (Boolean)((img->qp + img->bitdepth_luma_qp_scale)==0 && img->lossless_qpprime_flag==1);

  int residue_transform_flag = img->residue_transform_flag;


  for( i=0; i<8 && !lossless_qpprime; i++)
  {
    int a[8], b[8];
    a[0] = img->m7[ioff + 0][joff + i] + img->m7[ioff + 4][joff + i];
    a[4] = img->m7[ioff + 0][joff + i] - img->m7[ioff + 4][joff + i];
    a[2] = (img->m7[ioff + 2][joff + i]>>1) - img->m7[ioff + 6][joff + i];
    a[6] = img->m7[ioff + 2][joff + i] + (img->m7[ioff + 6][joff + i]>>1);

    b[0] = a[0] + a[6];
    b[2] = a[4] + a[2];
    b[4] = a[4] - a[2];
    b[6] = a[0] - a[6];

    a[1] = -img->m7[ioff + 3][joff + i] + img->m7[ioff + 5][joff + i] - img->m7[ioff + 7][joff + i] - (img->m7[ioff + 7][joff + i]>>1);
    a[3] = img->m7[ioff + 1][joff + i] + img->m7[ioff + 7][joff + i] - img->m7[ioff + 3][joff + i] - (img->m7[ioff + 3][joff + i]>>1);
    a[5] = -img->m7[ioff + 1][joff + i] + img->m7[ioff + 7][joff + i] + img->m7[ioff + 5][joff + i] + (img->m7[ioff + 5][joff + i]>>1);
    a[7] = img->m7[ioff + 3][joff + i] + img->m7[ioff + 5][joff + i] + img->m7[ioff + 1][joff + i] + (img->m7[ioff + 1][joff + i]>>1);

    b[1] = a[1] + (a[7]>>2);
//    b[7] = -(a[1]>>2 + 0) + a[7];  KS: do we need to add zero?
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
  for( i=0; i<8 && !lossless_qpprime; i++)
  {
    int a[8], b[8];
    a[0] = m6[i][0] + m6[i][4];
    a[4] = m6[i][0] - m6[i][4];
    a[2] = (m6[i][2]>>1) - m6[i][6];
    a[6] = m6[i][2] + (m6[i][6]>>1);

    b[0] = a[0] + a[6];
    b[2] = a[4] + a[2];
    b[4] = a[4] - a[2];
    b[6] = a[0] - a[6];

    a[1] = -m6[i][3] + m6[i][5] - m6[i][7] - (m6[i][7]>>1);
    a[3] = m6[i][1] + m6[i][7] - m6[i][3] - (m6[i][3]>>1);
    a[5] = -m6[i][1] + m6[i][7] + m6[i][5] + (m6[i][5]>>1);
    a[7] = m6[i][3] + m6[i][5] + m6[i][1] + (m6[i][1]>>1);

    b[1] = a[1] + (a[7]>>2);
    b[7] = -(a[1]>>2) + a[7];
    b[3] = a[3] + (a[5]>>2);
    b[5] = (a[3]>>2) - a[5];

    img->m7[ioff + i][joff + 0] = b[0] + b[7];
    img->m7[ioff + i][joff + 1] = b[2] + b[5];
    img->m7[ioff + i][joff + 2] = b[4] + b[3];
    img->m7[ioff + i][joff + 3] = b[6] + b[1];
    img->m7[ioff + i][joff + 4] = b[6] - b[1];
    img->m7[ioff + i][joff + 5] = b[4] - b[3];
    img->m7[ioff + i][joff + 6] = b[2] - b[5];
    img->m7[ioff + i][joff + 7] = b[0] - b[7];
  }
  for( i=0; i<8; i++)
  {
    for( j=0; j<8; j++)
    {
      // Residue Color Transform
      if(!residue_transform_flag)
      {
        if(lossless_qpprime)
          img->m7[i+ioff][j+joff] =min(img->max_imgpel_value,max(0,img->m7[ioff + i][joff + j]+(long)img->mpr[i+ioff][j+joff]));
        else
          img->m7[i+ioff][j+joff] =min(img->max_imgpel_value,max(0,(img->m7[ioff + i][joff + j]+((long)img->mpr[i+ioff][j+joff] << DQ_BITS_8)+DQ_ROUND_8)>>DQ_BITS_8));
      }
      else
      {
        if(lossless_qpprime)
          img->m7[i+ioff][j+joff] = img->m7[ioff + i][joff + j];
        else
          img->m7[i+ioff][j+joff] =(img->m7[ioff + i][joff + j]+DQ_ROUND_8)>>DQ_BITS_8;
      }
    }
  }
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


static int KLTCol8inProd[18][8], KLTRow8inProd[18][8];

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


  for(i = 0; i < 18; i++)
  {
    if(i == 2 || i == 9 || i > 11) // handled by DCT 
      continue;

    for(j = 0; j < 8; j++)
      for(k = 0; k < 8; k++)
        temp2D[j][k] = KLTCol8[i][k][j];
    compute_inner_product8x8(temp2D, 1, KLTCol8inProd[i]);

    for(j = 0; j < 8; j++)
      for(k = 0; k < 8; k++)
        temp2D[j][k] = KLTRow8[i][k][j];
    compute_inner_product8x8(temp2D, 0, KLTRow8inProd[i]);
  }
}

/*! 
*************************************************************************************
* \brief
*   Perform fast mode-dependent inverse 8x8 directional transform 
*
* \para itrans8x8klt_sep_fast()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
void itrans8x8klt_sep_fast( struct img_par *img, //!< image parameters
                            int ioff,            //!< index to 4x4 block
                            int joff,            //!<
                            int ipmode)
{
  int i, j, k; 
  int shift = 20, shiftby2 = 1<<(shift-1);
  int inprod[BLOCK_SIZE*2], 
      ilev  [BLOCK_SIZE*2][BLOCK_SIZE*2],
      temp2D[BLOCK_SIZE*2][BLOCK_SIZE*2], 
      y2D   [BLOCK_SIZE*2][BLOCK_SIZE*2];

  for(j = 0; j < 8; j++)
    for(i = 0; i < 8; i++)
    {
      ilev[j][i] = img->m7[ioff+i][joff+j];
    }

  compute_inner_product8x8(ilev, 0, inprod);

  // inverse transform, column first
  for (i=0;i<(BLOCK_SIZE*2);i++)
  {
    for (j=0;j<(BLOCK_SIZE*2);j++)
    {
      temp2D[i][j] = 0;
      for(k = 0; k < BLOCK_SIZE; k++)
        temp2D[i][j]+=
        (KLTCol8[ipmode][2*k][i]+ilev[2*k+1][j])*(KLTCol8[ipmode][2*k+1][i]+ilev[2*k][j]);
      temp2D[i][j] = temp2D[i][j] - KLTCol8inProd[ipmode][i] - inprod[j];
    }
  }

  compute_inner_product8x8(temp2D, 1, inprod);

  // inverse transform, row second 
  for (i=0;i<(BLOCK_SIZE*2);i++)
  {
    for (j=0;j<(BLOCK_SIZE*2);j++)
    {
      y2D[i][j] = 0;
      for(k = 0; k < BLOCK_SIZE; k++)
        y2D[i][j]+=
        (temp2D[i][2*k+1]+KLTRow8[ipmode][j][2*k])*(temp2D[i][2*k]+KLTRow8[ipmode][j][2*k+1]);
      y2D[i][j] = y2D[i][j] - inprod[i] - KLTRow8inProd[ipmode][j];
      
      img->m7[j+ioff][i+joff] =min(img->max_imgpel_value,max(0,(y2D[i][j]+((long)img->mpr[j+ioff][i+joff] << shift)+shiftby2)>>shift));
    }
  }
}

#endif
