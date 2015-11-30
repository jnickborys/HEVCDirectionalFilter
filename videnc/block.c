
/*!
 *************************************************************************************
 * \file block.c
 *
 * \brief
 *    Process one block
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
 *    - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *    - Stephan Wenger                  <stewe@cs.tu-berlin.de>
 *    - Jani Lainema                    <jani.lainema@nokia.com>
 *    - Detlev Marpe                    <marpe@hhi.de>
 *    - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *    - Ragip Kurceren                  <ragip.kurceren@nokia.com>
 *    - Greg Conklin                    <gregc@real.com>
 *************************************************************************************
 */

#include "contributors.h"


#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <string.h>
#include <math.h>

#include "global.h"

#include "image.h"
#include "mb_access.h"
#include "block.h"
#include "vlc.h"

#ifdef USE_INTRA_MDDT
long        quant_stat_rd[16];
long        quant_stat_rd16x16[256];
#endif 

const int quant_coef[6][4][4] = {
  {{13107, 8066,13107, 8066},{ 8066, 5243, 8066, 5243},{13107, 8066,13107, 8066},{ 8066, 5243, 8066, 5243}},
  {{11916, 7490,11916, 7490},{ 7490, 4660, 7490, 4660},{11916, 7490,11916, 7490},{ 7490, 4660, 7490, 4660}},
  {{10082, 6554,10082, 6554},{ 6554, 4194, 6554, 4194},{10082, 6554,10082, 6554},{ 6554, 4194, 6554, 4194}},
  {{ 9362, 5825, 9362, 5825},{ 5825, 3647, 5825, 3647},{ 9362, 5825, 9362, 5825},{ 5825, 3647, 5825, 3647}},
  {{ 8192, 5243, 8192, 5243},{ 5243, 3355, 5243, 3355},{ 8192, 5243, 8192, 5243},{ 5243, 3355, 5243, 3355}},
  {{ 7282, 4559, 7282, 4559},{ 4559, 2893, 4559, 2893},{ 7282, 4559, 7282, 4559},{ 4559, 2893, 4559, 2893}}
};

const int dequant_coef[6][4][4] = {
  {{10, 13, 10, 13},{ 13, 16, 13, 16},{10, 13, 10, 13},{ 13, 16, 13, 16}},
  {{11, 14, 11, 14},{ 14, 18, 14, 18},{11, 14, 11, 14},{ 14, 18, 14, 18}},
  {{13, 16, 13, 16},{ 16, 20, 16, 20},{13, 16, 13, 16},{ 16, 20, 16, 20}},
  {{14, 18, 14, 18},{ 18, 23, 18, 23},{14, 18, 14, 18},{ 18, 23, 18, 23}},
  {{16, 20, 16, 20},{ 20, 25, 20, 25},{16, 20, 16, 20},{ 20, 25, 20, 25}},
  {{18, 23, 18, 23},{ 23, 29, 23, 29},{18, 23, 18, 23},{ 23, 29, 23, 29}}
};
static const int A[4][4] = {
  { 16, 20, 16, 20},
  { 20, 25, 20, 25},
  { 16, 20, 16, 20},
  { 20, 25, 20, 25}
};


const byte QP_SCALE_CR[52]=
{
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,
  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,
  28,29,29,30,31,32,32,33,34,34,35,35,36,36,37,37,
  37,38,38,38,39,39,39,39
};

#ifdef RDO_Q
const int estErr4x4[6][4][4]=
{
  {
    {25600, 27040, 25600, 27040}, 
    {27040, 25600, 27040, 25600}, 
    {25600, 27040, 25600, 27040}, 
    {27040, 25600, 27040, 25600} 
  },
  {
    {30976, 31360, 30976, 31360}, 
    {31360, 32400, 31360, 32400}, 
    {30976, 31360, 30976, 31360}, 
    {31360, 32400, 31360, 32400} 
  },
  {
    {43264, 40960, 43264, 40960}, 
    {40960, 40000, 40960, 40000}, 
    {43264, 40960, 43264, 40960}, 
    {40960, 40000, 40960, 40000} 
  },
  {
    {50176, 51840, 50176, 51840}, 
    {51840, 52900, 51840, 52900}, 
    {50176, 51840, 50176, 51840}, 
    {51840, 52900, 51840, 52900} 
  },
  {
    {65536, 64000, 65536, 64000}, 
    {64000, 62500, 64000, 62500}, 
    {65536, 64000, 65536, 64000}, 
    {64000, 62500, 64000, 62500} 
  },
  {
    {82944, 84640, 82944, 84640}, 
    {84640, 84100, 84640, 84100}, 
    {82944, 84640, 82944, 84640}, 
    {84640, 84100, 84640, 84100} 
  }
};
#endif


// Notation for comments regarding prediction and predictors.
// The pels of the 4x4 block are labelled a..p. The predictor pels above
// are labelled A..H, from the left I..P, and from above left X, as follows:
//
//  X A B C D E F G H
//  I a b c d
//  J e f g h
//  K i j k l
//  L m n o p
//

// Predictor array index definitions
#define P_X (PredPel[0])
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

/*!
 ************************************************************************
 * \brief
 *    Make intra 4x4 prediction according to all 9 prediction modes.
 *    The routine uses left and upper neighbouring points from
 *    previous coded blocks to do this (if available). Notice that
 *    inaccessible neighbouring points are signalled with a negative
 *    value in the predmode array .
 *
 *  \par Input:
 *     Starting point of current 4x4 block image posision
 *
 *  \par Output:
 *      none
 ************************************************************************
 */
void intrapred_luma(int img_x,int img_y, int *left_available, int *up_available, int *all_available)
{
  int i,j;
  int s0;
  int PredPel[13];  // array of predictor pels
  imgpel **imgY = enc_picture->imgY;  // For MB level frame/field coding tools -- set default to imgY

  int ioff = (img_x & 15);
  int joff = (img_y & 15);
  int mb_nr=img->current_mb_nr;

  PixelPos pix_a[4];
  PixelPos pix_b, pix_c, pix_d;

  int block_available_up;
  int block_available_left;
  int block_available_up_left;
  int block_available_up_right;

  for (i=0;i<4;i++)
  {
    getNeighbour(mb_nr, ioff -1 , joff +i , 1, &pix_a[i]);
  }
    
  getNeighbour(mb_nr, ioff    , joff -1 , 1, &pix_b);
  getNeighbour(mb_nr, ioff +4 , joff -1 , 1, &pix_c);
  getNeighbour(mb_nr, ioff -1 , joff -1 , 1, &pix_d);

  pix_c.available = pix_c.available && !((ioff==4) && ((joff==4)||(joff==12)));

  if (input->UseConstrainedIntraPred)
  {
    for (i=0, block_available_left=1; i<4;i++)
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
  if (block_available_up)
  {
    P_A = imgY[pix_b.pos_y][pix_b.pos_x+0];
    P_B = imgY[pix_b.pos_y][pix_b.pos_x+1];
    P_C = imgY[pix_b.pos_y][pix_b.pos_x+2];
    P_D = imgY[pix_b.pos_y][pix_b.pos_x+3];

  }
  else
  {
    P_A = P_B = P_C = P_D = img->dc_pred_value_luma;
  }

  if (block_available_up_right)
  {
    P_E = imgY[pix_c.pos_y][pix_c.pos_x+0];
    P_F = imgY[pix_c.pos_y][pix_c.pos_x+1];
    P_G = imgY[pix_c.pos_y][pix_c.pos_x+2];
    P_H = imgY[pix_c.pos_y][pix_c.pos_x+3];
  }
  else
  {
    P_E = P_F = P_G = P_H = P_D;
  }

  if (block_available_left)
  {
    P_I = imgY[pix_a[0].pos_y][pix_a[0].pos_x];
    P_J = imgY[pix_a[1].pos_y][pix_a[1].pos_x];
    P_K = imgY[pix_a[2].pos_y][pix_a[2].pos_x];
    P_L = imgY[pix_a[3].pos_y][pix_a[3].pos_x];
  }
  else
  {
    P_I = P_J = P_K = P_L = img->dc_pred_value_luma;
  }

  if (block_available_up_left)
  {
    P_X = imgY[pix_d.pos_y][pix_d.pos_x];
  }
  else
  {
    P_X = img->dc_pred_value_luma;
  }

  for(i=0;i<9;i++)
    img->mprr[i][0][0]=-1;

  ///////////////////////////////
  // make DC prediction
  ///////////////////////////////
  s0 = 0;
  if (block_available_up && block_available_left)
  {   
    // no edge
    s0 = (P_A + P_B + P_C + P_D + P_I + P_J + P_K + P_L + 4) >> (BLOCK_SHIFT + 1);
  }
  else if (!block_available_up && block_available_left)
  {
    // upper edge
    s0 = (P_I + P_J + P_K + P_L + 2) >> BLOCK_SHIFT;;             
  }
  else if (block_available_up && !block_available_left)
  {
    // left edge
    s0 = (P_A + P_B + P_C + P_D + 2) >> BLOCK_SHIFT;             
  }
  else //if (!block_available_up && !block_available_left)
  {
    // top left corner, nothing to predict from
    s0 = img->dc_pred_value_luma;                           
  }

  // store DC prediction
  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < BLOCK_SIZE; i++)
      img->mprr[DC_PRED][j][i] = s0;
  }

  ///////////////////////////////
  // make horiz and vert prediction
  ///////////////////////////////

  for (i=0; i < BLOCK_SIZE; i++)
  {
    img->mprr[VERT_PRED][0][i] = 
    img->mprr[VERT_PRED][1][i] = 
    img->mprr[VERT_PRED][2][i] = 
    img->mprr[VERT_PRED][3][i] = (&P_A)[i];
    img->mprr[HOR_PRED][i][0]  = 
    img->mprr[HOR_PRED][i][1]  = 
    img->mprr[HOR_PRED][i][2]  = 
    img->mprr[HOR_PRED][i][3]  = (&P_I)[i];
  }

  if(!block_available_up)
    img->mprr[VERT_PRED][0][0]=-1;
  if(!block_available_left)
    img->mprr[HOR_PRED][0][0]=-1;

  if (block_available_up) 
  {
    // Mode DIAG_DOWN_LEFT_PRED
    img->mprr[DIAG_DOWN_LEFT_PRED][0][0] = (P_A + P_C + 2*(P_B) + 2) >> 2;
    img->mprr[DIAG_DOWN_LEFT_PRED][0][1] = 
    img->mprr[DIAG_DOWN_LEFT_PRED][1][0] = (P_B + P_D + 2*(P_C) + 2) >> 2;
    img->mprr[DIAG_DOWN_LEFT_PRED][0][2] =
    img->mprr[DIAG_DOWN_LEFT_PRED][1][1] =
    img->mprr[DIAG_DOWN_LEFT_PRED][2][0] = (P_C + P_E + 2*(P_D) + 2) >> 2;
    img->mprr[DIAG_DOWN_LEFT_PRED][0][3] = 
    img->mprr[DIAG_DOWN_LEFT_PRED][1][2] = 
    img->mprr[DIAG_DOWN_LEFT_PRED][2][1] = 
    img->mprr[DIAG_DOWN_LEFT_PRED][3][0] = (P_D + P_F + 2*(P_E) + 2) >> 2;
    img->mprr[DIAG_DOWN_LEFT_PRED][1][3] = 
    img->mprr[DIAG_DOWN_LEFT_PRED][2][2] = 
    img->mprr[DIAG_DOWN_LEFT_PRED][3][1] = (P_E + P_G + 2*(P_F) + 2) >> 2;
    img->mprr[DIAG_DOWN_LEFT_PRED][2][3] = 
    img->mprr[DIAG_DOWN_LEFT_PRED][3][2] = (P_F + P_H + 2*(P_G) + 2) >> 2;
    img->mprr[DIAG_DOWN_LEFT_PRED][3][3] = (P_G + 3*(P_H) + 2) >> 2;

    // Mode VERT_LEFT_PRED
    img->mprr[VERT_LEFT_PRED][0][0] = (P_A + P_B + 1) >> 1;
    img->mprr[VERT_LEFT_PRED][0][1] = 
    img->mprr[VERT_LEFT_PRED][2][0] = (P_B + P_C + 1) >> 1;
    img->mprr[VERT_LEFT_PRED][0][2] = 
    img->mprr[VERT_LEFT_PRED][2][1] = (P_C + P_D + 1) >> 1;
    img->mprr[VERT_LEFT_PRED][0][3] = 
    img->mprr[VERT_LEFT_PRED][2][2] = (P_D + P_E + 1) >> 1;
    img->mprr[VERT_LEFT_PRED][2][3] = (P_E + P_F + 1) >> 1;
    img->mprr[VERT_LEFT_PRED][1][0] = (P_A + 2*P_B + P_C + 2) >> 2;
    img->mprr[VERT_LEFT_PRED][1][1] = 
    img->mprr[VERT_LEFT_PRED][3][0] = (P_B + 2*P_C + P_D + 2) >> 2;
    img->mprr[VERT_LEFT_PRED][1][2] = 
    img->mprr[VERT_LEFT_PRED][3][1] = (P_C + 2*P_D + P_E + 2) >> 2;
    img->mprr[VERT_LEFT_PRED][1][3] = 
    img->mprr[VERT_LEFT_PRED][3][2] = (P_D + 2*P_E + P_F + 2) >> 2;
    img->mprr[VERT_LEFT_PRED][3][3] = (P_E + 2*P_F + P_G + 2) >> 2;

  }

  /*  Prediction according to 'diagonal' modes */
  if (block_available_left) 
  {
    // Mode HOR_UP_PRED
    img->mprr[HOR_UP_PRED][0][0] = (P_I + P_J + 1) >> 1;
    img->mprr[HOR_UP_PRED][0][1] = (P_I + 2*P_J + P_K + 2) >> 2;
    img->mprr[HOR_UP_PRED][0][2] = 
    img->mprr[HOR_UP_PRED][1][0] = (P_J + P_K + 1) >> 1;
    img->mprr[HOR_UP_PRED][0][3] = 
    img->mprr[HOR_UP_PRED][1][1] = (P_J + 2*P_K + P_L + 2) >> 2;
    img->mprr[HOR_UP_PRED][1][2] = 
    img->mprr[HOR_UP_PRED][2][0] = (P_K + P_L + 1) >> 1;
    img->mprr[HOR_UP_PRED][1][3] = 
    img->mprr[HOR_UP_PRED][2][1] = (P_K + 2*P_L + P_L + 2) >> 2;
    img->mprr[HOR_UP_PRED][3][0] = 
    img->mprr[HOR_UP_PRED][2][2] = 
    img->mprr[HOR_UP_PRED][2][3] = 
    img->mprr[HOR_UP_PRED][3][1] = 
    img->mprr[HOR_UP_PRED][3][2] = 
    img->mprr[HOR_UP_PRED][3][3] = P_L;
  }

  /*  Prediction according to 'diagonal' modes */
  if (block_available_up && block_available_left && block_available_up_left) 
  {
    // Mode DIAG_DOWN_RIGHT_PRED
    img->mprr[DIAG_DOWN_RIGHT_PRED][3][0] = (P_L + 2*P_K + P_J + 2) >> 2; 
    img->mprr[DIAG_DOWN_RIGHT_PRED][2][0] =
    img->mprr[DIAG_DOWN_RIGHT_PRED][3][1] = (P_K + 2*P_J + P_I + 2) >> 2; 
    img->mprr[DIAG_DOWN_RIGHT_PRED][1][0] =
    img->mprr[DIAG_DOWN_RIGHT_PRED][2][1] = 
    img->mprr[DIAG_DOWN_RIGHT_PRED][3][2] = (P_J + 2*P_I + P_X + 2) >> 2; 
    img->mprr[DIAG_DOWN_RIGHT_PRED][0][0] =
    img->mprr[DIAG_DOWN_RIGHT_PRED][1][1] =
    img->mprr[DIAG_DOWN_RIGHT_PRED][2][2] =
    img->mprr[DIAG_DOWN_RIGHT_PRED][3][3] = (P_I + 2*P_X + P_A + 2) >> 2; 
    img->mprr[DIAG_DOWN_RIGHT_PRED][0][1] =
    img->mprr[DIAG_DOWN_RIGHT_PRED][1][2] =
    img->mprr[DIAG_DOWN_RIGHT_PRED][2][3] = (P_X + 2*P_A + P_B + 2) >> 2;
    img->mprr[DIAG_DOWN_RIGHT_PRED][0][2] =
    img->mprr[DIAG_DOWN_RIGHT_PRED][1][3] = (P_A + 2*P_B + P_C + 2) >> 2;
    img->mprr[DIAG_DOWN_RIGHT_PRED][0][3] = (P_B + 2*P_C + P_D + 2) >> 2;

     // Mode VERT_RIGHT_PRED
    img->mprr[VERT_RIGHT_PRED][0][0] = 
    img->mprr[VERT_RIGHT_PRED][2][1] = (P_X + P_A + 1) >> 1;
    img->mprr[VERT_RIGHT_PRED][0][1] = 
    img->mprr[VERT_RIGHT_PRED][2][2] = (P_A + P_B + 1) >> 1;
    img->mprr[VERT_RIGHT_PRED][0][2] = 
    img->mprr[VERT_RIGHT_PRED][2][3] = (P_B + P_C + 1) >> 1;
    img->mprr[VERT_RIGHT_PRED][0][3] = (P_C + P_D + 1) >> 1;
    img->mprr[VERT_RIGHT_PRED][1][0] = 
    img->mprr[VERT_RIGHT_PRED][3][1] = (P_I + 2*P_X + P_A + 2) >> 2;
    img->mprr[VERT_RIGHT_PRED][1][1] = 
    img->mprr[VERT_RIGHT_PRED][3][2] = (P_X + 2*P_A + P_B + 2) >> 2;
    img->mprr[VERT_RIGHT_PRED][1][2] = 
    img->mprr[VERT_RIGHT_PRED][3][3] = (P_A + 2*P_B + P_C + 2) >> 2;
    img->mprr[VERT_RIGHT_PRED][1][3] = (P_B + 2*P_C + P_D + 2) >> 2;
    img->mprr[VERT_RIGHT_PRED][2][0] = (P_X + 2*P_I + P_J + 2) >> 2;
    img->mprr[VERT_RIGHT_PRED][3][0] = (P_I + 2*P_J + P_K + 2) >> 2;

    // Mode HOR_DOWN_PRED
    img->mprr[HOR_DOWN_PRED][0][0] = 
    img->mprr[HOR_DOWN_PRED][1][2] = (P_X + P_I + 1) >> 1;
    img->mprr[HOR_DOWN_PRED][0][1] = 
    img->mprr[HOR_DOWN_PRED][1][3] = (P_I + 2*P_X + P_A + 2) >> 2;
    img->mprr[HOR_DOWN_PRED][0][2] = (P_X + 2*P_A + P_B + 2) >> 2;
    img->mprr[HOR_DOWN_PRED][0][3] = (P_A + 2*P_B + P_C + 2) >> 2;
    img->mprr[HOR_DOWN_PRED][1][0] = 
    img->mprr[HOR_DOWN_PRED][2][2] = (P_I + P_J + 1) >> 1;
    img->mprr[HOR_DOWN_PRED][1][1] = 
    img->mprr[HOR_DOWN_PRED][2][3] = (P_X + 2*P_I + P_J + 2) >> 2;
    img->mprr[HOR_DOWN_PRED][2][0] = 
    img->mprr[HOR_DOWN_PRED][3][2] = (P_J + P_K + 1) >> 1;
    img->mprr[HOR_DOWN_PRED][2][1] = 
    img->mprr[HOR_DOWN_PRED][3][3] = (P_I + 2*P_J + P_K + 2) >> 2;
    img->mprr[HOR_DOWN_PRED][3][0] = (P_K + P_L + 1) >> 1;
    img->mprr[HOR_DOWN_PRED][3][1] = (P_J + 2*P_K + P_L + 2) >> 2;
  }
}

/*!
 ************************************************************************
 * \brief
 *    16x16 based luma prediction
 *
 * \par Input:
 *    Image parameters
 *
 * \par Output:
 *    none
 ************************************************************************
 */
void intrapred_luma_16x16()
{
  int s0=0,s1,s2;
  imgpel s[2][16];
  int i,j;

  int ih,iv;
  int ib,ic,iaa;

  imgpel   **imgY_pred = enc_picture->imgY;  // For Mb level field/frame coding tools -- default to frame pred
  int          mb_nr = img->current_mb_nr;

  PixelPos up;          //!< pixel position p(0,-1)
  PixelPos left[17];    //!< pixel positions p(-1, -1..15)

  int up_avail, left_avail, left_up_avail;

  for (i=0;i<17;i++)
  {
    getNeighbour(mb_nr, -1,  i-1, 1, &left[i]);
  }
  
  getNeighbour(mb_nr,    0,   -1, 1, &up);

  if (!(input->UseConstrainedIntraPred))
  {
    up_avail      = up.available;
    left_avail    = left[1].available;
    left_up_avail = left[0].available;
  }
  else
  {
    up_avail      = up.available ? img->intra_block[up.mb_addr] : 0;
    for (i=1, left_avail=1; i<17;i++)
      left_avail  &= left[i].available ? img->intra_block[left[i].mb_addr]: 0;
    left_up_avail = left[0].available ? img->intra_block[left[0].mb_addr]: 0;
  }

  s1=s2=0;
  // make DC prediction
  if (up_avail)
  {
    for (i=0; i < MB_BLOCK_SIZE; i++)
      s1 += imgY_pred[up.pos_y][up.pos_x+i];    // sum hor pix
  }

  if (left_avail)
  {
    for (i=0; i < MB_BLOCK_SIZE; i++)      
      s2 += imgY_pred[left[i+1].pos_y][left[i+1].pos_x];    // sum vert pix
  }

  if (up_avail && left_avail)
    s0=(s1+s2+16)/(2*MB_BLOCK_SIZE);             // no edge
  
  if (!up_avail && left_avail)
    s0=(s2+8)/MB_BLOCK_SIZE;                     // upper edge
  
  if (up_avail && !left_avail)
    s0=(s1+8)/MB_BLOCK_SIZE;                     // left edge
  
  if (!up_avail && !left_avail)
    s0=img->dc_pred_value_luma;                       // top left corner, nothing to predict from

  // vertical prediction
  if (up_avail)
    memcpy(s[0], &imgY_pred[up.pos_y][up.pos_x], MB_BLOCK_SIZE * sizeof(imgpel));
  
  // horizontal prediction
  if (left_avail)
  {
    for (i=0; i < MB_BLOCK_SIZE; i++)
      s[1][i]=imgY_pred[left[i+1].pos_y][left[i+1].pos_x];
  }

  for (j=0; j < MB_BLOCK_SIZE; j++)
  {
    memcpy(img->mprr_2[VERT_PRED_16][j], s[0], MB_BLOCK_SIZE * sizeof(imgpel)); // store vertical prediction
    for (i=0; i < MB_BLOCK_SIZE; i++)
    {
      img->mprr_2[HOR_PRED_16 ][j][i]=s[1][j]; // store horizontal prediction
      img->mprr_2[DC_PRED_16  ][j][i]=s0;      // store DC prediction
    }
  }
  if (!up_avail || !left_avail || !left_up_avail) // edge
    return;

  // 16 bit integer plan pred

  ih=0;
  iv=0;
  for (i=1;i<9;i++)
  {
    if (i<8)
      ih += i*(imgY_pred[up.pos_y][up.pos_x+7+i] - imgY_pred[up.pos_y][up.pos_x+7-i]);
    else
      ih += i*(imgY_pred[up.pos_y][up.pos_x+7+i] - imgY_pred[left[0].pos_y][left[0].pos_x]);
    
    iv += i*(imgY_pred[left[8+i].pos_y][left[8+i].pos_x] - imgY_pred[left[8-i].pos_y][left[8-i].pos_x]);
  }
  ib=(5*ih+32)>>6;
  ic=(5*iv+32)>>6;
  
  iaa=16*(imgY_pred[up.pos_y][up.pos_x+15]+imgY_pred[left[16].pos_y][left[16].pos_x]);

  for (j=0;j< MB_BLOCK_SIZE;j++)
  {
    for (i=0;i< MB_BLOCK_SIZE;i++)
    {
      img->mprr_2[PLANE_16][j][i]=max(0,min((int)img->max_imgpel_value,(iaa+(i-7)*ib +(j-7)*ic + 16)/32));// store plane prediction
    }
  }
}

#ifdef RDO_Q //TREL_CAVLC
typedef struct TRELLISNODE
{
  int level;
  int level_idx;
  struct TRELLISNODE *prev;

} TrellisNode;

int estCAVLCbits (int m4[4][4], int level_to_enc[16], int nnz, int block_type, int b8, int b4, int param);
int predict_nnz(int i,int j);
int predict_nnz_chroma(int i,int j);


int cmp(const void *arg1, const void *arg2)
{
  return (int)(((levelDataStruct *)arg2)->levelDouble - ((levelDataStruct *)arg1)->levelDouble);
}

void TrellisCAVLC_Q(int m4[4][4], levelDataStruct *levelData, int *levelTrellis, int block_type, int b8, int b4, int coeff_num, double lambda)
{
  int k, lastnonzero, coeff_ctr, dumb=0;
  int level_to_enc[16];
  int subblock_x, subblock_y, nnz;
  int cstat, bestcstat=0; 
  int nz_coeff=0;
  double lagr, lagrAcc, minlagr=0;


  if(block_type == CHROMA_AC)
  {
     
    static unsigned char chroma_ac_param[3][8][4] =
    {
      {{ 4, 20,  5, 21},
      {36, 52, 37, 53},
      { 0,  0,  0,  0},
      { 0,  0,  0,  0},
      { 0,  0,  0,  0},
      { 0,  0,  0,  0},
      { 0,  0,  0,  0},
      { 0,  0,  0,  0}},

      {{ 4, 20,  5, 21},
      { 6, 22,  7, 23},
      {36, 52, 37, 53},
      {38, 54, 39, 55},
      { 0,  0,  0,  0},
      { 0,  0,  0,  0},
      { 0,  0,  0,  0},
      { 0,  0,  0,  0}},

      {{ 4, 20,  5, 21},
      {36, 52, 37, 53},
      { 6, 22,  7, 23},
      {38, 54, 39, 55},
      { 8, 24,  9, 25},
      {40, 56, 41, 57},
      {10, 26, 11, 27},
      {42, 58, 43, 59}}
    };

    int param = chroma_ac_param[img->yuv_format-1][b8][b4];

    // chroma AC
    subblock_x = param >> 4;
    subblock_y = param & 15;
    nnz = predict_nnz_chroma(subblock_x,subblock_y);
  }
  else
  {
    subblock_x = ((b8&0x1)==0)?(((b4&0x1)==0)?0:1):(((b4&0x1)==0)?2:3); 
    // horiz. position for coeff_count context  
    subblock_y = (b8<2)?((b4<2)?0:1):((b4<2)?2:3); 
    // vert.  position for coeff_count context      
    nnz = predict_nnz(subblock_x,subblock_y);
  }

  lastnonzero = -1;
  lagrAcc=0;
  for (coeff_ctr=0;coeff_ctr < coeff_num;coeff_ctr++)
  {	 
    levelTrellis[coeff_ctr] = 0;

    for(k=0; k<levelData[coeff_ctr].noLevels; k++)
    {
      levelData[coeff_ctr].errLevel[k] /= 32768;
    }
           
    lagrAcc += levelData[coeff_ctr].errLevel[levelData[coeff_ctr].noLevels-1];

    level_to_enc[coeff_ctr] = levelData[coeff_ctr].pre_level;

    if(levelData[coeff_ctr].noLevels > 1)
    {
      levelData[coeff_ctr].coeff_ctr = coeff_ctr;
      lastnonzero = coeff_ctr;
    }
  }


  if(lastnonzero != -1)
  {

    //sort the coefficients based on their absolute value
    qsort(levelData, lastnonzero+1, sizeof(levelDataStruct), cmp);

    for(coeff_ctr=lastnonzero; coeff_ctr>=0; coeff_ctr--) // go over all coeff
    {
      if(levelData[coeff_ctr].noLevels == 1)
        continue;
		
      lagrAcc -= levelData[coeff_ctr].errLevel[levelData[coeff_ctr].noLevels-1];
      for(cstat=0; cstat<levelData[coeff_ctr].noLevels; cstat++) // go over all states of cur coeff k
      {		
        level_to_enc[levelData[coeff_ctr].coeff_ctr] = (int) levelData[coeff_ctr].level[cstat];
        lagr = lagrAcc + levelData[coeff_ctr].errLevel[cstat];
			                         
        lagr += lambda*estCAVLCbits(m4, level_to_enc, nnz, block_type, b8, b4, dumb);
			  
        if(cstat==0 || lagr<minlagr)
        {		
          minlagr = lagr;		
          bestcstat = cstat;
        }
      }
      lagrAcc += levelData[coeff_ctr].errLevel[bestcstat];
      level_to_enc[levelData[coeff_ctr].coeff_ctr] = (int)levelData[coeff_ctr].level[bestcstat];
    }

    for(coeff_ctr=0; coeff_ctr<=lastnonzero; coeff_ctr++)
    {
      levelTrellis[coeff_ctr] = level_to_enc[coeff_ctr];
      if(level_to_enc[coeff_ctr] != 0)
        nz_coeff++;
    }
  }
  img->nz_coeff [img->current_mb_nr ][subblock_x][subblock_y] = nz_coeff;
}

void TrellisCAVLC4x4(int m4[4][4], int q_bits, int qp_rem, int **levelscale, int **leveloffset, int *levelTrellis, int block_type, int b8, int b4, int coeff_num, double lambda)
{
  levelDataStruct levelData[16];
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
  double  normFact=pow(2,(2*DQ_BITS+19))*(1<<(2*img->BitDepthIncrease));
#else
  double  normFact=pow(2,(2*DQ_BITS+19));
#endif
  double err;
  int i, j, ii, coeff_ctr, lowerInt, noCoeff;
  int level;
  const byte (*pos_scan)[2] = SNGL_SCAN;

  noCoeff=0;
  for (coeff_ctr=0;coeff_ctr < coeff_num;coeff_ctr++)
  {
    if(block_type == LUMA_INTRA16x16AC || block_type == CHROMA_AC)
    {
      i=pos_scan[coeff_ctr+1][0]; // scan is shifted due to DC
      j=pos_scan[coeff_ctr+1][1]; // scan is shifted due to DC
    }
    else
    {
      i=pos_scan[coeff_ctr][0];
      j=pos_scan[coeff_ctr][1];
    }

    levelData[coeff_ctr].levelDouble=absm(m4[j][i]*levelscale[i][j]);
    level=(int)(levelData[coeff_ctr].levelDouble>>q_bits);

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
      levelData[coeff_ctr].errLevel[ii]=(err*err*(double)estErr4x4[qp_rem][i][j])/normFact; 
    }

    if(levelData[coeff_ctr].noLevels == 1)
      levelData[coeff_ctr].pre_level = 0;
    else
      levelData[coeff_ctr].pre_level = (absm (m4[j][i]) * levelscale[i][j] + leveloffset[i][j]) >> q_bits;
  }

  TrellisCAVLC_Q(m4, levelData, levelTrellis, block_type, b8, b4, coeff_num, lambda);
}
#endif

/*!
 ************************************************************************
 * \brief
 *    For new intra pred routines
 *
 * \par Input:
 *    Image par, 16x16 based intra mode
 *
 * \par Output:
 *    none
 ************************************************************************
 */
int dct_luma_16x16(int new_intra_mode)
{
  //int qp_const;
  int i,j;
  int ii,jj;
  int jdiv, jmod;
  int M1[16][16];
  int M4[4][4];
  int M5[4],M6[4];
  int M0[4][4][4][4];
  int run,scan_pos,coeff_ctr,level;
  int qp_per,qp_rem,q_bits;
  int ac_coef = 0;

  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  short is_field_mode = (img->field_picture || ( img->MbaffFrameFlag && currMB->mb_field));

  int   b8, b4;
  int*  DCLevel = img->cofDC[0][0];
  int*  DCRun   = img->cofDC[0][1];
  int*  ACLevel;
  int*  ACRun;
  int **levelscale,**leveloffset;
  int **invlevelscale;
  Boolean lossless_qpprime = (Boolean)((currMB->qp + img->bitdepth_luma_qp_scale)==0 && img->lossless_qpprime_flag==1);
  const byte (*pos_scan)[2] = is_field_mode ? FIELD_SCAN : SNGL_SCAN;
#ifdef RDO_Q
  levelDataStruct levelData[16];
  double  lambda_md=0;
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
  double normFact = pow(2, (2 * DQ_BITS + 19))*(1<<(2*img->BitDepthIncrease));
#else
  double normFact = pow(2, (2 * DQ_BITS + 19));
#endif
  double err;
  int lowerInt, levelTrellis[16], k, kStart, kStop, noCoeff, estBits;
#endif

  // Note that we could just use currMB->qp here
  qp_per    = qp_per_matrix[(currMB->qp + img->bitdepth_luma_qp_scale - MIN_QP)];
  qp_rem    = qp_rem_matrix[(currMB->qp + img->bitdepth_luma_qp_scale - MIN_QP)];

  q_bits    = Q_BITS+qp_per;

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
    levelscale    = LevelScale4x4Luma_IAQMS[img->mb_iaqms_idx][1][qp_rem];
    invlevelscale = InvLevelScale4x4Luma_IAQMS[img->mb_iaqms_idx][1][qp_rem];
  }
  else
  {
    levelscale    = LevelScale4x4Luma[1][qp_rem];
    invlevelscale = InvLevelScale4x4Luma[1][qp_rem];
  }
  leveloffset   = LevelOffset4x4Luma[1][qp_per];
#else
  levelscale    = LevelScale4x4Luma[1][qp_rem];
  leveloffset   = LevelOffset4x4Luma[1][qp_per];
  invlevelscale = InvLevelScale4x4Luma[1][qp_rem];
#endif

  for (j=0;j<16;j++)
  {
    jdiv = j >> 2;
    jmod = j & 0x03;
    jj = img->opix_y+j;
    for (i=0;i<16;i++)
    {
      // Residue Color Transform
      if(!img->residue_transform_flag)
        M1[j][i]=imgY_org[jj][img->opix_x+i]-img->mprr_2[new_intra_mode][j][i];
      else
        M1[j][i]=img->m7[j][i];

      M0[jdiv][i >> 2][jmod][i & 0x03]=M1[j][i];
    }
  }

  for (jj=0;jj<4 && !lossless_qpprime;jj++)
  {
    for (ii=0;ii<4;ii++)
    {
      for (j=0;j<4;j++)
      {
        M5[0] = M0[jj][ii][j][0] + M0[jj][ii][j][3];
        M5[1] = M0[jj][ii][j][1] + M0[jj][ii][j][2];
        M5[2] = M0[jj][ii][j][1] - M0[jj][ii][j][2];
        M5[3] = M0[jj][ii][j][0] - M0[jj][ii][j][3];

        M4[j][0] = M5[0]   + M5[1];
        M4[j][2] = M5[0]   - M5[1];
        M4[j][1] = M5[3]*2 + M5[2];
        M4[j][3] = M5[3]   - M5[2]*2;
      }
      // vertical
      for (i=0;i<4;i++)
      {
        M5[0] = M4[0][i] + M4[3][i];
        M5[1] = M4[1][i] + M4[2][i];
        M5[2] = M4[1][i] - M4[2][i];
        M5[3] = M4[0][i] - M4[3][i];
        
        M0[jj][ii][0][i] = M5[0]   + M5[1];
        M0[jj][ii][2][i] = M5[0]   - M5[1];
        M0[jj][ii][1][i] = M5[3]*2 + M5[2];
        M0[jj][ii][3][i] = M5[3]   - M5[2]*2;
      }
    }
  }

  // pick out DC coeff

  for (j=0;j<4;j++)
  {
    for (i=0;i<4;i++)
      M4[j][i]= M0[j][i][0][0];
  }

  if (!lossless_qpprime)
  {
    for (j=0;j<4;j++)
    {
      M5[0] = M4[j][0]+M4[j][3];
      M5[1] = M4[j][1]+M4[j][2];
      M5[2] = M4[j][1]-M4[j][2];
      M5[3] = M4[j][0]-M4[j][3];
      
      M4[j][0] = M5[0]+M5[1];
      M4[j][2] = M5[0]-M5[1];
      M4[j][1] = M5[3]+M5[2];
      M4[j][3] = M5[3]-M5[2];
    }
    
    // vertical
    
    for (i=0;i<4;i++)
    {    
      M5[0] = M4[0][i]+M4[3][i];
      M5[1] = M4[1][i]+M4[2][i];
      M5[2] = M4[1][i]-M4[2][i];
      M5[3] = M4[0][i]-M4[3][i];
      
      M4[0][i]=(M5[0]+M5[1])>>1;
      M4[2][i]=(M5[0]-M5[1])>>1;
      M4[1][i]=(M5[3]+M5[2])>>1;
      M4[3][i]=(M5[3]-M5[2])>>1;
    }
  }
  // quant

  run=-1;
  scan_pos=0;

  for (coeff_ctr=0;coeff_ctr<16;coeff_ctr++)
  {
    i=pos_scan[coeff_ctr][0];
    j=pos_scan[coeff_ctr][1];

    run++;

    if(lossless_qpprime)
      level= absm(M4[j][i]);
    else
      level= (absm(M4[j][i]) * levelscale[0][0] + (leveloffset[0][0]<<1)) >> (q_bits+1);

    if (input->symbol_mode == UVLC && img->qp < 10) 
    {
      if (level > CAVLC_LEVEL_LIMIT) 
        level = CAVLC_LEVEL_LIMIT;
    }

    if (level != 0)
    {
      DCLevel[scan_pos] = sign(level,M4[j][i]);
      DCRun  [scan_pos] = run;
      ++scan_pos;
      run=-1;
    }
    if(!lossless_qpprime)
      M4[j][i]=sign(level,M4[j][i]);
  }
  DCLevel[scan_pos]=0;

  // invers DC transform
  for (j=0;j<4 && !lossless_qpprime;j++)
  {
    M6[0]=M4[j][0]+M4[j][2];
    M6[1]=M4[j][0]-M4[j][2];
    M6[2]=M4[j][1]-M4[j][3];
    M6[3]=M4[j][1]+M4[j][3];

    M4[j][0] = M6[0]+M6[3];
    M4[j][1] = M6[1]+M6[2];
    M4[j][2] = M6[1]-M6[2];
    M4[j][3] = M6[0]-M6[3];
  }

  for (i=0;i<4 && !lossless_qpprime;i++)
  {
    
    M6[0]=M4[0][i]+M4[2][i];
    M6[1]=M4[0][i]-M4[2][i];
    M6[2]=M4[1][i]-M4[3][i];
    M6[3]=M4[1][i]+M4[3][i];
    
    if(qp_per<6)
    {
      M0[0][i][0][0] = ((M6[0]+M6[3])*invlevelscale[0][0]+(1<<(5-qp_per)))>>(6-qp_per);
      M0[1][i][0][0] = ((M6[1]+M6[2])*invlevelscale[0][0]+(1<<(5-qp_per)))>>(6-qp_per);
      M0[2][i][0][0] = ((M6[1]-M6[2])*invlevelscale[0][0]+(1<<(5-qp_per)))>>(6-qp_per);
      M0[3][i][0][0] = ((M6[0]-M6[3])*invlevelscale[0][0]+(1<<(5-qp_per)))>>(6-qp_per);
    }
    else
    {
      M0[0][i][0][0] = ((M6[0]+M6[3])*invlevelscale[0][0])<<(qp_per-6);
      M0[1][i][0][0] = ((M6[1]+M6[2])*invlevelscale[0][0])<<(qp_per-6);
      M0[2][i][0][0] = ((M6[1]-M6[2])*invlevelscale[0][0])<<(qp_per-6);
      M0[3][i][0][0] = ((M6[0]-M6[3])*invlevelscale[0][0])<<(qp_per-6);
    }   
  }

  // AC inverse trans/quant for MB
  for (jj=0;jj<4;jj++)
  {
    for (ii=0;ii<4;ii++)
    {
      for (j=0;j<4;j++)
      {
        memcpy(M4[j],M0[jj][ii][j], BLOCK_SIZE * sizeof(int));
      }

      run      = -1;
      scan_pos =  0;
      b8       = 2*(jj >> 1) + (ii >> 1);
      b4       = 2*(jj & 0x01) + (ii & 0x01);
      ACLevel  = img->cofAC [b8][b4][0];
      ACRun    = img->cofAC [b8][b4][1];

#ifdef RDO_Q
//#ifdef TREL_CAVLC
      if(input->UseRDO_Q && active_pps->entropy_coding_mode_flag == UVLC)
        TrellisCAVLC4x4(M4, q_bits, qp_rem, levelscale, leveloffset, levelTrellis, LUMA_INTRA16x16AC, b8, b4, 15, lambda_md);
//#endif
      if(input->UseRDO_Q && active_pps->entropy_coding_mode_flag == CABAC)
      {
        kStart=0; kStop=0; noCoeff=0;
        for (coeff_ctr=0;coeff_ctr < 15;coeff_ctr++)
        {
          i=pos_scan[coeff_ctr+1][0]; // scan is shifted due to DC
          j=pos_scan[coeff_ctr+1][1]; // scan is shifted due to DC

          levelData[coeff_ctr].levelDouble=absm(M4[j][i]*levelscale[i][j]);
          level = (int)(levelData[coeff_ctr].levelDouble >> q_bits);
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

          for (k=0; k<levelData[coeff_ctr].noLevels; k++)
          {
            err=(double)(levelData[coeff_ctr].level[k]<<q_bits)-(double)levelData[coeff_ctr].levelDouble;
            levelData[coeff_ctr].errLevel[k]=(err*err*(double)estErr4x4[qp_rem][i][j])/normFact; 
          }
        }
        estBits=est_write_and_store_CBP_block_bit(currMB, LUMA_16AC);

        est_writeRunLevel_CABAC(levelData, levelTrellis, LUMA_16AC, lambda_md, kStart, kStop, noCoeff, estBits);
      }
#endif

      for (coeff_ctr=1;coeff_ctr<16;coeff_ctr++) // set in AC coeff
      {

        i=pos_scan[coeff_ctr][0];
        j=pos_scan[coeff_ctr][1];

        run++;
#ifdef RDO_Q
        if(input->UseRDO_Q)
        {
          level=levelTrellis[coeff_ctr-1];
        }
        else
        {
          if(lossless_qpprime)
            level= absm( M4[j][i]);
          else          
            level= ( absm( M4[j][i]) * levelscale[i][j] + leveloffset[i][j]) >> q_bits;
        }
#else
        if(lossless_qpprime)
          level= absm( M4[j][i]);
        else          
          level= ( absm( M4[j][i]) * levelscale[i][j] + leveloffset[i][j]) >> q_bits;
#endif


        if (img->AdaptiveRounding)
        {
          if (lossless_qpprime || level == 0 )
          {
            img->fadjust4x4[2][jj*BLOCK_SIZE+j][ii*BLOCK_SIZE+i] = 0;
          }
          else
          {
            img->fadjust4x4[2][jj*BLOCK_SIZE+j][ii*BLOCK_SIZE+i] = 
              (AdaptRndWeight * (absm(M4[j][i]) * levelscale[i][j] - (level << q_bits)) + (1<< (q_bits))) >> (q_bits + 1);
          }
        }

        if (level != 0)
        {
          ac_coef = 15;
          ACLevel[scan_pos] = sign(level,M4[j][i]);
          ACRun  [scan_pos] = run;
          ++scan_pos;
          run=-1;
        }
        
        if(!lossless_qpprime)
        {
          level=sign(level, M4[j][i]);
          if(qp_per<4)
            M4[j][i]=(level*invlevelscale[i][j]+(1<<(3-qp_per)))>>(4-qp_per);
          else
            M4[j][i]=(level*invlevelscale[i][j])<<(qp_per-4);
        }
      }
      ACLevel[scan_pos] = 0;


      // IDCT horizontal
      for (j=0;j<4 && !lossless_qpprime;j++)
      {
        M6[0] = M4[j][0]     +  M4[j][2];
        M6[1] = M4[j][0]     -  M4[j][2];
        M6[2] =(M4[j][1]>>1) -  M4[j][3];
        M6[3] = M4[j][1]     + (M4[j][3]>>1);
        
        M4[j][0] = M6[0] + M6[3];
        M4[j][1] = M6[1] + M6[2];
        M4[j][2] = M6[1] - M6[2];
        M4[j][3] = M6[0] - M6[3];
      }

      // vert
      for (i=0;i<4 && !lossless_qpprime;i++)
      {
        M6[0]= M4[0][i]     +  M4[2][i];
        M6[1]= M4[0][i]     -  M4[2][i];
        M6[2]=(M4[1][i]>>1) -  M4[3][i];
        M6[3]= M4[1][i]     + (M4[3][i]>>1);
        
        M0[jj][ii][0][i] = M6[0] + M6[3];
        M0[jj][ii][1][i] = M6[1] + M6[2];
        M0[jj][ii][2][i] = M6[1] - M6[2];
        M0[jj][ii][3][i] = M6[0] - M6[3];
      }
    }
  }

  // Residue Color Transform
  if(!img->residue_transform_flag)
  {
    for (jj=0;jj<BLOCK_MULTIPLE; jj++)
      for (ii=0;ii<BLOCK_MULTIPLE; ii++)
        for (j=0;j<BLOCK_SIZE;j++)
        {
          memcpy(&M1[jj*BLOCK_SIZE + j][ii*BLOCK_SIZE], M0[jj][ii][j], BLOCK_SIZE * sizeof(int));
        }
  }
  else
  {
    if(lossless_qpprime)
    {
      for (j=0;j<MB_BLOCK_SIZE;j++)    
      {
        jdiv = j >> 2;
        jmod = j & 0x03;
        for (i=0;i<MB_BLOCK_SIZE;i++)
          img->m7[j][i]=M0[jdiv][i >> 2][jmod][i & 0x03];
      }
    }
    else
    {
      for (j=0;j<MB_BLOCK_SIZE;j++)    
      {
        jdiv = j >> 2;
        jmod = j & 0x03;
        for (i=0;i<MB_BLOCK_SIZE;i++)
          img->m7[j][i]=((M0[jdiv][i >> 2][jmod][i & 0x03]+DQ_ROUND)>>DQ_BITS);
      }
    }
  }

  if(!img->residue_transform_flag)
  {
    if(lossless_qpprime)
    {
      for (j=0;j<16;j++)
      {
        jj = img->pix_y+j;
        for (i=0;i<16;i++)
        {
          enc_picture->imgY[jj][img->pix_x+i]=(imgpel)(M1[j][i]+img->mprr_2[new_intra_mode][j][i]);
          if(img->type==SP_SLICE)
            lrec[jj][img->pix_x+i]=-16; //signals an I16 block in the SP frame 
        }
      }
    }
    else
    {
      for (j=0;j<16;j++)
      {
        jj = img->pix_y+j;
        for (i=0;i<16;i++)
        {
          enc_picture->imgY[jj][img->pix_x+i]=(imgpel)clip1a((M1[j][i]+((long)img->mprr_2[new_intra_mode][j][i]<<DQ_BITS)+DQ_ROUND)>>DQ_BITS);
          if(img->type==SP_SLICE)
            lrec[jj][img->pix_x+i]=-16; //signals an I16 block in the SP frame 
        }
      }
    }
  }
  return ac_coef;
}


/*!
************************************************************************
* \brief
*    The routine performs transform,quantization,inverse transform, adds the diff.
*    to the prediction and writes the result to the decoded luma frame. Includes the
*    RD constrained quantization also.
*
* \par Input:
*    block_x,block_y: Block position inside a macro block (0,4,8,12).
*
* \par Output_
*    nonzero: 0 if no levels are nonzero.  1 if there are nonzero levels.             \n
*    coeff_cost: Counter for nonzero coefficients, used to discard expensive levels.
************************************************************************
*/
#ifdef USE_INTRA_MDDT
int dct_luma(int block_x,int block_y,int *coeff_cost, int intra, int ipmode)
#else
int dct_luma(int block_x,int block_y,int *coeff_cost, int intra)
#endif
{
  int sign(int a,int b);

  int i,j,ilev, m4[4][4], m5[4],m6[4],coeff_ctr;
  int ii;
  //int qp_const;
  int level,scan_pos,run;
  int nonzero;
  int qp_per,qp_rem,q_bits;

  int   pos_x   = block_x >> BLOCK_SHIFT;
  int   pos_y   = block_y >> BLOCK_SHIFT;
  int   b8      = 2*(pos_y >> 1) + (pos_x >> 1);
  int   b4      = 2*(pos_y & 0x01) + (pos_x & 0x01);
  int*  ACLevel = img->cofAC[b8][b4][0];
  int*  ACRun   = img->cofAC[b8][b4][1];
  int   pix_y, pix_x;

  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  short is_field_mode = (img->field_picture || ( img->MbaffFrameFlag && currMB->mb_field));

  Boolean lossless_qpprime = (Boolean)((currMB->qp + img->bitdepth_luma_qp_scale)==0 && img->lossless_qpprime_flag==1);
  int **levelscale,**leveloffset;
  int **invlevelscale;
#ifdef USE_INTRA_MDDT
   byte pos_scan[16][2];
#else
  const byte (*pos_scan)[2] = is_field_mode ? FIELD_SCAN : SNGL_SCAN;
#endif

#ifdef RDO_Q
  levelDataStruct levelData[16];
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
  double  lambda_md=0, normFact=pow(2,(2*DQ_BITS+19))*(1<<(2*img->BitDepthIncrease));
#else
  double  lambda_md=0, normFact=pow(2,(2*DQ_BITS+19));
#endif
  double err;
  int lowerInt, levelTrellis[16], kStart=0, kStop=0, noCoeff, estBits;
#endif
#ifdef USE_INTRA_MDDT
  if(input->UseIntraMDDT && intra && !img->residue_transform_flag && !is_field_mode)
  {
    for(i = 0; i < 16; i++)
    {
      pos_scan[i][0] = img->scanOrder4x4[ipmode][i][0];
      pos_scan[i][1] = img->scanOrder4x4[ipmode][i][1];
    }
  }
  else if(!is_field_mode)
  {
    for(i = 0; i < 16; i++)
    {
      pos_scan[i][0] = SNGL_SCAN[i][0];
      pos_scan[i][1] = SNGL_SCAN[i][1];
    }
  }
  else
  {
    for(i = 0; i < 16; i++)
    {
      pos_scan[i][0] = FIELD_SCAN[i][0];
      pos_scan[i][1] = FIELD_SCAN[i][1];
    }
  }
#endif

  qp_per    = (currMB->qp + img->bitdepth_luma_qp_scale - MIN_QP)/6; 
  qp_rem    = (currMB->qp + img->bitdepth_luma_qp_scale - MIN_QP)%6; 
  q_bits    = Q_BITS+qp_per;

#ifdef RDO_Q
  if(input->UseRDO_Q)
  {
    if ((img->type==B_SLICE) && img->nal_reference_idc)
    {
      lambda_md = img->lambda_md[5][img->masterQP];  }
    else
    {
      lambda_md = img->lambda_md[img->type][img->masterQP];
    }
  }
#endif

#ifdef ADAPTIVE_QUANTIZATION
  if(img->slice_fractional_quant_flag)
  {
    levelscale    = LevelScale4x4Luma_IAQMS[img->mb_iaqms_idx][intra][qp_rem];
    invlevelscale = InvLevelScale4x4Luma_IAQMS[img->mb_iaqms_idx][intra][qp_rem];
  }
  else
  {
    levelscale    = LevelScale4x4Luma[intra][qp_rem];
    invlevelscale = InvLevelScale4x4Luma[intra][qp_rem];
  }
  leveloffset   = LevelOffset4x4Luma[intra][qp_per];
#else
  levelscale    = LevelScale4x4Luma[intra][qp_rem];
  leveloffset   = LevelOffset4x4Luma[intra][qp_per];
  invlevelscale = InvLevelScale4x4Luma[intra][qp_rem];
#endif

  //  Horizontal transform
  if (!lossless_qpprime)
  {
    for (j=0; j < BLOCK_SIZE; j++)
    {
      m5[0] = img->m7[j][0]+img->m7[j][3];
      m5[1] = img->m7[j][1]+img->m7[j][2];
      m5[2] = img->m7[j][1]-img->m7[j][2];
      m5[3] = img->m7[j][0]-img->m7[j][3];
      
      m4[j][0] = m5[0]   + m5[1];
      m4[j][2] = m5[0]   - m5[1];
      m4[j][1] = m5[3]*2 + m5[2];
      m4[j][3] = m5[3]   - m5[2]*2;
    }
    
    //  Vertical transform
    for (i=0; i < BLOCK_SIZE; i++)
    {    
      m5[0] = m4[0][i] + m4[3][i];
      m5[1] = m4[1][i] + m4[2][i];
      m5[2] = m4[1][i] - m4[2][i];
      m5[3] = m4[0][i] - m4[3][i];
      
      m4[0][i] = m5[0]   + m5[1];
      m4[2][i] = m5[0]   - m5[1];
      m4[1][i] = m5[3]*2 + m5[2];
      m4[3][i] = m5[3]   - m5[2]*2;
    }
  }
  // Quant

  nonzero=FALSE;

  run=-1;
  scan_pos=0;
#ifdef RDO_Q
//#ifdef TREL_CAVLC
  if(input->UseRDO_Q && active_pps->entropy_coding_mode_flag == UVLC)
    TrellisCAVLC4x4(m4, q_bits, qp_rem, levelscale, leveloffset, levelTrellis, LUMA, b8, b4, 16, lambda_md);
//#endif
  if(input->UseRDO_Q && active_pps->entropy_coding_mode_flag == CABAC)
  {
    noCoeff=0;
    for (coeff_ctr=0;coeff_ctr < 16;coeff_ctr++)
    {
      i=pos_scan[coeff_ctr][0];
      j=pos_scan[coeff_ctr][1];

      levelData[coeff_ctr].levelDouble = absm(m4[j][i] * levelscale[i][j]);
      level = (int)(levelData[coeff_ctr].levelDouble >> q_bits);
      lowerInt=(((int)levelData[coeff_ctr].levelDouble - (level << q_bits)) < (1 <<( q_bits - 1)))? 1: 0;

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
        levelData[coeff_ctr].errLevel[ii]=(err*err*(double)estErr4x4[qp_rem][i][j])/normFact; 
      }
    }

    estBits=est_write_and_store_CBP_block_bit(currMB, LUMA_4x4);
    est_writeRunLevel_CABAC(levelData, levelTrellis, LUMA_4x4, lambda_md, kStart, kStop, noCoeff, estBits);
  }
#endif

  for (coeff_ctr=0;coeff_ctr < 16;coeff_ctr++)
  {

    i=pos_scan[coeff_ctr][0];
    j=pos_scan[coeff_ctr][1];
    
    run++;
    ilev=0;

#ifdef RDO_Q
  if(input->UseRDO_Q)    
    level = levelTrellis[coeff_ctr];
  else
  {
    if(lossless_qpprime)
      level = absm (img->m7[j][i]);
    else
      level = (absm (m4[j][i]) * levelscale[i][j] + leveloffset[i][j]) >> q_bits;
  }
#else
    if(lossless_qpprime)
      level = absm (img->m7[j][i]);
    else
      level = (absm (m4[j][i]) * levelscale[i][j] + leveloffset[i][j]) >> q_bits;
#endif

    if (img->AdaptiveRounding)
    {
      if (lossless_qpprime || level == 0 )
      {
        img->fadjust4x4[intra][block_y+j][block_x+i] = 0;
      }
      else 
      {
        img->fadjust4x4[intra][block_y+j][block_x+i] = 
          (AdaptRndWeight * (absm(m4[j][i]) * levelscale[i][j] - (level << q_bits)) + (1<< (q_bits))) >> (q_bits + 1);         
      }
    }

    if (level != 0)
    {
      nonzero=TRUE;

      *coeff_cost += (level > 1 || lossless_qpprime) ? MAX_VALUE : COEFF_COST[input->disthres][run];

      if(lossless_qpprime)
        ACLevel[scan_pos] = sign(level,img->m7[j][i]);
      else
        ACLevel[scan_pos] = sign(level,m4[j][i]);

      ACRun  [scan_pos] = run;
      ++scan_pos;
      run=-1;                     // reset zero level counter

      level=sign(level, m4[j][i]);

      if(lossless_qpprime)
      {
        ilev=level;
      }
#if 0
      else if(qp_per<4)
      {
        ilev=(level*invlevelscale[i][j]+(1<<(3-qp_per)))>>(4-qp_per);
      }
      else
      {
        ilev=(level*invlevelscale[i][j])<<(qp_per-4);
      }
#else
      else
      {
        ilev=((((level*invlevelscale[i][j])<< qp_per) + 8 ) >> 4);
      }
#endif
    }
    if(!lossless_qpprime)
      m4[j][i]=ilev;
  }

  ACLevel[scan_pos] = 0;  
  
  //     IDCT.
  //     horizontal

  if (!lossless_qpprime)
  {
    for (j=0; j < BLOCK_SIZE; j++)
    {
      m6[0]=(m4[j][0]     +  m4[j][2]);
      m6[1]=(m4[j][0]     -  m4[j][2]);
      m6[2]=(m4[j][1]>>1) -  m4[j][3];
      m6[3]= m4[j][1]     + (m4[j][3]>>1);
      
      m4[j][0] = m6[0] + m6[3];
      m4[j][1] = m6[1] + m6[2];
      m4[j][2] = m6[1] - m6[2];
      m4[j][3] = m6[0] - m6[3];
    }
    
    //  vertical
    for (i=0; i < BLOCK_SIZE; i++)
    {
      
      m6[0]=(m4[0][i]     +  m4[2][i]);
      m6[1]=(m4[0][i]     -  m4[2][i]);
      m6[2]=(m4[1][i]>>1) -  m4[3][i];
      m6[3]= m4[1][i]     + (m4[3][i]>>1);
      
      ii = i + block_x;
      
      if (!img->residue_transform_flag)
      {
        img->m7[0][i] = min(img->max_imgpel_value,max(0,(m6[0]+m6[3]+((long)img->mpr[0 + block_y][ii] << DQ_BITS)+DQ_ROUND)>>DQ_BITS));
        img->m7[1][i] = min(img->max_imgpel_value,max(0,(m6[1]+m6[2]+((long)img->mpr[1 + block_y][ii] << DQ_BITS)+DQ_ROUND)>>DQ_BITS));
        img->m7[2][i] = min(img->max_imgpel_value,max(0,(m6[1]-m6[2]+((long)img->mpr[2 + block_y][ii] << DQ_BITS)+DQ_ROUND)>>DQ_BITS));
        img->m7[3][i] = min(img->max_imgpel_value,max(0,(m6[0]-m6[3]+((long)img->mpr[3 + block_y][ii] << DQ_BITS)+DQ_ROUND)>>DQ_BITS));
      } 
      else 
      {
        if(lossless_qpprime)
        {
          img->m7[0][i] = m6[0]+m6[3];
          img->m7[1][i] = m6[1]+m6[2];
          img->m7[2][i] = m6[1]-m6[2];
          img->m7[3][i] = m6[0]-m6[3];
        }
        else
        {
          img->m7[0][i] =(m6[0]+m6[3]+DQ_ROUND)>>DQ_BITS;
          img->m7[1][i] =(m6[1]+m6[2]+DQ_ROUND)>>DQ_BITS;
          img->m7[2][i] =(m6[1]-m6[2]+DQ_ROUND)>>DQ_BITS;
          img->m7[3][i] =(m6[0]-m6[3]+DQ_ROUND)>>DQ_BITS;
        }
      }
    }
  }
  //  Decoded block moved to frame memory
  if (!img->residue_transform_flag)
  {
    if(lossless_qpprime)
    {
      for (j=0; j < BLOCK_SIZE; j++)
      {
        pix_y = img->pix_y+block_y+j;
        for (i=0; i < BLOCK_SIZE; i++)
        {
         enc_picture->imgY[pix_y][img->pix_x+block_x+i]=img->m7[j][i]+img->mpr[j+block_y][i+block_x];
        }
      }
    }
    else
    {
      for (j=0; j < BLOCK_SIZE; j++)
      {
        pix_y = img->pix_y+block_y+j;
        pix_x = img->pix_x+block_x;
        for (i=0; i < BLOCK_SIZE; i++)
        {
          enc_picture->imgY[pix_y][pix_x + i]=img->m7[j][i];
        }
      }
    }
    
  }
  return nonzero;
}


/*!
 ************************************************************************
 * \brief
 *    Transform,quantization,inverse transform for chroma.
 *    The main reason why this is done in a separate routine is the
 *    additional 2x2 transform of DC-coeffs. This routine is called
 *    ones for each of the chroma components.
 *
 * \par Input:
 *    uv    : Make difference between the U and V chroma component  \n
 *    cr_cbp: chroma coded block pattern
 *
 * \par Output:
 *    cr_cbp: Updated chroma coded block pattern.
 ************************************************************************
 */
int dct_chroma(int uv,int cr_cbp)
{
  int i,j,i1,j2,ilev,n2,n1,j1,mb_y,coeff_ctr,level ,scan_pos,run;
  int m1[BLOCK_SIZE],m5[BLOCK_SIZE],m6[BLOCK_SIZE];
  int coeff_cost;
  int cr_cbp_tmp;
  int DCcoded=0 ;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
 
  int qp_per,qp_rem,q_bits;

  int   b4;
  int*  DCLevel = img->cofDC[uv+1][0];
  int*  DCRun   = img->cofDC[uv+1][1];
  int*  ACLevel;
  int*  ACRun;
  int   intra = IS_INTRA (currMB);
  int   uv_scale = uv*(img->num_blk8x8_uv >> 1);

  //FRExt
  int64 cbpblk_pattern[4]={0, 0xf0000, 0xff0000, 0xffff0000};
  int yuv = img->yuv_format;
  int b8;
  int m3[4][4];
  int m4[4][4];
  int qp_per_dc = 0;
  int qp_rem_dc = 0;
  int q_bits_422 = 0;  
  int ***levelscale, ***leveloffset;
  int ***invlevelscale;
  short pix_c_x, pix_c_y;
  short is_field_mode = (img->field_picture || ( img->MbaffFrameFlag && currMB->mb_field));
  const byte (*pos_scan)[2] = is_field_mode ? FIELD_SCAN : SNGL_SCAN;

  Boolean lossless_qpprime = (Boolean)((currMB->qp + img->bitdepth_luma_qp_scale)==0 && img->lossless_qpprime_flag==1);
#ifdef RDO_Q
  levelDataStruct levelData[16];
  double  lambda_md=0;
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
  double normFact = pow(2, (2 * DQ_BITS + 19))*(1<<(2*img->BitDepthIncrease));
#else
  double normFact = pow(2, (2 * DQ_BITS + 19));
#endif
  double err;
  int lowerInt, levelTrellis[16], k, kStart, kStop, noCoeff, estBits;
#endif

  qp_per    = qp_per_matrix[(currMB->qpc[uv] + img->bitdepth_chroma_qp_scale)];
  qp_rem    = qp_rem_matrix[(currMB->qpc[uv] + img->bitdepth_chroma_qp_scale)];

  q_bits    = Q_BITS+qp_per;

#ifdef RDO_Q
  //if (input->successive_Bframe>0)
  if(input->UseRDO_Q)
  {
    if ((img->type==B_SLICE) && img->nal_reference_idc)
    {
      lambda_md = img->lambda_md[5][img->masterQP];  
    }
    else
    {
      if (input->successive_Bframe==0 && img->type != B_SLICE)
      {
        int masterQP_C;
        masterQP_C = Clip3 ( -img->bitdepth_chroma_qp_scale, 51, img->masterQP + img->chroma_qp_offset[uv] );
        masterQP_C = (masterQP_C < 0 ? masterQP_C : QP_SCALE_CR[masterQP_C]) + img->bitdepth_chroma_qp_scale;

        if(masterQP_C < img->masterQP)
        {
          lambda_md = img->lambda_md[img->type][masterQP_C];
        }
        else
        {
          lambda_md = img->lambda_md[img->type][img->masterQP];
        }
      }
      else
      {
        lambda_md = img->lambda_md[img->type][img->masterQP];
      }
    }
  }
#endif

#ifdef ADAPTIVE_QUANTIZATION
  if(img->slice_fractional_quant_flag)
  {
    levelscale    = LevelScale4x4Chroma_IAQMS   [img->mb_iaqms_idx][uv][intra];
    invlevelscale = InvLevelScale4x4Chroma_IAQMS[img->mb_iaqms_idx][uv][intra];
  }
  else
  {
    levelscale = LevelScale4x4Chroma[uv][intra];
    invlevelscale = InvLevelScale4x4Chroma[uv][intra];
  }
  leveloffset = LevelOffset4x4Chroma[uv][intra];
#else
  levelscale = LevelScale4x4Chroma[uv][intra];
  leveloffset = LevelOffset4x4Chroma[uv][intra];
  invlevelscale = InvLevelScale4x4Chroma[uv][intra];
#endif

  if (img->yuv_format == YUV422)
  {
    //for YUV422 only
    qp_per_dc = (currMB->qpc[uv] + 3 + img->bitdepth_chroma_qp_scale)/6;
    qp_rem_dc = (currMB->qpc[uv] + 3 + img->bitdepth_chroma_qp_scale)%6;
    
    q_bits_422 = Q_BITS+qp_per_dc;  
  }

  
  //============= dct transform ===============  
  for (n2=0; n2 < img->mb_cr_size_y; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 < img->mb_cr_size_x; n1 += BLOCK_SIZE)
    {

      //  Horizontal transform.
      for (j=0; j < BLOCK_SIZE && !lossless_qpprime; j++)
      {
        mb_y=n2+j;
        
        m5[0]=img->m7[mb_y][n1  ]+img->m7[mb_y][n1+3];
        m5[1]=img->m7[mb_y][n1+1]+img->m7[mb_y][n1+2];
        m5[2]=img->m7[mb_y][n1+1]-img->m7[mb_y][n1+2];
        m5[3]=img->m7[mb_y][n1  ]-img->m7[mb_y][n1+3];
        
        img->m7[mb_y][n1  ] = (m5[0]   + m5[1]);
        img->m7[mb_y][n1+2] = (m5[0]   - m5[1]);
        img->m7[mb_y][n1+1] =  m5[3]*2 + m5[2];
        img->m7[mb_y][n1+3] =  m5[3]   - m5[2]*2;
      }

      //  Vertical transform.

      for (i=0; i < BLOCK_SIZE && !lossless_qpprime; i++)
      {
        j1=n1+i;
        m5[0] = img->m7[n2  ][j1] + img->m7[n2+3][j1];
        m5[1] = img->m7[n2+1][j1] + img->m7[n2+2][j1];
        m5[2] = img->m7[n2+1][j1] - img->m7[n2+2][j1];
        m5[3] = img->m7[n2  ][j1] - img->m7[n2+3][j1];

        img->m7[n2+0][j1] = (m5[0]   + m5[1]);
        img->m7[n2+2][j1] = (m5[0]   - m5[1]);
        img->m7[n2+1][j1] =  m5[3]*2 + m5[2];
        img->m7[n2+3][j1] =  m5[3]   - m5[2]*2;
      }
    }
  }
  
  if (yuv == YUV420)
  {
    //================== CHROMA DC YUV420 ===================
    //     2X2 transform of DC coeffs.
    if(lossless_qpprime)
    {
      m1[0]=img->m7[0][0];
      m1[1]=img->m7[0][4];
      m1[2]=img->m7[4][0];
      m1[3]=img->m7[4][4];
    }
    else 
    {
      m1[0]=(img->m7[0][0] + img->m7[0][4] + img->m7[4][0] + img->m7[4][4]);
      m1[1]=(img->m7[0][0] - img->m7[0][4] + img->m7[4][0] - img->m7[4][4]);
      m1[2]=(img->m7[0][0] + img->m7[0][4] - img->m7[4][0] - img->m7[4][4]);
      m1[3]=(img->m7[0][0] - img->m7[0][4] - img->m7[4][0] + img->m7[4][4]);
    }
    
    //     Quant of chroma 2X2 coeffs.
    run=-1;
    scan_pos=0;
    
    for (coeff_ctr=0; coeff_ctr < 4; coeff_ctr++)
    {
      run++;
      ilev=0;
      
      if(lossless_qpprime)
        level =absm(m1[coeff_ctr]);
      else 
        level =(absm(m1[coeff_ctr]) * levelscale[qp_rem][0][0] + (leveloffset[qp_per][0][0]<<1)) >> (q_bits+1);
      
      if (input->symbol_mode == UVLC && img->qp < 4) 
      {
        if (level > CAVLC_LEVEL_LIMIT) 
          level = CAVLC_LEVEL_LIMIT;
      }
      
      if (level  != 0)
      {
        currMB->cbp_blk |= 0xf0000 << (uv << 2) ;    // if one of the 2x2-DC levels is != 0 set the
        cr_cbp=max(1,cr_cbp);                     // coded-bit all 4 4x4 blocks (bit 16-19 or 20-23)
        DCcoded = 1 ;
        DCLevel[scan_pos] = sign(level ,m1[coeff_ctr]);
        DCRun  [scan_pos] = run;
        scan_pos++;
        run=-1;
        
        ilev=sign(level, m1[coeff_ctr]);
      }
      if(!lossless_qpprime)
        m1[coeff_ctr]=ilev;
    }
    DCLevel[scan_pos] = 0;
    
    //  Inverse transform of 2x2 DC levels
    if(!lossless_qpprime)
    {
      m5[0]=(m1[0] + m1[1] + m1[2] + m1[3]);
      m5[1]=(m1[0] - m1[1] + m1[2] - m1[3]);
      m5[2]=(m1[0] + m1[1] - m1[2] - m1[3]);
      m5[3]=(m1[0] - m1[1] - m1[2] + m1[3]);
      if(qp_per<5)
      {
        for(i=0; i<4; i++)
          m1[i]=(m5[i] * invlevelscale[qp_rem][0][0])>>(5-qp_per);
      }
      else
      {
        for(i=0; i<4; i++)
          m1[i]=(m5[i] * invlevelscale[qp_rem][0][0])<<(qp_per-5);
      }

      img->m7[0][0] = m1[0];
      img->m7[0][4] = m1[1];
      img->m7[4][0] = m1[2];
      img->m7[4][4] = m1[3];
    }
  }
  else if(yuv == YUV422)
  {
    //================== CHROMA DC YUV422 ===================
    //transform DC coeff
    //horizontal
    
    //pick out DC coeff
    for (j=0; j < img->mb_cr_size_y; j+=BLOCK_SIZE)
    {
      for (i=0; i < img->mb_cr_size_x; i+=BLOCK_SIZE)
        m3[i>>2][j>>2]= img->m7[j][i];
    } 
    //horizontal
    if(!lossless_qpprime)
    {
      m4[0][0] = m3[0][0] + m3[1][0];
      m4[0][1] = m3[0][1] + m3[1][1];
      m4[0][2] = m3[0][2] + m3[1][2];
      m4[0][3] = m3[0][3] + m3[1][3];
      
      m4[1][0] = m3[0][0] - m3[1][0];
      m4[1][1] = m3[0][1] - m3[1][1];
      m4[1][2] = m3[0][2] - m3[1][2];
      m4[1][3] = m3[0][3] - m3[1][3];
      
      // vertical
      for (i=0;i<2;i++)
      {
        m5[0] = m4[i][0] + m4[i][3];
        m5[1] = m4[i][1] + m4[i][2];
        m5[2] = m4[i][1] - m4[i][2];
        m5[3] = m4[i][0] - m4[i][3];
        
        m4[i][0] = (m5[0] + m5[1]);
        m4[i][2] = (m5[0] - m5[1]);
        m4[i][1] = (m5[3] + m5[2]);
        m4[i][3] = (m5[3] - m5[2]);
      }
    }
    
    run=-1;
    scan_pos=0;
    
    //quant of chroma DC-coeffs
    for (coeff_ctr=0;coeff_ctr<8;coeff_ctr++)
    {
      i=SCAN_YUV422[coeff_ctr][0];
      j=SCAN_YUV422[coeff_ctr][1];
      
      run++;

      if(lossless_qpprime)
      {
        level = absm(m3[i][j]);
        m4[i][j]=m3[i][j];
      }
      else 
        level =(absm(m4[i][j]) * levelscale[qp_rem_dc][0][0] + (leveloffset[qp_per_dc][0][0]*2)) >> (q_bits_422+1);

      if (level != 0)
      {
        //YUV422
        currMB->cbp_blk |= 0xff0000 << (uv << 3) ;   // if one of the DC levels is != 0 set the
        cr_cbp=max(1,cr_cbp);                           // coded-bit all 4 4x4 blocks (bit 16-31 or 32-47) //YUV444
        DCcoded = 1 ;
        
        DCLevel[scan_pos] = sign(level,m4[i][j]);
        DCRun  [scan_pos] = run;
        ++scan_pos;
        run=-1;
      }
      if(!lossless_qpprime)
        m3[i][j]=sign(level,m4[i][j]);
    }
    DCLevel[scan_pos]=0;

    //inverse DC transform
    //horizontal
    if(!lossless_qpprime)
    {
      m4[0][0] = m3[0][0] + m3[1][0];
      m4[0][1] = m3[0][1] + m3[1][1];
      m4[0][2] = m3[0][2] + m3[1][2];
      m4[0][3] = m3[0][3] + m3[1][3];
      
      m4[1][0] = m3[0][0] - m3[1][0];
      m4[1][1] = m3[0][1] - m3[1][1];
      m4[1][2] = m3[0][2] - m3[1][2];
      m4[1][3] = m3[0][3] - m3[1][3];      
      
      // vertical
      for (i=0;i<2;i++)
      {       
        m6[0]=m4[i][0]+m4[i][2];
        m6[1]=m4[i][0]-m4[i][2];
        m6[2]=m4[i][1]-m4[i][3];
        m6[3]=m4[i][1]+m4[i][3];
        
        if(qp_per_dc<4)
        {          
          img->m7[0 ][i*4]=((((m6[0]+m6[3])*invlevelscale[qp_rem_dc][0][0]+(1<<(3-qp_per_dc)))>>(4-qp_per_dc))+2)>>2;
          img->m7[4 ][i*4]=((((m6[1]+m6[2])*invlevelscale[qp_rem_dc][0][0]+(1<<(3-qp_per_dc)))>>(4-qp_per_dc))+2)>>2;
          img->m7[8 ][i*4]=((((m6[1]-m6[2])*invlevelscale[qp_rem_dc][0][0]+(1<<(3-qp_per_dc)))>>(4-qp_per_dc))+2)>>2;
          img->m7[12][i*4]=((((m6[0]-m6[3])*invlevelscale[qp_rem_dc][0][0]+(1<<(3-qp_per_dc)))>>(4-qp_per_dc))+2)>>2;
        }
        else
        {
          img->m7[0 ][i*4]=((((m6[0]+m6[3])*invlevelscale[qp_rem_dc][0][0])<<(qp_per_dc-4))+2)>>2;
          img->m7[4 ][i*4]=((((m6[1]+m6[2])*invlevelscale[qp_rem_dc][0][0])<<(qp_per_dc-4))+2)>>2;
          img->m7[8 ][i*4]=((((m6[1]-m6[2])*invlevelscale[qp_rem_dc][0][0])<<(qp_per_dc-4))+2)>>2;
          img->m7[12][i*4]=((((m6[0]-m6[3])*invlevelscale[qp_rem_dc][0][0])<<(qp_per_dc-4))+2)>>2;
        }
      }//for (i=0;i<2;i++)    
    }
  }
  else if(yuv == YUV444)
  {
    //================== CHROMA DC YUV444 ===================
    //transform DC coeff
    //pick out DC coeff
    for (j=0; j < img->mb_cr_size_y; j+=BLOCK_SIZE)
    {
      for (i=0; i < img->mb_cr_size_x; i+=BLOCK_SIZE)
        m4[i>>2][j>>2]= img->m7[j][i];
    }
    
    //horizontal
    for (j=0;j<4 && !lossless_qpprime;j++)
    {
      m5[0] = m4[0][j] + m4[3][j];
      m5[1] = m4[1][j] + m4[2][j];
      m5[2] = m4[1][j] - m4[2][j];
      m5[3] = m4[0][j] - m4[3][j];
      
      m4[0][j]=m5[0]+m5[1];
      m4[2][j]=m5[0]-m5[1];
      m4[1][j]=m5[3]+m5[2];
      m4[3][j]=m5[3]-m5[2];
    }
    // vertical
    for (i=0;i<4 && !lossless_qpprime;i++)
    {
      m5[0] = m4[i][0] + m4[i][3];
      m5[1] = m4[i][1] + m4[i][2];
      m5[2] = m4[i][1] - m4[i][2];
      m5[3] = m4[i][0] - m4[i][3];

      m4[i][0]=(m5[0]+m5[1])>>1;
      m4[i][2]=(m5[0]-m5[1])>>1;
      m4[i][1]=(m5[3]+m5[2])>>1;
      m4[i][3]=(m5[3]-m5[2])>>1;
    }

    run=-1;
    scan_pos=0;
    
    //quant of chroma DC-coeffs
    for (coeff_ctr=0;coeff_ctr<16;coeff_ctr++)
    {
      i=SNGL_SCAN[coeff_ctr][0];
      j=SNGL_SCAN[coeff_ctr][1];
      
      run++;
      
      if(lossless_qpprime)
        level = absm(m4[i][j]);
      else 
        level =(absm(m4[i][j]) * levelscale[qp_rem][0][0] + (leveloffset[qp_per][0][0]*2)) >> (q_bits+1);
      
      if (level != 0)
      {
        //YUV444
        currMB->cbp_blk |= ((int64)0xffff0000) << (uv << 4) ;   // if one of the DC levels is != 0 set the
        cr_cbp=max(1,cr_cbp);                           // coded-bit all 4 4x4 blocks (bit 16-31 or 32-47) //YUV444
        DCcoded = 1 ;
        
        DCLevel[scan_pos] = sign(level,m4[i][j]);
        DCRun  [scan_pos] = run;
        ++scan_pos;
        run=-1;
      }
      if(!lossless_qpprime)
        m4[i][j]=sign(level,m4[i][j]);
    }
    DCLevel[scan_pos]=0;

    // inverse DC transform
    //horizontal
    for (j=0;j<4 && !lossless_qpprime;j++)
    {     
      m6[0] = m4[0][j] + m4[2][j];
      m6[1] = m4[0][j] - m4[2][j];
      m6[2] = m4[1][j] - m4[3][j];
      m6[3] = m4[1][j] + m4[3][j];
      
      m4[0][j] = m6[0] + m6[3];
      m4[1][j] = m6[1] + m6[2];
      m4[2][j] = m6[1] - m6[2];
      m4[3][j] = m6[0] - m6[3];
    }
    
    //vertical
    for (i=0;i<4 && !lossless_qpprime;i++)
    {
      m6[0]=m4[i][0]+m4[i][2];
      m6[1]=m4[i][0]-m4[i][2];
      m6[2]=m4[i][1]-m4[i][3];
      m6[3]=m4[i][1]+m4[i][3];

      if(qp_per<4)
      {
        img->m7[0 ][i*4] = ((((m6[0] + m6[3])*invlevelscale[qp_rem][0][0]+(1<<(3-qp_per)))>>(4-qp_per))+2)>>2;
        img->m7[4 ][i*4] = ((((m6[1] + m6[2])*invlevelscale[qp_rem][0][0]+(1<<(3-qp_per)))>>(4-qp_per))+2)>>2;
        img->m7[8 ][i*4] = ((((m6[1] - m6[2])*invlevelscale[qp_rem][0][0]+(1<<(3-qp_per)))>>(4-qp_per))+2)>>2;
        img->m7[12][i*4] = ((((m6[0] - m6[3])*invlevelscale[qp_rem][0][0]+(1<<(3-qp_per)))>>(4-qp_per))+2)>>2;
      }
      else
      {
        img->m7[0 ][i*4] = ((((m6[0]+m6[3])*invlevelscale[qp_rem][0][0])<<(qp_per-4))+2)>>2;
        img->m7[4 ][i*4] = ((((m6[1]+m6[2])*invlevelscale[qp_rem][0][0])<<(qp_per-4))+2)>>2;
        img->m7[8 ][i*4] = ((((m6[1]-m6[2])*invlevelscale[qp_rem][0][0])<<(qp_per-4))+2)>>2;
        img->m7[12][i*4] = ((((m6[0]-m6[3])*invlevelscale[qp_rem][0][0])<<(qp_per-4))+2)>>2;
      }
    }
  }

  //     Quant of chroma AC-coeffs.
  coeff_cost=0;
  cr_cbp_tmp=0;

  for (b8=0; b8 < (img->num_blk8x8_uv >> 1); b8++)
  {
    for (b4=0; b4 < 4; b4++)
    {
      n1 = hor_offset[yuv][b8][b4];
      n2 = ver_offset[yuv][b8][b4];
      ACLevel = img->cofAC[4+b8+uv_scale][b4][0];
      ACRun   = img->cofAC[4+b8+uv_scale][b4][1];
      run=-1;
      scan_pos=0;



#ifdef RDO_Q
//#ifdef TREL_CAVLC
      if(input->UseRDO_Q && active_pps->entropy_coding_mode_flag == UVLC)
      {
        int M4[4][4], n, m;
        for(n=0; n<4; n++)
          for(m=0; m<4; m++)
            M4[n][m] = img->m7[n2+n][n1+m];

        TrellisCAVLC4x4(M4, q_bits, qp_rem, levelscale[qp_rem], leveloffset[qp_per], levelTrellis, CHROMA_AC, b8, b4, 15, lambda_md);
      }
//#endif
      if(input->UseRDO_Q && active_pps->entropy_coding_mode_flag == CABAC)
      {
        kStart=0; kStop=0; noCoeff=0;
        for (coeff_ctr=0;coeff_ctr < 15;coeff_ctr++)
        {
          i=pos_scan[coeff_ctr+1][0]; // scan is shifted due to DC
          j=pos_scan[coeff_ctr+1][1]; // scan is shifted due to DC

          levelData[coeff_ctr].levelDouble=absm(img->m7[n2+j][n1+i]*levelscale[qp_rem][i][j]);
          level = (int)(levelData[coeff_ctr].levelDouble >> q_bits);
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

          for (k=0; k<levelData[coeff_ctr].noLevels; k++)
          {
            err=(double)(levelData[coeff_ctr].level[k]<<q_bits)-(double)levelData[coeff_ctr].levelDouble;
            levelData[coeff_ctr].errLevel[k]=(err*err*(double)estErr4x4[qp_rem][i][j])/normFact; 
          }
        }
        estBits=est_write_and_store_CBP_block_bit(currMB, LUMA_16AC);

        est_writeRunLevel_CABAC(levelData, levelTrellis, LUMA_16AC, lambda_md, kStart, kStop, noCoeff, estBits);
      }
#endif


      for (coeff_ctr=1; coeff_ctr < 16; coeff_ctr++)// start change rd_quant
      {

        i=pos_scan[coeff_ctr][0];
        j=pos_scan[coeff_ctr][1];

        ++run;
        ilev=0;

#ifdef RDO_Q
        if(input->UseRDO_Q)
        {
          level=levelTrellis[coeff_ctr-1];
        }
        else
        {
          if(lossless_qpprime)
            level = absm(img->m7[n2+j][n1+i]);
          else 
            level=(absm(img->m7[n2+j][n1+i])*levelscale[qp_rem][i][j]+leveloffset[qp_per][i][j])>>q_bits;
        }
#else
        if(lossless_qpprime)
          level = absm(img->m7[n2+j][n1+i]);
        else 
          level=(absm(img->m7[n2+j][n1+i])*levelscale[qp_rem][i][j]+leveloffset[qp_per][i][j])>>q_bits;
#endif

        if (img->AdaptiveRounding)
        {
          if (lossless_qpprime || level == 0 )
          {
            img->fadjust4x4Cr[intra][uv][n2+j][n1+i] = 0;
          }
          else
          {
            img->fadjust4x4Cr[intra][uv][n2+j][n1+i] = 
              (AdaptRndWeight * (absm(img->m7[n2+j][n1+i]) * levelscale[qp_rem][i][j] - (level << q_bits)) + (1<< (q_bits))) >> (q_bits + 1); 
          }
        }

        if (level  != 0)
        {
          currMB->cbp_blk |= ((int64)1) << cbp_blk_chroma[b8 + uv_scale][b4];
          if (level > 1 || lossless_qpprime)
            coeff_cost += MAX_VALUE;                // set high cost, shall not be discarded
          else
            coeff_cost += COEFF_COST[input->disthres][run];

          cr_cbp_tmp=2;
          ACLevel[scan_pos] = sign(level,img->m7[n2+j][n1+i]);
          ACRun  [scan_pos] = run;
          ++scan_pos;
          run=-1;

          level=sign(level, img->m7[n2+j][n1+i]);
          if(lossless_qpprime)
          {
            ilev = level;
          }
          else if(qp_per<4)
          {
            ilev=(level*invlevelscale[qp_rem][i][j]+(1<<(3-qp_per)))>>(4-qp_per);
          }
          else
          {
            ilev=(level*invlevelscale[qp_rem][i][j])<<(qp_per-4);
          }
        }
        if(!lossless_qpprime)
          img->m7[n2+j][n1+i]=ilev;
      }
      ACLevel[scan_pos] = 0;
    }
  }

  // * reset chroma coeffs
  if(coeff_cost < _CHROMA_COEFF_COST_ && !lossless_qpprime)
  {
    cr_cbp_tmp = 0 ;
    
    for (b8=0; b8 < (img->num_blk8x8_uv >> 1); b8++)
    {
      for (b4=0; b4 < 4; b4++)
      {
        n1 = hor_offset[yuv][b8][b4];
        n2 = ver_offset[yuv][b8][b4];
        ACLevel = img->cofAC[4+b8+uv_scale][b4][0];
        ACRun   = img->cofAC[4+b8+uv_scale][b4][1];
        if( DCcoded == 0) 
          currMB->cbp_blk &= ~((int64)cbpblk_pattern[yuv] << (uv << (1+yuv)));  // if no chroma DC's: then reset coded-bits of this chroma subblock
        
        ACLevel[0] = 0;
        for (coeff_ctr=1; coeff_ctr < 16; coeff_ctr++)// ac coeff
        {
          i=pos_scan[coeff_ctr][0];
          j=pos_scan[coeff_ctr][1];
          
          img->m7[n2+j][n1+i]=0;
          ACLevel[coeff_ctr] = 0;
        }
      }
    }
  }

  if(cr_cbp_tmp==2)   
    cr_cbp = 2;
  
  //     IDCT.
  //     Horizontal.
  for (n2=0; n2 < img->mb_cr_size_y && !lossless_qpprime; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 < img->mb_cr_size_x; n1 += BLOCK_SIZE)
    {
      for (j=0; j < BLOCK_SIZE; j++)
      {
        j2 = n2 + j;
        for (i=0; i < BLOCK_SIZE; i++)
        {
          m5[i]=img->m7[j2][n1+i];
        }

        m6[0] = (m5[0]     +  m5[2]);
        m6[1] = (m5[0]     -  m5[2]);
        m6[2] = (m5[1]>>1) -  m5[3];
        m6[3] =  m5[1]     + (m5[3]>>1);

        img->m7[j2][n1  ] = m6[0] + m6[3];
        img->m7[j2][n1+1] = m6[1] + m6[2];
        img->m7[j2][n1+2] = m6[1] - m6[2];
        img->m7[j2][n1+3] = m6[0] - m6[3];
      }

      //     Vertical.
      for (i=0; i < BLOCK_SIZE && !lossless_qpprime; i++)
      {
        i1 = n1 + i;
        for (j=0; j < BLOCK_SIZE; j++)
        {
          m5[j]=img->m7[n2+j][i1];
        }
        m6[0]=(m5[0]+m5[2]);
        m6[1]=(m5[0]-m5[2]);
        m6[2]=(m5[1]>>1)-m5[3];
        m6[3]=m5[1]+(m5[3]>>1);

          // Residue Color Transform
        if (!img->residue_transform_flag)
        {
          img->m7[n2  ][i1] = min(img->max_imgpel_value_uv,max(0,(m6[0]+m6[3]+((long)img->mpr[n2  ][i1] << DQ_BITS)+DQ_ROUND)>>DQ_BITS));
          img->m7[n2+1][i1] = min(img->max_imgpel_value_uv,max(0,(m6[1]+m6[2]+((long)img->mpr[n2+1][i1] << DQ_BITS)+DQ_ROUND)>>DQ_BITS));
          img->m7[n2+2][i1] = min(img->max_imgpel_value_uv,max(0,(m6[1]-m6[2]+((long)img->mpr[n2+2][i1] << DQ_BITS)+DQ_ROUND)>>DQ_BITS));
          img->m7[n2+3][i1] = min(img->max_imgpel_value_uv,max(0,(m6[0]-m6[3]+((long)img->mpr[n2+3][i1] << DQ_BITS)+DQ_ROUND)>>DQ_BITS));
        } 
        else 
        {
          if(lossless_qpprime)
          {
            img->m7[n2  ][i1] = m6[0]+m6[3];
            img->m7[n2+1][i1] = m6[1]+m6[2];
            img->m7[n2+2][i1] = m6[1]-m6[2];
            img->m7[n2+3][i1] = m6[0]-m6[3];
          }
          else
          {
            img->m7[n2  ][i1] = (m6[0]+m6[3]+DQ_ROUND)>>DQ_BITS;
            img->m7[n2+1][i1] = (m6[1]+m6[2]+DQ_ROUND)>>DQ_BITS;
            img->m7[n2+2][i1] = (m6[1]-m6[2]+DQ_ROUND)>>DQ_BITS;
            img->m7[n2+3][i1] = (m6[0]-m6[3]+DQ_ROUND)>>DQ_BITS;
          }
        }
      }
    }
  }

  //  Decoded block moved to memory
  if (!img->residue_transform_flag)
  {
    for (j=0; j < img->mb_cr_size_y; j++)
    {
      pix_c_y = img->pix_c_y+j;
      for (i=0; i < img->mb_cr_size_x; i++)
      {
        pix_c_x = img->pix_c_x+i;
        if(lossless_qpprime)
          enc_picture->imgUV[uv][pix_c_y][pix_c_x]= img->m7[j][i]+img->mpr[j][i];
        else
          enc_picture->imgUV[uv][pix_c_y][pix_c_x]= img->m7[j][i];
      }
    }
  }
  return cr_cbp;
}


// Residue Color Transform
int dct_chroma4x4(int uv, int b8, int b4)
{
  int sign(int a,int b);

  int i,j,i1,j1,ilev,m5[4],m6[4],coeff_ctr;
  int level,scan_pos,run;
  int nonzeroAC;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int   intra = IS_INTRA (currMB);

  int qp_per,qp_rem,q_bits;

  int*  ACLevel = img->cofAC[b8][b4][0];
  int*  ACRun   = img->cofAC[b8][b4][1];

  int **levelscale, **leveloffset;
  int **invlevelscale;

  Boolean lossless_qpprime = (Boolean)((img->qp + img->bitdepth_luma_qp_scale)==0 && img->lossless_qpprime_flag==1);

  qp_per    = qp_per_matrix[(currMB->qpc[uv] + img->bitdepth_chroma_qp_scale)];
  qp_rem    = qp_rem_matrix[(currMB->qpc[uv] + img->bitdepth_chroma_qp_scale)];
  q_bits    = Q_BITS+qp_per;

#ifdef ADAPTIVE_QUANTIZATION
  if(img->slice_fractional_quant_flag)
  {
    levelscale    = LevelScale4x4Chroma_IAQMS   [img->mb_iaqms_idx][uv][intra][qp_rem];
    invlevelscale = InvLevelScale4x4Chroma_IAQMS[img->mb_iaqms_idx][uv][intra][qp_rem];
  }
  else
  {
    levelscale = LevelScale4x4Chroma[uv][intra][qp_rem];
    invlevelscale = InvLevelScale4x4Chroma[uv][intra][qp_rem];
  }
  leveloffset = LevelOffset4x4Chroma[uv][intra][qp_per];
#else
  levelscale = LevelScale4x4Chroma[uv][intra][qp_rem];
  leveloffset = LevelOffset4x4Chroma[uv][intra][qp_per];
  invlevelscale = InvLevelScale4x4Chroma[uv][intra][qp_rem];
#endif

  //  Horizontal transform
  if(!lossless_qpprime)
  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < 2; i++)
    {
      i1=3-i;
      m5[i]=img->m7[j][i]+img->m7[j][i1];
      m5[i1]=img->m7[j][i]-img->m7[j][i1];
    }
    img->m7[j][0]=(m5[0]+m5[1]);
    img->m7[j][2]=(m5[0]-m5[1]);
    img->m7[j][1]=m5[3]*2+m5[2];
    img->m7[j][3]=m5[3]-m5[2]*2;
  }

  //  Vertical transform
  if(!lossless_qpprime)
  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < 2; j++)
    {
      j1=3-j;
      m5[j]=img->m7[j][i]+img->m7[j1][i];
      m5[j1]=img->m7[j][i]-img->m7[j1][i];
    }
    img->m7[0][i]=(m5[0]+m5[1]);
    img->m7[2][i]=(m5[0]-m5[1]);
    img->m7[1][i]=m5[3]*2+m5[2];
    img->m7[3][i]=m5[3]-m5[2]*2;
  }

  // Quant

  nonzeroAC=FALSE;

  run=-1;
  scan_pos=0;

  if(lossless_qpprime)
    level = absm(img->m7[0][0]);
  else 
    level =(absm(img->m7[0][0]) * levelscale[0][0] + leveloffset[0][0]) >> q_bits;

  b8 -= 4*(uv+1);
  dc_level_temp[uv][2*(b8 & 0x01)+(b4 & 0x01)][2*(b8 >> 1)+(b4 >> 1)] = sign(level, img->m7[0][0]);

  /* Inverse Quantization */
  if(lossless_qpprime)
  {
    img->m7[0][0] = sign( level, img->m7[0][0]);
  }
  else
  {
    if(qp_per<4)
    {
      img->m7[0][0] = sign( ((level*invlevelscale[0][0]+(1<<(3-qp_per)))>>(4-qp_per)), img->m7[0][0]);
    }
    else
    {
      img->m7[0][0] = sign( ((level*invlevelscale[0][0])<<(qp_per-4)), img->m7[0][0]);
    }
  }

  for (coeff_ctr=1;coeff_ctr < 16;coeff_ctr++)
  {
    i=SNGL_SCAN[coeff_ctr][0];
    j=SNGL_SCAN[coeff_ctr][1];

    run++;
    ilev=0;

    if(lossless_qpprime)
      level = absm (img->m7[j][i]);
    else 
      level = (absm(img->m7[j][i])*levelscale[i][j]+leveloffset[i][j])>>q_bits;
    
    if (level != 0)
    {
      if(i||j) nonzeroAC=TRUE;
      
      ACLevel[scan_pos] = sign(level,img->m7[j][i]);
      ACRun  [scan_pos] = run;
      ++scan_pos;
      run=-1;                     // reset zero level counter
      
      level=sign(level, img->m7[j][i]);
      if(lossless_qpprime)
      {
        ilev=level;
      }
      else if(qp_per<4)
      {
        ilev=(level*invlevelscale[i][j]+(1<<(3-qp_per)))>>(4-qp_per);
      }
      else
      {
        ilev=(level*invlevelscale[i][j])<<(qp_per-4);
      }
    }
    if(!lossless_qpprime)
      img->m7[j][i]=ilev;
  }
  ACLevel[scan_pos] = 0;

  
  //     IDCT.
  //     horizontal
  if(!lossless_qpprime)
  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < BLOCK_SIZE; i++)
    {
      m5[i]=img->m7[j][i];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (i=0; i < 2; i++)
    {
      i1=3-i;
      img->m7[j][i]=m6[i]+m6[i1];
      img->m7[j][i1]=m6[i]-m6[i1];
    }
  }

  //  vertical
  if(!lossless_qpprime)
  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < BLOCK_SIZE; j++)
    {
      m5[j]=img->m7[j][i];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (j=0; j < 2; j++)
    {
      j1=3-j;
      img->m7[j][i] =(m6[j]+m6[j1]+DQ_ROUND)>>DQ_BITS;
      img->m7[j1][i]=(m6[j]-m6[j1]+DQ_ROUND)>>DQ_BITS;
    }
  }

  return nonzeroAC;
}

// Residue Color Transform
int dct_chroma_DC(int uv, int cr_cbp)
{
  int run, scan_pos, coeff_ctr, level, i, j;
  int*  DCLevel = img->cofDC[uv+1][0];
  int*  DCRun   = img->cofDC[uv+1][1];

  run=-1;
  scan_pos=0;

  for (coeff_ctr=0; coeff_ctr < 16; coeff_ctr++)
  {
    i=SNGL_SCAN[coeff_ctr][0];
    j=SNGL_SCAN[coeff_ctr][1];

    run++;

    level = absm(dc_level[uv][i][j]);

    if (level  != 0)
    {
      cr_cbp=max(1,cr_cbp);
      DCLevel[scan_pos] = sign(level ,dc_level[uv][i][j]);
      DCRun  [scan_pos] = run;
      scan_pos++;
      run=-1;
    }
  }
  DCLevel[scan_pos] = 0;

  return cr_cbp;
}


/*!
 ************************************************************************
 * \brief
 *    The routine performs transform,quantization,inverse transform, adds the diff.
 *    to the prediction and writes the result to the decoded luma frame. Includes the
 *    RD constrained quantization also.
 *
 * \par Input:
 *    block_x,block_y: Block position inside a macro block (0,4,8,12).
 *
 * \par Output:
 *    nonzero: 0 if no levels are nonzero.  1 if there are nonzero levels.              \n
 *    coeff_cost: Counter for nonzero coefficients, used to discard expensive levels.
 *
 *
 ************************************************************************
 */
int dct_luma_sp(int block_x,int block_y,int *coeff_cost)
{
  int sign(int a,int b);

  int i,j,i1,j1,ilev,m5[4],m6[4],coeff_ctr;
  int qp_const,level,scan_pos,run;
  int nonzero;

  int predicted_block[BLOCK_SIZE][BLOCK_SIZE],c_err,qp_const2;
  int qp_per,qp_rem,q_bits;
  int qp_per_sp,qp_rem_sp,q_bits_sp;

  int   pos_x   = block_x >> BLOCK_SHIFT;
  int   pos_y   = block_y >> BLOCK_SHIFT;
  int   b8      = 2*(pos_y >> 1) + (pos_x >> 1);
  int   b4      = 2*(pos_y & 0x01) + (pos_x & 0x01);
  int*  ACLevel = img->cofAC[b8][b4][0];
  int*  ACRun   = img->cofAC[b8][b4][1];
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  short is_field_mode = (img->field_picture || ( img->MbaffFrameFlag && currMB->mb_field));

  // For encoding optimization
  int c_err1, c_err2, level1, level2;
  double D_dis1, D_dis2;
  int len, info;
  double lambda_mode   = 0.85 * pow (2, (currMB->qp - SHIFT_QP)/3.0) * 4; 

  qp_per    = (currMB->qp-MIN_QP)/6;
  qp_rem    = (currMB->qp-MIN_QP)%6;
  q_bits    = Q_BITS+qp_per;
  qp_per_sp    = (currMB->qpsp-MIN_QP)/6;
  qp_rem_sp    = (currMB->qpsp-MIN_QP)%6;
  q_bits_sp    = Q_BITS+qp_per_sp;

  qp_const=(1<<q_bits)/6;    // inter
  qp_const2=(1<<q_bits_sp)/2;  //sp_pred

  //  Horizontal transform
  for (j=0; j< BLOCK_SIZE; j++)
    for (i=0; i< BLOCK_SIZE; i++)
    {
      img->m7[j][i]+=img->mpr[j+block_y][i+block_x];
      predicted_block[i][j]=img->mpr[j+block_y][i+block_x];
    }

  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < 2; i++)
    {
      i1=3-i;
      m5[i]=img->m7[j][i]+img->m7[j][i1];
      m5[i1]=img->m7[j][i]-img->m7[j][i1];
    }
    img->m7[j][0]=(m5[0]+m5[1]);
    img->m7[j][2]=(m5[0]-m5[1]);
    img->m7[j][1]=m5[3]*2+m5[2];
    img->m7[j][3]=m5[3]-m5[2]*2;
  }

  //  Vertical transform

  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < 2; j++)
    {
      j1=3-j;
      m5[j]=img->m7[j][i]+img->m7[j1][i];
      m5[j1]=img->m7[j][i]-img->m7[j1][i];
    }
    img->m7[0][i]=(m5[0]+m5[1]);
    img->m7[2][i]=(m5[0]-m5[1]);
    img->m7[1][i]=m5[3]*2+m5[2];
    img->m7[3][i]=m5[3]-m5[2]*2;
  }

  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < 2; i++)
    {
      i1=3-i;
      m5[i]=predicted_block[i][j]+predicted_block[i1][j];
      m5[i1]=predicted_block[i][j]-predicted_block[i1][j];
    }
    predicted_block[0][j]=(m5[0]+m5[1]);
    predicted_block[2][j]=(m5[0]-m5[1]);
    predicted_block[1][j]=m5[3]*2+m5[2];
    predicted_block[3][j]=m5[3]-m5[2]*2;
  }

  //  Vertical transform

  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < 2; j++)
    {
      j1=3-j;
      m5[j]=predicted_block[i][j]+predicted_block[i][j1];
      m5[j1]=predicted_block[i][j]-predicted_block[i][j1];
    }
    predicted_block[i][0]=(m5[0]+m5[1]);
    predicted_block[i][2]=(m5[0]-m5[1]);
    predicted_block[i][1]=m5[3]*2+m5[2];
    predicted_block[i][3]=m5[3]-m5[2]*2;
  }

  // Quant
  nonzero=FALSE;

  run=-1;
  scan_pos=0;
  
  for (coeff_ctr=0;coeff_ctr < 16;coeff_ctr++)     // 8 times if double scan, 16 normal scan
  {

    if (is_field_mode) 
    {  // Alternate scan for field coding
        i=FIELD_SCAN[coeff_ctr][0];
        j=FIELD_SCAN[coeff_ctr][1];
    }
    else 
    {
        i=SNGL_SCAN[coeff_ctr][0];
        j=SNGL_SCAN[coeff_ctr][1];
    }
    
    run++;
    ilev=0;
    
    // decide prediction
    
    // case 1
    level1 = (absm (predicted_block[i][j]) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp; 
    level1 = (level1 << q_bits_sp) / quant_coef[qp_rem_sp][i][j];                 
    c_err1 = img->m7[j][i]-sign(level1, predicted_block[i][j]);                   
    level1 = (absm (c_err1) * quant_coef[qp_rem][i][j] + qp_const) >> q_bits;
    
    // case 2
    c_err2=img->m7[j][i]-predicted_block[i][j];
    level2 = (absm (c_err2) * quant_coef[qp_rem][i][j] + qp_const) >> q_bits;
    
    // select prediction
    if ((level1 != level2) && (level1 != 0) && (level2 != 0))
    {
      D_dis1 = img->m7[j][i] - ((sign(level1,c_err1)*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6) - predicted_block[i][j]; 
      levrun_linfo_inter(level1, run, &len, &info);
      D_dis1 = D_dis1*D_dis1 + lambda_mode * len;
      
      D_dis2 = img->m7[j][i] - ((sign(level2,c_err2)*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6) - predicted_block[i][j]; 
      levrun_linfo_inter(level2, run, &len, &info);
      D_dis2 = D_dis2 * D_dis2 + lambda_mode * len;
      
      if (D_dis1 == D_dis2)
        level = (absm(level1) < absm(level2)) ? level1 : level2;
      else
      {
        if (D_dis1 < D_dis2)
          level = level1;
        else
          level = level2;
      }
      c_err = (level == level1) ? c_err1 : c_err2;
    }
    else if (level1 == level2)
    {
      level = level1;
      c_err = c_err1;
    }
    else
    {
      level = (level1 == 0) ? level1 : level2;
      c_err = (level1 == 0) ? c_err1 : c_err2;
    }
    
    if (level != 0)
    {
      nonzero=TRUE;
      if (level > 1)
        *coeff_cost += MAX_VALUE;                // set high cost, shall not be discarded
      else
        *coeff_cost += COEFF_COST[input->disthres][run];
      ACLevel[scan_pos] = sign(level,c_err);
      ACRun  [scan_pos] = run;
      ++scan_pos;
      run=-1;                     // reset zero level counter
      ilev=((sign(level,c_err)*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6);
    }
    ilev+=predicted_block[i][j] ; 
    if(!si_frame_indicator && !sp2_frame_indicator)//stores the SP frame coefficients in lrec, will be useful to encode these and create SI or SP switching frame
    {
      lrec[img->pix_y+block_y+j][img->pix_x+block_x+i]= 
        sign((abs(ilev) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp, ilev); 
    }
    img->m7[j][i] = sign((absm(ilev) * quant_coef[qp_rem_sp][i][j] + qp_const2)>> q_bits_sp, ilev) * dequant_coef[qp_rem_sp][i][j] << qp_per_sp;
  }
  ACLevel[scan_pos] = 0;
  
    
  //     IDCT.
  //     horizontal

  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < BLOCK_SIZE; i++)
    {
      m5[i]=img->m7[j][i];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (i=0; i < 2; i++)
    {
      i1=3-i;
      img->m7[j][i]=m6[i]+m6[i1];
      img->m7[j][i1]=m6[i]-m6[i1];
    }
  }

  //  vertical

  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < BLOCK_SIZE; j++)
    {
      m5[j]=img->m7[j][i];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (j=0; j < 2; j++)
    {
      j1=3-j;
      img->m7[j][i] =min(img->max_imgpel_value,max(0,(m6[j]+m6[j1]+DQ_ROUND)>>DQ_BITS));
      img->m7[j1][i]=min(img->max_imgpel_value,max(0,(m6[j]-m6[j1]+DQ_ROUND)>>DQ_BITS));
    }
  }

  //  Decoded block moved to frame memory

  for (j=0; j < BLOCK_SIZE; j++)
  for (i=0; i < BLOCK_SIZE; i++)
    enc_picture->imgY[img->pix_y+block_y+j][img->pix_x+block_x+i]=img->m7[j][i];

  return nonzero;
}

/*!
 ************************************************************************
 * \brief
 *    Transform,quantization,inverse transform for chroma.
 *    The main reason why this is done in a separate routine is the
 *    additional 2x2 transform of DC-coeffs. This routine is called
 *    ones for each of the chroma components.
 *
 * \par Input:
 *    uv    : Make difference between the U and V chroma component               \n
 *    cr_cbp: chroma coded block pattern
 *
 * \par Output:
 *    cr_cbp: Updated chroma coded block pattern.
 ************************************************************************
 */
int dct_chroma_sp(int uv,int cr_cbp)
{
  int i,j,i1,j2,ilev,n2,n1,j1,mb_y,coeff_ctr,qp_const,c_err,level ,scan_pos,run;
  int m1[BLOCK_SIZE],m5[BLOCK_SIZE],m6[BLOCK_SIZE];
  int coeff_cost;
  int cr_cbp_tmp;
  int predicted_chroma_block[MB_BLOCK_SIZE>>1][MB_BLOCK_SIZE>>1],qp_const2,mp1[BLOCK_SIZE];
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  short is_field_mode = (img->field_picture || ( img->MbaffFrameFlag && currMB->mb_field));

  int qp_per,qp_rem,q_bits;
  int qp_per_sp,qp_rem_sp,q_bits_sp;

  int   b4;
  int*  DCLevel = img->cofDC[uv+1][0];
  int*  DCRun   = img->cofDC[uv+1][1];
  int*  ACLevel;
  int*  ACRun;

  int c_err1, c_err2, level1, level2;
  int len, info;
  double D_dis1, D_dis2;
  double lambda_mode   = 0.85 * pow (2, (currMB->qp -SHIFT_QP)/3.0) * 4; 


  int qpChroma = Clip3(-img->bitdepth_chroma_qp_scale, 51, currMB->qp + active_pps->chroma_qp_index_offset);
  int qpChromaSP=Clip3(-img->bitdepth_chroma_qp_scale, 51, currMB->qpsp + active_pps->chroma_qp_index_offset);

  qp_per    = ((qpChroma<0?qpChroma:QP_SCALE_CR[qpChroma])-MIN_QP)/6;
  qp_rem    = ((qpChroma<0?qpChroma:QP_SCALE_CR[qpChroma])-MIN_QP)%6;
  q_bits    = Q_BITS+qp_per;
  qp_const=(1<<q_bits)/6;    // inter
  qp_per_sp    = ((qpChromaSP<0?currMB->qpsp:QP_SCALE_CR[qpChromaSP])-MIN_QP)/6;
  qp_rem_sp    = ((qpChromaSP<0?currMB->qpsp:QP_SCALE_CR[qpChromaSP])-MIN_QP)%6;
  q_bits_sp    = Q_BITS+qp_per_sp;
  qp_const2=(1<<q_bits_sp)/2;  //sp_pred


  for (j=0; j < MB_BLOCK_SIZE>>1; j++)
    for (i=0; i < MB_BLOCK_SIZE>>1; i++)
    {
      img->m7[j][i]+=img->mpr[j][i];
      predicted_chroma_block[i][j]=img->mpr[j][i];
    }

  for (n2=0; n2 <= BLOCK_SIZE; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 <= BLOCK_SIZE; n1 += BLOCK_SIZE)
    {

      //  Horizontal transform.
      for (j=0; j < BLOCK_SIZE; j++)
      {
        mb_y=n2+j;
        for (i=0; i < 2; i++)
        {
          i1=3-i;
          m5[i]=img->m7[mb_y][i+n1]+img->m7[mb_y][i1+n1];
          m5[i1]=img->m7[mb_y][i+n1]-img->m7[mb_y][i1+n1];
        }
        img->m7[mb_y][n1]  =(m5[0]+m5[1]);
        img->m7[mb_y][n1+2]=(m5[0]-m5[1]);
        img->m7[mb_y][n1+1]=m5[3]*2+m5[2];
        img->m7[mb_y][n1+3]=m5[3]-m5[2]*2;
      }

      //  Vertical transform.

      for (i=0; i < BLOCK_SIZE; i++)
      {
        j1=n1+i;
        for (j=0; j < 2; j++)
        {
          j2=3-j;
          m5[j]=img->m7[n2+j][j1]+img->m7[n2+j2][j1];
          m5[j2]=img->m7[n2+j][j1]-img->m7[n2+j2][j1];
        }
        img->m7[n2+0][j1]=(m5[0]+m5[1]);
        img->m7[n2+2][j1]=(m5[0]-m5[1]);
        img->m7[n2+1][j1]=m5[3]*2+m5[2];
        img->m7[n2+3][j1]=m5[3]-m5[2]*2;
      }
    }
  }
  for (n2=0; n2 <= BLOCK_SIZE; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 <= BLOCK_SIZE; n1 += BLOCK_SIZE)
    {

      //  Horizontal transform.
      for (j=0; j < BLOCK_SIZE; j++)
      {
        mb_y=n2+j;
        for (i=0; i < 2; i++)
        {
          i1=3-i;
          m5[i]=predicted_chroma_block[i+n1][mb_y]+predicted_chroma_block[i1+n1][mb_y];
          m5[i1]=predicted_chroma_block[i+n1][mb_y]-predicted_chroma_block[i1+n1][mb_y];
        }
        predicted_chroma_block[n1][mb_y]  =(m5[0]+m5[1]);
        predicted_chroma_block[n1+2][mb_y]=(m5[0]-m5[1]);
        predicted_chroma_block[n1+1][mb_y]=m5[3]*2+m5[2];
        predicted_chroma_block[n1+3][mb_y]=m5[3]-m5[2]*2;
      }

      //  Vertical transform.

      for (i=0; i < BLOCK_SIZE; i++)
      {
        j1=n1+i;
        for (j=0; j < 2; j++)
        {
          j2=3-j;
          m5[j]=predicted_chroma_block[j1][n2+j]+predicted_chroma_block[j1][n2+j2];
          m5[j2]=predicted_chroma_block[j1][n2+j]-predicted_chroma_block[j1][n2+j2];
        }
        predicted_chroma_block[j1][n2+0]=(m5[0]+m5[1]);
        predicted_chroma_block[j1][n2+2]=(m5[0]-m5[1]);
        predicted_chroma_block[j1][n2+1]=m5[3]*2+m5[2];
        predicted_chroma_block[j1][n2+3]=m5[3]-m5[2]*2;
      }
    }
  }

  //     2X2 transform of DC coeffs.
  m1[0]=(img->m7[0][0]+img->m7[0][4]+img->m7[4][0]+img->m7[4][4]);
  m1[1]=(img->m7[0][0]-img->m7[0][4]+img->m7[4][0]-img->m7[4][4]);
  m1[2]=(img->m7[0][0]+img->m7[0][4]-img->m7[4][0]-img->m7[4][4]);
  m1[3]=(img->m7[0][0]-img->m7[0][4]-img->m7[4][0]+img->m7[4][4]);

  //     2X2 transform of DC coeffs.
  mp1[0]=(predicted_chroma_block[0][0]+predicted_chroma_block[4][0]+predicted_chroma_block[0][4]+predicted_chroma_block[4][4]);
  mp1[1]=(predicted_chroma_block[0][0]-predicted_chroma_block[4][0]+predicted_chroma_block[0][4]-predicted_chroma_block[4][4]);
  mp1[2]=(predicted_chroma_block[0][0]+predicted_chroma_block[4][0]-predicted_chroma_block[0][4]-predicted_chroma_block[4][4]);
  mp1[3]=(predicted_chroma_block[0][0]-predicted_chroma_block[4][0]-predicted_chroma_block[0][4]+predicted_chroma_block[4][4]);

  run=-1;
  scan_pos=0;

  for (coeff_ctr=0; coeff_ctr < 4; coeff_ctr++)
  {
    run++;
    ilev=0;

  // case 1
    c_err1 = (absm (mp1[coeff_ctr]) * quant_coef[qp_rem_sp][0][0] + 2 * qp_const2) >> (q_bits_sp + 1);
    c_err1 = (c_err1 << (q_bits_sp + 1)) / quant_coef[qp_rem_sp][0][0];
    c_err1 = m1[coeff_ctr] - sign(c_err1, mp1[coeff_ctr]);
    level1 = (absm(c_err1) * quant_coef[qp_rem][0][0] + 2 * qp_const) >> (q_bits+1);

  // case 2
    c_err2 = m1[coeff_ctr] - mp1[coeff_ctr];
    level2 = (absm(c_err2) * quant_coef[qp_rem][0][0] + 2 * qp_const) >> (q_bits+1);

    if (level1 != level2 && level1 != 0 && level2 != 0)
    {
      D_dis1 = m1[coeff_ctr] - ((sign(level1,c_err1)*dequant_coef[qp_rem][0][0]*A[0][0]<< qp_per) >>5)- mp1[coeff_ctr];
      levrun_linfo_c2x2(level1, run, &len, &info);
      D_dis1 = D_dis1 * D_dis1 + lambda_mode * len;
      
      D_dis2 = m1[coeff_ctr] - ((sign(level2,c_err2)*dequant_coef[qp_rem][0][0]*A[0][0]<< qp_per) >>5)- mp1[coeff_ctr];
      levrun_linfo_c2x2(level2, run, &len, &info);
      D_dis2 = D_dis2 * D_dis2 + lambda_mode * len;
      
      if (D_dis1 == D_dis2)
        level = (absm(level1) < absm(level2)) ? level1 : level2;
      else
      {
        if (D_dis1 < D_dis2)
          level = level1;
        else
          level = level2;
      }
      c_err = (level == level1) ? c_err1 : c_err2;
    }
    else if (level1 == level2)
    {
      level = level1;
      c_err = c_err1;
    }
    else
    {
      level = (level1 == 0) ? level1 : level2;
      c_err = (level1 == 0) ? c_err1 : c_err2;
    }
    
    if (input->symbol_mode == UVLC && img->qp < 4) 
    {
      if (level > CAVLC_LEVEL_LIMIT) 
      {
        level = CAVLC_LEVEL_LIMIT;
      }
    }

    if (level  != 0)
    {
      currMB->cbp_blk |= 0xf0000 << (uv << 2) ;  // if one of the 2x2-DC levels is != 0 the coded-bit
      cr_cbp=max(1,cr_cbp);
      DCLevel[scan_pos] = sign(level ,c_err);
      DCRun  [scan_pos] = run;
      scan_pos++;
      run=-1;
      ilev=((sign(level,c_err)*dequant_coef[qp_rem][0][0]*A[0][0]<< qp_per) >>5);
    }
    ilev+= mp1[coeff_ctr];
    m1[coeff_ctr]=sign((absm(ilev)  * quant_coef[qp_rem_sp][0][0] + 2 * qp_const2) >> (q_bits_sp+1), ilev) * dequant_coef[qp_rem_sp][0][0] << qp_per_sp;
    if(!si_frame_indicator && !sp2_frame_indicator)
      lrec_uv[uv][img->pix_c_y+4*(coeff_ctr%2)][img->pix_c_x+4*(coeff_ctr/2)]=sign((absm(ilev)  * quant_coef[qp_rem_sp][0][0] + 2 * qp_const2) >> (q_bits_sp+1), ilev);// stores the SP frames coefficients, will be useful to encode SI or switching SP frame
  }
  DCLevel[scan_pos] = 0;

  //  Inverse transform of 2x2 DC levels

  img->m7[0][0]=(m1[0]+m1[1]+m1[2]+m1[3])/2;
  img->m7[0][4]=(m1[0]-m1[1]+m1[2]-m1[3])/2;
  img->m7[4][0]=(m1[0]+m1[1]-m1[2]-m1[3])/2;
  img->m7[4][4]=(m1[0]-m1[1]-m1[2]+m1[3])/2;

  //     Quant of chroma AC-coeffs.
  coeff_cost=0;
  cr_cbp_tmp=0;

  for (n2=0; n2 <= BLOCK_SIZE; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 <= BLOCK_SIZE; n1 += BLOCK_SIZE)
    {
      b4      = 2*(n2 >> 2) + (n1 >> 2);
      ACLevel = img->cofAC[uv+4][b4][0];
      ACRun   = img->cofAC[uv+4][b4][1];

      run      = -1;
      scan_pos =  0;

      for (coeff_ctr=1; coeff_ctr < 16; coeff_ctr++)// start change rd_quant
      {

        if (is_field_mode) 
        {  // Alternate scan for field coding
          i=FIELD_SCAN[coeff_ctr][0];
          j=FIELD_SCAN[coeff_ctr][1];
        }
        else 
        {
          i=SNGL_SCAN[coeff_ctr][0];
          j=SNGL_SCAN[coeff_ctr][1];
        }
        ++run;
        ilev=0;

    // quantization on prediction
    c_err1 = (absm(predicted_chroma_block[n1+i][n2+j]) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp;
    c_err1 = (c_err1 << q_bits_sp) / quant_coef[qp_rem_sp][i][j];
    c_err1 = img->m7[n2+j][n1+i] - sign(c_err1, predicted_chroma_block[n1+i][n2+j]);
    level1 = (absm(c_err1) * quant_coef[qp_rem][i][j] + qp_const) >> q_bits;

    // no quantization on prediction
    c_err2 = img->m7[n2+j][n1+i] - predicted_chroma_block[n1+i][n2+j];
    level2 = (absm(c_err2) * quant_coef[qp_rem][i][j] + qp_const) >> q_bits;

    if (level1 != level2 && level1 != 0 && level2 != 0)
    {
      D_dis1 = img->m7[n2+j][n1+i] - ((sign(level1,c_err1)*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6) - predicted_chroma_block[n1+i][n2+j]; 

      levrun_linfo_inter(level1, run, &len, &info);
      D_dis1 = D_dis1 * D_dis1 + lambda_mode * len;

      D_dis2 = img->m7[n2+j][n1+i] - ((sign(level2,c_err2)*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6) - predicted_chroma_block[n1+i][n2+j]; 
      levrun_linfo_inter(level2, run, &len, &info);
      D_dis2 = D_dis2 * D_dis2 + lambda_mode * len;
      
      if (D_dis1 == D_dis2)
        level = (absm(level1) < absm(level2)) ? level1 : level2;
      else
      {
        if (D_dis1 < D_dis2)
          level = level1;
        else
          level = level2;
      }
      c_err = (level == level1) ? c_err1 : c_err2;
    }
    else if (level1 == level2)
    {
      level = level1;
      c_err = c_err1;
    }
    else
    {
      level = (level1 == 0) ? level1 : level2;
      c_err = (level1 == 0) ? c_err1 : c_err2;
    }

        if (level  != 0)
        {
          currMB->cbp_blk |=  (int64)1 << (16 + (uv << 2) + ((n2 >> 1) + (n1 >> 2))) ;
          if (level > 1)
            coeff_cost += MAX_VALUE;                // set high cost, shall not be discarded
          else
            coeff_cost += COEFF_COST[input->disthres][run];

          cr_cbp_tmp=2;
          ACLevel[scan_pos] = sign(level,c_err);
          ACRun  [scan_pos] = run;
          ++scan_pos;
          run=-1;
          ilev=((sign(level,c_err)*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6);
        }
        ilev+=predicted_chroma_block[n1+i][n2+j];
        if(!si_frame_indicator && !sp2_frame_indicator)
          if(!( (n2+j) % 4==0 && (n1+i)%4 ==0 ))    
            lrec_uv[uv][img->pix_c_y+n1+j][img->pix_c_x+n2+i]=sign((absm(ilev) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp,ilev);//stores the SP frames coefficients, will be useful to encode SI or switching SP frame
        img->m7[n2+j][n1+i] = sign((absm(ilev) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp,ilev) * dequant_coef[qp_rem_sp][i][j] << qp_per_sp;
      }
      ACLevel[scan_pos] = 0;
    }
  }

  // * reset chroma coeffs

  if(cr_cbp_tmp==2)
      cr_cbp=2;
  //     IDCT.

      //     Horizontal.
  for (n2=0; n2 <= BLOCK_SIZE; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 <= BLOCK_SIZE; n1 += BLOCK_SIZE)
    {
      for (j=0; j < BLOCK_SIZE; j++)
      {
        for (i=0; i < BLOCK_SIZE; i++)
        {
          m5[i]=img->m7[n2+j][n1+i];
        }
        m6[0]=(m5[0]+m5[2]);
        m6[1]=(m5[0]-m5[2]);
        m6[2]=(m5[1]>>1)-m5[3];
        m6[3]=m5[1]+(m5[3]>>1);

        for (i=0; i < 2; i++)
        {
          i1=3-i;
          img->m7[n2+j][n1+i]=m6[i]+m6[i1];
          img->m7[n2+j][n1+i1]=m6[i]-m6[i1];
        }
      }

      //     Vertical.
      for (i=0; i < BLOCK_SIZE; i++)
      {
        for (j=0; j < BLOCK_SIZE; j++)
        {
          m5[j]=img->m7[n2+j][n1+i];
        }
        m6[0]=(m5[0]+m5[2]);
        m6[1]=(m5[0]-m5[2]);
        m6[2]=(m5[1]>>1)-m5[3];
        m6[3]=m5[1]+(m5[3]>>1);

        for (j=0; j < 2; j++)
        {
          j2=3-j;
          img->m7[n2+j][n1+i] =min(img->max_imgpel_value_uv,max(0,(m6[j]+m6[j2]+DQ_ROUND)>>DQ_BITS));
          img->m7[n2+j2][n1+i]=min(img->max_imgpel_value_uv,max(0,(m6[j]-m6[j2]+DQ_ROUND)>>DQ_BITS));
        }
      }
    }
  }

  //  Decoded block moved to memory
  for (j=0; j < BLOCK_SIZE*2; j++)
    for (i=0; i < BLOCK_SIZE*2; i++)
    {
      enc_picture->imgUV[uv][img->pix_c_y+j][img->pix_c_x+i]= img->m7[j][i];
    }

  return cr_cbp;
}

/*!
 ************************************************************************
 * \brief
 *    The routine performs transform,quantization,inverse transform, adds the diff.
 *    to the prediction and writes the result to the decoded luma frame. Includes the
 *    RD constrained quantization also.
 *
 * \par Input:
 *    block_x,block_y: Block position inside a macro block (0,4,8,12).
 *
 * \par Output:
 *    nonzero: 0 if no levels are nonzero.  1 if there are nonzero levels.            \n
 *    coeff_cost: Counter for nonzero coefficients, used to discard expencive levels.
 ************************************************************************
 */
void copyblock_sp(int block_x,int block_y)
{
  int sign(int a,int b);

  int i,j,i1,j1,m5[4],m6[4];

  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int predicted_block[BLOCK_SIZE][BLOCK_SIZE];
  int qp_per = (currMB->qpsp-MIN_QP)/6;
  int qp_rem = (currMB->qpsp-MIN_QP)%6;
  int q_bits    = Q_BITS+qp_per;
  int qp_const2=(1<<q_bits)/2;  //sp_pred

  //  Horizontal transform
  for (j=0; j< BLOCK_SIZE; j++)
    for (i=0; i< BLOCK_SIZE; i++)
    {
      predicted_block[i][j]=img->mpr[j+block_y][i+block_x];
    }

  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < 2; i++)
    {
      i1=3-i;
      m5[i]=predicted_block[i][j]+predicted_block[i1][j];
      m5[i1]=predicted_block[i][j]-predicted_block[i1][j];
    }
    predicted_block[0][j]=(m5[0]+m5[1]);
    predicted_block[2][j]=(m5[0]-m5[1]);
    predicted_block[1][j]=m5[3]*2+m5[2];
    predicted_block[3][j]=m5[3]-m5[2]*2;
  }

  //  Vertical transform

  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < 2; j++)
    {
      j1=3-j;
      m5[j]=predicted_block[i][j]+predicted_block[i][j1];
      m5[j1]=predicted_block[i][j]-predicted_block[i][j1];
    }
    predicted_block[i][0]=(m5[0]+m5[1]);
    predicted_block[i][2]=(m5[0]-m5[1]);
    predicted_block[i][1]=m5[3]*2+m5[2];
    predicted_block[i][3]=m5[3]-m5[2]*2;
  }

  // Quant
  for (j=0;j < BLOCK_SIZE; j++)
  {
    for (i=0; i < BLOCK_SIZE; i++)
    {
      img->m7[j][i]=sign((absm(predicted_block[i][j])* quant_coef[qp_rem][i][j]+qp_const2)>> q_bits,predicted_block[i][j])*dequant_coef[qp_rem][i][j]<<qp_per;
      if(!si_frame_indicator && !sp2_frame_indicator)
      {
        lrec[img->pix_y+block_y+j][img->pix_x+block_x+i] = 
        sign((abs(predicted_block[i][j]) * quant_coef[qp_rem][i][j] + qp_const2) >> q_bits, predicted_block[i][j]);// stores the SP frames coefficients, will be useful to encode SI or switching SP frame
      }
    }
  }

  //     IDCT.
  //     horizontal

  for (j=0;j<BLOCK_SIZE;j++)
  {
    for (i=0;i<BLOCK_SIZE;i++)
    {
      m5[i]=img->m7[j][i];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (i=0;i<2;i++)
    {
      i1=3-i;
      img->m7[j][i]=m6[i]+m6[i1];
      img->m7[j][i1]=m6[i]-m6[i1];
    }
  }
  // vertical
  for (i=0;i<BLOCK_SIZE;i++)
  {
    for (j=0;j<BLOCK_SIZE;j++)
      m5[j]=img->m7[j][i];

    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (j=0;j<2;j++)
    {
      j1=3-j;
      img->m7[j][i] =min(img->max_imgpel_value,max(0,(m6[j]+m6[j1]+DQ_ROUND)>>DQ_BITS));
      img->m7[j1][i]=min(img->max_imgpel_value,max(0,(m6[j]-m6[j1]+DQ_ROUND)>>DQ_BITS));
    }
  }

  //  Decoded block moved to frame memory

  for (j=0; j < BLOCK_SIZE; j++)
    for (i=0; i < BLOCK_SIZE; i++)
      enc_picture->imgY[img->pix_y+block_y+j][img->pix_x+block_x+i]=img->m7[j][i];
}



int writeIPCMBytes(Bitstream *currStream)
{
  int i,j, jj;
  int len = 0, uv;
  int             mb_nr     = img->current_mb_nr;
  Macroblock*     currMB    = &img->mb_data[mb_nr];
  SyntaxElement  *currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  
  
  for (j=0;j<16;j++)
  {
    jj = img->pix_y+j;
    for (i=0;i<16;i++)
    {
      currSE->len = img->bitdepth_luma;  
      len += currSE->len;
      currSE->bitpattern = enc_picture->imgY[jj][img->pix_x+i];
      writeSyntaxElement2Buf_Fixed(currSE, currStream);
    }
  }
  
  for (uv = 0; uv < 2; uv ++)
  {
    for (j=0;j<img->mb_cr_size_y;j++)
    {
      jj = img->pix_c_y+j;
      for (i=0;i<img->mb_cr_size_x;i++)
      {
        currSE->len = img->bitdepth_chroma;
        len += currSE->len;
        currSE->bitpattern = enc_picture->imgUV[uv][jj][img->pix_c_x+i];
        writeSyntaxElement2Buf_Fixed(currSE, currStream);
      }
    }
  }  
  return len;
}

int writePCMByteAlign(Bitstream *currStream)
{
  int len = 0;
  if (currStream->bits_to_go < 8)
  { // trailing bits to process
    len = 8 - currStream->bits_to_go;
    currStream->byte_buf = (currStream->byte_buf <<currStream->bits_to_go) | (0xff >> (8 - currStream->bits_to_go));
    stats->bit_use_stuffingBits[img->type]+=currStream->bits_to_go;
    currStream->streamBuffer[currStream->byte_pos++]=currStream->byte_buf;
    currStream->bits_to_go = 8;
  }
  return len;
}

/*!
 ************************************************************************
 * \brief Eric Setton 
 * Encoding of a secondary SP / SI frame.
 * For an SI frame the predicted block should only come from spatial pred.
 * The original image signal is the error coefficients of a primary SP in the raw data stream
 * the difference with the primary SP are :
 *  - the prediction signal is transformed and quantized (qpsp) but not dequantized 
 *  - only one kind of prediction is considered and not two 
 *  - the resulting error coefficients are not quantized before being sent to the VLC
 *
 * \para Input:
 *    block_x,block_y: Block position inside a macro block (0,4,8,12).
 *
 * \para Output:
 *    nonzero: 0 if no levels are nonzero.  1 if there are nonzero levels.              
 *    coeff_cost: Counter for nonzero coefficients, used to discard expencive levels.
 *
 *
 ************************************************************************
 */

int dct_luma_sp2(int block_x,int block_y,int *coeff_cost)
{
  int sign(int a,int b);

  int i,j,i1,j1,ilev,m5[4],m6[4],coeff_ctr;
  int qp_const,level,scan_pos,run;
  int nonzero;

  int predicted_block[BLOCK_SIZE][BLOCK_SIZE],c_err,qp_const2;
  int qp_per,qp_rem,q_bits;
  int qp_per_sp,qp_rem_sp,q_bits_sp;

  int   pos_x   = block_x >> BLOCK_SHIFT;
  int   pos_y   = block_y >> BLOCK_SHIFT;
  int   b8      = 2*(pos_y >> 1) + (pos_x >> 1);
  int   b4      = 2*(pos_y & 0x01) + (pos_x & 0x01);
  int*  ACLevel = img->cofAC[b8][b4][0];
  int*  ACRun   = img->cofAC[b8][b4][1];
  
  int level1;

  qp_per    = (img->qpsp-MIN_QP)/6 ;
  qp_rem    = (img->qpsp-MIN_QP)%6;
  q_bits    = Q_BITS+qp_per;
  qp_per_sp    = (img->qpsp-MIN_QP)/6;
  qp_rem_sp    = (img->qpsp-MIN_QP)%6;
  q_bits_sp    = Q_BITS+qp_per_sp;

  qp_const=(1<<q_bits)/6;    // inter
  qp_const2=(1<<q_bits_sp)/2;  //sp_pred

  //  Horizontal transform 
  for (j=0; j< BLOCK_SIZE; j++)
    for (i=0; i< BLOCK_SIZE; i++)
    {
      //Coefficients obtained from the prior encoding of the SP frame
      img->m7[j][i]=lrec[img->pix_y+block_y+j][img->pix_x+block_x+i]; 
      //Predicted block
      predicted_block[i][j]=img->mpr[j+block_y][i+block_x]; 
    }

  //Horizontal transform
  for (j=0; j < BLOCK_SIZE; j++)
    {
      for (i=0; i < 2; i++)
    {
      i1=3-i;
      m5[i]=predicted_block[i][j]+predicted_block[i1][j];
      m5[i1]=predicted_block[i][j]-predicted_block[i1][j];
    }
    predicted_block[0][j]=(m5[0]+m5[1]);
    predicted_block[2][j]=(m5[0]-m5[1]);
    predicted_block[1][j]=m5[3]*2+m5[2];
    predicted_block[3][j]=m5[3]-m5[2]*2;
  }

  //  Vertical transform of the predicted block 

  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < 2; j++)
    {
      j1=3-j;
      m5[j]=predicted_block[i][j]+predicted_block[i][j1];
      m5[j1]=predicted_block[i][j]-predicted_block[i][j1];
    }
    predicted_block[i][0]=(m5[0]+m5[1]);
    predicted_block[i][2]=(m5[0]-m5[1]);
    predicted_block[i][1]=m5[3]*2+m5[2];
    predicted_block[i][3]=m5[3]-m5[2]*2;
  }

  // Quant
  nonzero=FALSE;

  run=-1;
  scan_pos=0;
  
  for (coeff_ctr=0;coeff_ctr < 16;coeff_ctr++)     // 8 times if double scan, 16 normal scan
  {

    if (img->field_picture || ( mb_adaptive && img->field_mode )) 
    {  // Alternate scan for field coding
        i=FIELD_SCAN[coeff_ctr][0];
        j=FIELD_SCAN[coeff_ctr][1];
    }
    else 
    {
        i=SNGL_SCAN[coeff_ctr][0];
        j=SNGL_SCAN[coeff_ctr][1];
    }
    
    run++;
    ilev=0;
    
    //quantization of the predicted block
    level1 = (abs (predicted_block[i][j]) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp;  
    //substracted from lrec 
    c_err = img->m7[j][i]-sign(level1, predicted_block[i][j]);   //substracting the predicted block 

    
    level = abs(c_err);  
    if (level != 0)
    {
      nonzero=TRUE;
      if (level > 1)
        *coeff_cost += MAX_VALUE;                // set high cost, shall not be discarded
      else
        *coeff_cost += COEFF_COST[input->disthres][run];
      ACLevel[scan_pos] = sign(level,c_err);
      ACRun  [scan_pos] = run;
      ++scan_pos;
      run=-1;                     // reset zero level counter     
    }
    //from now on we are in decoder land
    ilev=c_err + sign(level1,predicted_block[i][j]) ;  // adding the quantized predicted block 
    img->m7[j][i] = ilev  *dequant_coef[qp_rem_sp][i][j] << qp_per_sp; 
    
  }
  ACLevel[scan_pos] = 0;
  
    
  //     IDCT.
  //     horizontal

  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < BLOCK_SIZE; i++)
    {
      m5[i]=img->m7[j][i];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (i=0; i < 2; i++)
    {
      i1=3-i;
      img->m7[j][i]=m6[i]+m6[i1];
      img->m7[j][i1]=m6[i]-m6[i1];
    }
  }

  //  vertical

  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < BLOCK_SIZE; j++)
    {
      m5[j]=img->m7[j][i];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (j=0; j < 2; j++)
    {
      j1=3-j;
      img->m7[j][i] =min(255,max(0,(m6[j]+m6[j1]+DQ_ROUND)>>DQ_BITS));
      img->m7[j1][i]=min(255,max(0,(m6[j]-m6[j1]+DQ_ROUND)>>DQ_BITS));
    }
  }

  //  Decoded block moved to frame memory
  for (j=0; j < BLOCK_SIZE; j++)
    for (i=0; i < BLOCK_SIZE; i++)
    {
      enc_picture->imgY[img->pix_y+block_y+i][img->pix_x+block_x+j]=img->m7[i][j]; 
    }
  return nonzero;
}


/*!
 ************************************************************************
 * \brief Eric Setton 
 * Encoding of the chroma of a  secondary SP / SI frame.
 * For an SI frame the predicted block should only come from spatial pred.
 * The original image signal is the error coefficients of a primary SP in the raw data stream
 * the difference with the primary SP are :
 *  - the prediction signal is transformed and quantized (qpsp) but not dequantized 
 *  - the resulting error coefficients are not quantized before being sent to the VLC
 *
 * \par Input:
 *    uv    : Make difference between the U and V chroma component             
 *    cr_cbp: chroma coded block pattern
 *
 * \par Output:
 *    cr_cbp: Updated chroma coded block pattern.
 *
 ************************************************************************
 */
int dct_chroma_sp2(int uv,int cr_cbp)
{
  int i,j,i1,j2,ilev,n2,n1,j1,mb_y,coeff_ctr,qp_const,c_err,level ,scan_pos,run;
  int m1[BLOCK_SIZE],m5[BLOCK_SIZE],m6[BLOCK_SIZE];
  int coeff_cost;
  int cr_cbp_tmp;
  int predicted_chroma_block[MB_BLOCK_SIZE/2][MB_BLOCK_SIZE/2],qp_const2,mp1[BLOCK_SIZE];
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int qp_per,qp_rem,q_bits;
  int qp_per_sp,qp_rem_sp,q_bits_sp;

  int   b4;
  int*  DCLevel = img->cofDC[uv+1][0];
  int*  DCRun   = img->cofDC[uv+1][1];
  int*  ACLevel;
  int*  ACRun;
  int  level1;
 
  qp_per    = ((img->qp<0?img->qp:QP_SCALE_CR[img->qp])-MIN_QP)/6;
  qp_rem    = ((img->qp<0?img->qp:QP_SCALE_CR[img->qp])-MIN_QP)%6;
  q_bits    = Q_BITS+qp_per;
  qp_const=(1<<q_bits)/6;    // inter

  qp_per_sp    = ((img->qpsp<0?img->qpsp:QP_SCALE_CR[img->qpsp])-MIN_QP)/6;
  qp_rem_sp    = ((img->qpsp<0?img->qpsp:QP_SCALE_CR[img->qpsp])-MIN_QP)%6;
  q_bits_sp    = Q_BITS+qp_per_sp;
  qp_const2=(1<<q_bits_sp)/2;  //sp_pred


  for (j=0; j < MB_BLOCK_SIZE>>1; j++)
    for (i=0; i < MB_BLOCK_SIZE>>1; i++)
    {
      predicted_chroma_block[i][j]=img->mpr[j][i];
      img->m7[j][i]=lrec_uv[uv][img->pix_c_y+j][img->pix_c_x+i];
    }

  for (n2=0; n2 <= BLOCK_SIZE; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 <= BLOCK_SIZE; n1 += BLOCK_SIZE)
    {

      //  Horizontal transform.
      for (j=0; j < BLOCK_SIZE; j++)
      {
        mb_y=n2+j;
        for (i=0; i < 2; i++)
        {
          i1=3-i;
          m5[i]=predicted_chroma_block[i+n1][mb_y]+predicted_chroma_block[i1+n1][mb_y];
          m5[i1]=predicted_chroma_block[i+n1][mb_y]-predicted_chroma_block[i1+n1][mb_y];
        }
        predicted_chroma_block[n1][mb_y]  =(m5[0]+m5[1]);
        predicted_chroma_block[n1+2][mb_y]=(m5[0]-m5[1]);
        predicted_chroma_block[n1+1][mb_y]=m5[3]*2+m5[2];
        predicted_chroma_block[n1+3][mb_y]=m5[3]-m5[2]*2;
      }

      //  Vertical transform.

      for (i=0; i < BLOCK_SIZE; i++)
      {
        j1=n1+i;
        for (j=0; j < 2; j++)
        {
          j2=3-j;
          m5[j]=predicted_chroma_block[j1][n2+j]+predicted_chroma_block[j1][n2+j2];
          m5[j2]=predicted_chroma_block[j1][n2+j]-predicted_chroma_block[j1][n2+j2];
        }
        predicted_chroma_block[j1][n2+0]=(m5[0]+m5[1]);
        predicted_chroma_block[j1][n2+2]=(m5[0]-m5[1]);
        predicted_chroma_block[j1][n2+1]=m5[3]*2+m5[2];
        predicted_chroma_block[j1][n2+3]=m5[3]-m5[2]*2;
      }
    }
  }

  //   DC coefficients already transformed and quantized 
  m1[0]= img->m7[0][0];
  m1[1]= img->m7[4][0];
  m1[2]= img->m7[0][4];
  m1[3]= img->m7[4][4];

  //     2X2 transform of predicted DC coeffs.
  mp1[0]=(predicted_chroma_block[0][0]+predicted_chroma_block[4][0]+predicted_chroma_block[0][4]+predicted_chroma_block[4][4]);
  mp1[1]=(predicted_chroma_block[0][0]-predicted_chroma_block[4][0]+predicted_chroma_block[0][4]-predicted_chroma_block[4][4]);
  mp1[2]=(predicted_chroma_block[0][0]+predicted_chroma_block[4][0]-predicted_chroma_block[0][4]-predicted_chroma_block[4][4]);
  mp1[3]=(predicted_chroma_block[0][0]-predicted_chroma_block[4][0]-predicted_chroma_block[0][4]+predicted_chroma_block[4][4]);

  run=-1;
  scan_pos=0;

  for (coeff_ctr=0; coeff_ctr < 4; coeff_ctr++)
  {
    run++;
    ilev=0;

    //quantization of predicted DC coeff  
    level1 = (abs (mp1[coeff_ctr]) * quant_coef[qp_rem_sp][0][0] + 2 * qp_const2) >> (q_bits_sp + 1);
    //substratcted from lrecUV
    c_err = m1[coeff_ctr] - sign(level1, mp1[coeff_ctr]);
    level = abs(c_err);

    if (level  != 0)
    {
      currMB->cbp_blk |= 0xf0000 << (uv << 2) ;  // if one of the 2x2-DC levels is != 0 the coded-bit
      cr_cbp=max(1,cr_cbp);
      DCLevel[scan_pos] = sign(level ,c_err);
      DCRun  [scan_pos] = run;
      scan_pos++;
      run=-1;
    }
    
    //from now on decoder world
    ilev = c_err + sign(level1,mp1[coeff_ctr]) ; // we have perfect reconstruction here

    m1[coeff_ctr]= ilev  * dequant_coef[qp_rem_sp][0][0] << qp_per_sp; 
  
  }
  DCLevel[scan_pos] = 0;

  //  Invers transform of 2x2 DC levels

  img->m7[0][0]=(m1[0]+m1[1]+m1[2]+m1[3])/2;
  img->m7[4][0]=(m1[0]-m1[1]+m1[2]-m1[3])/2; 
  img->m7[0][4]=(m1[0]+m1[1]-m1[2]-m1[3])/2; 
  img->m7[4][4]=(m1[0]-m1[1]-m1[2]+m1[3])/2;

  //     Quant of chroma AC-coeffs.
  coeff_cost=0;
  cr_cbp_tmp=0;

  for (n2=0; n2 <= BLOCK_SIZE; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 <= BLOCK_SIZE; n1 += BLOCK_SIZE)
    {
      b4      = 2*(n2/4) + (n1/4);
      ACLevel = img->cofAC[uv+4][b4][0];
      ACRun   = img->cofAC[uv+4][b4][1];

      run      = -1;
      scan_pos =  0;

      for (coeff_ctr=1; coeff_ctr < 16; coeff_ctr++)// start change rd_quant
      {

        if (img->field_picture || ( mb_adaptive && img->field_mode )) 
        {  // Alternate scan for field coding
          j=FIELD_SCAN[coeff_ctr][0];
          i=FIELD_SCAN[coeff_ctr][1];
        }
        else 
        {
          j=SNGL_SCAN[coeff_ctr][0];
          i=SNGL_SCAN[coeff_ctr][1];
        }
        ++run;
        ilev=0;
  // quantization on prediction
  level1 = (abs(predicted_chroma_block[n1+j][n2+i]) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp;
  //substracted from lrec 
  c_err  = img->m7[n1+i][n2+j] - sign(level1, predicted_chroma_block[n1+j][n2+i]); 
  level  = abs(c_err) ;

        if (level  != 0)
        {
          currMB->cbp_blk |=  (int64)1 << (16 + (uv << 2) + ((n2 >> 1) + (n1 >> 2))) ;
          if (level > 1)
            coeff_cost += MAX_VALUE;                // set high cost, shall not be discarded
          else
            coeff_cost += COEFF_COST[input->disthres][run];

          cr_cbp_tmp=2;
          ACLevel[scan_pos] = sign(level,c_err);
          ACRun  [scan_pos] = run;
          ++scan_pos;
          run=-1;
        }

  //from now on decoder land
        ilev=c_err + sign(level1,predicted_chroma_block[n1+j][n2+i]);
  img->m7[n1+i][n2+j] = ilev * dequant_coef[qp_rem_sp][i][j] << qp_per_sp;
      }
  ACLevel[scan_pos] = 0;      
    }
  }
  // * reset chroma coeffs

  if(cr_cbp_tmp==2)
      cr_cbp=2;
  //     IDCT.

      //     Horizontal.
  for (n2=0; n2 <= BLOCK_SIZE; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 <= BLOCK_SIZE; n1 += BLOCK_SIZE)
    {
      for (j=0; j < BLOCK_SIZE; j++)
      {
        for (i=0; i < BLOCK_SIZE; i++)
        {
          m5[i]=img->m7[n1+i][n2+j];
        }
        m6[0]=(m5[0]+m5[2]);
        m6[1]=(m5[0]-m5[2]);
        m6[2]=(m5[1]>>1)-m5[3];
        m6[3]=m5[1]+(m5[3]>>1);

        for (i=0; i < 2; i++)
        {
          i1=3-i;
          img->m7[n1+i][n2+j]=m6[i]+m6[i1];
          img->m7[n1+i1][n2+j]=m6[i]-m6[i1];
        }
      }

      //     Vertical.
      for (i=0; i < BLOCK_SIZE; i++)
      {
        for (j=0; j < BLOCK_SIZE; j++)
        {
          m5[j]=img->m7[n1+i][n2+j];
        }
        m6[0]=(m5[0]+m5[2]);
        m6[1]=(m5[0]-m5[2]);
        m6[2]=(m5[1]>>1)-m5[3];
        m6[3]=m5[1]+(m5[3]>>1);

        for (j=0; j < 2; j++)
        {
          j2=3-j;
          img->m7[n1+i][n2+j] =min(255,max(0,(m6[j]+m6[j2]+DQ_ROUND)>>DQ_BITS));
          img->m7[n1+i][n2+j2]=min(255,max(0,(m6[j]-m6[j2]+DQ_ROUND)>>DQ_BITS));
        }
      }
    }
  }

  //  Decoded block moved to memory
  for (j=0; j < BLOCK_SIZE; j++)
    for (i=0; i < BLOCK_SIZE; i++)
    {
  enc_picture->imgUV[uv][img->pix_c_y+i][img->pix_c_x+j]= img->m7[i][j];
  enc_picture->imgUV[uv][img->pix_c_y+i][img->pix_c_x+j+4]= img->m7[i+4][j];
  enc_picture->imgUV[uv][img->pix_c_y+i+4][img->pix_c_x+j]= img->m7[i][j+4];
  enc_picture->imgUV[uv][img->pix_c_y+i+4][img->pix_c_x+j+4]= img->m7[i+4][j+4];
    }
  
  return cr_cbp;
}

#ifdef USE_INTRA_MDDT

const int KLTCol[9][4][4]=
{
{ // 0
{ -42,  -61,  -73,  -74},
{  74,   65,  -16,  -81},
{ -80,   37,   73,  -57},
{ -53,   84,  -74,   33},
},
{ // 1
{ -35,  -62,  -79,  -71},
{  82,   65,  -24,  -70},
{ -78,   46,   59,  -68},
{ -49,   78,  -78,   42},
},
{ // 2
{ -47,  -63,  -72,  -71},
{  80,   57,  -25,  -78},
{ -75,   50,   67,  -62},
{ -46,   82,  -78,   37},
},
{ // 3
{ -30,  -60,  -79,  -75},
{  72,   73,  -10,  -76},
{ -84,   29,   69,  -61},
{ -56,   82,  -73,   34},
},
{ // 4
{ -30,  -61,  -80,  -74},
{  71,   73,  -11,  -77},
{ -85,   28,   68,  -62},
{ -57,   81,  -73,   35},
},
{ // 5
{ -29,  -57,  -79,  -78},
{  72,   74,   -7,  -75},
{ -83,   26,   71,  -61},
{ -58,   83,  -71,   32},
},
{ // 6
{ -32,  -61,  -79,  -74},
{  70,   73,  -14,  -77},
{ -85,   31,   67,  -61},
{ -56,   79,  -74,   38},
},
{ // 7
{ -34,  -61,  -77,  -75},
{  71,   71,  -10,  -79},
{  83,  -28,  -72,   59},
{  57,  -83,   72,  -32},
},
{ // 8
{ -45,  -61,  -71,  -74},
{  91,   51,  -31,  -68},
{ -68,   67,   54,  -66},
{ -38,   74,  -86,   45},
},
};

const int KLTRow[9][4][4]=
{
{ // 0
{ -41,   86,  -74,  -42},
{ -62,   56,   56,   79},
{ -75,  -26,   58,  -82},
{ -72,  -71,  -66,   41},
},
{ // 1
{ -37,   75,   79,   56},
{ -61,   67,  -33,  -84},
{ -74,  -13,  -74,   73},
{ -76,  -78,   60,  -30},
},
{ // 2
{ -44,   84,  -73,  -45},
{ -63,   56,   51,   82},
{ -73,  -25,   65,  -79},
{ -72,  -75,  -65,   36},
},
{ // 3
{ -47,   87,  -71,  -39},
{ -65,   50,   64,   75},
{ -73,  -34,   52,  -85},
{ -68,  -71,  -68,   46},
},
{ // 4
{ -31,   72,  -81,  -60},
{ -60,   73,   29,   81},
{ -78,  -13,   70,  -72},
{ -75,  -75,  -63,   34},
},
{ // 5
{ -34,   76,  -82,  -53},
{ -61,   70,   39,   79},
{ -78,  -18,   65,  -76},
{ -74,  -74,  -63,   39},
},
{ // 6
{ -31,   73,   83,   57},
{ -60,   71,  -28,  -83},
{ -77,   -9,  -71,   73},
{ -76,  -77,   60,  -31},
},
{ // 7
{ -46,   86,  -72,  -40},
{ -65,   52,   62,   75},
{ -73,  -34,   53,  -84},
{ -68,  -72,  -68,   46},
},
{ // 8
{ -34,   75,   81,   55},
{ -60,   69,  -33,  -83},
{ -76,  -12,  -72,   73},
{ -76,  -77,   61,  -32},
},
};

/*! 
*************************************************************************************
* \brief
*   Perform mode-dependent 4x4 directional transform 
*
* \para klt_luma_sep()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
int klt_luma_sep(int block_x,int block_y,int *coeff_cost, int ipmode)
{
  int sign(int a,int b);

  int i,j,k,coeff_ctr;
  int scan_pos,run;
  int nonzero;
  int qp_per,qp_rem, q_bits;

  int   pos_x   = block_x/BLOCK_SIZE;
  int   pos_y   = block_y/BLOCK_SIZE;
  int   b8      = 2*(pos_y/2) + (pos_x/2);
  int   b4      = 2*(pos_y%2) + (pos_x%2);
  int*  ACLevel = img->cofAC[b8][b4][0];
  int*  ACRun   = img->cofAC[b8][b4][1];
  byte (*pos_scan)[2] = img->scanOrder4x4[ipmode];
  int level;

  int coeff [BLOCK_SIZE][BLOCK_SIZE], ilev[BLOCK_SIZE][BLOCK_SIZE], 
      temp2D[BLOCK_SIZE][BLOCK_SIZE], y2D [BLOCK_SIZE][BLOCK_SIZE];

  int m=4, n=10, KLTprec=14;
  int qpOffset;
  int roundOffset;
  int R[6] = {40, 45, 50, 57, 63, 71};
  int Q[6] = {410,   364,   328,   287,   260,   231};

  long shift;
  int shiftBy2;

  static int first=1;

  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

#ifdef RDO_Q
  levelDataStruct levelData[16];
  double lambda_md = 0.0;
  double err;
  int lowerInt;
  int levelRDOQ[16];
  int kStart = 0;
  int kStop = 0;
  int noCoeff;
  int estBits;
  double normFact = pow(2, 2 * (KLTprec + m + n) - 15);
#endif

  pos_scan = img->scanOrder4x4[ipmode];

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
    img->delta_shift=1<<(n+qp_per+KLTprec-11); 
    img->QPVal=1<<(n+qp_per+KLTprec-3);  
    for (i=0; i<16; i++)
    {
      img->qp_shift[i]=(1<<(n+qp_per+KLTprec-2))/3;
    }
  }

  // column transform first
  for (i=0;i<(BLOCK_SIZE);i++)
  {
    for (j=0;j<(BLOCK_SIZE);j++)
    {
      temp2D[i][j] = 0;
      for(k = 0; k < BLOCK_SIZE; k++)
        temp2D[i][j]+=KLTCol[ipmode][i][k]*(long)(img->m7[k][j]);
    }
  }

  // row transform second
  for (i=0;i<(BLOCK_SIZE);i++)
  {
    for (j=0;j<(BLOCK_SIZE);j++)
    {
      coeff[i][j] = 0;
      for(k = 0; k < BLOCK_SIZE; k++)
        coeff[i][j]+=temp2D[i][k]*KLTRow[ipmode][k][j];
    }
  }

#ifdef RDO_Q
  if(input->UseRDO_Q)
  {
    noCoeff = 0;
    for (coeff_ctr = 0; coeff_ctr < 16; coeff_ctr++)
    {
      i = pos_scan[coeff_ctr][0];
      j = pos_scan[coeff_ctr][1];

      levelData[coeff_ctr].levelDouble = (int64)absm(coeff[j][i]) * Q[qp_rem];
      level = (int)(levelData[coeff_ctr].levelDouble >> qpOffset);
      lowerInt = ((levelData[coeff_ctr].levelDouble - (level << qpOffset)) < ((int64)1 << (qpOffset - 1)))? 1: 0;

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

    estBits = est_write_and_store_CBP_block_bit(currMB, LUMA_4x4);
    est_writeRunLevel_CABAC(levelData, levelRDOQ, LUMA_4x4, lambda_md, kStart, kStop, noCoeff, estBits);
  }
#endif

  // quantization and zig-zag scanning 
  nonzero=FALSE;

  run=-1;
  scan_pos=0;

  for (coeff_ctr=0;coeff_ctr < 16;coeff_ctr++)
  {
    i=pos_scan[coeff_ctr][0];
    j=pos_scan[coeff_ctr][1];
    
    run++;

#ifdef RDO_Q
    if(input->UseRDO_Q)
    {
      level = levelRDOQ[coeff_ctr];
    }
    else
#endif
    if (img->AdaptiveRounding)
    {
      level = (int)((absm((int64)coeff[j][i]) * Q[qp_rem] + img->qp_shift[coeff_ctr]) >> qpOffset);
      if (level != 0) 
      {
        quant_stat_rd[coeff_ctr] = (long)(absm((int64)coeff[j][i]) * Q[qp_rem] - ((int64)level << qpOffset));
      }
      else if (2 * absm((int64)coeff[j][i]) * Q[qp_rem] > ((int64)1 << qpOffset)) 
      {
        quant_stat_rd[coeff_ctr] = (long)(2 * absm((int64)coeff[j][i]) * Q[qp_rem] - ((int64)1 << qpOffset));
      }
    }    
    else
    {
      level = (int)((absm((int64)coeff[j][i]) * Q[qp_rem] + roundOffset) >> qpOffset);
    }

    if (level != 0)
    {
      nonzero=TRUE;

      *coeff_cost += (level > 1) ? MAX_VALUE : COEFF_COST[input->disthres][run];

      ACLevel[scan_pos] = sign(level, coeff[j][i]);
      ACRun  [scan_pos] = run;
      ++scan_pos;
      run=-1;                     // reset zero level counter

      ilev[j][i] = sign(level, coeff[j][i]) * (R[qp_rem] << qp_per);   
    }
    else 
      ilev[j][i]=0;
  }
  ACLevel[scan_pos] = 0;

  // inverse transform, column first
  for (i=0;i<(BLOCK_SIZE);i++)
  {
    for (j=0;j<(BLOCK_SIZE);j++)
    {
      temp2D[i][j] = 0;
      for(k = 0; k < BLOCK_SIZE; k++)
        temp2D[i][j]+=KLTCol[ipmode][k][i]*(ilev[k][j]);
    }
  }

  // inverse transform, row second 
  for (i=0;i<(BLOCK_SIZE);i++)
  {
    for (j=0;j<(BLOCK_SIZE);j++)
    {
      y2D[i][j] = 0;
      for(k = 0; k < BLOCK_SIZE; k++)
        y2D[i][j]+=temp2D[i][k]*KLTRow[ipmode][j][k];
    }
  }

  // store the encoded image
  shift=KLTprec+m+2; shiftBy2=1<<(shift-1);
  for (j=0; j < BLOCK_SIZE; j++)
    for (i=0; i < BLOCK_SIZE; i++)
    {
      short temp;

      temp=(short)((y2D[j][i]+((int64)img->mpr[j+block_y][i+block_x]<<shift)+shiftBy2)>>shift);
      temp=min(img->max_imgpel_value,max(0,temp));
      enc_picture->imgY[img->pix_y+block_y+j][img->pix_x+block_x+i]=temp;
    }

  return nonzero;
}

const int KLTCol16[4][16][16]=
{
{
{ -26,  -27,  -28,  -30,  -31,  -32,  -32,  -33,  -33,  -34,  -34,  -34,  -34,  -34,  -34,  -34},
{  39,   41,   40,   38,   32,   27,   21,   16,   -3,   -9,  -17,  -24,  -34,  -39,  -44,  -46},
{  36,   36,   29,   20,    0,  -12,  -24,  -32,  -49,  -45,  -36,  -22,    6,   24,   40,   49},
{ -40,  -36,  -22,   -3,   32,   46,   51,   45,  -21,  -33,  -37,  -33,  -13,    3,   19,   29},
{ -21,  -14,    3,   18,   29,   19,    0,  -19,  -52,  -32,    6,   43,   60,   33,  -17,  -58},
{  50,   26,  -18,  -53,  -51,  -13,   34,   57,  -12,  -22,  -16,    2,   28,   21,   -4,  -27},
{ -28,   -6,   24,   32,   -9,  -27,  -16,   15,   56,    6,  -48,  -49,   26,   48,   11,  -39},
{  49,   -2,  -52,  -36,   39,   46,   -7,  -52,   22,   20,  -12,  -31,    0,   29,   14,  -23},
{ -27,   14,   30,  -10,  -30,   12,   27,  -13,  -45,   33,   44,  -32,  -46,   36,   42,  -38},
{  41,  -37,  -39,   41,   34,  -46,  -30,   49,  -25,    2,   26,   -3,  -31,   13,   26,  -20},
{ -19,   25,    8,  -27,   18,    6,  -20,    6,   40,  -53,  -12,   57,  -49,    1,   51,  -31},
{  35,  -51,   -7,   52,  -55,   10,   51,  -40,   22,  -20,  -12,   27,  -16,   -4,   18,  -10},
{ -11,   24,  -21,    6,   16,  -31,   28,  -10,  -24,   55,  -53,   20,   22,  -52,   52,  -21},
{ -25,   55,  -53,   23,   18,  -48,   52,  -25,   13,  -22,   16,   -2,  -17,   33,  -30,   12},
{  14,  -32,   39,  -30,   26,  -29,   20,   -6,  -14,   35,  -47,   40,  -40,   48,  -37,   14},
{  15,  -33,   43,  -39,   41,  -49,   41,  -22,   18,  -32,   37,  -30,   26,  -28,   20,   -8},
},
{
{ -27,  -30,  -32,  -33,  -33,  -34,  -34,  -34,  -34,  -33,  -33,  -33,  -32,  -31,  -30,  -28},
{  47,   46,   41,   35,   28,   20,   11,    4,   -8,  -14,  -22,  -29,  -35,  -39,  -42,  -42},
{  53,   42,   22,    4,  -15,  -28,  -39,  -45,  -37,  -29,  -15,   -3,   15,   27,   38,   42},
{ -43,  -22,    5,   24,   37,   37,   30,   17,  -37,  -45,  -43,  -34,  -13,    9,   32,   44},
{  54,   11,  -34,  -53,  -33,    0,   34,   51,   22,    2,  -20,  -35,  -30,   -9,   17,   32},
{ -29,    2,   27,   28,   -1,  -22,  -30,  -20,   51,   42,    7,  -31,  -56,  -31,   16,   51},
{  33,  -16,  -46,  -22,   51,   46,   -4,  -49,  -24,   12,   35,   25,  -31,  -32,   -1,   27},
{ -24,   22,   32,   -7,  -38,  -15,   24,   31,  -49,  -19,   39,   51,  -20,  -48,  -14,   34},
{  35,  -42,  -28,   48,   20,  -39,  -25,   34,   20,  -27,  -23,   27,   28,  -34,  -30,   36},
{ -26,   37,   14,  -46,    2,   34,    0,  -28,   41,  -11,  -43,   11,   49,  -26,  -42,   38},
{  20,  -36,    6,   35,  -50,   15,   46,  -36,  -27,   39,    9,  -41,   35,   -3,  -37,   26},
{ -15,   28,  -10,  -22,   46,  -24,  -33,   41,  -45,   39,   25,  -51,   28,    9,  -35,   21},
{   7,  -18,   18,   -7,   -9,   20,  -18,    5,   22,  -45,   36,   -4,  -38,   69,  -64,   27},
{ -23,   59,  -65,   35,    9,  -42,   51,  -27,   10,   -9,    0,    7,  -17,   25,  -22,    8},
{ -11,   29,  -39,   33,  -35,   45,  -34,   11,   17,  -43,   52,  -37,   28,  -29,   21,   -8},
{   7,  -20,   30,  -30,   37,  -48,   44,  -25,   26,  -47,   49,  -36,   25,  -22,   15,   -6},
},
{
{ -25,  -28,  -31,  -32,  -33,  -34,  -35,  -36,  -35,  -34,  -34,  -33,  -32,  -31,  -29,  -27},
{  49,   47,   42,   36,   27,   19,   10,    3,   -8,  -14,  -22,  -28,  -34,  -38,  -41,  -40},
{  50,   41,   24,    6,  -15,  -28,  -38,  -43,  -37,  -28,  -15,   -2,   16,   30,   41,   46},
{ -38,  -23,    0,   20,   38,   40,   34,   20,  -34,  -45,  -45,  -37,  -12,   10,   31,   41},
{ -54,  -16,   32,   56,   37,    2,  -33,  -50,  -24,   -5,   17,   32,   30,   11,  -15,  -32},
{  31,    2,  -27,  -32,    0,   24,   31,   17,  -52,  -41,   -4,   36,   53,   27,  -17,  -49},
{ -30,   12,   42,   20,  -45,  -43,    3,   47,   27,  -13,  -41,  -27,   34,   38,    1,  -35},
{  33,  -20,  -41,    1,   43,   16,  -26,  -31,   47,   19,  -33,  -45,   16,   45,   13,  -35},
{ -34,   39,   31,  -45,  -25,   42,   28,  -42,  -19,   31,   20,  -32,  -22,   33,   24,  -30},
{  25,  -32,  -13,   40,   -2,  -34,    4,   31,  -48,   15,   47,  -18,  -47,   29,   38,  -36},
{ -26,   42,   -1,  -41,   51,  -11,  -46,   36,   17,  -34,    2,   33,  -38,    7,   35,  -25},
{  15,  -28,    7,   23,  -40,   20,   29,  -36,   44,  -41,  -20,   53,  -38,   -3,   41,  -25},
{ -22,   52,  -53,   24,   17,  -49,   52,  -23,  -11,   31,  -33,   15,    8,  -24,   26,  -11},
{  13,  -33,   37,  -22,    3,   13,  -22,   16,  -26,   50,  -42,   11,   27,  -56,   55,  -24},
{  13,  -33,   44,  -37,   38,  -46,   36,  -14,   -7,   24,  -34,   30,  -32,   41,  -34,   14},
{   7,  -19,   27,  -27,   33,  -42,   39,  -23,   25,  -44,   49,  -39,   35,  -36,   27,  -11},
},
{
{ -21,  -24,  -27,  -30,  -32,  -34,  -35,  -36,  -36,  -36,  -35,  -35,  -34,  -33,  -31,  -29},
{  42,   43,   42,   38,   32,   25,   16,    9,   -1,   -9,  -18,  -26,  -33,  -39,  -43,  -44},
{  45,   42,   31,   17,   -4,  -18,  -31,  -39,  -44,  -37,  -23,   -8,   11,   26,   40,   47},
{ -40,  -31,  -10,   12,   36,   44,   42,   31,  -28,  -40,  -43,  -36,  -16,    5,   25,   36},
{  47,   23,  -17,  -46,  -46,  -15,   23,   48,   34,   12,  -16,  -37,  -36,  -14,   16,   37},
{ -35,  -11,   22,   39,   16,  -14,  -35,  -33,   46,   39,    8,  -27,  -50,  -28,   15,   47},
{  30,    0,  -34,  -33,   30,   38,    7,  -31,  -39,    5,   46,   40,  -38,  -44,   -4,   35},
{ -39,    9,   47,   21,  -45,  -34,   18,   45,  -39,  -19,   27,   41,  -13,  -37,  -12,   26},
{  34,  -28,  -35,   30,   33,  -32,  -32,   33,   31,  -32,  -31,   34,   30,  -34,  -30,   33},
{ -32,   34,   27,  -41,  -17,   39,   13,  -39,   39,  -14,  -39,   17,   42,  -27,  -36,   33},
{  30,  -39,   -9,   42,  -40,    6,   39,  -27,  -28,   39,    8,  -40,   41,   -8,  -40,   30},
{ -24,   35,    0,  -35,   46,  -16,  -37,   38,  -38,   37,   17,  -46,   34,    1,  -36,   23},
{  19,  -42,   39,  -15,  -20,   47,  -47,   19,   17,  -43,   43,  -19,  -13,   35,  -37,   16},
{ -19,   44,  -45,   22,    7,  -29,   36,  -20,   22,  -41,   34,   -9,  -24,   49,  -47,   20},
{ -15,   36,  -43,   32,  -30,   37,  -28,    9,   10,  -31,   41,  -33,   35,  -47,   38,  -14},
{  11,  -28,   36,  -33,   37,  -47,   41,  -22,   21,  -39,   43,  -34,   30,  -32,   24,   -8},
},
};

const int KLTRow16[4][16][16]=
{
{
{ -23,   41,   40,  -47,  -44,   32,  -41,   37,  -38,   22,  -37,    6,  -28,    2,  -14,   -6},
{ -27,   44,   39,  -34,  -26,    8,    6,  -17,   35,  -26,   50,   -9,   66,   -5,   31,   14},
{ -29,   43,   29,   -5,   11,  -25,   46,  -44,   31,   -9,    6,    1,  -74,    5,  -33,  -17},
{ -31,   40,   16,   23,   41,  -38,   33,   -4,  -35,   29,  -61,   11,   52,    1,   12,    5},
{ -31,   34,    0,   41,   45,  -11,  -34,   44,  -30,   -3,   61,  -18,  -33,  -11,   26,   23},
{ -33,   24,  -20,   43,   15,   30,  -44,   11,   43,  -24,  -13,    8,   30,   18,  -54,  -50},
{ -34,   13,  -39,   32,  -24,   47,   -4,  -32,   30,    4,  -39,   13,  -25,  -17,   49,   52},
{ -36,    5,  -47,   12,  -42,   25,   39,  -24,  -54,   37,   40,  -32,   13,   13,  -24,  -22},
{ -35,   -2,  -45,  -18,  -31,  -41,   24,   42,  -12,  -56,   -2,   56,   -5,  -20,   18,  -24},
{ -35,  -11,  -36,  -35,   -6,  -50,  -15,   24,   34,   16,  -11,  -42,    4,   36,  -34,   58},
{ -35,  -21,  -18,  -42,   27,  -18,  -40,  -30,   16,   41,   -2,  -19,   -5,  -40,   40,  -60},
{ -34,  -30,    0,  -36,   45,   24,  -21,  -44,  -35,  -14,   19,   61,    4,   37,  -18,   29},
{ -34,  -36,   16,  -18,   33,   44,   36,   21,  -11,  -47,  -23,  -53,    2,  -49,  -19,    9},
{ -33,  -40,   30,    9,    5,   23,   37,   45,   29,   26,    3,    7,   -9,   67,   44,  -26},
{ -31,  -41,   40,   32,  -26,  -15,    0,   10,   19,   46,   22,   41,    9,  -58,  -39,   21},
{ -27,  -37,   39,   40,  -41,  -38,  -36,  -38,  -30,  -45,  -17,  -32,   -4,   24,   16,   -8},
},
{
{ -23,  -34,   37,   37,  -38,   41,  -38,   45,  -38,   23,  -32,   32,   -6,    1,   34,   13},
{ -26,  -39,   40,   35,  -29,   22,   -7,   -5,   23,  -21,   37,  -40,   10,   -2,  -72,  -28},
{ -28,  -41,   34,   19,    1,  -18,   34,  -49,   39,  -18,   15,  -10,   -7,    5,   74,   29},
{ -30,  -39,   23,   -4,   32,  -46,   37,  -22,  -21,   23,  -48,   51,    3,  -14,  -48,   -9},
{ -31,  -35,    5,  -31,   49,  -31,  -13,   44,  -40,   16,   26,  -55,   -2,   30,   23,  -30},
{ -33,  -28,  -15,  -45,   29,   12,  -38,   30,   31,  -35,   14,   20,    5,  -45,  -16,   62},
{ -34,  -20,  -31,  -45,   -7,   43,  -16,  -25,   37,   -9,  -20,   34,   -4,   45,   12,  -63},
{ -35,  -12,  -40,  -29,  -35,   36,   24,  -40,  -35,   46,  -12,  -50,    0,  -31,   -6,   28},
{ -35,   -2,  -47,   11,  -45,  -24,   42,   24,  -34,  -31,   51,   41,    3,   37,   -8,   17},
{ -35,    8,  -42,   33,  -23,  -42,    0,   35,   32,  -11,  -34,  -20,   -5,  -60,   21,  -44},
{ -35,   18,  -29,   45,   15,  -23,  -42,  -12,   30,   30,  -28,  -18,   10,   62,  -22,   45},
{ -35,   28,  -12,   41,   44,   16,  -34,  -44,  -29,   15,   54,   30,  -24,  -35,    9,  -23},
{ -34,   36,   10,   19,   41,   44,   34,    5,  -27,  -57,  -35,  -17,   56,    5,    7,   -1},
{ -34,   42,   27,   -7,   13,   27,   47,   38,   28,   19,   -2,   -4,  -82,   10,  -14,   10},
{ -33,   44,   39,  -29,  -21,  -12,    4,   14,   28,   55,   31,   18,   70,   -9,   10,   -8},
{ -30,   42,   41,  -38,  -40,  -40,  -43,  -32,  -32,  -46,  -19,  -11,  -28,    4,   -3,    2},
},
{
{ -23,   43,   45,  -42,  -43,   36,  -38,  -37,   38,  -28,   28,  -15,  -23,   21,   -2,   -4},
{ -28,   46,   41,  -30,  -22,    8,    6,   18,  -36,   35,  -38,   23,   54,  -49,    5,   10},
{ -31,   43,   28,   -3,   15,  -28,   43,   43,  -31,   11,   -4,   -3,  -60,   56,   -5,  -12},
{ -33,   38,   11,   24,   43,  -40,   29,    3,   37,  -41,   48,  -26,   37,  -37,   -7,    1},
{ -34,   30,   -8,   41,   43,   -8,  -35,  -44,   27,    9,  -50,   45,  -15,   14,   32,   27},
{ -35,   20,  -26,   42,   11,   30,  -41,  -15,  -40,   28,   12,  -27,    8,   -1,  -55,  -56},
{ -36,    9,  -39,   30,  -26,   44,   -2,   32,  -29,   -9,   33,  -25,   -5,   -2,   55,   57},
{ -36,    1,  -44,   10,  -45,   23,   36,   31,   43,  -34,  -24,   53,    1,   -1,  -35,  -24},
{ -36,   -6,  -41,  -24,  -35,  -39,   28,  -37,   25,   43,  -24,  -56,    5,    3,   33,  -23},
{ -35,  -14,  -31,  -40,   -8,  -47,  -13,  -30,  -35,   -3,   29,   29,  -11,   -2,  -52,   56},
{ -34,  -23,  -15,  -45,   26,  -17,  -41,   24,  -23,  -37,   11,   27,   12,    4,   54,  -58},
{ -33,  -32,    3,  -36,   45,   25,  -24,   46,   33,    2,  -43,  -48,  -17,  -16,  -31,   29},
{ -32,  -37,   19,  -14,   34,   45,   38,  -15,   20,   51,   43,   29,   35,   41,    4,    3},
{ -30,  -40,   32,   12,    5,   24,   41,  -45,  -29,  -22,   -7,    2,  -55,  -60,   11,  -16},
{ -27,  -41,   41,   34,  -25,  -15,   -1,  -12,  -24,  -48,  -36,  -25,   49,   52,  -10,   13},
{ -24,  -37,   41,   41,  -39,  -41,  -40,   37,   31,   44,   29,   16,  -20,  -21,    4,   -5},
},
{
{ -19,   37,   41,  -42,   45,  -39,   35,   39,   38,  -25,   33,  -25,   -8,   32,   -4,   -6},
{ -23,   42,   42,  -34,   28,  -13,    1,  -10,  -29,   27,  -39,   32,   18,  -73,   10,   15},
{ -27,   42,   34,  -12,   -8,   24,  -36,  -46,  -36,   16,  -15,    5,  -16,   78,  -10,  -17},
{ -30,   40,   21,   13,  -39,   42,  -35,  -17,   28,  -31,   51,  -41,    6,  -49,   -3,    3},
{ -32,   34,    4,   36,  -47,   17,   26,   45,   37,  -10,  -37,   52,   -1,   22,   25,   31},
{ -34,   26,  -16,   44,  -20,  -22,   38,   27,  -36,   33,    1,  -23,    3,  -14,  -45,  -64},
{ -36,   17,  -33,   39,   17,  -44,    9,  -26,  -36,    6,   24,  -31,   -3,   11,   45,   65},
{ -37,    8,  -41,   22,   41,  -31,  -28,  -42,   40,  -42,   -6,   50,   -1,   -7,  -28,  -28},
{ -37,    1,  -45,  -17,   38,   35,  -35,   32,   31,   38,  -41,  -47,    7,    3,   34,  -18},
{ -36,   -8,  -37,  -36,   16,   46,    5,   32,  -33,   -1,   32,   27,  -14,    1,  -61,   47},
{ -36,  -19,  -22,  -44,  -20,   20,   43,  -18,  -28,  -33,   23,   21,   17,    0,   63,  -48},
{ -35,  -29,   -4,  -39,  -43,  -20,   36,  -45,   30,   -2,  -50,  -40,  -24,   -5,  -33,   23},
{ -34,  -35,   13,  -20,  -37,  -45,  -39,   10,   23,   54,   39,   24,   51,   14,   -2,    3},
{ -32,  -41,   29,    6,   -9,  -26,  -45,   39,  -26,  -24,   -6,    3,  -79,  -17,   20,  -13},
{ -30,  -43,   40,   30,   23,   13,   -3,   13,  -26,  -51,  -31,  -23,   70,   13,  -17,   10},
{ -27,  -41,   43,   40,   40,   41,   40,  -30,   30,   47,   25,   15,  -28,   -5,    7,   -4},
},
};

static int KLTCol16inProd[4][16], KLTRow16inProd[4][16];

/*! 
*************************************************************************************
* \brief
*   Compute inner product for the input 16x16 matrix 
*
* \para compute_inner_product16x16()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
void compute_inner_product16x16(int matrix[16][16], int isLeft, int *inner_prod)
{
  int i, j;

  if(isLeft)
    for(i = 0; i < 16; i++)
    {
      inner_prod[i] = 0;
      for(j = 0; j < 8; j++)
      {
        inner_prod[i] += matrix[i][2*j+1]*matrix[i][2*j];
      }
    }
  else
    for(i = 0; i < 16; i++)
    {
      inner_prod[i] = 0;
      for(j = 0; j < 8; j++)
      {
        inner_prod[i] += matrix[2*j+1][i]*matrix[2*j][i];
      }
    }
}

/*! 
*************************************************************************************
* \brief
*   Precompute inner products for all 16x16 transform matrices
*
* \para precompute_all_inner_product16x16()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
void precompute_all_inner_product16x16()
{
  int i, j, k;
  int temp2D[16][16];


  for(i = 0; i < 4; i++)
  {
    for(j = 0; j < 16; j++)
      for(k = 0; k < 16; k++)
        temp2D[j][k] = KLTCol16[i][j][k];
    compute_inner_product16x16(temp2D, 1, KLTCol16inProd[i]);

    for(j = 0; j < 16; j++)
      for(k = 0; k < 16; k++)
        temp2D[j][k] = KLTRow16[i][j][k];
    compute_inner_product16x16(temp2D, 0, KLTRow16inProd[i]);
  }
}

/*! 
*************************************************************************************
* \brief
*   Perform fast mode-dependent 16x16 directional transform 
*
* \para klt_luma16x16_sep_fast()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
int klt_luma16x16_sep_fast(int ipmode)
{
  int sign(int a,int b);

  int i,j,k,coeff_ctr;
  int scan_pos,run;
  int nonzero;
  int qp_per,qp_rem, q_bits;

  int*  ACLevel = img->cofAC16[0];
  int*  ACRun   = img->cofAC16[1];
  byte (*pos_scan)[2] = img->scanOrder16x16[ipmode];
  int level;


  int coeff [MB_BLOCK_SIZE][MB_BLOCK_SIZE], ilev[MB_BLOCK_SIZE][MB_BLOCK_SIZE], 
      temp2D[MB_BLOCK_SIZE][MB_BLOCK_SIZE], y2D [MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  int M1[MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  int inprod[MB_BLOCK_SIZE];

  int m=4, n=10, KLTprec=14;
  int qpOffset;
  int roundOffset;
  int R[6] = {40, 45, 50, 57, 63, 71};
  int Q[6] = {410,   364,   328,   287,   260,   231};
  long shift;
  int shiftBy2;

  static int first=1;

  //Macroblock *currMB = &img->mb_data[img->current_mb_nr];

#ifdef RDO_Q
  levelDataStruct levelData[256];
  double lambda_md = 0.0;
  double err;
  int lowerInt;
  int levelRDOQ[256];
  int kStart = 0;
  int kStop = 0;
  int noCoeff;
  double normFact = pow(2, 2 * (KLTprec + m + n) - 15);
#endif

  pos_scan = img->scanOrder16x16[ipmode];

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
    img->delta_shift16x16=1<<(n+qp_per+KLTprec-11); 
    img->QPVal16x16=1<<(n+qp_per+KLTprec-3);  
    for (i=0; i<256; i++)
    {
      img->qp_shift16x16[i]=(1<<(n+qp_per+KLTprec-2))/3;
    }
  }

  for (j=0;j<16;j++)
  {
    for (i=0;i<16;i++)
    {
      M1[j][i]=imgY_org[img->opix_y+j][img->opix_x+i]-img->mprr_2[ipmode][j][i];
    }
  }

  compute_inner_product16x16(M1, 0, inprod);

  // column transform first
  for (i=0;i<(MB_BLOCK_SIZE);i++)
  {
    for (j=0;j<(MB_BLOCK_SIZE);j++)
    {
      temp2D[i][j] = 0;
      for(k = 0; k < MB_BLOCK_SIZE/2; k++)
        temp2D[i][j]+=
        (KLTCol16[ipmode][i][2*k]+M1[2*k+1][j])*(KLTCol16[ipmode][i][2*k+1]+M1[2*k][j]);
      temp2D[i][j] = temp2D[i][j] - KLTCol16inProd[ipmode][i] - inprod[j];
    }
  }

  compute_inner_product16x16(temp2D, 1, inprod);

  // row transform second
  for (i=0;i<(MB_BLOCK_SIZE);i++)
  {
    for (j=0;j<(MB_BLOCK_SIZE);j++)
    {
      coeff[i][j] = 0;
      for(k = 0; k < MB_BLOCK_SIZE/2; k++)
        coeff[i][j]+=
        (temp2D[i][2*k+1]+KLTRow16[ipmode][2*k][j])*(temp2D[i][2*k]+KLTRow16[ipmode][2*k+1][j]);
      coeff[i][j] = coeff[i][j] - inprod[i] - KLTRow16inProd[ipmode][j];
    }
  }

#ifdef RDO_Q
  if(input->UseRDO_Q)
  {
    noCoeff = 0;
    for (coeff_ctr = 0; coeff_ctr < 256; coeff_ctr++)
    {
      i = pos_scan[coeff_ctr][0];
      j = pos_scan[coeff_ctr][1];

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

    est_writeRunLevel_CABAC(levelData, levelRDOQ, LUMA_16x16, lambda_md, kStart, kStop, noCoeff, 0);
  }
#endif

  // quantization and zig-zag scanning 
  nonzero=FALSE;

  run=-1;
  scan_pos=0;

  for (coeff_ctr=0;coeff_ctr < 256;coeff_ctr++)
  {
    i=pos_scan[coeff_ctr][0];
    j=pos_scan[coeff_ctr][1];
    
    run++;

#ifdef RDO_Q
    if(input->UseRDO_Q)
    {
      level = levelRDOQ[coeff_ctr];
    }
    else
#endif
    if (img->AdaptiveRounding)
    {
      level = (int)((absm((int64)coeff[j][i]) * Q[qp_rem] + img->qp_shift16x16[coeff_ctr]) >> qpOffset);
      if (level != 0)
      {
        quant_stat_rd16x16[coeff_ctr] = (long)(absm((int64)coeff[j][i]) * Q[qp_rem] - ((int64)level << qpOffset));
      }
      else if (2 * absm((int64)coeff[j][i]) * Q[qp_rem] > ((int64)1 << qpOffset))
      {
        quant_stat_rd16x16[coeff_ctr] = (long)(2 * absm((int64)coeff[j][i]) * Q[qp_rem] - ((int64)1 << qpOffset));
      }
    }
    else
    {
      level = (int)((absm((int64)coeff[j][i]) * Q[qp_rem] + roundOffset) >> qpOffset);
    }

    if (level != 0)
    {
      nonzero=TRUE;

      ACLevel[scan_pos] = sign(level, coeff[j][i]);
      ACRun  [scan_pos] = run;
      ++scan_pos;
      run=-1;                     // reset zero level counter

      ilev[j][i] = sign(level, coeff[j][i]) * (R[qp_rem] << qp_per); 
    }
    else 
      ilev[j][i]=0;
  }
  ACLevel[scan_pos] = 0;

  // inverse transform, column first
  for (i=0;i<(MB_BLOCK_SIZE);i++)
  {
    for (j=0;j<(MB_BLOCK_SIZE);j++)
    {
      temp2D[i][j] = 0;
      for(k = 0; k < MB_BLOCK_SIZE; k++)
        temp2D[i][j]+=KLTCol16[ipmode][k][i]*(ilev[k][j]);
    }
  }

  // inverse transform, row second 
  for (i=0;i<(MB_BLOCK_SIZE);i++)
  {
    for (j=0;j<(MB_BLOCK_SIZE);j++)
    {
      y2D[i][j] = 0;
      for(k = 0; k < MB_BLOCK_SIZE; k++)
        y2D[i][j]+=temp2D[i][k]*KLTRow16[ipmode][j][k];
    }
  }

  // store the encoded image
  shift=KLTprec+m+2; shiftBy2=1<<(shift-1);
  for (j=0; j < MB_BLOCK_SIZE; j++)
    for (i=0; i < MB_BLOCK_SIZE; i++)
    {
      short temp;

      temp=(short)((y2D[j][i]+((int)img->mprr_2[ipmode][j][i]<<shift)+shiftBy2)>>shift);
      temp=min(img->max_imgpel_value,max(0,temp));
      enc_picture->imgY[img->pix_y+j][img->pix_x+i]=temp;
    }

  return nonzero;
}
#endif
