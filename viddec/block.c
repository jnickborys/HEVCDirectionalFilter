
/*!
 ***********************************************************************
 *  \file
 *      block.c
 *
 *  \brief
 *      Block functions
 *
 *  \author
 *      Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Inge Lille-Langoy          <inge.lille-langoy@telenor.com>
 *      - Rickard Sjoberg            <rickard.sjoberg@era.ericsson.se>
 ***********************************************************************
 */

#include "contributors.h"

#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "block.h"
#include "image.h"
#include "mb_access.h"
#define Q_BITS          15

#ifdef ADAPTIVE_QUANTIZATION
extern void AssignQuantParamForAQMS(pic_parameter_set_rbsp_t* pps);
extern void CalculateQuantParamForAQMS();
#endif

static const int quant_coef[6][4][4] = {
  {{13107, 8066,13107, 8066},{ 8066, 5243, 8066, 5243},{13107, 8066,13107, 8066},{ 8066, 5243, 8066, 5243}},
  {{11916, 7490,11916, 7490},{ 7490, 4660, 7490, 4660},{11916, 7490,11916, 7490},{ 7490, 4660, 7490, 4660}},
  {{10082, 6554,10082, 6554},{ 6554, 4194, 6554, 4194},{10082, 6554,10082, 6554},{ 6554, 4194, 6554, 4194}},
  {{ 9362, 5825, 9362, 5825},{ 5825, 3647, 5825, 3647},{ 9362, 5825, 9362, 5825},{ 5825, 3647, 5825, 3647}},
  {{ 8192, 5243, 8192, 5243},{ 5243, 3355, 5243, 3355},{ 8192, 5243, 8192, 5243},{ 5243, 3355, 5243, 3355}},
  {{ 7282, 4559, 7282, 4559},{ 4559, 2893, 4559, 2893},{ 7282, 4559, 7282, 4559},{ 4559, 2893, 4559, 2893}}
};
static const int A[4][4] = {
  { 16, 20, 16, 20},
  { 20, 25, 20, 25},
  { 16, 20, 16, 20},
  { 20, 25, 20, 25}
};

int quant_intra_default[16] = {
 6,13,20,28,
13,20,28,32,
20,28,32,37,
28,32,37,42
};

int quant_inter_default[16] = {
10,14,20,24,
14,20,24,27,
20,24,27,30,
24,27,30,34
};

int quant8_intra_default[64] = {
 6,10,13,16,18,23,25,27,
10,11,16,18,23,25,27,29,
13,16,18,23,25,27,29,31,
16,18,23,25,27,29,31,33,
18,23,25,27,29,31,33,36,
23,25,27,29,31,33,36,38,
25,27,29,31,33,36,38,40,
27,29,31,33,36,38,40,42
};

int quant8_inter_default[64] = {
 9,13,15,17,19,21,22,24,
13,13,17,19,21,22,24,25,
15,17,19,21,22,24,25,27,
17,19,21,22,24,25,27,28,
19,21,22,24,25,27,28,30,
21,22,24,25,27,28,30,32,
22,24,25,27,28,30,32,33,
24,25,27,28,30,32,33,35
};

int quant_org[16] = { //to be use if no q matrix is chosen
16,16,16,16,
16,16,16,16,
16,16,16,16,
16,16,16,16
};

int quant8_org[64] = { //to be use if no q matrix is chosen
16,16,16,16,16,16,16,16,
16,16,16,16,16,16,16,16,
16,16,16,16,16,16,16,16,
16,16,16,16,16,16,16,16,
16,16,16,16,16,16,16,16,
16,16,16,16,16,16,16,16,
16,16,16,16,16,16,16,16,
16,16,16,16,16,16,16,16
};

// Notation for comments regarding prediction and predictors.
// The pels of the 4x4 block are labelled a..p. The predictor pels above
// are labelled A..H, from the left I..L, and from above left X, as follows:
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
 ***********************************************************************
 * \brief
 *    makes and returns 4x4 blocks with all 5 intra prediction modes
 *
 * \return
 *    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
 *    SEARCH_SYNC   search next sync element as errors while decoding occured
 ***********************************************************************
 */

int intrapred( struct img_par *img,  //!< image parameters
               int ioff,             //!< pixel offset X within MB
               int joff,             //!< pixel offset Y within MB
               int img_block_x,      //!< location of block X, multiples of 4
               int img_block_y)      //!< location of block Y, multiples of 4
{
  int i,j;
  int s0;
  int img_y,img_x;
  int PredPel[13];  // array of predictor pels

  imgpel **imgY = dec_picture->imgY;

  PixelPos pix_a[4];
  PixelPos pix_b, pix_c, pix_d;

  int block_available_up;
  int block_available_left;
  int block_available_up_left;
  int block_available_up_right;

  int mb_nr=img->current_mb_nr;

  byte predmode = img->ipredmode[img_block_x][img_block_y];

  img_x=img_block_x*4;
  img_y=img_block_y*4;

  for (i=0;i<4;i++)
  {
    getNeighbour(mb_nr, ioff -1 , joff +i , 1, &pix_a[i]);
  }  
  
  getNeighbour(mb_nr, ioff    , joff -1 , 1, &pix_b);
  getNeighbour(mb_nr, ioff +4 , joff -1 , 1, &pix_c);
  getNeighbour(mb_nr, ioff -1 , joff -1 , 1, &pix_d);

  pix_c.available = pix_c.available && !((ioff==4) && ((joff==4)||(joff==12)));

  if (active_pps->constrained_intra_pred_flag)
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

  
  switch (predmode)
  {
  case DC_PRED:                         /* DC prediction */

    s0 = 0;
    if (block_available_up && block_available_left)
    {   
      // no edge
      s0 = (P_A + P_B + P_C + P_D + P_I + P_J + P_K + P_L + 4)/(2*BLOCK_SIZE);
    }
    else if (!block_available_up && block_available_left)
    {
      // upper edge
      s0 = (P_I + P_J + P_K + P_L + 2)/BLOCK_SIZE;             
    }
    else if (block_available_up && !block_available_left)
    {
      // left edge
      s0 = (P_A + P_B + P_C + P_D + 2)/BLOCK_SIZE;             
    }
    else //if (!block_available_up && !block_available_left)
    {
      // top left corner, nothing to predict from
      s0 = img->dc_pred_value_luma;
    }

    for (j=0; j < BLOCK_SIZE; j++)
    {
      for (i=0; i < BLOCK_SIZE; i++)
      {
        // store DC prediction
        img->mpr[i+ioff][j+joff] = s0;
      }
    }
    break;

  case VERT_PRED:                       /* vertical prediction from block above */
    if (!block_available_up)
      printf ("warning: Intra_4x4_Vertical prediction mode not allowed at mb %d\n",img->current_mb_nr);

    for(j=0;j<BLOCK_SIZE;j++)
      for(i=0;i<BLOCK_SIZE;i++)
        img->mpr[i+ioff][j+joff]=imgY[pix_b.pos_y][pix_b.pos_x+i];/* store predicted 4x4 block */
    break;

  case HOR_PRED:                        /* horizontal prediction from left block */
    if (!block_available_left)
      printf ("warning: Intra_4x4_Horizontal prediction mode not allowed at mb %d\n",img->current_mb_nr);

    for(j=0;j<BLOCK_SIZE;j++)
      for(i=0;i<BLOCK_SIZE;i++)
        img->mpr[i+ioff][j+joff]=imgY[pix_a[j].pos_y][pix_a[j].pos_x]; /* store predicted 4x4 block */
    break;

  case DIAG_DOWN_RIGHT_PRED:
    if ((!block_available_up)||(!block_available_left)||(!block_available_up_left))
      printf ("warning: Intra_4x4_Diagonal_Down_Right prediction mode not allowed at mb %d\n",img->current_mb_nr);

    img->mpr[0+ioff][3+joff] = (P_L + 2*P_K + P_J + 2) / 4; 
    img->mpr[0+ioff][2+joff] =
    img->mpr[1+ioff][3+joff] = (P_K + 2*P_J + P_I + 2) / 4; 
    img->mpr[0+ioff][1+joff] =
    img->mpr[1+ioff][2+joff] = 
    img->mpr[2+ioff][3+joff] = (P_J + 2*P_I + P_X + 2) / 4; 
    img->mpr[0+ioff][0+joff] =
    img->mpr[1+ioff][1+joff] =
    img->mpr[2+ioff][2+joff] =
    img->mpr[3+ioff][3+joff] = (P_I + 2*P_X + P_A + 2) / 4; 
    img->mpr[1+ioff][0+joff] =
    img->mpr[2+ioff][1+joff] =
    img->mpr[3+ioff][2+joff] = (P_X + 2*P_A + P_B + 2) / 4;
    img->mpr[2+ioff][0+joff] =
    img->mpr[3+ioff][1+joff] = (P_A + 2*P_B + P_C + 2) / 4;
    img->mpr[3+ioff][0+joff] = (P_B + 2*P_C + P_D + 2) / 4;
    break;

  case DIAG_DOWN_LEFT_PRED:
    if (!block_available_up)
      printf ("warning: Intra_4x4_Diagonal_Down_Left prediction mode not allowed at mb %d\n",img->current_mb_nr);

    img->mpr[0+ioff][0+joff] = (P_A + P_C + 2*(P_B) + 2) / 4;
    img->mpr[1+ioff][0+joff] = 
    img->mpr[0+ioff][1+joff] = (P_B + P_D + 2*(P_C) + 2) / 4;
    img->mpr[2+ioff][0+joff] =
    img->mpr[1+ioff][1+joff] =
    img->mpr[0+ioff][2+joff] = (P_C + P_E + 2*(P_D) + 2) / 4;
    img->mpr[3+ioff][0+joff] = 
    img->mpr[2+ioff][1+joff] = 
    img->mpr[1+ioff][2+joff] = 
    img->mpr[0+ioff][3+joff] = (P_D + P_F + 2*(P_E) + 2) / 4;
    img->mpr[3+ioff][1+joff] = 
    img->mpr[2+ioff][2+joff] = 
    img->mpr[1+ioff][3+joff] = (P_E + P_G + 2*(P_F) + 2) / 4;
    img->mpr[3+ioff][2+joff] = 
    img->mpr[2+ioff][3+joff] = (P_F + P_H + 2*(P_G) + 2) / 4;
    img->mpr[3+ioff][3+joff] = (P_G + 3*(P_H) + 2) / 4;
    break;

  case  VERT_RIGHT_PRED:/* diagonal prediction -22.5 deg to horizontal plane */
    if ((!block_available_up)||(!block_available_left)||(!block_available_up_left))
      printf ("warning: Intra_4x4_Vertical_Right prediction mode not allowed at mb %d\n",img->current_mb_nr);

    img->mpr[0+ioff][0+joff] = 
    img->mpr[1+ioff][2+joff] = (P_X + P_A + 1) / 2;
    img->mpr[1+ioff][0+joff] = 
    img->mpr[2+ioff][2+joff] = (P_A + P_B + 1) / 2;
    img->mpr[2+ioff][0+joff] = 
    img->mpr[3+ioff][2+joff] = (P_B + P_C + 1) / 2;
    img->mpr[3+ioff][0+joff] = (P_C + P_D + 1) / 2;
    img->mpr[0+ioff][1+joff] = 
    img->mpr[1+ioff][3+joff] = (P_I + 2*P_X + P_A + 2) / 4;
    img->mpr[1+ioff][1+joff] = 
    img->mpr[2+ioff][3+joff] = (P_X + 2*P_A + P_B + 2) / 4;
    img->mpr[2+ioff][1+joff] = 
    img->mpr[3+ioff][3+joff] = (P_A + 2*P_B + P_C + 2) / 4;
    img->mpr[3+ioff][1+joff] = (P_B + 2*P_C + P_D + 2) / 4;
    img->mpr[0+ioff][2+joff] = (P_X + 2*P_I + P_J + 2) / 4;
    img->mpr[0+ioff][3+joff] = (P_I + 2*P_J + P_K + 2) / 4;
    break;

  case  VERT_LEFT_PRED:/* diagonal prediction -22.5 deg to horizontal plane */
    if (!block_available_up)
      printf ("warning: Intra_4x4_Vertical_Left prediction mode not allowed at mb %d\n",img->current_mb_nr);
    
    img->mpr[0+ioff][0+joff] = (P_A + P_B + 1) / 2;
    img->mpr[1+ioff][0+joff] = 
    img->mpr[0+ioff][2+joff] = (P_B + P_C + 1) / 2;
    img->mpr[2+ioff][0+joff] = 
    img->mpr[1+ioff][2+joff] = (P_C + P_D + 1) / 2;
    img->mpr[3+ioff][0+joff] = 
    img->mpr[2+ioff][2+joff] = (P_D + P_E + 1) / 2;
    img->mpr[3+ioff][2+joff] = (P_E + P_F + 1) / 2;
    img->mpr[0+ioff][1+joff] = (P_A + 2*P_B + P_C + 2) / 4;
    img->mpr[1+ioff][1+joff] = 
    img->mpr[0+ioff][3+joff] = (P_B + 2*P_C + P_D + 2) / 4;
    img->mpr[2+ioff][1+joff] = 
    img->mpr[1+ioff][3+joff] = (P_C + 2*P_D + P_E + 2) / 4;
    img->mpr[3+ioff][1+joff] = 
    img->mpr[2+ioff][3+joff] = (P_D + 2*P_E + P_F + 2) / 4;
    img->mpr[3+ioff][3+joff] = (P_E + 2*P_F + P_G + 2) / 4;
    break;

  case  HOR_UP_PRED:/* diagonal prediction -22.5 deg to horizontal plane */
    if (!block_available_left)
      printf ("warning: Intra_4x4_Horizontal_Up prediction mode not allowed at mb %d\n",img->current_mb_nr);
    
    img->mpr[0+ioff][0+joff] = (P_I + P_J + 1) / 2;
    img->mpr[1+ioff][0+joff] = (P_I + 2*P_J + P_K + 2) / 4;
    img->mpr[2+ioff][0+joff] = 
    img->mpr[0+ioff][1+joff] = (P_J + P_K + 1) / 2;
    img->mpr[3+ioff][0+joff] = 
    img->mpr[1+ioff][1+joff] = (P_J + 2*P_K + P_L + 2) / 4;
    img->mpr[2+ioff][1+joff] = 
    img->mpr[0+ioff][2+joff] = (P_K + P_L + 1) / 2;
    img->mpr[3+ioff][1+joff] = 
    img->mpr[1+ioff][2+joff] = (P_K + 2*P_L + P_L + 2) / 4;
    img->mpr[3+ioff][2+joff] = 
    img->mpr[1+ioff][3+joff] = 
    img->mpr[0+ioff][3+joff] = 
    img->mpr[2+ioff][2+joff] = 
    img->mpr[2+ioff][3+joff] = 
    img->mpr[3+ioff][3+joff] = P_L;
    break;

  case  HOR_DOWN_PRED:/* diagonal prediction -22.5 deg to horizontal plane */
    if ((!block_available_up)||(!block_available_left)||(!block_available_up_left))
      printf ("warning: Intra_4x4_Horizontal_Down prediction mode not allowed at mb %d\n",img->current_mb_nr);

    img->mpr[0+ioff][0+joff] = 
    img->mpr[2+ioff][1+joff] = (P_X + P_I + 1) / 2;
    img->mpr[1+ioff][0+joff] = 
    img->mpr[3+ioff][1+joff] = (P_I + 2*P_X + P_A + 2) / 4;
    img->mpr[2+ioff][0+joff] = (P_X + 2*P_A + P_B + 2) / 4;
    img->mpr[3+ioff][0+joff] = (P_A + 2*P_B + P_C + 2) / 4;
    img->mpr[0+ioff][1+joff] = 
    img->mpr[2+ioff][2+joff] = (P_I + P_J + 1) / 2;
    img->mpr[1+ioff][1+joff] = 
    img->mpr[3+ioff][2+joff] = (P_X + 2*P_I + P_J + 2) / 4;
    img->mpr[0+ioff][2+joff] = 
    img->mpr[2+ioff][3+joff] = (P_J + P_K + 1) / 2;
    img->mpr[1+ioff][2+joff] = 
    img->mpr[3+ioff][3+joff] = (P_I + 2*P_J + P_K + 2) / 4;
    img->mpr[0+ioff][3+joff] = (P_K + P_L + 1) / 2;
    img->mpr[1+ioff][3+joff] = (P_J + 2*P_K + P_L + 2) / 4;
    break;

  default:
    printf("Error: illegal intra_4x4 prediction mode: %d\n",predmode);
    return SEARCH_SYNC;
    break;
  }

  return DECODING_OK;
}


/*!
 ***********************************************************************
 * \return
 *    best SAD
 ***********************************************************************
 */
int intrapred_luma_16x16(struct img_par *img, //!< image parameters
                         int predmode)        //!< prediction mode
{
  int s0=0,s1,s2;

  int i,j;

  int ih,iv;
  int ib,ic,iaa;

  imgpel **imgY=dec_picture->imgY;

  int mb_nr=img->current_mb_nr;

  PixelPos up;          //!< pixel position p(0,-1)
  PixelPos left[17];    //!< pixel positions p(-1, -1..15)

  int up_avail, left_avail, left_up_avail;

  s1=s2=0;

  for (i=0;i<17;i++)
  {
    getNeighbour(mb_nr, -1 ,  i-1 , 1, &left[i]);
  }
  
  getNeighbour(mb_nr, 0     ,  -1 , 1, &up);

  if (!active_pps->constrained_intra_pred_flag)
  {
    up_avail   = up.available;
    left_avail = left[1].available;
    left_up_avail = left[0].available;
  }
  else
  {
    up_avail      = up.available ? img->intra_block[up.mb_addr] : 0;
    for (i=1, left_avail=1; i<17;i++)
      left_avail  &= left[i].available ? img->intra_block[left[i].mb_addr]: 0;
    left_up_avail = left[0].available ? img->intra_block[left[0].mb_addr]: 0;
  }

  switch (predmode)
  {
  case VERT_PRED_16:                       // vertical prediction from block above
    if (!up_avail)
      error ("invalid 16x16 intra pred Mode VERT_PRED_16",500);
    for(j=0;j<MB_BLOCK_SIZE;j++)
      for(i=0;i<MB_BLOCK_SIZE;i++)
        img->mpr[i][j]=imgY[up.pos_y][up.pos_x+i];// store predicted 16x16 block
    break;

  case HOR_PRED_16:                        // horizontal prediction from left block
    if (!left_avail)
      error ("invalid 16x16 intra pred Mode HOR_PRED_16",500);
    for(j=0;j<MB_BLOCK_SIZE;j++)
      for(i=0;i<MB_BLOCK_SIZE;i++)
        img->mpr[i][j]=imgY[left[j+1].pos_y][left[j+1].pos_x]; // store predicted 16x16 block
    break;

  case DC_PRED_16:                         // DC prediction
    s1=s2=0;
    for (i=0; i < MB_BLOCK_SIZE; i++)
    {
      if (up_avail)
        s1 += imgY[up.pos_y][up.pos_x+i];    // sum hor pix
      if (left_avail)
        s2 += imgY[left[i+1].pos_y][left[i+1].pos_x];    // sum vert pix
    }
    if (up_avail && left_avail)
      s0=(s1+s2+16)>>5;       // no edge
    if (!up_avail && left_avail)
      s0=(s2+8)>>4;              // upper edge
    if (up_avail && !left_avail)
      s0=(s1+8)>>4;              // left edge
    if (!up_avail && !left_avail)
      s0=img->dc_pred_value_luma;                            // top left corner, nothing to predict from
    for(i=0;i<MB_BLOCK_SIZE;i++)
      for(j=0;j<MB_BLOCK_SIZE;j++)
      {
        img->mpr[i][j]=s0;
      }
    break;
  case PLANE_16:// 16 bit integer plan pred
    if (!up_avail || !left_up_avail  || !left_avail)
      error ("invalid 16x16 intra pred Mode PLANE_16",500);

    ih=0;
    iv=0;
    for (i=1;i<9;i++)
    {
      if (i<8)
        ih += i*(imgY[up.pos_y][up.pos_x+7+i] - imgY[up.pos_y][up.pos_x+7-i]);
      else
        ih += i*(imgY[up.pos_y][up.pos_x+7+i] - imgY[left[0].pos_y][left[0].pos_x]);

      iv += i*(imgY[left[8+i].pos_y][left[8+i].pos_x] - imgY[left[8-i].pos_y][left[8-i].pos_x]);
    }
    ib=(5*ih+32)>>6;
    ic=(5*iv+32)>>6;

    iaa=16*(imgY[up.pos_y][up.pos_x+15]+imgY[left[16].pos_y][left[16].pos_x]);
    for (j=0;j< MB_BLOCK_SIZE;j++)
    {
      for (i=0;i< MB_BLOCK_SIZE;i++)
      {
        img->mpr[i][j]=max(0,min((iaa+(i-7)*ib +(j-7)*ic + 16)>>5, img->max_imgpel_value));
      }
    }// store plane prediction
    break;
    
  default:
    {                                    // indication of fault in bitstream,exit
      printf("illegal 16x16 intra prediction mode input: %d\n",predmode);
      return SEARCH_SYNC;
    }
  }
  
  return DECODING_OK;
}


void intrapred_chroma(struct img_par *img, int uv)
{
  int i,j, ii, jj, ioff, joff;
  
  imgpel ***imgUV = dec_picture->imgUV;
  
  int js[4][4];
  
  int pred;
  int ih, iv, ib, ic, iaa;
  
  int      b8, b4;
  int      yuv = dec_picture->chroma_format_idc - 1;
  int      blk_x, blk_y;
  int      block_pos[3][4][4]= //[yuv][b8][b4]
  {
    { {0, 1, 2, 3},{0, 0, 0, 0},{0, 0, 0, 0},{0, 0, 0, 0}},
    { {0, 1, 2, 3},{2, 3, 2, 3},{0, 0, 0, 0},{0, 0, 0, 0}},
    { {0, 1, 2, 3},{1, 1, 3, 3},{2, 3, 2, 3},{3, 3, 3, 3}}
  };
  int      s0, s1, s2, s3;

  int mb_nr=img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  PixelPos up;        //!< pixel position  p(0,-1)
  PixelPos left[17];  //!< pixel positions p(-1, -1..16)

  int up_avail, left_avail[2], left_up_avail;

  int cr_MB_x = img->mb_cr_size_x;
  int cr_MB_y = img->mb_cr_size_y;

  for (i=0;i<cr_MB_y+1;i++)
  {
    getNeighbour(mb_nr, -1, i-1, 0, &left[i]);
  }
  
  getNeighbour(mb_nr, 0, -1, 0, &up);

  if (!active_pps->constrained_intra_pred_flag)
  {
    up_avail      = up.available;
    left_avail[0] = left_avail[1] = left[1].available;
    left_up_avail = left[0].available;
  }
  else
  {
    up_avail = up.available ? img->intra_block[up.mb_addr] : 0;
    for (i=0, left_avail[0]=1; i<cr_MB_y/2;i++)
      left_avail[0]  &= left[i+1].available ? img->intra_block[left[i+1].mb_addr]: 0;
    for (i=cr_MB_y/2, left_avail[1]=1; i<cr_MB_y;i++)
      left_avail[1]  &= left[i+1].available ? img->intra_block[left[i+1].mb_addr]: 0;
    left_up_avail = left[0].available ? img->intra_block[left[0].mb_addr]: 0;
  }

  
  if (currMB->c_ipred_mode == DC_PRED_8)
  {
    // DC prediction
    for(b8=0; b8<img->num_blk8x8_uv/2;b8++)
    {
      for (b4=0; b4<4; b4++)
      {
        blk_y = subblk_offset_y[yuv][b8][b4] + 1;
        blk_x = subblk_offset_x[yuv][b8][b4]; 
        
        s0=s1=s2=s3=0;
        js[b8][b4]=img->dc_pred_value_chroma;
        
        //===== get prediction value =====
        switch (block_pos[yuv][b8][b4])
        {
        case 0:  //===== TOP LEFT =====
          if      (up_avail)       for (i=blk_x;i<(blk_x+4);i++)  s0 += imgUV[uv][up.pos_y][up.pos_x + i];
          if      (left_avail[0])  for (i=blk_y;i<(blk_y+4);i++)  s2 += imgUV[uv][left[i].pos_y][left[i].pos_x];
          if      (up_avail && left_avail[0])         js[b8][b4]  = (s0+s2+4) >> 3;
          else if (up_avail)                          js[b8][b4]  = (s0   +2) >> 2;
          else if (left_avail[0])                     js[b8][b4]  = (s2   +2) >> 2;
          break;
        case 1: //===== TOP RIGHT =====
          if      (up_avail)       for (i=blk_x;i<(blk_x+4);i++)  s1 += imgUV[uv][up.pos_y][up.pos_x + i];
          else if (left_avail[0])  for (i=blk_y;i<(blk_y+4);i++)  s2 += imgUV[uv][left[i].pos_y][left[i].pos_x];
          if      (up_avail)                          js[b8][b4]  = (s1   +2) >> 2;
          else if (left_avail[0])                     js[b8][b4]  = (s2   +2) >> 2;
          break;
        case 2: //===== BOTTOM LEFT =====
          if      (left_avail[1])  for (i=blk_y;i<(blk_y+4);i++)  s3 += imgUV[uv][left[i].pos_y][left[i].pos_x];
          else if (up_avail)       for (i=blk_x;i<(blk_x+4);i++)  s0 += imgUV[uv][up.pos_y][up.pos_x + i];
          if      (left_avail[1])                     js[b8][b4]  = (s3   +2) >> 2;
          else if (up_avail)                          js[b8][b4]  = (s0   +2) >> 2;
          break;
        case 3: //===== BOTTOM RIGHT =====
          if      (up_avail)       for (i=blk_x;i<(blk_x+4);i++)  s1 += imgUV[uv][up.pos_y][up.pos_x + i];
          if      (left_avail[1])  for (i=blk_y;i<(blk_y+4);i++)  s3 += imgUV[uv][left[i].pos_y][left[i].pos_x];
          if      (up_avail && left_avail[1])         js[b8][b4]  = (s1+s3+4) >> 3;
          else if (up_avail)                          js[b8][b4]  = (s1   +2) >> 2;
          else if (left_avail[1])                     js[b8][b4]  = (s3   +2) >> 2;
          break;
        }
      }
    }
  }
  if (PLANE_8 == currMB->c_ipred_mode)
  {
    // plane prediction
    if (!left_up_avail || !left_avail[0] || !left_avail[1] || !up_avail)
      error("unexpected PLANE_8 chroma intra prediction mode",-1);
    
    ih = cr_MB_x/2*(imgUV[uv][up.pos_y][up.pos_x+cr_MB_x-1] - imgUV[uv][left[0].pos_y][left[0].pos_x]);
    for (i=0;i<cr_MB_x/2-1;i++)
      ih += (i+1)*(imgUV[uv][up.pos_y][up.pos_x+cr_MB_x/2  +i] -
      imgUV[uv][up.pos_y][up.pos_x+cr_MB_x/2-2-i]);
    
    iv = cr_MB_y/2*(imgUV[uv][left[cr_MB_y].pos_y][left[cr_MB_y].pos_x] - imgUV[uv][left[0].pos_y][left[0].pos_x]);
    for (i=0;i<cr_MB_y/2-1;i++)
      iv += (i+1)*(imgUV[uv][left[cr_MB_y/2+1+i].pos_y][left[cr_MB_y/2+1+i].pos_x] -
      imgUV[uv][left[cr_MB_y/2-1-i].pos_y][left[cr_MB_y/2-1-i].pos_x]);
    
    ib= ((cr_MB_x == 8?17:5)*ih+2*cr_MB_x)>>(cr_MB_x == 8?5:6);
    ic= ((cr_MB_y == 8?17:5)*iv+2*cr_MB_y)>>(cr_MB_y == 8?5:6);
    
    iaa=16*(imgUV[uv][left[cr_MB_y].pos_y][left[cr_MB_y].pos_x] +
            imgUV[uv][up.pos_y][up.pos_x+cr_MB_x-1]);
    
    for (j=0; j<cr_MB_y; j++)
      for (i=0; i<cr_MB_x; i++)
        img->mpr[i][j]=max(0,min(img->max_imgpel_value_uv,(iaa+(i-cr_MB_x/2+1)*ib+(j-cr_MB_y/2+1)*ic+16)>>5));
  }
  else
  {
    switch (currMB->c_ipred_mode)
    {
    case DC_PRED_8:
      for (b8=0;b8<img->num_blk8x8_uv/2;b8++)
      {
        for (b4=0;b4<4;b4++)
        {
          joff = subblk_offset_y[yuv][b8][b4];
          ioff = subblk_offset_x[yuv][b8][b4];
          for (ii=0; ii<BLOCK_SIZE; ii++)
          for (jj=0; jj<BLOCK_SIZE; jj++)
          {
            img->mpr[ii+ioff][jj+joff]=js[b8][b4];
          }
        }
      }
      break;
    case HOR_PRED_8:
      if (!left_avail[0] || !left_avail[1])
        error("unexpected HOR_PRED_8 chroma intra prediction mode",-1);

      for (j=0;j<2;j++)
      {
        joff=j*cr_MB_y/2;
        for(i=0;i<2;i++)
        {
          ioff=i*cr_MB_x/2;
          for (jj=0; jj<cr_MB_y/2; jj++)
          {
            pred = imgUV[uv][left[1+jj+joff].pos_y][left[1+jj+joff].pos_x];
            for (ii=0; ii<cr_MB_x/2; ii++)
              img->mpr[ii+ioff][jj+joff]=pred;
          }
        }
      }
      break;
    case VERT_PRED_8:
      if (!up_avail)
        error("unexpected VERT_PRED_8 chroma intra prediction mode",-1);

      for (j=0;j<2;j++)
      {
        joff=j*cr_MB_y/2;
        for(i=0;i<2;i++)
        {
          ioff=i*cr_MB_x/2;
          for (ii=0; ii<cr_MB_x/2; ii++)
          {
            pred = imgUV[uv][up.pos_y][up.pos_x+ii+ioff];
            for (jj=0; jj<cr_MB_y/2; jj++)
              img->mpr[ii+ioff][jj+joff]=pred;
          }
        }
      }
      break;
    default:
      error("illegal chroma intra prediction mode", 600);
      break;
    }
  }
}

/*!
 ***********************************************************************
 * \brief
 *    Inverse 4x4 transformation, transforms cof to m7
 ***********************************************************************
 */
void itrans(struct img_par *img, //!< image parameters
            int ioff,            //!< index to 4x4 block
            int joff,            //!<
            int i0,              //!<
            int j0,
            int chroma)
{
  int i,j,i1,j1;
  int m5[4];
  int m6[4];

  Boolean lossless_qpprime = (Boolean)((img->qp + img->bitdepth_luma_qp_scale)==0 && img->lossless_qpprime_flag==1);

  // Residue Color Transform
  int residue_transform_flag = img->residue_transform_flag;

  // horizontal
  for (j=0;j<BLOCK_SIZE && !lossless_qpprime;j++)
  {
    for (i=0;i<BLOCK_SIZE;i++)
    {
      m5[i]=img->cof[i0][j0][i][j];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (i=0;i<2;i++)
    {
      i1=3-i;
      img->m7[i][j]=m6[i]+m6[i1];
      img->m7[i1][j]=m6[i]-m6[i1];
    }
  }
  // vertical
  for (i=0;i<BLOCK_SIZE && !lossless_qpprime;i++)
  {
    for (j=0;j<BLOCK_SIZE;j++)
      m5[j]=img->m7[i][j];

    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (j=0;j<2;j++)
    {
      j1=3-j;
      // Residue Color Transform
      if(!residue_transform_flag)
      {
        if(!chroma)
        {
          img->m7[i][j] =max(0,min(img->max_imgpel_value,(m6[j]+m6[j1]+((long)img->mpr[i+ioff][j+joff] <<DQ_BITS)+DQ_ROUND)>>DQ_BITS));
          img->m7[i][j1]=max(0,min(img->max_imgpel_value,(m6[j]-m6[j1]+((long)img->mpr[i+ioff][j1+joff]<<DQ_BITS)+DQ_ROUND)>>DQ_BITS));
        } else 
        {
          img->m7[i][j] =max(0,min(img->max_imgpel_value_uv,(m6[j]+m6[j1]+((long)img->mpr[i+ioff][j+joff] <<DQ_BITS)+DQ_ROUND)>>DQ_BITS));
          img->m7[i][j1]=max(0,min(img->max_imgpel_value_uv,(m6[j]-m6[j1]+((long)img->mpr[i+ioff][j1+joff]<<DQ_BITS)+DQ_ROUND)>>DQ_BITS));
        }
      }
      else
      {
        img->m7[i][j] =(m6[j]+m6[j1]+DQ_ROUND)>>DQ_BITS;
        img->m7[i][j1]=(m6[j]-m6[j1]+DQ_ROUND)>>DQ_BITS;
      }
    }
  }

  // Residue Color Transform
  if(!residue_transform_flag)
  {
    for (i=0;i<BLOCK_SIZE && lossless_qpprime;i++)
      for (j=0;j<BLOCK_SIZE;j++)
        if(!chroma)
          img->m7[i][j] = max(0,min(img->max_imgpel_value,img->cof[i0][j0][i][j]+(long)img->mpr[i+ioff][j+joff]));
        else
          img->m7[i][j] = max(0,min(img->max_imgpel_value_uv,img->cof[i0][j0][i][j]+(long)img->mpr[i+ioff][j+joff]));
  }
  else
  {
    for (i=0;i<BLOCK_SIZE && lossless_qpprime;i++)
      for (j=0;j<BLOCK_SIZE;j++)
        img->m7[i][j] = img->cof[i0][j0][i][j];
  }
}

/*!
 ************************************************************************
 * \brief
 *    For mapping the q-matrix to the active id and calculate quantisation values
 *
 * \param pps
 *    Picture parameter set
 * \param sps
 *    Sequence parameter set
 *
 ************************************************************************
 */
void AssignQuantParam(pic_parameter_set_rbsp_t* pps, seq_parameter_set_rbsp_t* sps)
{
  int i;
  
  if(!pps->pic_scaling_matrix_present_flag && !sps->seq_scaling_matrix_present_flag)
  {
    for(i=0; i<8; i++)
      qmatrix[i] = (i<6) ? quant_org:quant8_org;
  }
  else
  {
    if(sps->seq_scaling_matrix_present_flag) // check sps first
    {
      for(i=0; i<8; i++)
      {
        if(i<6)
        {
          if(!sps->seq_scaling_list_present_flag[i]) // fall-back rule A
          {
            if((i==0) || (i==3))
              qmatrix[i] = (i==0) ? quant_intra_default:quant_inter_default;
            else
              qmatrix[i] = qmatrix[i-1];
          }
          else
          {
            if(sps->UseDefaultScalingMatrix4x4Flag[i])
              qmatrix[i] = (i<3) ? quant_intra_default:quant_inter_default;
            else
              qmatrix[i] = sps->ScalingList4x4[i];
          }
        }
        else
        {
          if(!sps->seq_scaling_list_present_flag[i] || sps->UseDefaultScalingMatrix8x8Flag[i-6]) // fall-back rule A
            qmatrix[i] = (i==6) ? quant8_intra_default:quant8_inter_default;
          else
            qmatrix[i] = sps->ScalingList8x8[i-6];
        }
      }
    }
    
    if(pps->pic_scaling_matrix_present_flag) // then check pps
    {
      for(i=0; i<8; i++)
      {
        if(i<6)
        {
          if(!pps->pic_scaling_list_present_flag[i]) // fall-back rule B
          {
            if((i==0) || (i==3))
            {
              if(!sps->seq_scaling_matrix_present_flag)
                qmatrix[i] = (i==0) ? quant_intra_default:quant_inter_default;
            }
            else
              qmatrix[i] = qmatrix[i-1];
          }
          else
          {
            if(pps->UseDefaultScalingMatrix4x4Flag[i])
              qmatrix[i] = (i<3) ? quant_intra_default:quant_inter_default;
            else
              qmatrix[i] = pps->ScalingList4x4[i];
          }
        }
        else
        {
          if(!pps->pic_scaling_list_present_flag[i]) // fall-back rule B
          {
            if(!sps->seq_scaling_matrix_present_flag)
              qmatrix[i] = (i==6) ? quant8_intra_default:quant8_inter_default;
          }
          else if(pps->UseDefaultScalingMatrix8x8Flag[i-6])
            qmatrix[i] = (i==6) ? quant8_intra_default:quant8_inter_default;
          else
            qmatrix[i] = pps->ScalingList8x8[i-6];
        }
      }
    }
  }
#ifdef ADAPTIVE_QUANTIZATION
  AssignQuantParamForAQMS(pps);
#endif
  CalculateQuantParam();
  if(pps->transform_8x8_mode_flag)
    CalculateQuant8Param();
}

/*!
 ************************************************************************
 * \brief
 *    For calculating the quantisation values at frame level
 *
 ************************************************************************
 */
void CalculateQuantParam()
{
  int i, j, k, temp;

  for(k=0; k<6; k++)
    for(j=0; j<4; j++)
      for(i=0; i<4; i++)
      {
        temp = (i<<2)+j;
        InvLevelScale4x4Luma_Intra[k][j][i]      = dequant_coef[k][j][i]*qmatrix[0][temp];
        InvLevelScale4x4Chroma_Intra[0][k][j][i] = dequant_coef[k][j][i]*qmatrix[1][temp];
        InvLevelScale4x4Chroma_Intra[1][k][j][i] = dequant_coef[k][j][i]*qmatrix[2][temp];

        InvLevelScale4x4Luma_Inter[k][j][i]      = dequant_coef[k][j][i]*qmatrix[3][temp];
        InvLevelScale4x4Chroma_Inter[0][k][j][i] = dequant_coef[k][j][i]*qmatrix[4][temp];
        InvLevelScale4x4Chroma_Inter[1][k][j][i] = dequant_coef[k][j][i]*qmatrix[5][temp];
      }
#ifdef ADAPTIVE_QUANTIZATION
  CalculateQuantParamForAQMS();
#endif
}

/*!
 ***********************************************************************
 * \brief
 *    Luma DC inverse transform
 ***********************************************************************
 */
void itrans_2(struct img_par *img) //!< image parameters
{
  int i,j,i1,j1;
  int M5[4];
  int M6[4];

  int qp_per = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)/6;
  int qp_rem = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)%6;

  int qp_const = 1<<(5-qp_per);
#ifdef ADAPTIVE_QUANTIZATION
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
#endif
  // horizontal
  for (j=0;j<4;j++)
  {
    for (i=0;i<4;i++)
      M5[i]=img->cof[i][j][0][0];

    M6[0]=M5[0]+M5[2];
    M6[1]=M5[0]-M5[2];
    M6[2]=M5[1]-M5[3];
    M6[3]=M5[1]+M5[3];

    for (i=0;i<2;i++)
    {
      i1=3-i;
      img->cof[i ][j][0][0]= M6[i]+M6[i1];
      img->cof[i1][j][0][0]=M6[i]-M6[i1];
    }
  }

  // vertical
  for (i=0;i<4;i++)
  {
    for (j=0;j<4;j++)
      M5[j]=img->cof[i][j][0][0];

    M6[0]=M5[0]+M5[2];
    M6[1]=M5[0]-M5[2];
    M6[2]=M5[1]-M5[3];
    M6[3]=M5[1]+M5[3];

    for (j=0;j<2;j++)
    {
      j1=3-j;
#ifdef ADAPTIVE_QUANTIZATION
      if(img->slice_fractional_quant_flag)
      {
        if(qp_per<6)
        {
          img->cof[i][j][0][0] =((M6[j]+M6[j1])*InvLevelScale4x4Luma_Intra_IAQMS[currMB->mb_iaqms_idx][qp_rem][0][0]+qp_const)>>(6-qp_per);
          img->cof[i][j1][0][0]=((M6[j]-M6[j1])*InvLevelScale4x4Luma_Intra_IAQMS[currMB->mb_iaqms_idx][qp_rem][0][0]+qp_const)>>(6-qp_per);
        }
        else
        {
          img->cof[i][j][0][0] =((M6[j]+M6[j1])*InvLevelScale4x4Luma_Intra_IAQMS[currMB->mb_iaqms_idx][qp_rem][0][0])<<(qp_per-6);
          img->cof[i][j1][0][0]=((M6[j]-M6[j1])*InvLevelScale4x4Luma_Intra_IAQMS[currMB->mb_iaqms_idx][qp_rem][0][0])<<(qp_per-6);
        }
      }
      else
      {
#endif
        if(qp_per<6)
        {
          img->cof[i][j][0][0] =((M6[j]+M6[j1])*InvLevelScale4x4Luma_Intra[qp_rem][0][0]+qp_const)>>(6-qp_per);
          img->cof[i][j1][0][0]=((M6[j]-M6[j1])*InvLevelScale4x4Luma_Intra[qp_rem][0][0]+qp_const)>>(6-qp_per);
        }
        else
        {
          img->cof[i][j][0][0] =((M6[j]+M6[j1])*InvLevelScale4x4Luma_Intra[qp_rem][0][0])<<(qp_per-6);
          img->cof[i][j1][0][0]=((M6[j]-M6[j1])*InvLevelScale4x4Luma_Intra[qp_rem][0][0])<<(qp_per-6);
        }
#ifdef ADAPTIVE_QUANTIZATION
      }
#endif
    }
  }
}


void itrans_sp(struct img_par *img,  //!< image parameters
               int ioff,             //!< index to 4x4 block
               int joff,             //!<
               int i0,               //!<
               int j0)               //!<
{
  int i,j,i1,j1;
  int m5[4];
  int m6[4];
  int predicted_block[BLOCK_SIZE][BLOCK_SIZE],ilev;
  
  int qp_per = (img->qp-MIN_QP)/6;
  int qp_rem = (img->qp-MIN_QP)%6;
  int q_bits    = Q_BITS+qp_per;

  int qp_per_sp = (img->qpsp-MIN_QP)/6;
  int qp_rem_sp = (img->qpsp-MIN_QP)%6;
  int q_bits_sp    = Q_BITS+qp_per_sp;
  int qp_const2=(1<<q_bits_sp)/2;  //sp_pred
  if (img->type == SI_SLICE) //ES modified
  {
    qp_per = (img->qpsp-MIN_QP)/6;
    qp_rem = (img->qpsp-MIN_QP)%6;
    q_bits = Q_BITS+qp_per;
  }

  for (j=0; j< BLOCK_SIZE; j++)
  for (i=0; i< BLOCK_SIZE; i++)
      predicted_block[i][j]=img->mpr[i+ioff][j+joff];
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

  //  Vertival transform

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

  for (j=0;j<BLOCK_SIZE;j++)
  for (i=0;i<BLOCK_SIZE;i++)
  {
    // recovering coefficient since they are already dequantized earlier
    img->cof[i0][j0][i][j]=(img->cof[i0][j0][i][j] >> qp_per) / dequant_coef[qp_rem][i][j]; 
    if(img->sp_switch || img->type==SI_SLICE)  //M.W. patched for SI
    {
      ilev=(abs(predicted_block[i][j]) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp; //ES added
      ilev= sign(ilev,predicted_block[i][j])+ img->cof[i0][j0][i][j];                           //ES added
      img->cof[i0][j0][i][j] = sign(abs(ilev) * dequant_coef[qp_rem_sp][i][j] << qp_per_sp ,ilev) ; //ES added 
    }                                                                                             //ES added
    else
    {                                                                                          //ES added
      ilev=((img->cof[i0][j0][i][j]*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6)+predicted_block[i][j] ;
      img->cof[i0][j0][i][j]=sign((abs(ilev) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp, ilev) * dequant_coef[qp_rem_sp][i][j] << qp_per_sp;
    }
  }
  // horizontal
  for (j=0;j<BLOCK_SIZE;j++)
  {
    for (i=0;i<BLOCK_SIZE;i++)
    {
      m5[i]=img->cof[i0][j0][i][j];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (i=0;i<2;i++)
    {
      i1=3-i;
      img->m7[i][j]=m6[i]+m6[i1];
      img->m7[i1][j]=m6[i]-m6[i1];
    }
  }
  // vertical
  for (i=0;i<BLOCK_SIZE;i++)
  {
    for (j=0;j<BLOCK_SIZE;j++)
      m5[j]=img->m7[i][j];

    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (j=0;j<2;j++)
    {
      j1=3-j;
      img->m7[i][j] =max(0,min(img->max_imgpel_value,(m6[j]+m6[j1]+DQ_ROUND)>>DQ_BITS));
      img->m7[i][j1]=max(0,min(img->max_imgpel_value,(m6[j]-m6[j1]+DQ_ROUND)>>DQ_BITS));
    }
  }
}

/*!
 ***********************************************************************
 * \brief
 *    The routine performs transform,quantization,inverse transform, adds the diff.
 *    to the prediction and writes the result to the decoded luma frame. Includes the
 *    RD constrained quantization also.
 *
 * \par Input:
 *    block_x,block_y: Block position inside a macro block (0,4,8,12).
 *
 * \par Output:
 *    nonzero: 0 if no levels are nonzero.  1 if there are nonzero levels. \n
 *    coeff_cost: Counter for nonzero coefficients, used to discard expencive levels.
 ************************************************************************
 */
void copyblock_sp(struct img_par *img,int block_x,int block_y)
{
  int sign(int a,int b);

  int i,j,i1,j1,m5[4],m6[4];

  int predicted_block[BLOCK_SIZE][BLOCK_SIZE];
  int qp_per = (img->qpsp-MIN_QP)/6;
  int qp_rem = (img->qpsp-MIN_QP)%6;
  int q_bits    = Q_BITS+qp_per;
  int qp_const2=(1<<q_bits)/2;  //sp_pred


  //  Horizontal transform
  for (j=0; j< BLOCK_SIZE; j++)
  for (i=0; i< BLOCK_SIZE; i++)
    predicted_block[i][j]=img->mpr[i+block_x][j+block_y];

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

  //  Vertival transform

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
  for (i=0; i < BLOCK_SIZE; i++)
    img->m7[i][j]=sign((abs(predicted_block[i][j])* quant_coef[qp_rem][i][j]+qp_const2)>> q_bits,predicted_block[i][j])*dequant_coef[qp_rem][i][j]<<qp_per;

  //     IDCT.
  //     horizontal

  for (j=0;j<BLOCK_SIZE;j++)
  {
    for (i=0;i<BLOCK_SIZE;i++)
    {
      m5[i]=img->m7[i][j];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (i=0;i<2;i++)
    {
      i1=3-i;
      img->m7[i][j]=m6[i]+m6[i1];
      img->m7[i1][j]=m6[i]-m6[i1];
    }
  }
  // vertical
  for (i=0;i<BLOCK_SIZE;i++)
  {
    for (j=0;j<BLOCK_SIZE;j++)
      m5[j]=img->m7[i][j];

    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (j=0;j<2;j++)
    {
      j1=3-j;
      img->m7[i][j] =max(0,min(img->max_imgpel_value,(m6[j]+m6[j1]+DQ_ROUND)>>DQ_BITS));
      img->m7[i][j1]=max(0,min(img->max_imgpel_value,(m6[j]-m6[j1]+DQ_ROUND)>>DQ_BITS));
    }
  }

  //  Decoded block moved to frame memory

  for (j=0; j < BLOCK_SIZE; j++)
    for (i=0; i < BLOCK_SIZE; i++)
      dec_picture->imgY[img->pix_y+block_y+j][img->pix_x+block_x+i]=img->m7[i][j];

}

void itrans_sp_chroma(struct img_par *img,int ll)
{
  int i,j,i1,j2,ilev,n2,n1,j1,mb_y;
  int m5[BLOCK_SIZE];
  int predicted_chroma_block[MB_BLOCK_SIZE/2][MB_BLOCK_SIZE/2],mp1[BLOCK_SIZE];
  int qp_per,qp_rem,q_bits;
  int qp_per_sp,qp_rem_sp,q_bits_sp,qp_const2;

  qp_per    = ((img->qp<0?img->qp:QP_SCALE_CR[img->qp])-MIN_QP)/6;
  qp_rem    = ((img->qp<0?img->qp:QP_SCALE_CR[img->qp])-MIN_QP)%6;
  q_bits    = Q_BITS+qp_per;

  qp_per_sp    = ((img->qpsp<0?img->qpsp:QP_SCALE_CR[img->qpsp])-MIN_QP)/6;
  qp_rem_sp    = ((img->qpsp<0?img->qpsp:QP_SCALE_CR[img->qpsp])-MIN_QP)%6;
  q_bits_sp    = Q_BITS+qp_per_sp;
  qp_const2=(1<<q_bits_sp)/2;  //sp_pred

  if (img->type == SI_SLICE)
  {
    qp_per    = ((img->qpsp < 0 ? img->qpsp : QP_SCALE_CR[img->qpsp]) - MIN_QP) / 6;
    qp_rem    = ((img->qpsp < 0 ? img->qpsp : QP_SCALE_CR[img->qpsp]) - MIN_QP) % 6;
    q_bits    = Q_BITS + qp_per;
  }

  for (j=0; j < MB_BLOCK_SIZE/2; j++)
  for (i=0; i < MB_BLOCK_SIZE/2; i++)
  {
    predicted_chroma_block[i][j]=img->mpr[i][j];
    img->mpr[i][j]=0;
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
  mp1[0]=(predicted_chroma_block[0][0]+predicted_chroma_block[4][0]+predicted_chroma_block[0][4]+predicted_chroma_block[4][4]);
  mp1[1]=(predicted_chroma_block[0][0]-predicted_chroma_block[4][0]+predicted_chroma_block[0][4]-predicted_chroma_block[4][4]);
  mp1[2]=(predicted_chroma_block[0][0]+predicted_chroma_block[4][0]-predicted_chroma_block[0][4]-predicted_chroma_block[4][4]);
  mp1[3]=(predicted_chroma_block[0][0]-predicted_chroma_block[4][0]-predicted_chroma_block[0][4]+predicted_chroma_block[4][4]);

  for (n1=0; n1 < 2; n1 ++)
  for (n2=0; n2 < 2; n2 ++)
  {
    if (img->sp_switch || img->type==SI_SLICE)  //M.W. patched for SI
    {
      //quantization fo predicted block
      ilev=(abs (mp1[n1+n2*2]) * quant_coef[qp_rem_sp][0][0] + 2 * qp_const2) >> (q_bits_sp + 1); 
      //addition   
      ilev=img->cof[n1+ll][4+n2][0][0]+sign(ilev,mp1[n1+n2*2]);                                   
      //dequantization
      mp1[n1+n2*2] =ilev*dequant_coef[qp_rem_sp][0][0]<<qp_per_sp;                                
    }   
    else
    {
      ilev=((img->cof[n1+ll][4+n2][0][0]*dequant_coef[qp_rem][0][0]*A[0][0]<< qp_per) >>5)+mp1[n1+n2*2] ;
      mp1[n1+n2*2]=sign((abs(ilev)* quant_coef[qp_rem_sp][0][0]+ 2 * qp_const2)>> (q_bits_sp+1),ilev)*dequant_coef[qp_rem_sp][0][0]<<qp_per_sp;
    }
  }


  for (n2=0; n2 < 2; n2 ++)
  for (n1=0; n1 < 2; n1 ++)
  for (i=0;i< BLOCK_SIZE; i++)
  for (j=0;j< BLOCK_SIZE; j++)
  {
  // recovering coefficient since they are already dequantized earlier
    img->cof[n1+ll][4+n2][i][j] = (img->cof[n1+ll][4+n2][i][j] >> qp_per) / dequant_coef[qp_rem][i][j];

    if (img->sp_switch || img->type==SI_SLICE)  //M.W. patched for SI
    {
      //quantization of the predicted block
      ilev =  (abs(predicted_chroma_block[n1*BLOCK_SIZE+i][n2*BLOCK_SIZE+j]) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp;
      //addition of the residual
      ilev = sign(ilev,predicted_chroma_block[n1*BLOCK_SIZE+i][n2*BLOCK_SIZE+j]) + img->cof[n1+ll][4+n2][i][j];
      // Inverse quantization 
      img->cof[n1+ll][4+n2][i][j] = ilev * dequant_coef[qp_rem_sp][i][j] << qp_per_sp  ;
    }
    else
    {
      //dequantization and addition of the predicted block
      ilev=((img->cof[n1+ll][4+n2][i][j]*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6)+predicted_chroma_block[n1*BLOCK_SIZE+i][n2*BLOCK_SIZE+j] ;
      //quantization and dequantization
      img->cof[n1+ll][4+n2][i][j] = sign((abs(ilev) * quant_coef[qp_rem_sp][i][j] + qp_const2)>> q_bits_sp,ilev)*dequant_coef[qp_rem_sp][i][j]<<qp_per_sp;
    }
  }
  img->cof[0+ll][4][0][0]=(mp1[0]+mp1[1]+mp1[2]+mp1[3])>>1;
  img->cof[1+ll][4][0][0]=(mp1[0]-mp1[1]+mp1[2]-mp1[3])>>1;
  img->cof[0+ll][5][0][0]=(mp1[0]+mp1[1]-mp1[2]-mp1[3])>>1;
  img->cof[1+ll][5][0][0]=(mp1[0]-mp1[1]-mp1[2]+mp1[3])>>1;
}

int sign(int a , int b)
{
  int x;

  x=abs(a);
  if (b>0)
    return(x);
  else return(-x);
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
*   Perform mode-dependent inverse 4x4 directional transform 
*
* \para itransklt_sep()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
void itransklt_sep( struct img_par *img, //!< image parameters
                    int ioff,            //!< index to 4x4 block
                    int joff,            //!<
                    int i0,              //!<
                    int j0,
                    int ipmode)
{
  int i, j, k; 
  int shift = 20, shiftby2 = 1<<(shift-1);
  int ilev[BLOCK_SIZE][BLOCK_SIZE], 
      temp2D[BLOCK_SIZE][BLOCK_SIZE], y2D [BLOCK_SIZE][BLOCK_SIZE];

  for(j = 0; j < BLOCK_SIZE; j++)
    for(i = 0; i < BLOCK_SIZE; i++)
    {
      ilev[j][i] = img->cof[i0][j0][i][j];
    }

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

  for(j = 0; j < BLOCK_SIZE; j++)
    for(i = 0; i < BLOCK_SIZE; i++)
    {
      int temp = max(0,min(img->max_imgpel_value,(y2D[j][i]+((long)img->mpr[i+ioff][j+joff] <<shift)+shiftby2)>>shift));
      img->m7[i][j] = temp;
    }

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
        temp2D[j][k] = KLTCol16[i][k][j];
    compute_inner_product16x16(temp2D, 1, KLTCol16inProd[i]);

    for(j = 0; j < 16; j++)
      for(k = 0; k < 16; k++)
        temp2D[j][k] = KLTRow16[i][k][j];
    compute_inner_product16x16(temp2D, 0, KLTRow16inProd[i]);
  }
}

/*! 
*************************************************************************************
* \brief
*   Perform fast mode-dependent inverse 16x16 directional transform 
*
* \para itransklt16x16_sep_fast()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
void itransklt16x16_sep_fast( struct img_par *img, //!< image parameters
                              int ipmode)
{
  int i, j, k; 
  int temp;
  int shift = 20, shiftby2 = 1<<(shift-1);
  int inprod[MB_BLOCK_SIZE], 
      temp2D[MB_BLOCK_SIZE][MB_BLOCK_SIZE], 
      y2D   [MB_BLOCK_SIZE][MB_BLOCK_SIZE];

  compute_inner_product16x16(img->cof16x16, 0, inprod);

  // inverse transform, column first
  for (i=0;i<(MB_BLOCK_SIZE);i++)
  {
    for (j=0;j<(MB_BLOCK_SIZE);j++)
    {
      temp2D[i][j] = 0;
      for(k = 0; k < MB_BLOCK_SIZE/2; k++)
        temp2D[i][j]+=
        (KLTCol16[ipmode][2*k][i]+img->cof16x16[2*k+1][j])*(KLTCol16[ipmode][2*k+1][i]+img->cof16x16[2*k][j]);
      temp2D[i][j] = temp2D[i][j] - KLTCol16inProd[ipmode][i] - inprod[j];
    }
  }

  compute_inner_product16x16(temp2D, 1, inprod);

  // inverse transform, row second 
  for (i=0;i<(MB_BLOCK_SIZE);i++)
  {
    for (j=0;j<(MB_BLOCK_SIZE);j++)
    {
      y2D[i][j] = 0;
      for(k = 0; k < MB_BLOCK_SIZE/2; k++)
        y2D[i][j]+=
        (temp2D[i][2*k+1]+KLTRow16[ipmode][j][2*k])*(temp2D[i][2*k]+KLTRow16[ipmode][j][2*k+1]);
      y2D[i][j] = y2D[i][j] - inprod[i] - KLTRow16inProd[ipmode][j];
      
      temp = max(0,min(img->max_imgpel_value,(y2D[i][j]+((long)img->mpr[j][i] <<shift)+shiftby2)>>shift));
      img->m7[j][i] = temp;
    }
  }
}

#endif

