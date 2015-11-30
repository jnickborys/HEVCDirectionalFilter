//
//  Decides which adaptive filters should be used
//
//
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "defines.h"
#include "adaptive_filter.h"
#include "header.h"
#include "global.h" 
#include "image.h"

#ifdef E_DAIF
#define int_pos             -1            // This should go where the defines for a_pos - o_pos are

#define STD_INTERPOLATION   0             // Interpolation selection in GetAIFPixel()
#define AIF_INTERPOLATION   1

#define IS_REF_VALID(F)     ((F) != -1)   // Checks reference

#ifndef EAIF
#define D_SIZE              6             // Size of the decision vector (one integer + 5 sub-pel position)
#endif

#ifdef EAIF
#define MAX_D_SIZE              12             // Size of the decision vector (one integer + 5 sub-pel position)
#endif

#define USE_POC_WEIGHTS_IN_RDDECISION

extern int FILTER_OFFSET_INT;
#ifndef EDAIF2
extern int FILTER_SIZE;
#endif
extern double FilterCoefInt[SQR_FILTER_INT];
extern double FilterCoef[15][SQR_FILTER];
extern double STANDARD_2D_FILTER_orig[15][SQR_FILTER];
extern double STANDARD_2D_FILTER[15][SQR_FILTER];

//
//  Maps an integer or a sub-pixel position to an index in the range [0,...,d_size-1] 
//
static int SpPos2Idx(int position, int d_size)
{
  switch(d_size)
  {
  case 6:
    {
      static int map[16] = {0, 1, 2, 1,         // I, a, b, a
        1, 3, 4, 3,         // a, e, f, e 
        2, 4, 5, 4,         // b, f, j, f
        1, 3, 4, 3};        // a, e, f, e
      return map[position + 1];
    }
  case 9:
    {
      static int map[16] = {0, 1, 2, 1,         // I, a, b, a
        3, 4, 5, 4,         // d, e, f, e 
        6, 7, 8, 7,         // h, i, j, i
        3, 4, 5, 4};        // d, e, f, e
      return map[position + 1];
    }
  case 16:
    {
      static int map[16] = {0 , 1,  2,  3,       // I, a, b, c
        4,  5,  6,  7,       // d, e, f, g 
        8 , 9, 10, 11,       // h, i, j, k
        12, 13, 14, 15};      // l, m, n, o
      return map[position + 1];
    }
  default:
    {
      assert(0);
      return -1;
    }
  }
}

//
//  Maps an index in [0,...,d_size-1] to the corresponding integer or sub-pixel position 
//
static int Idx2SpPos(int position, int d_size)
{
  switch(d_size)
  {
  case 6:
    {
      static int map[6] = {int_pos, a_pos, b_pos, 
        e_pos, f_pos, j_pos};
      return map[position];
    }
  case 9:
    {
      static int map[9] = {int_pos, a_pos, b_pos, 
        d_pos, e_pos, f_pos, 
        h_pos, i_pos, j_pos};
      return map[position];
    }
  case 16:
    {
      static int map[16] = {int_pos, a_pos, b_pos, c_pos,
        d_pos, e_pos, f_pos, g_pos,
        h_pos, i_pos, j_pos, k_pos,
        l_pos, m_pos, n_pos, o_pos};
      return map[position];
    }
  default:
    {
      assert(0);
      return -1;
    }
  }
}


//
//  Return pixel after interpolation with AIF or standard filter
//
static double GetAIFPixel(imgpel **refY, int i, int j, int sub_pos, int mvx, int mvy, int method)
{
  int ii, jj;
  double pix = 0.0;
  int pos_x, pos_y;
  int img_width  = dpb.fs_ref[0]->frame->size_x;
  int img_height = dpb.fs_ref[0]->frame->size_y;

  if(method == STD_INTERPOLATION)
  {
    if(sub_pos == int_pos)
    {
      // Integer position
      pos_y = FindPositionInt(img_height, i, FILTER_OFFSET_INT, mvy);
      pos_x = FindPositionInt(img_width, j, FILTER_OFFSET_INT, mvx);
      pix = refY[pos_y][pos_x];
    }
    else
    {
      // Sub-pel position
      for(ii = 0; ii < FILTER_SIZE; ++ii)
      {
        for(jj = 0; jj < FILTER_SIZE; ++jj)
        {
          pos_y = FindPosition(img_height, i + ii, 0, mvy);
          pos_x = FindPosition(img_width,  j + jj, 0, mvx);
          pix += refY[pos_y][pos_x] * (input->UseAdaptiveFilter == FILTER_TYPE_EDAIF ? 
            STANDARD_2D_FILTER_orig[sub_pos][FILTER_SIZE * ii + jj]:STANDARD_2D_FILTER[sub_pos][FILTER_SIZE * ii + jj]);
        }
      }
    }
  }
  else if(method == AIF_INTERPOLATION)
  {
    if(sub_pos == int_pos)
    {
      // Integer position
      for(ii = 0; ii < FILTER_SIZE_INT; ++ii)
      {
        for(jj = 0; jj < FILTER_SIZE_INT; ++jj)
        {
          pos_y = FindPositionInt(img_height, i + ii, 0, mvy);
          pos_x = FindPositionInt(img_width, j + jj, 0, mvx);
          pix += refY[pos_y][pos_x] * FilterCoefInt[FILTER_SIZE_INT * ii + jj];
        }
      }
      pix += FilterCoefInt[SQR(FILTER_SIZE_INT)];
    }
    else
    {
      // Sub-pel position
      for(ii = 0; ii < FILTER_SIZE; ++ii)
      {
        for(jj = 0; jj < FILTER_SIZE; ++jj)
        {
          pos_y = FindPosition(img_height, i + ii, 0, mvy);
          pos_x = FindPosition(img_width, j + jj, 0, mvx);
          pix += refY[pos_y][pos_x] * FilterCoef[sub_pos][FILTER_SIZE * ii + jj];
        }
      }
      pix += FilterCoef[sub_pos][SQR(FILTER_SIZE)];
    }
  }
  else
  {
    assert(0);
  }
  return pix;
}


static double GetAIFPixel2(imgpel **refY, int i, int j, int sub_pos, int mvx, int mvy, int method, int filterID)
{
  FilterEntity *p_Filter = &SetAIF[filterID];
  int ii, jj;
  double pix = 0.0;
  int pos_x, pos_y;
  int img_width  = dpb.fs_ref[0]->frame->size_x;
  int img_height = dpb.fs_ref[0]->frame->size_y;

  if(method == STD_INTERPOLATION)
  {
    if(sub_pos == int_pos)
    {
      // Integer position
      pos_y = FindPositionInt(img_height, i, FILTER_OFFSET_INT, mvy);
      pos_x = FindPositionInt(img_width, j, FILTER_OFFSET_INT, mvx);
      pix = refY[pos_y][pos_x];
    }
    else
    {
      // Sub-pel position
      for(ii = 0; ii < FILTER_SIZE; ++ii)
      {
        for(jj = 0; jj < FILTER_SIZE; ++jj)
        {
          pos_y = FindPosition(img_height, i + ii, 0, mvy);
          pos_x = FindPosition(img_width,  j + jj, 0, mvx);
          pix += refY[pos_y][pos_x] * (p_Filter->STANDARD_2D_FILTER[sub_pos][FILTER_SIZE * ii + jj]);
        }
      }
    }
  }
  else if(method == AIF_INTERPOLATION)
  {
    if(sub_pos == int_pos)
    {
      // Integer position
      for(ii = 0; ii < FILTER_SIZE_INT; ++ii)
      {
        for(jj = 0; jj < FILTER_SIZE_INT; ++jj)
        {
          pos_y = FindPositionInt(img_height, i + ii, 0, mvy);
          pos_x = FindPositionInt(img_width, j + jj, 0, mvx);
          pix += refY[pos_y][pos_x] * FilterCoefInt[FILTER_SIZE_INT * ii + jj];
        }
      }
      pix += FilterCoefInt[SQR(FILTER_SIZE_INT)];
    }
    else
    {
      // Sub-pel position
      for(ii = 0; ii < FILTER_SIZE; ++ii)
      {
        for(jj = 0; jj < FILTER_SIZE; ++jj)
        {
          pos_y = FindPosition(img_height, i + ii, 0, mvy);
          pos_x = FindPosition(img_width, j + jj, 0, mvx);
          pix += refY[pos_y][pos_x] * p_Filter->FilterCoef[sub_pos][FILTER_SIZE * ii + jj];
        }
      }
      pix += p_Filter->FilterCoef[sub_pos][SQR(FILTER_SIZE)];
    }
  }
  else
  {
    assert(0);
  }
  return pix;
}

//
//  Returns the sub-pel position or -1 if integer
//
static int GetSubPos(int mvx, int mvy)
{
  int mvx_sub = (mvx >= 0)? mvx % 4: (4 - abs(mvx) % 4) % 4;
  int mvy_sub = (mvy >= 0)? mvy % 4: (4 - abs(mvy) % 4) % 4;
  return(mvx_sub + 4 * mvy_sub - 1);          // pos 0..14 in a 4x4 block
}

//
//  Return the cost in bits of the filter corresponding to D[pos]
//  Returned value takes into account whether the filter has already been sent
//
#ifdef EAIF
static int BitCost(int pos, int nBits[MAX_D_SIZE], int D[MAX_D_SIZE])
#else
static int BitCost(int pos, int nBits[D_SIZE], int D[D_SIZE])
#endif
{
#ifdef EAIF
  int D_SIZE = input->UseAdaptiveFilter == FILTER_TYPE_EAIF ? 9:6;
#endif

  if((D_SIZE == 6) || (D_SIZE == 9))
  {
    return nBits[pos];
  }
  else if(D_SIZE == 16)
  {
    switch(pos)
    {
    case  0:    // int_pos
    case  1:    // a_pos
    case  2:    // b_pos
    case  4:    // d_pos
    case  5:    // e_pos
    case  6:    // f_pos
    case  8:    // h_pos
    case  9:    // i_pos
    case 10:    // j_pos
      return nBits[pos];              // Primary filter, send cost
    case  3:    // c_pos
      return ((D[1])? 0: nBits[1]);   // If filter has not been sent then return cost
    case  7:    // g_pos
      return ((D[5])? 0: nBits[5]);
    case 11:    // k_pos
      return ((D[9])? 0: nBits[9]);
    case 12:    // l_pos
      return ((D[4])? 0: nBits[4]);
    case 13:    // m_pos
      return ((D[5] || D[7])? 0: nBits[5]);
    case 14:    // n_pos
      return ((D[6])? 0: nBits[6]);
    case 15:    // o_pos
      return ((D[5] || D[7] || D[13])? 0: nBits[5]);
    default:
      {
        assert(0);
        return -1;
      }
    }
  }
  else
  {
    assert(0);
    return -1;
  }
}

//
//  Decide filter usage by inspecting the diagonal of the cost matrix
//
#ifdef EAIF
static void Diagonal(double C[2][2][MAX_D_SIZE][MAX_D_SIZE], int AvailFilter[MAX_D_SIZE], int nBits[MAX_D_SIZE], int D[MAX_D_SIZE])
#else
static void Diagonal(double C[2][2][D_SIZE][D_SIZE], int AvailFilter[D_SIZE], int nBits[D_SIZE], int D[D_SIZE])
#endif
{
  int i;
  double adaptiveCost;
  double standardCost;
  double lambda; 
  int qp = input->UseRDO_Q ? img->masterQP:img->qp;
#ifdef EAIF
  int D_SIZE = input->UseAdaptiveFilter == FILTER_TYPE_EAIF ? 9:6;
#endif

  if ((img->type==B_SLICE) && img->nal_reference_idc)
  {
    lambda = img->lambda_me[5][qp];  
  }
  else
  {
    lambda = img->lambda_me[img->type][qp];
  }

  for(i = 0; i < D_SIZE; ++i)
  {
    if(AvailFilter[i])
    {
      standardCost = C[0][0][i][i];
      adaptiveCost = C[1][1][i][i] + BitCost(i, nBits, D) * lambda;
      D[i] = (standardCost <= adaptiveCost)? 0: 1;
    }
    else
    {
      D[i] = 0;
    }
  }
}

//
//  Check whether a decision is incompatible with the filters available
//  If it is, converts the decision to a vector D[]
//
#ifdef EAIF
static int is_incompatible(int decision, int AvailFilter[MAX_D_SIZE], int D[MAX_D_SIZE])
#else
static int is_incompatible(int decision, int AvailFilter[D_SIZE], int D[D_SIZE])
#endif
{
  int i;
#ifdef EAIF
  int D_SIZE = input->UseAdaptiveFilter == FILTER_TYPE_EAIF ? 9:6;
#endif

  for(i = 0; i < D_SIZE; ++i)
  {
    D[i] = decision & 0x1;
    if(D[i] && !AvailFilter[i])     // No corresponding filter
      return 1;
    else
      decision >>= 1;
  }
  return 0;
}

//
//  Decide filter usage by exhaustive search
//
#ifdef EAIF
static void Exhaustive(double C[2][2][MAX_D_SIZE][MAX_D_SIZE], int AvailFilter[MAX_D_SIZE], int nBits[MAX_D_SIZE], int D[MAX_D_SIZE])
#else
static void Exhaustive(double C[2][2][D_SIZE][D_SIZE], int AvailFilter[D_SIZE], int nBits[D_SIZE], int D[D_SIZE])
#endif
{
  int i, j;
  int decision;
  double BestCost = DBL_MAX;
  double TmpCost; 
#ifdef EAIF
  int TmpD[MAX_D_SIZE] = {0};
  int D_SIZE = input->UseAdaptiveFilter == FILTER_TYPE_EAIF ? 9:6;
#else
  int TmpD[D_SIZE] = {0};
#endif
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

  for(decision = 0; decision < (1 << D_SIZE); ++decision)
  {
    if(is_incompatible(decision, AvailFilter, TmpD))
      continue;

    TmpCost = 0.0;

    for(i = 0; i < D_SIZE; ++i)
      for(j = 0; j < D_SIZE; ++j)
        TmpCost += C[TmpD[i]][TmpD[j]][i][j];

    for(i = 0; i < D_SIZE; ++i)
      if(TmpD[i])
        TmpCost += BitCost(i, nBits, TmpD) * lambda;

    if(TmpCost < BestCost)
    {
      BestCost = TmpCost;
      for(i = 0; i < D_SIZE; ++i)
        D[i] = TmpD[i];
    }
  }
}

//
//  SAD-based optimal decision on the integer and sub-pel positions
//
void AIFDecision(int *FilterFlagInt, int *FilterFlag, int *SymmetryPosition)
{
  int img_width, img_height;
  int mbx, mby; 
  int x_orig, y_orig;
  int curr_mb_nr;
  int subblock;
  int max_mb_nr = (img->width * img->height) / 256;
#ifdef EAIF
  double C[2][2][MAX_D_SIZE][MAX_D_SIZE] = {{{{0.0}}}};
  int nBits[MAX_D_SIZE] = {0};                // Bits required to encode offset and filter coefficients
  int AvailFilter[MAX_D_SIZE] = {0};          // 1 if the corresponding filter is available
  int D[MAX_D_SIZE] = {0};                    // Decision vector (0/1 don't use/use corresponding filter)
#else
  double C[2][2][D_SIZE][D_SIZE] = {{{{0.0}}}};
  int nBits[D_SIZE] = {0};                // Bits required to encode offset and filter coefficients
  int AvailFilter[D_SIZE] = {0};          // 1 if the corresponding filter is available
  int D[D_SIZE] = {0};                    // Decision vector (0/1 don't use/use corresponding filter)
#endif
  int x, y;
  int ref_frameF, ref_frameB;
  int F_nrefs = 0;                        // Counts number of F and B references
  int B_nrefs = 0;
  int mvxF = 0;
  int mvyF = 0;
  int mvxB = 0;
  int mvyB = 0;
  int sub_posF = 0;
  int sub_posB = 0;
  imgpel **refYF = NULL;
  imgpel **refYB = NULL;
  int sub_pos;
  float b_scaling = 1.0;
  double pix_00, pix_01, pix_10, pix_11;
  imgpel rec_00, rec_01, rec_10, rec_11;
  imgpel orig;
  int filterF, filterB;
  int i, j;
  double pixAifF, pixAifB;
  double pixStdF, pixStdB;

#ifdef  INTERNAL_BIT_DEPTH_INCREASE
  imgpel max_pixel_value = (1 << img->bitdepth_luma) - 1;
#else
  imgpel max_pixel_value = 255;
#endif

#ifdef EAIF
  int D_SIZE = input->UseAdaptiveFilter == FILTER_TYPE_EAIF ? 9:6;
#endif

  img_width  = dpb.fs_ref[0]->frame->size_x;
  img_height = dpb.fs_ref[0]->frame->size_y;

  // Cost of encoding offset and coefficients
#ifdef EAIF
  if(*FilterFlagInt && (input->UseAdaptiveFilter == FILTER_TYPE_EDAIF || input->UseAdaptiveFilter == FILTER_TYPE_EAIF))
#else
  if(*FilterFlagInt && input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)
#endif
  {
    nBits[SpPos2Idx(int_pos, D_SIZE)] = sendAIFInteger(NULL) + sendAIFOffset(-1, NULL);
    AvailFilter[SpPos2Idx(int_pos, D_SIZE)] = 1;
  }
  for(sub_pos = a_pos; sub_pos <= o_pos; ++sub_pos)
  {
    if(SymmetryPosition[sub_pos] && FilterFlag[sub_pos])
    {
#ifdef EAIF
      nBits[SpPos2Idx(sub_pos, D_SIZE)] = 
        input->UseAdaptiveFilter != FILTER_TYPE_EAIF ? estimateCostOfCoefs_EDAIF(sub_pos):estimateCostOfCoefs_EAIF(sub_pos);
#else
      nBits[SpPos2Idx(sub_pos, D_SIZE)] = estimateCostOfCoefs_EDAIF(sub_pos);
#endif
#ifdef EAIF
      if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF || input->UseAdaptiveFilter == FILTER_TYPE_EAIF)
#else
      if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)
#endif
        nBits[SpPos2Idx(sub_pos, D_SIZE)] += sendAIFOffset(sub_pos, NULL);
      AvailFilter[SpPos2Idx(sub_pos, D_SIZE)] |= FilterFlag[sub_pos];
    }
  }

  // SAD computation
  for(curr_mb_nr = 0; curr_mb_nr < max_mb_nr; ++curr_mb_nr)
  {
    mbx = curr_mb_nr % (img_width / MB_BLOCK_SIZE);
    mby = curr_mb_nr / (img_width / MB_BLOCK_SIZE); 
    x_orig = mbx * MB_BLOCK_SIZE;
    y_orig = mby * MB_BLOCK_SIZE;

    if(IS_INTRA(&img->mb_data[curr_mb_nr]))
      continue; 

    for(subblock = 0; subblock < 16; ++subblock)
    {
      x = x_orig + 4 * (subblock % 4);
      y = y_orig + 4 * (subblock / 4);

      ref_frameF = enc_picture->ref_idx[LIST_0][y / 4][x / 4];
      ref_frameB = enc_picture->ref_idx[LIST_1][y / 4][x / 4];

      if(IS_REF_VALID(ref_frameF))
      {
        mvxF = enc_picture->mv[LIST_0][y / 4][x / 4][0];
        mvyF = enc_picture->mv[LIST_0][y / 4][x / 4][1];
        sub_posF = GetSubPos(mvxF, mvyF);
        refYF = listX[LIST_0][ref_frameF]->imgY;
        ++F_nrefs;
      }

      if(IS_REF_VALID(ref_frameB))
      {
        mvxB = enc_picture->mv[LIST_1][y / 4][x / 4][0];
        mvyB = enc_picture->mv[LIST_1][y / 4][x / 4][1];
        sub_posB = GetSubPos(mvxB, mvyB);
        refYB = listX[LIST_1][ref_frameB]->imgY; 
        ++B_nrefs;
      }

      // b_scaling computation
      if(IS_REF_VALID(ref_frameF) && IS_REF_VALID(ref_frameB))
      {
#ifdef USE_POC_WEIGHTS_IN_RDDECISION
        int td, tb;
        td = Clip3(-128, 127, (listX[LIST_1][ref_frameB]->poc - listX[LIST_0][ref_frameF]->poc));
        tb = Clip3(-128, 127, (enc_picture->poc - listX[LIST_0][ref_frameF]->poc));
        if (td != 0)
          b_scaling = (float) tb / (float) td;      // POC-based
        else
#endif
          b_scaling = 0.5;                          // 1/2
      }
      else if(IS_REF_VALID(ref_frameF) || IS_REF_VALID(ref_frameB))
      {
        b_scaling = 1.0; 
      }

      // Check if the filters exist and is primary
      filterF = IS_REF_VALID(ref_frameF)? ((sub_posF == int_pos)? *FilterFlagInt: FilterFlag[sub_posF]): 0;
      filterB = IS_REF_VALID(ref_frameB)? ((sub_posB == int_pos)? *FilterFlagInt: FilterFlag[sub_posB]): 0;
      if(!filterF && !filterB)
        continue;


      // Bidirectional
      if(IS_REF_VALID(ref_frameF) && IS_REF_VALID(ref_frameB))
      {
        for(i = y; i < y + 4; ++i)
        {
          for(j = x; j < x + 4; ++j)
          {
            if(filterF)
              pixAifF = GetAIFPixel(refYF, i, j, sub_posF, mvxF, mvyF, AIF_INTERPOLATION);
            else
              pixAifF = GetAIFPixel(refYF, i, j, sub_posF, mvxF, mvyF, STD_INTERPOLATION);

            if(filterB)
              pixAifB = GetAIFPixel(refYB, i, j, sub_posB, mvxB, mvyB, AIF_INTERPOLATION);
            else
              pixAifB = GetAIFPixel(refYB, i, j, sub_posB, mvxB, mvyB, STD_INTERPOLATION);

            pixStdF = GetAIFPixel(refYF, i, j, sub_posF, mvxF, mvyF, STD_INTERPOLATION);
            pixStdB = GetAIFPixel(refYB, i, j, sub_posB, mvxB, mvyB, STD_INTERPOLATION);

            pix_00 = (1 - b_scaling) * pixStdF + b_scaling * pixStdB;
            pix_01 = (1 - b_scaling) * pixStdF + b_scaling * pixAifB;
            pix_10 = (1 - b_scaling) * pixAifF + b_scaling * pixStdB;
            pix_11 = (1 - b_scaling) * pixAifF + b_scaling * pixAifB;

            rec_00 = max(0, min(max_pixel_value, (int)(pix_00 + 0.5)));
            rec_01 = max(0, min(max_pixel_value, (int)(pix_01 + 0.5)));
            rec_10 = max(0, min(max_pixel_value, (int)(pix_10 + 0.5)));
            rec_11 = max(0, min(max_pixel_value, (int)(pix_11 + 0.5)));

            orig = (int)(imgY_org[i][j]);

            if(filterF && filterB)
            {
              // Both filters are valid
              C[0][0][SpPos2Idx(sub_posF, D_SIZE)][SpPos2Idx(sub_posB, D_SIZE)] += (double)abs(orig - rec_00);
              C[1][1][SpPos2Idx(sub_posF, D_SIZE)][SpPos2Idx(sub_posB, D_SIZE)] += (double)abs(orig - rec_11);

              if(sub_posF != sub_posB)
              {
                // Sub-pel positions are not using the same filter
                C[0][1][SpPos2Idx(sub_posF, D_SIZE)][SpPos2Idx(sub_posB, D_SIZE)] += (double)abs(orig - rec_01);
                C[1][0][SpPos2Idx(sub_posF, D_SIZE)][SpPos2Idx(sub_posB, D_SIZE)] += (double)abs(orig - rec_10);
              }
            } 
            else if(filterF)
            {
              // Filter F is valid
              C[0][0][SpPos2Idx(sub_posF, D_SIZE)][SpPos2Idx(sub_posB, D_SIZE)] += (double)abs(orig - rec_00);
              C[1][0][SpPos2Idx(sub_posF, D_SIZE)][SpPos2Idx(sub_posB, D_SIZE)] += (double)abs(orig - rec_10);
            }
            else
            {
              // Filter B is valid
              C[0][0][SpPos2Idx(sub_posF, D_SIZE)][SpPos2Idx(sub_posB, D_SIZE)] += (double)abs(orig - rec_00);
              C[0][1][SpPos2Idx(sub_posF, D_SIZE)][SpPos2Idx(sub_posB, D_SIZE)] += (double)abs(orig - rec_01);
            }
          }
        }
      }
      else
      {
        if(IS_REF_VALID(ref_frameB))
        {
          if(!filterB)
          {
            continue;
          }
          else
          {
            sub_posF = sub_posB;
            mvxF = mvxB; 
            mvyF = mvyB; 
            refYF = refYB;
          }
        }
        else if(!filterF)
        {
          continue;
        }

        for(i = y; i < y + 4; ++i)
        {
          for(j = x; j < x + 4; ++j)
          {            
            pix_00 = GetAIFPixel(refYF, i, j, sub_posF, mvxF, mvyF, STD_INTERPOLATION);
            pix_11 = GetAIFPixel(refYF, i, j, sub_posF, mvxF, mvyF, AIF_INTERPOLATION);

            rec_00 = max(0, min(max_pixel_value, (int)(pix_00 + 0.5)));
            rec_11 = max(0, min(max_pixel_value, (int)(pix_11 + 0.5)));

            orig   = (int)(imgY_org[i][j]);

            C[0][0][SpPos2Idx(sub_posF, D_SIZE)][SpPos2Idx(sub_posF, D_SIZE)] += (double)abs(orig - rec_00);
            C[1][1][SpPos2Idx(sub_posF, D_SIZE)][SpPos2Idx(sub_posF, D_SIZE)] += (double)abs(orig - rec_11);
          }
        }
      }     // Bidirectional if
    }     // Subblock loop
  }     // Macroblock loop

  // Decision
  if((F_nrefs == 0) || (B_nrefs == 0))
    Diagonal(C, AvailFilter, nBits, D);            // All references from the same direction
  else
    Exhaustive(C, AvailFilter, nBits, D);

  // Output
  *FilterFlagInt = D[0];
  for(i = 1; i < D_SIZE; ++i)
    FilterFlag[Idx2SpPos(i, D_SIZE)] = D[i];
}


#ifdef EDAIF2
void AIFDecision2(int *FilterFlagInt, int *FilterFlag, int *SymmetryPosition, int *FILTER_NEW_SUB_POS)
{
  //extern int FilterFlagInt; 
  int img_width, img_height;
  int mbx, mby; 
  int x_orig, y_orig;
  int curr_mb_nr;
  int subblock;
  int max_mb_nr = (img->width * img->height) / 256;
  double C[2][2][MAX_NUM_AIF+1][MAX_NUM_AIF+1] = {{{{0.0}}}};
  int nBits[MAX_NUM_AIF+1] = {0};                // Bits required to encode offset and filter coefficients
  int AvailFilter[MAX_NUM_AIF+1] = {0};          // 1 if the corresponding filter is available

  int x, y;
  int ref_frameF, ref_frameB;
  int F_nrefs = 0;                        // Counts number of F and B references
  int B_nrefs = 0;
  int mvxF = 0;
  int mvyF = 0;
  int mvxB = 0;
  int mvyB = 0;
  int sub_posF = 0;
  int sub_posB = 0;
  imgpel **refYF = NULL;
  imgpel **refYB = NULL;
  int sub_pos;
  float b_scaling = 1.0;
  double pix_00, pix_01, pix_10, pix_11;
  imgpel rec_00, rec_01, rec_10, rec_11;
  imgpel orig;
  int filterF, filterB;
  int i, j;
  double pixAifF, pixAifB;
  double pixStdF, pixStdB;

  double costMCP[2][MAX_NUM_AIF+1];

#ifdef  INTERNAL_BIT_DEPTH_INCREASE
  imgpel max_pixel_value = (1 << img->bitdepth_luma) - 1;
#else
  imgpel max_pixel_value = 255;
#endif

  img_width  = dpb.fs_ref[0]->frame->size_x;
  img_height = dpb.fs_ref[0]->frame->size_y;

  // Cost of encoding offset and coefficients

  nBits[MAX_NUM_AIF] = sendAIFInteger(NULL) + sendAIFOffset(-1, NULL);
  AvailFilter[MAX_NUM_AIF] = 0;
  for(sub_pos = a_pos; sub_pos < MAX_NUM_AIF; sub_pos++)
  {
    nBits[sub_pos] = 0;
    costMCP[0][sub_pos] = 4161601.0;; // initialize
    costMCP[1][sub_pos] = 4161601.0;;// initialize

    AvailFilter[sub_pos] = 0;
    if(SymmetryPosition[sub_pos] && (FilterFlag[sub_pos]%2!=0))// if filter is adaptive and is primary estimate bit-overhead
    {
      nBits[sub_pos] = estimateCostOfCoefs_EDAIF(sub_pos);
      nBits[sub_pos] += sendAIFOffset(sub_pos, NULL);
    }
  }


  // SAD computation
  for(curr_mb_nr = 0; curr_mb_nr < max_mb_nr; ++curr_mb_nr)
  {
    mbx = curr_mb_nr % (img_width / MB_BLOCK_SIZE);
    mby = curr_mb_nr / (img_width / MB_BLOCK_SIZE); 
    x_orig = mbx * MB_BLOCK_SIZE;
    y_orig = mby * MB_BLOCK_SIZE;

    if(IS_INTRA(&img->mb_data[curr_mb_nr]))
      continue; 

    for(subblock = 0; subblock < 16; ++subblock)
    {
      x = x_orig + 4 * (subblock % 4);
      y = y_orig + 4 * (subblock / 4);

      ref_frameF = enc_picture->ref_idx[LIST_0][y / 4][x / 4];
      ref_frameB = enc_picture->ref_idx[LIST_1][y / 4][x / 4];

      if(IS_REF_VALID(ref_frameF))
      {
        mvxF = enc_picture->mv[LIST_0][y / 4][x / 4][0];
        mvyF = enc_picture->mv[LIST_0][y / 4][x / 4][1];
        sub_posF = GetSubPos(mvxF, mvyF);
        refYF = listX[LIST_0][ref_frameF]->imgY;
        ++F_nrefs;
      }

      if(IS_REF_VALID(ref_frameB))
      {
        mvxB = enc_picture->mv[LIST_1][y / 4][x / 4][0];
        mvyB = enc_picture->mv[LIST_1][y / 4][x / 4][1];
        sub_posB = GetSubPos(mvxB, mvyB);
        refYB = listX[LIST_1][ref_frameB]->imgY; 
        ++B_nrefs;
      }

      // b_scaling computation
      if(IS_REF_VALID(ref_frameF) && IS_REF_VALID(ref_frameB))
      {
#ifdef USE_POC_WEIGHTS_IN_RDDECISION
        int td, tb;
        td = Clip3(-128, 127, (listX[LIST_1][ref_frameB]->poc - listX[LIST_0][ref_frameF]->poc));
        tb = Clip3(-128, 127, (enc_picture->poc - listX[LIST_0][ref_frameF]->poc));
        if (td != 0)
          b_scaling = (float) tb / (float) td;      // POC-based
        else
#endif
          b_scaling = 0.5;                          // 1/2
      }
      else if(IS_REF_VALID(ref_frameF) || IS_REF_VALID(ref_frameB))
      {
        b_scaling = 1.0; 
      }

      // Check if the filters exist and is primary
      filterF = IS_REF_VALID(ref_frameF)? ((sub_posF == int_pos)? *FilterFlagInt: (FilterFlag[sub_posF]%2!=0)): 0;
      filterB = IS_REF_VALID(ref_frameB)? ((sub_posB == int_pos)? *FilterFlagInt: (FilterFlag[sub_posB]%2!=0)): 0;
      if(!filterF && !filterB)
        continue;


      // Bidirectional
      if(IS_REF_VALID(ref_frameF) && IS_REF_VALID(ref_frameB))
      {
        int mother_pelF = (sub_posF!=-1) ? FILTER_NEW_SUB_POS[sub_posF]:MAX_NUM_AIF;
        int mother_pelB = (sub_posB!=-1) ? FILTER_NEW_SUB_POS[sub_posB]:MAX_NUM_AIF;
        AvailFilter[mother_pelF] = 1;// it filter was used
        AvailFilter[mother_pelB] = 1;// it filter was used
        for(i = y; i < y + 4; ++i)
        {
          for(j = x; j < x + 4; ++j)
          {
            if(filterF)
              pixAifF = GetAIFPixel(refYF, i, j, sub_posF, mvxF, mvyF, AIF_INTERPOLATION);
            else
              pixAifF = GetAIFPixel(refYF, i, j, sub_posF, mvxF, mvyF, STD_INTERPOLATION);

            if(filterB)
              pixAifB = GetAIFPixel(refYB, i, j, sub_posB, mvxB, mvyB, AIF_INTERPOLATION);
            else
              pixAifB = GetAIFPixel(refYB, i, j, sub_posB, mvxB, mvyB, STD_INTERPOLATION);

            pixStdF = GetAIFPixel(refYF, i, j, sub_posF, mvxF, mvyF, STD_INTERPOLATION);
            pixStdB = GetAIFPixel(refYB, i, j, sub_posB, mvxB, mvyB, STD_INTERPOLATION);

            pix_00 = (1 - b_scaling) * pixStdF + b_scaling * pixStdB;
            pix_01 = (1 - b_scaling) * pixStdF + b_scaling * pixAifB;
            pix_10 = (1 - b_scaling) * pixAifF + b_scaling * pixStdB;
            pix_11 = (1 - b_scaling) * pixAifF + b_scaling * pixAifB;

            rec_00 = max(0, min(max_pixel_value, (int)(pix_00 + 0.5)));
            rec_01 = max(0, min(max_pixel_value, (int)(pix_01 + 0.5)));
            rec_10 = max(0, min(max_pixel_value, (int)(pix_10 + 0.5)));
            rec_11 = max(0, min(max_pixel_value, (int)(pix_11 + 0.5)));

            orig = (imgpel)(imgY_org[i][j]);

            if(filterF && filterB)
            {
              // Both filters are valid
              C[0][0][mother_pelF][mother_pelB] += (double)abs(orig - rec_00);
              C[1][1][mother_pelF][mother_pelB] += (double)abs(orig - rec_11);

              if(sub_posF != sub_posB)
              {
                // Sub-pel positions are not using the same filter
                C[0][1][mother_pelF][mother_pelB] += (double)abs(orig - rec_01);
                C[1][0][mother_pelF][mother_pelB] += (double)abs(orig - rec_10);
              }
            } 
            else if(filterF)
            {
              // Filter F is valid
              C[0][0][mother_pelF][mother_pelB] += (double)abs(orig - rec_00);
              C[1][0][mother_pelF][mother_pelB] += (double)abs(orig - rec_10);
            }
            else
            {
              // Filter B is valid
              C[0][0][mother_pelF][mother_pelB] += (double)abs(orig - rec_00);
              C[0][1][mother_pelF][mother_pelB] += (double)abs(orig - rec_01);
            }
          }
        }
      }
      else
      {
        int mother_pel = (sub_posF!=-1) ? FILTER_NEW_SUB_POS[sub_posF]:MAX_NUM_AIF;
        AvailFilter[mother_pel] = 1;// it filter was used
        if(IS_REF_VALID(ref_frameB))
        {
          if(!filterB)
          {
            continue;
          }
          else
          {
            sub_posF = sub_posB;
            mvxF = mvxB; 
            mvyF = mvyB; 
            refYF = refYB;
          }
        }
        else if(!filterF)
        {
          continue;
        }

        for(i = y; i < y + 4; ++i)
        {
          for(j = x; j < x + 4; ++j)
          {            
            pix_00 = GetAIFPixel(refYF, i, j, sub_posF, mvxF, mvyF, STD_INTERPOLATION);
            pix_11 = GetAIFPixel(refYF, i, j, sub_posF, mvxF, mvyF, AIF_INTERPOLATION);

            rec_00 = max(0, min(max_pixel_value, (int)(pix_00 + 0.5)));
            rec_11 = max(0, min(max_pixel_value, (int)(pix_11 + 0.5)));

            orig   = (int)(imgY_org[i][j]);

            C[0][0][0][mother_pel] += (double)abs(orig - rec_00);
            C[1][1][0][mother_pel] += (double)abs(orig - rec_11);
          }
        }
      }     // Bidirectional if
    }     // Subblock loop
  }     // Macroblock loop

  {// decission on filters
    double lambda; 
    int qp = input->UseRDO_Q ? img->masterQP:img->qp;
    if ((img->type==B_SLICE) && img->nal_reference_idc)
    {
      lambda = img->lambda_me[5][qp];  
    }
    else
    {
      lambda = img->lambda_me[img->type][qp];
    }
    {
      for(sub_pos = a_pos; sub_pos < MAX_NUM_AIF; sub_pos++)
      {
        costMCP[0][sub_pos] = 0.0;
        costMCP[1][sub_pos] = lambda*nBits[sub_pos];  // filters budget
      }
      //if(input->checkAIF_MCP)
      {
        for(sub_pos = a_pos; sub_pos < MAX_NUM_AIF; sub_pos++)
        {
          int mother_pel = FILTER_NEW_SUB_POS[FILTER_NEW_SUB_POS[sub_pos]];
          if((AvailFilter[sub_pos]))
          {
            costMCP[0][mother_pel] += C[0][0][0][sub_pos];
            costMCP[1][mother_pel] += C[1][1][0][sub_pos];  // filters budget
          }
        }
        for(sub_pos = a_pos; sub_pos < MAX_NUM_AIF; sub_pos++)
        {
          if((SymmetryPosition[sub_pos])&&(AvailFilter[sub_pos]))
          {
            if(costMCP[0][sub_pos]<costMCP[1][sub_pos])
            {// Assign the static filter of the same fitler structure (DAIF or RAIF)
              FilterFlag[sub_pos] = (FilterFlag[sub_pos]%2!=0) ? FilterFlag[sub_pos]-1: FilterFlag[sub_pos];
            }
          }
        }
      }
      sub_pos = MAX_NUM_AIF;
      if(AvailFilter[sub_pos])
      {
        costMCP[0][sub_pos] = C[0][0][0][sub_pos];
        costMCP[1][sub_pos] = C[1][1][0][sub_pos] + lambda*nBits[sub_pos];  // filters budget
      }
      *FilterFlagInt = (costMCP[0][sub_pos]<costMCP[1][sub_pos])? 0: 1;
    }
  }
}

void GetErrorMCP(int *FilterFlagInt,int filterID)
{
  int img_width, img_height;
  int mbx, mby; 
  int x_orig, y_orig;
  int curr_mb_nr;
  int subblock;
  int max_mb_nr = (img->width * img->height) / 256;
  double C[2][2][MAX_NUM_AIF+1][MAX_NUM_AIF+1] = {{{{0.0}}}};
  int nBits[MAX_NUM_AIF+1] = {0};                // Bits required to encode offset and filter coefficients
  int AvailFilter[MAX_NUM_AIF+1] = {0};          // 1 if the corresponding filter is available

  int x, y;
  int ref_frameF, ref_frameB;
  int F_nrefs = 0;                        // Counts number of F and B references
  int B_nrefs = 0;
  int mvxF = 0;
  int mvyF = 0;
  int mvxB = 0;
  int mvyB = 0;
  int sub_posF = 0;
  int sub_posB = 0;
  imgpel **refYF = NULL;
  imgpel **refYB = NULL;
  int sub_pos;
  float b_scaling = 1.0;
  double pix_00, pix_01, pix_10, pix_11;
  imgpel rec_00, rec_01, rec_10, rec_11;
  imgpel orig;
  int filterF, filterB;
  int i, j;
  double pixAifF, pixAifB;
  double pixStdF, pixStdB;

   FilterEntity *p_Filter = &SetAIF[filterID];

#ifdef  INTERNAL_BIT_DEPTH_INCREASE
  imgpel max_pixel_value = (1 << img->bitdepth_luma) - 1;
#else
  imgpel max_pixel_value = 255;
#endif

  img_width  = dpb.fs_ref[0]->frame->size_x;
  img_height = dpb.fs_ref[0]->frame->size_y;

  // Cost of encoding offset and coefficients

  nBits[MAX_NUM_AIF] = 0;
  AvailFilter[MAX_NUM_AIF] = 0;
  for(sub_pos = a_pos; sub_pos < MAX_NUM_AIF; sub_pos++)
  {
    nBits[sub_pos] = 0;
    costMCP[2*(filterID)][sub_pos] = 4161601.0;; // initialize
    costMCP[2*(filterID)+1][sub_pos] = 4161601.0;;// initialize
    AvailFilter[sub_pos] = 0;
    if(p_Filter->SymmetryPosition[sub_pos] && (p_Filter->FilterFlag[sub_pos]%2!=0))
    {// if filter is adaptive and is primary estimate bit-overhead
      nBits[sub_pos] = 0;
      nBits[sub_pos] += 0;
    }
  }

  // SAD computation
  for(curr_mb_nr = 0; curr_mb_nr < max_mb_nr; ++curr_mb_nr)
  {
    mbx = curr_mb_nr % (img_width / MB_BLOCK_SIZE);
    mby = curr_mb_nr / (img_width / MB_BLOCK_SIZE); 
    x_orig = mbx * MB_BLOCK_SIZE;
    y_orig = mby * MB_BLOCK_SIZE;

    if(IS_INTRA(&img->mb_data[curr_mb_nr]))
      continue; 

    for(subblock = 0; subblock < 16; ++subblock)
    {
      x = x_orig + 4 * (subblock % 4);
      y = y_orig + 4 * (subblock / 4);

      ref_frameF = enc_picture->ref_idx[LIST_0][y / 4][x / 4];
      ref_frameB = enc_picture->ref_idx[LIST_1][y / 4][x / 4];

      if(IS_REF_VALID(ref_frameF))
      {
        mvxF = enc_picture->mv[LIST_0][y / 4][x / 4][0];
        mvyF = enc_picture->mv[LIST_0][y / 4][x / 4][1];
        sub_posF = GetSubPos(mvxF, mvyF);
        refYF = listX[LIST_0][ref_frameF]->imgY;
        ++F_nrefs;
      }

      if(IS_REF_VALID(ref_frameB))
      {
        mvxB = enc_picture->mv[LIST_1][y / 4][x / 4][0];
        mvyB = enc_picture->mv[LIST_1][y / 4][x / 4][1];
        sub_posB = GetSubPos(mvxB, mvyB);
        refYB = listX[LIST_1][ref_frameB]->imgY; 
        ++B_nrefs;
      }

      // b_scaling computation
      if(IS_REF_VALID(ref_frameF) && IS_REF_VALID(ref_frameB))
      {
#ifdef USE_POC_WEIGHTS_IN_RDDECISION
        int td, tb;
        td = Clip3(-128, 127, (listX[LIST_1][ref_frameB]->poc - listX[LIST_0][ref_frameF]->poc));
        tb = Clip3(-128, 127, (enc_picture->poc - listX[LIST_0][ref_frameF]->poc));
        if (td != 0)
          b_scaling = (float) tb / (float) td;      // POC-based
        else
#endif
          b_scaling = 0.5;                          // 1/2
      }
      else if(IS_REF_VALID(ref_frameF) || IS_REF_VALID(ref_frameB))
      {
        b_scaling = 1.0; 
      }

      // Check if the filters exist and is primary
      filterF = IS_REF_VALID(ref_frameF)? ((sub_posF == int_pos)? *FilterFlagInt: (p_Filter->FilterFlag[sub_posF]%2!=0)): 0;
      filterB = IS_REF_VALID(ref_frameB)? ((sub_posB == int_pos)? *FilterFlagInt: (p_Filter->FilterFlag[sub_posB]%2!=0)): 0;
      if(!filterF && !filterB)
        continue;


      // Bidirectional
      if(IS_REF_VALID(ref_frameF) && IS_REF_VALID(ref_frameB))
      {
        int mother_pelF = (sub_posF!=-1) ? p_Filter->FILTER_NEW_SUB_POS[sub_posF]:MAX_NUM_AIF;
        int mother_pelB = (sub_posB!=-1) ? p_Filter->FILTER_NEW_SUB_POS[sub_posB]:MAX_NUM_AIF;
        AvailFilter[mother_pelF] = 1;// it filter was used
        AvailFilter[mother_pelB] = 1;// it filter was used
        for(i = y; i < y + 4; ++i)
        {
          for(j = x; j < x + 4; ++j)
          {
            if(filterF)
              pixAifF = GetAIFPixel2(refYF, i, j, sub_posF, mvxF, mvyF, AIF_INTERPOLATION, filterID);
            else
              pixAifF = GetAIFPixel2(refYF, i, j, sub_posF, mvxF, mvyF, STD_INTERPOLATION, filterID);

            if(filterB)
              pixAifB = GetAIFPixel2(refYB, i, j, sub_posB, mvxB, mvyB, AIF_INTERPOLATION, filterID);
            else
              pixAifB = GetAIFPixel2(refYB, i, j, sub_posB, mvxB, mvyB, STD_INTERPOLATION, filterID);

            pixStdF = GetAIFPixel2(refYF, i, j, sub_posF, mvxF, mvyF, STD_INTERPOLATION, filterID);
            pixStdB = GetAIFPixel2(refYB, i, j, sub_posB, mvxB, mvyB, STD_INTERPOLATION, filterID);

            pix_00 = (1 - b_scaling) * pixStdF + b_scaling * pixStdB;
            pix_01 = (1 - b_scaling) * pixStdF + b_scaling * pixAifB;
            pix_10 = (1 - b_scaling) * pixAifF + b_scaling * pixStdB;
            pix_11 = (1 - b_scaling) * pixAifF + b_scaling * pixAifB;

            rec_00 = max(0, min(max_pixel_value, (int)(pix_00 + 0.5)));
            rec_01 = max(0, min(max_pixel_value, (int)(pix_01 + 0.5)));
            rec_10 = max(0, min(max_pixel_value, (int)(pix_10 + 0.5)));
            rec_11 = max(0, min(max_pixel_value, (int)(pix_11 + 0.5)));

            orig = (int)(imgY_org[i][j]);

            if(filterF && filterB)
            {
              // Both filters are valid
              C[0][0][mother_pelF][mother_pelB] += (double)abs(orig - rec_00);
              C[1][1][mother_pelF][mother_pelB] += (double)abs(orig - rec_11);

              if(sub_posF != sub_posB)
              {
                // Sub-pel positions are not using the same filter
                C[0][1][mother_pelF][mother_pelB] += (double)abs(orig - rec_01);
                C[1][0][mother_pelF][mother_pelB] += (double)abs(orig - rec_10);
              }
            } 
            else if(filterF)
            {
              // Filter F is valid
              C[0][0][mother_pelF][mother_pelB] += (double)abs(orig - rec_00);
              C[1][0][mother_pelF][mother_pelB] += (double)abs(orig - rec_10);
            }
            else
            {
              // Filter B is valid
              C[0][0][mother_pelF][mother_pelB] += (double)abs(orig - rec_00);
              C[0][1][mother_pelF][mother_pelB] += (double)abs(orig - rec_01);
            }
          }
        }
      }
      else
      {
        int mother_pel = (sub_posF!=-1) ? p_Filter->FILTER_NEW_SUB_POS[sub_posF]:MAX_NUM_AIF;
        AvailFilter[mother_pel] = 1;// it filter was used
        if(IS_REF_VALID(ref_frameB))
        {
          if(!filterB)
          {
            continue;
          }
          else
          {
            sub_posF = sub_posB;
            mvxF = mvxB; 
            mvyF = mvyB; 
            refYF = refYB;
          }
        }
        else if(!filterF)
        {
          continue;
        }

        for(i = y; i < y + 4; ++i)
        {
          for(j = x; j < x + 4; ++j)
          {            
            pix_00 = GetAIFPixel2(refYF, i, j, sub_posF, mvxF, mvyF, STD_INTERPOLATION,filterID);
            pix_11 = GetAIFPixel2(refYF, i, j, sub_posF, mvxF, mvyF, AIF_INTERPOLATION,filterID);

            rec_00 = max(0, min(max_pixel_value, (int)(pix_00 + 0.5)));
            rec_11 = max(0, min(max_pixel_value, (int)(pix_11 + 0.5)));

            orig   = (int)(imgY_org[i][j]);

            C[0][0][0][mother_pel] += (double)abs(orig - rec_00);
            C[1][1][0][mother_pel] += (double)abs(orig - rec_11);
          }
        }
      }     // Bidirectional if
    }     // Subblock loop
  }     // Macroblock loop


  {
    for(sub_pos = a_pos; sub_pos < MAX_NUM_AIF; sub_pos++)
    {
      if(AvailFilter[sub_pos])
      {
        costMCP[2*(filterID)][sub_pos] = C[0][0][0][sub_pos];
        costMCP[2*(filterID)+1][sub_pos] = C[1][1][0][sub_pos];
      }

    }
  }
  //printf("Done!\n");
}

#endif
#endif // E_DAIF
