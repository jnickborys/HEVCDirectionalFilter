/*!
************************************************************************
* \file adaptive_quantization.c
*
* \brief
*    adaptive quantization matrix selection routines
*
* \author
*  Akiyuki Tanizawa       <akiyuki.tanizawa@toshiba.co.jp>           \n
*  Takeshi Chujoh         <takeshi.chujoh@toshiba.co.jp>             \n
*
* This software may be used for research purposes only.
* TOSHIBA would also appreciate that any technical report, conference or
* journal paper which uses the software includes one of the following references:
*    - A. Tanizawa, T. Chujoh, "Adaptive Quantization Matrix Selection",
*                              VCEG contribution D.266, Geneva, November 2006.
*    - A. Tanizawa, T. Chujoh, "Simulation results of AQMS on KTA software",
*                              VCEG contribution VCEG-AC07, Klagenfurt, July 2006.
*    - A. Tanizawa, T. Chujoh, "Adaptive Quantization Matrix Selection on KTA software",
*                              VCEG contribution VCEG-AD06, Hangzhou, Oct. 2006.
*    - A. Tanizawa, T. Chujoh, "Simulation results of AQMS on KTA software version 1.3",
*                              VCEG contribution VCEG-AF08, San Jose, April 2007.
*    - A. Tanizawa, T. Chujoh, T. Yamakage, 
*                              "Encoding complexity reduction of AQMS",
*                              VCEG contribution C.403, Geneve, April 2008.
*    - A. Tanizawa, T. Chujoh, T. Yamakage, 
*                              "Improvement of Adaptive Quantization Matrix Selection",
*                              VCEG contribution VCEG-AI19, Berlin, July 2008.
************************************************************************
*/

/*
* CONFIGURATION:
*
* Configuration of the .Cfg file:
* In section KTA STUFF:
*
* // AQMS command:
* UseAdaptiveQuantMatrix    = 0 (0:disable, 1:enable (IAQMS))
*
* NOTE:
* IAQMS can perform on VCEG-AA10 conditions.
*    - TK. Tan, G. Sullivan, T. Wedi, "Recommended Simulation Common Conditions for Coding Efficiency Experiments",
*                              VCEG contribution VCEG-AA10, Nice, Oct. 2005.
*    - TK. Tan, G. Sullivan, T. Wedi, "Recommended Simulation Common Conditions for Coding Efficiency Experiments",
*                              VCEG contribution VCEG-AE10r1, Marrakech, Jan. 2007.
********************
* KNOWN LIMITATIONS
********************
* Not compatible with:                
*  - interlace coding
*  - RdOpt != 1   
*  - RateControl
*
********************
* VERSIONS
********************
* 1.0 KTA version 1.2  - First Implementation for AQMS
* 1.1 KTA version 1.3  - Modification of Compatibility with AIF, APEC, 
* 1.2 KTA version 1.9  - Addition of Fast AQMS and Modification of Compatibility with RDOQ, 
* 1.3 KTA version 2.1  - Replacement of IAQMS (Improved AQMS) and Deletion of Adaptive QP
*
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include "global.h"
#include "parset.h"
#include "adaptive_quantization.h"
#include "memalloc.h"
#include "malloc.h"
#include "defines.h"
#include "mbuffer.h"
#include "image.h"

#ifdef ADAPTIVE_QUANTIZATION

extern const int dequant_coef[6][4][4];
static const int dequant_coef8[6][8][8] = 
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

/*!
 ************************************************************************
 * \brief
 *    For assigning the modeling quantization matrix
 *
 ************************************************************************
 */
void AssignQuantParamForAQMS(pic_parameter_set_rbsp_t* pps)
{
  int i;

  if( img->slice_fractional_quant_flag )
  {
    int j;
    for(i=0; i<8; i++)
    {
      if(i<6)
        qmatrix[i] = pps->ScalingList4x4[i];
      else
        qmatrix[i] = pps->ScalingList8x8[i-6];
    }

    for(j=0;j<NUM_OF_IAQMS_MAT;j++)
    {
      for(i=0; i<8; i++)
      {
        qmatrixIAQMS[j][i] = (i<6) ? pps->ScalingList4x4_IAQMS[j][i] : pps->ScalingList8x8_IAQMS[j][i-6];
      }
    }
  }
}
/*!
 ************************************************************************
 * \brief
 *    For calculating inverse-level scale for AQMS
 *
 ************************************************************************
 */
void CalculateQuantParamForAQMS()
{
  int i,j,k,temp;

  if(img->slice_fractional_quant_flag)
  {
    int l;
    for(l=0; l<NUM_OF_IAQMS_MAT; l++)
      for(k=0; k<6; k++)
        for(j=0; j<4; j++)
          for(i=0; i<4; i++)
          {
            temp = (i<<2)+j;
            InvLevelScale4x4Luma_Intra_IAQMS[l][k][j][i]      = dequant_coef[k][j][i]*qmatrixIAQMS[l][0][temp];
            InvLevelScale4x4Chroma_Intra_IAQMS[l][0][k][j][i] = dequant_coef[k][j][i]*qmatrixIAQMS[l][1][temp];
            InvLevelScale4x4Chroma_Intra_IAQMS[l][1][k][j][i] = dequant_coef[k][j][i]*qmatrixIAQMS[l][2][temp];

            InvLevelScale4x4Luma_Inter_IAQMS[l][k][j][i]      = dequant_coef[k][j][i]*qmatrixIAQMS[l][3][temp];
            InvLevelScale4x4Chroma_Inter_IAQMS[l][0][k][j][i] = dequant_coef[k][j][i]*qmatrixIAQMS[l][4][temp];
            InvLevelScale4x4Chroma_Inter_IAQMS[l][1][k][j][i] = dequant_coef[k][j][i]*qmatrixIAQMS[l][5][temp];
          }
  }
}

/*!
 ************************************************************************
 * \brief
 *    For calculating the inverse level scale
 *
 ************************************************************************
 */
void CalculateQuant8ParamForAQMS()
{
  int i,j,k,temp;

  if(img->slice_fractional_quant_flag)
  {
    int l;
    for(l=0; l<NUM_OF_IAQMS_MAT;l++)
      for(k=0; k<6; k++)
        for(j=0; j<8; j++)
          for(i=0; i<8; i++)
          {
            temp = (i<<3)+j;
            InvLevelScale8x8Luma_Intra_IAQMS[l][k][j][i] = dequant_coef8[k][j][i]*qmatrixIAQMS[l][6][temp];
            InvLevelScale8x8Luma_Inter_IAQMS[l][k][j][i] = dequant_coef8[k][j][i]*qmatrixIAQMS[l][7][temp];
          }
  }
}

/*!
 ************************************************************************
 * \brief
 *    For generating the modeling quantization matrix
 *
 ************************************************************************
 */
void SetAutoScalingListDataForIAQMS()
{
  int ListIdx=0;
  int *GradientTab;
  int cGradientTab_HR[6][NUM_OF_IAQMS_MAT] ={{  0,  -8,256,255}, // P-slice
                                            {  0, 256,256,256}, // B-slice
                                            {  0,  -8,-16,255}, // I-slice
                                            {  0,   4, -4,  8},
                                            {  0,   4, -4,  8},
                                            {  0,   4, -4,  8}};
  int cGradientTab_BR[6][NUM_OF_IAQMS_MAT] ={{  0,256,256,256}, // P-slice
                                            {  0, 256,256,256}, // B-slice
                                            {  0,  -8,-16,255}, // I-slice
                                            {  0,   4, -4,  8},
                                            {  0,   4, -4,  8},
                                            {  0,   4, -4,  8}};

  for(ListIdx=0;ListIdx<NUM_OF_IAQMS_MAT;ListIdx++)
  {
    if(active_sps->profile_idc < FREXT_HP)
      GradientTab = cGradientTab_BR[img->type];
    else
      GradientTab = cGradientTab_HR[img->type];

    // 4x4 block
    CreateScalingList4x4ForIAQMS(GradientTab[ListIdx], ListIdx);
    // 8x8 block
    CreateScalingList8x8ForIAQMS(GradientTab[ListIdx], ListIdx);
  }
}
/*!
 ************************************************************************
 * \brief
 *    For generating the modeling quantization matrix at 4x4 block
 *
 ************************************************************************
 */
void CreateScalingList4x4ForIAQMS(int Gradient, int ListIdx)
{
  int i,j,r=0;
  int QM4[16];

  if(Gradient==255)
  {
    for(i=0;i<16;i++)
      for(j=0;j<6;j++)
        ScalingList4x4_IAQMS[ListIdx][j][i] = Clip3(1,255,(ScalingList4x4[j][i]-1));
  }
  else if(Gradient==256)
  {
    for(i=0;i<16;i++)
      for(j=0;j<6;j++)
        ScalingList4x4_IAQMS[ListIdx][j][i] = Clip3(1,255,(ScalingList4x4[j][i]+1));
  }
  else
  {
    // 4x4 Matrix
    for(j=0;j<4;j++)
      for(i=0;i<4;i++)
      {
        r=i+j;
        QM4[i+(j<<2)]=(((Gradient*r)+8)>>4);
      }

    for(i=0;i<16;i++)
      for(j=0;j<6;j++)
        ScalingList4x4_IAQMS[ListIdx][j][i] = Clip3(1,255,(QM4[i]+ScalingList4x4[j][i]));
  }
}
/*!
 ************************************************************************
 * \brief
 *    For generating the modeling quantization matrix at 8x8 block
 *
 ************************************************************************
 */
void CreateScalingList8x8ForIAQMS(int Gradient, int ListIdx)
{
  int i,j,r=0;
  int QM[64];

  if(Gradient==255)
  {
    for(i=0;i<64;i++)
      for(j=0;j<2;j++)
        ScalingList8x8_IAQMS[ListIdx][j][i] = Clip3(1,255,(ScalingList8x8[j][i]-1));
  }
  else if(Gradient==256)
  {
    for(i=0;i<64;i++)
      for(j=0;j<2;j++)
        ScalingList8x8_IAQMS[ListIdx][j][i] = Clip3(1,255,(ScalingList8x8[j][i]+1));
  }
  else
  {
    // 8x8 Matrix
    for(j=0;j<8;j++)
      for(i=0;i<8;i++)
      {
        r=i+j;
        QM[i+(j<<3)]=(((Gradient*r)+16)>>5);
      }

    for(i=0;i<64;i++)
      for(j=0;j<2;j++)
        ScalingList8x8_IAQMS[ListIdx][j][i] = Clip3(1,255,(QM[i]+ScalingList8x8[j][i]));
  }
}

/*!
 ************************************************************************
 * \brief
 *    For generating the modeling quantization matrix
 *
 ************************************************************************
 */
void CreateModelingMatrixForIAQMS()
{
  int i,j;

  for(i=0;i<6;i++)
    ScalingList4x4[i] = active_pps->ScalingList4x4[i];
  for(i=0;i<6;i++)
    ScalingList8x8[i] = active_pps->ScalingList8x8[i];

  for(j=0;j<NUM_OF_IAQMS_MAT;j++)
  {
    for(i=0;i<6;i++)
      ScalingList4x4_IAQMS[j][i] = active_pps->ScalingList4x4_IAQMS[j][i];
    for(i=0;i<6;i++)
      ScalingList8x8_IAQMS[j][i] = active_pps->ScalingList8x8_IAQMS[j][i];
  }

  SetAutoScalingListDataForIAQMS();
}
/*!
************************************************************************
* \brief
*    For generating the modeling quantization matrix
*
************************************************************************
*/
void SetAutoScalingListDataForMQM(int param0, int param1)
{
  
  int i, j,  pA = 0, pB = 0;
  
  for(i=0;i<6;i++)
    ScalingList4x4[i] = active_pps->ScalingList4x4[i];
  for(i=0;i<6;i++)
    ScalingList8x8[i] = active_pps->ScalingList8x8[i];

  for(j=0;j<NUM_OF_IAQMS_MAT;j++)
  {
    for(i=0;i<6;i++)
      ScalingList4x4_IAQMS[j][i] = active_pps->ScalingList4x4_IAQMS[j][i];
    for(i=0;i<6;i++)
      ScalingList8x8_IAQMS[j][i] = active_pps->ScalingList8x8_IAQMS[j][i];
  }

  pA = (param0-32);
  pB = (param1);
  
  // 4x4 block
  CreateScalingList4x4ForMQM(pA, pB);
  // 8x8 block
  CreateScalingList8x8ForMQM(pA, pB);
}
/*!
************************************************************************
* \brief
*    For generating the modeling quantization matrix at 4x4 block
*
************************************************************************
*/
void CreateScalingList4x4ForMQM(int pA, int pB)
{
  int i,j,r=0;
  unsigned int QM4[16];
  
  // 4x4 Matrix
  for(j=0;j<4;j++)
  {
    for(i=0;i<4;i++)
    {
      r=i+j;
      QM4[i+(j<<2)]=Clip3(1,255,(((pA*r)+8)>>4) + pB);
    }
  }
  
  for(i=0;i<16;i++)
    for(j=0;j<6;j++)
      ScalingList4x4[j][i] = Clip3(1,255,QM4[i]);
}
/*!
************************************************************************
* \brief
*    For generating the modeling quantization matrix at 8x8 block
*
************************************************************************
*/
void CreateScalingList8x8ForMQM(int pA, int pB)
{
  int i,j,r=0;
  unsigned int QM[64];
  
  // 8x8 Matrix
  for(j=0;j<8;j++)
  {
    for(i=0;i<8;i++)
    {
      r=i+j;
      QM[i+j*8]=Clip3(1,255,( ((pA*r+16)>>5 ) + pB) );
    }
  }
  
  for(i=0;i<64;i++)
    for(j=0;j<2;j++)
      ScalingList8x8[j][i] = Clip3(1,255,QM[i]);
}
/*!
 ************************************************************************
 * \brief
 *    initialize scalinglists for SPS
 *
 ************************************************************************
 */
void InitScalingListSPS(seq_parameter_set_rbsp_t *sps)
{
  int i, j;
  for(j=0;j<6;j++)
    for(i=0;i<16;i++)
      sps->ScalingList4x4[j][i] = 16;

  for(j=0;j<2;j++)
    for(i=0;i<64;i++)
      sps->ScalingList8x8[j][i] = 16;
}
/*!
 ************************************************************************
 * \brief
 *    initialize scalinglists for PPS
 *
 ************************************************************************
 */
void InitScalingListPPS(pic_parameter_set_rbsp_t *pps)
{
  int i, j;
  for(j=0;j<6;j++)
    for(i=0;i<16;i++)
      pps->ScalingList4x4[j][i] = 16;

  for(j=0;j<2;j++)
    for(i=0;i<64;i++)
      pps->ScalingList8x8[j][i] = 16;
}

#endif
