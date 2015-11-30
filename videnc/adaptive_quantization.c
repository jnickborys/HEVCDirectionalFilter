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
#include "rdopt_coding_state.h"
#include "cabac.h"

#ifdef ADAPTIVE_QUANTIZATION

extern pic_parameter_set_rbsp_t *PicParSet[MAXPPS];
extern void frame_picture (Picture *frame, int method);
const int aqms_rd_pass_id = 4;

/*!
*************************************************************************************
* \brief
*    computing the same lambda for different QP values
*************************************************************************************
*/
void CalculateLambdaForMatrix(RD_PARAMS   *enc_mb)
{
  short bslice      = (img->type==B_SLICE);
  
  // Anchor QP setting
  if(input->UseAdaptiveQuantMatrix)
  {
    if (bslice && img->nal_reference_idc)
    {
      enc_mb->lambda_md = img->lambda_md[5][gdwBaseQp];
      enc_mb->lambda_me = img->lambda_me[5][gdwBaseQp];
      enc_mb->lambda_mf = img->lambda_mf[5][gdwBaseQp];
    }
    else
    {
      enc_mb->lambda_md = img->lambda_md[img->type][gdwBaseQp];
      enc_mb->lambda_me = img->lambda_me[img->type][gdwBaseQp];
      enc_mb->lambda_mf = img->lambda_mf[img->type][gdwBaseQp];
    }
  }
}

/*!
***********************************************************************
* \brief
*   storing the input parameters (for speed up)
***********************************************************************
*/
void SaveInputParam(int AQMS_Skip_flag, sDefaultInputParam_t *defInputParam)
{
  defInputParam->default_me_modeflag  =input->FMEnable;
  defInputParam->default_ChromaIntraDisable = input->ChromaIntraDisable;
  defInputParam->default_SelectiveIntraEnable = input->SelectiveIntraEnable;
  defInputParam->default_EarlySkipEnable = input->EarlySkipEnable;
  defInputParam->default_IAQMS_flag = img->slice_fractional_quant_flag;

  // storing the additional command options for speed up encoding.
  // Please add the new command option if you make.
#ifdef ADAPTIVE_FD_SD_CODING
  defInputParam->default_APEC_in_FD_and_SD = input->APEC_in_FD_and_SD;
#endif
#ifdef MV_COMPETITION
  defInputParam->default_MVCompetition = input->mv_competition;
#endif

  // swtiching off the additional command options for speed up encoding.
  // Please add the new command option if you make.
  img->slice_fractional_quant_flag = 0;
#ifdef ADAPTIVE_FD_SD_CODING
  input->APEC_in_FD_and_SD = 0;
#endif
#ifdef MV_COMPETITION
  input->mv_competition = 0;
#endif
}
/*!
***********************************************************************
* \brief
*   loading the input parameters (for speed up)
***********************************************************************
*/
void LoadInputParam(sDefaultInputParam_t *defInputParam)
{
  input->FMEnable             = defInputParam->default_me_modeflag;
  input->ChromaIntraDisable = defInputParam->default_ChromaIntraDisable;
  input->SelectiveIntraEnable = defInputParam->default_SelectiveIntraEnable;
  input->EarlySkipEnable = defInputParam->default_EarlySkipEnable;
  img->slice_fractional_quant_flag = defInputParam->default_IAQMS_flag;

  // loading the additional command option for speed up encoding.
  // Please add the new command option if you make.
#ifdef ADAPTIVE_FD_SD_CODING
  input->APEC_in_FD_and_SD    = defInputParam->default_APEC_in_FD_and_SD;
#endif
#ifdef MV_COMPETITION
  input->mv_competition       = defInputParam->default_MVCompetition;
#endif
}

/*!
***********************************************************************
* \brief
*   checking RD cost for previous flat quantization matrix
***********************************************************************
*/
int ChooseBestPicture(int rd_qp)
{
  int mincost_flg=0;
  img->write_macroblock = 0;
  
#ifdef ADAPTIVE_FILTER
  if(input->RDPictureDecision && img->rd_pass == 3)
  {
    mincost_flg=picture_coding_decision_for_AQMS(frame_pic_aqms, frame_pic_aif, rd_qp);
    if(mincost_flg)
    {
      active_pps = PicParSet[active_pps->pic_parameter_set_id];
      enc_picture = enc_frame_picture_aif;
      frame_pic   = frame_pic_aif;
    }
  }
  else
#endif
    if(input->RDPictureDecision && img->rd_pass == 2)
    {
      mincost_flg=picture_coding_decision_for_AQMS(frame_pic_aqms, frame_pic_3, rd_qp);
      if(mincost_flg)
      {
        active_pps = PicParSet[active_pps->pic_parameter_set_id];
        enc_picture = enc_frame_picture3;
        frame_pic   = frame_pic_3;
      }
    }
    else if(input->RDPictureDecision && img->rd_pass == 1)
    {
      mincost_flg=picture_coding_decision_for_AQMS(frame_pic_aqms, frame_pic_2, rd_qp);
      if(mincost_flg)
      {
        active_pps = PicParSet[active_pps->pic_parameter_set_id];
        enc_picture = enc_frame_picture2;
        frame_pic   = frame_pic_2;
      }
    }
    else
    {
      mincost_flg=picture_coding_decision_for_AQMS(frame_pic_aqms, frame_pic_1, rd_qp);
      if(mincost_flg)
      {
        active_pps = PicParSet[active_pps->pic_parameter_set_id];
        enc_picture = enc_frame_picture;
        frame_pic   = frame_pic_1;
      }
    }
    
    return mincost_flg;
}
#ifdef MV_COMPETITION
/*! 
*************************************************************************************
* \brief
*   reencoding with mv_competition
*
* \para 
*
*************************************************************************************
*/
void ReEncodeMVCompetition(int rd_pass)
{
  img->slice_fractional_quant_flag = 0;

  if(rd_pass==2)
  {
    if(enc_frame_picture3)
    {
      free_storable_picture (enc_frame_picture3);
      enc_frame_picture3 = NULL;
    }
    if (frame_pic_3)
    {
      free_slice_list(frame_pic_3);
    }
    frame_picture (frame_pic_3,img->rd_pass);
  }
  else if(rd_pass==1)
  {
    if(enc_frame_picture2)
    {
      free_storable_picture (enc_frame_picture2);
      enc_frame_picture2 = NULL;
    }
    if (frame_pic_2)
    {
      free_slice_list(frame_pic_2);
    }
    frame_picture (frame_pic_2,img->rd_pass);
  }
  else
  {
    if(enc_frame_picture)
    {
      free_storable_picture (enc_frame_picture);
      enc_frame_picture = NULL;
    }
    if (frame_pic_1)
    {
      free_slice_list(frame_pic_1);
    }
    frame_picture (frame_pic_1,img->rd_pass);
  }
}
#endif

/*!
 ***********************************************************************
 * \brief
 *   fractional quantization
 ***********************************************************************
 */
int rdPictureCodingForIAQMS(int FAQ_Skip_flag, sDefaultInputParam_t *defInputParam, int *mincost_flg, int rd_qp)
{
  int idr_flag = ((!IMG_NUMBER) && (!(img->structure==BOTTOM_FIELD))) || (input->idr_enable && (img->type == I_SLICE || img->type==SI_SLICE)&& (!(img->structure==BOTTOM_FIELD)));
  best_rd_pass = img->rd_pass;

  if( (input->UseAdaptiveQuantMatrix && img->type == I_SLICE && !idr_flag) || (FAQ_Skip_flag))
  {
    *mincost_flg = 1;
    return 0;
  }

  // set slice modeling quantization matrix
  img->slice_modeling_qm_param0 = VALUE_OF_QM_PARAM0;
  img->slice_modeling_qm_param1 = VALUE_OF_QM_PARAM1;
  SetAutoScalingListDataForMQM();
  // set macroblock quantization matrix
  SetAutoScalingListDataForIAQMS();

  ReGenerateParameterSets(active_pps->pic_parameter_set_id);
  active_pps  = PicParSet[PPSID_AQMS];
  // return default input parameter
  LoadInputParam(defInputParam);
  SetFastEncodingForIAQMS(defInputParam);

  img->write_macroblock = 0;
  {
    enc_picture = enc_frame_picture_aqms;
    frame_pic   = frame_pic_aqms;
  }

  // Final encode
  frame_picture(frame_pic_aqms, gdwAQMS_frame_num);
  *mincost_flg = ChooseBestPicture(rd_qp);

  // return default input parameter
  LoadInputParam(defInputParam);

  if(*mincost_flg==0)
    img->rd_pass = gdwAQMS_frame_num;
  else
    img->rd_pass = best_rd_pass;

  SkipRatio[(img->type)+(img->type==B_SLICE? (img->nal_reference_idc<<2) : 0)] = ((double)dwSkipMbCountForSlice[(img->type)+(img->type==B_SLICE? (img->nal_reference_idc<<2) : 0)]/(double)dwSkipEnableMbCountForSlice[(img->type)+(img->type==B_SLICE? (img->nal_reference_idc<<2) : 0)]*100.0);

#ifdef MV_COMPETITION // MVC bug fix by AT
  if(*mincost_flg && input->mv_competition)
    ReEncodeMVCompetition(img->rd_pass);
#endif

  return 0;
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

  if(Gradient==255) // Flat base value minus 1
  {
    for(i=0;i<16;i++)
      for(j=0;j<6;j++)
        ScalingList4x4_IAQMS[ListIdx][j][i] = Clip3(1,255,(ScalingList4x4[j][i]-1));
  }
  else if(Gradient==256)  // Flat base value plus 1
  {
    for(i=0;i<16;i++)
      for(j=0;j<6;j++)
        ScalingList4x4_IAQMS[ListIdx][j][i] = Clip3(1,255,(ScalingList4x4[j][i]+1));
  }
  else  // slope matrix
  {
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
  int QM8[64];

  if(Gradient==255) // Flat base value minus 1
  {
    for(i=0;i<64;i++)
      for(j=0;j<2;j++)
        ScalingList8x8_IAQMS[ListIdx][j][i] = Clip3(1,255,(ScalingList8x8[j][i]-1));
  }
  else if(Gradient==256)  // Flat base value plus 1
  {
    for(i=0;i<64;i++)
      for(j=0;j<2;j++)
        ScalingList8x8_IAQMS[ListIdx][j][i] = Clip3(1,255,(ScalingList8x8[j][i]+1));
  }
  else  // slope matrix
  {
    for(j=0;j<8;j++)
      for(i=0;i<8;i++)
      {
        r=i+j;
        QM8[i+(j<<3)]=(((Gradient*r)+16)>>5);
      }

    for(i=0;i<64;i++)
      for(j=0;j<2;j++)
        ScalingList8x8_IAQMS[ListIdx][j][i] = Clip3(1,255,(QM8[i]+ScalingList8x8[j][i]));
  }
}

/*!
************************************************************************
* \brief
*    For generating the modeling quantization matrix at slice level
*
************************************************************************
*/
void SetAutoScalingListDataForMQM()
{
  int pA = 0, pB = 0;
  
  pA = ((img->slice_modeling_qm_param0)-32);
  pB = (img->slice_modeling_qm_param1);
  
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
  unsigned int QM8[64];
  
  // 8x8 Matrix
  for(j=0;j<8;j++)
  {
    for(i=0;i<8;i++)
    {
      r=i+j;
      QM8[i+(j<<3)]=Clip3(1,255,( ((pA*r+16)>>5 ) + pB) );
    }
  }
  
  for(i=0;i<64;i++)
    for(j=0;j<2;j++)
      ScalingList8x8[j][i] = Clip3(1,255,QM8[i]);
}

/*!
 ***********************************************************************
 * \brief
 *   storing the input parameters (for speed up)
 ***********************************************************************
 */
void SetFastEncodingForIAQMS(sDefaultInputParam_t *defInputParam)
{
  // for fast encoding
  input->SelectiveIntraEnable = 1;
  input->EarlySkipEnable = 1;
  input->FMEnable = 3;

  if(img->type!=I_SLICE)
  {
    input->ChromaIntraDisable     = 1;
  }
}
/*
 ***********************************************************************
 * \brief
 *   initialize scaling lists 
 ***********************************************************************
 */
void InitScalingList()
{
  int i, j;
  for(j=0;j<6;j++)
    for(i=0;i<16;i++)
      ScalingList4x4[j][i] = 16;

  for(j=0;j<2;j++)
    for(i=0;i<64;i++)
      ScalingList8x8[j][i] = 16;
}

/*
 ***********************************************************************
 * \brief
 *   set mode availability
 ***********************************************************************
 */
void SetModeAvalability(RD_PARAMS *enc_mb)
{
  int i;

  if(final_mb_encoding == 1)
  {
    if(saved_best_mode<=P8x8)
    {
      for(i=P8x8+1; i<MAXMODE; i++)
        enc_mb->valid[i] = 0;
    }
  }
}

#endif
