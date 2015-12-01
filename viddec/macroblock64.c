
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "global.h"
#include "image.h"
#include "elements.h"
#include "errorconcealment.h"
//#include "macroblock.h"
#include "fmo.h"
#include "cabac.h"
#include "vlc.h"
#include "image.h"
#include "mb_access.h"
#include "biaridecod.h"
#include "transform8x8.h"

#ifdef MV_COMPETITION
#include "mv_competition.h"
extern MV_Competition mv_comp;                  // to each possible predictor
#endif

#ifdef ADAPTIVE_FILTER
#include "adaptive_filter.h"
#endif

#ifdef MB32X32
#define DCT16PREC   7

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#include <math.h>
#if TRACE
#define TRACE_STRING(s) strncpy(currSE.tracestring, s, TRACESTRING_SIZE)
#else
#define TRACE_STRING(s) // do nothing
#endif


int mb32_mvd[2][BLOCK_MULTIPLE<<MAX_MB_EXT_LEVEL][BLOCK_MULTIPLE<<MAX_MB_EXT_LEVEL][2];          //!< indices correspond to [forw,backw][block_y][block_x][x,y]
extern ColocatedParams *Co_located;
extern const int BLOCK_STEP[8][2];
extern const byte QP_SCALE_CR[52];

extern int last_dquant;

extern int delta_qp_sent;
extern Macroblock  MB32;
extern Macroblock  MB64;//MB64X64

int mb_is_available(int mbAddr, int currMbAddr);
void readMB_typeInfo_CABAC32(SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readCBP_CABAC32(SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readRefFrame_CABAC32( SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readMVD_CABAC32( SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void interpret_mb_mode_P(struct img_par *img);
void interpret_mb_mode_I(struct img_par *img);
void interpret_mb_mode_SI(struct img_par *img);
void interpret_mb_mode_B(struct img_par *img);
void SetB8Mode (struct img_par* img, Macroblock* currMB, int value, int i);
void init_macroblock(struct img_par *img);
void field_flag_inference();
int BType2CtxRef (int btype);
void reset_coeffs();
void set_chroma_qp(Macroblock* currMB);

int Not_MB16_3(int mb_nr)
{
  int mb_nr_x, mb_nr_y;
  mb_nr_x = mb_nr % img->PicWidthInMbs;
  mb_nr_y = mb_nr / img->PicWidthInMbs;

  if(mb_nr_x%2==1 && mb_nr_y%2==1)
  {
    if((unsigned)mb_nr_x < img->PicWidthInMbs-1)
    {
      return 0;
    }
  }
  return 1;
}

/*!
 ************************************************************************
 * \brief
 *    Checks the availability of neighboring macroblocks of
 *    the current macroblock for prediction and context determination;
 ************************************************************************
 */
void CheckAvailabilityOfNeighbors_Z()
{
  const int mb_nr = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];

  // mark all neighbors as unavailable
  currMB->mb_available_up   = NULL;
  currMB->mb_available_left = NULL;

  if (dec_picture->MbaffFrameFlag)
  {
    currMB->mbAddrA = 2 * (mb_nr/2 - 1);
    currMB->mbAddrB = 2 * (mb_nr/2 - dec_picture->PicWidthInMbs);
    currMB->mbAddrC = 2 * (mb_nr/2 - dec_picture->PicWidthInMbs + 1);
    currMB->mbAddrD = 2 * (mb_nr/2 - dec_picture->PicWidthInMbs - 1);
    
    currMB->mbAvailA = mb_is_available(currMB->mbAddrA, mb_nr) && (((mb_nr/2) % dec_picture->PicWidthInMbs)!=0);
    currMB->mbAvailB = mb_is_available(currMB->mbAddrB, mb_nr);
    currMB->mbAvailC = mb_is_available(currMB->mbAddrC, mb_nr) && (((mb_nr/2 +1) % dec_picture->PicWidthInMbs)!=0);
    currMB->mbAvailD = mb_is_available(currMB->mbAddrD, mb_nr) && (((mb_nr/2) % dec_picture->PicWidthInMbs)!=0);
  }
  else
  {
    currMB->mbAddrA = mb_nr - 1;
    currMB->mbAddrB = mb_nr - dec_picture->PicWidthInMbs;
    currMB->mbAddrC = mb_nr - dec_picture->PicWidthInMbs + 1;
    currMB->mbAddrD = mb_nr - dec_picture->PicWidthInMbs - 1;

    currMB->mbAvailA = mb_is_available(currMB->mbAddrA, mb_nr) && ((mb_nr % dec_picture->PicWidthInMbs)!=0);
    currMB->mbAvailB = mb_is_available(currMB->mbAddrB, mb_nr);
    currMB->mbAvailC =  Not_MB16_3(mb_nr) && mb_is_available(currMB->mbAddrC, mb_nr) && (((mb_nr+1) % dec_picture->PicWidthInMbs)!=0);
    currMB->mbAvailD = mb_is_available(currMB->mbAddrD, mb_nr) && ((mb_nr % dec_picture->PicWidthInMbs)!=0);
  }
}

int readCBP32(struct img_par *img,struct inp_par *inp, int dquant)
{
  //int i,j,k;
  //int level;
  int mb_nr = img->current_mb_nr;
  //int ii,jj;
  //int m2,jg2;// i1,j1;
  Macroblock *currMB = &img->mb_data[mb_nr];
  int cbp;
  SyntaxElement currSE;
  Slice *currSlice = img->currentSlice;
  DataPartition *dP;
  int *partMap = assignSE2partition[currSlice->dp_mode];
  //int coef_ctr, i0, j0, b8;
  //int ll;
  //int block_x,block_y;
  //int start_scan;
  //int run, len;
  //int levarr[16], runarr[16], numcoeff;
  //
  //int qp_const;
  //int qp_per    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)/6;
  //int qp_rem    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)%6;
  //int smb       = ((img->type==SP_SLICE) && IS_INTER (currMB)) || (img->type == SI_SLICE && currMB->mb_type == SI4MB);
  //
  //int uv;
  //int qp_const_uv[2];
  //int qp_per_uv[2];
  //int qp_rem_uv[2];
  //
  //int intra     = IS_INTRA (currMB);
  //int temp[4];
  //
  //int b4;
  //int yuv = dec_picture->chroma_format_idc-1;
  //int m5[4];
  //int m6[4];
  
  //int need_transform_size_flag;
  //Boolean lossless_qpprime = (Boolean)((img->qp + img->bitdepth_luma_qp_scale)==0 && img->lossless_qpprime_flag==1);
  
  // Residue Color Transform
  //Boolean residual_transform_dc = (Boolean)((img->residue_transform_flag==1) && (IS_OLDINTRA(currMB)||currMB->mb_type==I8MB) );
//#ifdef USE_INTRA_MDDT
//  int R[6] = {40, 45, 50, 57, 63, 71};
//#endif
  
  //if(img->type==SP_SLICE  && currMB->mb_type!=I16MB )
  //  smb=1;
  //
  // QPI
  //init constants for every chroma qp offset
  //if (dec_picture->chroma_format_idc != YUV400)
  //{
  //  for (i=0; i<2; i++)
  //  {
  //    qp_per_uv[i] = (currMB->qpc[i] + img->bitdepth_chroma_qp_scale)/6;
  //    qp_rem_uv[i] = (currMB->qpc[i] + img->bitdepth_chroma_qp_scale)%6;
  //  }
  //}
  
//#ifdef ADAPTIVE_FD_SD_CODING
//  currMB->SD_Coding_on_off=0;
//#endif
  
  // read CBP if not new intra mode
  //if (!IS_NEWINTRA (currMB))
  {
    //=====   C B P   =====
    //---------------------
    if (IS_OLDINTRA (currMB) || currMB->mb_type == SI4MB || currMB->mb_type == I8MB)   currSE.type = SE_CBP_INTRA;
    else                        currSE.type = SE_CBP_INTER;
    
    dP = &(currSlice->partArr[partMap[currSE.type]]);
    
    if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
    {
      if (IS_OLDINTRA (currMB) || currMB->mb_type == SI4MB || currMB->mb_type == I8MB)  currSE.mapping = linfo_cbp_intra;
      else                        currSE.mapping = linfo_cbp_inter;
    }
    else
    {
      currSE.reading = readCBP_CABAC32;
    }
    
    TRACE_STRING("coded_block_pattern");
    
    dP->readSyntaxElement(&currSE,img,inp,dP);
    currMB->cbp = cbp = currSE.value1;
 
    if (cbp !=0 && dquant)
    {
      if (IS_INTER (currMB))  currSE.type = SE_DELTA_QUANT_INTER;
      else                    currSE.type = SE_DELTA_QUANT_INTRA;
      
      dP = &(currSlice->partArr[partMap[currSE.type]]);
      
      if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
      {
        currSE.mapping = linfo_se;
      }
      else
        currSE.reading= readDquant_CABAC; //gabi
      
      TRACE_STRING("mb_qp_delta");
      dP->readSyntaxElement(&currSE,img,inp,dP);
      currMB->delta_quant = currSE.value1;
      if ((currMB->delta_quant < -(26 + img->bitdepth_luma_qp_scale/2)) || (currMB->delta_quant > (25 + img->bitdepth_luma_qp_scale/2)))
        error ("mb_qp_delta is out of range", 500);
      
      img->qp= ((img->qp + currMB->delta_quant + 52 + 2*img->bitdepth_luma_qp_scale)%(52+img->bitdepth_luma_qp_scale)) -
        img->bitdepth_luma_qp_scale;
    }
 
    return cbp;   
  }
}
/*!
 ************************************************************************
 * \brief
 *    Checks the availability of neighboring macroblocks of
 *    the current macroblock for prediction and context determination;
 *    mb_ext_level: 2-this is 64x64 block; 1-this is 32x32 block; 0-a 16x16 block inside a 32x32 block
 ************************************************************************
 */
void CheckAvailabilityOfNeighbors32(int mb_ext_level)
{
  const int mb_nr = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];

  // mark all neighbors as unavailable
  currMB->mb_available_up   = NULL;
  currMB->mb_available_left = NULL;

  if (img->MbaffFrameFlag)
  {
    currMB->mbAddrA = 2 * (mb_nr/2 - 1);
    currMB->mbAddrB = 2 * (mb_nr/2 - img->PicWidthInMbs);
    currMB->mbAddrC = 2 * (mb_nr/2 - img->PicWidthInMbs + 1);
    currMB->mbAddrD = 2 * (mb_nr/2 - img->PicWidthInMbs - 1);
    
    currMB->mbAvailA = mb_is_available(currMB->mbAddrA, mb_nr) && (((mb_nr/2) % img->PicWidthInMbs)!=0);
    currMB->mbAvailB = mb_is_available(currMB->mbAddrB, mb_nr);
    currMB->mbAvailC = mb_is_available(currMB->mbAddrC, mb_nr) && (((mb_nr/2 +1) % img->PicWidthInMbs)!=0);
    currMB->mbAvailD = mb_is_available(currMB->mbAddrD, mb_nr) && (((mb_nr/2) % img->PicWidthInMbs)!=0);
  }
  else
  {
    currMB->mbAddrA = mb_nr - 1;
    currMB->mbAddrB = mb_nr - img->PicWidthInMbs;
    currMB->mbAddrC = mb_nr - img->PicWidthInMbs + (1<<mb_ext_level);
    currMB->mbAddrD = mb_nr - img->PicWidthInMbs - 1;

    currMB->mbAvailA = mb_is_available(currMB->mbAddrA, mb_nr) && ((mb_nr % img->PicWidthInMbs)!=0);
    currMB->mbAvailB = mb_is_available(currMB->mbAddrB, mb_nr);
    if(img->type != I_SLICE)
    {
      if(mb_ext_level)
        currMB->mbAvailC = mb_is_available(currMB->mbAddrC, mb_nr) && (((mb_nr+(1<<mb_ext_level)) % img->PicWidthInMbs)!=0);
      else
        currMB->mbAvailC = Not_MB16_3(mb_nr) && mb_is_available(currMB->mbAddrC, mb_nr) && (((mb_nr+1) % img->PicWidthInMbs)!=0);
    }
    else
    currMB->mbAvailC = mb_is_available(currMB->mbAddrC, mb_nr) && (((mb_nr+1) % img->PicWidthInMbs)!=0);

    currMB->mbAvailD = mb_is_available(currMB->mbAddrD, mb_nr) && ((mb_nr % img->PicWidthInMbs)!=0);
  }

  if (currMB->mbAvailA) currMB->mb_available_left = &(img->mb_data[currMB->mbAddrA]);
  if (currMB->mbAvailB) currMB->mb_available_up   = &(img->mb_data[currMB->mbAddrB]);
}

/*!
************************************************************************
* \brief
*    based on the subpartition of 32x32 block decide the 16x16 block on which the neighborhoold will be decided 
*  \param mb_x: x-coordinate inside 32x32   (in the unit of pel)
*  \param mb_y: y-coordinate inside 32x32   (in the unit of pel)
*  \param reset_block_x
*      input x block position
*  \param reset_block_y
*      input y block position
* \return
*    16x16 block number in the raster scan order
************************************************************************
*/

int ResetAvailabilityOfNeighbors32(int mb_x, int mb_y, int *reset_block_x, int *reset_block_y, int blockshape_x, int blockshape_y)
{
  int reset_mb_nr=img->current_mb_nr;

  if(mb_x==0 && mb_y==0 && blockshape_x==64 && blockshape_y==64) // 64x64 partition
    CheckAvailabilityOfNeighbors32(2);

  if(mb_x==0 && mb_y==0 && blockshape_x==32 && blockshape_y==32) // 32x32 partition
    CheckAvailabilityOfNeighbors32(1);
  
  if(mb_x==0 && mb_y==0 && blockshape_x==64 && blockshape_y==32) // 1st partition 64X32
    CheckAvailabilityOfNeighbors32(2);
  if(mb_x==0 && mb_y==0 && blockshape_x==32 && blockshape_y==16) // 1st partition 32x16
    CheckAvailabilityOfNeighbors32(1);


  if(mb_x==0 && mb_y==32 && blockshape_x==64 && blockshape_y==32)// 2nd partition 64X32
  {
    int saved_mb_nr;
    saved_mb_nr = img->current_mb_nr;
    reset_mb_nr = img->current_mb_nr + img->PicWidthInMbs * 2;
    img->current_mb_nr = reset_mb_nr;
    CheckAvailabilityOfNeighbors32(0);
    img->current_mb_nr = saved_mb_nr;
    *reset_block_y = 0;
  }
  if(mb_x==0 && mb_y==16 && blockshape_x==32 && blockshape_y==16)// 2nd partition 32x16
  {
    int saved_mb_nr;
    saved_mb_nr = img->current_mb_nr;
    reset_mb_nr = img->current_mb_nr + img->PicWidthInMbs;
    img->current_mb_nr = reset_mb_nr;
    CheckAvailabilityOfNeighbors32(0);
    img->current_mb_nr = saved_mb_nr;
    *reset_block_y = 0;
  }

  if(mb_x==0 && mb_y==0 && blockshape_x==32 && blockshape_y==64) // 1st partition 32x64
    CheckAvailabilityOfNeighbors32(1);
  if(mb_x==0 && mb_y==0 && blockshape_x==16 && blockshape_y==32) // 1st partition 16x32
    CheckAvailabilityOfNeighbors32(0);

  if(mb_x==32 && mb_y==0 && blockshape_x==32 && blockshape_y==64)// 2nd partition 32X64
  {
    int saved_mb_nr;
    saved_mb_nr = img->current_mb_nr;
    reset_mb_nr = img->current_mb_nr + 2;
    img->current_mb_nr = reset_mb_nr;
    CheckAvailabilityOfNeighbors32(1);
    img->current_mb_nr = saved_mb_nr;
    *reset_block_x = 0; 
  }
  if(mb_x==16 && mb_y==0 && blockshape_x==16 && blockshape_y==32)// 2nd partition 16x32
  {
    int saved_mb_nr;
    saved_mb_nr = img->current_mb_nr;
    reset_mb_nr = img->current_mb_nr + 1;
    img->current_mb_nr = reset_mb_nr;
    CheckAvailabilityOfNeighbors32(0);
    img->current_mb_nr = saved_mb_nr;
    *reset_block_x = 0; 
  }
  return reset_mb_nr;
}

/*!
************************************************************************
* \brief
*    Set motion vector predictor
************************************************************************
*/
void SetMotionVectorPredictor32 (struct img_par  *img,
                               short           *pmv_x,
                               short           *pmv_y,
                               char            ref_frame,
                               byte            list,
                               char            ***refPic,
                               short           ****tmp_mv,
                               int             block_x,
                               int             block_y,
                               int             blockshape_x,
                               int             blockshape_y,
                               int             mb_ext_level)
{
  int mb_x                 = BLOCK_SIZE*block_x;
  int mb_y                 = BLOCK_SIZE*block_y;
  
  int mv_a, mv_b, mv_c, pred_vec=0;
  int mvPredType, rFrameL, rFrameU, rFrameUR;
  int hv;
  
  
  PixelPos block_a, block_b, block_c, block_d;
  {
    int reset_mb_nr, reset_block_x=block_x, reset_block_y=block_y;

    reset_mb_nr = ResetAvailabilityOfNeighbors32(mb_x, mb_y, &reset_block_x, &reset_block_y, blockshape_x, blockshape_y);

    getLuma4x4Neighbour(reset_mb_nr, reset_block_x, reset_block_y,           -1,  0, &block_a);
    getLuma4x4Neighbour(reset_mb_nr, reset_block_x, reset_block_y,            0, -1, &block_b);
    getLuma4x4Neighbour(reset_mb_nr, reset_block_x, reset_block_y, blockshape_x, -1, &block_c);
    getLuma4x4Neighbour(reset_mb_nr, reset_block_x, reset_block_y,           -1, -1, &block_d);
  }

  if (mb_y > 0)
  {
    if (mb_x < (8<<mb_ext_level))  // first column of 8x8 blocks //MB32X32
    {
      if (mb_y==(8<<mb_ext_level)) //MB32X32
      {
        if (blockshape_x == (16<<mb_ext_level))      block_c.available  = 0; //MB32X32
      }
      else
      {
        error("never reach here!\n", -1);//MB32X32
        if (mb_x+blockshape_x == 8)  block_c.available  = 0; 
      }
    }
    else
    {
      error("never reach here!\n", -1);
      if (mb_x+blockshape_x == 16)   block_c.available  = 0;
    }
  }
  
  if (!block_c.available)
  {
    block_c=block_d;
  }
  
  mvPredType = MVPRED_MEDIAN;
  
  if (!img->MbaffFrameFlag)
  {
    rFrameL    = block_a.available    ? refPic[list][block_a.pos_y][block_a.pos_x] : -1;
    rFrameU    = block_b.available    ? refPic[list][block_b.pos_y][block_b.pos_x] : -1;
    rFrameUR   = block_c.available    ? refPic[list][block_c.pos_y][block_c.pos_x] : -1;
  }
  else
  {
    if (img->mb_data[img->current_mb_nr].mb_field)
    {
      rFrameL    = block_a.available    ? 
        img->mb_data[block_a.mb_addr].mb_field ? 
        refPic[list][block_a.pos_y][block_a.pos_x]:
      refPic[list][block_a.pos_y][block_a.pos_x] * 2: 
      -1;
      rFrameU    = block_b.available    ? 
        img->mb_data[block_b.mb_addr].mb_field ? 
        refPic[list][block_b.pos_y][block_b.pos_x]:
      refPic[list][block_b.pos_y][block_b.pos_x] * 2: 
      -1;
      rFrameUR    = block_c.available    ? 
        img->mb_data[block_c.mb_addr].mb_field ? 
        refPic[list][block_c.pos_y][block_c.pos_x]:
      refPic[list][block_c.pos_y][block_c.pos_x] * 2: 
      -1;
    }
    else
    {
      rFrameL    = block_a.available    ? 
        img->mb_data[block_a.mb_addr].mb_field ? 
        refPic[list][block_a.pos_y][block_a.pos_x] >>1:
      refPic[list][block_a.pos_y][block_a.pos_x] : 
      -1;
      rFrameU    = block_b.available    ? 
        img->mb_data[block_b.mb_addr].mb_field ? 
        refPic[list][block_b.pos_y][block_b.pos_x] >>1:
      refPic[list][block_b.pos_y][block_b.pos_x] : 
      -1;
      rFrameUR    = block_c.available    ? 
        img->mb_data[block_c.mb_addr].mb_field ? 
        refPic[list][block_c.pos_y][block_c.pos_x] >>1:
      refPic[list][block_c.pos_y][block_c.pos_x] : 
      -1;
    }
  }
  
  
  /* Prediction if only one of the neighbors uses the reference frame
  * we are checking
  */
  if(rFrameL == ref_frame && rFrameU != ref_frame && rFrameUR != ref_frame)       mvPredType = MVPRED_L;
  else if(rFrameL != ref_frame && rFrameU == ref_frame && rFrameUR != ref_frame)  mvPredType = MVPRED_U;
  else if(rFrameL != ref_frame && rFrameU != ref_frame && rFrameUR == ref_frame)  mvPredType = MVPRED_UR;
  // Directional predictions 
  if(blockshape_x == (8<<mb_ext_level) && blockshape_y == (16<<mb_ext_level)) //MB32X32
  {
    if(mb_x == 0)
    {
      if(rFrameL == ref_frame)
        mvPredType = MVPRED_L;
    }
    else
    {
      if( rFrameUR == ref_frame)
        mvPredType = MVPRED_UR;
    }
  }
  else if(blockshape_x == (16<<mb_ext_level) && blockshape_y == (8<<mb_ext_level))//MB32X32
  {
    if(mb_y == 0)
    {
      if(rFrameU == ref_frame)
        mvPredType = MVPRED_U;
    }
    else
    {
      if(rFrameL == ref_frame)
        mvPredType = MVPRED_L;
    }
  }
  
  for (hv=0; hv < 2; hv++)
  {
    if (!img->MbaffFrameFlag || hv==0)
    {
      mv_a = block_a.available  ? tmp_mv[list][block_a.pos_y][block_a.pos_x][hv] : 0;
      mv_b = block_b.available  ? tmp_mv[list][block_b.pos_y][block_b.pos_x][hv] : 0;
      mv_c = block_c.available  ? tmp_mv[list][block_c.pos_y][block_c.pos_x][hv] : 0;
    }
    else
    {
      if (img->mb_data[img->current_mb_nr].mb_field)
      {
        mv_a = block_a.available  ? img->mb_data[block_a.mb_addr].mb_field?
          tmp_mv[list][block_a.pos_y][block_a.pos_x][hv]:
        tmp_mv[list][block_a.pos_y][block_a.pos_x][hv] / 2: 
        0;
        mv_b = block_b.available  ? img->mb_data[block_b.mb_addr].mb_field?
          tmp_mv[list][block_b.pos_y][block_b.pos_x][hv]:
        tmp_mv[list][block_b.pos_y][block_b.pos_x][hv] / 2: 
        0;
        mv_c = block_c.available  ? img->mb_data[block_c.mb_addr].mb_field?
          tmp_mv[list][block_c.pos_y][block_c.pos_x][hv]:
        tmp_mv[list][block_c.pos_y][block_c.pos_x][hv] / 2: 
        0;
      }
      else
      {
        mv_a = block_a.available  ? img->mb_data[block_a.mb_addr].mb_field?
          tmp_mv[list][block_a.pos_y][block_a.pos_x][hv] * 2:
        tmp_mv[list][block_a.pos_y][block_a.pos_x][hv]: 
        0;
        mv_b = block_b.available  ? img->mb_data[block_b.mb_addr].mb_field?
          tmp_mv[list][block_b.pos_y][block_b.pos_x][hv] * 2:
        tmp_mv[list][block_b.pos_y][block_b.pos_x][hv]: 
        0;
        mv_c = block_c.available  ? img->mb_data[block_c.mb_addr].mb_field?
          tmp_mv[list][block_c.pos_y][block_c.pos_x][hv] * 2:
        tmp_mv[list][block_c.pos_y][block_c.pos_x][hv]: 
        0;
      }
    }
    
    switch (mvPredType)
    {
    case MVPRED_MEDIAN:
      if(!(block_b.available || block_c.available))
        pred_vec = mv_a;
      else
        pred_vec = mv_a+mv_b+mv_c-min(mv_a,min(mv_b,mv_c))-max(mv_a,max(mv_b,mv_c));
      break;
    case MVPRED_L:
      pred_vec = mv_a;
      break;
    case MVPRED_U:
      pred_vec = mv_b;
      break;
    case MVPRED_UR:
      pred_vec = mv_c;
      break;
    default:
      break;
    }
    
    if (hv==0)  *pmv_x = pred_vec;
    else        *pmv_y = pred_vec;
    
  }
}


/*!
************************************************************************
* \brief
*    init macroblock I and P frames
************************************************************************
*/
void init_macroblock32(struct img_par *img, int mb_ext_level)
{
  int i,j;
  for (i=0;i<(BLOCK_SIZE<<mb_ext_level);i++)
  {                           // reset vectors and pred. modes
    for(j=0;j<(BLOCK_SIZE<<mb_ext_level);j++)
    {
      dec_picture->mv[LIST_0][img->block_y+j][img->block_x+i][0]=0;
      dec_picture->mv[LIST_0][img->block_y+j][img->block_x+i][1]=0;
      dec_picture->mv[LIST_1][img->block_y+j][img->block_x+i][0]=0;
      dec_picture->mv[LIST_1][img->block_y+j][img->block_x+i][1]=0;
      
      img->ipredmode[img->block_x+i][img->block_y+j] = DC_PRED;
    }
  }
  
  for (j=0; j<(BLOCK_SIZE<<mb_ext_level); j++)
    for (i=0; i<(BLOCK_SIZE<<mb_ext_level); i++)
    {
      dec_picture->ref_idx[LIST_0][img->block_y+j][img->block_x+i] = -1;
      dec_picture->ref_idx[LIST_1][img->block_y+j][img->block_x+i] = -1;
      dec_picture->ref_pic_id[LIST_0][img->block_y+j][img->block_x+i] = INT64_MIN;
      dec_picture->ref_pic_id[LIST_1][img->block_y+j][img->block_x+i] = INT64_MIN;
    }
}

/*!
************************************************************************
* \brief
*    Read motion info
************************************************************************
*/
void readMotionInfoFromNAL32 (struct img_par *img, struct inp_par *inp, int mb_ext_level)
{
  int i,j,k;
  int step_h,step_v;
  int curr_mvd;
  Macroblock *currMB  = &img->mb_data[img->current_mb_nr];
  SyntaxElement currSE;
  Slice *currSlice    = img->currentSlice;
  DataPartition *dP;
  int *partMap        = assignSE2partition[currSlice->dp_mode];
  int bframe          = (img->type==B_SLICE);
  int partmode        = (IS_P8x8(currMB)?4:currMB->mb_type);
  int step_h0         = BLOCK_STEP [partmode][0]<<mb_ext_level;
  int step_v0         = BLOCK_STEP [partmode][1]<<mb_ext_level;
  int mv_mode, i0, j0;
  char refframe;
  short pmv[2];
  int j4, i4, ii,jj;
  
  int motion_vector;
  int mv_scale = 0;
  
  int flag_mode;
  
  int list_offset = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;
  
  byte  **    moving_block;
  short ****  co_located_mv;
  char  ***   co_located_ref_idx;
  int64 ***   co_located_ref_id;
  
  if ((img->MbaffFrameFlag)&&(currMB->mb_field))
  {
    if(img->current_mb_nr%2)
    {
      moving_block = Co_located->bottom_moving_block;
      co_located_mv = Co_located->bottom_mv;
      co_located_ref_idx = Co_located->bottom_ref_idx;
      co_located_ref_id = Co_located->bottom_ref_pic_id;
    }
    else
    {
      moving_block = Co_located->top_moving_block;
      co_located_mv = Co_located->top_mv;
      co_located_ref_idx = Co_located->top_ref_idx;
      co_located_ref_id = Co_located->top_ref_pic_id;
    }
  }
  else
  {
    moving_block = Co_located->moving_block;
    co_located_mv = Co_located->mv;
    co_located_ref_idx = Co_located->ref_idx;
    co_located_ref_id = Co_located->ref_pic_id;
  }
  
  if (bframe && IS_P8x8 (currMB))
  {
    if (img->direct_spatial_mv_pred_flag)
    {
      int imgblock_y= ((img->MbaffFrameFlag)&&(currMB->mb_field))? (img->current_mb_nr%2) ? (img->block_y-4)/2:img->block_y/2: img->block_y;
      int fw_rFrameL, fw_rFrameU, fw_rFrameUL, fw_rFrameUR;
      int bw_rFrameL, bw_rFrameU, bw_rFrameUL, bw_rFrameUR;
      
      PixelPos mb_left, mb_up, mb_upleft, mb_upright;
      
      char  fw_rFrame,bw_rFrame;
      short pmvfw[2]={0,0},
        pmvbw[2]={0,0};
      
      
      getLuma4x4Neighbour(img->current_mb_nr, 0, 0, -1,  0, &mb_left);
      getLuma4x4Neighbour(img->current_mb_nr, 0, 0,  0, -1, &mb_up);
      getLuma4x4Neighbour(img->current_mb_nr, 0, 0, 16, -1, &mb_upright);
      getLuma4x4Neighbour(img->current_mb_nr, 0, 0, -1, -1, &mb_upleft);
      
      if (!img->MbaffFrameFlag)
      {
        fw_rFrameL = mb_left.available ? dec_picture->ref_idx[LIST_0][mb_left.pos_y][mb_left.pos_x] : -1;
        fw_rFrameU = mb_up.available ? dec_picture->ref_idx[LIST_0][mb_up.pos_y][mb_up.pos_x] : -1;
        fw_rFrameUL = mb_upleft.available ? dec_picture->ref_idx[LIST_0][mb_upleft.pos_y][mb_upleft.pos_x] : -1;
        fw_rFrameUR = mb_upright.available ? dec_picture->ref_idx[LIST_0][mb_upright.pos_y][mb_upright.pos_x] : fw_rFrameUL;      
        
        bw_rFrameL = mb_left.available ? dec_picture->ref_idx[LIST_1][mb_left.pos_y][mb_left.pos_x] : -1;
        bw_rFrameU = mb_up.available ? dec_picture->ref_idx[LIST_1][mb_up.pos_y][mb_up.pos_x] : -1;
        bw_rFrameUL = mb_upleft.available ? dec_picture->ref_idx[LIST_1][mb_upleft.pos_y][mb_upleft.pos_x] : -1;
        bw_rFrameUR = mb_upright.available ? dec_picture->ref_idx[LIST_1][mb_upright.pos_y][mb_upright.pos_x] : bw_rFrameUL;      
      }
      else
      {
        if (img->mb_data[img->current_mb_nr].mb_field)
        {
          fw_rFrameL = mb_left.available ? 
            img->mb_data[mb_left.mb_addr].mb_field  || dec_picture->ref_idx[LIST_0][mb_left.pos_y][mb_left.pos_x] < 0? 
            dec_picture->ref_idx[LIST_0][mb_left.pos_y][mb_left.pos_x] : 
          dec_picture->ref_idx[LIST_0][mb_left.pos_y][mb_left.pos_x] * 2: -1;
          fw_rFrameU = mb_up.available ? 
            img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_up.pos_y][mb_up.pos_x] < 0? 
            dec_picture->ref_idx[LIST_0][mb_up.pos_y][mb_up.pos_x] : 
          dec_picture->ref_idx[LIST_0][mb_up.pos_y][mb_up.pos_x] * 2: -1;
          
          fw_rFrameUL = mb_upleft.available ? 
            img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upleft.pos_y][mb_upleft.pos_x] < 0?
            dec_picture->ref_idx[LIST_0][mb_upleft.pos_y][mb_upleft.pos_x] : 
          dec_picture->ref_idx[LIST_0][mb_upleft.pos_y][mb_upleft.pos_x] *2: -1;      
          
          fw_rFrameUR = mb_upright.available ? 
            img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upright.pos_y][mb_upright.pos_x] < 0 ?
            dec_picture->ref_idx[LIST_0][mb_upright.pos_y][mb_upright.pos_x] : 
          dec_picture->ref_idx[LIST_0][mb_upright.pos_y][mb_upright.pos_x] * 2: fw_rFrameUL;      
          
          bw_rFrameL = mb_left.available ? 
            img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_left.pos_y][mb_left.pos_x]  < 0 ?
            dec_picture->ref_idx[LIST_1][mb_left.pos_y][mb_left.pos_x] : 
          dec_picture->ref_idx[LIST_1][mb_left.pos_y][mb_left.pos_x] * 2: -1;
          
          bw_rFrameU = mb_up.available ? 
            img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_up.pos_y][mb_up.pos_x]  < 0 ?
            dec_picture->ref_idx[LIST_1][mb_up.pos_y][mb_up.pos_x] : 
          dec_picture->ref_idx[LIST_1][mb_up.pos_y][mb_up.pos_x] * 2: -1;
          
          bw_rFrameUL = mb_upleft.available ? 
            img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upleft.pos_y][mb_upleft.pos_x]  < 0 ?
            dec_picture->ref_idx[LIST_1][mb_upleft.pos_y][mb_upleft.pos_x] : 
          dec_picture->ref_idx[LIST_1][mb_upleft.pos_y][mb_upleft.pos_x] *2: -1;      
          
          bw_rFrameUR = mb_upright.available ? 
            img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upright.pos_y][mb_upright.pos_x]  < 0 ?
            dec_picture->ref_idx[LIST_1][mb_upright.pos_y][mb_upright.pos_x] : 
          dec_picture->ref_idx[LIST_1][mb_upright.pos_y][mb_upright.pos_x] * 2: bw_rFrameUL;      
          
        }
        else
        {
          fw_rFrameL = mb_left.available ? 
            img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_left.pos_y][mb_left.pos_x]  < 0 ?
            dec_picture->ref_idx[LIST_0][mb_left.pos_y][mb_left.pos_x] >> 1 : 
          dec_picture->ref_idx[LIST_0][mb_left.pos_y][mb_left.pos_x]: -1;
          
          fw_rFrameU = mb_up.available ? 
            img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_up.pos_y][mb_up.pos_x]  < 0 ?
            dec_picture->ref_idx[LIST_0][mb_up.pos_y][mb_up.pos_x] >> 1 :  
          dec_picture->ref_idx[LIST_0][mb_up.pos_y][mb_up.pos_x] : -1;
          
          fw_rFrameUL = mb_upleft.available ? 
            img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upleft.pos_y][mb_upleft.pos_x] < 0 ?
            dec_picture->ref_idx[LIST_0][mb_upleft.pos_y][mb_upleft.pos_x]>> 1 : 
          dec_picture->ref_idx[LIST_0][mb_upleft.pos_y][mb_upleft.pos_x] : -1;      
          
          fw_rFrameUR = mb_upright.available ? 
            img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upright.pos_y][mb_upright.pos_x]  < 0 ?
            dec_picture->ref_idx[LIST_0][mb_upright.pos_y][mb_upright.pos_x] >> 1 :  
          dec_picture->ref_idx[LIST_0][mb_upright.pos_y][mb_upright.pos_x] : fw_rFrameUL;      
          
          bw_rFrameL = mb_left.available ? 
            img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_left.pos_y][mb_left.pos_x] < 0 ?
            dec_picture->ref_idx[LIST_1][mb_left.pos_y][mb_left.pos_x] >> 1 :  
          dec_picture->ref_idx[LIST_1][mb_left.pos_y][mb_left.pos_x] : -1;
          bw_rFrameU = mb_up.available ? 
            img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_up.pos_y][mb_up.pos_x] < 0 ?
            dec_picture->ref_idx[LIST_1][mb_up.pos_y][mb_up.pos_x] >> 1 : 
          dec_picture->ref_idx[LIST_1][mb_up.pos_y][mb_up.pos_x] : -1;
          
          bw_rFrameUL = mb_upleft.available ? 
            img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upleft.pos_y][mb_upleft.pos_x]  < 0 ?
            dec_picture->ref_idx[LIST_1][mb_upleft.pos_y][mb_upleft.pos_x] >> 1 : 
          dec_picture->ref_idx[LIST_1][mb_upleft.pos_y][mb_upleft.pos_x] : -1;      
          
          bw_rFrameUR = mb_upright.available ? 
            img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upright.pos_y][mb_upright.pos_x] < 0 ?
            dec_picture->ref_idx[LIST_1][mb_upright.pos_y][mb_upright.pos_x] >> 1: 
          dec_picture->ref_idx[LIST_1][mb_upright.pos_y][mb_upright.pos_x] : bw_rFrameUL;      
        }
      }
      
      fw_rFrame = (fw_rFrameL >= 0 && fw_rFrameU >= 0) ? min(fw_rFrameL,fw_rFrameU): max(fw_rFrameL,fw_rFrameU);
      fw_rFrame = (fw_rFrame >= 0 && fw_rFrameUR >= 0) ? min(fw_rFrame,fw_rFrameUR): max(fw_rFrame,fw_rFrameUR);
      
      bw_rFrame = (bw_rFrameL >= 0 && bw_rFrameU >= 0) ? min(bw_rFrameL,bw_rFrameU): max(bw_rFrameL,bw_rFrameU);
      bw_rFrame = (bw_rFrame >= 0 && bw_rFrameUR >= 0) ? min(bw_rFrame,bw_rFrameUR): max(bw_rFrame,bw_rFrameUR);
      
      
      if (fw_rFrame >=0)
        SetMotionVectorPredictor (img, pmvfw, pmvfw+1, fw_rFrame, LIST_0, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);
      
      if (bw_rFrame >=0)
        SetMotionVectorPredictor (img, pmvbw, pmvbw+1, bw_rFrame, LIST_1, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);
      
      
      for (i=0;i<4;i++)
      {
        if (currMB->b8mode[i] == 0)
          for(j=2*(i/2);j<2*(i/2)+2;j++)
            for(k=2*(i%2);k<2*(i%2)+2;k++)
            {
              int j6 = imgblock_y+j;
              j4 = img->block_y+j;
              i4 = img->block_x+k;
              
              
              if (fw_rFrame >= 0)
              {
                
                if  (!fw_rFrame  && ((!moving_block[j6][i4]) && (!listX[1+list_offset][0]->is_long_term)))
                {                    
                  dec_picture->mv  [LIST_0][j4][i4][0] = 0;
                  dec_picture->mv  [LIST_0][j4][i4][1] = 0;
                  dec_picture->ref_idx[LIST_0][j4][i4] = 0;                    
                }
                else
                {
                  
                  dec_picture->mv  [LIST_0][j4][i4][0] = pmvfw[0];
                  dec_picture->mv  [LIST_0][j4][i4][1] = pmvfw[1];
                  dec_picture->ref_idx[LIST_0][j4][i4] = fw_rFrame;
                }
              }
              else
              {
                dec_picture->mv  [LIST_0][j4][i4][0] = 0;
                dec_picture->mv  [LIST_0][j4][i4][1] = 0;
                dec_picture->ref_idx[LIST_0][j4][i4] = -1;
              }
              if (bw_rFrame >= 0)
              {
                if  (bw_rFrame==0 && ((!moving_block[j6][i4])&& (!listX[1+list_offset][0]->is_long_term)))
                {
                  dec_picture->mv  [LIST_1][j4][i4][0] = 0;
                  dec_picture->mv  [LIST_1][j4][i4][1] = 0;
                  dec_picture->ref_idx[LIST_1][j4][i4] = 0;
                }
                else
                {
                  dec_picture->mv  [LIST_1][j4][i4][0] = pmvbw[0];
                  dec_picture->mv  [LIST_1][j4][i4][1] = pmvbw[1];
                  dec_picture->ref_idx[LIST_1][j4][i4] = bw_rFrame;
                }
              }
              else
              {
                dec_picture->mv  [LIST_1][j4][i4][0] = 0;
                dec_picture->mv  [LIST_1][j4][i4][1] = 0;
                dec_picture->ref_idx[LIST_1][j4][i4] = -1;                               
              }
              
              if (fw_rFrame <0 && bw_rFrame <0)
              {
                dec_picture->ref_idx[LIST_0][j4][i4] = 0;
                dec_picture->ref_idx[LIST_1][j4][i4] = 0;                  
              }
            }
      }
    }
    else
    {
      for (i=0;i<4;i++)
      {
        if (currMB->b8mode[i] == 0)
        {
          for(j=2*(i/2);j<2*(i/2)+2;j++)
          {
            for(k=2*(i%2);k<2*(i%2)+2;k++)
            {
              
              int list_offset = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;
              int imgblock_y= ((img->MbaffFrameFlag)&&(currMB->mb_field))? (img->current_mb_nr%2) ? (img->block_y-4)/2 : img->block_y/2 : img->block_y;
              int refList = co_located_ref_idx[LIST_0 ][imgblock_y+j][img->block_x+k]== -1 ? LIST_1 : LIST_0;
              int ref_idx = co_located_ref_idx[refList][imgblock_y + j][img->block_x + k];
              int mapped_idx=-1, iref;                             
              
              if (ref_idx == -1)
              {
                dec_picture->ref_idx [LIST_0][img->block_y + j][img->block_x + k] = 0;
                dec_picture->ref_idx [LIST_1][img->block_y + j][img->block_x + k] = 0;                
              }
              else
              {
                for (iref=0;iref<min(img->num_ref_idx_l0_active,listXsize[LIST_0 + list_offset]);iref++)
                {
#if 1
                  int curr_mb_field = ((img->MbaffFrameFlag)&&(currMB->mb_field));
                  
                  if(img->structure==0 && curr_mb_field==0)
                  {
                    // If the current MB is a frame MB and the colocated is from a field picture, 
                    // then the co_located_ref_id may have been generated from the wrong value of 
                    // frame_poc if it references it's complementary field, so test both POC values
                    if(listX[0][iref]->top_poc*2 == co_located_ref_id[refList][imgblock_y + j][img->block_x + k] 
                      || listX[0][iref]->bottom_poc*2 == co_located_ref_id[refList][imgblock_y + j][img->block_x + k])
                    {
                      mapped_idx=iref;
                      break;
                    }
                    else //! invalid index. Default to zero even though this case should not happen
                      mapped_idx=INVALIDINDEX;
                    continue;
                  }    
#endif                                        
                  if (dec_picture->ref_pic_num[img->current_slice_nr][LIST_0 + list_offset][iref]==co_located_ref_id[refList][imgblock_y + j][img->block_x + k])
                  {
                    mapped_idx=iref;
                    break;
                  }
                  else //! invalid index. Default to zero even though this case should not happen
                    mapped_idx=INVALIDINDEX;
                }
                if (INVALIDINDEX == mapped_idx)
                {
                  error("temporal direct error\ncolocated block has ref that is unavailable",-1111);
                }
                dec_picture->ref_idx [LIST_0][img->block_y + j][img->block_x + k] = mapped_idx;
                dec_picture->ref_idx [LIST_1][img->block_y + j][img->block_x + k] = 0;                
              }
            }
          }
        }
      }
    }
  } 
  
  //  If multiple ref. frames, read reference frame for the MB *********************************
  if(img->num_ref_idx_l0_active>1) 
  {
    flag_mode = ( img->num_ref_idx_l0_active == 2 ? 1 : 0);
    
    currSE.type = SE_REFFRAME;
    dP = &(currSlice->partArr[partMap[SE_REFFRAME]]);
    
    if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)   currSE.mapping = linfo_ue;
    else                              
    {
      currSE.mb_ext_level = mb_ext_level;
      currSE.reading = readRefFrame_CABAC32;
    }

    for (j0=0; j0<(4<<mb_ext_level); j0+=step_v0)
    {
      for (i0=0; i0<(4<<mb_ext_level); i0+=step_h0)
      {
        k=2*(j0>>(mb_ext_level+1))+(i0>>(mb_ext_level+1));
        if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)
        {
          TRACE_STRING("ref_idx_l0");
          
          img->subblock_x = i0;
          img->subblock_y = j0;
          
          if (!IS_P8x8 (currMB) || bframe || (!bframe && !img->allrefzero))
          {
            currSE.context = BType2CtxRef (currMB->b8mode[k]);
            if( (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag) && flag_mode )
            {
              currSE.len = 1;
              readSyntaxElement_FLC(&currSE, dP->bitstream);
              currSE.value1 = 1 - currSE.value1;
            }
            else
            {
              currSE.value2 = LIST_0;
              dP->readSyntaxElement (&currSE,img,inp,dP);
            }
            refframe = currSE.value1;            
          }
          else
          {
            refframe = 0;
          }
          
          /*
          if (bframe && refframe>img->buf_cycle)    // img->buf_cycle should be correct for field MBs now
          {
          set_ec_flag(SE_REFFRAME);
          refframe = 1;
          }
          */
          
          for (j=j0; j<j0+step_v0;j++)
            for (i=i0; i<i0+step_h0;i++)
            {
              dec_picture->ref_idx[LIST_0][img->block_y + j][img->block_x + i] = refframe;
            }
            
        }
      }
    }
  }
  else
  {
    for (j0=0; j0<(4<<mb_ext_level); j0+=step_v0)
    {
      for (i0=0; i0<(4<<mb_ext_level); i0+=step_h0)
      {
        k=2*(j0>>(mb_ext_level+1))+(i0>>(mb_ext_level+1));
        if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)
        {
          for (j=j0; j<j0+step_v0;j++)
            for (i=i0; i<i0+step_h0;i++)
            {
              dec_picture->ref_idx[LIST_0][img->block_y + j][img->block_x + i] = 0;
            }
        }
      }
    }
  }
  
  //  If backward multiple ref. frames, read backward reference frame for the MB *********************************
  if(img->num_ref_idx_l1_active>1)
  {
    flag_mode = ( img->num_ref_idx_l1_active == 2 ? 1 : 0);
    
    currSE.type = SE_REFFRAME;
    dP = &(currSlice->partArr[partMap[SE_REFFRAME]]);
    if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)   currSE.mapping = linfo_ue;
    else
    {
      currSE.mb_ext_level = mb_ext_level;
      currSE.reading = readRefFrame_CABAC32;
    }

    for (j0=0; j0<(4<<mb_ext_level); j0+=step_v0)
    {
      for (i0=0; i0<(4<<mb_ext_level); i0+=step_h0)
      {
        k=2*(j0>>(mb_ext_level+1))+(i0>>(mb_ext_level+1));
        if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)
        {
          TRACE_STRING("ref_idx_l1");
          
          img->subblock_x = i0;
          img->subblock_y = j0;
          
          currSE.context = BType2CtxRef (currMB->b8mode[k]);
          if( (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag) && flag_mode )
          {
            currSE.len = 1;
            readSyntaxElement_FLC(&currSE, dP->bitstream);
            currSE.value1 = 1-currSE.value1;
          }
          else
          {
            currSE.value2 = LIST_1;
            dP->readSyntaxElement (&currSE,img,inp,dP);
          }
          refframe = currSE.value1;

          for (j=j0; j<j0+step_v0;j++)
          {
            for (i=i0; i<i0+step_h0;i++)
            {
              dec_picture->ref_idx[LIST_1][img->block_y + j][img->block_x + i] = refframe;
            }
          }
        }
      }
    }
  }
  else
  {
    for (j0=0; j0<(4<<mb_ext_level); j0+=step_v0)
    {
      for (i0=0; i0<(4<<mb_ext_level); i0+=step_h0)
      {
        k=2*(j0>>(mb_ext_level+1))+(i0>>(mb_ext_level+1));
        if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)
        {
          for (j=j0; j<j0+step_v0;j++)
            for (i=i0; i<i0+step_h0;i++)
            {
              dec_picture->ref_idx[LIST_1][img->block_y + j][img->block_x + i] = 0;
            }
        }
      }
    }
  }
  
  //=====  READ FORWARD MOTION VECTORS =====
  currSE.type = SE_MVD;
  dP = &(currSlice->partArr[partMap[SE_MVD]]);
  
  if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag) currSE.mapping = linfo_se;
  else                                               
  {
    currSE.mb_ext_level = mb_ext_level;
    currSE.reading = readMVD_CABAC32;
  }

  for (j0=0; j0<(4<<mb_ext_level); j0+=step_v0)
    for (i0=0; i0<(4<<mb_ext_level); i0+=step_h0)
    {
      k=2*(j0>>(mb_ext_level+1))+(i0>>(mb_ext_level+1));
      
      if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && (currMB->b8mode[k] !=0))//has forward vector
      {
        mv_mode  = currMB->b8mode[k];
        step_h   = BLOCK_STEP [mv_mode][0]<<mb_ext_level;
        step_v   = BLOCK_STEP [mv_mode][1]<<mb_ext_level;
        refframe = dec_picture->ref_idx[LIST_0][img->block_y+j0][img->block_x+i0];
        
        for (j=j0; j<j0+step_v0; j+=step_v)
        {
          for (i=i0; i<i0+step_h0; i+=step_h)
          {
            j4 = img->block_y+j;
            i4 = img->block_x+i;
            
            // first make mv-prediction
#ifdef MB32X32_MVC
            if (mv_comp.mv_competition > 0)
            {
              int mvd[2];

              for (k=0; k < 2; k++) 
              {
                TRACE_STRING("mvd_l0");
                
                img->subblock_x = i; // position used for context determination
                img->subblock_y = j; // position used for context determination
                currSE.value2 = k<<1; // identifies the component; only used for context determination
                dP->readSyntaxElement(&currSE,img,inp,dP);
                curr_mvd = currSE.value1; 
                
                mvd[k]=curr_mvd;
              }
              
              SetMotionVectorPredictor_Competition32 (img,inp, pmv, pmv+1, refframe, LIST_0, dec_picture->ref_idx, dec_picture->mv, i, j, 4*step_h, 4*step_v, mb_ext_level);
              
              for (k=0; k < 2; k++) 
              {
                for(ii=0;ii<step_h;ii++)
                {
                  for(jj=0;jj<step_v;jj++)
                  {
                    dec_picture->mv        [LIST_0][j4+jj][i4+ii][k] = mvd[k] + pmv[k];
                                //the mvd for the whole 32x32 is saved in mb32_mvd and will be copied to each 16x16 block's mvd in write_macroblock_cluster
                                 mb32_mvd  [LIST_0][j +jj][i +ii][k] = mvd[k];
                  }
                }

                //only save mvd for the upper left 16x16 block for the context of the 2nd partition of 32x16 and 16x32
                if(j==0 && i==0)
                {
                  for(ii=0;ii<4;ii++)
                  {
                    for(jj=0;jj<4;jj++)
                    {
                      currMB->mvd      [LIST_0][j+jj] [i+ii] [k] = mvd[k];
                    }
                  }
                }
              }
            }
            
            else
#endif
            {
              
              SetMotionVectorPredictor32 (img, pmv, pmv+1, refframe, LIST_0, dec_picture->ref_idx, dec_picture->mv, i, j, 4*step_h, 4*step_v, mb_ext_level);
              
              for (k=0; k < 2; k++) 
              {
                TRACE_STRING("mvd_l0");
                
                img->subblock_x = i; // position used for context determination
                img->subblock_y = j; // position used for context determination
                currSE.value2 = k<<1; // identifies the component; only used for context determination
                dP->readSyntaxElement(&currSE,img,inp,dP);
                curr_mvd = currSE.value1; 
                
                motion_vector=curr_mvd+pmv[k];           /* find motion vector */
                

                for(ii=0;ii<step_h;ii++)
                {
                  for(jj=0;jj<step_v;jj++)
                  {
                    dec_picture->mv  [LIST_0][j4+jj][i4+ii][k] = motion_vector;
                  }
                }

                //only save mvd for the upper left 16x16 block for the context of the 2nd partition of 32x16 and 16x32
                if(j==0 && i==0)
                {
                  for(ii=0;ii<4;ii++)
                  {
                    for(jj=0;jj<4;jj++)
                    {
                      currMB->mvd      [LIST_0][j+jj] [i+ii] [k] = curr_mvd;
                    }
                  }
                }
                //the mvd for the whole 32x32 is saved in mb32_mvd and will be copied to each 16x16 block's mvd in write_macroblock_cluster
                for(ii=0;ii<step_h;ii++)
                {
                  for(jj=0;jj<step_v;jj++)
                  {
                    mb32_mvd      [LIST_0][j+jj] [i+ii] [k] = curr_mvd;
                  }
                }


              }
            }
          }
        }
      }
      else if (currMB->b8mode[k=2*(j0>>(mb_ext_level+1))+(i0>>(mb_ext_level+1))]==0)      
      {  
        if (!img->direct_spatial_mv_pred_flag)
        {
          int list_offset = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;
          int imgblock_y= ((img->MbaffFrameFlag)&&(currMB->mb_field))? (img->current_mb_nr%2) ? (img->block_y-4)/2:img->block_y/2 : img->block_y;
          
          int refList = (co_located_ref_idx[LIST_0][imgblock_y+j0][img->block_x+i0]== -1 ? LIST_1 : LIST_0);
          int ref_idx =  co_located_ref_idx[refList][imgblock_y+j0][img->block_x+i0];          
          
          if (ref_idx==-1)
          {
            for (j=j0; j<j0+step_v0; j++)
              for (i=i0; i<i0+step_h0; i++)
              {            
                dec_picture->ref_idx [LIST_1][img->block_y+j][img->block_x+i]=0;
                dec_picture->ref_idx [LIST_0][img->block_y+j][img->block_x+i]=0; 
                j4 = img->block_y+j;
                i4 = img->block_x+i;            
                for (ii=0; ii < 2; ii++) 
                {                                    
                  dec_picture->mv [LIST_0][j4][i4][ii]=0;
                  dec_picture->mv [LIST_1][j4][i4][ii]=0;                  
                }
              }
          }
          else 
          {        
            int mapped_idx=-1, iref;                             
            int j6;
            
            for (iref=0;iref<min(img->num_ref_idx_l0_active,listXsize[LIST_0 + list_offset]);iref++)
            {
              
#if 1
              int curr_mb_field = ((img->MbaffFrameFlag)&&(currMB->mb_field));
              
              if(img->structure==0 && curr_mb_field==0)
              {
                // If the current MB is a frame MB and the colocated is from a field picture, 
                // then the co_located_ref_id may have been generated from the wrong value of 
                // frame_poc if it references it's complementary field, so test both POC values
                if(listX[0][iref]->top_poc*2 == co_located_ref_id[refList][imgblock_y + j0][img->block_x + i0] 
                  || listX[0][iref]->bottom_poc*2 == co_located_ref_id[refList][imgblock_y + j0][img->block_x + i0])
                {
                  mapped_idx=iref;
                  break;
                }
                else //! invalid index. Default to zero even though this case should not happen
                  mapped_idx=INVALIDINDEX;
                continue;
              }    
#endif                                        
              if (dec_picture->ref_pic_num[img->current_slice_nr][LIST_0 + list_offset][iref]==co_located_ref_id[refList][imgblock_y+j0][img->block_x+i0])
              {
                mapped_idx=iref;
                break;
              }
              else //! invalid index. Default to zero even though this case should not happen
                mapped_idx=INVALIDINDEX;
            }
            
            if (INVALIDINDEX == mapped_idx)
            {
              error("temporal direct error\ncolocated block has ref that is unavailable",-1111);
            }
            
            
            for (j=j0; j<j0+step_v0; j++)
              for (i=i0; i<i0+step_h0; i++)
              {
                {
                  mv_scale = img->mvscale[LIST_0 + list_offset][mapped_idx];
                  
                  dec_picture->ref_idx [LIST_0][img->block_y+j][img->block_x+i] = mapped_idx;
                  dec_picture->ref_idx [LIST_1][img->block_y+j][img->block_x+i] = 0;
                  
                  j4 = img->block_y+j;
                  j6 = imgblock_y+j;
                  i4 = img->block_x+i;
                  
                  for (ii=0; ii < 2; ii++) 
                  {              
                    //if (iTRp==0)
                    if (mv_scale == 9999 || listX[LIST_0+list_offset][mapped_idx]->is_long_term)
                      //                    if (mv_scale==9999 || Co_located->is_long_term)
                    {                      
                      dec_picture->mv  [LIST_0][j4][i4][ii]=co_located_mv[refList][j6][i4][ii];
                      dec_picture->mv  [LIST_1][j4][i4][ii]=0;
                    }
                    else
                    {
                      dec_picture->mv  [LIST_0][j4][i4][ii]=(mv_scale * co_located_mv[refList][j6][i4][ii] + 128 ) >> 8;
                      dec_picture->mv  [LIST_1][j4][i4][ii]=dec_picture->mv[LIST_0][j4][i4][ii] - co_located_mv[refList][j6][i4][ii];
                    }
                  }
                } 
              }
          }  
        } 
      }
  }
  
  //=====  READ BACKWARD MOTION VECTORS =====
  currSE.type = SE_MVD;
  dP          = &(currSlice->partArr[partMap[SE_MVD]]);
  
  if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag) currSE.mapping = linfo_se;
  else                                       
  {
    currSE.mb_ext_level = mb_ext_level;
    currSE.reading = readMVD_CABAC32;
  }

  for (j0=0; j0<(4<<mb_ext_level); j0+=step_v0)
  {
    for (i0=0; i0<(4<<mb_ext_level); i0+=step_h0)
    {
      k=2*(j0>>(mb_ext_level+1))+(i0>>(mb_ext_level+1));
      if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && (currMB->b8mode[k]!=0))//has backward vector
      {
        mv_mode  = currMB->b8mode[k];
        step_h   = BLOCK_STEP [mv_mode][0]<<mb_ext_level;
        step_v   = BLOCK_STEP [mv_mode][1]<<mb_ext_level;
        refframe = dec_picture->ref_idx[LIST_1][img->block_y+j0][img->block_x+i0];
        
        for (j=j0; j<j0+step_v0; j+=step_v)
        {
          for (i=i0; i<i0+step_h0; i+=step_h)
          {
            j4 = img->block_y+j;
            i4 = img->block_x+i;
            
            // first make mv-prediction
            
#ifdef MB32X32_MVC
            if (mv_comp.mv_competition > 0)
            {
              int mvd[2];
              for (k=0; k < 2; k++) 
              {
                TRACE_STRING("mvd_l1");
                
                img->subblock_x = i; // position used for context determination
                img->subblock_y = j; // position used for context determination
                currSE.value2   = (k<<1) +1; // identifies the component; only used for context determination
                dP->readSyntaxElement(&currSE,img,inp,dP);
                curr_mvd = currSE.value1; 
                
                mvd[k] = curr_mvd;
              }
 
              SetMotionVectorPredictor_Competition32 (img,inp, pmv, pmv+1, refframe, LIST_1, dec_picture->ref_idx, dec_picture->mv, i, j, 4*step_h, 4*step_v, mb_ext_level);
              
              for (k=0; k < 2; k++) 
              {
                for(ii=0;ii<step_h;ii++)
                {
                  for(jj=0;jj<step_v;jj++)
                  {
                    dec_picture->mv        [LIST_1][j4+jj][i4+ii][k] = mvd[k] + pmv[k];
                                //the mvd for the whole 32x32 is saved in mb32_mvd and will be copied to each 16x16 block's mvd in write_macroblock_cluster
                                 mb32_mvd  [LIST_1][j +jj][i +ii][k] = mvd[k];
                  }
                }

                //only save mvd for the upper left 16x16 block for the context of the 2nd partition of 32x16 and 16x32
                if(j==0 && i==0)
                {
                  for(ii=0;ii<4;ii++)
                  {
                    for(jj=0;jj<4;jj++)
                    {
                      currMB->mvd      [LIST_1][j+jj] [i+ii] [k] = mvd[k];
                    }
                  }
                }
              }
            }
            else
              
#endif
            {
              SetMotionVectorPredictor32 (img, pmv, pmv+1, refframe, LIST_1, dec_picture->ref_idx, dec_picture->mv, i, j, 4*step_h, 4*step_v, mb_ext_level);

              for (k=0; k < 2; k++) 
              {
                TRACE_STRING("mvd_l1");
                
                img->subblock_x = i; // position used for context determination
                img->subblock_y = j; // position used for context determination
                currSE.value2   = (k<<1) +1; // identifies the component; only used for context determination
                dP->readSyntaxElement(&currSE,img,inp,dP);
                curr_mvd = currSE.value1; 
                
                motion_vector=curr_mvd+pmv[k];           /* find motion vector */
                
                for(ii=0;ii<step_h;ii++)
                {
                  for(jj=0;jj<step_v;jj++)
                  {
                    dec_picture->mv  [LIST_1][j4+jj][i4+ii][k] = motion_vector;
                  }
                }

                //only save mvd for the upper left 16x16 block for the context of the 2nd partition of 32x16 and 16x32
                if(j==0 && i==0)
                {
                  for(ii=0;ii<4;ii++)
                  {
                    for(jj=0;jj<4;jj++)
                    {
                      currMB->mvd      [LIST_1][j+jj] [i+ii] [k] = curr_mvd;
                    }
                  }
                }
                //the mvd for the whole 32x32 is saved in mb32_mvd and will be copied to each 16x16 block's mvd in write_macroblock_cluster
                for(ii=0;ii<step_h;ii++)
                {
                  for(jj=0;jj<step_v;jj++)
                  {
                    mb32_mvd      [LIST_1][j+jj] [i+ii] [k] = curr_mvd;
                  }
                }


              }
            }
          }
        }
      }
    }
  }
  // record reference picture Ids for deblocking decisions
  for(i4=img->block_x;i4<(img->block_x+(4<<mb_ext_level));i4++)
    for(j4=img->block_y;j4<(img->block_y+(4<<mb_ext_level));j4++)
    {
      if(dec_picture->ref_idx[LIST_0][j4][i4]>=0)
        dec_picture->ref_pic_id[LIST_0][j4][i4] = dec_picture->ref_pic_num[img->current_slice_nr][LIST_0 + list_offset][(short)dec_picture->ref_idx[LIST_0][j4][i4]];
      else
        dec_picture->ref_pic_id[LIST_0][j4][i4] = INT64_MIN;
      if(dec_picture->ref_idx[LIST_1][j4][i4]>=0)
        dec_picture->ref_pic_id[LIST_1][j4][i4] = dec_picture->ref_pic_num[img->current_slice_nr][LIST_1 + list_offset][(short)dec_picture->ref_idx[LIST_1][j4][i4]];  
      else
        dec_picture->ref_pic_id[LIST_1][j4][i4] = INT64_MIN;  
    }
}

/*!
************************************************************************
* \brief
*    Get the syntax elements from the NAL
* \param (mb32,mbi)
*        (1,0) - 32x32, 32x16, 16x32 mode info, mv info, coeff of upper left 16x16, dequant
*        (1,1~)-                                         coeff of other      16x16
*        (0,0) - 16x16 mode, upper left 16x16 partition mode, mv, coeff, dequant
*        (0,1~)-                        16x16 partition mode, mv, coeff
************************************************************************
*/
int read_one_macroblock32(struct img_par *img,struct inp_par *inp, int mbi)
{
  int i;
  int mb32=0;

  SyntaxElement currSE;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  
  Slice *currSlice = img->currentSlice;
  DataPartition *dP;
  int *partMap = assignSE2partition[currSlice->dp_mode];
  Macroblock *topMB = NULL;
  int  prevMbSkipped = 0;
  int  img_block_y;
  int  check_bottom, read_bottom, read_top;
  
#ifdef ADAPTIVE_FD_SD_CODING
  int           j;
  for (i=0;i<16;i++)
  {
    for (j=0;j<16;j++)
    {
      currMB->quantizer_indices[j][i]=0;
    }
  }
  for (i=0;i<2;i++)
  {
    for (j=0;j<2;j++)
    {
      currMB->SD_or_FD     [j][i]=0;
    }
  }
  currMB->SD_or_FD_t8x8=0;
  currMB->written_SD_Coding_on_off=0;
#endif

  if (img->MbaffFrameFlag)
  {
    if (img->current_mb_nr%2)
    {
      topMB= &img->mb_data[img->current_mb_nr-1];
      if(!(img->type == B_SLICE))
        prevMbSkipped = (topMB->mb_type == 0);
      else 
        prevMbSkipped = topMB->skip_flag;
    }
    else 
      prevMbSkipped = 0;
  }


  if (img->current_mb_nr%2 == 0)
    currMB->mb_field = 0;
  else
    currMB->mb_field = img->mb_data[img->current_mb_nr-1].mb_field;
  
  
  currMB->qp = img->qp ;
  for (i=0; i<2; i++)
  {
    currMB->qpc[i] = Clip3 ( -img->bitdepth_chroma_qp_scale, 51, img->qp + dec_picture->chroma_qp_offset[i] );
    currMB->qpc[i] = currMB->qpc[i] < 0 ? currMB->qpc[i] : QP_SCALE_CR[currMB->qpc[i]];
  }
  
  currSE.type = SE_MBTYPE;
  
  //  read MB mode *****************************************************************
  dP = &(currSlice->partArr[partMap[currSE.type]]);
  
  if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)   currSE.mapping = linfo_ue;
  
  if(img->type == I_SLICE || img->type == SI_SLICE)
  {
    // read MB aff
    if (img->MbaffFrameFlag && img->current_mb_nr%2==0)
    {
      TRACE_STRING("mb_field_decoding_flag");
      if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
      {
        currSE.len = 1;
        readSyntaxElement_FLC(&currSE, dP->bitstream);
      }
      else
      {
        currSE.reading = readFieldModeInfo_CABAC;
        dP->readSyntaxElement(&currSE,img,inp,dP);
      }
      currMB->mb_field = currSE.value1;
    }
    if(active_pps->entropy_coding_mode_flag  == CABAC)
      CheckAvailabilityOfNeighborsCABAC();
    
    //  read MB type
    TRACE_STRING("mb_type");
    currSE.reading = readMB_typeInfo_CABAC;
    dP->readSyntaxElement(&currSE,img,inp,dP);
    
    currMB->mb_type = currSE.value1;
    if(!dP->bitstream->ei_flag)
      currMB->ei_flag = 0;
  }
  // non I/SI-slice CABAC
  else if (active_pps->entropy_coding_mode_flag == CABAC)
  {


    /********************************************
     *MB type for mode 32x32, 32x16, 16x32, 16x16
     *******************************************/
    if(mbi==0)
    {
      // read skip flag
      CheckAvailabilityOfNeighborsCABAC();
      TRACE_STRING("mb_skip_flag");
      currSE.reading = readMB_skip_flagInfo_CABAC;
      dP->readSyntaxElement(&currSE,img,inp,dP);
      
      currMB->mb_type   = currSE.value1;
      currMB->skip_flag = !(currSE.value1);
      
      if (img->type==B_SLICE)
        currMB->cbp = currSE.value2;
      
      if(!dP->bitstream->ei_flag)
        currMB->ei_flag = 0;
      
      if ((img->type==B_SLICE) && currSE.value1==0 && currSE.value2==0)
        img->cod_counter=0;

      // read MB AFF
      if (img->MbaffFrameFlag) 
      {
        check_bottom=read_bottom=read_top=0;
        if (img->current_mb_nr%2==0)
        {
          check_bottom =  currMB->skip_flag;
          read_top = !check_bottom;
        }
        else
        {
          read_bottom = (topMB->skip_flag && (!currMB->skip_flag));
        }
        
        if (read_bottom || read_top)
        {
          TRACE_STRING("mb_field_decoding_flag");
          currSE.reading = readFieldModeInfo_CABAC;
          dP->readSyntaxElement(&currSE,img,inp,dP);
          currMB->mb_field = currSE.value1;
        }
        if (check_bottom)
          check_next_mb_and_get_field_mode_CABAC(&currSE,img,inp,dP);
        
      }

      CheckAvailabilityOfNeighborsCABAC();

      // read MB type
      if (currMB->mb_type != 0 )
      {
        currSE.reading = readMB_typeInfo_CABAC32;
        TRACE_STRING("mb_type");
        dP->readSyntaxElement(&currSE,img,inp,dP);
        currMB->mb_type = currSE.value1;
        if(!dP->bitstream->ei_flag)
          currMB->ei_flag = 0;
      }

      if ((img->type==P_SLICE ))    // inter frame
        interpret_mb_mode_P(img);
      else if (img->type==I_SLICE)                                  // intra frame
        interpret_mb_mode_I(img);
      else if ((img->type==B_SLICE))       // B frame
        interpret_mb_mode_B(img);
      else if ((img->type==SP_SLICE))     // SP frame
        interpret_mb_mode_P(img);
      else if (img->type==SI_SLICE)     // SI frame
        interpret_mb_mode_SI(img);

      MB32.mb_type = currMB->mb_type; //untill this line, MB32.mb_type=1 could mean mode 1 of P or direct mode of B
      MB32.cbp     = currMB->cbp;
      memcpy(MB32.b8mode, currMB->b8mode, sizeof(int)*4);
      memcpy(MB32.b8pdir, currMB->b8pdir, sizeof(int)*4);

      assert((MB32.mb_type >=0 && MB32.mb_type<=3)||MB32.mb_type==P8x8);

      
      if(MB32.mb_type == P8x8)
        mb32=0;
      else
        mb32=1;
    }
    else
    {
      if( MB32.mb_type == P8x8)
        mb32=0;
      else
        mb32=1;
    }





    // read MB skip_flag
    if (img->MbaffFrameFlag && (img->current_mb_nr%2 == 0||prevMbSkipped))
      field_flag_inference();
    
    if(mb32==0)
    {
      CheckAvailabilityOfNeighborsCABAC();
      TRACE_STRING("mb_skip_flag");
      currSE.reading = readMB_skip_flagInfo_CABAC;
      dP->readSyntaxElement(&currSE,img,inp,dP);
      
      currMB->mb_type   = currSE.value1;
      currMB->skip_flag = !(currSE.value1);
      
      if (img->type==B_SLICE)
        currMB->cbp = currSE.value2;
      
      if(!dP->bitstream->ei_flag)
        currMB->ei_flag = 0;
      
      if ((img->type==B_SLICE) && currSE.value1==0 && currSE.value2==0)
        img->cod_counter=0;
      
      // read MB AFF
      if (img->MbaffFrameFlag) 
      {
        check_bottom=read_bottom=read_top=0;
        if (img->current_mb_nr%2==0)
        {
          check_bottom =  currMB->skip_flag;
          read_top = !check_bottom;
        }
        else
        {
          read_bottom = (topMB->skip_flag && (!currMB->skip_flag));
        }
        
        if (read_bottom || read_top)
        {
          TRACE_STRING("mb_field_decoding_flag");
          currSE.reading = readFieldModeInfo_CABAC;
          dP->readSyntaxElement(&currSE,img,inp,dP);
          currMB->mb_field = currSE.value1;
        }
        if (check_bottom)
          check_next_mb_and_get_field_mode_CABAC(&currSE,img,inp,dP);
        
      }
    }

    CheckAvailabilityOfNeighborsCABAC();
    if(mb32==0)
    {    
      // read MB type
      if (currMB->mb_type != 0 )
      {
        currSE.reading = readMB_typeInfo_CABAC;
        TRACE_STRING("mb_type");
        dP->readSyntaxElement(&currSE,img,inp,dP);
        currMB->mb_type = currSE.value1;
        if(!dP->bitstream->ei_flag)
          currMB->ei_flag = 0;
      }
    }
  }
  // VLC Non-Intra
  else
  {
    if(img->cod_counter == -1)
    {
      TRACE_STRING("mb_skip_run");
      dP->readSyntaxElement(&currSE,img,inp,dP);
      img->cod_counter = currSE.value1;
    }
    if (img->cod_counter==0)
    {
      // read MB aff
      if ((img->MbaffFrameFlag) && ((img->current_mb_nr%2==0) || ((img->current_mb_nr%2) && prevMbSkipped)))
      {
        TRACE_STRING("mb_field_decoding_flag");
        currSE.len = 1;
        readSyntaxElement_FLC(&currSE, dP->bitstream);
        currMB->mb_field = currSE.value1;
      }
      
      // read MB type
      TRACE_STRING("mb_type");
      dP->readSyntaxElement(&currSE,img,inp,dP);
      if(img->type == P_SLICE || img->type == SP_SLICE)
        currSE.value1++;
      currMB->mb_type = currSE.value1;
      if(!dP->bitstream->ei_flag)
        currMB->ei_flag = 0;
      img->cod_counter--;
      currMB->skip_flag = 0;
    } 
    else
    {
      img->cod_counter--;
      currMB->mb_type = 0;
      currMB->ei_flag = 0;
      currMB->skip_flag = 1;
      
      // read field flag of bottom block
      if(img->MbaffFrameFlag)
      {
        if(img->cod_counter == 0 && (img->current_mb_nr%2 == 0))
        {
          TRACE_STRING("mb_field_decoding_flag (of coded bottom mb)");
          currSE.len = 1;
          readSyntaxElement_FLC(&currSE, dP->bitstream);
          dP->bitstream->frame_bitoffset--;
          currMB->mb_field = currSE.value1;
        }
        else if(img->cod_counter > 0 && (img->current_mb_nr%2 == 0))
        {
          // check left macroblock pair first
          if (mb_is_available(img->current_mb_nr-2, img->current_mb_nr)&&((img->current_mb_nr%(img->PicWidthInMbs*2))!=0))
          {
            currMB->mb_field = img->mb_data[img->current_mb_nr-2].mb_field;
          }
          else
          {
            // check top macroblock pair
            if (mb_is_available(img->current_mb_nr-2*img->PicWidthInMbs, img->current_mb_nr))
            {
              currMB->mb_field = img->mb_data[img->current_mb_nr-2*img->PicWidthInMbs].mb_field;
            }
            else
              currMB->mb_field = 0;
          }
        }
      }
    }
  }
  
  dec_picture->mb_field[img->current_mb_nr] = currMB->mb_field;
  
  img->siblock[img->mb_x][img->mb_y]=0;
  
  if(mb32==0)
  {

  if ((img->type==P_SLICE ))    // inter frame
    interpret_mb_mode_P(img);
  else if (img->type==I_SLICE)                                  // intra frame
    interpret_mb_mode_I(img);
  else if ((img->type==B_SLICE))       // B frame
    interpret_mb_mode_B(img);
  else if ((img->type==SP_SLICE))     // SP frame
    interpret_mb_mode_P(img);
  else if (img->type==SI_SLICE)     // SI frame
    interpret_mb_mode_SI(img);
  }

#ifdef ADAPTIVE_QUANTIZATION
  if( currMB->cbp==0 || currMB->mb_type==IPCM) currMB->mb_iaqms_idx = 0;
#endif  
  if(img->MbaffFrameFlag)
  {
    if(currMB->mb_field)
    {
      img->num_ref_idx_l0_active <<=1;
      img->num_ref_idx_l1_active <<=1;
    }
  }
  
  //init NoMbPartLessThan8x8Flag
  currMB->NoMbPartLessThan8x8Flag = (IS_DIRECT(currMB) && !(active_sps->direct_8x8_inference_flag))? 0: 1;
  
  //====== READ 8x8 SUB-PARTITION MODES (modes of 8x8 blocks) and Intra VBST block modes ======
  if (IS_P8x8 (currMB))
  {
    currSE.type    = SE_MBTYPE;
    dP = &(currSlice->partArr[partMap[SE_MBTYPE]]);
    
    for (i=0; i<4; i++)
    {
      if (active_pps->entropy_coding_mode_flag ==UVLC || dP->bitstream->ei_flag) currSE.mapping = linfo_ue;
      else                                                  currSE.reading = readB8_typeInfo_CABAC;
      
      TRACE_STRING("sub_mb_type");
      dP->readSyntaxElement (&currSE, img, inp, dP);
      SetB8Mode (img, currMB, currSE.value1, i);
      
      //set NoMbPartLessThan8x8Flag for P8x8 mode
      currMB->NoMbPartLessThan8x8Flag &= (currMB->b8mode[i]==0 && active_sps->direct_8x8_inference_flag) || 
        (currMB->b8mode[i]==4);
    }
    //--- init macroblock data ---
    init_macroblock       (img);
    readMotionInfoFromNAL (img, inp);
  }
  
  
  //============= Transform Size Flag for INTRA MBs =============
  //-------------------------------------------------------------
  //transform size flag for INTRA_4x4 and INTRA_8x8 modes
  if (currMB->mb_type == I4MB && img->Transform8x8Mode)
  {
    currSE.type   =  SE_HEADER;
    dP = &(currSlice->partArr[partMap[SE_HEADER]]);
    currSE.reading = readMB_transform_size_flag_CABAC;
    TRACE_STRING("transform size 8x8 flag");
    
    // read UVLC transform_size_8x8_flag
    if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
    {
      currSE.len = 1;
      readSyntaxElement_FLC(&currSE, dP->bitstream);
    } 
    else 
    {
      dP->readSyntaxElement(&currSE,img,inp,dP);
    } 
    
    currMB->luma_transform_size_8x8_flag = currSE.value1;
    
    if (currMB->luma_transform_size_8x8_flag)
    {
      currMB->mb_type = I8MB;
      for (i=0;i<4;i++) 
      {
        currMB->b8mode[i]=I8MB;
        currMB->b8pdir[i]=-1;
      }
    }
  }
  else
  {
    currMB->luma_transform_size_8x8_flag = 0;
  }
  
  if(active_pps->constrained_intra_pred_flag && (img->type==P_SLICE|| img->type==B_SLICE))        // inter frame
  {
    if( !IS_INTRA(currMB) )
    {
      img->intra_block[img->current_mb_nr] = 0;
    }
  }
  
  //! TO for error concealment
  //! If we have an INTRA Macroblock and we lost the partition
  //! which contains the intra coefficients Copy MB would be better 
  //! than just a gray block.
  //! Seems to be a bit at the wrong place to do this right here, but for this case 
  //! up to now there is no other way.
  dP = &(currSlice->partArr[partMap[SE_CBP_INTRA]]);
  if(IS_INTRA (currMB) && dP->bitstream->ei_flag && img->number)
  {
    currMB->mb_type = 0;
    currMB->ei_flag = 1;
    for (i=0;i<4;i++)
    {
      currMB->b8mode[i]=currMB->b8pdir[i]=0;
    }
  }
  dP = &(currSlice->partArr[partMap[currSE.type]]);
  //! End TO
  
  
  //--- init macroblock data ---
  if (!IS_P8x8 (currMB))
  {
    if(mb32==0)
    {
      init_macroblock       (img);
    }
    else if(mb32==1 && mbi==0)
    {
      init_macroblock32     (img, 1);
    }
  }
  
  if (IS_DIRECT (currMB) && img->cod_counter >= 0)
  {
    currMB->cbp = 0;
    if(mb32==1)
    {
      if(mbi==0)
        reset_coeffs(); //should we have a specific function for 32x32 ???
    }
    else
    reset_coeffs();
    
    if (active_pps->entropy_coding_mode_flag ==CABAC)
      img->cod_counter=-1;
    
    return DECODE_MB;
  }
  if (((mb32 && mbi==0) || !mb32) && IS_COPY (currMB) )
  {
    int i, j, k;
    short pmv[2];
    int zeroMotionAbove;
    int zeroMotionLeft;
    PixelPos mb_a, mb_b;
    int      a_mv_y = 0;
    int      a_ref_idx = 0;
    int      b_mv_y = 0;
    int      b_ref_idx = 0;
    int      list_offset = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;
    
    getLuma4x4Neighbour(img->current_mb_nr,0,0,-1, 0,&mb_a);
    getLuma4x4Neighbour(img->current_mb_nr,0,0, 0,-1,&mb_b);
    
    if (mb_a.available)
    {
      a_mv_y    = dec_picture->mv[LIST_0][mb_a.pos_y][mb_a.pos_x][1];
      a_ref_idx = dec_picture->ref_idx[LIST_0][mb_a.pos_y][mb_a.pos_x];
      
      if (currMB->mb_field && !img->mb_data[mb_a.mb_addr].mb_field)
      {
        a_mv_y    /=2;
        a_ref_idx *=2;
      }
      if (!currMB->mb_field && img->mb_data[mb_a.mb_addr].mb_field)
      {
        a_mv_y    *=2;
        a_ref_idx >>=1;
      }
    }
    
    if (mb_b.available)
    {
      b_mv_y    = dec_picture->mv[LIST_0][mb_b.pos_y][mb_b.pos_x][1];
      b_ref_idx = dec_picture->ref_idx[LIST_0][mb_b.pos_y][mb_b.pos_x];
      
      if (currMB->mb_field && !img->mb_data[mb_b.mb_addr].mb_field)
      {
        b_mv_y    /=2;
        b_ref_idx *=2;
      }
      if (!currMB->mb_field && img->mb_data[mb_b.mb_addr].mb_field)
      {
        b_mv_y    *=2;
        b_ref_idx >>=1;
      }
    }
    
    zeroMotionLeft  = !mb_a.available ? 1 : a_ref_idx==0 && dec_picture->mv[LIST_0][mb_a.pos_y][mb_a.pos_x][0]==0 && a_mv_y==0 ? 1 : 0;
    zeroMotionAbove = !mb_b.available ? 1 : b_ref_idx==0 && dec_picture->mv[LIST_0][mb_b.pos_y][mb_b.pos_x][0]==0 && b_mv_y==0 ? 1 : 0;
    
    currMB->cbp = 0;
    reset_coeffs();
    
    img_block_y   = img->block_y;
    
#ifdef MV_COMPETITION
#ifdef MB32X32_MVC
    if (mv_comp.mv_competition > 0)
    {
      short prediction_mode_for_skip;
      short predictor;

//    32x32 SKIP   ---rahul  
      if(mb32 && mbi==0)
      {
        for (prediction_mode_for_skip = 0; prediction_mode_for_skip<mv_comp.nb_mode_for_skip; prediction_mode_for_skip++)
        {
  //      SetMotionVectorPredictor_Skip  (pmv, pmv+1, dec_picture->ref_idx, dec_picture->mv, 0, LIST_0, 0, 0, 16, 16,prediction_mode_for_skip);
          SetMotionVectorPredictor_Skip32(pmv, pmv+1, dec_picture->ref_idx, dec_picture->mv, 0, LIST_0, 0, 0, 32, 32,prediction_mode_for_skip, 1);
          mv_comp.mv_pred_skip[prediction_mode_for_skip][0] = pmv[0];
          mv_comp.mv_pred_skip[prediction_mode_for_skip][1] = pmv[1];
        }
        if ((mv_comp.nb_mode_for_skip> 1) && (read_index_for_skip_mode() == TRUE))
          predictor = readPredictorForSkip(img,inp);
        else 
          predictor = 0;

        for(i=0;i<BLOCK_SIZE*2;i++)
          for(j=0;j<BLOCK_SIZE*2;j++)
            for (k=0;k<2;k++)
              dec_picture->mv[LIST_0][img_block_y+j][img->block_x+i][k] = mv_comp.mv_pred_skip[predictor][k];
      }
//    16x16 skip inside 32x32  ---rahul
      if(mb32==0)
      {
        for (prediction_mode_for_skip = 0; prediction_mode_for_skip<mv_comp.nb_mode_for_skip; prediction_mode_for_skip++)
        {
          SetMotionVectorPredictor_Skip  (pmv, pmv+1, dec_picture->ref_idx, dec_picture->mv, 0, LIST_0, 0, 0, 16, 16,prediction_mode_for_skip);
          mv_comp.mv_pred_skip[prediction_mode_for_skip][0] = pmv[0];
          mv_comp.mv_pred_skip[prediction_mode_for_skip][1] = pmv[1];
        }

        if ((mv_comp.nb_mode_for_skip> 1) && (read_index_for_skip_mode() == TRUE))
          predictor = readPredictorForSkip(img,inp);
        else 
          predictor = 0;

        for(i=0;i<BLOCK_SIZE;i++)
          for(j=0;j<BLOCK_SIZE;j++)
            for (k=0;k<2;k++)
              dec_picture->mv[LIST_0][img_block_y+j][img->block_x+i][k] = mv_comp.mv_pred_skip[predictor][k];
      }
    }
#else
    if (mv_comp.mv_competition > 0 && !mb32)
    {
      short prediction_mode_for_skip;
      short predictor;
      
      for (prediction_mode_for_skip = 0; prediction_mode_for_skip<mv_comp.nb_mode_for_skip; prediction_mode_for_skip++)
      {
        SetMotionVectorPredictor_Skip(pmv, pmv+1, dec_picture->ref_idx, dec_picture->mv, 0, LIST_0, 0, 0, 16, 16,prediction_mode_for_skip);
        mv_comp.mv_pred_skip[prediction_mode_for_skip][0] = pmv[0];
        mv_comp.mv_pred_skip[prediction_mode_for_skip][1] = pmv[1];
      }
      
      if ((mv_comp.nb_mode_for_skip> 1) && (read_index_for_skip_mode() == TRUE))
        predictor = readPredictorForSkip(img,inp);
      else 
        predictor = 0;
      
      for(i=0;i<BLOCK_SIZE;i++)
        for(j=0;j<BLOCK_SIZE;j++)
          for (k=0;k<2;k++)
            dec_picture->mv[LIST_0][img_block_y+j][img->block_x+i][k] = mv_comp.mv_pred_skip[predictor][k];
    }
#endif
    else
#endif
      
    {
      
      
      if (zeroMotionAbove || zeroMotionLeft)
      {
        if(mb32)
        {
          if(mbi==0)
          {
            for(i=0;i<BLOCK_SIZE*2;i++)
              for(j=0;j<BLOCK_SIZE*2;j++)
                for (k=0;k<2;k++)
                  dec_picture->mv[LIST_0][img->block_y+j][img->block_x+i][k] = 0;
          }
        }
        else
        {
        for(i=0;i<BLOCK_SIZE;i++)
          for(j=0;j<BLOCK_SIZE;j++)
            for (k=0;k<2;k++)
              dec_picture->mv[LIST_0][img->block_y+j][img->block_x+i][k] = 0;
        }
      }
      else
      {
        if(mb32)
        {
          if(mbi==0)
          {
            SetMotionVectorPredictor32 (img, pmv, pmv+1, 0, LIST_0, dec_picture->ref_idx, dec_picture->mv, 0, 0, 32, 32, 1);
            for(i=0;i<BLOCK_SIZE*2;i++)
              for(j=0;j<BLOCK_SIZE*2;j++)
                for (k=0;k<2;k++)
                {
                  dec_picture->mv[LIST_0][img_block_y+j][img->block_x+i][k] = pmv[k];
                }
            }
        }
        else
        {
        SetMotionVectorPredictor (img, pmv, pmv+1, 0, LIST_0, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);
        
        for(i=0;i<BLOCK_SIZE;i++)
          for(j=0;j<BLOCK_SIZE;j++)
            for (k=0;k<2;k++)
            {
              dec_picture->mv[LIST_0][img_block_y+j][img->block_x+i][k] = pmv[k];
            }
        }
      }
      
      
    }
    if(mb32) 
    {
      if(mbi==0)
      {
        for(i=0;i<BLOCK_SIZE*2;i++)
          for(j=0;j<BLOCK_SIZE*2;j++)
          {
            dec_picture->ref_idx[LIST_0][img_block_y+j][img->block_x+i] = 0;
            dec_picture->ref_pic_id[LIST_0][img_block_y+j][img->block_x+i] = 
              dec_picture->ref_pic_num[img->current_slice_nr][LIST_0 + list_offset][(short)dec_picture->ref_idx[LIST_0][img_block_y+j][img->block_x+i]];
          }        
      }
    }
    else
    {
    for(i=0;i<BLOCK_SIZE;i++)
      for(j=0;j<BLOCK_SIZE;j++)
      {
        dec_picture->ref_idx[LIST_0][img_block_y+j][img->block_x+i] = 0;
        dec_picture->ref_pic_id[LIST_0][img_block_y+j][img->block_x+i] = 
          dec_picture->ref_pic_num[img->current_slice_nr][LIST_0 + list_offset][(short)dec_picture->ref_idx[LIST_0][img_block_y+j][img->block_x+i]];
      }        
    }
      return DECODE_MB;
  }
  if(currMB->mb_type!=IPCM)
  {
    
    // intra prediction modes for a macroblock 4x4 **********************************************
    read_ipred_modes(img,inp);
    
    // read inter frame vector data *********************************************************
    if(mb32 && mbi==0)
    {
      if (MB32.mb_type != 0)
         readMotionInfoFromNAL32 (img, inp, 1);
    }
    else if(mb32==0 && IS_INTERMV (currMB) && (!IS_P8x8(currMB)))
    {
      readMotionInfoFromNAL (img, inp);
    }
    // read CBP and Coeffs  ***************************************************************

    if(mb32 && mbi==0)
    {      
      MB32.cbp=readCBP32(img,inp,!delta_qp_sent);//with dquant
    }
    if(!mb32 || (mb32 && MB32.cbp !=0 ) )
      readCBPandCoeffsFromNAL (img,inp, mb32,mbi);
  }
  else
  {
    //read pcm_alignment_zero_bit and pcm_byte[i] 
    
    // here dP is assigned with the same dP as SE_MBTYPE, because IPCM syntax is in the 
    // same category as MBTYPE
    dP = &(currSlice->partArr[partMap[SE_MBTYPE]]);
    readIPCMcoeffsFromNAL(img,inp,dP);
  }
  
#ifdef ADAPTIVE_FD_SD_CODING
  if (currMB->written_SD_Coding_on_off==0)
  {
    currMB->SD_Coding_on_off=0;
  }
#endif
  return DECODE_MB;
}

void read_mbtype(struct img_par *img, struct inp_par *inp, int mb_ext_level)//MB64X64
{
  Slice *currSlice = img->currentSlice;
  SyntaxElement currSE;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int *partMap = assignSE2partition[currSlice->dp_mode];
  DataPartition *dP;
  currSE.type = SE_MBTYPE;

  //  read MB mode *****************************************************************
  dP = &(currSlice->partArr[partMap[currSE.type]]);

      // read skip flag
      CheckAvailabilityOfNeighborsCABAC();
      TRACE_STRING("mb_skip_flag");
      currSE.reading = readMB_skip_flagInfo_CABAC;
      dP->readSyntaxElement(&currSE,img,inp,dP);
      
      currMB->mb_type   = currSE.value1;
      currMB->skip_flag = !(currSE.value1);
      
      if (img->type==B_SLICE)
        currMB->cbp = currSE.value2;
      
      if(!dP->bitstream->ei_flag)
        currMB->ei_flag = 0;
      
      if ((img->type==B_SLICE) && currSE.value1==0 && currSE.value2==0)
        img->cod_counter=0;


      CheckAvailabilityOfNeighborsCABAC();

      // read MB type
      if (currMB->mb_type != 0 )
      {
        if(mb_ext_level != 0)
          currSE.reading = readMB_typeInfo_CABAC32;
        else
          currSE.reading = readMB_typeInfo_CABAC;

        TRACE_STRING("mb_type");
        dP->readSyntaxElement(&currSE,img,inp,dP);
        currMB->mb_type = currSE.value1;
        if(!dP->bitstream->ei_flag)
          currMB->ei_flag = 0;
      }

      if ((img->type==P_SLICE ))    // inter frame
        interpret_mb_mode_P(img);
      else if (img->type==I_SLICE)                                  // intra frame
        interpret_mb_mode_I(img);
      else if ((img->type==B_SLICE))       // B frame
        interpret_mb_mode_B(img);
      else if ((img->type==SP_SLICE))     // SP frame
        interpret_mb_mode_P(img);
      else if (img->type==SI_SLICE)     // SI frame
        interpret_mb_mode_SI(img);

}
int read_one_macroblock64(struct img_par *img,struct inp_par *inp, int mb32count, int mb16count)
{
  int i;
  int mb64=0, mb32=0; // mb64:0- 32x32 mode, 1-the others; mb32:0-16x16 mode, 1-the others

  SyntaxElement currSE;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  
  Slice *currSlice = img->currentSlice;
  DataPartition *dP;
  int *partMap = assignSE2partition[currSlice->dp_mode];
  Macroblock *topMB = NULL;
  int  prevMbSkipped = 0;
  int  img_block_y;  

#ifdef ADAPTIVE_FD_SD_CODING
  int           j;
  for (i=0;i<16;i++)
  {
    for (j=0;j<16;j++)
    {
      currMB->quantizer_indices[j][i]=0;
    }
  }
  for (i=0;i<2;i++)
  {
    for (j=0;j<2;j++)
    {
      currMB->SD_or_FD     [j][i]=0;
    }
  }
  currMB->SD_or_FD_t8x8=0;
  currMB->written_SD_Coding_on_off=0;
#endif

  if (img->MbaffFrameFlag)
  {
    if (img->current_mb_nr%2)
    {
      topMB= &img->mb_data[img->current_mb_nr-1];
      if(!(img->type == B_SLICE))
        prevMbSkipped = (topMB->mb_type == 0);
      else 
        prevMbSkipped = topMB->skip_flag;
    }
    else 
      prevMbSkipped = 0;
  }

  if (img->current_mb_nr%2 == 0)
    currMB->mb_field = 0;
  else
    currMB->mb_field = img->mb_data[img->current_mb_nr-1].mb_field;
  
  
  currMB->qp = img->qp ;
  for (i=0; i<2; i++)
  {
    currMB->qpc[i] = Clip3 ( -img->bitdepth_chroma_qp_scale, 51, img->qp + dec_picture->chroma_qp_offset[i] );
    currMB->qpc[i] = currMB->qpc[i] < 0 ? currMB->qpc[i] : QP_SCALE_CR[currMB->qpc[i]];
  }
  
  currSE.type = SE_MBTYPE;
  
  //  read MB mode *****************************************************************
  dP = &(currSlice->partArr[partMap[currSE.type]]);
  
  if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)   currSE.mapping = linfo_ue;
  
  if(img->type == I_SLICE || img->type == SI_SLICE)
  {
    // read MB aff
    if (img->MbaffFrameFlag && img->current_mb_nr%2==0)
    {
      TRACE_STRING("mb_field_decoding_flag");
      if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
      {
        currSE.len = 1;
        readSyntaxElement_FLC(&currSE, dP->bitstream);
      }
      else
      {
        currSE.reading = readFieldModeInfo_CABAC;
        dP->readSyntaxElement(&currSE,img,inp,dP);
      }
      currMB->mb_field = currSE.value1;
    }
    if(active_pps->entropy_coding_mode_flag  == CABAC)
      CheckAvailabilityOfNeighborsCABAC();
    
    //  read MB type
    TRACE_STRING("mb_type");
    currSE.reading = readMB_typeInfo_CABAC;
    dP->readSyntaxElement(&currSE,img,inp,dP);
    
    currMB->mb_type = currSE.value1;
    if(!dP->bitstream->ei_flag)
      currMB->ei_flag = 0;
  }
  // non I/SI-slice CABAC
  else if (active_pps->entropy_coding_mode_flag == CABAC)
  {
    CheckAvailabilityOfNeighborsCABAC();


    /********************************************
     *MB type for mode 64x64, 64x32, 32x64, 32x32
     *******************************************/
    if(img->mb_ext_level == 2)//MB64X64
    {
      if(mb32count==0 && mb16count==0)
      {
        read_mbtype(img, inp, 2);

        MB64.mb_type = currMB->mb_type; //untill this line, MB32.mb_type=1 could mean mode 1 of P or direct mode of B
        MB64.cbp     = currMB->cbp;
        memcpy(MB64.b8mode, currMB->b8mode, sizeof(int)*4);
        memcpy(MB64.b8pdir, currMB->b8pdir, sizeof(int)*4);

        assert((MB64.mb_type >=0 && MB64.mb_type<=3)||MB64.mb_type==P8x8);      
      }

      if( MB64.mb_type == P8x8)
        mb64=0;
      else
        mb64=1;
    }

    /********************************************
     *MB type for mode 32x32, 32x16, 16x32, 16x16
     *******************************************/
    if(mb64==0)
    {
      if((mb16count&3)==0)
      {
        read_mbtype(img, inp, 1);

        MB32.mb_type = currMB->mb_type; //untill this line, MB32.mb_type=1 could mean mode 1 of P or direct mode of B
        MB32.cbp     = currMB->cbp;
        memcpy(MB32.b8mode, currMB->b8mode, sizeof(int)*4);
        memcpy(MB32.b8pdir, currMB->b8pdir, sizeof(int)*4);

        assert((MB32.mb_type >=0 && MB32.mb_type<=3)||MB32.mb_type==P8x8);        
      }
      if( MB32.mb_type == P8x8)
        mb32=0;
      else
        mb32=1;
    }

    
    if(mb64==0 && mb32==0)
    {
      read_mbtype(img, inp, 0);
    }
  }
  // VLC Non-Intra
  else
  {
    if(img->cod_counter == -1)
    {
      TRACE_STRING("mb_skip_run");
      dP->readSyntaxElement(&currSE,img,inp,dP);
      img->cod_counter = currSE.value1;
    }
    if (img->cod_counter==0)
    {
      // read MB aff
      if ((img->MbaffFrameFlag) && ((img->current_mb_nr%2==0) || ((img->current_mb_nr%2) && prevMbSkipped)))
      {
        TRACE_STRING("mb_field_decoding_flag");
        currSE.len = 1;
        readSyntaxElement_FLC(&currSE, dP->bitstream);
        currMB->mb_field = currSE.value1;
      }
      
      // read MB type
      TRACE_STRING("mb_type");
      dP->readSyntaxElement(&currSE,img,inp,dP);
      if(img->type == P_SLICE || img->type == SP_SLICE)
        currSE.value1++;
      currMB->mb_type = currSE.value1;
      if(!dP->bitstream->ei_flag)
        currMB->ei_flag = 0;
      img->cod_counter--;
      currMB->skip_flag = 0;
    } 
    else
    {
      img->cod_counter--;
      currMB->mb_type = 0;
      currMB->ei_flag = 0;
      currMB->skip_flag = 1;
      
      // read field flag of bottom block
      if(img->MbaffFrameFlag)
      {
        if(img->cod_counter == 0 && (img->current_mb_nr%2 == 0))
        {
          TRACE_STRING("mb_field_decoding_flag (of coded bottom mb)");
          currSE.len = 1;
          readSyntaxElement_FLC(&currSE, dP->bitstream);
          dP->bitstream->frame_bitoffset--;
          currMB->mb_field = currSE.value1;
        }
        else if(img->cod_counter > 0 && (img->current_mb_nr%2 == 0))
        {
          // check left macroblock pair first
          if (mb_is_available(img->current_mb_nr-2, img->current_mb_nr)&&((img->current_mb_nr%(img->PicWidthInMbs*2))!=0))
          {
            currMB->mb_field = img->mb_data[img->current_mb_nr-2].mb_field;
          }
          else
          {
            // check top macroblock pair
            if (mb_is_available(img->current_mb_nr-2*img->PicWidthInMbs, img->current_mb_nr))
            {
              currMB->mb_field = img->mb_data[img->current_mb_nr-2*img->PicWidthInMbs].mb_field;
            }
            else
              currMB->mb_field = 0;
          }
        }
      }
    }
      if ((img->type==P_SLICE ))    // inter frame
        interpret_mb_mode_P(img);
      else if (img->type==I_SLICE)                                  // intra frame
        interpret_mb_mode_I(img);
      else if ((img->type==B_SLICE))       // B frame
        interpret_mb_mode_B(img);
      else if ((img->type==SP_SLICE))     // SP frame
        interpret_mb_mode_P(img);
      else if (img->type==SI_SLICE)     // SI frame
        interpret_mb_mode_SI(img);
    
  }
  dec_picture->mb_field[img->current_mb_nr] = currMB->mb_field;
    
  img->siblock[img->mb_x][img->mb_y]=0;    
  

#ifdef ADAPTIVE_QUANTIZATION
  if( currMB->cbp==0 || currMB->mb_type==IPCM) currMB->mb_iaqms_idx = 0;
#endif  
  if(img->MbaffFrameFlag)
  {
    if(currMB->mb_field)
    {
      img->num_ref_idx_l0_active <<=1;
      img->num_ref_idx_l1_active <<=1;
    }
  }
  
  //init NoMbPartLessThan8x8Flag
  currMB->NoMbPartLessThan8x8Flag = (IS_DIRECT(currMB) && !(active_sps->direct_8x8_inference_flag))? 0: 1;
  
  //====== READ 8x8 SUB-PARTITION MODES (modes of 8x8 blocks) and Intra VBST block modes ======
  if (IS_P8x8 (currMB))
  {
    currSE.type    = SE_MBTYPE;
    dP = &(currSlice->partArr[partMap[SE_MBTYPE]]);
    
    for (i=0; i<4; i++)
    {
      if (active_pps->entropy_coding_mode_flag ==UVLC || dP->bitstream->ei_flag) currSE.mapping = linfo_ue;
      else                                                  currSE.reading = readB8_typeInfo_CABAC;
      
      TRACE_STRING("sub_mb_type");
      dP->readSyntaxElement (&currSE, img, inp, dP);
      SetB8Mode (img, currMB, currSE.value1, i);
      
      //set NoMbPartLessThan8x8Flag for P8x8 mode
      currMB->NoMbPartLessThan8x8Flag &= (currMB->b8mode[i]==0 && active_sps->direct_8x8_inference_flag) || 
        (currMB->b8mode[i]==4);
    }
    //--- init macroblock data ---
    init_macroblock       (img);
    readMotionInfoFromNAL (img, inp);
  }
  
  
  //============= Transform Size Flag for INTRA MBs =============
  //-------------------------------------------------------------
  //transform size flag for INTRA_4x4 and INTRA_8x8 modes
  if (currMB->mb_type == I4MB && img->Transform8x8Mode)
  {
    currSE.type   =  SE_HEADER;
    dP = &(currSlice->partArr[partMap[SE_HEADER]]);
    currSE.reading = readMB_transform_size_flag_CABAC;
    TRACE_STRING("transform size 8x8 flag");
    
    // read UVLC transform_size_8x8_flag
    if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
    {
      currSE.len = 1;
      readSyntaxElement_FLC(&currSE, dP->bitstream);
    } 
    else 
    {
      dP->readSyntaxElement(&currSE,img,inp,dP);
    } 
    
    currMB->luma_transform_size_8x8_flag = currSE.value1;
    
    if (currMB->luma_transform_size_8x8_flag)
    {
      currMB->mb_type = I8MB;
      for (i=0;i<4;i++) 
      {
        currMB->b8mode[i]=I8MB;
        currMB->b8pdir[i]=-1;
      }
    }
  }
  else
  {
    currMB->luma_transform_size_8x8_flag = 0;
  }
  
  if(active_pps->constrained_intra_pred_flag && (img->type==P_SLICE|| img->type==B_SLICE))        // inter frame
  {
    if( !IS_INTRA(currMB) )
    {
      img->intra_block[img->current_mb_nr] = 0;
    }
  }
  
  //! TO for error concealment
  //! If we have an INTRA Macroblock and we lost the partition
  //! which contains the intra coefficients Copy MB would be better 
  //! than just a gray block.
  //! Seems to be a bit at the wrong place to do this right here, but for this case 
  //! up to now there is no other way.
  dP = &(currSlice->partArr[partMap[SE_CBP_INTRA]]);
  if(IS_INTRA (currMB) && dP->bitstream->ei_flag && img->number)
  {
    currMB->mb_type = 0;
    currMB->ei_flag = 1;
    for (i=0;i<4;i++)
    {
      currMB->b8mode[i]=currMB->b8pdir[i]=0;
    }
  }
  dP = &(currSlice->partArr[partMap[currSE.type]]);
  //! End TO
  
  
  //--- init macroblock data ---
  if (!IS_P8x8 (currMB))
  {
    if(mb64==0 && mb32==0)
    {
      init_macroblock       (img);
    }
    else if(mb64==0 && mb32==1 && (mb16count&3)==0)
    {
      init_macroblock32     (img, 1);
    }
    else if(mb64==1 && mb32count == 0 && mb16count==0)
    {
      init_macroblock32     (img, 2);
    }
  }
  
  if (IS_DIRECT (currMB) && img->cod_counter >= 0)
  {
    currMB->cbp = 0;
    reset_coeffs();//reset_coeffs has been called in image.c 
    
    if (active_pps->entropy_coding_mode_flag ==CABAC)
      img->cod_counter=-1;
    
    return DECODE_MB;
  }
  if ( ((mb64 && mb32count==0 && mb16count==0) || (!mb64 && mb32 && (mb16count&3)==0) || (!mb64 &&!mb32)) 
    && IS_COPY (currMB) )
  {
    int i, j, k;
    short pmv[2];
    int zeroMotionAbove;
    int zeroMotionLeft;
    PixelPos mb_a, mb_b;
    int      a_mv_y = 0;
    int      a_ref_idx = 0;
    int      b_mv_y = 0;
    int      b_ref_idx = 0;
    int      list_offset = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;
    
    getLuma4x4Neighbour(img->current_mb_nr,0,0,-1, 0,&mb_a);
    getLuma4x4Neighbour(img->current_mb_nr,0,0, 0,-1,&mb_b);
    
    if (mb_a.available)
    {
      a_mv_y    = dec_picture->mv[LIST_0][mb_a.pos_y][mb_a.pos_x][1];
      a_ref_idx = dec_picture->ref_idx[LIST_0][mb_a.pos_y][mb_a.pos_x];
      
      if (currMB->mb_field && !img->mb_data[mb_a.mb_addr].mb_field)
      {
        a_mv_y    /=2;
        a_ref_idx *=2;
      }
      if (!currMB->mb_field && img->mb_data[mb_a.mb_addr].mb_field)
      {
        a_mv_y    *=2;
        a_ref_idx >>=1;
      }
    }
    
    if (mb_b.available)
    {
      b_mv_y    = dec_picture->mv[LIST_0][mb_b.pos_y][mb_b.pos_x][1];
      b_ref_idx = dec_picture->ref_idx[LIST_0][mb_b.pos_y][mb_b.pos_x];
      
      if (currMB->mb_field && !img->mb_data[mb_b.mb_addr].mb_field)
      {
        b_mv_y    /=2;
        b_ref_idx *=2;
      }
      if (!currMB->mb_field && img->mb_data[mb_b.mb_addr].mb_field)
      {
        b_mv_y    *=2;
        b_ref_idx >>=1;
      }
    }
    
    zeroMotionLeft  = !mb_a.available ? 1 : a_ref_idx==0 && dec_picture->mv[LIST_0][mb_a.pos_y][mb_a.pos_x][0]==0 && a_mv_y==0 ? 1 : 0;
    zeroMotionAbove = !mb_b.available ? 1 : b_ref_idx==0 && dec_picture->mv[LIST_0][mb_b.pos_y][mb_b.pos_x][0]==0 && b_mv_y==0 ? 1 : 0;
    
    currMB->cbp = 0;
    reset_coeffs();
    
    img_block_y   = img->block_y;
    
#ifdef MV_COMPETITION
#ifdef MB32X32_MVC
    if (mv_comp.mv_competition > 0)
    {
      short prediction_mode_for_skip;
      short predictor;

      if(((mb64 && mb32count==0 && mb16count==0) || (!mb64 && mb32 && (mb16count&3)==0) || (!mb64 && !mb32)))
      {

      for (prediction_mode_for_skip = 0; prediction_mode_for_skip<mv_comp.nb_mode_for_skip; prediction_mode_for_skip++)
      {
//      SetMotionVectorPredictor_Skip  (pmv, pmv+1, dec_picture->ref_idx, dec_picture->mv, 0, LIST_0, 0, 0, 16, 16,prediction_mode_for_skip);
        if(mb32 || mb64)
          SetMotionVectorPredictor_Skip32(pmv, pmv+1, dec_picture->ref_idx, dec_picture->mv, 0, LIST_0, 0, 0, 16<<(mb64*2+mb32), 16<<(mb64*2+mb32), prediction_mode_for_skip, (mb64*2+mb32));
        else
          SetMotionVectorPredictor_Skip  (pmv, pmv+1, dec_picture->ref_idx, dec_picture->mv, 0, LIST_0, 0, 0, 16, 16,prediction_mode_for_skip);

        mv_comp.mv_pred_skip[prediction_mode_for_skip][0] = pmv[0];
        mv_comp.mv_pred_skip[prediction_mode_for_skip][1] = pmv[1];
      }

        if ((mv_comp.nb_mode_for_skip> 1) && (read_index_for_skip_mode() == TRUE))
          predictor = readPredictorForSkip(img,inp);
        else 
          predictor = 0;

        if(mb32 || mb64)
        {
          for(i=0;i<(BLOCK_SIZE<<(mb64*2+mb32));i++)
            for(j=0;j<(BLOCK_SIZE<<(mb64*2+mb32));j++)
              for (k=0;k<2;k++)
                dec_picture->mv[LIST_0][img_block_y+j][img->block_x+i][k] = mv_comp.mv_pred_skip[predictor][k];
        }
        else
        {
          for(i=0;i<BLOCK_SIZE;i++)
            for(j=0;j<BLOCK_SIZE;j++)
              for (k=0;k<2;k++)
                dec_picture->mv[LIST_0][img_block_y+j][img->block_x+i][k] = mv_comp.mv_pred_skip[predictor][k];
        }

      }
    }
#else
    
    if (mv_comp.mv_competition > 0 && (!mb64 &&!mb32))
      
    {
      short prediction_mode_for_skip;
      short predictor;
      
      for (prediction_mode_for_skip = 0; prediction_mode_for_skip<mv_comp.nb_mode_for_skip; prediction_mode_for_skip++)
      {
        SetMotionVectorPredictor_Skip(pmv, pmv+1, dec_picture->ref_idx, dec_picture->mv, 0, LIST_0, 0, 0, 16, 16,prediction_mode_for_skip);
        mv_comp.mv_pred_skip[prediction_mode_for_skip][0] = pmv[0];
        mv_comp.mv_pred_skip[prediction_mode_for_skip][1] = pmv[1];
      }
      
      if ((mv_comp.nb_mode_for_skip> 1) && (read_index_for_skip_mode() == TRUE))
        predictor = readPredictorForSkip(img,inp);
      else 
        predictor = 0;
      
      for(i=0;i<BLOCK_SIZE;i++)
        for(j=0;j<BLOCK_SIZE;j++)
          for (k=0;k<2;k++)
            dec_picture->mv[LIST_0][img_block_y+j][img->block_x+i][k] = mv_comp.mv_pred_skip[predictor][k];
    }
#endif
    else
#endif
      
    {
      
      
      if (zeroMotionAbove || zeroMotionLeft)
      {
        for(i=0;i<(BLOCK_SIZE<<(mb64*2+mb32));i++)
          for(j=0;j<(BLOCK_SIZE<<(mb64*2+mb32));j++)
            for (k=0;k<2;k++)
              dec_picture->mv[LIST_0][img->block_y+j][img->block_x+i][k] = 0;
      }
      else
      {
        if(mb32 || mb64)
        {
          SetMotionVectorPredictor32 (img, pmv, pmv+1, 0, LIST_0, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16<<(mb64*2+mb32), 16<<(mb64*2+mb32), (mb64*2+mb32));
          for(i=0;i<(BLOCK_SIZE<<(mb64*2+mb32));i++)
            for(j=0;j<(BLOCK_SIZE<<(mb64*2+mb32));j++)
              for (k=0;k<2;k++)
              {
                dec_picture->mv[LIST_0][img_block_y+j][img->block_x+i][k] = pmv[k];
              }
        }
        else
        {
        SetMotionVectorPredictor (img, pmv, pmv+1, 0, LIST_0, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);
        
        for(i=0;i<BLOCK_SIZE;i++)
          for(j=0;j<BLOCK_SIZE;j++)
            for (k=0;k<2;k++)
            {
              dec_picture->mv[LIST_0][img_block_y+j][img->block_x+i][k] = pmv[k];
            }
        }
      }
      
      
    }
    for(i=0;i<(BLOCK_SIZE<<(mb64*2+mb32));i++)
      for(j=0;j<(BLOCK_SIZE<<(mb64*2+mb32));j++)
      {
        dec_picture->ref_idx[LIST_0][img_block_y+j][img->block_x+i] = 0;
        dec_picture->ref_pic_id[LIST_0][img_block_y+j][img->block_x+i] = 
          dec_picture->ref_pic_num[img->current_slice_nr][LIST_0 + list_offset][(short)dec_picture->ref_idx[LIST_0][img_block_y+j][img->block_x+i]];
      }              
      return DECODE_MB;
  }
  if(currMB->mb_type!=IPCM)
  {
    
    // intra prediction modes for a macroblock 4x4 **********************************************
    read_ipred_modes(img,inp);
    
    // read inter frame vector data *********************************************************
    if(mb64 && mb32count==0 && mb16count==0)
    {
      if (MB64.mb_type != 0)
         readMotionInfoFromNAL32 (img, inp, 2);
    }
    else if(!mb64 && mb32 && (mb16count&3)==0)
    {
      if (MB32.mb_type != 0)
         readMotionInfoFromNAL32 (img, inp, 1);
    }
    else if(mb64==0 && mb32==0 && IS_INTERMV (currMB) && (!IS_P8x8(currMB)))
    {
      readMotionInfoFromNAL (img, inp);
    }
    // read CBP and Coeffs  ***************************************************************

    if(mb64 && mb32count==0 && mb16count==0)
    {
      MB64.cbp=readCBP32(img,inp,1);//with dquant
    }


    if(mb64 && MB64.cbp !=0 && (mb16count&3)==0)
    {      
      MB32.cbp=readCBP32(img,inp,0);//with dquant
    }

    if(!mb64 && mb32 && (mb16count&3)==0)
    {      
      MB32.cbp=readCBP32(img,inp,!delta_qp_sent);//with dquant
    }

    if( (mb64 && MB64.cbp!=0 && MB32.cbp !=0)
      ||(!mb64 && mb32 && MB32.cbp !=0 )
      ||(!mb64 && !mb32) )
    {
        readCBPandCoeffsFromNAL (img,inp, (mb64|mb32), 0);
    }
  }
  else
  {
    //read pcm_alignment_zero_bit and pcm_byte[i] 
    
    // here dP is assigned with the same dP as SE_MBTYPE, because IPCM syntax is in the 
    // same category as MBTYPE
    dP = &(currSlice->partArr[partMap[SE_MBTYPE]]);
    readIPCMcoeffsFromNAL(img,inp,dP);
  }
  
#ifdef ADAPTIVE_FD_SD_CODING
  if (currMB->written_SD_Coding_on_off==0)
  {
    currMB->SD_Coding_on_off=0;
  }
#endif
  return DECODE_MB;
}


void propagate_mb_info(Macroblock* currMB32, int CurrentMbAddr, int CurrentMbAddr_x, int CurrentMbAddr_y, int extend_x, int extend_y,
                       int saved_cluster_qp)      
{
  int i, j, encodeMbAddr;
        //set 16x16 block type 
        for(i=CurrentMbAddr_y; i<=extend_y; i++)
        {
          for(j=CurrentMbAddr_x; j<=extend_x; j++)
          {
            encodeMbAddr = i*img->PicWidthInMbs + j;

            if(img->type==P_SLICE && currMB32->mb_type == 0)
              assert(currMB32->cbp==0);

            if(currMB32->mb_type != P8x8)
            {
              if(currMB32->mb_type == 0 && currMB32->cbp==0)
                img->mb_data[encodeMbAddr].skip_flag = 1;
              else
                img->mb_data[encodeMbAddr].skip_flag = 0;

              if(img->mb_data[CurrentMbAddr].mb_type == 0)
                img->mb_data[encodeMbAddr].mb_type = 0;
              else
              {
                if( CurrentMbAddr != encodeMbAddr)
                {
                  img->mb_data[encodeMbAddr].mb_type = 1;
                }
              }

            }

          }//j        
        }//i


      //qp handling
      for(i=CurrentMbAddr_y; i<=extend_y; i++)
      {
        for(j=CurrentMbAddr_x; j<=extend_x; j++)
        {
          //int b4_x,b4_y, list, k;

          encodeMbAddr = i*img->PicWidthInMbs + j;
          img->mb_data[encodeMbAddr].qp = saved_cluster_qp;
          set_chroma_qp(&img->mb_data[encodeMbAddr]);
          //img->mb_data[encodeMbAddr].delta_quant = saved_delta_qp;
        }
      }
      
}
void propagate_qp_info(int CurrentMbAddr_x, int CurrentMbAddr_y, int extend_x, int extend_y,
                       int saved_cluster_qp, int saved_delta_qp)      
{
  int i, j, encodeMbAddr;
      //qp handling
      for(i=CurrentMbAddr_y; i<=extend_y; i++)
      {
        for(j=CurrentMbAddr_x; j<=extend_x; j++)
        {
          //int b4_x,b4_y, list, k;

          encodeMbAddr = i*img->PicWidthInMbs + j;
          img->mb_data[encodeMbAddr].qp = saved_cluster_qp;
          set_chroma_qp(&img->mb_data[encodeMbAddr]);
          img->mb_data[encodeMbAddr].delta_quant = saved_delta_qp;
        }
      }
      last_dquant = saved_delta_qp;

}
/*!
 ************************************************************************
 * \brief
 *    CurrentMbAddr   :    16x16MB index in the raster scan order starting from 0
 *    CurrentMbAddr_x : 16x16MB x-axis index
 *    CurrentMbAddr_y : 16x16MB y_axis index
 *    extend_x        : the horizonal next 16x16MB x-axis index in the current cluster
 *    extend_y        : the vertical  next 16x16MB y-axis index in the current cluster
 *    num_MB          : number of 16x16 MB in the 32x32 cluster
 *
 ************************************************************************
 */
void find_Mb_cluster(int level, int CurrentMbAddr, int *CurrentMbAddr_x, int *CurrentMbAddr_y, int *extend_x, int *extend_y, int *num_MB, int *nextMbAddr)
{
  int extend = 1<<level;

  *CurrentMbAddr_x = CurrentMbAddr % img->PicWidthInMbs;
  *CurrentMbAddr_y = CurrentMbAddr / img->PicWidthInMbs;

  *extend_x = min(*CurrentMbAddr_x+extend-1, (int) img->PicWidthInMbs-1);
  *extend_y = min(*CurrentMbAddr_y+extend-1, (int) img->FrameHeightInMbs-1);

  *num_MB = (*extend_x - *CurrentMbAddr_x + 1) * (*extend_y - *CurrentMbAddr_y + 1);

  if(*extend_x+1 <= (int) img->PicWidthInMbs-1)
    *nextMbAddr = CurrentMbAddr + extend;
  else if(*extend_y+1 <= (int) img->FrameHeightInMbs-1)
    *nextMbAddr = (*extend_y+1) * img->PicWidthInMbs;
  else
    *nextMbAddr = -1;

}
Boolean read_macroblock_cluster32(struct inp_par *inp, int saved_last_dquant, int *delta_qp_sent)
{
  Boolean end_of_slice = FALSE;
  int read_flag;

      int i, j, mb_count;
      int CurrentMbAddr = img->current_mb_nr, encodeMbAddr;
      int CurrentMbAddr_x, CurrentMbAddr_y, extend_x, extend_y, num_MB, nextMbAddr;
      int saved_delta_qp = 0, saved_cluster_qp=0;
      find_Mb_cluster(1, CurrentMbAddr, &CurrentMbAddr_x, &CurrentMbAddr_y, &extend_x, &extend_y, &num_MB, &nextMbAddr);

      saved_delta_qp = 0;      
      last_dquant = saved_last_dquant;

      if(num_MB == 4)
      {
        mb_count=0;

        for(i=CurrentMbAddr_y; i<=extend_y; i++)
        {
          for(j=CurrentMbAddr_x; j<=extend_x; j++)
          {
            //int b4_x,b4_y, list, k;

            img->current_mb_nr = encodeMbAddr = i*img->PicWidthInMbs + j;
            
            //start_macroblock32 (encodeMbAddr, FALSE, 0);//reset img->current_mb_nr and other global variables
            // Initializes the current macroblock
            start_macroblock(img,inp, img->current_mb_nr);
          }
        }

        memset(mb32_mvd, 0, 2*(BLOCK_MULTIPLE<<MAX_MB_EXT_LEVEL)*(BLOCK_MULTIPLE<<MAX_MB_EXT_LEVEL)*2 * sizeof(int));

        for(i=CurrentMbAddr_y; i<=extend_y; i++)
        {
          for(j=CurrentMbAddr_x; j<=extend_x; j++)
          {
            int b4_x,b4_y, list, k;

            img->current_mb_nr = encodeMbAddr = i*img->PicWidthInMbs + j;
            
            //start_macroblock32 (encodeMbAddr, FALSE, 0);//reset img->current_mb_nr and other global variables
            // Initializes the current macroblock
            start_macroblock(img,inp, img->current_mb_nr);
            
            reset_coeffs();
            
            // Get the syntax elements from the NAL

            if(mb_count !=0 && MB32.mb_type != P8x8)
            {
              int ii;
              img->mb_data[encodeMbAddr].mb_type = MB32.mb_type;
              for(ii=0; ii<4; ii++)
              {
                img->mb_data[encodeMbAddr].b8mode[ii] = MB32.mb_type;
                img->mb_data[encodeMbAddr].b8pdir[ii] = MB32.b8pdir[mb_count];
              }

            }
            read_flag = read_one_macroblock32(img,inp, mb_count);   

            if(mb_count == 0 && MB32.mb_type != P8x8)
            {
              int ii;
              img->mb_data[encodeMbAddr].mb_type = MB32.mb_type;
              for(ii=0; ii<4; ii++)
              {
                img->mb_data[encodeMbAddr].b8mode[ii] = MB32.mb_type;
                img->mb_data[encodeMbAddr].b8pdir[ii] = MB32.b8pdir[mb_count];
              }
              saved_cluster_qp = img->qp;        
              saved_delta_qp = img->mb_data[encodeMbAddr].delta_quant;
              if(MB32.cbp !=0)
                *delta_qp_sent = 1;
            }

            if(MB32.mb_type == P8x8)
            {
              if( (img->mb_data[encodeMbAddr].cbp > 0 || img->mb_data[encodeMbAddr].mb_type==I16MB) && *delta_qp_sent==0 )
              {
                *delta_qp_sent = 1;
                saved_delta_qp = img->mb_data[encodeMbAddr].delta_quant;
                saved_cluster_qp = img->qp;
              }
              //printf("mb %d saved_delta_qp %d\n", img->current_mb_nr, saved_delta_qp);
            }

            if(MB32.mb_type == P8x8 && *delta_qp_sent == 0)
              last_dquant = saved_last_dquant;



#ifdef ADAPTIVE_FILTER  
            if(UseAdaptiveFilterForCurrentFrame() && img->AdaptiveFilterFlag == 2) // separable aif
              decode_one_macroblock_sepAIF(img,inp);
            else
#endif

            decode_one_macroblock(img,inp, 1, mb_count);


            if(MB32.mb_type != P8x8 && MB32.mb_type != 0)
            {
              for(list=0; list<2; list++)            
                for(b4_y=0; b4_y<4; b4_y++)
                  for(b4_x=0; b4_x<4; b4_x++)
                    for(k=0; k<2; k++)
                    {                    
                      img->mb_data[encodeMbAddr].mvd[list][b4_y][b4_x][k] = mb32_mvd[list][(mb_count>>1)*4+b4_y][(mb_count&1)*4+b4_x][k];
                    }
            }

            mb_count++;
           }//j
        }//i
     
        img->num_dec_mb += num_MB-1; //img->num_dec_mb will be increased by 1 in exit_macroblock()
        img->current_mb_nr = nextMbAddr-1;//img->current_mb_nr will be increased by 1 in exit_macroblock()

        
        end_of_slice=(Boolean)exit_macroblock(img,inp,(!img->MbaffFrameFlag||img->current_mb_nr%2)); //exit_macroblock will access eos bit

        if(MB32.mb_type == P8x8 && *delta_qp_sent==0)
          saved_cluster_qp = img->qp;

        propagate_mb_info(&MB32, CurrentMbAddr, CurrentMbAddr_x,CurrentMbAddr_y,extend_x, extend_y, 
           saved_cluster_qp);
        
        last_dquant = saved_delta_qp;

      }//num_MB=4
      else
      {
        mb_count = 0;
        for(i=CurrentMbAddr_y; i<=extend_y; i++)
        {
          for(j=CurrentMbAddr_x; j<=extend_x; j++)
          {
            //int b4_x,b4_y, list, k;

            img->current_mb_nr = encodeMbAddr = i*img->PicWidthInMbs + j;

            // Initializes the current macroblock
            start_macroblock(img,inp, img->current_mb_nr);
            // Get the syntax elements from the NAL
            read_flag = read_one_macroblock(img,inp);  
            
            //qp handling
            if((img->mb_data[encodeMbAddr].cbp > 0 || img->mb_data[encodeMbAddr].mb_type==I16MB) && *delta_qp_sent==0 )
            {
              *delta_qp_sent = 1;
              saved_cluster_qp = img->qp;
              saved_delta_qp = img->mb_data[encodeMbAddr].delta_quant;
            }
            if(*delta_qp_sent == 0)
              last_dquant = saved_last_dquant;

#ifdef ADAPTIVE_FILTER  
            if(UseAdaptiveFilterForCurrentFrame() && img->AdaptiveFilterFlag == 2) // separable aif
              decode_one_macroblock_sepAIF(img,inp);
            else
#endif

            decode_one_macroblock(img,inp, 0, mb_count);
            end_of_slice=(Boolean)exit_macroblock(img,inp,(!img->MbaffFrameFlag||img->current_mb_nr%2));

            mb_count ++;
          }//j
        }//i

        if(*delta_qp_sent==0)
          saved_cluster_qp = img->qp;

      propagate_qp_info(CurrentMbAddr_x,CurrentMbAddr_y,extend_x, extend_y, 
        saved_cluster_qp, saved_delta_qp);

      }//num_MB!=4

      return end_of_slice;
}

Boolean read_macroblock_cluster64(struct inp_par *inp, int *delta_qp_sent)//MB64X64
{
  Boolean end_of_slice = FALSE;
  int read_flag;

  int i, j, ii, jj, mb32count, mb16count;
  int CurrentMbAddr = img->current_mb_nr, encodeMbAddr;
  int CurrentMbAddr_x, CurrentMbAddr_y, extend_x, extend_y, num_MB, nextMbAddr;
  int saved_delta_qp = 0, saved_cluster_qp=0, saved_last_dquant;

  find_Mb_cluster(2, CurrentMbAddr, &CurrentMbAddr_x, &CurrentMbAddr_y, &extend_x, &extend_y, &num_MB, &nextMbAddr);

      saved_delta_qp = 0;
      saved_last_dquant = last_dquant;

      if(num_MB == 16)
      {
        for(i=CurrentMbAddr_y; i<=extend_y; i+=2)
        {
          for(j=CurrentMbAddr_x; j<=extend_x; j+=2)
          {
            for(ii=0; ii<2; ii++)
            {
              for(jj=0; jj<2; jj++)
              {
                //int b4_x,b4_y, list, k;

                img->current_mb_nr = encodeMbAddr = (i+ii)*img->PicWidthInMbs + (j+jj);
                
                //start_macroblock32 (encodeMbAddr, FALSE, 0);//reset img->current_mb_nr and other global variables
                // Initializes the current macroblock
                start_macroblock(img,inp, img->current_mb_nr);
              }
            }
          }
        }

        memset(mb32_mvd, 0, 2*(BLOCK_MULTIPLE<<MAX_MB_EXT_LEVEL)*(BLOCK_MULTIPLE<<MAX_MB_EXT_LEVEL)*2 * sizeof(int));

        mb16count=0;
        mb32count=0;
        for(i=CurrentMbAddr_y; i<=extend_y; i+=2)
        {
          for(j=CurrentMbAddr_x; j<=extend_x; j+=2)
          {
            if(MB64.mb_type == P8x8 && mb32count > 0)
            {
              memset(mb32_mvd, 0, 2*(BLOCK_MULTIPLE<<MAX_MB_EXT_LEVEL)*(BLOCK_MULTIPLE<<MAX_MB_EXT_LEVEL)*2 * sizeof(int));
            }

            for(ii=0; ii<2; ii++)
            {
              for(jj=0; jj<2; jj++)
              {
                int b4_x,b4_y, list, k;

                img->current_mb_nr = encodeMbAddr = (i+ii)*img->PicWidthInMbs + (j+jj);
                
                //start_macroblock32 (encodeMbAddr, FALSE, 0);//reset img->current_mb_nr and other global variables
                // Initializes the current macroblock
                start_macroblock(img,inp, img->current_mb_nr);
                
                reset_coeffs();
                
                // Get the syntax elements from the NAL

                if(MB64.mb_type != P8x8 && mb16count !=0)
                {
                  int k;
                  img->mb_data[encodeMbAddr].mb_type = MB64.mb_type;
                  for(k=0; k<4; k++)
                  {
                    img->mb_data[encodeMbAddr].b8mode[k] = MB64.mb_type;
                    img->mb_data[encodeMbAddr].b8pdir[k] = MB64.b8pdir[mb32count];
                  }
                }
                if(MB64.mb_type == P8x8 && (mb16count&3) !=0 && MB32.mb_type != P8x8)
                {
                  int k;
                  img->mb_data[encodeMbAddr].mb_type = MB32.mb_type;
                  for(k=0; k<4; k++)
                  {
                    img->mb_data[encodeMbAddr].b8mode[k] = MB32.mb_type;
                    img->mb_data[encodeMbAddr].b8pdir[k] = MB32.b8pdir[mb16count&3];
                  }

                }

                read_flag = read_one_macroblock64(img,inp, mb32count, mb16count);   

                if(mb16count == 0 && MB64.mb_type != P8x8)
                {
                  int k;
                  img->mb_data[encodeMbAddr].mb_type = MB64.mb_type;
                  for(k=0; k<4; k++)
                  {
                    img->mb_data[encodeMbAddr].b8mode[k] = MB64.mb_type;
                    img->mb_data[encodeMbAddr].b8pdir[k] = MB64.b8pdir[mb32count];
                  }
                //qp handling
                  saved_cluster_qp = img->qp;        
                  saved_delta_qp = img->mb_data[encodeMbAddr].delta_quant;
                  if(MB64.cbp !=0)
                    *delta_qp_sent = 1;
                }

                if((mb16count&3) == 0 && MB64.mb_type == P8x8 && MB32.mb_type != P8x8)
                {
                  int k;
                  img->mb_data[encodeMbAddr].mb_type = MB32.mb_type;
                  for(k=0; k<4; k++)
                  {
                    img->mb_data[encodeMbAddr].b8mode[k] = MB32.mb_type;
                    img->mb_data[encodeMbAddr].b8pdir[k] = MB32.b8pdir[mb16count&3];
                  }
                  if(MB32.cbp !=0 && *delta_qp_sent==0)
                  {
                    *delta_qp_sent = 1;                                        

                    saved_cluster_qp = img->qp;                          
                    saved_delta_qp = img->mb_data[encodeMbAddr].delta_quant;
                  }

                  if(*delta_qp_sent == 0)
                    last_dquant = saved_last_dquant;
                }



                if(MB64.mb_type == P8x8 && MB32.mb_type == P8x8)
                {
                  if( (img->mb_data[encodeMbAddr].cbp > 0 || img->mb_data[encodeMbAddr].mb_type==I16MB) && *delta_qp_sent==0 )
                  {
                    *delta_qp_sent = 1;
                    saved_delta_qp = img->mb_data[encodeMbAddr].delta_quant;
                    saved_cluster_qp = img->qp;
                  }
                  //printf("mb %d saved_delta_qp %d\n", img->current_mb_nr, saved_delta_qp);

                  if(*delta_qp_sent == 0)
                    last_dquant = saved_last_dquant;
                }



#ifdef ADAPTIVE_FILTER  
                if(UseAdaptiveFilterForCurrentFrame() && img->AdaptiveFilterFlag == 2) // separable aif
                  decode_one_macroblock_sepAIF(img,inp);
                else
#endif

                  if(MB64.mb_type != P8x8)
                    decode_one_macroblock(img,inp, 2, mb16count);
                  else
                    decode_one_macroblock(img,inp, 1, mb16count);


                if( MB64.mb_type == P8x8 && (MB32.mb_type != P8x8 && MB32.mb_type != 0) )
                {
                  int mb_x, mb_y;
                  mb_x = (mb16count&3)&1;
                  mb_y = (mb16count&3)>>1;
                  for(list=0; list<2; list++)            
                    for(b4_y=0; b4_y<4; b4_y++)
                      for(b4_x=0; b4_x<4; b4_x++)
                        for(k=0; k<2; k++)
                        {                    
                          img->mb_data[encodeMbAddr].mvd[list][b4_y][b4_x][k] = mb32_mvd[list][mb_y*4+b4_y][mb_x*4+b4_x][k];
                        }
                }

                if(MB64.mb_type != P8x8 && MB64.mb_type != 0)
                {
                  int mb_x, mb_y;
                  mb_x = (mb32count&1)*2 + ((mb16count&3)&1);
                  mb_y = (mb32count>>1)*2 + ((mb16count&3)>>1);
                  for(list=0; list<2; list++)            
                    for(b4_y=0; b4_y<4; b4_y++)
                      for(b4_x=0; b4_x<4; b4_x++)
                        for(k=0; k<2; k++)
                        {                    
                          img->mb_data[encodeMbAddr].mvd[list][b4_y][b4_x][k] = mb32_mvd[list][mb_y*4+b4_y][mb_x*4+b4_x][k];
                        }
                }


                mb16count++;
                }                 
             }

             if(MB64.mb_type == P8x8)
             {               
               propagate_mb_info(&MB32, (i*img->PicWidthInMbs+j), j,i,j+1, i+1,                
                   saved_cluster_qp);
             }


             mb32count++;
           }//j
        }//i
     
        img->num_dec_mb += num_MB-1; //img->num_dec_mb will be increased by 1 in exit_macroblock()
        img->current_mb_nr = nextMbAddr-1;//img->current_mb_nr will be increased by 1 in exit_macroblock()

        end_of_slice=(Boolean)exit_macroblock(img,inp,(!img->MbaffFrameFlag||img->current_mb_nr%2)); //exit_macroblock will access eos bit

        if(*delta_qp_sent==0)
          saved_cluster_qp = img->qp;

        if(MB64.mb_type != P8x8)
        {
          propagate_mb_info(&MB64, CurrentMbAddr, CurrentMbAddr_x,CurrentMbAddr_y,extend_x, extend_y, 
           saved_cluster_qp);         
        }
        else
        {
          propagate_qp_info(CurrentMbAddr_x,CurrentMbAddr_y,extend_x, extend_y,         
            saved_cluster_qp, saved_delta_qp);
        }
        last_dquant = saved_delta_qp;

      }//num_MB=16
      else
      {
        int prev_qp = img->qp;
        mb32count = 0;
        for(i=CurrentMbAddr_y; i<=extend_y; i+=2)
        {
          for(j=CurrentMbAddr_x; j<=extend_x; j+=2)
          {
            //int b4_x,b4_y, list, k;

            img->current_mb_nr = encodeMbAddr = i*img->PicWidthInMbs + j;

            // Initializes the current macroblock
            start_macroblock(img,inp, img->current_mb_nr);

            end_of_slice = read_macroblock_cluster32(inp, saved_last_dquant, delta_qp_sent);
           
            mb32count ++;
          }//j
        }//i

        img->current_mb_nr = nextMbAddr;

        if(*delta_qp_sent)
        {       
          propagate_qp_info(CurrentMbAddr_x,CurrentMbAddr_y,extend_x, extend_y,         
            img->qp, img->qp-prev_qp);
        }
        else
        {
          last_dquant = 0;
        }

      }//num_MB!=16

      return end_of_slice;
}

void readRunLevel_CABAC16x16(SyntaxElement *se, struct inp_par *inp, struct img_par *img,  DecodingEnvironmentPtr dep_dp);

const byte SNGL_SCAN16x16[256][2] = 
{
  { 0, 0}, { 1, 0}, { 0, 1}, { 0, 2}, { 1, 1}, { 2, 0}, { 3, 0}, { 2, 1}, { 1, 2}, { 0, 3}, { 0, 4}, { 1, 3}, { 2, 2}, { 3, 1}, { 4, 0}, { 5, 0},
  { 4, 1}, { 3, 2}, { 2, 3}, { 1, 4}, { 0, 5}, { 0, 6}, { 1, 5}, { 2, 4}, { 3, 3}, { 4, 2}, { 5, 1}, { 6, 0}, { 7, 0}, { 6, 1}, { 5, 2}, { 4, 3},
  { 3, 4}, { 2, 5}, { 1, 6}, { 0, 7}, { 0, 8}, { 1, 7}, { 2, 6}, { 3, 5}, { 4, 4}, { 5, 3}, { 6, 2}, { 7, 1}, { 8, 0}, { 9, 0}, { 8, 1}, { 7, 2}, 
  { 6, 3}, { 5, 4}, { 4, 5}, { 3, 6}, { 2, 7}, { 1, 8}, { 0, 9}, { 0,10}, { 1, 9}, { 2, 8}, { 3, 7}, { 4, 6}, { 5, 5}, { 6, 4}, { 7, 3}, { 8, 2},
  { 9, 1}, {10, 0}, {11, 0}, {10, 1}, { 9, 2}, { 8, 3}, { 7, 4}, { 6, 5}, { 5, 6}, { 4, 7}, { 3, 8}, { 2, 9}, { 1,10}, { 0,11}, { 0,12}, { 1,11},
  { 2,10}, { 3, 9}, { 4, 8}, { 5, 7}, { 6, 6}, { 7, 5}, { 8, 4}, { 9, 3}, {10, 2}, {11, 1}, {12, 0}, {13, 0}, {12, 1}, {11, 2}, {10, 3}, { 9, 4},
  { 8, 5}, { 7, 6}, { 6, 7}, { 5, 8}, { 4, 9}, { 3,10}, { 2,11}, { 1,12}, { 0,13}, { 0,14}, { 1,13}, { 2,12}, { 3,11}, { 4,10}, { 5, 9}, { 6, 8},
  { 7, 7}, { 8, 6}, { 9, 5}, {10, 4}, {11, 3}, {12, 2}, {13, 1}, {14, 0}, {15, 0}, {14, 1}, {13, 2}, {12, 3}, {11, 4}, {10, 5}, { 9, 6}, { 8, 7},
  { 7, 8}, { 6, 9}, { 5,10}, { 4,11}, { 3,12}, { 2,13}, { 1,14}, { 0,15}, { 1,15}, { 2,14}, { 3,13}, { 4,12}, { 5,11}, { 6,10}, { 7, 9}, { 8, 8},
  { 9, 7}, {10, 6}, {11, 5}, {12, 4}, {13, 3}, {14, 2}, {15, 1}, {15, 2}, {14, 3}, {13, 4}, {12, 5}, {11, 6}, {10, 7}, { 9, 8}, { 8, 9}, { 7,10},
  { 6,11}, { 5,12}, { 4,13}, { 3,14}, { 2,15}, { 3,15}, { 4,14}, { 5,13}, { 6,12}, { 7,11}, { 8,10}, { 9, 9}, {10, 8}, {11, 7}, {12, 6}, {13, 5},
  {14, 4}, {15, 3}, {15, 4}, {14, 5}, {13, 6}, {12, 7}, {11, 8}, {10, 9}, { 9,10}, { 8,11}, { 7,12}, { 6,13}, { 5,14}, { 4,15}, { 5,15}, { 6,14},
  { 7,13}, { 8,12}, { 9,11}, {10,10}, {11, 9}, {12, 8}, {13, 7}, {14, 6}, {15, 5}, {15, 6}, {14, 7}, {13, 8}, {12, 9}, {11,10}, {10,11}, { 9,12},
  { 8,13}, { 7,14}, { 6,15}, { 7,15}, { 8,14}, { 9,13}, {10,12}, {11,11}, {12,10}, {13, 9}, {14, 8}, {15, 7}, {15, 8}, {14, 9}, {13,10}, {12,11},
  {11,12}, {10,13}, { 9,14}, { 8,15}, { 9,15}, {10,14}, {11,13}, {12,12}, {13,11}, {14,10}, {15, 9}, {15,10}, {14,11}, {13,12}, {12,13}, {11,14},
  {10,15}, {11,15}, {12,14}, {13,13}, {14,12}, {15,11}, {15,12}, {14,13}, {13,14}, {12,15}, {13,15}, {14,14}, {15,13}, {15,14}, {14,15}, {15,15}
};

const byte SNGL_SCAN16x8[128][2] = 
{
  { 0, 0}, { 0, 1}, { 1, 0}, { 2, 0}, { 1, 1}, { 0, 2}, { 0, 3}, { 1, 2}, { 2, 1}, { 3, 0}, { 4, 0}, { 3, 1}, { 2, 2}, { 1, 3}, { 0, 4}, { 0, 5}, 
  { 1, 4}, { 2, 3}, { 3, 2}, { 4, 1}, { 5, 0}, { 6, 0}, { 5, 1}, { 4, 2}, { 3, 3}, { 2, 4}, { 1, 5}, { 0, 6}, { 0, 7}, { 1, 6}, { 2, 5}, { 3, 4}, 
  { 4, 3}, { 5, 2}, { 6, 1}, { 7, 0}, { 8, 0}, { 7, 1}, { 6, 2}, { 5, 3}, { 4, 4}, { 3, 5}, { 2, 6}, { 1, 7}, { 2, 7}, { 3, 6}, { 4, 5}, { 5, 4}, 
  { 6, 3}, { 7, 2}, { 8, 1}, { 9, 0}, {10, 0}, { 9, 1}, { 8, 2}, { 7, 3}, { 6, 4}, { 5, 5}, { 4, 6}, { 3, 7}, { 4, 7}, { 5, 6}, { 6, 5}, { 7, 4}, 
  { 8, 3}, { 9, 2}, {10, 1}, {11, 0}, {12, 0}, {11, 1}, {10, 2}, { 9, 3}, { 8, 4}, { 7, 5}, { 6, 6}, { 5, 7}, { 6, 7}, { 7, 6}, { 8, 5}, { 9, 4}, 
  {10, 3}, {11, 2}, {12, 1}, {13, 0}, {14, 0}, {13, 1}, {12, 2}, {11, 3}, {10, 4}, { 9, 5}, { 8, 6}, { 7, 7}, { 8, 7}, { 9, 6}, {10, 5}, {11, 4}, 
  {12, 3}, {13, 2}, {14, 1}, {15, 0}, {15, 1}, {14, 2}, {13, 3}, {12, 4}, {11, 5}, {10, 6}, { 9, 7}, {10, 7}, {11, 6}, {12, 5}, {13, 4}, {14, 3},
  {15, 2}, {15, 3}, {14, 4}, {13, 5}, {12, 6}, {11, 7}, {12, 7}, {13, 6}, {14, 5}, {15, 4}, {15, 5}, {14, 6}, {13, 7}, {14, 7}, {15, 6}, {15, 7}
};
static void readLumaCoeff16x16_CABAC (struct img_par *img,struct inp_par *inp, int mb32)
{
  int mb_nr = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];
  int cbp = currMB->cbp;
  SyntaxElement currSE;
  Slice *currSlice = img->currentSlice;
  DataPartition *dP;
  int *partMap = assignSE2partition[currSlice->dp_mode];
  int coef_ctr;// i0, j0;
  int start_scan;
  int any_coeff;
  int numCoeff = MB_BLOCK_SIZE*MB_BLOCK_SIZE+1;

  int R[6] = {40, 45, 50, 57, 63, 71};
  int level;
  int run, len;

  int qp_per    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)/6;
  int qp_rem    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)%6;

  int i, j, k;

  for(i = 0; i < MB_BLOCK_SIZE; i++)
    memset(img->cof16x16P[i], 0, sizeof(int)*MB_BLOCK_SIZE);

  if (cbp & 0xf)  // are there any coeff in current block at all
  {
    start_scan = 0; // take all coeffs
    coef_ctr = start_scan-1;
    level    = 1;

    for(k=start_scan;(k < numCoeff) && (level != 0);k++)
    {
      //============ read =============
      /*
      * make distinction between INTRA and INTER coded
      * luminance coefficients
      */
      currSE.context      = LUMA_16x16P;
      currSE.type         = (k==0 ? SE_LUM_DC_INTER : SE_LUM_AC_INTER);

#if TRACE
      sprintf(currSE.tracestring, "Luma16x16 sng ");
#endif
      dP = &(currSlice->partArr[partMap[currSE.type]]);
      currSE.reading = mb32 ? readRunLevel_CABAC16x16 : readRunLevel_CABAC;
      dP->readSyntaxElement(&currSE,img,inp,dP);
      level = currSE.value1;
      run   = currSE.value2;
      len   = currSE.len;

      //============ decode =============
      if (level != 0)    /* leave if len=1 */
      {
        any_coeff=1;
        coef_ctr += run+1;

        i=SNGL_SCAN16x16[coef_ctr][0];
        j=SNGL_SCAN16x16[coef_ctr][1];

        currMB->cbp_blk |= 0xffffffff;

        img->cof16x16P[j][i] = level*(R[qp_rem] << qp_per);; // dequantization

        //if(i!=0 || j != 0) img->cof16x16P[j][i]=0;
      }
    }
  }
}


static void readLumaCoeff16x8_CABAC (struct img_par *img,struct inp_par *inp, int blk)
{
  int mb_nr = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];
  int cbp = currMB->cbp;
  SyntaxElement currSE;
  Slice *currSlice = img->currentSlice;
  DataPartition *dP;
  int *partMap = assignSE2partition[currSlice->dp_mode];
  int coef_ctr;// i0, j0;
  int start_scan;
  int any_coeff;
  int numCoeff = MB_BLOCK_SIZE*MB_BLOCK_SIZE/2+1;

  int R[6] = {40, 45, 50, 57, 63, 71};
  int level;
  int run, len;

  int qp_per    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)/6;
  int qp_rem    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)%6;

  int i, j, k;

  for(i = 0; i < MB_BLOCK_SIZE; i++)
    memset(img->cof16x8[blk][i], 0, sizeof(int)*MB_BLOCK_SIZE);

  img->subblock_x = 0;               // horiz. position for coeff_count context
  img->subblock_y = (blk == 0)?0:2;  // vert.  position for coeff_count context

  if (cbp & (0x3 << blk*2))  // are there any coeff in current block at all
  {
    start_scan = 0; // take all coeffs
    coef_ctr = start_scan-1;
    level    = 1;

    for(k=start_scan;(k < numCoeff) && (level != 0);k++)
    {
      //============ read =============
      /*
      * make distinction between INTRA and INTER coded
      * luminance coefficients
      */
      currSE.context      = LUMA_16x8P;
      currSE.type         = (k==0 ? SE_LUM_DC_INTER : SE_LUM_AC_INTER);

#if TRACE
      sprintf(currSE.tracestring, "Luma16x16 sng ");
#endif
      dP = &(currSlice->partArr[partMap[currSE.type]]);
      currSE.reading = readRunLevel_CABAC;

      dP->readSyntaxElement(&currSE,img,inp,dP);
      level = currSE.value1;
      run   = currSE.value2;
      len   = currSE.len;


      //============ decode =============
      if (level != 0)    /* leave if len=1 */
      {
        any_coeff=1;
        coef_ctr += run+1;

        i=SNGL_SCAN16x8[coef_ctr][0];
        j=SNGL_SCAN16x8[coef_ctr][1];

        currMB->cbp_blk |= (0xff<<blk*8);

        img->cof16x8[blk][j][i] = level*(R[qp_rem] << qp_per);; // dequantization
      }
    }
  }
}


static void readLumaCoeff8x16_CABAC (struct img_par *img,struct inp_par *inp, int blk)
{
  int mb_nr = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];
  int cbp = currMB->cbp;
  SyntaxElement currSE;
  Slice *currSlice = img->currentSlice;
  DataPartition *dP;
  int *partMap = assignSE2partition[currSlice->dp_mode];
  int coef_ctr;// i0, j0;
  int start_scan;
  int any_coeff;
  int numCoeff = MB_BLOCK_SIZE*MB_BLOCK_SIZE/2+1;

  int R[6] = {40, 45, 50, 57, 63, 71};
  int level;
  int run, len;

  int qp_per    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)/6;
  int qp_rem    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)%6;

  int i, j, k;

  for(i = 0; i < MB_BLOCK_SIZE; i++)
    memset(img->cof8x16[blk][i], 0, sizeof(int)*MB_BLOCK_SIZE);

  img->subblock_x = (blk == 0)?0:2;     // horiz. position for coeff_count context
  img->subblock_y = 0;                  // vert.  position for coeff_count context

  if (cbp & (0x5<<blk))  // are there any coeff in current block at all
  {
    start_scan = 0; // take all coeffs
    coef_ctr = start_scan-1;
    level    = 1;

    for(k=start_scan;(k < numCoeff) && (level != 0);k++)
    {
      //============ read =============
      /*
      * make distinction between INTRA and INTER coded
      * luminance coefficients
      */
      currSE.context      = LUMA_8x16P;
      currSE.type         = (k==0 ? SE_LUM_DC_INTER : SE_LUM_AC_INTER);

#if TRACE
      sprintf(currSE.tracestring, "Luma16x16 sng ");
#endif
      dP = &(currSlice->partArr[partMap[currSE.type]]);
      currSE.reading = readRunLevel_CABAC;

      dP->readSyntaxElement(&currSE,img,inp,dP);
      level = currSE.value1;
      run   = currSE.value2;
      len   = currSE.len;

      //============ decode =============
      if (level != 0)    /* leave if len=1 */
      {
        any_coeff=1;
        coef_ctr += run+1;

        j=SNGL_SCAN16x8[coef_ctr][0];
        i=SNGL_SCAN16x8[coef_ctr][1];

        currMB->cbp_blk |= (0x3333 << blk*2);

        img->cof8x16[blk][j][i] = level*(R[qp_rem] << qp_per);; // dequantization
      }
    }
  }
}
void readLumaCoeffBigBlock_CABAC (struct img_par *img,struct inp_par *inp, int mb32)
{
  int mb_nr = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];
  int i;

  if(currMB->mb_type == P16x16 || currMB->mb_type == 0 || mb32)
  {
    if(currMB->cbp & 0xf)
      readLumaCoeff16x16_CABAC(img, inp, mb32);
    else
    {
      // set coefficients to all zero
      for(i = 0; i < MB_BLOCK_SIZE; i++)
        memset(img->cof16x16P[i], 0, sizeof(int)*MB_BLOCK_SIZE);
    }
  }
  else if(currMB->mb_type == P16x8)
  {
    readLumaCoeff16x8_CABAC(img, inp, 0);
    readLumaCoeff16x8_CABAC(img, inp, 1);
  }
  else //if(currMB->mb_type == P8x16)
  {
    readLumaCoeff8x16_CABAC(img, inp, 0);
    readLumaCoeff8x16_CABAC(img, inp, 1);
  }
}

// integer DCT16 basis
void calcCoeff16x16I (int coeffI[16][16])
{
  int i, j;
  double a = M_PI/32.;
  double b = 1./(sqrt(2.));
  double c = 1./sqrt(8.);
  double coeff[16][16];
  int factor = 1 << DCT16PREC;

  for(j = 0; j < 16; j++)
    for(i = 0; i < 16; i++)
      coeff[j][i] = cos((2*i+1)*j*a)*c; 

  for(i = 0; i < 16; i++)
    coeff[0][i] *= b; 

  for(j = 0; j < 16; j++)
    for(i = 0; i < 16; i++)
    {
      coeffI[j][i] = (int)(fabs(coeff[j][i])*factor+0.5); 
      if(coeff[j][i] < 0.)
        coeffI[j][i] = - coeffI[j][i];
    }
}

// integer DCT8 basis
void calcCoeff8x8I(int coeffI[16][16])
{
  int i, j;
  double a = M_PI/16.;
  double b = 1./(sqrt(2.));
  double c = 1./sqrt(4.);
  double coeff[16][16];
  int factor = 1 << DCT16PREC;

  for(j = 0; j < 8; j++)
    for(i = 0; i < 8; i++)
      coeff[j][i] = cos((2*i+1)*j*a)*c; 

  for(i = 0; i < 8; i++)
    coeff[0][i] *= b; 

  for(j = 0; j < 8; j++)
    for(i = 0; i < 8; i++)
    {
      coeffI[j][i] = (int)(fabs(coeff[j][i])*factor+0.5); 
      if(coeff[j][i] < 0.)
        coeffI[j][i] = - coeffI[j][i];
    }
}

// inverse 16x16 transform 
void itrans16x16(struct img_par *img)
{
  int x, y, k;
  int shift = DCT16PREC*2+6, shiftby2 = 1<<(shift-1);
  int y2D[MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  int temp2D[MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  int temp;

  // horizontal 
  for(y = 0; y < MB_BLOCK_SIZE; y++)
  {
    for(x = 0; x < MB_BLOCK_SIZE; x++)
    {
      temp2D[y][x] = 0;
      for(k = 0; k < MB_BLOCK_SIZE; k++)
        temp2D[y][x] += img->dct_coeff16I[k][y]*img->cof16x16P[x][k];
    }
  }
  
  // vertical 
  for(y = 0; y < MB_BLOCK_SIZE; y++)
  {
    for(x = 0; x < MB_BLOCK_SIZE; x++)
    {
      y2D[y][x] = 0;
      for(k = 0; k < MB_BLOCK_SIZE; k++)
        y2D[y][x] += img->dct_coeff16I[k][y]*temp2D[x][k];

      temp = max(0,min(img->max_imgpel_value,(y2D[y][x]+((long)img->mpr[x][y] <<shift)+shiftby2)>>shift));
      img->m7[x][y] = temp;

      //printf("%d(%d) ", img->mpr[x][y], y2D[y][x]);

    }
    //printf("\n");
      
  }
}

// inverse 16x8 transform 
void itrans16x8(struct img_par *img, int blk)
{
  int x, y, k;
  int shift = DCT16PREC*2+6, shiftby2 = 1<<(shift-1);
  int y2D[MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  int temp2D[MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  int temp;
  int bw, bh;

  bw = MB_BLOCK_SIZE; bh = MB_BLOCK_SIZE/2;
  // horizontal
  for(y = 0; y < bw; y++)
  {
    for(x = 0; x < bh; x++)
    {
      temp2D[x][y] = 0;
      for(k = 0; k < bw; k++)
        temp2D[x][y] += img->dct_coeff16I[k][y]*img->cof16x8[blk][x][k];
    }
  }
    
  // vertical 
  for(y = 0; y < bh; y++)
  {
    for(x = 0; x < bw; x++)
    {
      y2D[y][x] = 0;
      for(k = 0; k < bh; k++)
        y2D[y][x] += img->dct_coeff8I[k][y]*temp2D[k][x];

      temp = max(0,min(img->max_imgpel_value,(y2D[y][x]+((long)img->mpr[x][y+blk*MB_BLOCK_SIZE/2] <<shift)+shiftby2)>>shift));
      img->m7[x][y+blk*MB_BLOCK_SIZE/2] = temp;
    }
  }
}


// inverse 8x16 transform 
void itrans8x16(struct img_par *img, int blk)
{
  int x, y, k;
  int shift = DCT16PREC*2+6, shiftby2 = 1<<(shift-1);
  int y2D[MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  int temp2D[MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  int temp;
  int bw, bh;

  bh = MB_BLOCK_SIZE; bw = MB_BLOCK_SIZE/2;
  // horizontal
  for(y = 0; y < bw; y++)
  {
    for(x = 0; x < bh; x++)
    {
      temp2D[x][y] = 0;
      for(k = 0; k < bw; k++)
        temp2D[x][y] += img->dct_coeff8I[k][y]*img->cof8x16[blk][x][k];
    }
  }
  
  // vertical 
  for(y = 0; y < bh; y++)
  {
    for(x = 0; x < bw; x++)
    {
      y2D[y][x] = 0;
      for(k = 0; k < bh; k++)
        y2D[y][x] += img->dct_coeff16I[k][y]*temp2D[k][x];

      temp = max(0,min(img->max_imgpel_value,(y2D[y][x]+((long)img->mpr[x+blk*MB_BLOCK_SIZE/2][y] <<shift)+shiftby2)>>shift));
      img->m7[x+blk*MB_BLOCK_SIZE/2][y] = temp;
    }
  }
#ifdef PROFILE_DEC
  if(img->current_mb_nr == 278)
  {
    FILE *fp = fopen("mb278-dec.txt", "a");
    int x, y;

    for(y = 0; y < 16; y++)
    {
      for(x = 0; x < 8; x++)
        fprintf(fp, "%d  ", img->m7[x+blk*MB_BLOCK_SIZE/2][y]);
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n\n");

    for(y = 0; y < 16; y++)
    {
      for(x = 0; x < 8; x++)
        fprintf(fp, "%d  ", y2D[y][x]);
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n\n");
    fclose(fp);
  }
#endif
}

#endif

