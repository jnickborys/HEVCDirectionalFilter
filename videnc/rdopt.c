
/*!
***************************************************************************
* \file rdopt.c
*
* \brief
*    Rate-Distortion optimized mode decision
*
* \author
*    - Heiko Schwarz              <hschwarz@hhi.de>
*    - Valeri George              <george@hhi.de>
*    - Lowell Winger              <lwinger@lsil.com>
*    - Alexis Michael Tourapis    <alexismt@ieee.org>
* \date
*    12. April 2001
**************************************************************************
*/

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <memory.h>
#include <string.h>

#include "global.h"

#include "rdopt_coding_state.h"
#include "memalloc.h"
#include "mb_access.h"
#include "elements.h"
#include "intrarefresh.h"
#include "image.h"
#include "transform8x8.h"
#include "cabac.h"   
#include "vlc.h"
#include "fast_me.h"
#include "ratectl.h"            // head file for rate control
#include "mode_decision.h"
#include "fmo.h"
#define KS_MV
#ifdef ADAPTIVE_FILTER
#include "adaptive_filter.h"
#endif

#ifdef MV_COMPETITION
#include "mv_competition.h"
extern MV_Competition mv_comp;
extern int skip_mode;
#endif

#ifdef RDO_Q
int saved_cbp, saved_I16cbp, saved_I8cbp;
int saved_i16mode;
int saved_cofAC[4][4][2][65];
int saved_I16cofAC[4][4][2][65], saved_I8cofAC[4][4][2][65];
int saved_cofDC[3][2][65];
imgpel saved_enc_picture_imgY[16][16], saved_I16enc_picture_imgY[16][16], saved_I8enc_picture_imgY[16][16];
imgpel saved_img_mpr[16][16], saved_I8img_mpr[16][16];
#endif
#ifdef USE_INTRA_MDDT
extern long quant_stat_rd[16];
extern long quant_stat_rd16x16[256];

int   bestCofAC16[2][16*17];
#endif 


//Rate control

int QP,QP2;
int DELTA_QP,DELTA_QP2;

imgpel pred[16][16];

#define FASTMODE 1
//#define RESET_STATE

extern const int LEVELMVLIMIT[17][6];
extern int   QP2QUANT[40];

extern short OffsetList4x4[15][16];
extern short OffsetList8x8[5][64];
extern const int OffsetBits;

#ifdef MB32X32
int   bestCofAC16x16 [2][257];
int   bestCofAC16x8  [2][2][129];
int   bestCofAC8x16  [2][2][129];
extern int   ****cofAC32[1<<(2*MAX_MB_EXT_LEVEL)], ****cofAC32_RDCost[1<<(2*MAX_MB_EXT_LEVEL)];// ****cofAC8x8=NULL;        // [8x8block][4x4block][level/run][scan_pos]
extern int   ***cofDC32[1<<(2*MAX_MB_EXT_LEVEL)], ***cofDC32_RDCost[1<<(2*MAX_MB_EXT_LEVEL)];                       // [yuv][level/run][scan_pos]
extern int    **cof32AC16x16[1<<(2*MAX_MB_EXT_LEVEL)], **cof32AC16x16_RDCost[1<<(2*MAX_MB_EXT_LEVEL)];//[2][257];
extern int    ***cof32AC16x8[1<<(2*MAX_MB_EXT_LEVEL)], ***cof32AC16x8_RDCost[1<<(2*MAX_MB_EXT_LEVEL)];//[2][2][129];
extern int    ***cof32AC8x16[1<<(2*MAX_MB_EXT_LEVEL)], ***cof32AC8x16_RDCost[1<<(2*MAX_MB_EXT_LEVEL)];//[2][2][129];
imgpel   rec_mbY[16<<MAX_MB_EXT_LEVEL][16<<MAX_MB_EXT_LEVEL], rec_mbU[16<<MAX_MB_EXT_LEVEL][16<<MAX_MB_EXT_LEVEL], rec_mbV[16<<MAX_MB_EXT_LEVEL][16<<MAX_MB_EXT_LEVEL];    // reconstruction values
#else
imgpel   rec_mbY[16][16], rec_mbU[16][16], rec_mbV[16][16];    // reconstruction values
#endif

int lrec_rec[16][16],lrec_rec_U[16][16],lrec_rec_V[16][16]; // store the transf. and quantized coefficients for SP frames

RD_8x8DATA tr4x4, tr8x8;

#ifdef ADAPTIVE_FD_SD_CODING
int   best_SD_or_FD           [2][2];
int   best_SD_or_FD_t8x8;    
int   best_quantizer_indices  [16][16];
int   best_SD_Coding_on_off;
#endif

int   bestInterFAdjust4x4[16][16], bestIntraFAdjust4x4[16][16];
int   bestInterFAdjust8x8[16][16], bestIntraFAdjust8x8[16][16];
int   bestInterFAdjust4x4Cr[2][16][16], bestIntraFAdjust4x4Cr[2][16][16];
int   fadjust8x8[16][16], fadjust4x4[16][16], fadjust4x4Cr[2][16][16], fadjust8x8Cr[2][16][16];    

#ifdef ADAPTIVE_FD_SD_CODING
float adjust_adaptive_f_spatial_domain_8x8;
float adjust_adaptive_f_spatial_domain_4x4;
float best_adjust_adaptive_f_spatial_domain_8x8;
float best_adjust_adaptive_f_spatial_domain_4x4;
#endif

int   ****cofAC=NULL, ****cofAC8x8=NULL;        // [8x8block][4x4block][level/run][scan_pos]
int   ***cofDC=NULL;                       // [yuv][level/run][scan_pos]
int   **cofAC4x4=NULL, ****cofAC4x4intern=NULL; // [level/run][scan_pos]
int   cbp, cbp8x8, cnt_nonz_8x8;
int64 cbp_blk;
int   cbp_blk8x8;
char  frefframe[4][4], brefframe[4][4];
int   b8mode[4], b8pdir[4];
short best8x8mode [4];                // [block]
short best8x8pdir  [MAXMODE][4];       // [mode][block]
short best8x8fwref [MAXMODE][4];       // [mode][block]
short best8x8bwref [MAXMODE][4];       // [mode][block]

#ifdef ADAPTIVE_QUANTIZATION
int   best_mb_iaqms_idx;
#endif
CSptr cs_mb=NULL, cs_b8=NULL, cs_cm=NULL, cs_imb=NULL, cs_ib8=NULL, cs_ib4=NULL, cs_pc=NULL;
#ifdef USE_INTRA_MDDT
CSptr cs_i16 = NULL;
#endif 
int   best_c_imode;
int   best_i16offset;
short best_mode;
short  bi_pred_me;

//mixed transform sizes definitions
int   luma_transform_size_8x8_flag;

short all_mv8x8[2][2][4][4][2];       //[8x8_data/temp_data][LIST][block_x][block_y][MVx/MVy]
short pred_mv8x8[2][2][4][4][2];

int   ****cofAC_8x8ts = NULL;        // [8x8block][4x4block][level/run][scan_pos]

int64    cbp_blk8_8x8ts;
int      cbp8_8x8ts;
int      cost8_8x8ts;
int      cnt_nonz8_8x8ts;

#ifdef RDO_Q
imgpel rec_mbY_intra[4][4][16][16], rec_mbU_intra[4][4][16][16], rec_mbV_intra[4][4][16][16];    // reconstruction values
int   lrec_rec_intra[4][4][16][16],lrec_rec_U_intra[4][4][16][16],lrec_rec_V_intra[4][4][16][16]; // store the transf. and quantized coefficients for SP frames
int   i16offset_intra[4][4];
int   cbp_intra[4][4];
int64 cbp_blk_intra[4][4];
double rdcost_intra[4][4];
int   cofAC_intra[4][4][12][4][2][65];
int   cofDC_intra[4][4][3][2][65];                       // [yuv][level/run][scan_pos]
char   b4_ipredmode_intra[4][4][16];
char   b8_ipredmode8x8_intra[4][4][4][4];
char   b4_ipredmode_mb_intra[4][4][16];
char   b8_ipredmode8x8_mb_intra[4][4][16];
#define IS_INTRA_MODE(mode)  ((mode==I4MB)||(mode==I8MB)||(mode==I16MB)||(mode==IPCM))
#define IS_INTER_SLICE(type) ((type) == B_SLICE || (type) == P_SLICE)
#endif

void StoreMV8x8(int dir);
void RestoreMV8x8(int dir);
// end of mixed transform sizes definitions

//Adaptive Rounding update function
void update_offset_params(int mode, int luma_transform_size_8x8_flag);

// Residue Color Transform
int   cofAC4x4_chroma[2][2][18];
int   rec_resG_8x8[16][16], resTrans_R_8x8[16][16], resTrans_B_8x8[16][16];
int   rec_resG_8x8ts[16][16], resTrans_R_8x8ts[16][16], resTrans_B_8x8ts[16][16];
int   mprRGB_8x8[3][16][16], mprRGB_8x8ts[3][16][16];
char  b4_ipredmode[16], b4_intra_pred_modes[16];

/*!
************************************************************************
* \brief
*    delete structure for RD-optimized mode decision
************************************************************************
*/
void clear_rdopt ()
{
#ifdef MB32X32
  if(input->UseExtMB>0)
  {
    int i;
    for(i=0; i<(1<<(2*MAX_MB_EXT_LEVEL)); i++)
    {
      free_mem_DCcoeff (cofDC32[i]);
      free_mem_ACcoeff (cofAC32[i]);
      free_mem_DCcoeff (cofDC32_RDCost[i]);
      free_mem_ACcoeff (cofAC32_RDCost[i]);
      free_mem2Dint(cof32AC16x16[i]);
      free_mem3Dint(cof32AC16x8[i], 2);
      free_mem3Dint(cof32AC8x16[i], 2);
      free_mem2Dint(cof32AC16x16_RDCost[i]);
      free_mem3Dint(cof32AC16x8_RDCost[i], 2);
      free_mem3Dint(cof32AC8x16_RDCost[i], 2);
    }
  }
#endif
  free_mem_DCcoeff (cofDC);
  free_mem_ACcoeff (cofAC);
  free_mem_ACcoeff (cofAC8x8);
  free_mem_ACcoeff (cofAC4x4intern);
  
  if (input->Transform8x8Mode)
  {
    free_mem_ACcoeff (cofAC_8x8ts);
  }
  
  // structure for saving the coding state
  delete_coding_state (cs_mb);
  delete_coding_state (cs_b8);
  delete_coding_state (cs_cm);
  delete_coding_state (cs_imb);
  delete_coding_state (cs_ib8);
  delete_coding_state (cs_ib4);
  delete_coding_state (cs_pc);
#ifdef USE_INTRA_MDDT
  delete_coding_state (cs_i16);
#endif 
  
#ifdef ADAPTIVE_FD_SD_CODING
  delete_coding_state (cs_spatial_domain_coding);
  delete_coding_state (cs_spatial_domain_coding_rd_opt);
  delete_coding_state (cs_slice_coding);
  delete_coding_state (cs_slice_coding1);
#endif
  
}


/*!
************************************************************************
* \brief
*    create structure for RD-optimized mode decision
************************************************************************
*/
void init_rdopt ()
{
#ifdef MB32X32
  if(input->UseExtMB>0)
  {
    int i;
    for(i=0; i<(1<<(2*MAX_MB_EXT_LEVEL)); i++)
    {
      get_mem_DCcoeff (&cofDC32[i]);
      get_mem_ACcoeff (&cofAC32[i]);
      get_mem_DCcoeff (&cofDC32_RDCost[i]);
      get_mem_ACcoeff (&cofAC32_RDCost[i]);

      get_mem2Dint(&cof32AC16x16[i],2,257);
      get_mem3Dint(&cof32AC16x8[i], 2, 2, 129);
      get_mem3Dint(&cof32AC8x16[i], 2, 2, 129);
      get_mem2Dint(&cof32AC16x16_RDCost[i],2,257);
      get_mem3Dint(&cof32AC16x8_RDCost[i], 2, 2, 129);
      get_mem3Dint(&cof32AC8x16_RDCost[i], 2, 2, 129);
    }
  }
#endif

  rdopt = NULL;
  
  get_mem_DCcoeff (&cofDC);
  get_mem_ACcoeff (&cofAC);
  get_mem_ACcoeff (&cofAC8x8);
  get_mem_ACcoeff (&cofAC4x4intern);
  cofAC4x4 = cofAC4x4intern[0][0];
  
  if (input->Transform8x8Mode)
  {
    get_mem_ACcoeff (&cofAC_8x8ts);
  }
  
  // structure for saving the coding state
  cs_mb  = create_coding_state ();
  cs_b8  = create_coding_state ();
  cs_cm  = create_coding_state ();
  cs_imb = create_coding_state ();
  cs_ib8 = create_coding_state ();
  cs_ib4 = create_coding_state ();
  cs_pc  = create_coding_state ();
#ifdef USE_INTRA_MDDT
  cs_i16  = create_coding_state ();
#endif 
  
#ifdef ADAPTIVE_FD_SD_CODING
  cs_spatial_domain_coding                = create_coding_state ();
  cs_spatial_domain_coding_rd_opt         = create_coding_state ();
  cs_slice_coding                         = create_coding_state ();
  cs_slice_coding1                        = create_coding_state ();
#endif
  
}



/*!
*************************************************************************************
* \brief
*    Updates the pixel map that shows, which reference frames are reliable for
*    each MB-area of the picture.
*
* \note
*    The new values of the pixel_map are taken from the temporary buffer refresh_map
*
*************************************************************************************
*/
void UpdatePixelMap()
{
  int mx,my,y,x,i,j;
  if (img->type==I_SLICE)
  {
    for (y=0; y<img->height; y++)
      for (x=0; x<img->width; x++)
      {
        pixel_map[y][x]=1;
      }
  }
  else
  {
    for (my=0; my<img->height >> 3; my++)
      for (mx=0; mx<img->width >> 3;  mx++)
      {
        j = my*8 + 8;
        i = mx*8 + 8;
        if (refresh_map[my][mx])
        {
          for (y=my*8; y<j; y++)
            for (x=mx*8; x<i; x++)  
              pixel_map[y][x] = 1;
        }
        else
        {
          for (y=my*8; y<j; y++)
            for (x=mx*8; x<i; x++)  
              pixel_map[y][x] = min(pixel_map[y][x] + 1, input->num_ref_frames+1);
        }
      }
  }
}

/*!
*************************************************************************************
* \brief
*    Checks if a given reference frame is reliable for the current
*    macroblock, given the motion vectors that the motion search has
*    returned.
*
* \return
*    If the return value is 1, the reference frame is reliable. If it
*    is 0, then it is not reliable.
*
* \note
*    A specific area in each reference frame is assumed to be unreliable
*    if the same area has been intra-refreshed in a subsequent frame.
*    The information about intra-refreshed areas is kept in the pixel_map.
*
*************************************************************************************
*/
int CheckReliabilityOfRef (int block, int list_idx, int ref, int mode)
{
  int y,x, block_y, block_x, dy, dx, y_pos, x_pos, yy, xx, pres_x, pres_y;
  int maxold_x  = img->width-1;
  int maxold_y  = img->height-1;
  int ref_frame = ref+1;
  
  int by0 = (mode>=4?2*(block >> 1):mode==2?2*block:0);
  int by1 = by0 + (mode>=4||mode==2?2:4);
  int bx0 = (mode>=4?2*(block & 0x01):mode==3?2*block:0);
  int bx1 = bx0 + (mode>=4||mode==3?2:4);
  
  for (block_y=by0; block_y<by1; block_y++)
  {
    for (block_x=bx0; block_x<bx1; block_x++)
    {
      y_pos  = img->all_mv[block_y][block_x][list_idx][ref][mode][1];
      y_pos += (img->block_y + block_y) * BLOCK_SIZE * 4;
      x_pos  = img->all_mv[block_y][block_x][list_idx][ref][mode][0];
      x_pos += (img->block_x + block_x) * BLOCK_SIZE * 4;
      
      /* Here we specify which pixels of the reference frame influence
      the reference values and check their reliability. This is
      based on the function Get_Reference_Pixel */
      
      dy = y_pos & 3;
      dx = x_pos & 3;
      
      y_pos = (y_pos-dy) >> 2;
      x_pos = (x_pos-dx) >> 2;
      
      if (dy==0 && dx==0) //full-pel
      {
        for (y=y_pos ; y < y_pos + BLOCK_SIZE ; y++)
          for (x=x_pos ; x < x_pos + BLOCK_SIZE ; x++)
            if (pixel_map[max(0,min(maxold_y,y))][max(0,min(maxold_x,x))] < ref_frame)
              return 0;
      }
      else  /* other positions */
      {
        if (dy == 0)
        {
          for (y = y_pos ; y < y_pos + BLOCK_SIZE ; y++)
          {
            pres_y = max(0,min(maxold_y,y));
            for (x = x_pos ; x < x_pos + BLOCK_SIZE ; x++)
            {
              for(xx = -2 ; xx < 4 ; xx++)
              {
                pres_x = max(0, min( maxold_x, x + xx));
                if (pixel_map[pres_y][pres_x] < ref_frame)
                  return 0;
              }
            }
          }
        }        
        else if (dx == 0)
        {
          for (y = y_pos ; y < y_pos + BLOCK_SIZE ; y++)
            for (x=x_pos ; x < x_pos + BLOCK_SIZE ; x++)
            {
              pres_x = max(0,min(maxold_x,x));
              for(yy=-2;yy<4;yy++) 
              {
                pres_y = max(0,min(maxold_y, yy + y));
                if (pixel_map[pres_y][pres_x] < ref_frame)
                  return 0;
              }
            }
        }
        else if (dx == 2)
        {
          for (y = y_pos ; y < y_pos + BLOCK_SIZE ; y++)
            for (x = x_pos ; x < x_pos + BLOCK_SIZE ; x++)
            {
              for(yy=-2;yy<4;yy++) 
              {
                pres_y = max(0,min(maxold_y, yy + y));
                for(xx=-2;xx<4;xx++) 
                {
                  pres_x = max(0,min(maxold_x, xx + x));
                  if (pixel_map[pres_y][pres_x] < ref_frame)
                    return 0;
                }
              }
            }
        }
        else if (dy == 2)
        {
          for (y = y_pos ; y < y_pos + BLOCK_SIZE ; y++)
            for (x = x_pos ; x < x_pos + BLOCK_SIZE ; x++)
            {
              for(xx=-2;xx<4;xx++) 
              {
                pres_x = max(0,min(maxold_x, xx + x));
                for(yy=-2;yy<4;yy++) 
                {
                  pres_y = max(0,min(maxold_y, yy + y));
                  if (pixel_map[pres_y][pres_x] < ref_frame)
                    return 0;
                }
              }
            }
        }
        else
        {
          for (y = y_pos ; y < y_pos + BLOCK_SIZE ; y++)
          {
            for (x = x_pos ; x < x_pos + BLOCK_SIZE ; x++)
            {
              pres_y = dy == 1 ? y : y + 1;
              pres_y = max(0,min(maxold_y,pres_y));
              
              for(xx=-2;xx<4;xx++) 
              {
                pres_x = max(0,min(maxold_x,xx + x));
                if (pixel_map[pres_y][pres_x] < ref_frame)
                  return 0;
              }
              
              pres_x = dx == 1 ? x : x + 1;
              pres_x = max(0,min(maxold_x,pres_x));
              
              for(yy=-2;yy<4;yy++) 
              {
                pres_y = max(0,min(maxold_y, yy + y));
                if (pixel_map[pres_y][pres_x] < ref_frame)
                  return 0;
              }
            }
          }
        }        
      }
    }
  }
  return 1;
}



/*!
*************************************************************************************
* \brief
*    R-D Cost for an 4x4 Intra block
*************************************************************************************
*/
double RDCost_for_4x4IntraBlocks (int*    nonzero,
                                  int     b8,
                                  int     b4,
                                  int     ipmode,
                                  double  lambda,
                                  double  min_rdcost,
                                  int mostProbableMode)
{
  double  rdcost;
  int     dummy, x, y, rate;
  int64   distortion  = 0;
  int     block_x     = 8*(b8 & 0x01)+4*(b4 & 0x01);
  int     block_y     = 8*(b8 >> 1)+4*(b4 >> 1);
  int     pic_pix_x   = img->pix_x+block_x;
  int     pic_pix_y   = img->pix_y+block_y;
  int     pic_opix_y  = img->opix_y+block_y;
  imgpel  **imgY      = enc_picture->imgY;
  
  Slice          *currSlice    =  img->currentSlice;
  Macroblock     *currMB       = &img->mb_data[img->current_mb_nr];
  SyntaxElement  *currSE       = &img->MB_SyntaxElements[currMB->currSEnr];
  const int      *partMap      = assignSE2partition[input->partition_mode];
  DataPartition  *dataPart;
  
  //===== perform DCT, Q, IQ, IDCT, Reconstruction =====
  dummy = 0;
  
  if(img->type!=SP_SLICE)
#ifdef USE_INTRA_MDDT
    if(input->UseIntraMDDT)
    {
      *nonzero = (ipmode == 2) ? 
        dct_luma (block_x, block_y, &dummy, 1, ipmode): 
      klt_luma_sep (block_x, block_y, &dummy, ipmode);
    }
    else 
      *nonzero = dct_luma (block_x, block_y, &dummy, 1, ipmode);
#else
    *nonzero = dct_luma (block_x, block_y, &dummy, 1);
#endif
    else if(!si_frame_indicator && !sp2_frame_indicator)
    {
      *nonzero = dct_luma_sp(block_x, block_y, &dummy);
    }
    else
    {   
      *nonzero = dct_luma_sp2(block_x, block_y, &dummy);
    }
    
    //===== get distortion (SSD) of 4x4 block =====
    if(!img->residue_transform_flag)
    {
      for (y=0; y<4; y++)
      {
        for (x=pic_pix_x; x<pic_pix_x+4; x++)
        {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          distortion += SQR_DEPTH(imgY_org[pic_opix_y+y][x], imgY[pic_pix_y+y][x], input->BitDepthLuma, img->BitDepthIncrease);
#else
          distortion += img->quad [imgY_org[pic_opix_y+y][x] - imgY[pic_pix_y+y][x]];
#endif
        }
      }
    }
    
    //===== RATE for INTRA PREDICTION MODE  (SYMBOL MODE MUST BE SET TO UVLC) =====
    currSE->value1 = (mostProbableMode == ipmode) ? -1 : ipmode < mostProbableMode ? ipmode : ipmode-1;
    
    //--- set position and type ---
    currSE->context = 4*b8 + b4;
    currSE->type    = SE_INTRAPREDMODE;
    
    //--- choose data partition ---
    dataPart = &(currSlice->partArr[partMap[SE_INTRAPREDMODE]]);
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
      rate  += writeCoeff4x4_CAVLC (LUMA, b8, b4, 0);
    }
    else
    {
      rate  += writeLumaCoeff4x4_CABAC (b8, b4, 1);
    }
    //reset_coding_state (cs_cm);
    rdcost = (double)distortion + lambda*(double)rate;
    
    if(img->residue_transform_flag)
      return (double)rate;
    else
      return rdcost;
                                  }
                                  
                                  
                                  // Residue Color Transform
int RDCost_for_4x4Blocks_Chroma (int     b8,
                                 int     b4,
                                 int  chroma)
{
  int     rate=0;
  
  Slice          *currSlice    =  img->currentSlice;
  Macroblock     *currMB       = &img->mb_data[img->current_mb_nr];
  SyntaxElement  *currSE       = &img->MB_SyntaxElements[currMB->currSEnr];
  const int      *partMap      = assignSE2partition[input->partition_mode];
  int uv;
  
  //===== perform DCT, Q, IQ, IDCT, Reconstruction =====
  if(b8 > 7) 
    uv = 1;
  else 
    uv = 0;
  
  cbp_chroma_block_temp[uv][2*((b8-4*(uv+1)) & 0x01)+(b4 & 0x01)][2*((b8-4*(uv+1)) >> 1)+(b4 >> 1)] = dct_chroma4x4 (chroma, b8, b4);
  
  //===== RATE for LUMINANCE COEFFICIENTS =====
  if (input->symbol_mode == UVLC)
  {
    rate  = writeCoeff4x4_CAVLC (CHROMA_AC, b8, b4, ((2*(b8 & 0x01)+ (b4 & 0x01))<<4) | (2*(b8 >> 1)+(b4 >> 1)));
  }
  else
  {
    int * ACLevel, * ACRun;
    int level, run, k;
    DataPartition*  dataPart;
    int*            bitCount  = currMB->bitcounter;
    ACLevel = img->cofAC[b8][b4][0];
    ACRun   = img->cofAC[b8][b4][1];
    
    level=1;
    
    img->subblock_y = b4 >> 1;
    img->subblock_x = b4 & 0x01;
    
    for (k=0; k < 17 && level != 0; k++)
    {
      level = currSE->value1 = ACLevel[k]; // level
      run   = currSE->value2 = ACRun  [k]; // run
      
      if (input->symbol_mode == UVLC)   
        currSE->mapping = levrun_linfo_inter;
      else                              
        currSE->writing = writeRunLevel_CABAC;
      
      currSE->context     = CHROMA_AC;
      currSE->type        = SE_CHR_AC_INTRA;
      
      img->is_intra_block =  IS_INTRA(currMB);
      img->is_v_block     = uv;
      
      // choose the appropriate data partition
      dataPart = &(currSlice->partArr[partMap[currSE->type]]);
      dataPart->writeSyntaxElement (currSE, dataPart);
      bitCount[BITS_COEFF_UV_MB] += currSE->len;
      rate                       += currSE->len;
      
      // proceed to next SE
      currSE++;
      currMB->currSEnr++;
    }
  }
  reset_coding_state (cs_cm);
  
  return rate;
}


/*!
*************************************************************************************
* \brief
*    Mode Decision for an 4x4 Intra block
*************************************************************************************
*/
int Mode_Decision_for_4x4IntraBlocks (int  b8,  int  b4,  double  lambda,  int*  min_cost)
{
  int     ipmode, best_ipmode = 0, i, j, k, x, y, cost, dummy;
  int     c_nz, nonzero = 0, diff[16];
  imgpel  rec4x4[4][4];
  double  rdcost;
  int     block_x     = 8*(b8 & 0x01)+4*(b4 & 0x01);
  int     block_y     = 8*(b8 >> 1)+4*(b4 >> 1);
  int     block_x4    = block_x >> 2;
  int     block_y4    = block_y >> 2;
  int     pic_pix_x   = img->pix_x+block_x;
  int     pic_pix_y   = img->pix_y+block_y;
  int     pic_opix_x   = img->opix_x+block_x;
  int     pic_opix_y   = img->opix_y+block_y;
  int     pic_block_x = pic_pix_x >> 2;
  int     pic_block_y = pic_pix_y >> 2;
  double  min_rdcost  = 1e30;
  
  int left_available, up_available, all_available;
  
  char   upMode;
  char   leftMode;
  int     mostProbableMode;
  
  PixelPos left_block;
  PixelPos top_block;
  
  int     lrec4x4[4][4]; 
  int     lrec4x4_c[2][4][4];
  // Residue Color Transform
  int residue_R, residue_G, residue_B;
  int rate, temp;
  int64 distortion;
  int c_ipmode = img->mb_data[img->current_mb_nr].c_ipred_mode;
  int fixedcost = (int) floor(4 * lambda );
  imgpel rec4x4_c[2][4][4];
#ifdef USE_INTRA_MDDT
  long quant_stat_best[16];
#endif
#ifdef BEST_NZ_COEFF
  int best_nz_coeff = 0;
  int best_coded_block_flag = 0;
#ifdef ADAPTIVE_FD_SD_CODING
  int best_SD_or_FD_flag = 0;
#endif
  int bit_pos = 1 + ((((b8>>1)<<1)+(b4>>1))<<2) + (((b8&1)<<1)+(b4&1));
  static int64 cbp_bits;
#ifdef ADAPTIVE_FD_SD_CODING
  static int64 FD_or_SD_bits;
#endif
  
#ifdef ADAPTIVE_FD_SD_CODING
  if (b8==0 && b4==0)
  {
#else
    if (b8==0 && b4==0)
#endif
      cbp_bits = 0;
#ifdef ADAPTIVE_FD_SD_CODING
    FD_or_SD_bits = 0;
  }
#endif
  
#endif
  
  getLuma4x4Neighbour(img->current_mb_nr, block_x4, block_y4, -1,  0, &left_block);
  getLuma4x4Neighbour(img->current_mb_nr, block_x4, block_y4,  0, -1, &top_block);
  
  // constrained intra pred
  if (input->UseConstrainedIntraPred)
  {
    left_block.available = left_block.available ? img->intra_block[left_block.mb_addr] : 0;
    top_block.available  = top_block.available  ? img->intra_block[top_block.mb_addr]  : 0;
  }
  
  upMode            =  top_block.available ? img->ipredmode[top_block.pos_y ][top_block.pos_x ] : -1;
  leftMode          = left_block.available ? img->ipredmode[left_block.pos_y][left_block.pos_x] : -1;
  
  mostProbableMode  = (upMode < 0 || leftMode < 0) ? DC_PRED : upMode < leftMode ? upMode : leftMode;
  
  *min_cost = INT_MAX;
  
  //===== INTRA PREDICTION FOR 4x4 BLOCK =====
  intrapred_luma (pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);
  
  //===== LOOP OVER ALL 4x4 INTRA PREDICTION MODES =====
  for (ipmode=0; ipmode<NO_INTRA_PMODE; ipmode++)
  {
    int available_mode =  (ipmode==DC_PRED) ||
      ((ipmode==VERT_PRED||ipmode==VERT_LEFT_PRED||ipmode==DIAG_DOWN_LEFT_PRED) && up_available ) ||
      ((ipmode==HOR_PRED||ipmode==HOR_UP_PRED) && left_available ) ||(all_available);
    
    if (input->IntraDisableInterOnly==0 || img->type != I_SLICE)
    {
      if (input->Intra4x4ParDisable && (ipmode==VERT_PRED||ipmode==HOR_PRED))
        continue;      
      
      if (input->Intra4x4DiagDisable && (ipmode==DIAG_DOWN_LEFT_PRED||ipmode==DIAG_DOWN_RIGHT_PRED))
        continue;      
      
      if (input->Intra4x4DirDisable && ipmode>=VERT_RIGHT_PRED)
        continue;
    }
    
    if( available_mode)
    {
      if (!input->rdopt)
      {
        for (k=j=0; j<4; j++)
        {
          int jj = pic_opix_y+j;
          for (i=0; i<4; i++, k++)
          {
            diff[k] = imgY_org[jj][pic_opix_x+i] - img->mprr[ipmode][j][i];
          }
        }
        //cost  = (ipmode == mostProbableMode) ? 0 : (int)floor(4 * lambda );
        cost  = (ipmode == mostProbableMode) ? 0 : fixedcost;
        cost += SATD (diff, input->hadamard);
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
          for (j=0; j<4; j++)
          {
            memcpy(&img->mpr[block_y+j][block_x], img->mprr[ipmode][j], BLOCK_SIZE * sizeof(imgpel));
            for (i=0; i<4; i++)
            {
              img->m7[j][i] = (int) (imgY_org[pic_opix_y+j][pic_opix_x+i] - img->mprr[ipmode][j][i]);
            }
          }
          
          //===== store the coding state =====
          //store_coding_state (cs_cm);
          // get and check rate-distortion cost
#ifdef BEST_NZ_COEFF
          img->mb_data[img->current_mb_nr].cbp_bits = cbp_bits;
#ifdef ADAPTIVE_FD_SD_CODING
          img->mb_data[img->current_mb_nr].FD_or_SD_bits = FD_or_SD_bits;
#endif
#endif
          if ((rdcost = RDCost_for_4x4IntraBlocks (&c_nz, b8, b4, ipmode, lambda, min_rdcost, mostProbableMode)) < min_rdcost)
          {
#ifdef USE_INTRA_MDDT
            if(input->UseIntraMDDT)
              memcpy(quant_stat_best, quant_stat_rd, 16*sizeof(long)); 
#endif
            //--- set coefficients ---
            memcpy(cofAC4x4[0],img->cofAC[b8][b4][0], 18 * sizeof(int));
            memcpy(cofAC4x4[1],img->cofAC[b8][b4][1], 18 * sizeof(int));
            
            //--- set reconstruction ---
            for (y=0; y<4; y++)
            {
              memcpy(rec4x4[y],&enc_picture->imgY[pic_pix_y+y][pic_pix_x], BLOCK_SIZE * sizeof(imgpel));
              if(img->type==SP_SLICE &&(!si_frame_indicator && !sp2_frame_indicator))
                memcpy(lrec4x4[y],&lrec[pic_pix_y+y][pic_pix_x], BLOCK_SIZE * sizeof(int));// stores the mode coefficients
            } 
            //--- flag if dct-coefficients must be coded ---
            nonzero = c_nz;
            
            //--- set best mode update minimum cost ---
            min_rdcost    = rdcost;
            best_ipmode   = ipmode;
#ifdef BEST_NZ_COEFF
            best_nz_coeff = img->nz_coeff [img->current_mb_nr][block_x4][block_y4];
            best_coded_block_flag = (int)((img->mb_data[img->current_mb_nr].cbp_bits>>bit_pos)&(int64)(1));
#ifdef ADAPTIVE_FD_SD_CODING
            best_SD_or_FD_flag    = (int)((img->mb_data[img->current_mb_nr].FD_or_SD_bits>>bit_pos)&(int64)(1));
#endif
#endif            
            //store_coding_state (cs_ib4);
            if (img->AdaptiveRounding)
            {
              for (j=0; j<4; j++)
                memcpy(&fadjust4x4[block_y+j][block_x],&img->fadjust4x4[1][block_y+j][block_x], BLOCK_SIZE * sizeof(int));
            }
          }
          
#ifndef RESET_STATE
          reset_coding_state (cs_cm);
#endif
        }
        else 
        {
          for (j=0; j<4; j++)
          {
            for (i=0; i<4; i++)
            {
              residue_B = imgUV_org[0][pic_opix_y+j][pic_opix_x+i] - img->mprr_c[0][c_ipmode][block_y+j][block_x+i];
              residue_G = imgY_org[pic_opix_y+j][pic_opix_x+i] - img->mprr[ipmode][j][i];
              residue_R = imgUV_org[1][pic_opix_y+j][pic_opix_x+i] - img->mprr_c[1][c_ipmode][block_y+j][block_x+i];
              
              /* Foward Residue Transform */
              resTrans_R[j][i] = residue_R-residue_B;
              temp = residue_B+(resTrans_R[j][i]>>1);
              resTrans_B[j][i] = residue_G-temp;
              resTrans_G[j][i] = temp+(resTrans_B[j][i]>>1);
            }
          }
          
          for (j=0; j<4; j++)
          {
            for (i=0; i<4; i++)
            {
              img->m7[j][i]  = resTrans_G[j][i];
            }
          }
          
          store_coding_state (cs_cm);
          // yuwen 2005.11.17 => should be an error here, requiring process for img->mb_data[img->current_mb_nr].cbp_bits
          rate = (int) RDCost_for_4x4IntraBlocks (&c_nz, b8, b4, ipmode, lambda, min_rdcost, mostProbableMode);
          reset_coding_state (cs_cm);
          
          for (j=0; j<4; j++)
          {
            for (i=0; i<4; i++)
            {
              rec_resG[j][i] = img->m7[j][i];
              img->m7[j][i]  = resTrans_B[j][i];
            }
          }
          //store_coding_state (cs_cm);
          // yuwen 2005.11.17 => should be an error here, requiring process for img->mb_data[img->current_mb_nr].cbp_bits
          rate += RDCost_for_4x4Blocks_Chroma (b8+4, b4, 0);
          for (j=0; j<4; j++)
          {
            for (i=0; i<4; i++)
            {
              rec_resB[j][i] = img->m7[j][i];
              img->m7[j][i]  = resTrans_R[j][i];
            }
          }
          // yuwen 2005.11.17 => should be an error here, requiring process for img->mb_data[img->current_mb_nr].cbp_bits
          rate += RDCost_for_4x4Blocks_Chroma (b8+8, b4, 1);
          
          reset_coding_state (cs_cm);
          for (j=0; j<4; j++)
          {
            for (i=0; i<4; i++)
            {
              rec_resR[j][i] = img->m7[j][i];
            }
          }
          
          for (j=0; j<4; j++)
          {
            for (i=0; i<4; i++)
            {
              /* Inverse Residue Transform */
              temp      = rec_resG[j][i]-(rec_resB[j][i]>>1);
              residue_G = rec_resB[j][i]+temp;
              residue_B = temp - (rec_resR[j][i]>>1);
              residue_R = residue_B+rec_resR[j][i];
              enc_picture->imgUV[0][pic_pix_y+j][pic_pix_x+i] = min(img->max_imgpel_value_uv,max(0,residue_B+(int)img->mprr_c[0][c_ipmode][block_y+j][block_x+i]));
              enc_picture->imgY[pic_pix_y+j][pic_pix_x+i]     = min(img->max_imgpel_value,max(0,residue_G+(int)img->mprr[ipmode][j][i]));
              enc_picture->imgUV[1][pic_pix_y+j][pic_pix_x+i] = min(img->max_imgpel_value_uv,max(0,residue_R+(int)img->mprr_c[1][c_ipmode][block_y+j][block_x+i]));
            }
          } 
          //===== get distortion (SSD) of 4x4 block =====
          distortion = 0;
          for (y=0; y<4; y++)
          {
            for (x=pic_pix_x; x<pic_pix_x+4; x++)
            {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
              distortion += SQR_DEPTH(imgY_org    [pic_pix_y+y][x], enc_picture->imgY    [pic_pix_y+y][x], input->BitDepthLuma,   img->BitDepthIncrease);
              distortion += SQR_DEPTH(imgUV_org[0][pic_pix_y+y][x], enc_picture->imgUV[0][pic_pix_y+y][x], input->BitDepthChroma, img->BitDepthIncreaseChroma);
              distortion += SQR_DEPTH(imgUV_org[1][pic_pix_y+y][x], enc_picture->imgUV[1][pic_pix_y+y][x], input->BitDepthChroma, img->BitDepthIncreaseChroma);
#else
              distortion += img->quad[imgY_org    [pic_pix_y+y][x] - enc_picture->imgY    [pic_pix_y+y][x]];
              distortion += img->quad[imgUV_org[0][pic_pix_y+y][x] - enc_picture->imgUV[0][pic_pix_y+y][x]];
              distortion += img->quad[imgUV_org[1][pic_pix_y+y][x] - enc_picture->imgUV[1][pic_pix_y+y][x]];
#endif
            }
          }
          rdcost = (double)distortion + lambda*(double)rate;
          
          if (rdcost < min_rdcost)
          {
            //--- set coefficients ---
            for (j=0; j<2; j++)
            {
              for (i=0; i<18;i++)  
                cofAC4x4[j][i]=img->cofAC[b8][b4][j][i];
              for (i=0; i<18;i++)  
                cofAC4x4_chroma[0][j][i]=img->cofAC[b8+4][b4][j][i];
              for (i=0; i<18;i++)
                cofAC4x4_chroma[1][j][i]=img->cofAC[b8+8][b4][j][i];
            }
            
            for (i=0; i<2; i++)
            { //uv
              dc_level        [i][2*(b8 & 0x01)+(b4 & 0x01)][2*(b8 >> 1)+(b4 >> 1)] = dc_level_temp        [i][2*(b8 & 0x01)+(b4 & 0x01)][2*(b8 >> 1)+(b4 >> 1)];
              cbp_chroma_block[i][2*(b8 & 0x01)+(b4 & 0x01)][2*(b8 >> 1)+(b4 >> 1)] = cbp_chroma_block_temp[i][2*(b8 & 0x01)+(b4 & 0x01)][2*(b8 >> 1)+(b4 >> 1)];
              //--- set reconstruction ---
              for (y=0; y<BLOCK_SIZE; y++)
              {
                memcpy(rec4x4_c[i][y],&enc_picture->imgUV[i][pic_pix_y+y][pic_pix_x], BLOCK_SIZE * sizeof(imgpel)); 
                if(img->type==SP_SLICE &&(!si_frame_indicator && !sp2_frame_indicator))
                  memcpy(lrec4x4[y],&lrec[pic_pix_y+y][pic_pix_x], BLOCK_SIZE * sizeof(int));//stores the coefficients for the mode
              }
            }
            
            //--- set reconstruction ---
            for (y=0; y<BLOCK_SIZE; y++)
            {
              memcpy(rec4x4[y],&enc_picture->imgY[pic_pix_y+y][pic_pix_x], BLOCK_SIZE * sizeof(imgpel));
              if(img->type==SP_SLICE &&(!si_frame_indicator && !sp2_frame_indicator))
                memcpy(lrec4x4[y],&lrec[pic_pix_y+y][pic_pix_x], BLOCK_SIZE * sizeof(int));
            }
            
            //--- flag if dct-coefficients must be coded ---
            nonzero = c_nz;
            
            //--- set best mode update minimum cost ---
            min_rdcost  = rdcost;
            best_ipmode = ipmode;
#ifdef BEST_NZ_COEFF
            best_nz_coeff = img->nz_coeff [img->current_mb_nr][block_x4][block_y4];
            // yuwen 2005.11.17 => should be an error here, requiring process for img->mb_data[img->current_mb_nr].cbp_bits
#endif
          }
        }
      }
    }
  }
  
#ifdef BEST_NZ_COEFF
  img->nz_coeff [img->current_mb_nr][block_x4][block_y4] = best_nz_coeff;
  cbp_bits &= (~(int64)(1<<bit_pos));
  cbp_bits |= (int64)(best_coded_block_flag<<bit_pos);
#ifdef ADAPTIVE_FD_SD_CODING
  FD_or_SD_bits &= (~(int64)(1<<bit_pos));
  FD_or_SD_bits |= (int64)(best_SD_or_FD_flag<<bit_pos);
#endif
  
#endif
  
#ifdef USE_INTRA_MDDT
  if(input->UseIntraMDDT)
    for (j=0; j<16; j++)
    {
      img->quant_stat[j]+=(quant_stat_best[j]>>5);
    }
#endif
    
    //===== set intra mode prediction =====
    img->ipredmode[pic_block_y][pic_block_x] = best_ipmode;
    img->mb_data[img->current_mb_nr].intra_pred_modes[4*b8+b4] = mostProbableMode == best_ipmode ? -1 : best_ipmode < mostProbableMode ? best_ipmode : best_ipmode-1;
    
    if (!input->rdopt)
    {
      // Residue Color Transform
      if(!img->residue_transform_flag)
      {
        // get prediction and prediction error
        for (j=0; j<4; j++)
        {
          int jj = pic_opix_y+j;
          for (i=0; i<4; i++)
          {
            img->mpr[block_y+j][block_x+i]  = img->mprr[best_ipmode][j][i];
            img->m7[j][i]                   = imgY_org[jj][pic_opix_x+i] - img->mprr[best_ipmode][j][i];
          }
        }
#ifdef USE_INTRA_MDDT
        if(input->UseIntraMDDT)
        {
          nonzero = (best_ipmode == 2) ? 
            dct_luma     (block_x, block_y, &dummy, 1, best_ipmode) : 
          klt_luma_sep (block_x, block_y, &dummy, best_ipmode);
        }
        else 
          nonzero = dct_luma     (block_x, block_y, &dummy, 1, best_ipmode);
#else
        nonzero = dct_luma (block_x, block_y, &dummy, 1);
#endif
      } 
      else 
      {
        int y_pos = 2*(b8 & 0x01)+(b4 & 0x01);
        int x_pos = 2*(b8 >> 1)+(b4 >> 1);
        for (j=0; j<4; j++)
        {
          for (i=0; i<4; i++)
          {
            residue_B = imgUV_org[0][pic_opix_y+j][pic_opix_x+i] - img->mprr_c[0][c_ipmode][block_y+j][block_x+i];
            residue_G = imgY_org[pic_opix_y+j][pic_opix_x+i] - img->mprr[best_ipmode][j][i];
            residue_R = imgUV_org[1][pic_opix_y+j][pic_opix_x+i] - img->mprr_c[1][c_ipmode][block_y+j][block_x+i];
            
            /* Forward Residue Transform */
            resTrans_R[j][i] = residue_R-residue_B;
            temp = residue_B+(resTrans_R[j][i]>>1);
            resTrans_B[j][i] = residue_G-temp;
            resTrans_G[j][i] = temp+(resTrans_B[j][i]>>1);
          }
        }
        
        for (j=0; j<4; j++)
        {
          for (i=0; i<4; i++)
          {
            img->m7[j][i]  = resTrans_G[j][i];
          }
        }
#ifdef USE_INTRA_MDDT
        nonzero = dct_luma (block_x, block_y, &dummy, 1, c_ipmode);
#else
        nonzero = dct_luma (block_x, block_y, &dummy, 1);
#endif
        for (j=0; j<4; j++)
        {
          for (i=0; i<4; i++)
          {
            rec_resG[j][i] = img->m7[j][i];
            img->m7[j][i]  = resTrans_B[j][i];
          }
        }
        cbp_chroma_block[0][y_pos][x_pos] = dct_chroma4x4 (0, b8+4, b4);
        dc_level        [0][y_pos][x_pos] = dc_level_temp[0][y_pos][x_pos];
        for (j=0; j<4; j++)
        {
          for (i=0; i<4; i++)
          {
            rec_resB[j][i] = img->m7[j][i];
            img->m7[j][i]  = resTrans_R[j][i];
          }
        }
        cbp_chroma_block[1][y_pos][x_pos] = dct_chroma4x4 (1, b8+8, b4);
        dc_level        [1][y_pos][x_pos] = dc_level_temp[1][y_pos][x_pos];
        for (j=0; j<4; j++)
        {
          for (i=0; i<4; i++)
          {
            rec_resR[j][i] = img->m7[j][i];
          }
        }
        for (j=0; j<4; j++)
        {
          for (i=0; i<4; i++)
          {
            /* Inverse Residue Transform */
            temp      = rec_resG[j][i]-(rec_resB[j][i]>>1);
            residue_G = rec_resB[j][i]+temp;
            residue_B = temp - (rec_resR[j][i]>>1);
            residue_R = residue_B+rec_resR[j][i];
            enc_picture->imgUV[0][pic_pix_y+j][pic_pix_x+i] = min(img->max_imgpel_value_uv,max(0,residue_B+(int)img->mprr_c[0][c_ipmode][block_y+j][block_x+i]));
            enc_picture->imgY[pic_pix_y+j][pic_pix_x+i]     = min(img->max_imgpel_value,max(0,residue_G+(int)img->mprr[best_ipmode][j][i]));
            enc_picture->imgUV[1][pic_pix_y+j][pic_pix_x+i] = min(img->max_imgpel_value_uv,max(0,residue_R+(int)img->mprr_c[1][c_ipmode][block_y+j][block_x+i]));
          }
        }
      }
  }
  else
  {
    //===== restore coefficients =====
    for (j=0; j<2; j++)
    {
      memcpy (img->cofAC[b8][b4][j],cofAC4x4[j], 18 * sizeof(int));
    } 
    // Residue Color Transform
    if(img->residue_transform_flag)
    {
      for (j=0; j<2; j++)
      {
        memcpy (img->cofAC[b8+4][b4][j],cofAC4x4_chroma[0][j], 18 * sizeof(int));            
        memcpy (img->cofAC[b8+8][b4][j],cofAC4x4_chroma[1][j], 18 * sizeof(int));            
      }
    }
    
    //===== restore reconstruction and prediction (needed if single coeffs are removed) =====
    for (y=0; y<BLOCK_SIZE; y++)
    {
      memcpy (&enc_picture->imgY[pic_pix_y+y][pic_pix_x],rec4x4[y],    BLOCK_SIZE * sizeof(imgpel));
      memcpy (&img->mpr[block_y+y][block_x],img->mprr[best_ipmode][y], BLOCK_SIZE * sizeof(imgpel));
      if(img->type==SP_SLICE &&(!si_frame_indicator && !sp2_frame_indicator))
        memcpy (&lrec[pic_pix_y+y][pic_pix_x],lrec4x4[y], BLOCK_SIZE * sizeof(int));//restore coefficients when encoding primary SP frame
    }
    
    if (img->AdaptiveRounding)
    {
      for (j=0; j<BLOCK_SIZE; j++)
        memcpy (&img->fadjust4x4[1][block_y+j][block_x],&fadjust4x4[block_y+j][block_x], BLOCK_SIZE * sizeof(int)); 
    }
    
    // Residue Color Transform
    if(img->residue_transform_flag)
    {
      for (i=0; i<2; i++)
      { //uv
        //--- set reconstruction ---
        for (y=0; y<4; y++)
        {
          memcpy(&enc_picture->imgUV[i][pic_pix_y+y][pic_pix_x],rec4x4_c[i][y], BLOCK_SIZE * sizeof(imgpel));
          if(img->type==SP_SLICE &&(!si_frame_indicator && !sp2_frame_indicator))
            memcpy (&lrec_uv[i][pic_pix_y+y][pic_pix_x],lrec4x4_c[i][y], BLOCK_SIZE * sizeof(int));//restore coefficients for primary SP frame encoding
        }
      }
    }    
  }  
  return nonzero;
}


/*!
*************************************************************************************
* \brief
*    Mode Decision for an 8x8 Intra block
*************************************************************************************
*/
int Mode_Decision_for_8x8IntraBlocks(int b8,double lambda,int *cost)
{
  int  nonzero=0, b4;
  int  cost4x4;
  
  *cost = (int)floor(6.0 * lambda + 0.4999);
  
  for (b4=0; b4<4; b4++)
  {
    if (Mode_Decision_for_4x4IntraBlocks (b8, b4, lambda, &cost4x4))
    {
      nonzero        = 1;
    }
    *cost += cost4x4;
  }
#ifdef RESET_STATE
  reset_coding_state (cs_cm);
#endif
  
  return nonzero;
}

/*!
*************************************************************************************
* \brief
*    4x4 Intra mode decision for an macroblock
*************************************************************************************
*/
int Mode_Decision_for_Intra4x4Macroblock (double lambda,  int* cost)
{
  int  cbp=0, b8, cost8x8;
  
  for (*cost=0, b8=0; b8<4; b8++)
  {
    if (Mode_Decision_for_8x8IntraBlocks (b8, lambda, &cost8x8))
    {
      cbp |= (1<<b8);
    }
    *cost += cost8x8;
  }
  
  return cbp;
}


/*!
*************************************************************************************
* \brief
*    R-D Cost for an 8x8 Partition
*************************************************************************************
*/
double RDCost_for_8x8blocks (int*    cnt_nonz,   // --> number of nonzero coefficients
                             int64*  cbp_blk,    // --> cbp blk
                             double  lambda,     // <-- lagrange multiplier
                             int     block,      // <-- 8x8 block number
                             int     mode,       // <-- partitioning mode
                             short   pdir,       // <-- prediction direction
                             short   ref,        // <-- reference frame
                             short   bwd_ref)    // <-- abp type
{
  int  i, j, k;
  int  rate=0;
  int64 distortion=0;
  int  dummy = 0, mrate;
  int  fw_mode, bw_mode;
  int  cbp     = 0;
  int  pax     = 8*(block & 0x01);
  int  pay     = 8*(block >> 1);
  int  i0      = pax >> 2;
  int  j0      = pay >> 2;
  int  bframe  = (img->type==B_SLICE);
  int  direct  = (bframe && mode==0);
  int  b8value = B8Mode2Value (mode, pdir);
  
  Macroblock    *currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement *currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice         *currSlice = img->currentSlice;
  DataPartition *dataPart;
  const int     *partMap   = assignSE2partition[input->partition_mode];
  
  EncodingEnvironmentPtr eep_dp;
  
  // Residue Color Transform
  int residue_R, residue_G, residue_B, temp, b4;
  int b4_x, b4_y;
  
  //=====
  //=====  GET COEFFICIENTS, RECONSTRUCTIONS, CBP
  //=====
  currMB->bi_pred_me=0;
  
  if (direct)
  {
    if (direct_pdir[img->block_y+j0][img->block_x+i0]<0) // mode not allowed
      return (1e20);
    else
      *cnt_nonz = LumaResidualCoding8x8 (&cbp, cbp_blk, block, direct_pdir[img->block_y+j0][img->block_x+i0], 0, 0, 
      (short)max(0,direct_ref_idx[LIST_0][img->block_y+j0][img->block_x+i0]), direct_ref_idx[LIST_1][img->block_y+j0][img->block_x+i0]);
  }
  else
  {
    fw_mode   = (pdir==0||pdir==2 ? mode : 0);
    bw_mode   = (pdir==1||pdir==2 ? mode : 0);
    *cnt_nonz = LumaResidualCoding8x8 (&cbp, cbp_blk, block, pdir, fw_mode, bw_mode, ref, bwd_ref);
  }
  
  // Residue Color Transform
  if(img->residue_transform_flag)
  {
    for(b4 = 0; b4 < 4; b4++)
    {
      b4_x = pax+(b4 & 0x01)*4;
      b4_y = pay+(b4 >> 1  )*4;
      for (j=0; j<4; j++)
      {
        for (i=0; i<4; i++)
          img->m7[j][i] = resTrans_B[j+b4_y][i+b4_x];
      }
      rate += RDCost_for_4x4Blocks_Chroma (block+4, b4, 0);
      
      for (j=0; j<4; j++)
      {
        for (i=0; i<4; i++)
        {
          rec_resB[j+b4_y][i+b4_x] = img->m7[j][i];
          img->m7[j][i]  = resTrans_R[j+b4_y][i+b4_x];
        }
      }
      rate += RDCost_for_4x4Blocks_Chroma (block+8, b4, 1);
      
      for (j=0; j<4; j++)
      {
        for (i=0; i<4; i++)
        {
          rec_resR[j+b4_y][i+b4_x] = img->m7[j][i];
        }
      }
    }
    reset_coding_state (cs_cm);  
    /* Inverse Residue Transform */
    for (j=pay; j<pay+8; j++)
      for (i=pax; i<pax+8; i++)
      {
        /* YCoCg-R */
        temp      = rec_resG[j][i]-(rec_resB[j][i]>>1);
        residue_G = rec_resB[j][i]+temp;
        residue_B = temp - (rec_resR[j][i]>>1);
        residue_R = residue_B+rec_resR[j][i];
        
        enc_picture->imgUV[0][img->pix_y+j][img->pix_x+i] = min(img->max_imgpel_value_uv,max(0,residue_B+mprRGB[1][j][i]));
        enc_picture->imgY[img->pix_y+j][img->pix_x+i]     = min(img->max_imgpel_value,max(0,residue_G+mprRGB[0][j][i]));
        enc_picture->imgUV[1][img->pix_y+j][img->pix_x+i] = min(img->max_imgpel_value_uv,max(0,residue_R+mprRGB[2][j][i]));
      }
  }
  
  //===== get residue =====
  if (input->rdopt==3 && img->type!=B_SLICE)
  {
    // We need the reconstructed prediction residue for the simulated decoders.
    compute_residue_b8block (block, -1);
  }
  
  //=====
  //=====   GET DISTORTION
  //=====
  if (input->rdopt==3 && img->type!=B_SLICE)
  {
    for (k=0; k<input->NoOfDecoders ;k++)
    {
      decode_one_b8block (k, P8x8, block, mode, ref);
      for (j=img->opix_y+pay; j<img->opix_y+pay+8; j++)
        for (i=img->opix_x+pax; i<img->opix_x+pax+8; i++)
        {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          distortion += SQR_DEPTH(imgY_org[j][i], decs->decY[k][j][i], input->BitDepthLuma, img->BitDepthIncrease);
#else
          distortion += img->quad[imgY_org[j][i] - decs->decY[k][j][i]];
#endif
        }
    }
    distortion /= input->NoOfDecoders;
  }
  else
  {
    for (j=pay; j<pay+8; j++)
      for (i=img->pix_x+pax; i<img->pix_x+pax+8; i++)
      {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        distortion += SQR_DEPTH(imgY_org[img->opix_y+j][i], enc_picture->imgY[img->pix_y+j][i], input->BitDepthLuma, img->BitDepthIncrease);
#else
        distortion += img->quad [imgY_org[img->opix_y+j][i] - enc_picture->imgY[img->pix_y+j][i]];
#endif
        // Residue Color Transform
        if(img->residue_transform_flag)
        {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          distortion += SQR_DEPTH(imgUV_org[0][img->opix_y+j][i], enc_picture->imgUV[0][img->pix_y+j][i], input->BitDepthChroma, img->BitDepthIncreaseChroma);
          distortion += SQR_DEPTH(imgUV_org[1][img->opix_y+j][i], enc_picture->imgUV[1][img->pix_y+j][i], input->BitDepthChroma, img->BitDepthIncreaseChroma);
#else
          distortion += img->quad [imgUV_org[0][img->opix_y+j][i] - enc_picture->imgUV[0][img->pix_y+j][i]];
          distortion += img->quad [imgUV_org[1][img->opix_y+j][i] - enc_picture->imgUV[1][img->pix_y+j][i]];
#endif
        }
      }
  }
  
  //=====
  //=====   GET RATE
  //=====
  //----- block 8x8 mode -----
  if (input->symbol_mode == UVLC)
  {
    ue_linfo (b8value, dummy, &mrate, &dummy);
    rate += mrate;
  }
  else
  {
    currSE->value1  = b8value;
    currSE->writing = writeB8_typeInfo_CABAC;
    currSE->type    = SE_MBTYPE;
    dataPart = &(currSlice->partArr[partMap[currSE->type]]);
    dataPart->writeSyntaxElement (currSE, dataPart);
    rate += currSE->len;
    currSE++;
    currMB->currSEnr++;
  }
  
  //----- motion information -----
  if (!direct)
  {
    if ((img->num_ref_idx_l0_active > 1 ) && (pdir==0 || pdir==2))
      rate  += writeReferenceFrame (mode, i0, j0, 1, ref);
    if(img->num_ref_idx_l1_active > 1 && img->type== B_SLICE)
    {
      if (pdir==1 || pdir==2)
      {
        rate  += writeReferenceFrame (mode, i0, j0, 0, bwd_ref);
      }
    }
#ifdef MV_COMPETITION
    if (input->mv_competition)
      pass_with_writing = 0;
#endif
    if (pdir==0 || pdir==2)
    {
      rate  += writeMotionVector8x8 (i0, j0, i0+2, j0+2, ref,LIST_0, mode);
    }
    if (pdir==1 || pdir==2)
    {
      rate  += writeMotionVector8x8 (i0, j0, i0+2, j0+2, bwd_ref, LIST_1, mode);
    }
#ifdef MV_COMPETITION
    if (input->mv_competition)
      pass_with_writing = 1;
#endif
  }
  
  //----- coded block pattern (for CABAC only) -----
  if (input->symbol_mode == CABAC)
  {
    dataPart = &(currSlice->partArr[partMap[SE_CBP_INTER]]);
    eep_dp   = &(dataPart->ee_cabac);
    mrate    = arienco_bits_written (eep_dp);
    writeCBP_BIT_CABAC (block, ((*cnt_nonz>0)?1:0), cbp8x8, currMB, 1, eep_dp);
    mrate    = arienco_bits_written (eep_dp) - mrate;
    rate    += mrate;
  }
  
  //----- luminance coefficients -----
  if (*cnt_nonz)
  {
    rate += writeLumaCoeff8x8 (block, mode, currMB->luma_transform_size_8x8_flag);
  }
  
  return (double)distortion + lambda * (double)rate;
}


/*!
*************************************************************************************
* \brief
*    Gets mode offset for intra16x16 mode
*************************************************************************************
*/
int I16Offset (int cbp, int i16mode)
{
  return (cbp&15?13:1) + i16mode + ((cbp&0x30)>>2);
}


/*!
*************************************************************************************
* \brief
*    Sets modes and reference frames for a macroblock
*************************************************************************************
*/
void SetModesAndRefframeForBlocks (int mode)
{
  int i,j,k,l;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int  bframe  = (img->type==B_SLICE);
  int  block_x, block_y;
  int  cur_ref[2];  
  
  //--- macroblock type ---
  currMB->mb_type = mode;    
  currMB->bi_pred_me= (mode == 1 ? img->bi_pred_me[mode] : 0);  
  
  //--- block 8x8 mode and prediction direction ---
  switch (mode)
  {
  case 0:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] = 0;
      currMB->b8pdir[i] = (bframe ? direct_pdir[img->block_y + (i >> 1)*2][img->block_x + (i & 0x01)*2] : 0);
    }
    break;
  case 1:
  case 2:
  case 3:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] = mode;
      currMB->b8pdir[i] = best8x8pdir[mode][i];
    }
    break;
  case P8x8:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i]   = best8x8mode[i];
      currMB->b8pdir[i]   = best8x8pdir[mode][i];
    }
    break;
  case I4MB:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] = IBLOCK; 
      currMB->b8pdir[i] = -1;
    }
    break;
  case I16MB:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] =  0;
      currMB->b8pdir[i] = -1;
    }
    break;
  case I8MB:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] = I8MB;
      currMB->b8pdir[i] = -1;
    }
    //switch to 8x8 transform
    currMB->luma_transform_size_8x8_flag = 1;
    break;
  case IPCM:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] = IPCM;
      currMB->b8pdir[i] = -1;
    }
    currMB->luma_transform_size_8x8_flag = 0;
    break;
  default:
    printf ("Unsupported mode in SetModesAndRefframeForBlocks!\n");
    exit (1);
  }
  
#define IS_FW ((best8x8pdir[mode][k]==0 || best8x8pdir[mode][k]==2) && (mode!=P8x8 || best8x8mode[k]!=0 || !bframe))
#define IS_BW ((best8x8pdir[mode][k]==1 || best8x8pdir[mode][k]==2) && (mode!=P8x8 || best8x8mode[k]!=0))
  //--- reference frame arrays ---
  if (mode==0 || mode==I4MB || mode==I16MB || mode==I8MB)
  {
    if (bframe)
    {
      if (!mode)
      {
        for (j = img->block_y; j < img->block_y + 4; j++)
        {
          memcpy(&enc_picture->ref_idx[LIST_0][j][img->block_x],&direct_ref_idx[LIST_0][j][img->block_x], 4 * sizeof(char));
          memcpy(&enc_picture->ref_idx[LIST_1][j][img->block_x],&direct_ref_idx[LIST_1][j][img->block_x], 4 * sizeof(char));
        }
      }
      else
      {
        for (j = img->block_y; j < img->block_y + 4; j++)
        {
          memset(&enc_picture->ref_idx[LIST_0][j][img->block_x],-1, 4 * sizeof(char));
          memset(&enc_picture->ref_idx[LIST_1][j][img->block_x],-1, 4 * sizeof(char));
        }
      }
    }
    else
    {
      if (!mode)
      {
        for (j = img->block_y; j < img->block_y + 4; j++)
          memset(&enc_picture->ref_idx[LIST_0][j][img->block_x],0, 4 * sizeof(char));
      }
      else
      {
        for (j = img->block_y; j < img->block_y + 4; j++)
          memset(&enc_picture->ref_idx[LIST_0][j][img->block_x],-1, 4 * sizeof(char));
      }
    }
  }
  else
  {
    if (bframe)
    {
      for (j=0;j<4;j++)
      {
        block_y = img->block_y + j;
        for (i=0;i<4;i++)
        {
          block_x = img->block_x + i;
          k = 2*(j >> 1) + (i >> 1);
          l = 2*(j & 0x01) + (i & 0x01);
          
          if(mode == P8x8 && best8x8mode[k]==0)
          {
            enc_picture->ref_idx[LIST_0][block_y][block_x] = direct_ref_idx[LIST_0][block_y][block_x];
            enc_picture->ref_idx[LIST_1][block_y][block_x] = direct_ref_idx[LIST_1][block_y][block_x];
          }
          else if (mode ==1 && currMB->bi_pred_me && IS_FW && IS_BW)
          {
            enc_picture->ref_idx[LIST_0][block_y][block_x] = 0;
            enc_picture->ref_idx[LIST_1][block_y][block_x] = 0;
          }
          else
          {
            enc_picture->ref_idx[LIST_0][block_y][block_x] = (IS_FW ? best8x8fwref[mode][k] : -1);
            enc_picture->ref_idx[LIST_1][block_y][block_x] = (IS_BW ? best8x8bwref[mode][k] : -1);
          }
        }
      }
    }
    else
    {
      for (j=0;j<4;j++)
      {
        block_y = img->block_y + j;
        for (i=0;i<4;i++)
        {
          block_x = img->block_x + i;
          k = 2*(j >> 1) + (i >> 1);
          l = 2*(j & 0x01) + (i & 0x01);
          enc_picture->ref_idx[LIST_0][block_y][block_x] = (IS_FW ? best8x8fwref[mode][k] : -1);
        }
      }
    }
  }
  
  if (bframe)
  {
    
    for (j = img->block_y; j < img->block_y + 4; j++)
      for (i = img->block_x; i < img->block_x + 4;i++)
      {
        cur_ref[LIST_0] = (int) enc_picture->ref_idx[LIST_0][j][i];
        cur_ref[LIST_1] = (int) enc_picture->ref_idx[LIST_1][j][i];
        
        enc_picture->ref_pic_id [LIST_0][j][i] = (cur_ref[LIST_0]>=0 
          ? enc_picture->ref_pic_num[LIST_0 + currMB->list_offset][cur_ref[LIST_0]]
          : -1);
        enc_picture->ref_pic_id [LIST_1][j][i] = (cur_ref[LIST_1]>=0 
          ? enc_picture->ref_pic_num[LIST_1 + currMB->list_offset][cur_ref[LIST_1]]
          : -1);
      }
  }
  else
  {  
    for (j = img->block_y; j < img->block_y + 4; j++)
      for (i = img->block_x; i < img->block_x + 4;i++)
      {
        cur_ref[LIST_0] = (int) enc_picture->ref_idx[LIST_0][j][i];
        enc_picture->ref_pic_id [LIST_0][j][i] = (cur_ref[LIST_0]>=0 
          ? enc_picture->ref_pic_num[LIST_0 + currMB->list_offset][cur_ref[LIST_0]]
          : -1);
      }
  }
  
#undef IS_FW
#undef IS_BW
}


#ifdef USE_INTRA_MDDT
/*! 
*************************************************************************************
* \brief
*   RD-optimized mode decision for INTRA16x16 mode
*
* \para Intra16x16_Mode_Decision_RDopt()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
void Intra16x16_Mode_Decision_RDopt (Macroblock* currMB, int* i16mode, double lambda)
{
  extern int MBType2Value (Macroblock* currMB);
  int i, j;
  
  int mode16, best_mode16;
  int max_mode = 4;
  double rdcost, min_rdcost = 1e30;
  int distortion, rate; 
  imgpel **imgY = enc_picture->imgY; 
  int x, y;
  extern CSptr cs_i16;
  SyntaxElement  *currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
  DataPartition*  dataPart;
  const int*      partMap    = assignSE2partition[input->partition_mode];
  Slice*          currSlice  = img->currentSlice;
  long quant_stat_best16x16[256];
  int  bestCofAC16[2][16*17];
  int  best_cbp=0;
  imgpel bestRec[MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  
  intrapred_luma_16x16 ();   /* make intra pred for all 4 new modes */
  
  //if(!img->residue_transform_flag)  // assume this flag is always OFF
  {
    PixelPos up;          //!< pixel position p(0,-1)
    PixelPos left[17];    //!< pixel positions p(-1, -1..15)
    
    int up_avail, left_avail, left_up_avail;
    
    for (i=0;i<17;i++)
    {
      getNeighbour(img->current_mb_nr, -1 ,  i-1 , 1, &left[i]);
    }
    
    getNeighbour(img->current_mb_nr, 0     ,  -1 , 1, &up);
    
    if (!(input->UseConstrainedIntraPred))
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
    
    
    best_mode16 = DC_PRED_16;
    for(mode16 = 0; mode16 < max_mode; mode16++)
    {
      if (input->IntraDisableInterOnly == 0 || img->type != I_SLICE)
      {
        if (input->Intra16x16ParDisable && (mode16==VERT_PRED_16||mode16==HOR_PRED_16))
          continue;
        
        if (input->Intra16x16PlaneDisable && mode16==PLANE_16)
          continue;
      }
      
      //check if there are neighbours to predict from
      if ((mode16==VERT_PRED_16 && !up_avail  ) || 
        (mode16==HOR_PRED_16  && !left_avail) || 
        (mode16==PLANE_16     && (!left_avail || !up_avail || !left_up_avail)))
      {
        continue; // edge, do nothing
      }
      
      if(input->symbol_mode == UVLC)
      {    
        currMB->cbp = dct_luma_16x16 (mode16);
        
      }   // if UVLC coding
      else
      {
        currMB->cbp = klt_luma16x16_sep_fast (mode16); 
        //===== get distortion (SSD) of 4x4 block =====
        distortion = 0;
        for (y=0; y<16; y++)
        {
          for (x=0; x<16; x++)
          {
            distortion += 
              img->quad [imgY_org[img->opix_y+y][img->opix_x+x] - imgY[img->pix_y+y][img->pix_x+x]];
          }
        }
        
        store_coding_state (cs_i16);
        
        // rate for MB_TYPE 
        img->i16offset = I16Offset  (currMB->cbp, mode16);
        currMB->mb_type = I16MB; 
        currSE->value1  = MBType2Value (currMB);
        currSE->value2  = 0;
        currSE->type    = SE_MBTYPE;
        
        if (input->symbol_mode == UVLC)  
          currSE->mapping = ue_linfo;
        else
          currSE->writing = writeMB_typeInfo_CABAC;
        dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);
        dataPart->writeSyntaxElement( currSE, dataPart);
        rate                = currSE->len;
        currSE++;
        currMB->currSEnr++;
        
        //===== RATE for LUMINANCE COEFFICIENTS =====
#ifdef MB32X32
        rate += writeCBPandLumaCoeff(1,0,0);
#else
        rate += writeCBPandLumaCoeff(1);
#endif
        reset_coding_state (cs_i16);
        
        rdcost = (double)distortion + lambda*(double)rate;
        if(rdcost < min_rdcost)
        {
          best_mode16 = mode16;
          min_rdcost = rdcost;
          best_cbp = currMB->cbp;
          for (j=0; j<256; j++)
          {
            quant_stat_best16x16[j]=quant_stat_rd16x16[j];
          }
          memcpy(bestCofAC16[0], img->cofAC16[0], sizeof(int)*16*17);
          memcpy(bestCofAC16[1], img->cofAC16[1], sizeof(int)*16*17);
          for(j = 0; j < MB_BLOCK_SIZE; j++)
            memcpy(bestRec[j], &imgY[img->pix_y+j][img->pix_x], sizeof(imgpel)*MB_BLOCK_SIZE);
        }
      }
    }
    
    *i16mode = best_mode16;
  }
  
  for (j=0; j<256; j++)
  {
    img->quant_stat16x16[j]+=(quant_stat_best16x16[j]>>5);
  }
  
  currMB->cbp = best_cbp;
  if(currMB->cbp)
    currMB->cbp = 15;
  memcpy(img->cofAC16[0], bestCofAC16[0], sizeof(int)*16*17);
  memcpy(img->cofAC16[1], bestCofAC16[1], sizeof(int)*16*17);
  for(y = 0; y < MB_BLOCK_SIZE; y++)
  {
    memcpy(&imgY[img->pix_y+y][img->pix_x], &bestRec[y][0], MB_BLOCK_SIZE*sizeof(imgpel));
  }
}
#endif

/*!
*************************************************************************************
* \brief
*    Intra 16x16 mode decision
*************************************************************************************
*/
void Intra16x16_Mode_Decision (Macroblock* currMB, int* i16mode)
{
  // Residue Color Transform
  int residue_R, residue_G, residue_B;
  int c_ipmode = img->mb_data[img->current_mb_nr].c_ipred_mode;
  int i, j, temp;
  int pic_pix_x   = img->pix_x;
  int pic_pix_y   = img->pix_y;
  pel_t   **imgY_orig  = imgY_org;
  pel_t   ***imgUV_orig  = imgUV_org;
  int cr_cbp;
  
  intrapred_luma_16x16 ();   /* make intra pred for all 4 new modes */
  
  if(!img->residue_transform_flag)
    find_sad_16x16 (i16mode);   /* get best new intra mode */
  
  // Residue Color Transform
  if(img->residue_transform_flag)
  {
    for (j=0; j < MB_BLOCK_SIZE; j++)
      for (i=0; i < MB_BLOCK_SIZE; i++)
      {
        residue_B = imgUV_orig[0][pic_pix_y+j][pic_pix_x+i] - img->mprr_c[0][c_ipmode][j][i];
        residue_G = imgY_orig[pic_pix_y+j][pic_pix_x+i] - img->mprr_2[*i16mode][j][i];
        residue_R = imgUV_orig[1][pic_pix_y+j][pic_pix_x+i] - img->mprr_c[1][c_ipmode][j][i];
        
        /* Forward Residue Transform */
        resTrans_R[j][i] = residue_R-residue_B;
        temp = residue_B+(resTrans_R[j][i]>>1);
        resTrans_B[j][i] = residue_G-temp;
        resTrans_G[j][i] = temp+(resTrans_B[j][i]>>1);
        
        img->m7[j][i]  = resTrans_G[j][i];
      }
  }
  
  currMB->cbp = dct_luma_16x16 (*i16mode);
  
  // Residue Color Transform
  if(img->residue_transform_flag)
  {
    for (j=0; j < MB_BLOCK_SIZE; j++)
    {
      for (i=0; i < MB_BLOCK_SIZE; i++)
      {
        rec_resG[j][i] = img->m7[j][i];
        img->m7[j][i]  = resTrans_B[j][i];
      }
    }
    cr_cbp = dct_chroma(0, 0);
    
    for (j=0; j < MB_BLOCK_SIZE; j++)
    {
      for (i=0; i < MB_BLOCK_SIZE; i++)
      {
        rec_resB[j][i] = img->m7[j][i];
        img->m7[j][i]  = resTrans_R[j][i];
      }
    }
    cr_cbp = dct_chroma(1, cr_cbp);
    
    for (j = 0; j < MB_BLOCK_SIZE; j++)
    {
      for (i = 0; i < MB_BLOCK_SIZE; i++)
        
        rec_resR[j][i] = img->m7[j][i];
    }
    
    currMB->cbp += (cr_cbp<<4);
    
    /* Inverse Residue Transform */
    for (j = 0; j < MB_BLOCK_SIZE; j++)
    {
      for (i = 0; i < MB_BLOCK_SIZE; i++)
      {
        temp      = rec_resG[j][i]-(rec_resB[j][i]>>1);
        residue_G = rec_resB[j][i]+temp;
        residue_B = temp - (rec_resR[j][i]>>1);
        residue_R = residue_B+rec_resR[j][i];
        
        enc_picture->imgUV[0][pic_pix_y+j][pic_pix_x+i] = min(img->max_imgpel_value_uv,max(0,residue_B+(int)img->mprr_c[0][c_ipmode][j][i]));
        enc_picture->imgY[pic_pix_y+j][pic_pix_x+i]     = min(img->max_imgpel_value,max(0,residue_G+(int)img->mprr_2[*i16mode][j][i]));
        enc_picture->imgUV[1][pic_pix_y+j][pic_pix_x+i] = min(img->max_imgpel_value_uv,max(0,residue_R+(int)img->mprr_c[1][c_ipmode][j][i]));
      }
    }
  }
}


/*!
*************************************************************************************
* \brief
*    Sets Coefficients and reconstruction for an 8x8 block
*************************************************************************************
*/
void SetCoeffAndReconstruction8x8 (Macroblock* currMB)
{
  int block, k, j, i;
  int cur_ref[2];
  
  //============= MIXED TRANSFORM SIZES FOR 8x8 PARTITION ==============
  //--------------------------------------------------------------------
  int l;
  int bframe = img->type==B_SLICE; 
  
  if (currMB->luma_transform_size_8x8_flag)
  {
    
    //============= set mode and ref. frames ==============
    for(i = 0;i<4;i++)
    {
      currMB->b8mode[i]   = tr8x8.part8x8mode[i];
      currMB->b8pdir[i]   = tr8x8.part8x8pdir[i];
    }
    
    if (bframe)
    {
      for (j = 0;j<4;j++)
        for (i = 0;i<4;i++)
        {
          k = 2*(j >> 1)+(i >> 1);
          l = 2*(j & 0x01)+(i & 0x01);
          enc_picture->ref_idx[LIST_0][img->block_y+j][img->block_x+i] = ((currMB->b8pdir[k] & 0x01) == 0) ? tr8x8.part8x8fwref[k] : - 1;
          enc_picture->ref_idx[LIST_1][img->block_y+j][img->block_x+i] = (currMB->b8pdir[k] > 0) ? tr8x8.part8x8bwref[k] : - 1;
        }
    }
    else
    {
      for (j = 0;j<4;j++)
        for (i = 0;i<4;i++)
        {
          k = 2*(j >> 1)+(i >> 1);
          l = 2*(j & 0x01)+(i & 0x01);
          enc_picture->ref_idx[LIST_0][img->block_y+j][img->block_x+i] = tr8x8.part8x8fwref[k];
        }
    }
    
    
    for (j = img->block_y;j<img->block_y + BLOCK_MULTIPLE;j++)
    {
      for (i = img->block_x;i<img->block_x + BLOCK_MULTIPLE;i++)
      {
        cur_ref[LIST_0] = (int) enc_picture->ref_idx[LIST_0][j][i];
        
        enc_picture->ref_pic_id [LIST_0][j][i] =(cur_ref[LIST_0]>=0 
          ? enc_picture->ref_pic_num[LIST_0 + currMB->list_offset][cur_ref[LIST_0]]
          : -1);
      }
    }
    
    if (bframe)
    {
      for (j = img->block_y; j < img->block_y + BLOCK_MULTIPLE; j++)
      {
        for (i = img->block_x;i<img->block_x + BLOCK_MULTIPLE;i++)
        {
          cur_ref[LIST_1] = (int) enc_picture->ref_idx[LIST_1][j][i];
          
          enc_picture->ref_pic_id [LIST_1][j][i] = (cur_ref[LIST_1]>=0 
            ? enc_picture->ref_pic_num[LIST_1 + currMB->list_offset][cur_ref[LIST_1]]
            : -1);
        }
        
      }
    }
    
    //====== set the mv's for 8x8 partition with transform size 8x8 ======
    //save the mv data for 4x4 transform
    StoreMV8x8(1);
    //set new mv data for 8x8 transform
    RestoreMV8x8(0);
    
    //============= get pre-calculated data ==============
    //restore coefficients from 8x8 transform
    
    for (block = 0; block<4; block++)
    {
      for (k = 0; k<4; k++)
        for (j = 0; j<2; j++)
          memcpy (img->cofAC[block][k][j],cofAC_8x8ts[block][k][j], 65 * sizeof(int));            
    }     
    //restore reconstruction 
    if (cnt_nonz8_8x8ts <= _LUMA_8x8_COEFF_COST_ && 
      ((img->qp + img->bitdepth_luma_qp_scale)!=0 || img->lossless_qpprime_flag==0) &&
      (img->type!=SP_SLICE))// modif ES added last condition (we probably never go there so is the next modification useful ? check)
    {
      currMB->cbp     = 0;
      currMB->cbp_blk = 0;
      if(!img->residue_transform_flag) // Residue Color Transform
      {
        for (j = 0; j < MB_BLOCK_SIZE; j++)
        {
          memcpy(&enc_picture->imgY[img->pix_y+j][img->pix_x], tr8x8.mpr8x8[j], MB_BLOCK_SIZE * sizeof(imgpel));
          if(img->type==SP_SLICE &&(!si_frame_indicator && !sp2_frame_indicator ))
            memcpy(&lrec[img->pix_y+j][img->pix_x],tr8x8.lrec[j], MB_BLOCK_SIZE * sizeof(int));
        }
      }
      else
      {
        for (j = 0; j < MB_BLOCK_SIZE; j++)
          memset(rec_resG[j], 0, MB_BLOCK_SIZE * sizeof(int));
      }                  
    }
    else
    {
      currMB->cbp     = cbp8_8x8ts;
      currMB->cbp_blk = cbp_blk8_8x8ts;
      
#ifdef ADAPTIVE_FD_SD_CODING
      currMB->SD_Coding_on_off=tr8x8.SD_Coding_on_off8x8;
      
      memcpy(currMB->quantizer_indices,tr8x8.quantizer_indices8x8,16*16 * sizeof(int));
      
      memcpy(currMB->SD_or_FD,tr8x8.SD_or_FD8x8,2*2 * sizeof(int));
      
      currMB->SD_or_FD_t8x8=tr8x8.SD_or_FD_t8x88x8;
#endif
      
      if(!img->residue_transform_flag)                // Residue Color Transform
      {
        for (j = 0; j < MB_BLOCK_SIZE; j++)
        {
          memcpy (&enc_picture->imgY[img->pix_y+j][img->pix_x],tr8x8.rec_mbY8x8[j], MB_BLOCK_SIZE * sizeof(imgpel));   
          if(img->type==SP_SLICE &&(!si_frame_indicator && !sp2_frame_indicator))
            memcpy (&lrec[img->pix_y+j][img->pix_x],tr8x8.lrec[j], MB_BLOCK_SIZE * sizeof(int)); 
        }
      }
      else
      {
        for (j = 0; j < MB_BLOCK_SIZE; j++)
          memcpy (rec_resG[j], tr8x8.rec_resG_8x8[j], MB_BLOCK_SIZE * sizeof(int));
      }
    }
    
    if(img->residue_transform_flag)                // Residue Color Transform
    {
      // Residue Color Transform
      for (j = 0; j < MB_BLOCK_SIZE; j++)
        for (i = 0; i < MB_BLOCK_SIZE; i++)
        {
          mprRGB[0][j][i]  = tr8x8.mprRGB_8x8[0][j][i];
          mprRGB[1][j][i]  = tr8x8.mprRGB_8x8[1][j][i];
          mprRGB[2][j][i]  = tr8x8.mprRGB_8x8[2][j][i];
          resTrans_R[j][i] = tr8x8.resTrans_R_8x8[j][i];
          resTrans_B[j][i] = tr8x8.resTrans_B_8x8[j][i];
        }
    }
  }
  else
  {
    //============= get pre-calculated data ==============
    //---------------------------------------------------
    //--- restore coefficients ---
    for (block = 0; block<4+img->num_blk8x8_uv; block++)
    {
      for (k = 0; k<4; k++)
        for (j = 0; j<2; j++)
          memcpy (img->cofAC[block][k][j],cofAC8x8[block][k][j], 65 * sizeof(int));            
    }
    
    if (cnt_nonz_8x8<=5 && img->type!=SP_SLICE &&
      ((img->qp + img->bitdepth_luma_qp_scale)!=0 || img->lossless_qpprime_flag==0))
    {
      currMB->cbp     = 0;
      currMB->cbp_blk = 0;
      if(!img->residue_transform_flag) // Residue Color Transform
      {
        for (j = 0; j < MB_BLOCK_SIZE; j++)
        {
          memcpy (&enc_picture->imgY[img->pix_y+j][img->pix_x],tr4x4.mpr8x8[j], MB_BLOCK_SIZE * sizeof(imgpel));            
          if(img->type ==SP_SLICE &&(!si_frame_indicator && !sp2_frame_indicator))
            memcpy (&lrec[img->pix_y+j][img->pix_x],tr4x4.lrec[j], MB_BLOCK_SIZE * sizeof(int)); // restore coeff. SP frame
        }
      }
      else
      {
        for (j = 0; j < MB_BLOCK_SIZE; j++)
          memset(rec_resG[j], 0, MB_BLOCK_SIZE * sizeof(int));
      }
    }
    else
    {
      currMB->cbp     = cbp8x8;
      currMB->cbp_blk = cbp_blk8x8;
      
#ifdef ADAPTIVE_FD_SD_CODING
      currMB->SD_Coding_on_off=tr4x4.SD_Coding_on_off8x8;
      
      
      memcpy(currMB->quantizer_indices,tr4x4.quantizer_indices8x8  ,16*16*sizeof(int));
      
      for (i=0;i<2;i++)
        for (j=0;j<2;j++)
        {
          currMB->SD_or_FD[j][i]=tr4x4.SD_or_FD8x8[j][i];
        }
        
        currMB->SD_or_FD_t8x8=tr4x4.SD_or_FD_t8x88x8;
#endif
        
        if(!img->residue_transform_flag)           // Residue Color Transform
        {
          for (j = 0; j < MB_BLOCK_SIZE; j++)
          {
            memcpy (&enc_picture->imgY[img->pix_y+j][img->pix_x],tr4x4.rec_mbY8x8[j], MB_BLOCK_SIZE * sizeof(imgpel));            
            if(img->type==SP_SLICE &&(!si_frame_indicator && !sp2_frame_indicator))
              memcpy (&lrec[img->pix_y+j][img->pix_x],tr4x4.lrec[j], MB_BLOCK_SIZE * sizeof(int));
          }
        }
        else
        {
          for (j = 0; j < MB_BLOCK_SIZE; j++)
            memcpy(rec_resG[j], tr4x4.rec_resG_8x8[j], MB_BLOCK_SIZE * sizeof(int));
        }
    }
    
    if(img->residue_transform_flag)           // Residue Color Transform
    {
      // Residue Color Transform
      for (j = 0; j < MB_BLOCK_SIZE; j++)
      {
        memcpy(mprRGB[0][j], tr4x4.mprRGB_8x8[0][j], MB_BLOCK_SIZE * sizeof(int));
        memcpy(mprRGB[1][j], tr4x4.mprRGB_8x8[1][j], MB_BLOCK_SIZE * sizeof(int));
        memcpy(mprRGB[2][j], tr4x4.mprRGB_8x8[2][j], MB_BLOCK_SIZE * sizeof(int));
        memcpy(resTrans_R[j], tr4x4.resTrans_R_8x8[j], MB_BLOCK_SIZE * sizeof(int));
        memcpy(resTrans_B[j], tr4x4.resTrans_B_8x8[j], MB_BLOCK_SIZE * sizeof(int));
      }
    }
  }
}


/*!
*************************************************************************************
* \brief
*    Sets motion vectors for a macroblock
*************************************************************************************
*/
void SetMotionVectorsMB (Macroblock* currMB, int bframe)
{
  int i, j, k, l, m, mode8, pdir8, ref, by, bx;
  short ******all_mv  = img->all_mv;
  short ******pred_mv = img->pred_mv;
  int  bw_ref;
  int jdiv, jmod;
  
  if (!bframe)
  {
    for (j = 0; j<4; j++)
    {
      jmod = j & 0x01;
      jdiv = j >>   1;
      by    = img->block_y+j;
      for (i = 0; i<4; i++)
      {
        mode8 = currMB->b8mode[k=2*jdiv+(i>>1)];
        l     = 2*jmod + (i & 0x01);
        
        bx   = img->block_x+i;
        
        pdir8 = currMB->b8pdir[k];
        ref    = enc_picture->ref_idx[LIST_0][by][bx];
        
        if (pdir8>=0) 
        {
          enc_picture->mv[LIST_0][by][bx][0] = all_mv [j][i][LIST_0][ ref][mode8][0];
          enc_picture->mv[LIST_0][by][bx][1] = all_mv [j][i][LIST_0][ ref][mode8][1];
        }
        else
        {
          enc_picture->mv[LIST_0][by][bx][0] = 0;
          enc_picture->mv[LIST_0][by][bx][1] = 0;
        }
      }
    }
  }
  else
  {
    for (j = 0; j<4; j++)
    {
      jmod = j & 0x01;
      jdiv = j >>   1;
      by    = img->block_y+j;
      for (i = 0; i<4; i++)
      {
        mode8 = currMB->b8mode[k=2*jdiv+(i>>1)];
        l     = 2*jmod + (i & 0x01);
        
        bx    = img->block_x+i;
        
        pdir8 = currMB->b8pdir[k];
        ref    = enc_picture->ref_idx[LIST_0][by][bx];
        bw_ref = enc_picture->ref_idx[LIST_1][by][bx];
        
        if (currMB->bi_pred_me && (pdir8 == 2) && currMB->mb_type==1)
        {
          all_mv  = currMB->bi_pred_me == 1 ? img->bipred_mv1 : img->bipred_mv2;
          ref = 0;
          bw_ref = 0;
        }
        
        if (pdir8==-1) // intra
        {
          enc_picture->mv[LIST_0][by][bx][0] = 0;
          enc_picture->mv[LIST_0][by][bx][1] = 0;
          enc_picture->mv[LIST_1][by][bx][0] = 0;
          enc_picture->mv[LIST_1][by][bx][1] = 0;
        }
        else if (pdir8==0) // list 0
        {
          enc_picture->mv[LIST_0][by][bx][0] = all_mv [j][i][LIST_0][ ref][mode8][0];
          enc_picture->mv[LIST_0][by][bx][1] = all_mv [j][i][LIST_0][ ref][mode8][1];
          enc_picture->mv[LIST_1][by][bx][0] = 0;
          enc_picture->mv[LIST_1][by][bx][1] = 0;
          enc_picture->ref_idx[LIST_1][by][bx] = -1;
        }
        else if (pdir8==1) // list 1
        {
          enc_picture->mv[LIST_0][by][bx][0] = 0;
          enc_picture->mv[LIST_0][by][bx][1] = 0;          
          enc_picture->ref_idx[LIST_0][by][bx] = -1;
          enc_picture->mv[LIST_1][by][bx][0] = all_mv [j][i][LIST_1][bw_ref][mode8][0];
          enc_picture->mv[LIST_1][by][bx][1] = all_mv [j][i][LIST_1][bw_ref][mode8][1];
        }
        else if (pdir8==2) // bipredictive
        {
          enc_picture->mv[LIST_0][by][bx][0] = all_mv [j][i][LIST_0][ ref][mode8][0];
          enc_picture->mv[LIST_0][by][bx][1] = all_mv [j][i][LIST_0][ ref][mode8][1];                   
          enc_picture->mv[LIST_1][by][bx][0] = all_mv [j][i][LIST_1][bw_ref][mode8][0];
          enc_picture->mv[LIST_1][by][bx][1] = all_mv [j][i][LIST_1][bw_ref][mode8][1];
        }
        else
        {
          error("invalid direction mode", 255);
        }
        
      }
    }
  }
  
  // copy all the motion vectors into rdopt structure
  // Can simplify this by copying the MV's of the best mode (TBD)
  if(img->MbaffFrameFlag)
  {
    for(i = 0;i<4;i++)
    {
      for(j = 0;j<4;j++)
      {
        for (k = 0;k<2;k++)
        {
          for(l = 0;l<img->max_num_references;l++)
          {
            for(m = 0;m<9;m++)
            {
              rdopt->all_mv [j][i][k][l][m][0]  = all_mv [j][i][k][l][m][0];
              rdopt->pred_mv[j][i][k][l][m][0]  = pred_mv[j][i][k][l][m][0];
              
              rdopt->all_mv [j][i][k][l][m][1]  = all_mv [j][i][k][l][m][1];
              rdopt->pred_mv[j][i][k][l][m][1]  = pred_mv[j][i][k][l][m][1];
            }
          }
        }
      }
    }
  }
}



/*!
*************************************************************************************
* \brief
*    R-D Cost for a macroblock
*************************************************************************************
*/
int RDCost_for_macroblocks (double   lambda,       // <-- lagrange multiplier
                            int      mode,         // <-- modus (0-COPY/DIRECT, 1-16x16, 2-16x8, 3-8x16, 4-8x8(+), 5-Intra4x4, 6-Intra16x16)
                            double*  min_rdcost,   // <-> minimum rate-distortion cost
                            double*  min_rate,     // --> bitrate of mode which has minimum rate-distortion cost. 
                            int i16mode )
{
  
  int         i, j, k; 
  int         j1, j2;
  int         rate = 0, coeff_rate = 0;
  int64       distortion = 0;
  double      rdcost;
  Macroblock  *currMB   = &img->mb_data[img->current_mb_nr];
  Macroblock  *prevMB   = img->current_mb_nr ? &img->mb_data[img->current_mb_nr-1] : NULL;
  int         bframe    = (img->type==B_SLICE);
  int         use_of_cc =  (img->type!=I_SLICE &&  input->symbol_mode!=CABAC);
  int         dummy;
  int         cr_cbp = 0, uv;
  int      return_val = 0;
  int tmp_cc;
  int cc_rate;
  
#ifdef MV_COMPETITION
  short iter, nbiter;
  short ******all_mv = img->all_mv;
#endif
  
  //=====
  //=====  SET REFERENCE FRAMES AND BLOCK MODES
  //=====
  SetModesAndRefframeForBlocks (mode);
  
  
#ifdef MV_COMPETITION
  if ((input->mv_competition > 0) && (mode == 0) && (img->type==P_SLICE))
    nbiter = mv_comp.nb_mode_for_skip + 1;
  else
    nbiter = 1;
  
  // Let's iterate on each possible mode for the skip
  for (iter= 0; iter < nbiter; iter++)
  {
    if ((input->mv_competition > 0) && (mode == 0) && (img->type==P_SLICE))  
    {
      distortion = 0;
      coeff_rate = 0;
      rate = 0;
      cr_cbp = 0;
    }
    
#endif
    //=====
    //=====  GET COEFFICIENTS, RECONSTRUCTIONS, CBP
    //=====
    if (bframe && mode==0)
    {
      int block_x=img->pix_x>>2;
      int block_y=img->pix_y>>2;    
      for (j = block_y;j< block_y + 4;j++)
        for (i = block_x;i<block_x + 4;i++)
          if (direct_pdir[j][i] < 0)
            return 0;
    }
    
    // Test MV limits for Skip Mode. This could be necessary for MBAFF case Frame MBs. 
    if ((img->MbaffFrameFlag) && (!currMB->mb_field) && (img->type==P_SLICE) && (mode==0) )
    {
      if ( img->all_mv[0][0][0][0][0][0] < - 8192 
        || img->all_mv[0][0][0][0][0][0] > 8191 
        || img->all_mv[0][0][0][0][0][1] < LEVELMVLIMIT[img->LevelIndex][4] 
        || img->all_mv[0][0][0][0][0][1] > LEVELMVLIMIT[img->LevelIndex][5])
        return 0;
    }
    
    if (img->AdaptiveRounding)
    {
#ifdef ADAPTIVE_FD_SD_CODING
      img->adjust_adaptive_f_spatial_domain_4x4=0;
      img->adjust_adaptive_f_spatial_domain_8x8=0;
#endif
      for (j = 0;j < MB_BLOCK_SIZE;j++)
      {
        memset(img->fadjust4x4[0][j], 0, MB_BLOCK_SIZE * sizeof(int));;
        memset(img->fadjust8x8[0][j], 0, MB_BLOCK_SIZE * sizeof(int));
        memset(img->fadjust4x4Cr[0][0][j], 0, MB_BLOCK_SIZE * sizeof(int));
        memset(img->fadjust4x4Cr[0][1][j], 0, MB_BLOCK_SIZE * sizeof(int));
      }
    }
    
#ifdef MV_COMPETITION
    if ((input->mv_competition > 0) && (mode == 0) && (img->type==P_SLICE))
    {
      // The last iteration consist of re-launching with the bestmode,
      // to set needed variables for the JM
      if (iter == (nbiter - 1))
        mv_comp.current_predictor_for_skip = currMB->best_predictor_for_skip;
      else
      {
        mv_comp.current_predictor_for_skip = iter;
        if ((mv_comp.mv_pred_skip[mv_comp.current_predictor_for_skip][0] == NOT_AVAILABLE) || (mv_comp.mv_pred_skip[mv_comp.current_predictor_for_skip][1] == NOT_AVAILABLE))
        {
          continue;
        }
      }
      
      LumaResidualCoding ();
      
      if ((!(img->residue_transform_flag && (mode==I4MB || mode==I16MB || mode==I8MB))) && img->yuv_format!=YUV400)
        ChromaResidualCoding (&dummy);
    }
    
    else 
    {
#endif
      if (mode<P8x8)
      {
        LumaResidualCoding ();
        
        if(mode==0 && currMB->cbp!=0 && (img->type != B_SLICE || img->NoResidueDirect==1)) 
          return 0;
        if(mode==0 && currMB->cbp==0 && currMB->luma_transform_size_8x8_flag==1) //for B_skip, luma_transform_size_8x8_flag=0 only
          return 0;
      }
      else if (mode==P8x8)
      {
        SetCoeffAndReconstruction8x8 (currMB);
      }
      else if (mode==I4MB)
      {
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
        if(input->UseRDO_Q || img->slice_fractional_quant_flag)
#else
        if(input->UseRDO_Q)
#endif
        {
#ifdef MB32X32
          if(Motion_Selected && IS_INTER_SLICE(img->type) && input->UseExtMB==0)
#else
          if(Motion_Selected && IS_INTER_SLICE(img->type))
#endif
          { 
            set_stored_macroblock_parameters_intra(mode); 
          }
          else
          {
            if (currMB->c_ipred_mode == DC_PRED_8)
              currMB->cbp = Mode_Decision_for_Intra4x4Macroblock (lambda, &dummy);
            
            set_I4MB_parameters();

            // Residue Color Transform
            if(img->residue_transform_flag)
            {
              for(i = 0; i<2; i++)
              {
                for(j = 0; j<4; j++)
                  for(k = 0; k<4; k++)
                    if(cbp_chroma_block[i][j][k])
                      cr_cbp = 2;
              }
              for(uv=0; uv<2; uv++)
                cr_cbp = dct_chroma_DC(uv, cr_cbp);
              
              currMB->cbp += (cr_cbp<<4);
            }
          }
        }
        else
        {
          currMB->cbp = Mode_Decision_for_Intra4x4Macroblock (lambda, &dummy);
          // Residue Color Transform
          if(img->residue_transform_flag)
          {
            for(i = 0; i<2; i++)
            {
              for(j = 0; j<4; j++)
                for(k = 0; k<4; k++)
                  if(cbp_chroma_block[i][j][k])
                    cr_cbp = 2;
            }
            for(uv=0; uv<2; uv++)
              cr_cbp = dct_chroma_DC(uv, cr_cbp);
            
            currMB->cbp += (cr_cbp<<4);
          }
        }
#else
        currMB->cbp = Mode_Decision_for_Intra4x4Macroblock (lambda, &dummy);
        // Residue Color Transform
        if(img->residue_transform_flag)
        {
          for(i = 0; i<2; i++)
          {
            for(j = 0; j<4; j++)
              for(k = 0; k<4; k++)
                if(cbp_chroma_block[i][j][k])
                  cr_cbp = 2;
          }
          for(uv=0; uv<2; uv++)
            cr_cbp = dct_chroma_DC(uv, cr_cbp);
          
          currMB->cbp += (cr_cbp<<4);
        }
#endif
        
        
#ifdef USE_INTRA_MDDT
        if(input->UseIntraMDDT)
        {
          for (i=0; i<16; i++)
          {
            if (img->quant_stat[i]>0)
            {
              img->qp_shift[i]+=(img->delta_shift);
              img->quant_stat[i]=(img->quant_stat[i]+1)/2;
            }
            else
            {
              img->qp_shift[i]-=(img->delta_shift);
              img->quant_stat[i]=(img->quant_stat[i]-1)/2;
            }
            img->qp_shift[i] = img->qp_shift[i]<0? 0 : img->qp_shift[i];
            img->qp_shift[i] = img->qp_shift[i]>img->QPVal? img->QPVal : img->qp_shift[i]; 
          }
        }
#endif
      }
      else if (mode == I16MB)
      {
#ifdef RDO_Q    
#ifdef ADAPTIVE_QUANTIZATION
        if(input->UseRDO_Q || img->slice_fractional_quant_flag) 
#else
        if(input->UseRDO_Q)
#endif
        {          
#ifdef MB32X32
          if(Motion_Selected && IS_INTER_SLICE(img->type) && input->UseExtMB==0)
#else
          if(Motion_Selected && IS_INTER_SLICE(img->type))
#endif
          {        
            set_stored_macroblock_parameters_intra(mode); 
          }
          else
          {          
              if (currMB->c_ipred_mode == DC_PRED_8)
              {
#ifdef USE_INTRA_MDDT         
                if(input->symbol_mode == CABAC && input->UseIntraMDDT)
                {                
                  Intra16x16_Mode_Decision_RDopt(currMB, &i16mode, lambda);
                }                 
                else
#endif
                {                
                  Intra16x16_Mode_Decision(currMB, &i16mode);
                }
              }
              
              set_I16MB_parameters(&i16mode);
              
          }
        }
        else
        {
          
#ifdef USE_INTRA_MDDT
          if(input->symbol_mode == CABAC && input->UseIntraMDDT)
          {
            Intra16x16_Mode_Decision_RDopt(currMB, &i16mode, lambda);
          }
          else
#endif
          {
            Intra16x16_Mode_Decision(currMB, &i16mode);
          }
        }
#else
#ifdef USE_INTRA_MDDT
        if(input->symbol_mode == CABAC && input->UseIntraMDDT)
        {
          Intra16x16_Mode_Decision_RDopt(currMB, &i16mode, lambda);
        }
        else
#endif
        {
          Intra16x16_Mode_Decision(currMB, &i16mode);
        }
#endif
        
#ifdef USE_INTRA_MDDT
        if(input->symbol_mode == CABAC && input->UseIntraMDDT)
        {
          for (i=0; i<256; i++)
          {
            if (img->quant_stat16x16[i]>0)
            {
              img->qp_shift16x16[i]+=(img->delta_shift16x16);
              img->quant_stat16x16[i]=(img->quant_stat16x16[i]+1)/2;
            }
            else
            {
              img->qp_shift16x16[i]-=(img->delta_shift16x16);
              img->quant_stat16x16[i]=(img->quant_stat16x16[i]-1)/2;
            }
            img->qp_shift16x16[i] = img->qp_shift16x16[i]<0? 0 : img->qp_shift16x16[i];
            img->qp_shift16x16[i] = img->qp_shift16x16[i]>img->QPVal16x16? img->QPVal16x16 : img->qp_shift16x16[i]; 
          }
        }
#endif
      }
      else if(mode==I8MB)
      {
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
        if(input->UseRDO_Q || img->slice_fractional_quant_flag)
#else
        if(input->UseRDO_Q)
#endif
        {
#ifdef MB32X32
          if(Motion_Selected && IS_INTER_SLICE(img->type) && input->UseExtMB==0)
#else
          if(Motion_Selected && IS_INTER_SLICE(img->type))
#endif

          {             
            set_stored_macroblock_parameters_intra(mode); 
          }          
          else
          { 
            if (currMB->c_ipred_mode == DC_PRED_8)
              currMB->cbp = Mode_Decision_for_new_Intra8x8Macroblock(lambda, &dummy);    
            
            set_I8MB_parameters();
          }
          // Residue Color Transform
          if(img->residue_transform_flag)
          {
            for(i = 0; i<2; i++)
            {
              for(j = 0; j<4; j++)
                for(k = 0; k<4; k++)
                  if(cbp_chroma_block[i][j][k])
                    cr_cbp = 2;
            }     
            for(uv = 0; uv<2; uv++)
              cr_cbp = dct_chroma_DC(uv, cr_cbp);
            
            currMB->cbp += (cr_cbp<<4);
          }
        }
        else
        {
          currMB->cbp = Mode_Decision_for_new_Intra8x8Macroblock(lambda, &dummy);
          // Residue Color Transform
          if(img->residue_transform_flag)
          {
            for(i = 0; i<2; i++)
            {
              for(j = 0; j<4; j++)
                for(k = 0; k<4; k++)
                  if(cbp_chroma_block[i][j][k])
                    cr_cbp = 2;
            }     
            for(uv = 0; uv<2; uv++)
              cr_cbp = dct_chroma_DC(uv, cr_cbp);
            
            currMB->cbp += (cr_cbp<<4);
          }
        }
#else
        currMB->cbp = Mode_Decision_for_new_Intra8x8Macroblock(lambda, &dummy);
        // Residue Color Transform
        if(img->residue_transform_flag)
        {
          for(i = 0; i<2; i++)
          {
            for(j = 0; j<4; j++)
              for(k = 0; k<4; k++)
                if(cbp_chroma_block[i][j][k])
                  cr_cbp = 2;
          }     
          for(uv = 0; uv<2; uv++)
            cr_cbp = dct_chroma_DC(uv, cr_cbp);
          
          currMB->cbp += (cr_cbp<<4);
        }
#endif
        
#ifdef USE_INTRA_MDDT
        if(input->UseIntraMDDT)
        {
          for (i=0; i<64; i++)
          {
            if (img->quant_stat8x8[i]>0)
            {
              img->qp_shift8x8[i]+=(img->delta_shift8x8);
              img->quant_stat8x8[i]=(img->quant_stat8x8[i]+1)/2;
            }
            else
            {
              img->qp_shift8x8[i]-=(img->delta_shift8x8);
              img->quant_stat8x8[i]=(img->quant_stat8x8[i]-1)/2;
            }
            img->qp_shift8x8[i] = img->qp_shift8x8[i]<0? 0 : img->qp_shift8x8[i];
            img->qp_shift8x8[i] = img->qp_shift8x8[i]>img->QPVal8x8? img->QPVal8x8 : img->qp_shift8x8[i]; 
          }
        }
#endif
      }
      else if(mode==IPCM)
      {
        for (j = 0; j < MB_BLOCK_SIZE; j++)
        {
          j1 = j + img->opix_y;
          j2 = j + img->pix_y;
          for (i=img->opix_x; i<img->opix_x+MB_BLOCK_SIZE; i++)        
            enc_picture->imgY[j2][i] = imgY_org[j1][i];
        }
        if (img->yuv_format != YUV400)
        {
          // CHROMA
          for (j = 0; j<img->mb_cr_size_y; j++)
          {
            j1 = j + img->opix_c_y;
            j2 = j + img->pix_c_y;
            for (i=img->opix_c_x; i<img->opix_c_x+img->mb_cr_size_x; i++)
            {
              enc_picture->imgUV[0][j2][i] = imgUV_org[0][j1][i];
              enc_picture->imgUV[1][j2][i] = imgUV_org[1][j1][i];
            }
          }
        }  
        for (j=0;j<4;j++)
          for (i=0; i<(4+img->num_blk8x8_uv); i++)
            img->nz_coeff[img->current_mb_nr][j][i] = 16;
          
      }
      
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
      if(input->UseRDO_Q || img->slice_fractional_quant_flag)
#else
      if(input->UseRDO_Q)
#endif
      {
#ifdef MB32X32
        if (img->type == I_SLICE  || !IS_INTRA_MODE(mode) || (IS_INTRA_MODE(mode) && (!Motion_Selected || input->UseExtMB!= 0) && IS_INTER_SLICE(img->type)))
#else
        if (img->type == I_SLICE  || !IS_INTRA_MODE(mode) || (IS_INTRA_MODE(mode) && !Motion_Selected && IS_INTER_SLICE(img->type)))
#endif
        { 
          if (input->rdopt==3 && img->type!=B_SLICE)
          {
            // We need the reconstructed prediction residue for the simulated decoders.
            compute_residue_mb (mode==I16MB?i16mode:-1);
          }
          
          //Rate control
          if (input->RCEnable)
          {
            if (mode == I16MB)
              memcpy(pred,img->mprr_2[i16mode],MB_PIXELS * sizeof(imgpel));
            else
              memcpy(pred,img->mpr,MB_PIXELS * sizeof(imgpel));
          }
          
          img->i16offset = 0;
          dummy = 0;
          if ((!(img->residue_transform_flag && (mode==I4MB || mode==I16MB || mode==I8MB)))
            && (img->yuv_format!=YUV400) && (mode != IPCM))
            ChromaResidualCoding (&dummy);
          
          if (mode==I16MB)     
            img->i16offset = I16Offset  (currMB->cbp, i16mode);
          
#ifdef MB32X32
          if(IS_INTER_SLICE(img->type) && IS_INTRA_MODE(mode) && !Motion_Selected && input->UseExtMB==0)
#else
          if(IS_INTER_SLICE(img->type) && IS_INTRA_MODE(mode) && !Motion_Selected)
#endif
            store_macroblock_parameters_intra(mode); 
        }
      }
      else
      {
        if (input->rdopt==3 && img->type!=B_SLICE)
        {
          // We need the reconstructed prediction residue for the simulated decoders.
          compute_residue_mb (mode==I16MB?i16mode:-1);
        }
        
        //Rate control
        if (input->RCEnable)
        {
          if (mode == I16MB)
            memcpy(pred,img->mprr_2[i16mode],MB_PIXELS * sizeof(imgpel));
          else
            memcpy(pred,img->mpr,MB_PIXELS * sizeof(imgpel));
        }
        
        img->i16offset = 0;
        dummy = 0;
        if ((!(img->residue_transform_flag && (mode==I4MB || mode==I16MB || mode==I8MB)))
          && (img->yuv_format!=YUV400) && (mode != IPCM))
          ChromaResidualCoding (&dummy);
        
        if (mode==I16MB)     
          img->i16offset = I16Offset  (currMB->cbp, i16mode);
      }
#else
      if (input->rdopt==3 && img->type!=B_SLICE)
      {
        // We need the reconstructed prediction residue for the simulated decoders.
        compute_residue_mb (mode==I16MB?i16mode:-1);
      }
      
      //Rate control
      if (input->RCEnable)
      {
        if (mode == I16MB)
          memcpy(pred,img->mprr_2[i16mode],MB_PIXELS * sizeof(imgpel));
        else
          memcpy(pred,img->mpr,MB_PIXELS * sizeof(imgpel));
      }
      
      img->i16offset = 0;
      dummy = 0;
      if ((!(img->residue_transform_flag && (mode==I4MB || mode==I16MB || mode==I8MB)))
        && (img->yuv_format!=YUV400) && (mode != IPCM))
        ChromaResidualCoding (&dummy);
      
      if (mode==I16MB)     
        img->i16offset = I16Offset  (currMB->cbp, i16mode);
      
#endif
      
#ifdef MV_COMPETITION
    }
#endif
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
#ifdef MB32X32
    if( (input->UseRDO_Q || img->slice_fractional_quant_flag) && IS_INTRA_MODE(mode) && Motion_Selected && IS_INTER_SLICE(img->type) && input->UseExtMB==0)
#else
    if( (input->UseRDO_Q || img->slice_fractional_quant_flag) && IS_INTRA_MODE(mode) && Motion_Selected && IS_INTER_SLICE(img->type))
#endif
#else
    if(input->UseRDO_Q && IS_INTRA_MODE(mode) && Motion_Selected && IS_INTER_SLICE(img->type))
#endif
    {
      int         intra_mode_mapping_table[] = 
      {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, -1, -1, 2, 3};
      int         luma_mode = intra_mode_mapping_table[mode];
      rdcost = rdcost_intra[luma_mode][currMB->c_ipred_mode];
    }
    else
    {
#endif
      
      //=====
      //=====   GET DISTORTION
      //=====
      // LUMA
      if (input->rdopt==3 && img->type!=B_SLICE)
      {
        for (k = 0; k<input->NoOfDecoders ;k++)
        {
          decode_one_mb (k, currMB);
          for (j = 0; j<MB_BLOCK_SIZE; j++)
          {
            for (i=img->opix_x; i<img->opix_x+MB_BLOCK_SIZE; i++)
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
            {
              distortion += (img->BitDepthIncrease)? img->quad [(imgY_org[img->opix_y+j][i]>>img->BitDepthIncrease)-Clip3(0,((1<<input->InputBitDepth)-1),(decs->decY[k][img->opix_y+j][i]+(1<<(img->BitDepthIncrease-1)))>>img->BitDepthIncrease)]
                : img->quad [imgY_org[img->opix_y+j][i] - decs->decY[k][img->opix_y+j][i]];
            }
#else
            distortion += img->quad [imgY_org[img->opix_y+j][i] - decs->decY[k][img->opix_y+j][i]];
#endif
          }
        }
        distortion /= input->NoOfDecoders;
      }
      else
      {
        for (j = 0; j < MB_BLOCK_SIZE; j++)
        {
          j1 = j + img->opix_y;
          j2 = j + img->pix_y;
          for (i=img->opix_x; i<img->opix_x+MB_BLOCK_SIZE; i++)        
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
            distortion += SQR_DEPTH(imgY_org[j1][i], enc_picture->imgY[j2][i], input->BitDepthLuma, img->BitDepthIncrease);
#else
          distortion += img->quad [imgY_org[j1][i] - enc_picture->imgY[j2][i]];
#endif
        }
      }
      
      if (img->yuv_format != YUV400)
      {
        // CHROMA
        for (j = 0; j<img->mb_cr_size_y; j++)
        {
          j1 = j + img->opix_c_y;
          j2 = j + img->pix_c_y;
          for (i=img->opix_c_x; i<img->opix_c_x+img->mb_cr_size_x; i++)
          {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
            distortion += SQR_DEPTH(imgUV_org[0][j1][i], enc_picture->imgUV[0][j2][i], input->BitDepthChroma, img->BitDepthIncreaseChroma);
            distortion += SQR_DEPTH(imgUV_org[1][j1][i], enc_picture->imgUV[1][j2][i], input->BitDepthChroma, img->BitDepthIncreaseChroma);
#else
            distortion += img->quad [imgUV_org[0][j1][i] - enc_picture->imgUV[0][j2][i]];
            distortion += img->quad [imgUV_org[1][j1][i] - enc_picture->imgUV[1][j2][i]];
#endif
          }
        }
      }  
      
      //=====   S T O R E   C O D I N G   S T A T E   =====
      //---------------------------------------------------
      store_coding_state (cs_cm);
      
      
      //=====
      //=====   GET RATE
      //=====
      //----- macroblock header -----
      if (use_of_cc)
      {
        if (currMB->mb_type!=0 || (bframe && currMB->cbp!=0))
        {
          // cod counter and macroblock mode are written ==> do not consider code counter
          tmp_cc = img->cod_counter;
          rate   = writeMBLayer (1, &coeff_rate);
          ue_linfo (tmp_cc, dummy, &cc_rate, &dummy);
          rate  -= cc_rate;
          img->cod_counter = tmp_cc;
        }
        else
        {
          // cod counter is just increased  ==> get additional rate
          ue_linfo (img->cod_counter+1, dummy, &rate,    &dummy);
          ue_linfo (img->cod_counter,   dummy, &cc_rate, &dummy);
          rate -= cc_rate;
        }
      }
      else
      {
        rate = writeMBLayer (1, &coeff_rate);
      }
      
      //=====   R E S T O R E   C O D I N G   S T A T E   =====
      //-------------------------------------------------------
      reset_coding_state (cs_cm);
      
      rdcost = (double)distortion + lambda * max(0.5,(double)rate);
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
#ifdef MB32X32
      if(input->UseExtMB==0 && (input->UseRDO_Q || img->slice_fractional_quant_flag) && IS_INTRA_MODE(mode) && IS_INTER_SLICE(img->type))
#else
      if( (input->UseRDO_Q || img->slice_fractional_quant_flag) && IS_INTRA_MODE(mode) && IS_INTER_SLICE(img->type))
#endif
#else
      if(input->UseRDO_Q && IS_INTRA_MODE(mode) && IS_INTER_SLICE(img->type))
#endif
      {
        int         intra_mode_mapping_table[] = 
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, -1, -1, 2, 3};
        int         luma_mode = intra_mode_mapping_table[mode];
        
        rdcost_intra[luma_mode][currMB->c_ipred_mode] = rdcost; 
      }
    }
#endif
    
    if (!(rdcost >= *min_rdcost ||
      ((img->qp + img->bitdepth_luma_qp_scale)==0 && img->lossless_qpprime_flag==1 && distortion!=0)))
    {                
      if (((img->MbaffFrameFlag) && (mode ? 0: ((img->type == B_SLICE) ? !currMB->cbp:1)))  // AFF and current is skip
        && (img->current_mb_nr & 0x01) //bottom
        && (prevMB->mb_type ? 0:((img->type == B_SLICE) ? !prevMB->cbp:1)) //top is skip
        && (!(field_flag_inference() == currMB->mb_field)))
        return_val = 0;
      else      
      {
        //=====   U P D A T E   M I N I M U M   C O S T   =====
        //-----------------------------------------------------
        *min_rdcost = rdcost;
        *min_rate = lambda * (double)coeff_rate;
        
#ifdef BEST_NZ_COEFF
        for (j=0;j<4;j++)
          for (i=0; i<(4+img->num_blk8x8_uv); i++)
            gaaiMBAFF_NZCoeff[j][i] = img->nz_coeff[img->current_mb_nr][j][i]; 
#endif
          
#ifdef MV_COMPETITION      
          if ((input->mv_competition > 0) && (mode == 0) && (img->type==P_SLICE))
          {
            int by, bx;
            
            currMB->best_predictor_for_skip = mv_comp.current_predictor_for_skip;
            
            for (by = 0;by < 4;by++)
            {
              for (bx = 0;bx < 4;bx++)
              {
                all_mv [by][bx][0][0][0][0] = mv_comp.mv_pred_skip[currMB->best_predictor_for_skip][0];
                all_mv [by][bx][0][0][0][1] = mv_comp.mv_pred_skip[currMB->best_predictor_for_skip][1];
              }
            }
          }
#endif
          return_val = 1;
      }
      
      
    }
    else 
    {
#ifdef MV_COMPETITION
      if (!((input->mv_competition > 0) && (mode == 0) && (img->type==P_SLICE)))
      {
#endif
        
#if FASTMODE
        // Reordering RDCost comparison order of mode 0 and mode 1 in P_SLICE
        // if RDcost of mode 0 and mode 1 is same, we choose best_mode is 0
        // This might not always be good since mode 0 is more biased towards rate than quality.
        if((img->type!=P_SLICE || mode != 0 || rdcost != *min_rdcost) || input->ProfileIDC>=FREXT_HP)
#endif
          return_val = 0;
#ifdef MV_COMPETITION
      }
#endif
    }
#ifdef MV_COMPETITION
  }  //For
#endif  
  
  return (return_val);
}

/*!
*************************************************************************************
* \brief
*    Store adaptive rounding parameters
*************************************************************************************
*/
void store_adaptive_rounding_parameters (int mode, Macroblock *currMB)
{
  int i,j;
#ifdef ADAPTIVE_FD_SD_CODING
  best_adjust_adaptive_f_spatial_domain_4x4=0.0;
  best_adjust_adaptive_f_spatial_domain_8x8=0.0;
#endif
  
  for (j = 0; j < MB_BLOCK_SIZE; j++)
  {
    for (i = 0; i < MB_BLOCK_SIZE; i++)
    {      
      if ((mode == P8x8))
      {
        if (currMB->luma_transform_size_8x8_flag)
        {
          bestInterFAdjust8x8[j][i]=img->fadjust8x8[2][j][i]; 
#ifdef ADAPTIVE_FD_SD_CODING
          best_adjust_adaptive_f_spatial_domain_8x8=img->adjust_adaptive_f_spatial_domain_8x8;
#endif
          if ((i >> 3) ==0 && (j >> 3) == 0)
          {
            bestInterFAdjust4x4Cr[0][j][i]=img->fadjust8x8Cr[0][0][j][i]; 
            bestInterFAdjust4x4Cr[1][j][i]=img->fadjust8x8Cr[0][1][j][i]; 
          }
        }
        else
        {
          bestInterFAdjust4x4[j][i]=img->fadjust4x4[3][j][i]; 
#ifdef ADAPTIVE_FD_SD_CODING
          best_adjust_adaptive_f_spatial_domain_4x4=img->adjust_adaptive_f_spatial_domain_4x4;
#endif
          if ((i >> 3) ==0 && (j >> 3) == 0)
          {
            bestInterFAdjust4x4Cr[0][j][i]=img->fadjust4x4Cr[2][0][j][i]; 
            bestInterFAdjust4x4Cr[1][j][i]=img->fadjust4x4Cr[2][1][j][i]; 
          }
        }
      }
      else if ((mode != I4MB)&&(mode != I16MB)&&(mode != I8MB))
      {
        if (currMB->luma_transform_size_8x8_flag)
        {
          bestInterFAdjust8x8[j][i]=img->fadjust8x8[0][j][i]; 
#ifdef ADAPTIVE_FD_SD_CODING
          best_adjust_adaptive_f_spatial_domain_8x8=img->adjust_adaptive_f_spatial_domain_8x8;
#endif
        }
        else
        {
          bestInterFAdjust4x4[j][i]=img->fadjust4x4[0][j][i]; 
#ifdef ADAPTIVE_FD_SD_CODING
          best_adjust_adaptive_f_spatial_domain_4x4=img->adjust_adaptive_f_spatial_domain_4x4;
#endif
        }        
        if ((i >> 3) ==0 && (j >> 3) == 0)
        {
          bestInterFAdjust4x4Cr[0][j][i]=img->fadjust4x4Cr[0][0][j][i]; 
          bestInterFAdjust4x4Cr[1][j][i]=img->fadjust4x4Cr[0][1][j][i]; 
        }        
      }
      else if (mode != I8MB)
      {
        bestIntraFAdjust4x4[j][i]=img->fadjust4x4[1 + mode == I16MB][j][i]; 
        
        if ((i >> 3) ==0 && (j >> 3) == 0)
        {
          bestIntraFAdjust4x4Cr[0][j][i]=img->fadjust4x4Cr[1][0][j][i]; 
          bestIntraFAdjust4x4Cr[1][j][i]=img->fadjust4x4Cr[1][1][j][i]; 
        } 
      }
      else
      {
        bestIntraFAdjust8x8[j][i]=img->fadjust8x8[1][j][i]; 
        
        if ((i >> 3) ==0 && (j >> 3) == 0)
        {
          bestIntraFAdjust4x4Cr[0][j][i]=img->fadjust4x4Cr[0][1][j][i]; 
          bestIntraFAdjust4x4Cr[1][j][i]=img->fadjust4x4Cr[1][1][j][i]; 
        } 
      }
    }
  }
}

#ifdef RDO_Q
/*!
****************************************************************************
* \brief
*    save encoding information for intra MB
****************************************************************************
*/
void store_macroblock_parameters_intra (int mode)
{
  int        i, j, k;
  int        num_blk8x8 = 4 + img->num_blk8x8_uv;
  Macroblock *currMB  = &img->mb_data[img->current_mb_nr];
  int        intra_mode_mapping_table[] = 
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, -1, -1, 2, 3};
  int        luma_mode = intra_mode_mapping_table[mode];
  int        chroma_mode = currMB->c_ipred_mode; 
  
  assert(luma_mode != -1);
  
  if(currMB->mb_type == I8MB)
  {
    for(j = 0; j < BLOCK_MULTIPLE; j++)
    {
      //memcpy(b8_ipredmode8x8_intra[luma_mode][chroma_mode][j], &img->ipredmode[img->block_y+j][img->block_x],BLOCK_MULTIPLE * sizeof(char));
      memcpy(b8_ipredmode8x8_intra[luma_mode][chroma_mode][j], &img->ipredmode8x8[img->block_y+j][img->block_x], BLOCK_MULTIPLE * sizeof(char));
    }
    memcpy(b8_ipredmode8x8_mb_intra[luma_mode][chroma_mode], currMB->intra_pred_modes8x8, BLOCK_MULTIPLE*BLOCK_MULTIPLE*sizeof(char));
  }
  else if (currMB->mb_type == I4MB)
  {
    for(j = 0; j < BLOCK_MULTIPLE; j++)
      memcpy(&b4_ipredmode_intra[luma_mode][chroma_mode][BLOCK_MULTIPLE * j], &img->ipredmode[img->block_y + j][img->block_x],BLOCK_MULTIPLE * sizeof(char));
    memcpy(b4_ipredmode_mb_intra[luma_mode][chroma_mode], currMB->intra_pred_modes, BLOCK_MULTIPLE*BLOCK_MULTIPLE*sizeof(char));
  }
  
  //--- reconstructed blocks ----
  for (j = 0; j < MB_BLOCK_SIZE; j++)
  {
    memcpy(rec_mbY_intra[luma_mode][chroma_mode][j],&enc_picture->imgY[img->pix_y+j][img->pix_x], MB_BLOCK_SIZE * sizeof(imgpel));
  }
  if (img->yuv_format != YUV400)
  {
    for (j = 0; j<img->mb_cr_size_y; j++)
    {
      memcpy(rec_mbU_intra[luma_mode][chroma_mode][j],&enc_picture->imgUV[0][img->pix_c_y+j][img->pix_c_x], img->mb_cr_size_x * sizeof(imgpel));
      memcpy(rec_mbV_intra[luma_mode][chroma_mode][j],&enc_picture->imgUV[1][img->pix_c_y+j][img->pix_c_x], img->mb_cr_size_x * sizeof(imgpel));
    }
  }
  
  
  //--- coeff, cbp, kac ---
  for(i = 0; i < num_blk8x8; i++)
    for(j = 0; j < 4; j++)
      for(k = 0; k < 2; k++)
        memcpy(cofAC_intra[luma_mode][chroma_mode][i][j][k], img->cofAC[i][j][k], sizeof(int)*65);
      for(i = 0; i < 3; i++)
        for(j = 0; j < 2; j++)
          memcpy(cofDC_intra[luma_mode][chroma_mode][i][j], img->cofDC[i][j], sizeof(int)*65);
        cbp_intra[luma_mode][chroma_mode]     = currMB->cbp;
        cbp_blk_intra[luma_mode][chroma_mode] = currMB->cbp_blk;
        i16offset_intra[luma_mode][chroma_mode] = img->i16offset;
}

/*!
****************************************************************************
* \brief
*    assign saved encoding information for intra MB
****************************************************************************
*/
void set_stored_macroblock_parameters_intra (int mode)
{
  int  i, j, k;
  int num_blk8x8 = 4 + img->num_blk8x8_uv;
  Macroblock  *currMB  = &img->mb_data[img->current_mb_nr];
  int         intra_mode_mapping_table[] = 
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, -1, -1, 2, 3};
  int         luma_mode = intra_mode_mapping_table[mode];
  int        chroma_mode = currMB->c_ipred_mode; 
  int         bframe   = (img->type==B_SLICE);  
  int        block_x, block_y;  
  
  assert(luma_mode != -1);
  
  //===== reconstruction values =====
  for (j = 0; j < MB_BLOCK_SIZE; j++)
  {
    memcpy(&enc_picture->imgY[img->pix_y+j][img->pix_x],rec_mbY_intra[luma_mode][chroma_mode][j], MB_BLOCK_SIZE * sizeof(imgpel));
    if(img->MbaffFrameFlag)
      memcpy(rdopt->rec_mbY[j],rec_mbY_intra[luma_mode][chroma_mode][j], MB_BLOCK_SIZE * sizeof(imgpel));
  }
  if (img->yuv_format != YUV400)
  {
    for (j = 0; j<img->mb_cr_size_y; j++)
    {
      memcpy(&enc_picture->imgUV[0][img->pix_c_y+j][img->pix_c_x],rec_mbU_intra[luma_mode][chroma_mode][j], img->mb_cr_size_x * sizeof(imgpel));
      memcpy(&enc_picture->imgUV[1][img->pix_c_y+j][img->pix_c_x],rec_mbV_intra[luma_mode][chroma_mode][j], img->mb_cr_size_x * sizeof(imgpel));
      if(img->MbaffFrameFlag)
      {
        memcpy(rdopt->rec_mbU[j],rec_mbU_intra[luma_mode][chroma_mode][j], img->mb_cr_size_x * sizeof(imgpel));
        memcpy(rdopt->rec_mbV[j],rec_mbV_intra[luma_mode][chroma_mode][j], img->mb_cr_size_x * sizeof(imgpel));
      }
    }
  }
  
  //===== coefficients and cbp =====
  for(i = 0; i < num_blk8x8; i++)
    for(j = 0; j < 4; j++)
      for(k = 0; k < 2; k++)
        memcpy(img->cofAC[i][j][k], cofAC_intra[luma_mode][chroma_mode][i][j][k], sizeof(int)*65);
      for(i = 0; i < 3; i++)
        for(j = 0; j < 2; j++)
          memcpy(img->cofDC[i][j], cofDC_intra[luma_mode][chroma_mode][i][j], sizeof(int)*65);
        currMB->cbp      = cbp_intra[luma_mode][chroma_mode];
        currMB->cbp_blk  = cbp_blk_intra[luma_mode][chroma_mode];
        //==== macroblock type ====
        currMB->mb_type = mode;
        img->i16offset = i16offset_intra[luma_mode][chroma_mode];
        
        //==== transform size flag ====
        currMB->luma_transform_size_8x8_flag = (currMB->mb_type == I8MB);
        
        rdopt->luma_transform_size_8x8_flag  = currMB->luma_transform_size_8x8_flag;
        
        if (input->rdopt==3 && img->type!=B_SLICE)
        {
          //! save the MB Mode of every macroblock
          decs->dec_mb_mode[img->mb_x][img->mb_y] = mode;
        }
        
        if (bframe)
        {
          for (j=0; j<4; j++)
          {
            block_y = img->block_y + j;
            for (i=0; i<4; i++)
            {          
              block_x = img->block_x + i;
              k = 2*(j >> 1)+(i >> 1);
              
              // forward
              enc_picture->ref_idx    [LIST_1][block_y][block_x]    = -1;
              enc_picture->ref_pic_id [LIST_1][block_y][block_x]    = -1;
              enc_picture->mv         [LIST_1][block_y][block_x][0] = 0;
              enc_picture->mv         [LIST_1][block_y][block_x][1] = 0;
              if(img->MbaffFrameFlag)
                rdopt->refar[LIST_1][j][i] = -1;
            }
          }
        }
        
        if(currMB->mb_type == I8MB)
        {
          for(j = 0; j < BLOCK_MULTIPLE; j++)
          {
            memcpy(&img->ipredmode[img->block_y+j][img->block_x],b8_ipredmode8x8_intra[luma_mode][chroma_mode][j], BLOCK_MULTIPLE * sizeof(char));
            memcpy(&img->ipredmode8x8[img->block_y+j][img->block_x], b8_ipredmode8x8_intra[luma_mode][chroma_mode][j], BLOCK_MULTIPLE * sizeof(char));
          }
          memcpy(currMB->intra_pred_modes8x8, b8_ipredmode8x8_mb_intra[luma_mode][chroma_mode], BLOCK_MULTIPLE*BLOCK_MULTIPLE*sizeof(char));
        }
        else if (currMB->mb_type == I4MB)
        {
          for(j = 0; j < BLOCK_MULTIPLE; j++)
            memcpy(&img->ipredmode[img->block_y + j][img->block_x],&b4_ipredmode_intra[luma_mode][chroma_mode][BLOCK_MULTIPLE * j], BLOCK_MULTIPLE * sizeof(char));
          memcpy(currMB->intra_pred_modes, b4_ipredmode_mb_intra[luma_mode][chroma_mode], BLOCK_MULTIPLE*BLOCK_MULTIPLE*sizeof(char));
        }
}

/*!
****************************************************************************
* \brief
*    save/assign encoding information for Intra4x4 MB
****************************************************************************
*/
void set_I4MB_parameters()
{
  Macroblock  *currMB   = &img->mb_data[img->current_mb_nr];
  int *src1, *dest1;
  imgpel *src2, *dest2, *src3, *dest3;
  int b8, b4,j, y;
  
  for (b8=0; b8<4; b8++)
  {
    for (b4=0; b4<4; b4++)
    {
      for (j=0; j<2; j++)    
      {
        if (currMB->c_ipred_mode == DC_PRED_8)
        {
          src1 = img->cofAC[b8][b4][j];
          dest1 = saved_cofAC[b8][b4][j];
        }
        else
        {
          dest1 = img->cofAC[b8][b4][j];
          src1 = saved_cofAC[b8][b4][j];
        }
        memcpy (dest1, src1, 18 * sizeof(int));
      }
    }
  }
  
  for (y=0; y<MB_BLOCK_SIZE; y++)
  {                
    if (currMB->c_ipred_mode == DC_PRED_8)
    {
      src2 = &enc_picture->imgY[img->pix_y+y][img->pix_x];
      dest2 = &saved_enc_picture_imgY[y][0];
      src3 = &img->mpr[y][0];
      dest3 = &saved_img_mpr[y][0];
    }
    else
    {
      dest2 = &enc_picture->imgY[img->pix_y+y][img->pix_x];
      src2 = &saved_enc_picture_imgY[y][0];
      dest3 = &img->mpr[y][0];
      src3 = &saved_img_mpr[y][0];
    }
    
    memcpy (dest2, src2, MB_BLOCK_SIZE * sizeof(imgpel));                
    memcpy (dest3, src3, MB_BLOCK_SIZE * sizeof(imgpel)); 
  }
  
  if(currMB->c_ipred_mode == DC_PRED_8)
  {
    saved_cbp = currMB->cbp;
  }
  else
  {
    currMB->cbp = saved_cbp;
  }
}

/*!
****************************************************************************
* \brief
*    save/assign encoding information for Intra16x16 MB
****************************************************************************
*/
void set_I16MB_parameters(int *i16mode)
{
  Macroblock  *currMB   = &img->mb_data[img->current_mb_nr];
  int *src1, *dest1;
  imgpel *src2, *dest2;
  int b8, b4,j, y;
  
  for (b8=0; b8<4; b8++)
  {
    for (b4=0; b4<4; b4++)
    {
      for (j=0; j<2; j++)    
      {
        if (currMB->c_ipred_mode == DC_PRED_8)
        {
          src1 = img->cofAC[b8][b4][j];
          dest1 = saved_I16cofAC[b8][b4][j];
        }
        else
        {
          dest1 = img->cofAC[b8][b4][j];
          src1 = saved_I16cofAC[b8][b4][j];
        }
        memcpy (dest1, src1, 18 * sizeof(int));
      }
    }
  }
  
  for(j=0; j<2; j++)
  {
    if (currMB->c_ipred_mode == DC_PRED_8)
    {
      src1 = img->cofDC[0][j];
      dest1 = saved_cofDC[0][j];
    }
    else
    {
      dest1 = img->cofDC[0][j];
      src1 = saved_cofDC[0][j];
    }
    memcpy(dest1, src1, 18 * sizeof(int));
  }
  
  for (y=0; y<MB_BLOCK_SIZE; y++)
  {                
    if (currMB->c_ipred_mode == DC_PRED_8)
    {
      src2 = &enc_picture->imgY[img->pix_y+y][img->pix_x];
      dest2 = &saved_I16enc_picture_imgY[y][0];
    }
    else
    {
      dest2 = &enc_picture->imgY[img->pix_y+y][img->pix_x];
      src2 = &saved_I16enc_picture_imgY[y][0];
    }
    memcpy (dest2, src2, MB_BLOCK_SIZE * sizeof(imgpel));                
  }
  
  if(currMB->c_ipred_mode == DC_PRED_8)
  {
    saved_I16cbp = currMB->cbp;
    saved_i16mode = *i16mode;
  }
  else
  {
    currMB->cbp = saved_I16cbp;
    *i16mode = saved_i16mode;
  }
}

/*!
****************************************************************************
* \brief
*    save/assign encoding information for Intra8x8 MB
****************************************************************************
*/
void set_I8MB_parameters()
{
  Macroblock  *currMB   = &img->mb_data[img->current_mb_nr];
  int *src1, *dest1;
  imgpel *src2, *dest2, *src3, *dest3;
  int b8, b4,j, y;
  
  for (b8=0; b8<4; b8++)
  {
    for (b4=0; b4<4; b4++)
    {
      for (j=0; j<2; j++)    
      {
        if (currMB->c_ipred_mode == DC_PRED_8)
        {
          src1 = img->cofAC[b8][b4][j];
          dest1 = saved_I8cofAC[b8][b4][j];
        }    
        else
        {
          dest1 = img->cofAC[b8][b4][j];
          src1 = saved_I8cofAC[b8][b4][j];
        }
        memcpy (dest1, src1, 65 * sizeof(int));
      }
    }
  }
  
  for (y=0; y<MB_BLOCK_SIZE; y++)
  {                
    if (currMB->c_ipred_mode == DC_PRED_8)
    {
      src2 = &enc_picture->imgY[img->pix_y+y][img->pix_x];
      dest2 = &saved_I8enc_picture_imgY[y][0];
      src3 = &img->mpr[y][0];
      dest3 = &saved_I8img_mpr[y][0];
    }
    else
    {
      dest2 = &enc_picture->imgY[img->pix_y+y][img->pix_x];
      src2 = &saved_I8enc_picture_imgY[y][0];
      dest3 = &img->mpr[y][0];
      src3 = &saved_I8img_mpr[y][0];
    }
    memcpy (dest2, src2, MB_BLOCK_SIZE * sizeof(imgpel));                
    memcpy (dest3, src3, MB_BLOCK_SIZE * sizeof(imgpel)); 
  }
  
  if(currMB->c_ipred_mode == DC_PRED_8)
  {
    saved_I8cbp = currMB->cbp;
  }
  else
  {
    currMB->cbp = saved_I8cbp;
  }
}

#endif

/*!
*************************************************************************************
* \brief
*    Store macroblock parameters
*************************************************************************************
*/
void store_macroblock_parameters (int mode)
{
  int  i, j, k, ****i4p, ***i3p;
  Macroblock *currMB  = &img->mb_data[img->current_mb_nr];
  int        bframe   = (img->type==B_SLICE);
  
  //--- store best mode ---
  best_mode = mode;
  best_c_imode = currMB->c_ipred_mode;
  best_i16offset = img->i16offset;
  
  // If condition is not really necessary.
  bi_pred_me = (mode == 1) ? currMB->bi_pred_me : 0;  
  
  memcpy(b8mode, currMB->b8mode, BLOCK_MULTIPLE * sizeof(int));
  memcpy(b8pdir, currMB->b8pdir, BLOCK_MULTIPLE * sizeof(int));
  
  // Residue Color Transform
  //for (k = 0, j=img->block_y; j<img->block_y+BLOCK_MULTIPLE; j++, k+=BLOCK_MULTIPLE)
  memcpy(b4_intra_pred_modes,currMB->intra_pred_modes, MB_BLOCK_PARTITIONS * sizeof(char));
  memcpy(b8_intra_pred_modes8x8,currMB->intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));
  for (j = 0 ; j<BLOCK_MULTIPLE; j++)
  {
    memcpy(&b4_ipredmode[j * BLOCK_MULTIPLE],&img->ipredmode[img->block_y + j][img->block_x],BLOCK_MULTIPLE * sizeof(char));
    memcpy(b8_ipredmode8x8[j],&img->ipredmode8x8[img->block_y + j][img->block_x],BLOCK_MULTIPLE * sizeof(char));
  }  
  //--- reconstructed blocks ----
  for (j = 0; j < MB_BLOCK_SIZE; j++)
  {
    memcpy(rec_mbY[j],&enc_picture->imgY[img->pix_y+j][img->pix_x], MB_BLOCK_SIZE * sizeof(imgpel));
    if((img->type==SP_SLICE) && (si_frame_indicator==0 && sp2_frame_indicator==0))
      memcpy(lrec_rec[j],&lrec[img->pix_y+j][img->pix_x], MB_BLOCK_SIZE * sizeof(int));//store coefficients SP frame
  }
  
  if (img->AdaptiveRounding)
    store_adaptive_rounding_parameters (mode, currMB);
  
  if (img->yuv_format != YUV400)
  {
    for (j = 0; j<img->mb_cr_size_y; j++)
    {
      memcpy(rec_mbU[j],&enc_picture->imgUV[0][img->pix_c_y+j][img->pix_c_x], img->mb_cr_size_x * sizeof(imgpel));
      memcpy(rec_mbV[j],&enc_picture->imgUV[1][img->pix_c_y+j][img->pix_c_x], img->mb_cr_size_x * sizeof(imgpel));
      if((img->type==SP_SLICE) && (si_frame_indicator==0 && sp2_frame_indicator==0))
      {//store uv coefficients SP frame
        memcpy(lrec_rec_U[j],&lrec_uv[0][img->pix_c_y+j][img->pix_c_x], img->mb_cr_size_x * sizeof(int));
        memcpy(lrec_rec_V[j],&lrec_uv[1][img->pix_c_y+j][img->pix_c_x], img->mb_cr_size_x * sizeof(int));
      }
    }
  }
  
  
  //--- store results of decoders ---
  if (input->rdopt==3 && img->type!=B_SLICE)
  {
    for (k = 0; k<input->NoOfDecoders; k++)
    {
      for (j=img->pix_y; j<img->pix_y+16; j++)
        for (i=img->pix_x; i<img->pix_x+16; i++)
        {
          // Keep the decoded values of each MB for updating the ref frames
          decs->decY_best[k][j][i] = decs->decY[k][j][i];
        }
    }
  }
  
  //--- coeff, cbp, kac ---
  if (mode || bframe)
  {
    i4p=cofAC; cofAC=img->cofAC; img->cofAC=i4p;
    i3p=cofDC; cofDC=img->cofDC; img->cofDC=i3p;
    cbp     = currMB->cbp;
    cbp_blk = currMB->cbp_blk;
  }
  else
  {
    cbp_blk = cbp = 0;
  }
  
  //--- store transform size ---
  luma_transform_size_8x8_flag = currMB->luma_transform_size_8x8_flag;
  
  for (j = 0; j<4; j++)
    memcpy(frefframe[j],&enc_picture->ref_idx[LIST_0][img->block_y+j][img->block_x], BLOCK_MULTIPLE * sizeof(char));
  
  if (bframe)
  {
    for (j = 0; j<4; j++)
      memcpy(brefframe[j],&enc_picture->ref_idx[LIST_1][img->block_y+j][img->block_x], BLOCK_MULTIPLE * sizeof(char));
  }  
  
  
  
#ifdef ADAPTIVE_FD_SD_CODING
  best_SD_Coding_on_off=currMB->SD_Coding_on_off;
  memcpy(best_quantizer_indices,currMB->quantizer_indices,256*sizeof(int));
  memcpy(best_SD_or_FD,currMB->SD_or_FD,4*sizeof(int));
  best_SD_or_FD_t8x8=currMB->SD_or_FD_t8x8;
#endif
#ifdef ADAPTIVE_QUANTIZATION
  best_mb_iaqms_idx = currMB->mb_iaqms_idx;
#endif
#ifdef USE_INTRA_MDDT
  if(input->UseIntraMDDT)
  {
    memcpy(bestCofAC16[0], img->cofAC16[0], sizeof(int)*16*17);
    memcpy(bestCofAC16[1], img->cofAC16[1], sizeof(int)*16*17);
  }
#endif 
#ifdef MB32X32
  if(input->UseExtMB > 0)
  {
    memcpy(bestCofAC16x16[0],    img->cofAC16x16[0],    sizeof(int)*257);
    memcpy(bestCofAC16x16[1],    img->cofAC16x16[1],    sizeof(int)*257);
    memcpy(bestCofAC16x8 [0][0], img->cofAC16x8 [0][0], sizeof(int)*129);
    memcpy(bestCofAC16x8 [0][1], img->cofAC16x8 [0][1], sizeof(int)*129);
    memcpy(bestCofAC16x8 [1][0], img->cofAC16x8 [1][0], sizeof(int)*129);
    memcpy(bestCofAC16x8 [1][1], img->cofAC16x8 [1][1], sizeof(int)*129);
    memcpy(bestCofAC8x16 [0][0], img->cofAC8x16 [0][0], sizeof(int)*129);
    memcpy(bestCofAC8x16 [0][1], img->cofAC8x16 [0][1], sizeof(int)*129);
    memcpy(bestCofAC8x16 [1][0], img->cofAC8x16 [1][0], sizeof(int)*129);
    memcpy(bestCofAC8x16 [1][1], img->cofAC8x16 [1][1], sizeof(int)*129);
  }
#endif

}


/*!
*************************************************************************************
* \brief
*    Set stored macroblock parameters
*************************************************************************************
*/
void set_stored_macroblock_parameters ()
{
  int  i, j, k, ****i4p, ***i3p;
  Macroblock  *currMB  = &img->mb_data[img->current_mb_nr];
  int         mode     = best_mode;
  int         bframe   = (img->type==B_SLICE);
  char    **ipredmodes = img->ipredmode;
  
  imgpel        **imgY  = enc_picture->imgY;
  imgpel       ***imgUV = enc_picture->imgUV;
  int        block_x, block_y;  
  short   *cur_mv;
  
#ifdef ADAPTIVE_FD_SD_CODING
  currMB->SD_Coding_on_off=best_SD_Coding_on_off;
  memcpy(currMB->quantizer_indices,best_quantizer_indices,256*sizeof(int));
  memcpy(currMB->SD_or_FD,best_SD_or_FD,4*sizeof(int));
  currMB->SD_or_FD_t8x8=best_SD_or_FD_t8x8;
#endif
#ifdef ADAPTIVE_QUANTIZATION
  currMB->mb_iaqms_idx = best_mb_iaqms_idx;
#endif

  //===== reconstruction values =====
  for (j = 0; j < MB_BLOCK_SIZE; j++)
  {
    memcpy(&imgY[img->pix_y+j][img->pix_x],rec_mbY[j], MB_BLOCK_SIZE * sizeof(imgpel));
    if((img->type==SP_SLICE) &&(si_frame_indicator==0 && sp2_frame_indicator==0 ))
      memcpy(&lrec[img->pix_y+j][img->pix_x],lrec_rec[j], MB_BLOCK_SIZE * sizeof(int)); //restore coeff SP frame
    if(img->MbaffFrameFlag)
      memcpy(rdopt->rec_mbY[j],rec_mbY[j], MB_BLOCK_SIZE * sizeof(imgpel));
  }
#ifdef ADAPTIVE_FD_SD_CODING
  if(img->MbaffFrameFlag)
  {
    rdopt->best_SD_Coding_on_off=best_SD_Coding_on_off;
    memcpy(rdopt->best_quantizer_indices,best_quantizer_indices,256*sizeof(int));
    memcpy(rdopt->best_SD_or_FD,best_SD_or_FD,4*sizeof(int));
    rdopt->best_SD_or_FD_t8x8=best_SD_or_FD_t8x8;
  }
#endif
  
  if (img->AdaptiveRounding)
  {
    update_offset_params(mode,luma_transform_size_8x8_flag);
  }
  
  if (img->yuv_format != YUV400)
  {
    for (j = 0; j<img->mb_cr_size_y; j++)
    {
      memcpy(&imgUV[0][img->pix_c_y+j][img->pix_c_x],rec_mbU[j], img->mb_cr_size_x * sizeof(imgpel));
      memcpy(&imgUV[1][img->pix_c_y+j][img->pix_c_x],rec_mbV[j], img->mb_cr_size_x * sizeof(imgpel));
      if((img->type==SP_SLICE) &&(!si_frame_indicator && !sp2_frame_indicator))
      {
        memcpy(&lrec_uv[0][img->pix_c_y+j][img->pix_c_x],lrec_rec_U[j], img->mb_cr_size_x * sizeof(int));
        memcpy(&lrec_uv[1][img->pix_c_y+j][img->pix_c_x],lrec_rec_V[j], img->mb_cr_size_x * sizeof(int));
      }
      if(img->MbaffFrameFlag)
      {
        memcpy(rdopt->rec_mbU[j],rec_mbU[j], img->mb_cr_size_x * sizeof(imgpel));
        memcpy(rdopt->rec_mbV[j],rec_mbV[j], img->mb_cr_size_x * sizeof(imgpel));
      }
    }
  }
  
  //===== coefficients and cbp =====
  i4p=cofAC; cofAC=img->cofAC; img->cofAC=i4p;
  i3p=cofDC; cofDC=img->cofDC; img->cofDC=i3p;
  currMB->cbp      = cbp;
  currMB->cbp_blk = cbp_blk;
  //==== macroblock type ====
  currMB->mb_type = mode;
  
  if(img->MbaffFrameFlag)
  {
    rdopt->mode = mode;
    rdopt->i16offset = img->i16offset;
    rdopt->cbp = cbp;
    rdopt->cbp_blk = cbp_blk;
    rdopt->mb_type  = mode;
    
    rdopt->prev_qp=currMB->prev_qp;
    rdopt->prev_delta_qp=currMB->prev_delta_qp;
    rdopt->delta_qp = currMB->delta_qp;
    rdopt->qp=currMB->qp;
    rdopt->prev_cbp=currMB->prev_cbp;

    for(i = 0;i<4+img->num_blk8x8_uv;i++)
    {
      for(j = 0;j<4;j++)
        for(k = 0;k<2;k++)
          memcpy(rdopt->cofAC[i][j][k], img->cofAC[i][j][k], 65 * sizeof(int));
    }     
    for(i = 0;i<3;i++)
      for(k = 0;k<2;k++)
        memcpy(rdopt->cofDC[i][k], img->cofDC[i][k], 18 * sizeof(int));
  }
  
  
  memcpy(currMB->b8mode,b8mode, BLOCK_MULTIPLE * sizeof(int));
  memcpy(currMB->b8pdir,b8pdir, BLOCK_MULTIPLE * sizeof(int));
  if(img->MbaffFrameFlag)
  {
    memcpy(rdopt->b8mode,b8mode, BLOCK_MULTIPLE * sizeof(int));
    memcpy(rdopt->b8pdir,b8pdir, BLOCK_MULTIPLE * sizeof(int));    
  }
  
  currMB->bi_pred_me = currMB->mb_type == 1 ? bi_pred_me : 0;  
  
  
  //if P8x8 mode and transform size 4x4 choosen, restore motion vector data for this transform size 
  if (mode == P8x8 && !luma_transform_size_8x8_flag && input->Transform8x8Mode)
    RestoreMV8x8(1);
  
  //==== transform size flag ====
  if (((currMB->cbp & 15) == 0) && !(IS_OLDINTRA(currMB) || currMB->mb_type == I8MB))
    currMB->luma_transform_size_8x8_flag = 0;
  else
    currMB->luma_transform_size_8x8_flag = luma_transform_size_8x8_flag;
  
  rdopt->luma_transform_size_8x8_flag  = currMB->luma_transform_size_8x8_flag;
  
  if (input->rdopt==3 && img->type!=B_SLICE)
  {
    //! save the MB Mode of every macroblock
    decs->dec_mb_mode[img->mb_x][img->mb_y] = mode;
  }
  
  //==== reference frames =====
  for (j = 0; j < 4; j++)
  {
    block_y = img->block_y + j;
    for (i = 0; i < 4; i++)
    {
      block_x = img->block_x + i;
      k = 2*(j >> 1)+(i >> 1);
      
      // backward prediction or intra
      if ((currMB->b8pdir[k] == 1) || IS_INTRA(currMB))
      {
        enc_picture->ref_idx    [LIST_0][block_y][block_x]    = -1;
        enc_picture->ref_pic_id [LIST_0][block_y][block_x]    = -1;          
        enc_picture->mv         [LIST_0][block_y][block_x][0] = 0;
        enc_picture->mv         [LIST_0][block_y][block_x][1] = 0;
        if(img->MbaffFrameFlag)
          rdopt->refar[LIST_0][j][i] = -1;
      }
      else
      {
        if (currMB->bi_pred_me && (currMB->b8pdir[k] == 2) && currMB->mb_type==1)
        {
          cur_mv = currMB->bi_pred_me == 1 
            ? img->bipred_mv1[j][i][LIST_0][0][currMB->b8mode[k]] 
            : img->bipred_mv2[j][i][LIST_0][0][currMB->b8mode[k]];
          
          enc_picture->ref_idx    [LIST_0][block_y][block_x] = 0;                         
          enc_picture->ref_pic_id [LIST_0][block_y][block_x] = enc_picture->ref_pic_num[LIST_0 + currMB->list_offset][0];
          enc_picture->mv         [LIST_0][block_y][block_x][0] = cur_mv[0];
          enc_picture->mv         [LIST_0][block_y][block_x][1] = cur_mv[1];
          if(img->MbaffFrameFlag)
            rdopt->refar[LIST_0][j][i] = 0;        
        }
        else
        {
          cur_mv = img->all_mv[j][i][LIST_0][(short)frefframe[j][i]][currMB->b8mode[k]];
          
          enc_picture->ref_idx    [LIST_0][block_y][block_x]    = frefframe[j][i];
          enc_picture->ref_pic_id [LIST_0][block_y][block_x]    = enc_picture->ref_pic_num[LIST_0 + currMB->list_offset][(short)frefframe[j][i]];
          enc_picture->mv         [LIST_0][block_y][block_x][0] = cur_mv[0];
          enc_picture->mv         [LIST_0][block_y][block_x][1] = cur_mv[1];
          if(img->MbaffFrameFlag)
            rdopt->refar[LIST_0][j][i] = frefframe[j][i];
        }
      }
      
      // forward prediction or intra
      if ((currMB->b8pdir[k] == 0) || IS_INTRA(currMB))
      {
        enc_picture->ref_idx    [LIST_1][block_y][block_x]    = -1;
        enc_picture->ref_pic_id [LIST_1][block_y][block_x]    = -1;
        enc_picture->mv         [LIST_1][block_y][block_x][0] = 0;
        enc_picture->mv         [LIST_1][block_y][block_x][1] = 0;
        if(img->MbaffFrameFlag)
          rdopt->refar[LIST_1][j][i] = -1;
      }
    }
  }
  
  if (bframe)
  {
    for (j=0; j<4; j++)
    {
      block_y = img->block_y + j;
      for (i=0; i<4; i++)
      {          
        block_x = img->block_x + i;
        k = 2*(j >> 1)+(i >> 1);
        
        // forward
        if (IS_INTRA(currMB)||(currMB->b8pdir[k] == 0))
        {
          enc_picture->ref_idx    [LIST_1][block_y][block_x]    = -1;
          enc_picture->ref_pic_id [LIST_1][block_y][block_x]    = -1;
          enc_picture->mv         [LIST_1][block_y][block_x][0] = 0;
          enc_picture->mv         [LIST_1][block_y][block_x][1] = 0;
          if(img->MbaffFrameFlag)
            rdopt->refar[LIST_1][j][i] = -1;
        }
        else
        {
          if (currMB->bi_pred_me && (currMB->b8pdir[k] == 2) && currMB->mb_type==1)
          {
            cur_mv = currMB->bi_pred_me == 1 
              ? img->bipred_mv1[j][i][LIST_1][0][currMB->b8mode[k]] 
              : img->bipred_mv2[j][i][LIST_1][0][currMB->b8mode[k]];
            
            enc_picture->ref_idx    [LIST_1][block_y][block_x] = 0; 
            enc_picture->ref_pic_id [LIST_1][block_y][block_x] = enc_picture->ref_pic_num[LIST_1 + currMB->list_offset][0];
            enc_picture->mv         [LIST_1][block_y][block_x][0] = cur_mv[0];
            enc_picture->mv         [LIST_1][block_y][block_x][1] = cur_mv[1];        
            if(img->MbaffFrameFlag)
              rdopt->refar[LIST_1][j][i] = 0;        
          }
          else
          {
            cur_mv = img->all_mv[j][i][LIST_1][(short)brefframe[j][i]][currMB->b8mode[k]];
            
            enc_picture->ref_idx    [LIST_1][block_y][block_x] = brefframe[j][i];
            enc_picture->ref_pic_id [LIST_1][block_y][block_x] = enc_picture->ref_pic_num[LIST_1 + currMB->list_offset][(short)brefframe[j][i]];
            enc_picture->mv         [LIST_1][block_y][block_x][0] = cur_mv[0];
            enc_picture->mv         [LIST_1][block_y][block_x][1] = cur_mv[1];
            if(img->MbaffFrameFlag)
              rdopt->refar[LIST_1][j][i] = brefframe[j][i];
          }
        }
      }
    }
  }
  
  //==== intra prediction modes ====
  currMB->c_ipred_mode = best_c_imode;
  img->i16offset = best_i16offset;
  
  if(currMB->mb_type == I8MB)
  {
    memcpy(currMB->intra_pred_modes8x8,b8_intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));            
    memcpy(currMB->intra_pred_modes,b8_intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));            
    for(j = 0; j < BLOCK_MULTIPLE; j++)
    {
      memcpy(&img->ipredmode[img->block_y+j][img->block_x],b8_ipredmode8x8[j], BLOCK_MULTIPLE * sizeof(char));
      memcpy(&img->ipredmode8x8[img->block_y+j][img->block_x], b8_ipredmode8x8[j], BLOCK_MULTIPLE * sizeof(char));
    }
  }
  else if (mode!=I4MB && mode!=I8MB)
  {
    memset(currMB->intra_pred_modes,DC_PRED, MB_BLOCK_PARTITIONS * sizeof(char));
    for(j = img->block_y; j < img->block_y + BLOCK_MULTIPLE; j++)
      memset(&img->ipredmode[j][img->block_x], DC_PRED, BLOCK_MULTIPLE * sizeof(char));
  }
  // Residue Color Transform
  else if (mode == I4MB)
  {
    memcpy(currMB->intra_pred_modes,b4_intra_pred_modes, MB_BLOCK_PARTITIONS * sizeof(char));
    for(j = 0; j < BLOCK_MULTIPLE; j++)
      memcpy(&img->ipredmode[img->block_y + j][img->block_x],&b4_ipredmode[BLOCK_MULTIPLE * j], BLOCK_MULTIPLE * sizeof(char));
  }
  
  if(img->MbaffFrameFlag)
  {
    rdopt->c_ipred_mode = currMB->c_ipred_mode;
    rdopt->i16offset = img->i16offset;  
    memcpy(rdopt->intra_pred_modes,currMB->intra_pred_modes, MB_BLOCK_PARTITIONS * sizeof(char));
    memcpy(rdopt->intra_pred_modes8x8,currMB->intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));
    for(j = img->block_y; j < img->block_y +BLOCK_MULTIPLE; j++)
      memcpy(&rdopt->ipredmode[j][img->block_x],&ipredmodes[j][img->block_x], BLOCK_MULTIPLE * sizeof(char));
  }
  
  //==== motion vectors =====
  SetMotionVectorsMB (currMB, bframe);
  
#ifdef USE_INTRA_MDDT
  if(input->UseIntraMDDT)
  {
    memcpy(img->cofAC16[0], bestCofAC16[0], sizeof(int)*16*17);
    memcpy(img->cofAC16[1], bestCofAC16[1], sizeof(int)*16*17);
  }
#endif 
#ifdef MB32X32
  if(input->UseExtMB > 0)
  {
    memcpy(img->cofAC16x16[0],    bestCofAC16x16[0],    sizeof(int)*257);
    memcpy(img->cofAC16x16[1],    bestCofAC16x16[1],    sizeof(int)*257);
    memcpy(img->cofAC16x8 [0][0], bestCofAC16x8 [0][0], sizeof(int)*129);
    memcpy(img->cofAC16x8 [0][1], bestCofAC16x8 [0][1], sizeof(int)*129);
    memcpy(img->cofAC16x8 [1][0], bestCofAC16x8 [1][0], sizeof(int)*129);
    memcpy(img->cofAC16x8 [1][1], bestCofAC16x8 [1][1], sizeof(int)*129);
    memcpy(img->cofAC8x16 [0][0], bestCofAC8x16 [0][0], sizeof(int)*129);
    memcpy(img->cofAC8x16 [0][1], bestCofAC8x16 [0][1], sizeof(int)*129);
    memcpy(img->cofAC8x16 [1][0], bestCofAC8x16 [1][0], sizeof(int)*129);
    memcpy(img->cofAC8x16 [1][1], bestCofAC8x16 [1][1], sizeof(int)*129);
  }
#endif
}



/*!
*************************************************************************************
* \brief
*    Set reference frames and motion vectors
*************************************************************************************
*/
void SetRefAndMotionVectors (int block, int mode, int pdir, int fwref, int bwref)
{
  int     i, j=0;
  int     bslice  = (img->type==B_SLICE);
  int     pmode   = (mode==1||mode==2||mode==3?mode:4);
  int     j0      = ((block >> 1)<<1);
  int     i0      = ((block & 0x01)<<1);
  int     j1      = j0 + (input->part_size[pmode][1]);
  int     i1      = i0 + (input->part_size[pmode][0]);
  int     block_x, block_y;
  short   *cur_mv;
  Macroblock  *currMB  = &img->mb_data[img->current_mb_nr];
  
  if (pdir<0)
  {
    for (j = img->block_y + j0; j < img->block_y + j1; j++)
    {
      for (i=img->block_x + i0; i<img->block_x +i1; i++)
      {
        enc_picture->ref_pic_id[LIST_0][j][i] = -1;
        enc_picture->ref_pic_id[LIST_1][j][i] = -1;
      }
      memset(&enc_picture->ref_idx[LIST_0][j][img->block_x + i0], -1, (input->part_size[pmode][0]) * sizeof(char));
      memset(&enc_picture->ref_idx[LIST_1][j][img->block_x + i0], -1, (input->part_size[pmode][0]) * sizeof(char));
      memset(enc_picture->mv[LIST_0][j][img->block_x + i0], 0, 2*(input->part_size[pmode][0]) * sizeof(short));
      memset(enc_picture->mv[LIST_1][j][img->block_x + i0], 0, 2*(input->part_size[pmode][0]) * sizeof(short));
    }
    return;
  }
  
  if (!bslice)
  {
    for (j=j0; j<j1; j++)
    {
      block_y = img->block_y + j;
      memset(&enc_picture->ref_idx   [LIST_0][block_y][img->block_x + i0], fwref, (input->part_size[pmode][0]) * sizeof(char));
      for (i=i0; i<i1; i++)
      {
        block_x = img->block_x + i;        
        cur_mv = img->all_mv[j][i][LIST_0][fwref][mode];                
        enc_picture->mv        [LIST_0][block_y][block_x][0] = cur_mv[0];
        enc_picture->mv        [LIST_0][block_y][block_x][1] = cur_mv[1];
        enc_picture->ref_pic_id[LIST_0][block_y][block_x] = enc_picture->ref_pic_num[LIST_0+currMB->list_offset][fwref];
      }
    }
    return;
  }
  else
  {
    for (j=j0; j<j1; j++)
    {
      block_y = img->block_y + j;
      for (i=i0; i<i1; i++)
      {
        block_x = img->block_x + i;
        if (mode==0)
        {
          pdir  = direct_pdir[block_y][block_x];
          fwref = direct_ref_idx[LIST_0][block_y][block_x];
          bwref = direct_ref_idx[LIST_1][block_y][block_x];
        }
        
        if ((pdir==0 || pdir==2))
        {
          if (currMB->bi_pred_me && (pdir == 2) && mode == 1)
          {
            cur_mv = currMB->bi_pred_me == 1 
              ? img->bipred_mv1[j][i][LIST_0][0][mode]
              : img->bipred_mv2[j][i][LIST_0][0][mode];
            
            enc_picture->mv        [LIST_0][block_y][block_x][0] = cur_mv[0];
            enc_picture->mv        [LIST_0][block_y][block_x][1] = cur_mv[1];
            enc_picture->ref_idx   [LIST_0][block_y][block_x]    = 0;            
            enc_picture->ref_pic_id[LIST_0][block_y][block_x]    = enc_picture->ref_pic_num[LIST_0+currMB->list_offset][0];
          }
          else
          {
            cur_mv = img->all_mv[j][i][LIST_0][fwref][mode];
            
            enc_picture->mv        [LIST_0][block_y][block_x][0] = cur_mv[0];
            enc_picture->mv        [LIST_0][block_y][block_x][1] = cur_mv[1];
            enc_picture->ref_idx   [LIST_0][block_y][block_x] = fwref;
            enc_picture->ref_pic_id[LIST_0][block_y][block_x] = 
              enc_picture->ref_pic_num[LIST_0+currMB->list_offset][(short)enc_picture->ref_idx[LIST_0][block_y][block_x]];
          }
        }
        else
        {
          enc_picture->mv        [LIST_0][block_y][block_x][0] = 0;
          enc_picture->mv        [LIST_0][block_y][block_x][1] = 0;
          enc_picture->ref_idx   [LIST_0][block_y][block_x]    = -1;
          enc_picture->ref_pic_id[LIST_0][block_y][block_x]    = -1;
        }
        
        if ((pdir==1 || pdir==2))
        {
          if (currMB->bi_pred_me && (pdir == 2) && mode == 1)
          {
            cur_mv = currMB->bi_pred_me == 1 
              ? img->bipred_mv1[j][i][LIST_1][0][mode]
              : img->bipred_mv2[j][i][LIST_1][0][mode];
            
            enc_picture->mv        [LIST_1][block_y][block_x][0] = cur_mv[0];
            enc_picture->mv        [LIST_1][block_y][block_x][1] = cur_mv[1];
            enc_picture->ref_idx   [LIST_1][block_y][block_x]    = 0;            
            enc_picture->ref_pic_id[LIST_1][block_y][block_x]    = enc_picture->ref_pic_num[LIST_1+currMB->list_offset][0];
          }
          else
          {
            cur_mv = img->all_mv[j][i][LIST_1][bwref][mode];
            
            enc_picture->mv        [LIST_1][block_y][block_x][0] = cur_mv[0];
            enc_picture->mv        [LIST_1][block_y][block_x][1] = cur_mv[1];
            enc_picture->ref_idx   [LIST_1][block_y][block_x] = bwref;
            enc_picture->ref_pic_id[LIST_1][block_y][block_x] = 
              enc_picture->ref_pic_num[LIST_1+currMB->list_offset][(short)enc_picture->ref_idx[LIST_1][block_y][block_x]];
          }
        }
        else
        {
          enc_picture->mv        [LIST_1][block_y][block_x][0] = 0;
          enc_picture->mv        [LIST_1][block_y][block_x][1] = 0;
          enc_picture->ref_idx   [LIST_1][block_y][block_x]    = -1;
          enc_picture->ref_pic_id[LIST_1][block_y][block_x]    = -1;
        }
      }
    }
  }
}

/*!
*************************************************************************************
* \brief
*    skip macroblock field inference
* \return
*    inferred field flag
*************************************************************************************
*/
int field_flag_inference()
{
  int mb_field;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  
  if (currMB->mbAvailA)
  {
    mb_field = img->mb_data[currMB->mbAddrA].mb_field;
  }
  else
  {
    // check top macroblock pair
    if (currMB->mbAvailB)
      mb_field = img->mb_data[currMB->mbAddrB].mb_field;
    else
      mb_field = 0;
  }
  
  return mb_field;
}

/*!
*************************************************************************************
* \brief
*    Store motion vectors for 8x8 partition
*************************************************************************************
*/

void StoreMVBlock8x8(int dir, int block8x8, int mode, int ref, int bw_ref, int pdir8, int bframe)
{
  int i, j, i0, j0, ii, jj;
  short ******all_mv  = img->all_mv;
  short ******pred_mv = img->pred_mv;
  
  
  i0 = (block8x8 & 0x01) << 1;
  j0 = (block8x8 >> 1) << 1;
  ii = i0+2;
  jj = j0+2;
  
  if (!bframe)
  {
    if (pdir8>=0) //(mode8!=IBLOCK)&&(mode8!=I16MB))  // && ref != -1)
    {
      for (j=j0; j<jj; j++)
        for (i=i0; i<ii; i++)
        {
          all_mv8x8 [dir][LIST_0][j][i][0] = all_mv [j][i][LIST_0][ref][4][0];
          all_mv8x8 [dir][LIST_0][j][i][1] = all_mv [j][i][LIST_0][ref][4][1];
          pred_mv8x8[dir][LIST_0][j][i][0] = pred_mv[j][i][LIST_0][ref][4][0];
          pred_mv8x8[dir][LIST_0][j][i][1] = pred_mv[j][i][LIST_0][ref][4][1];
        }
    }
  }
  else
  {
    if (pdir8==0) // forward
    {
      for (j=j0; j<jj; j++)
        for (i=i0; i<ii; i++)
        {
          all_mv8x8 [dir][LIST_0][j][i][0] = all_mv [j][i][LIST_0][ref][mode][0];
          all_mv8x8 [dir][LIST_0][j][i][1] = all_mv [j][i][LIST_0][ref][mode][1];
          pred_mv8x8[dir][LIST_0][j][i][0] = pred_mv[j][i][LIST_0][ref][mode][0];
          pred_mv8x8[dir][LIST_0][j][i][1] = pred_mv[j][i][LIST_0][ref][mode][1];
        }
    }
    else if (pdir8==1) // backward
    {
      for (j=j0; j<jj; j++)
        for (i=i0; i<ii; i++)
        {
          all_mv8x8 [dir][LIST_1][j][i][0] = all_mv [j][i][LIST_1][bw_ref][mode][0];
          all_mv8x8 [dir][LIST_1][j][i][1] = all_mv [j][i][LIST_1][bw_ref][mode][1];
          pred_mv8x8[dir][LIST_1][j][i][0] = pred_mv[j][i][LIST_1][bw_ref][mode][0];
          pred_mv8x8[dir][LIST_1][j][i][1] = pred_mv[j][i][LIST_1][bw_ref][mode][1];
        }
    }
    else if (pdir8==2) // bidir
    {
      for (j=j0; j<jj; j++)
        for (i=i0; i<ii; i++)
        {
          all_mv8x8 [dir][LIST_0][j][i][0] = all_mv [j][i][LIST_0][ref][mode][0];
          all_mv8x8 [dir][LIST_0][j][i][1] = all_mv [j][i][LIST_0][ref][mode][1];
          pred_mv8x8[dir][LIST_0][j][i][0] = pred_mv[j][i][LIST_0][ref][mode][0];
          pred_mv8x8[dir][LIST_0][j][i][1] = pred_mv[j][i][LIST_0][ref][mode][1];
          
          all_mv8x8 [dir][LIST_1][j][i][0] = all_mv [j][i][LIST_1][bw_ref][mode][0];
          all_mv8x8 [dir][LIST_1][j][i][1] = all_mv [j][i][LIST_1][bw_ref][mode][1];
          pred_mv8x8[dir][LIST_1][j][i][0] = pred_mv[j][i][LIST_1][bw_ref][mode][0];
          pred_mv8x8[dir][LIST_1][j][i][1] = pred_mv[j][i][LIST_1][bw_ref][mode][1];
        }
    }
    else
    {
      error("invalid direction mode", 255);
    }
  }
}



/*!
*************************************************************************************
* \brief
*    Store motion vectors of 8x8 partitions of one macroblock
*************************************************************************************
*/
void StoreMV8x8(int dir)
{
  int block8x8;
  
  int bframe = (img->type == B_SLICE);
  
  for (block8x8=0; block8x8<4; block8x8++)
    StoreMVBlock8x8(dir, block8x8, tr8x8.part8x8mode[block8x8], tr8x8.part8x8fwref[block8x8], 
    tr8x8.part8x8bwref[block8x8], tr8x8.part8x8pdir[block8x8], bframe);
}

/*!
*************************************************************************************
* \brief
*    Restore motion vectors for 8x8 partition
*************************************************************************************
*/
void RestoreMVBlock8x8(int dir, int block8x8, RD_8x8DATA tr, int bframe)
{
  int i, j, i0, j0, ii, jj;
  short ******all_mv  = img->all_mv;
  short ******pred_mv = img->pred_mv;
  short pdir8  = tr.part8x8pdir [block8x8];
  short mode   = tr.part8x8mode [block8x8];
  short ref    = tr.part8x8fwref[block8x8];
  short bw_ref = tr.part8x8bwref[block8x8];
  
  i0 = (block8x8 & 0x01) << 1;
  j0 = (block8x8 >> 1) << 1;
  ii = i0+2;
  jj = j0+2;
  
  if (!bframe)
  {
    if (pdir8>=0) //(mode8!=IBLOCK)&&(mode8!=I16MB))  // && ref != -1)
    {
      for (j=j0; j<jj; j++)
        for (i=i0; i<ii; i++)
        {
          all_mv [j][i][LIST_0][ref][4][0] = all_mv8x8 [dir][LIST_0][j][i][0] ;
          all_mv [j][i][LIST_0][ref][4][1] = all_mv8x8 [dir][LIST_0][j][i][1] ;
          pred_mv[j][i][LIST_0][ref][4][0] = pred_mv8x8[dir][LIST_0][j][i][0];
          pred_mv[j][i][LIST_0][ref][4][1] = pred_mv8x8[dir][LIST_0][j][i][1];
        }
    }
  }
  else
  {
    if (pdir8==0) // forward
    {
      for (j=j0; j<jj; j++)
        for (i=i0; i<ii; i++)
        {
          all_mv [j][i][LIST_0][ref][mode][0] = all_mv8x8 [dir][LIST_0][j][i][0] ;
          all_mv [j][i][LIST_0][ref][mode][1] = all_mv8x8 [dir][LIST_0][j][i][1] ;
          pred_mv[j][i][LIST_0][ref][mode][0] = pred_mv8x8[dir][LIST_0][j][i][0];
          pred_mv[j][i][LIST_0][ref][mode][1] = pred_mv8x8[dir][LIST_0][j][i][1];
        }
    }
    else if (pdir8==1) // backward
    {
      for (j=j0; j<jj; j++)
        for (i=i0; i<ii; i++)
        {
          all_mv [j][i][LIST_1][bw_ref][mode][0] = all_mv8x8 [dir][LIST_1][j][i][0] ;
          all_mv [j][i][LIST_1][bw_ref][mode][1] = all_mv8x8 [dir][LIST_1][j][i][1] ;
          pred_mv[j][i][LIST_1][bw_ref][mode][0] = pred_mv8x8[dir][LIST_1][j][i][0];
          pred_mv[j][i][LIST_1][bw_ref][mode][1] = pred_mv8x8[dir][LIST_1][j][i][1];
        }
    }
    else if (pdir8==2) // bidir
    {
      for (j=j0; j<jj; j++)
        for (i=i0; i<ii; i++)
        {
          all_mv [j][i][LIST_0][ref][mode][0] = all_mv8x8 [dir][LIST_0][j][i][0] ;
          all_mv [j][i][LIST_0][ref][mode][1] = all_mv8x8 [dir][LIST_0][j][i][1] ;
          pred_mv[j][i][LIST_0][ref][mode][0] = pred_mv8x8[dir][LIST_0][j][i][0];
          pred_mv[j][i][LIST_0][ref][mode][1] = pred_mv8x8[dir][LIST_0][j][i][1];
          
          all_mv [j][i][LIST_1][bw_ref][mode][0] = all_mv8x8 [dir][LIST_1][j][i][0] ;
          all_mv [j][i][LIST_1][bw_ref][mode][1] = all_mv8x8 [dir][LIST_1][j][i][1] ;
          pred_mv[j][i][LIST_1][bw_ref][mode][0] = pred_mv8x8[dir][LIST_1][j][i][0];
          pred_mv[j][i][LIST_1][bw_ref][mode][1] = pred_mv8x8[dir][LIST_1][j][i][1];
        }
    }
    else
    {
      error("invalid direction mode", 255);
    }
  }
}

/*!
*************************************************************************************
* \brief
*    Restore motion vectors of 8x8 partitions of one macroblock
*************************************************************************************
*/
void RestoreMV8x8(int dir)
{
  int block8x8;
  
  int bframe = (img->type == B_SLICE);
  
  for (block8x8=0; block8x8<4; block8x8++)
    RestoreMVBlock8x8(dir, block8x8, tr8x8, bframe);   
}


/*!
*************************************************************************************
* \brief
*    Store predictors for 8x8 partition
*************************************************************************************
*/

void StoreNewMotionVectorsBlock8x8(int dir, int block8x8, int mode, int fw_ref, int bw_ref, int pdir8, int bframe)
{
  int i, j, i0, j0, ii, jj;
  short ******all_mv  = img->all_mv;
  short ******pred_mv = img->pred_mv;
  
  
  i0 = (block8x8 & 0x01) << 1;
  j0 = (block8x8 >> 1) << 1;
  ii = i0+2;
  jj = j0+2;
  
  if (pdir8<0)
  {
    for (j=j0; j<jj; j++)
    {
      memset(&all_mv8x8[dir][LIST_0][j][i0], 0, 2 * 2 * sizeof(short));
      memset(&all_mv8x8[dir][LIST_1][j][i0], 0, 2 * 2 * sizeof(short));
    }
    return;
  }
  
  if (!bframe)
  {
    for (j=j0; j<jj; j++)
    {
      for (i=i0; i<ii; i++)
      {
        all_mv8x8 [dir][LIST_0][j][i][0] = all_mv [j][i][LIST_0][fw_ref][4][0];
        all_mv8x8 [dir][LIST_0][j][i][1] = all_mv [j][i][LIST_0][fw_ref][4][1];
        pred_mv8x8[dir][LIST_0][j][i][0] = pred_mv[j][i][LIST_0][fw_ref][4][0];
        pred_mv8x8[dir][LIST_0][j][i][1] = pred_mv[j][i][LIST_0][fw_ref][4][1];
      }
    }
    return;
  }
  else
  {
    if ((pdir8==0 || pdir8==2))
    {
      for (j=j0; j<jj; j++)
        for (i=i0; i<ii; i++)
        {
          all_mv8x8 [dir][LIST_0][j][i][0] = all_mv [j][i][LIST_0][fw_ref][mode][0];
          all_mv8x8 [dir][LIST_0][j][i][1] = all_mv [j][i][LIST_0][fw_ref][mode][1];
          pred_mv8x8[dir][LIST_0][j][i][0] = pred_mv[j][i][LIST_0][fw_ref][mode][0];
          pred_mv8x8[dir][LIST_0][j][i][1] = pred_mv[j][i][LIST_0][fw_ref][mode][1];
        }
    }
    else
    {
      for (j=j0; j<jj; j++)
        memset(&all_mv8x8[dir][LIST_0][j][i0], 0, 2 * 2 * sizeof(short));
    }
    
    if ((pdir8==1 || pdir8==2))
    {
      for (j=j0; j<jj; j++)
        for (i=i0; i<ii; i++)
        {
          all_mv8x8 [dir][LIST_1][j][i][0] = all_mv [j][i][LIST_1][bw_ref][mode][0];
          all_mv8x8 [dir][LIST_1][j][i][1] = all_mv [j][i][LIST_1][bw_ref][mode][1];
          pred_mv8x8[dir][LIST_1][j][i][0] = pred_mv[j][i][LIST_1][bw_ref][mode][0];
          pred_mv8x8[dir][LIST_1][j][i][1] = pred_mv[j][i][LIST_1][bw_ref][mode][1];
        }
    }
    else
    {
      for (j=j0; j<jj; j++)
        memset(&all_mv8x8[dir][LIST_1][j][i0], 0, 2 * 2 * sizeof(short));
    }
  }
}

/*!
************************************************************************
* \brief
*    Makes the decision if 8x8 tranform will be used (for RD-off)
************************************************************************
*/
int GetBestTransformP8x8()
{
  int    block_y, block_x, pic_pix_y, pic_pix_x, i, j, k;
  int    mb_y, mb_x, block8x8;
  int    cost8x8=0, cost4x4=0;
  int    diff4x4[64], *diff_ptr;
  int    diff8x8[64];
  
  if(input->Transform8x8Mode==2) //always allow 8x8 transform
    return 1;
  
  for (block8x8=0; block8x8<4; block8x8++)
  {
    mb_y = (block8x8 >>   1) << 3;
    mb_x = (block8x8 & 0x01) << 3;
    //===== loop over 4x4 blocks =====
    k=0;
    for (block_y=mb_y; block_y<mb_y+8; block_y+=4)
    {
      pic_pix_y = img->opix_y + block_y;
      
      //get cost for transform size 4x4
      for (block_x=mb_x; block_x<mb_x+8; block_x+=4)
      {
        pic_pix_x = img->opix_x + block_x;
        
        //===== get displaced frame difference ======
        diff_ptr=&diff4x4[k];
        for (j=0; j<4; j++)
        {
          for (i=0; i<4; i++, k++)
          {
            //4x4 transform size
            diff4x4[k] = imgY_org[pic_pix_y+j][pic_pix_x+i] - tr4x4.mpr8x8[j+block_y][i+block_x];
            //8x8 transform size
            diff8x8[k] = imgY_org[pic_pix_y+j][pic_pix_x+i] - tr8x8.mpr8x8[j+block_y][i+block_x];
          }
        } 
        
        cost4x4 += SATD (diff_ptr, input->hadamard);
      }
    }    
    cost8x8 += SATD8X8 (diff8x8, input->hadamard);
  }  
  return (cost8x8 < cost4x4);
}

/*!
************************************************************************
* \brief
*    Sets MBAFF RD parameters
************************************************************************
*/
void set_mbaff_parameters()
{
  int  i, j, k;
  Macroblock  *currMB  = &img->mb_data[img->current_mb_nr];
  int         mode     = best_mode;
  int         bframe   = (img->type==B_SLICE);
  char    **ipredmodes = img->ipredmode;
  
  
#ifdef ADAPTIVE_FD_SD_CODING
  rdopt->best_SD_Coding_on_off = currMB->SD_Coding_on_off;
  memcpy(rdopt->best_quantizer_indices,currMB->quantizer_indices,256*sizeof(int));
  memcpy(rdopt->best_SD_or_FD,currMB->SD_or_FD,4*sizeof(int));
  rdopt->best_SD_or_FD_t8x8=currMB->SD_or_FD_t8x8;
#endif
  
  
  //===== reconstruction values =====
  for (j=0; j < MB_BLOCK_SIZE; j++)
    memcpy(rdopt->rec_mbY[j],&enc_picture->imgY[img->pix_y + j][img->pix_x], MB_BLOCK_SIZE * sizeof(imgpel));
  
  if (img->yuv_format != YUV400)
  {
    for (j=0; j<img->mb_cr_size_y; j++)
    {
      memcpy(rdopt->rec_mbU[j],&enc_picture->imgUV[0][img->pix_c_y + j][img->pix_c_x], img->mb_cr_size_x * sizeof(imgpel));
      memcpy(rdopt->rec_mbV[j],&enc_picture->imgUV[1][img->pix_c_y + j][img->pix_c_x], img->mb_cr_size_x * sizeof(imgpel));
    }
  }
  
  //===== coefficients and cbp =====
  rdopt->mode      = mode;
  rdopt->i16offset = img->i16offset;
  rdopt->cbp       = currMB->cbp;
  rdopt->cbp_blk   = currMB->cbp_blk;
  rdopt->mb_type   = currMB->mb_type;
  
  rdopt->luma_transform_size_8x8_flag = currMB->luma_transform_size_8x8_flag;
  if(rdopt->mb_type == 0 && mode != 0)
  {
    mode=0;
    rdopt->mode=0;
  }
  
  for(i=0;i<4+img->num_blk8x8_uv;i++)
  {
    for(j=0;j<4;j++)
      for(k=0;k<2;k++)
        memcpy(rdopt->cofAC[i][j][k], img->cofAC[i][j][k], 65 * sizeof(int));
  }
  
  for(i=0;i<3;i++)
  {
    for(k=0;k<2;k++)
      memcpy(rdopt->cofDC[i][k], img->cofDC[i][k], 18 * sizeof(int));
  }   
  
  memcpy(rdopt->b8mode,currMB->b8mode, BLOCK_MULTIPLE * sizeof(int));
  memcpy(rdopt->b8pdir,currMB->b8pdir, BLOCK_MULTIPLE * sizeof(int));
  
  //==== reference frames =====
  if (bframe)
  {
    for (j = 0; j < BLOCK_MULTIPLE; j++)
    {
      memcpy(rdopt->refar[LIST_0][j],&enc_picture->ref_idx[LIST_0][img->block_y + j][img->block_x] , BLOCK_MULTIPLE * sizeof(char));
      memcpy(rdopt->refar[LIST_1][j],&enc_picture->ref_idx[LIST_1][img->block_y + j][img->block_x] , BLOCK_MULTIPLE * sizeof(char));
    }
    rdopt->bi_pred_me = currMB->bi_pred_me;
  }
  else
  {
    for (j = 0; j < BLOCK_MULTIPLE; j++)
      memcpy(rdopt->refar[LIST_0][j],&enc_picture->ref_idx[LIST_0][img->block_y + j][img->block_x] , BLOCK_MULTIPLE * sizeof(char));
  }      
  
  memcpy(rdopt->intra_pred_modes,currMB->intra_pred_modes, MB_BLOCK_PARTITIONS * sizeof(char));
  memcpy(rdopt->intra_pred_modes8x8,currMB->intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));
  for (j = img->block_y; j < img->block_y + 4; j++)
  {
    memcpy(&rdopt->ipredmode[j][img->block_x],&ipredmodes[j][img->block_x], BLOCK_MULTIPLE * sizeof(char));
  }
}

/*!
************************************************************************
* \brief
*    store coding state (for rd-optimized mode decision), used for 8x8 transformation
************************************************************************
*/
void store_coding_state_cs_cm()
{
  store_coding_state(cs_cm);
}

/*!
************************************************************************
* \brief
*    restore coding state (for rd-optimized mode decision), used for 8x8 transformation
************************************************************************
*/
void reset_coding_state_cs_cm()
{
  reset_coding_state(cs_cm);
}

/*!
************************************************************************
* \brief
*    update rounding offsets based on JVT-N011
************************************************************************
*/
void update_offset_params(int mode, int luma_transform_size_8x8_flag)
{
  int i,j;
  int temp = 0;
  int offsetRange = 1 << (OffsetBits - 1);
  
  for (j=0; j < MB_BLOCK_SIZE; j++)
    for (i=0; i < MB_BLOCK_SIZE; i++)
    {     
      if ((mode != I4MB)&&(mode != I16MB) && (mode != I8MB) )
      {
        if (img->type == B_SLICE)
        {
          if (!luma_transform_size_8x8_flag )
          {
            temp = ((j & 0x03)<<2)+ (i & 0x03);    
            OffsetList4x4[12][temp] += bestInterFAdjust4x4[j][i];
            OffsetList4x4[12][temp] = Clip3(0,offsetRange,OffsetList4x4[12][temp]);
#ifdef ADAPTIVE_FD_SD_CODING
            if (i==0 && j==0)
              adaptive_f_spatial_domain_4x4[1]+=best_adjust_adaptive_f_spatial_domain_4x4;
#endif
          }
          else
          {
            temp = ((j & 0x07)<<3)+ (i & 0x07);
            OffsetList8x8[4][temp] += bestInterFAdjust8x8[j][i];
            OffsetList8x8[4][temp] = Clip3(0,offsetRange,OffsetList8x8[4][temp]);
#ifdef ADAPTIVE_FD_SD_CODING
            if (i==0 && j==0)
              adaptive_f_spatial_domain_8x8[1]+=best_adjust_adaptive_f_spatial_domain_8x8;
#endif
          }
          
          if (input->AdaptRndChroma && (i >> 3) ==0 && (j >> 3) == 0)
          {
            temp = ((j & 0x03)<<2)+ (i & 0x03);
            OffsetList4x4[13][temp] += bestInterFAdjust4x4Cr[0][j][i];
            OffsetList4x4[13][temp] = Clip3(0,offsetRange,OffsetList4x4[13][temp]);
            OffsetList4x4[14][temp] += bestInterFAdjust4x4Cr[1][j][i];
            OffsetList4x4[14][temp] = Clip3(0,offsetRange,OffsetList4x4[14][temp]);
          }
        }
        else
        {
          if (!luma_transform_size_8x8_flag )
          {
            temp = ((j & 0x03)<<2)+ (i & 0x03);    
            OffsetList4x4[9][temp] += bestInterFAdjust4x4[j][i];
            OffsetList4x4[9][temp] = Clip3(0,offsetRange,OffsetList4x4[9][temp]);
#ifdef ADAPTIVE_FD_SD_CODING
            if (i==0 && j==0)
              adaptive_f_spatial_domain_4x4[0]+=best_adjust_adaptive_f_spatial_domain_4x4;
#endif
          }
          else
          {
            temp = ((j & 0x07)<<3)+ (i & 0x07);
            OffsetList8x8[3][temp] += bestInterFAdjust8x8[j][i];
            OffsetList8x8[3][temp] = Clip3(0,offsetRange,OffsetList8x8[3][temp]);
#ifdef ADAPTIVE_FD_SD_CODING
            if (i==0 && j==0)
              adaptive_f_spatial_domain_8x8[0]+=best_adjust_adaptive_f_spatial_domain_8x8;
#endif
          }
          
          if (input->AdaptRndChroma && (i >> 3) ==0 && (j >> 3) == 0)
          {
            temp = ((j & 0x03)<<2)+ (i & 0x03);
            OffsetList4x4[10][temp] += bestInterFAdjust4x4Cr[0][j][i];
            OffsetList4x4[10][temp] = Clip3(0,offsetRange,OffsetList4x4[10][temp]);
            OffsetList4x4[11][temp] += bestInterFAdjust4x4Cr[1][j][i];
            OffsetList4x4[11][temp] = Clip3(0,offsetRange,OffsetList4x4[11][temp]);
          }
        }
      }
      else if (mode != I8MB)
      {
        if (img->type == I_SLICE)
        {
          
          temp = ((j & 0x03)<<2)+ (i & 0x03);
          OffsetList4x4[0][temp] += bestIntraFAdjust4x4[j][i];
          OffsetList4x4[0][temp] = Clip3(0,offsetRange,OffsetList4x4[0][temp]);
          if (input->AdaptRndChroma && (i >> 3) ==0 && (j >> 3) == 0)
          {
            temp = ((j & 0x03)<<2)+ (i & 0x03);
            OffsetList4x4[1][temp] += bestIntraFAdjust4x4Cr[0][j][i];
            OffsetList4x4[1][temp] = Clip3(0,offsetRange,OffsetList4x4[1][temp]);
            OffsetList4x4[2][temp] += bestIntraFAdjust4x4Cr[1][j][i];
            OffsetList4x4[2][temp] = Clip3(0,offsetRange,OffsetList4x4[2][temp]);
          }
        }
        else if (img->type == B_SLICE)
        {
          temp = ((j & 0x03)<<2)+ (i & 0x03);
          OffsetList4x4[6][temp] += bestIntraFAdjust4x4[j][i];
          OffsetList4x4[6][temp] = Clip3(0,offsetRange,OffsetList4x4[6][temp]);
          if (input->AdaptRndChroma && (i >> 3) ==0 && (j >> 3) == 0)
          {
            temp = ((j & 0x03)<<2)+ (i & 0x03);
            OffsetList4x4[7][temp] += bestIntraFAdjust4x4Cr[0][j][i];
            OffsetList4x4[7][temp] = Clip3(0,offsetRange,OffsetList4x4[7][temp]);
            OffsetList4x4[8][temp] += bestIntraFAdjust4x4Cr[1][j][i];
            OffsetList4x4[8][temp] = Clip3(0,offsetRange,OffsetList4x4[8][temp]);
          }
        }
        else
        {
          temp = ((j & 0x03)<<2)+ (i & 0x03);
          OffsetList4x4[3][temp] += bestIntraFAdjust4x4[j][i];
          OffsetList4x4[3][temp] = Clip3(0,offsetRange,OffsetList4x4[3][temp]);
          if (input->AdaptRndChroma && (i >> 3) ==0 && (j >> 3) == 0)
          {
            temp = ((j & 0x03)<<2)+ (i & 0x03);
            OffsetList4x4[4][temp] += bestIntraFAdjust4x4Cr[0][j][i];
            OffsetList4x4[4][temp] = Clip3(0,offsetRange,OffsetList4x4[4][temp]);
            OffsetList4x4[5][temp] += bestIntraFAdjust4x4Cr[1][j][i];
            OffsetList4x4[5][temp] = Clip3(0,offsetRange,OffsetList4x4[5][temp]);
          }
        }
        
      }
      else
      {
        if (img->type == I_SLICE)
        {
          temp = ((j & 0x07)<<3)+ (i & 0x07);
          OffsetList8x8[0][temp] += bestIntraFAdjust8x8[j][i];
          OffsetList8x8[0][temp] = Clip3(0,offsetRange,OffsetList8x8[0][temp]);
          if (input->AdaptRndChroma && (i >> 3) ==0 && (j >> 3) == 0)
          {
            temp = ((j & 0x03)<<2)+ (i & 0x03);
            OffsetList4x4[1][temp] += bestIntraFAdjust4x4Cr[0][j][i];
            OffsetList4x4[1][temp] = Clip3(0,offsetRange,OffsetList4x4[1][temp]);
            OffsetList4x4[2][temp] += bestIntraFAdjust4x4Cr[1][j][i];
            OffsetList4x4[2][temp] = Clip3(0,offsetRange,OffsetList4x4[2][temp]);
          }
        }
        else if (img->type == B_SLICE)
        {
          temp = ((j & 0x07)<<3)+ (i & 0x07);
          OffsetList8x8[2][temp] += bestIntraFAdjust8x8[j][i];
          OffsetList8x8[2][temp] = Clip3(0,offsetRange,OffsetList8x8[2][temp]);
          if (input->AdaptRndChroma && (i >> 3) ==0 && (j >> 3) == 0)
          {
            temp = ((j & 0x03)<<2)+ (i & 0x03);
            OffsetList4x4[7][temp] += bestIntraFAdjust4x4Cr[0][j][i];
            OffsetList4x4[7][temp] = Clip3(0,offsetRange,OffsetList4x4[7][temp]);
            OffsetList4x4[8][temp] += bestIntraFAdjust4x4Cr[1][j][i];
            OffsetList4x4[8][temp] = Clip3(0,offsetRange,OffsetList4x4[8][temp]);
          }
        }
        else
        {
          temp = ((j & 0x07)<<3)+ (i & 0x07);
          OffsetList8x8[1][temp] += bestIntraFAdjust8x8[j][i];
          OffsetList8x8[1][temp] = Clip3(0,offsetRange,OffsetList8x8[1][temp]);
          if (input->AdaptRndChroma && (i >> 3) ==0 && (j >> 3) == 0)
          {
            temp = ((j & 0x03)<<2)+ (i & 0x03);
            OffsetList4x4[4][temp] += bestIntraFAdjust4x4Cr[0][j][i];
            OffsetList4x4[4][temp] = Clip3(0,offsetRange,OffsetList4x4[4][temp]);
            OffsetList4x4[5][temp] += bestIntraFAdjust4x4Cr[1][j][i];
            OffsetList4x4[5][temp] = Clip3(0,offsetRange,OffsetList4x4[5][temp]);
          }
        }
      }
    }
}

void assign_enc_picture_params(int mode, int best_pdir, int block, int list_offset, int best_fw_ref, int best_bw_ref, int bframe)
{
  int i,j;
  int block_x, block_y;
  short *cur_mv;
  
  if (mode==1)
  {
    if (best_pdir==1)
    {
      for (j=img->block_y+(block&2); j<img->block_y+(block&2) + BLOCK_MULTIPLE; j++)
      {
        block_x = img->block_x+(block&1)*2;
        
        memset(&enc_picture->ref_idx[LIST_0][j][block_x], -1 ,     BLOCK_MULTIPLE * sizeof(char));
        memset(enc_picture->mv      [LIST_0][j][block_x],  0 , 2 * BLOCK_MULTIPLE * sizeof(short));
        for (i=block_x; i<block_x + BLOCK_MULTIPLE; i++)
        {
          enc_picture->ref_pic_id [LIST_0][j][i]    = -1;
        }
      }
    }
    else if (img->bi_pred_me[mode])
    {
      for (j=0; j<BLOCK_MULTIPLE; j++)
      {
        block_y = img->block_y+(block&2)+j;
        block_x = img->block_x+(block&1)*2;
        memset(&enc_picture->ref_idx[LIST_0][block_y][block_x], 0, BLOCK_MULTIPLE * sizeof(char));
        for (i=0; i<BLOCK_MULTIPLE; i++)
        {            
          cur_mv = img->bi_pred_me[mode] == 1 
            ? img->bipred_mv1[i][j][LIST_0][0][mode] 
            : img->bipred_mv2[i][j][LIST_0][0][mode];
          
          enc_picture->ref_pic_id [LIST_0][block_y][block_x + i]    = enc_picture->ref_pic_num[LIST_0 + list_offset][0];  
          enc_picture->mv         [LIST_0][block_y][block_x + i][0] = cur_mv[0];
          enc_picture->mv         [LIST_0][block_y][block_x + i][1] = cur_mv[1];
        }
      }
    }
    else 
    {
      for (j=0; j<BLOCK_MULTIPLE; j++)
      {
        block_y = img->block_y+(block&2)+j;
        block_x = img->block_x+(block&1)*2;
        memset(&enc_picture->ref_idx[LIST_0][block_y][block_x], best_fw_ref , BLOCK_MULTIPLE * sizeof(char));
        for (i=0; i<BLOCK_MULTIPLE; i++)
        {                                
          cur_mv = img->all_mv[j][i][LIST_0][best_fw_ref][mode];
          
          enc_picture->ref_pic_id [LIST_0][block_y][block_x + i]    = enc_picture->ref_pic_num[LIST_0 + list_offset][best_fw_ref];  
          enc_picture->mv         [LIST_0][block_y][block_x + i][0] = cur_mv[0];
          enc_picture->mv         [LIST_0][block_y][block_x + i][1] = cur_mv[1];
        }          
      }
    }
    
    if (bframe)
    {
      if (best_pdir==0)
      {
        for (j=img->block_y+(block&2); j<img->block_y+(block&2) + BLOCK_MULTIPLE; j++)
        {
          block_x = img->block_x+(block&1)*2;
          memset(&enc_picture->ref_idx[LIST_1][j][block_x], -1 , BLOCK_MULTIPLE * sizeof(char));
          memset(enc_picture->mv[LIST_1][j][block_x], 0 , 2 * BLOCK_MULTIPLE * sizeof(short));
          for (i=block_x; i<block_x + BLOCK_MULTIPLE; i++)
          {
            enc_picture->ref_pic_id [LIST_1][j][i] = -1;
          }
        }
      }
      else
      {
        if (img->bi_pred_me[mode])
        {
          for (j=0; j<BLOCK_MULTIPLE; j++)
          {
            block_y = img->block_y+(block&2)+j;
            block_x = img->block_x+(block&1)*2;
            memset(&enc_picture->ref_idx[LIST_1][block_y][block_x], 0, BLOCK_MULTIPLE * sizeof(char));
            for (i=0; i<BLOCK_MULTIPLE; i++)
            {                     
              cur_mv = img->bi_pred_me[mode] == 1 
                ? img->bipred_mv1[i][j][LIST_1][0][mode] 
                : img->bipred_mv2[i][j][LIST_1][0][mode];
              
              enc_picture->ref_pic_id [LIST_1][block_y][block_x + i] = 
                enc_picture->ref_pic_num[LIST_1 + list_offset][0];
              enc_picture->mv         [LIST_1][block_y][block_x + i][0] = cur_mv[0];
              enc_picture->mv         [LIST_1][block_y][block_x + i][1] = cur_mv[1];
            }
          }
        }
        else 
        {
          for (j=0; j<BLOCK_MULTIPLE; j++)
          {
            block_y = img->block_y+(block&2)+j;
            block_x = img->block_x+(block&1)*2;
            memset(&enc_picture->ref_idx[LIST_1][block_y][block_x], best_bw_ref, BLOCK_MULTIPLE * sizeof(char));
            for (i=0; i<BLOCK_MULTIPLE; i++)
            {                     
              
              enc_picture->ref_pic_id [LIST_1][block_y][block_x + i] = 
                enc_picture->ref_pic_num[LIST_1 + list_offset][best_bw_ref];
              if(best_bw_ref>=0)
              {
                cur_mv = img->all_mv[j][i][LIST_1][best_bw_ref][mode];
                enc_picture->mv[LIST_1][block_y][block_x + i][0] = cur_mv[0];
                enc_picture->mv[LIST_1][block_y][block_x + i][1] = cur_mv[1];
              }
            }            
          }
        }
      }
    }
  }
  else if (mode==2)
  {
    for (j=0; j<2; j++)
    {
      block_y = img->block_y + block * 2 + j;
      for (i=0; i<BLOCK_MULTIPLE; i++)
      {
        block_x = img->block_x + i;
        if (best_pdir==1)
        {
          enc_picture->ref_idx    [LIST_0][block_y][block_x]    = -1;
          enc_picture->ref_pic_id [LIST_0][block_y][block_x]    = -1;
          enc_picture->mv         [LIST_0][block_y][block_x][0] = 0;
          enc_picture->mv         [LIST_0][block_y][block_x][1] = 0;
        }
        else
        {                     
          cur_mv = img->all_mv[j+block*2][i][LIST_0][best_fw_ref][mode];
          
          enc_picture->ref_idx    [LIST_0][block_y][block_x]    = best_fw_ref;
          enc_picture->ref_pic_id [LIST_0][block_y][block_x]    = 
            enc_picture->ref_pic_num[LIST_0 + list_offset][best_fw_ref];
          enc_picture->mv         [LIST_0][block_y][block_x][0] = cur_mv[0];
          enc_picture->mv         [LIST_0][block_y][block_x][1] = cur_mv[1];
        }
        
        if (bframe)
        {
          if (best_pdir==0)
          {
            enc_picture->ref_idx    [LIST_1][block_y][block_x]    = -1;
            enc_picture->ref_pic_id [LIST_1][block_y][block_x]    = -1;
            enc_picture->mv         [LIST_1][block_y][block_x][0] = 0;
            enc_picture->mv         [LIST_1][block_y][block_x][1] = 0;
          }
          else
          {
            enc_picture->ref_idx[LIST_1][block_y][block_x] = best_bw_ref;
            if(best_bw_ref>=0)
            {
              cur_mv = img->all_mv[j+ block*2][i][LIST_1][best_bw_ref][mode];
              
              enc_picture->ref_pic_id [LIST_1][block_y][block_x] = 
                enc_picture->ref_pic_num[LIST_1 + list_offset][best_bw_ref];
              enc_picture->mv[LIST_1][block_y][block_x][0] = cur_mv[0];
              enc_picture->mv[LIST_1][block_y][block_x][1] = cur_mv[1];
            }                       
          }
        }
      }
    }
  }
  else
  {
    for (j=0; j<BLOCK_MULTIPLE; j++)
    {
      block_y = img->block_y+j;
      for (i=0; i<2; i++)
      {
        block_x = img->block_x + block*2 + i;
        if (best_pdir==1)
        {
          enc_picture->ref_idx    [LIST_0][block_y][block_x]    = -1;
          enc_picture->ref_pic_id [LIST_0][block_y][block_x]    = -1;
          enc_picture->mv         [LIST_0][block_y][block_x][0] = 0;
          enc_picture->mv         [LIST_0][block_y][block_x][1] = 0;
        }
        else
        {
          cur_mv = img->all_mv[j][block*2+i][LIST_0][best_fw_ref][mode];
          
          enc_picture->ref_idx    [LIST_0][block_y][block_x] = best_fw_ref;
          enc_picture->ref_pic_id [LIST_0][block_y][block_x] = 
            enc_picture->ref_pic_num[LIST_0 + list_offset][best_fw_ref];          
          enc_picture->mv[LIST_0][block_y][block_x][0] = cur_mv[0];
          enc_picture->mv[LIST_0][block_y][block_x][1] = cur_mv[1];
        }
        
        if (bframe)
        {
          if (best_pdir==0)
          {
            enc_picture->ref_idx    [LIST_1][block_y][block_x]    = -1;
            enc_picture->ref_pic_id [LIST_1][block_y][block_x]    = -1;
            enc_picture->mv         [LIST_1][block_y][block_x][0] = 0;
            enc_picture->mv         [LIST_1][block_y][block_x][1] = 0;
          }
          else
          {
            enc_picture->ref_idx[LIST_1][block_y][block_x] = best_bw_ref;
            if(best_bw_ref>=0)
            {
              cur_mv = img->all_mv[j][block*2+i][LIST_1][best_bw_ref][mode];
              enc_picture->ref_pic_id [LIST_1][block_y][block_x] = 
                enc_picture->ref_pic_num[LIST_1 + list_offset][best_bw_ref];
              
              enc_picture->mv[LIST_1][block_y][block_x][0] = cur_mv[0];
              enc_picture->mv[LIST_1][block_y][block_x][1] = cur_mv[1];
            }
          }
        }
      }
    }
  }
}

void update_refresh_map(int intra, int intra1, Macroblock *currMB)
{
  if (input->RestrictRef==1)
  {
    // Modified for Fast Mode Decision. Inchoon Choi, SungKyunKwan Univ.
    if (input->rdopt<2)
    {
      refresh_map[2*img->mb_y  ][2*img->mb_x  ] = (intra ? 1 : 0);
      refresh_map[2*img->mb_y  ][2*img->mb_x+1] = (intra ? 1 : 0);
      refresh_map[2*img->mb_y+1][2*img->mb_x  ] = (intra ? 1 : 0);
      refresh_map[2*img->mb_y+1][2*img->mb_x+1] = (intra ? 1 : 0);
    }
    else if (input->rdopt==3)
    {
      refresh_map[2*img->mb_y  ][2*img->mb_x  ] = (intra1==0 && (currMB->mb_type==I16MB || currMB->mb_type==I4MB) ? 1 : 0);
      refresh_map[2*img->mb_y  ][2*img->mb_x+1] = (intra1==0 && (currMB->mb_type==I16MB || currMB->mb_type==I4MB) ? 1 : 0);
      refresh_map[2*img->mb_y+1][2*img->mb_x  ] = (intra1==0 && (currMB->mb_type==I16MB || currMB->mb_type==I4MB) ? 1 : 0);
      refresh_map[2*img->mb_y+1][2*img->mb_x+1] = (intra1==0 && (currMB->mb_type==I16MB || currMB->mb_type==I4MB) ? 1 : 0);
    }
  }
  else if (input->RestrictRef==2)
  {
    refresh_map[2*img->mb_y  ][2*img->mb_x  ] = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
    refresh_map[2*img->mb_y  ][2*img->mb_x+1] = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
    refresh_map[2*img->mb_y+1][2*img->mb_x  ] = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
    refresh_map[2*img->mb_y+1][2*img->mb_x+1] = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
  }
}  

void set_stored_macroblock_parameters1 ()
{
  Macroblock  *currMB  = &img->mb_data[img->current_mb_nr];
  memcpy(currMB->b8mode,b8mode, BLOCK_MULTIPLE * sizeof(int));
  memcpy(currMB->b8pdir,b8pdir, BLOCK_MULTIPLE * sizeof(int));
}
