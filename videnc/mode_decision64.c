#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include <memory.h>

#include "global.h"
#include "rdopt_coding_state.h"
#include "mb_access.h"
#include "intrarefresh.h"
#include "image.h"
#include "transform8x8.h"
#include "fast_me.h"
#include "simplified_fast_me.h"
#include "ratectl.h"            
#include "mode_decision.h"
#include "fmo.h"

#ifdef MV_COMPETITION
#include "mv_competition.h"
extern MV_Competition mv_comp;
extern int skip_mode;
#endif

#define SIMPLIFY_CODE

#ifdef MB32X32

#define IS_INTER_SLICE(type) ((type) == B_SLICE || (type) == P_SLICE) //also defined in rdopt.c

void
PartitionMotionSearch32 (int    blocktype,
                       int    block8x8, //MB32x32: block8x8=32x32, 32x16, 16x32, 16x16; MB16x16: block8x8=16x16, 16x8, 8x16, 8x8 
                       int    lambda_factor,
#ifdef RDO_Q
                       int   transform8x8,
#endif
                       int mb_ext_level
                       );
void FindSkipModeMotionVector32 ();
void Get_Direct_Motion_Vectors32 (int mb_ext_level);
void CheckAvailabilityOfNeighborsCABAC(void);
int BPredPartitionCost32 (int   blocktype,
                        int   block8x8,
                        short fw_ref,
                        short bw_ref,
                        int   lambda_factor,
                        int   list,
                        int   mb_ext_level);
int BIDPartitionCost32 (int   blocktype,
                      int   block8x8,
                      short fw_ref,
                      short bw_ref,
                      int   lambda_factor,
                      int   mb_ext_level);
void find_Mb_cluster(int level, int CurrentMbAddr, int *CurrentMbAddr_x, int *CurrentMbAddr_y, int *extend_x, int *extend_y, int *num_MB, int *nextMbAddr);
void assign_qp_info(int mbcount, int encodeMbAddr, int saved_cluster_qp, int saved_prev_qp, int saved_prev_delta_qp);
int writeMBLayer32 (int rdopt, int *coeff_rate, int mb32, int mbi);
void store_adaptive_rounding_parameters (int mode, Macroblock *currMB);
void update_offset_params(int mode, int luma_transform_size_8x8_flag);


#ifdef ADAPTIVE_FD_SD_CODING
extern int   best_SD_or_FD           [2][2];
extern int   best_SD_or_FD_t8x8;    
extern int   best_quantizer_indices  [16][16];
extern int   best_SD_Coding_on_off;
#endif
#ifdef USE_INTRA_MDDT
extern long quant_stat_rd[16];
extern long quant_stat_rd16x16[256];

extern int   bestCofAC16[2][16*17];
#endif 
#ifdef ADAPTIVE_QUANTIZATION
extern int   best_mb_iaqms_idx;
#endif


int   cbp16[1<<(2*MAX_MB_EXT_LEVEL)], cbp16_RDCost[1<<(2*MAX_MB_EXT_LEVEL)];
int   cbp_32[1<<(2*(MAX_MB_EXT_LEVEL-1))], cbp_32_RDCost[1<<(2*(MAX_MB_EXT_LEVEL-1))];
int64 cbp_blk16[1<<(2*MAX_MB_EXT_LEVEL)], cbp_blk16_RDCost[1<<(2*MAX_MB_EXT_LEVEL)];
int   ****cofAC32[1<<(2*MAX_MB_EXT_LEVEL)], ****cofAC32_RDCost[1<<(2*MAX_MB_EXT_LEVEL)];// ****cofAC8x8=NULL;        // [8x8block][4x4block][level/run][scan_pos]
int   ***cofDC32[1<<(2*MAX_MB_EXT_LEVEL)], ***cofDC32_RDCost[1<<(2*MAX_MB_EXT_LEVEL)];                       // [yuv][level/run][scan_pos]
int    **cof32AC16x16[1<<(2*MAX_MB_EXT_LEVEL)], **cof32AC16x16_RDCost[1<<(2*MAX_MB_EXT_LEVEL)];//[2][257];
int    ***cof32AC16x8[1<<(2*MAX_MB_EXT_LEVEL)], ***cof32AC16x8_RDCost[1<<(2*MAX_MB_EXT_LEVEL)];//[2][2][129];
int    ***cof32AC8x16[1<<(2*MAX_MB_EXT_LEVEL)], ***cof32AC8x16_RDCost[1<<(2*MAX_MB_EXT_LEVEL)];//[2][2][129];
int only_16x16_transform=0;
int   luma_transform_size_8x8_flag32[1<<(2*MAX_MB_EXT_LEVEL)], luma_transform_size_8x8_flag32_RDCost[1<<(2*MAX_MB_EXT_LEVEL)];
char  frefframe32[4<<MAX_MB_EXT_LEVEL][4<<MAX_MB_EXT_LEVEL], brefframe32[4<<MAX_MB_EXT_LEVEL][4<<MAX_MB_EXT_LEVEL];
extern char  frefframe[4][4], brefframe[4][4];
int   best16x16pdir[4];
int bestmode32;
extern int best_c_imode;
extern int best_i16offset;
//extern int   ****cofAC=NULL, ****cofAC8x8=NULL;        // [8x8block][4x4block][level/run][scan_pos]
//extern int   ***cofDC=NULL;                       // [yuv][level/run][scan_pos]
extern int   b8mode[4], b8pdir[4];
extern int   luma_transform_size_8x8_flag;
extern  int ****cofAC;               //!< AC coefficients [8x8block][4x4block][level/run][scan_pos]
extern  int ***cofDC;                //!< DC coefficients [yuv][level/run][scan_pos]
extern int write_mb32_mbtype;
Macroblock  MB32;

//#ifdef MB64X64
Macroblock  MB64;
int bestmode64;
int mb32type[1<<(2*(MAX_MB_EXT_LEVEL-1))];
int write_mb64_mbtype;
int RDCost_for_macroblocks64 (double   lambda,       // <-- lagrange multiplier
                            int      mode,         // <-- modus (0-COPY/DIRECT, 1-16x16, 2-16x8, 3-8x16, 4-8x8(+), 5-Intra4x4, 6-Intra16x16)
                            double*  min_rdcost,   // <-> minimum rate-distortion cost
                            double*  min_rate,     // --> bitrate of mode which has minimum rate-distortion cost. 
                            int i16mode);
int writeMBLayer64 (int rdopt, int *coeff_rate, int mb32, int mb16);
//#endif

Macroblock MB32_16x16;//32x32 block, 16x16 mb type
int saved_mb32_nr;

#ifdef RDO_Q
extern int curr_mbi, curr_mb32i;
#endif
int Not_MB16_3(int mb_nr)
{
  int mb_nr_x, mb_nr_y;
  mb_nr_x = mb_nr % img->PicWidthInMbs;
  mb_nr_y = mb_nr / img->PicWidthInMbs;

  if(mb_nr_x%2==1 && mb_nr_y%2==1)
  {
    if(mb_nr_x < (signed)(img->PicWidthInMbs-1))
    {
      return 0;
    }
  }
  return 1;
}
void CheckAvailabilityOfNeighbors_Z(void)
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
    currMB->mbAddrC = mb_nr - img->PicWidthInMbs + 1;
    currMB->mbAddrD = mb_nr - img->PicWidthInMbs - 1;

    currMB->mbAvailA = mb_is_available(currMB->mbAddrA, mb_nr) && ((mb_nr % img->PicWidthInMbs)!=0);
    currMB->mbAvailB = mb_is_available(currMB->mbAddrB, mb_nr);
    currMB->mbAvailC =  Not_MB16_3(mb_nr) && mb_is_available(currMB->mbAddrC, mb_nr) && (((mb_nr+1) % img->PicWidthInMbs)!=0);
    currMB->mbAvailD = mb_is_available(currMB->mbAddrD, mb_nr) && ((mb_nr % img->PicWidthInMbs)!=0);
  }

  if (currMB->mbAvailA) currMB->mb_available_left = &(img->mb_data[currMB->mbAddrA]);
  if (currMB->mbAvailB) currMB->mb_available_up   = &(img->mb_data[currMB->mbAddrB]);
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
*    initializes the current macroblock
*
* \param 
*    mb32: 
*       2-initialization for 64x64 block;
*       1-initialization for 32x32 block;
*       0-initialization for a 16x16 block inside a 32x32 block
*    mb_addr: macroblock address in scan order
*    if mb32 is 1, mb_addr is the upper left 16x16 block index in scan order
*    if mb32 is 0, mb_addr is the 16x16 block index in scan order
* \param mb_field
*    true for field macroblock coding
************************************************************************
*/
void start_macroblock32(int mb_addr, int mb_field, int mb32)
{
  int i,j,l;
  Macroblock *currMB = &img->mb_data[mb_addr];
  Slice *curr_slice = img->currentSlice;
  //DataPartition *dataPart;
  //Bitstream *currStream;
  //EncodingEnvironmentPtr eep;
  //int max_qp_delta = 25 + img->bitdepth_luma_qp_scale/2;
  //int min_qp_delta = (26 + img->bitdepth_luma_qp_scale/2);
  int prev_mb;
  
  currMB->mb_field = mb_field;
  
  enc_picture->mb_field[mb_addr] = mb_field;
  
  set_MB_parameters (mb_addr); //get the upper left corner coordinators
  
  prev_mb = FmoGetPreviousMBNr(img->current_mb_nr);
  
  
  // Save the slice number of this macroblock. When the macroblock below
  // is coded it will use this to decide if prediction for above is possible
  currMB->slice_nr = img->current_slice_nr;
  

  
  if(mb32)
  {
    if (prev_mb>-1)
    {
      currMB->prev_qp = img->mb_data[prev_mb].qp;
      if (img->mb_data[prev_mb].slice_nr == img->current_slice_nr)
      {
        currMB->prev_delta_qp = img->mb_data[prev_mb].delta_qp; //for delta_qp encoding in cabac
      }
      else
      {
        currMB->prev_delta_qp = 0;
      }
    }
    else
    {
      currMB->prev_qp = curr_slice->qp;
      currMB->prev_delta_qp = 0;
    }
    
#ifndef RDO_Q
    currMB->qp = currSlice->qp;
#else
    if(input->UseRDO_Q)
      currMB->qp = img->qp ; 
    else
      currMB->qp = curr_slice->qp;    
#endif
    
  }
  else
  {
    //the following will be overwritten after this function exits when MB32_DELTA_QP is defined
    currMB->prev_qp = MB32.prev_qp;
    currMB->prev_delta_qp = MB32.prev_delta_qp; //for delta_qp encoding in cabac
    currMB->qp = MB32.qp;
  }
  currMB->delta_qp = currMB->qp - currMB->prev_qp;
  //DELTA_QP = DELTA_QP2 = currMB->delta_qp;
  //QP = QP2 = currMB->qp;
  
  set_chroma_qp (currMB);
  
  // Initialize counter for MB symbols
  currMB->currSEnr=0;
  
  // loop filter parameter
  if (active_pps->deblocking_filter_control_present_flag)
  {
    currMB->LFDisableIdc    = img->LFDisableIdc;
    currMB->LFAlphaC0Offset = img->LFAlphaC0Offset;
    currMB->LFBetaOffset    = img->LFBetaOffset;
  }
  else
  {
    currMB->LFDisableIdc    = 0;
    currMB->LFAlphaC0Offset = 0;
    currMB->LFBetaOffset    = 0;
  }
  
  // If MB is next to a slice boundary, mark neighboring blocks unavailable for prediction
  CheckAvailabilityOfNeighbors32(mb32);
  
  if (input->symbol_mode == CABAC)
    CheckAvailabilityOfNeighborsCABAC();
  
  if(mb32)
  {
    // Reset vectors and reference indices before doing motion search in motion_search().
    for (l=0; l<2; l++)
    {
      for (j=img->block_y; j < img->block_y + (BLOCK_MULTIPLE<<mb32); j++)
      {
        memset(&enc_picture->ref_idx[l][j][img->block_x], -1, (BLOCK_MULTIPLE<<mb32) * sizeof(char));     
        memset(enc_picture->mv [l][j][img->block_x], 0, 2 * (BLOCK_MULTIPLE<<mb32) * sizeof(short));
        for (i=img->block_x; i < img->block_x + (BLOCK_MULTIPLE<<mb32); i++)
          enc_picture->ref_pic_id[l][j][i]= -1;
      }
    }
  }
  // Reset syntax element entries in MB struct
  currMB->mb_type      = 0;
  currMB->cbp_blk      = 0;
  currMB->cbp          = 0;
  currMB->cbp_bits     = 0;
#ifdef ADAPTIVE_FD_SD_CODING
  currMB->FD_or_SD_bits= 0;
#endif
#ifdef ADAPTIVE_QUANTIZATION
  currMB->mb_iaqms_idx = img->mb_iaqms_idx;
#endif
  currMB->c_ipred_mode = DC_PRED_8;
  
  memset (currMB->mvd, 0, BLOCK_CONTEXT * sizeof(int));  
  memset (currMB->intra_pred_modes, DC_PRED, MB_BLOCK_PARTITIONS * sizeof(char)); // changing this to char would allow us to use memset
  memset (currMB->intra_pred_modes8x8, DC_PRED, MB_BLOCK_PARTITIONS * sizeof(char));
  
  
  // Initialize bitcounters for this macroblock
  if(prev_mb < 0) // No slice header to account for
  {
    currMB->bitcounter[BITS_HEADER] = 0;
  }
  else if (currMB->slice_nr == img->mb_data[prev_mb].slice_nr) // current MB belongs to the
    // same slice as the last MB
  {
    currMB->bitcounter[BITS_HEADER] = 0;
  }
  
  currMB->bitcounter[BITS_MB_MODE       ] = 0;
  currMB->bitcounter[BITS_COEFF_Y_MB    ] = 0;
  currMB->bitcounter[BITS_INTER_MB      ] = 0;
  currMB->bitcounter[BITS_CBP_MB        ] = 0;
  currMB->bitcounter[BITS_DELTA_QUANT_MB] = 0;
  currMB->bitcounter[BITS_COEFF_UV_MB   ] = 0;
#ifdef MB32X32_MVC
  currMB->bitcounter[BITS_SKIP_PRED     ] = 0;
  currMB->bitcounter[BITS_MV_PRED       ] = 0;
#endif
  
#ifdef _FAST_FULL_ME_
  if(!input->FMEnable)
    ResetFastFullIntegerSearch ();
#endif
}
/*!
*************************************************************************************
* \brief
*    Sets modes and reference frames for a 16x16 macroblock
*************************************************************************************
*/
void SetModesAndRefframeForBlocks32 (int mode, short best16x16pdir, short best16x16fwref, short best16x16bwref)
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
      if(i>0)
      {
        assert(currMB->b8pdir[i] == currMB->b8pdir[i-1]);
      }
    }
    break;
  case 1:
  case 2:
  case 3:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] = mode;
      currMB->b8pdir[i] = best16x16pdir;
    }
    break;

  default:
    printf ("Unsupported mode in SetModesAndRefframeForBlocks32!\n");
    exit (1);
  }
#define IS_FW ((best16x16pdir==0 || best16x16pdir==2))
#define IS_BW ((best16x16pdir==1 || best16x16pdir==2))
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
#ifndef SIMPLIFY_CODE
      else
      {
        for (j = img->block_y; j < img->block_y + 4; j++)
        {
          memset(&enc_picture->ref_idx[LIST_0][j][img->block_x],-1, 4 * sizeof(char));
          memset(&enc_picture->ref_idx[LIST_1][j][img->block_x],-1, 4 * sizeof(char));
        }
      }
#endif
    }
    else
    {
      if (!mode)
      {
        for (j = img->block_y; j < img->block_y + 4; j++)
          memset(&enc_picture->ref_idx[LIST_0][j][img->block_x],0, 4 * sizeof(char));
      }
#ifndef SIMPLIFY_CODE
      else
      {
        for (j = img->block_y; j < img->block_y + 4; j++)
          memset(&enc_picture->ref_idx[LIST_0][j][img->block_x],-1, 4 * sizeof(char));
      }
#endif
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
          
#ifndef SIMPLIFY_CODE
          if(mode == P8x8 && best8x8mode[k]==0)
          {
            enc_picture->ref_idx[LIST_0][block_y][block_x] = direct_ref_idx[LIST_0][block_y][block_x];
            enc_picture->ref_idx[LIST_1][block_y][block_x] = direct_ref_idx[LIST_1][block_y][block_x];
          }
          else
#endif
            if (mode ==1 && currMB->bi_pred_me && IS_FW && IS_BW)
          {
            enc_picture->ref_idx[LIST_0][block_y][block_x] = 0;
            enc_picture->ref_idx[LIST_1][block_y][block_x] = 0;
          }
          else
          {
            enc_picture->ref_idx[LIST_0][block_y][block_x] = (IS_FW ? best16x16fwref : -1);
            enc_picture->ref_idx[LIST_1][block_y][block_x] = (IS_BW ? best16x16bwref : -1);
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
          enc_picture->ref_idx[LIST_0][block_y][block_x] = (IS_FW ? best16x16fwref : -1);
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

void SetModesAndRefframeForBlocks32_notused (int mode)
{
  int i,j,k;
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
      currMB->b8pdir[i] = (bframe ? direct_pdir[img->block_y + (i >> 1)*4][img->block_x + (i & 0x01)*4] : 0);
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
#ifndef SIMPLIFY_CODE
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
#endif
  default:
    printf ("Unsupported mode in SetModesAndRefframeForBlocks32!\n");
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
        for (j = img->block_y; j < img->block_y + 8; j++)
        {
          memcpy(&enc_picture->ref_idx[LIST_0][j][img->block_x],&direct_ref_idx[LIST_0][j][img->block_x], 8 * sizeof(char));
          memcpy(&enc_picture->ref_idx[LIST_1][j][img->block_x],&direct_ref_idx[LIST_1][j][img->block_x], 8 * sizeof(char));
        }
      }
#ifndef SIMPLIFY_CODE
      else
      {
        for (j = img->block_y; j < img->block_y + 8; j++)
        {
          memset(&enc_picture->ref_idx[LIST_0][j][img->block_x],-1, 8 * sizeof(char));
          memset(&enc_picture->ref_idx[LIST_1][j][img->block_x],-1, 8 * sizeof(char));
        }
      }
#endif
    }
    else
    {
      if (!mode)
      {
        for (j = img->block_y; j < img->block_y + 8; j++)
          memset(&enc_picture->ref_idx[LIST_0][j][img->block_x],0, 8 * sizeof(char));
      }
      else
      { error("should not get here\n", -1);
        for (j = img->block_y; j < img->block_y + 4; j++)
          memset(&enc_picture->ref_idx[LIST_0][j][img->block_x],-1, 4 * sizeof(char));
      }
    }
  }
  else
  {
    if (bframe)
    {
      for (j=0;j<8;j++)
      {
        block_y = img->block_y + j;
        for (i=0;i<8;i++)
        {
          block_x = img->block_x + i;
          k = 2*(j >> 2) + (i >> 2);
#ifndef SIMPLIFY_CODE
          l = 2*(j & 0x01) + (i & 0x01);

          if(mode == P8x8 && best8x8mode[k]==0)
          {
            enc_picture->ref_idx[LIST_0][block_y][block_x] = direct_ref_idx[LIST_0][block_y][block_x];
            enc_picture->ref_idx[LIST_1][block_y][block_x] = direct_ref_idx[LIST_1][block_y][block_x];
          }
          else 
#endif
            if (mode ==1 && currMB->bi_pred_me && IS_FW && IS_BW)
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
      for (j=0;j<8;j++)
      {
        block_y = img->block_y + j;
        for (i=0;i<8;i++)
        {
          block_x = img->block_x + i;
          k = 2*(j >> 2) + (i >> 2);

#ifndef SIMPLIFY_CODE
          l = 2*(j & 0x01) + (i & 0x01);
#endif
          enc_picture->ref_idx[LIST_0][block_y][block_x] = (IS_FW ? best8x8fwref[mode][k] : -1);
        }
      }
    }
  }
  
  if (bframe)
  {    
    for (j = img->block_y; j < img->block_y + 8; j++)
      for (i = img->block_x; i < img->block_x + 8;i++)
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
    for (j = img->block_y; j < img->block_y + 8; j++)
      for (i = img->block_x; i < img->block_x + 8;i++)
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

void switchback_data_after_write(int encodeMbAddr, int mbcount, int mode)
{
  int ****i4p, ***i3p;
  int **i2p;
  i4p=cofAC32_RDCost[mbcount]; cofAC32_RDCost[mbcount]=img->cofAC; img->cofAC=i4p;
  i3p=cofDC32_RDCost[mbcount]; cofDC32_RDCost[mbcount]=img->cofDC; img->cofDC=i3p;

  i2p=cof32AC16x16_RDCost[mbcount]; cof32AC16x16_RDCost[mbcount]=img->cofAC16x16; img->cofAC16x16=i2p;
}

void prepare_data_for_write_RDCost(int encodeMbAddr, int mbcount, int mode)
{
  int ****i4p, ***i3p;
  Macroblock  *currMB  = &img->mb_data[encodeMbAddr];
  int **i2p;

  currMB->cbp = cbp16_RDCost[mbcount];
  currMB->cbp_blk = cbp_blk16_RDCost[mbcount];
  currMB->mb_type = mode;
  //currMB->luma_transform_size_8x8_flag = luma_transform_size_8x8_flag32_RDCost[mbcount];

  i4p=cofAC32_RDCost[mbcount]; cofAC32_RDCost[mbcount]=img->cofAC; img->cofAC=i4p;
  i3p=cofDC32_RDCost[mbcount]; cofDC32_RDCost[mbcount]=img->cofDC; img->cofDC=i3p;

  i2p=cof32AC16x16_RDCost[mbcount]; cof32AC16x16_RDCost[mbcount]=img->cofAC16x16; img->cofAC16x16=i2p;
}

void prepare_mv(int mbi_x, int mbi_y)
{
  int i, j, list, ref, mode,k;
  short******i6p;

  for(i=0; i<4; i++)
    for(j=0; j<4; j++)
      for(list=0; list<2; list++)
        for(ref=0; ref<img->max_num_references; ref++)
          for(mode=0; mode<9; mode++)
            for(k=0; k<2; k++)                    
              img->all_mv32[i][j][list][ref][mode][k] = img->all_mv[mbi_y*4+i][mbi_x*4+j][list][ref][mode][k];

  i6p=img->all_mv; img->all_mv=img->all_mv32; img->all_mv32=i6p;

}
void switch_back_mv()
{
  short******i6p;
  i6p=img->all_mv; img->all_mv=img->all_mv32; img->all_mv32=i6p;
}

/*!
*************************************************************************************
* \brief
*    R-D Cost for a macroblock
*************************************************************************************
*/
int RDCost_for_macroblocks32 (double   lambda,       // <-- lagrange multiplier
                            int      mode,         // <-- modus (0-COPY/DIRECT, 1-16x16, 2-16x8, 3-8x16, 4-8x8(+), 5-Intra4x4, 6-Intra16x16)
                            double*  min_rdcost,   // <-> minimum rate-distortion cost
                            double*  min_rate,     // --> bitrate of mode which has minimum rate-distortion cost. 
                            int i16mode,
                            int mb_ext_level)
{
  
  int         i, j, mbi, mbj; 
  int         j1, j2;
  int         rate = 0, coeff_rate = 0;
  int64       distortion = 0;
  double      rdcost;
  Macroblock  *currMB;
  //Macroblock  *prevMB;
  //int         bframe    = (img->type==B_SLICE);
  //int         use_of_cc =  (img->type!=I_SLICE &&  input->symbol_mode!=CABAC);
  int         dummy;
  int      return_val = 0;
  //int tmp_cc;
  //int cc_rate;
  int saved_luma_transform_size_8x8_flag;
  int encodeMbAddr, CurrentMbAddr_x, CurrentMbAddr_y, extend_x, extend_y, num_MB, nextMbAddr, mb_count=0;
  int ****i4p, ***i3p;
  int **i2p;
#ifdef MB32X32_MVC
  int         cr_cbp = 0;
  short iter, nbiter;
  short ******all_mv = img->all_mv;
  short loopbreak=0;
#endif

  if(mb_ext_level == 2)//MB64X64
  {    
    return (RDCost_for_macroblocks64 (lambda, mode, min_rdcost, min_rate, i16mode));
  }

#ifdef MB32X32_MVC
  if ((input->mv_competition > 0) && (mode == 0) && (img->type==P_SLICE))
    nbiter = mv_comp.nb_mode_for_skip + 1;
  else
    nbiter = 1;
  
  // Let's iterate on each possible mode for the SKIP
  for (iter= 0; iter < nbiter; iter++)
  {
    loopbreak = 0;
    if ((input->mv_competition > 0) && (mode == 0) && (img->type==P_SLICE))  
    {
      distortion = 0;
      coeff_rate = 0;
      rate = 0;
      cr_cbp = 0;
    }
#endif

  /*******************************
   * restore some global variables
   *******************************/
  set_MB_parameters(saved_mb32_nr);
  saved_luma_transform_size_8x8_flag = img->mb_data[img->current_mb_nr].luma_transform_size_8x8_flag;


  //=====   S T O R E   C O D I N G   S T A T E   =====
  //---------------------------------------------------
  store_coding_state (cs_cm);
    
  
  find_Mb_cluster(mb_ext_level, img->current_mb_nr, &CurrentMbAddr_x, &CurrentMbAddr_y, &extend_x, &extend_y, &num_MB, &nextMbAddr);

  mb_count=0;
  for(mbi=CurrentMbAddr_y; mbi<=extend_y; mbi++)  
  {      
    for(mbj=CurrentMbAddr_x; mbj<=extend_x; mbj++)   
    {             
      encodeMbAddr = mbi*img->PicWidthInMbs + mbj;
      
      start_macroblock32 (encodeMbAddr, FALSE, 0);//reset img->current_mb_nr and other global variables

      currMB   = &img->mb_data[img->current_mb_nr];   
      //prevMB   = img->current_mb_nr ? &img->mb_data[img->current_mb_nr-1] : NULL;
      //=====
      //=====  SET REFERENCE FRAMES AND BLOCK MODES
      //=====    
      SetModesAndRefframeForBlocks32 (mode, best8x8pdir[mode][mb_count], best8x8fwref[mode][mb_count], best8x8bwref[mode][mb_count]);
      best16x16pdir[mb_count] = currMB->b8pdir[0];

      mb_count++;
    }
  }

  mb_count=0;
  for(mbi=CurrentMbAddr_y; mbi<=extend_y; mbi++)  
  {      
    for(mbj=CurrentMbAddr_x; mbj<=extend_x; mbj++)   
    {             
      encodeMbAddr = mbi*img->PicWidthInMbs + mbj;
      
      start_macroblock32 (encodeMbAddr, FALSE, 0);//reset img->current_mb_nr and other global variables

      currMB   = &img->mb_data[img->current_mb_nr];   
      //prevMB   = img->current_mb_nr ? &img->mb_data[img->current_mb_nr-1] : NULL;
      //=====
      //=====  SET REFERENCE FRAMES AND BLOCK MODES
      //=====    
      if(mb_count>0)
      {
        img->mb_data[encodeMbAddr].luma_transform_size_8x8_flag = saved_luma_transform_size_8x8_flag;
      }
      SetModesAndRefframeForBlocks32 (mode, best8x8pdir[mode][mb_count], best8x8fwref[mode][mb_count], best8x8bwref[mode][mb_count]);

    

      //=====
      //=====  GET COEFFICIENTS, RECONSTRUCTIONS, CBP
      //=====
#ifndef SIMPLIFY_CODE
      if (bframe && mode==0)
      {
        int block_x=img->pix_x>>2;
        int block_y=img->pix_y>>2;    
        for (j = block_y;j< block_y + 4;j++)
          for (i = block_x;i<block_x + 4;i++)
            if (direct_pdir[j][i] < 0) // not necessary in the spatial direct mode
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
#endif
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
      

#ifdef MB32X32_MVC
      if ((input->mv_competition > 0) && (mode == 0) && (img->type==P_SLICE))
      {
        // The last iteration consist of re-launching with the bestmode,
        // to set needed variables for the JM
        if (iter == (nbiter - 1))
          mv_comp.current_predictor_for_skip = MB32.best_predictor_for_skip;
        else
        {
          mv_comp.current_predictor_for_skip = iter;
          if ((mv_comp.mv_pred_skip[mv_comp.current_predictor_for_skip][0] == NOT_AVAILABLE) || (mv_comp.mv_pred_skip[mv_comp.current_predictor_for_skip][1] == NOT_AVAILABLE))
          {
            loopbreak=1;
            goto nextinloop;
            //continue;
          }
        }

        prepare_mv((mb_count&1), (mb_count>>1));
        only_16x16_transform = (currMB->luma_transform_size_8x8_flag == 2);
        LumaResidualCoding ();
        only_16x16_transform = 0;
        switch_back_mv();

        prepare_mv((mb_count&1), (mb_count>>1));
        img->i16offset = 0;
        dummy = 0;
        if (img->yuv_format!=YUV400)
          ChromaResidualCoding (&dummy);
        switch_back_mv();      
      }

      else //if((input->mv_competition > 0) && (mode == 0) && (img->type==P_SLICE))
      {
#endif
        if (mode<P8x8)
        {
          prepare_mv((mb_count&1), (mb_count>>1));

          only_16x16_transform = currMB->luma_transform_size_8x8_flag == 2;

          LumaResidualCoding ();

          only_16x16_transform = 0;

          switch_back_mv();

        }
#ifndef SIMPLIFY_CODE
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
            if(Motion_Selected && IS_INTER_SLICE(img->type))
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
            if(Motion_Selected && IS_INTER_SLICE(img->type))
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
            if(Motion_Selected && IS_INTER_SLICE(img->type))
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
#endif


#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
#if 1 //def TRELLIS_ONLY
        if(0)
#else
        if(input->UseRDO_Q || img->slice_fractional_quant_flag)
#endif
#else
        if(input->UseRDO_Q)
#endif
        {
          if (img->type == I_SLICE  || !IS_INTRA_MODE(mode) || (IS_INTRA_MODE(mode) && !Motion_Selected && IS_INTER_SLICE(img->type)))
          { 
#ifndef SIMPLIFY_CODE
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
#endif          



            img->i16offset = 0;
            dummy = 0;
            if ((!(img->residue_transform_flag && (mode==I4MB || mode==I16MB || mode==I8MB)))
              && (img->yuv_format!=YUV400) && (mode != IPCM))
              ChromaResidualCoding (&dummy);
            
#ifndef SIMPLIFY_CODE
            if (mode==I16MB)     
              img->i16offset = I16Offset  (currMB->cbp, i16mode);
            
            if(IS_INTER_SLICE(img->type) && IS_INTRA_MODE(mode) && !Motion_Selected)
              store_macroblock_parameters_intra(mode); 
#endif
          }
        }
        else
        {
#ifndef SIMPLIFY_CODE
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
#endif

          prepare_mv((mb_count&1), (mb_count>>1));

          img->i16offset = 0;
          dummy = 0;
          if ((!(img->residue_transform_flag && (mode==I4MB || mode==I16MB || mode==I8MB)))
            && (img->yuv_format!=YUV400) && (mode != IPCM))
            ChromaResidualCoding (&dummy);

          switch_back_mv();


#ifndef SIMPLIFY_CODE        
          if (mode==I16MB)     
            img->i16offset = I16Offset  (currMB->cbp, i16mode);
#endif
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
#ifdef MB32X32_MVC
      }
#endif
        

#ifndef SIMPLIFY_CODE
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
      if( (input->UseRDO_Q || img->slice_fractional_quant_flag) && IS_INTRA_MODE(mode) && Motion_Selected && IS_INTER_SLICE(img->type))
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
#endif      
        //=====
        //=====   GET DISTORTION
        //=====
        // LUMA
#ifndef SIMPLIFY_CODE
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
#endif
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
#ifndef SIMPLIFY_CODE
        //=====
        //=====   GET RATE
        //=====
        //----- macroblock header -----
        if (use_of_cc) //use_of_cc =  (img->type!=I_SLICE &&  input->symbol_mode!=CABAC);
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
#endif

        /******************************************
         *for the upper left 16x16 block, set b8mode
         *b8pdir for the whole 32x32 block in order
         *to write motion info in 1 shot
         *****************************************/
        //if(mb_count == 0)
        //{
        //  int i;
        //  for(i=0;i<4;i++)
        //  {
        //    currMB->b8mode[i] = mode;
        //    currMB->b8pdir[i] = best16x16pdir[i];
        //  }    
        //}
//      SetModesAndRefframeForBlocks32 (mode, best8x8pdir[mode][mb_count], best8x8fwref[mode][mb_count], best8x8bwref[mode][mb_count]);
        //rate += writeMBLayer32 (1, &coeff_rate, 1, mb_count);


        /******************************************
         *save macroblock info for each 16x16 block
         *****************************************/
        i4p=cofAC32_RDCost[mb_count]; cofAC32_RDCost[mb_count]=img->cofAC; img->cofAC=i4p;  
        i3p=cofDC32_RDCost[mb_count]; cofDC32_RDCost[mb_count]=img->cofDC; img->cofDC=i3p;  
        i2p=cof32AC16x16_RDCost[mb_count]; cof32AC16x16_RDCost[mb_count]=img->cofAC16x16; img->cofAC16x16=i2p;
        cbp16_RDCost[mb_count] = currMB->cbp;  
        cbp_blk16_RDCost[mb_count] = currMB->cbp_blk;  
  
        if (((currMB->cbp & 15) == 0) && !(IS_OLDINTRA(currMB) || currMB->mb_type == I8MB))//originally done in set_stored_macroblock_parameters
          luma_transform_size_8x8_flag32_RDCost[mb_count] = 0;
        else
          luma_transform_size_8x8_flag32_RDCost[mb_count] = currMB->luma_transform_size_8x8_flag;
        
        mb_count++;
      } // j             
   }//i
      
#ifdef MB32X32_MVC
nextinloop:
  if(loopbreak==1)
    continue;   //for (iter= 0; iter < nbiter; iter++)
#endif

  //prepare 32x32 block cbp     
   MB32.cbp=0;
  for(mbi=0; mbi<4; mbi++)
  {
    if(cbp16_RDCost[mbi] != 0)
    {
      MB32.cbp = 1;
      break;
    }
  }      
  MB32.mb_type = mode;

  mb_count=0;
  for(mbi=CurrentMbAddr_y; mbi<=extend_y; mbi++)  
  {      
    for(mbj=CurrentMbAddr_x; mbj<=extend_x; mbj++)   
    {             
      encodeMbAddr = mbi*img->PicWidthInMbs + mbj;
      
      start_macroblock32 (encodeMbAddr, FALSE, 0);//reset img->current_mb_nr and other global variables

      currMB   = &img->mb_data[img->current_mb_nr];   
 
        /******************************************
         *for the upper left 16x16 block, set b8mode
         *b8pdir for the whole 32x32 block in order
         *to write motion info in 1 shot
         * mv and ref_idx are set in SetModesAndRefframeForBlocks32
         *****************************************/
        if(mb_count == 0)
        {
          int i;
          for(i=0;i<4;i++)
          {
            currMB->b8mode[i] = mode;
            currMB->b8pdir[i] = best16x16pdir[i];
          }    
        }

        prepare_data_for_write_RDCost(encodeMbAddr, mb_count, mode);
        img->mb_data[encodeMbAddr].luma_transform_size_8x8_flag = saved_luma_transform_size_8x8_flag;

        rate += writeMBLayer32 (1, &coeff_rate, 1, mb_count);
        switchback_data_after_write(encodeMbAddr, mb_count, mode);

        mb_count++;
      } // j       
    }//i
      
      
      //=====   R E S T O R E   C O D I N G   S T A T E   =====
      //-------------------------------------------------------
      //restore some global variables
      set_MB_parameters(saved_mb32_nr);
      //restore some global variables

      reset_coding_state (cs_cm);

      
      rdcost = (double)distortion + lambda * max(0.5,(double)rate);


      //if(mb_ext_level==1 && img->current_mb_nr==1836)
      //{
      //  int dist = (int)distortion;
      //  printf("mb_nr %d mode %4d dist %6d, rate %3d, cost %6d\n", saved_mb32_nr, mode, dist, rate, (int)rdcost);
      //}

#ifndef SIMPLIFY_CODE
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
      if( (input->UseRDO_Q || img->slice_fractional_quant_flag) && IS_INTRA_MODE(mode) && IS_INTER_SLICE(img->type))
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
#endif    
    if (!(rdcost >= *min_rdcost ||
      ((img->qp + img->bitdepth_luma_qp_scale)==0 && img->lossless_qpprime_flag==1 && distortion!=0)))
    {               
#ifndef SIMPLIFY_CODE
      if (((img->MbaffFrameFlag) && (mode ? 0: ((img->type == B_SLICE) ? !currMB->cbp:1)))  // AFF and current is skip
        && (img->current_mb_nr & 0x01) //bottom
        && (prevMB->mb_type ? 0:((img->type == B_SLICE) ? !prevMB->cbp:1)) //top is skip
        && (!(field_flag_inference() == currMB->mb_field)))
        return_val = 0;
      else      
#endif
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
#ifdef MB32X32_MVC
          if ((input->mv_competition > 0) && (mode == 0) && (img->type==P_SLICE))
          {
            int by, bx;
            
            MB32.best_predictor_for_skip = mv_comp.current_predictor_for_skip;
            
            for (by = 0;by < 4<<mb_ext_level;by++)
            {
              for (bx = 0;bx < 4<<mb_ext_level;bx++)
              {
                all_mv [by][bx][0][0][0][0] = mv_comp.mv_pred_skip[MB32.best_predictor_for_skip][0];
                all_mv [by][bx][0][0][0][1] = mv_comp.mv_pred_skip[MB32.best_predictor_for_skip][1];
              }
            }
          }
#endif
          
          return_val = 1;
      }
      
      
    }
    else 
    {
#ifdef MB32X32_MVC
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
#ifdef MB32X32_MVC
      }
#endif
    }
#ifdef MB32X32_MVC
  }  //for (iter= 0; iter < nbiter; iter++)
#endif 
  
  return (return_val);
}

 
int RDCost_for_macroblocks64 (double   lambda,       // <-- lagrange multiplier
                            int      mode,         // <-- modus (0-COPY/DIRECT, 1-16x16, 2-16x8, 3-8x16, 4-8x8(+), 5-Intra4x4, 6-Intra16x16)
                            double*  min_rdcost,   // <-> minimum rate-distortion cost
                            double*  min_rate,     // --> bitrate of mode which has minimum rate-distortion cost. 
                            int i16mode)//MB64X64
{
  
  int         i, j, k, mbi, mbj, ii, jj; 
  int         j1, j2;
  int         rate = 0, coeff_rate = 0;
  int64       distortion = 0;
  double      rdcost;
  Macroblock  *currMB;
  //Macroblock  *prevMB;
  //int         bframe    = (img->type==B_SLICE);
  //int         use_of_cc =  (img->type!=I_SLICE &&  input->symbol_mode!=CABAC);
  int         dummy;
  int      return_val = 0;
  //int tmp_cc;
  //int cc_rate;
  int saved_luma_transform_size_8x8_flag;
  int encodeMbAddr, CurrentMbAddr_x, CurrentMbAddr_y, extend_x, extend_y, num_MB, nextMbAddr, mb16_count=0, mb32_count=0;
  int ****i4p, ***i3p;
  //int mb32_cbp[4];//cbp of 4 32x32 blocks
  int **i2p;

#ifdef MB32X32_MVC
  int         cr_cbp = 0;
  short iter, nbiter;
  short ******all_mv = img->all_mv;
  short loopbreak=0;

  if ((input->mv_competition > 0) && (mode == 0) && (img->type==P_SLICE))
    nbiter = mv_comp.nb_mode_for_skip + 1;
  else
    nbiter = 1;
  
  // Let's iterate on each possible mode for the SKIP
  for (iter= 0; iter < nbiter; iter++)
  {
    loopbreak = 0;
    if ((input->mv_competition > 0) && (mode == 0) && (img->type==P_SLICE))  
    {
      distortion = 0;
      coeff_rate = 0;
      rate = 0;
      cr_cbp = 0;
    }
#endif

  /*******************************
   * restore some global variables
   *******************************/
  set_MB_parameters(saved_mb32_nr);
  saved_luma_transform_size_8x8_flag = img->mb_data[img->current_mb_nr].luma_transform_size_8x8_flag;


  //=====   S T O R E   C O D I N G   S T A T E   =====
  //---------------------------------------------------
  store_coding_state (cs_cm);
    
  
  find_Mb_cluster(2, img->current_mb_nr, &CurrentMbAddr_x, &CurrentMbAddr_y, &extend_x, &extend_y, &num_MB, &nextMbAddr);

  mb32_count=0;
  for(mbi=CurrentMbAddr_y; mbi<=extend_y; mbi+=2)  
  {      
    for(mbj=CurrentMbAddr_x; mbj<=extend_x; mbj+=2)   
    {             
      for(ii=0; ii<2; ii++)
      {
        for(jj=0; jj<2; jj++)
        {
          encodeMbAddr = (mbi+ii)*img->PicWidthInMbs + (mbj+jj);
          
          start_macroblock32 (encodeMbAddr, FALSE, 0);//reset img->current_mb_nr and other global variables

          currMB   = &img->mb_data[img->current_mb_nr];   
          //prevMB   = img->current_mb_nr ? &img->mb_data[img->current_mb_nr-1] : NULL;
          //=====
          //=====  SET REFERENCE FRAMES AND BLOCK MODES
          //=====    
          SetModesAndRefframeForBlocks32 (mode, best8x8pdir[mode][mb32_count], best8x8fwref[mode][mb32_count], best8x8bwref[mode][mb32_count]);
        }
      }
      best16x16pdir[mb32_count] = currMB->b8pdir[0];
      mb32_count++;
    }
  }

  mb16_count=0;
  mb32_count=0;
  for(mbi=CurrentMbAddr_y; mbi<=extend_y; mbi+=2)  
  {      
  for(mbj=CurrentMbAddr_x; mbj<=extend_x; mbj+=2)    
  {             
    for(ii=0; ii<2; ii++)    
    {
    for(jj=0; jj<2; jj++)
    {
      encodeMbAddr = (mbi+ii)*img->PicWidthInMbs + (mbj+jj);
      
      start_macroblock32 (encodeMbAddr, FALSE, 0);//reset img->current_mb_nr and other global variables

      currMB   = &img->mb_data[img->current_mb_nr];   
      //prevMB   = img->current_mb_nr ? &img->mb_data[img->current_mb_nr-1] : NULL;
      //=====
      //=====  SET REFERENCE FRAMES AND BLOCK MODES
      //=====    
      img->mb_data[encodeMbAddr].luma_transform_size_8x8_flag = saved_luma_transform_size_8x8_flag;
     
      SetModesAndRefframeForBlocks32 (mode, best8x8pdir[mode][mb32_count], best8x8fwref[mode][mb32_count], best8x8bwref[mode][mb32_count]);

    
      //=====
      //=====  GET COEFFICIENTS, RECONSTRUCTIONS, CBP
      //=====
#ifndef SIMPLIFY_CODE
      if (bframe && mode==0)
      {
        int block_x=img->pix_x>>2;
        int block_y=img->pix_y>>2;    
        for (j = block_y;j< block_y + 4;j++)
          for (i = block_x;i<block_x + 4;i++)
            if (direct_pdir[j][i] < 0) // not necessary in the spatial direct mode
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
#endif     

#ifdef MB32X32_MVC
      if ((input->mv_competition > 0) && (mode == 0) && (img->type==P_SLICE))
      {
        // The last iteration consist of re-launching with the bestmode,
        // to set needed variables for the JM
        if (iter == (nbiter - 1))
          mv_comp.current_predictor_for_skip = MB64.best_predictor_for_skip;
        else
        {
          mv_comp.current_predictor_for_skip = iter;
          if ((mv_comp.mv_pred_skip[mv_comp.current_predictor_for_skip][0] == NOT_AVAILABLE) || (mv_comp.mv_pred_skip[mv_comp.current_predictor_for_skip][1] == NOT_AVAILABLE))
          {
            loopbreak=1;
            goto nextinloop;
            //continue;
          }
        }

        prepare_mv(mbj+jj-CurrentMbAddr_x, mbi+ii-CurrentMbAddr_y);
        only_16x16_transform = (currMB->luma_transform_size_8x8_flag == 2);
        LumaResidualCoding ();
        only_16x16_transform = 0;
        switch_back_mv();

        prepare_mv(mbj+jj-CurrentMbAddr_x, mbi+ii-CurrentMbAddr_y);
        img->i16offset = 0;
        dummy = 0;
        if ((img->yuv_format!=YUV400) && (mode != IPCM))
          ChromaResidualCoding (&dummy);
        switch_back_mv();
      }

      else //if((input->mv_competition > 0) && (mode == 0) && (img->type==P_SLICE))
      {
#endif
        if (mode<P8x8)
        {
          prepare_mv(mbj+jj-CurrentMbAddr_x, mbi+ii-CurrentMbAddr_y);

          only_16x16_transform = currMB->luma_transform_size_8x8_flag == 2;

          LumaResidualCoding ();

          only_16x16_transform = 0;
          switch_back_mv();

        }
#ifndef SIMPLIFY_CODE
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
            if(Motion_Selected && IS_INTER_SLICE(img->type))
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
            if(Motion_Selected && IS_INTER_SLICE(img->type))
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
            if(Motion_Selected && IS_INTER_SLICE(img->type))
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
#endif

#ifdef MB32X32_MVC      
      }
#endif

#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
#if 1 //def TRELLIS_ONLY
        if(0)
#else
        if(input->UseRDO_Q || img->slice_fractional_quant_flag)
#endif
#else
        if(input->UseRDO_Q)
#endif
        {
          if (img->type == I_SLICE  || !IS_INTRA_MODE(mode) || (IS_INTRA_MODE(mode) && !Motion_Selected && IS_INTER_SLICE(img->type)))
          { 
#ifndef SIMPLIFY_CODE
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
#endif          



            img->i16offset = 0;
            dummy = 0;
            if ((!(img->residue_transform_flag && (mode==I4MB || mode==I16MB || mode==I8MB)))
              && (img->yuv_format!=YUV400) && (mode != IPCM))
              ChromaResidualCoding (&dummy);
            
#ifndef SIMPLIFY_CODE
            if (mode==I16MB)     
              img->i16offset = I16Offset  (currMB->cbp, i16mode);
            
            if(IS_INTER_SLICE(img->type) && IS_INTRA_MODE(mode) && !Motion_Selected)
              store_macroblock_parameters_intra(mode); 
#endif
          }
        }
        else
        {
#ifndef SIMPLIFY_CODE
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
#endif

          prepare_mv(mbj+jj-CurrentMbAddr_x, mbi+ii-CurrentMbAddr_y);

          img->i16offset = 0;
          dummy = 0;
          if ((!(img->residue_transform_flag && (mode==I4MB || mode==I16MB || mode==I8MB)))
            && (img->yuv_format!=YUV400) && (mode != IPCM))
            ChromaResidualCoding (&dummy);

          switch_back_mv();


#ifndef SIMPLIFY_CODE        
          if (mode==I16MB)     
            img->i16offset = I16Offset  (currMB->cbp, i16mode);
#endif
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
        

#ifndef SIMPLIFY_CODE
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
      if( (input->UseRDO_Q || img->slice_fractional_quant_flag) && IS_INTRA_MODE(mode) && Motion_Selected && IS_INTER_SLICE(img->type))
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
#endif      
        //=====
        //=====   GET DISTORTION
        //=====
        // LUMA
#ifndef SIMPLIFY_CODE
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
#endif
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
#ifndef SIMPLIFY_CODE
        //=====
        //=====   GET RATE
        //=====
        //----- macroblock header -----
        if (use_of_cc) //use_of_cc =  (img->type!=I_SLICE &&  input->symbol_mode!=CABAC);
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
#endif

        /******************************************
         *for the upper left 16x16 block, set b8mode
         *b8pdir for the whole 32x32 block in order
         *to write motion info in 1 shot
         *****************************************/
        //if(mb_count == 0)
        //{
        //  int i;
        //  for(i=0;i<4;i++)
        //  {
        //    currMB->b8mode[i] = mode;
        //    currMB->b8pdir[i] = best16x16pdir[i];
        //  }    
        //}
//      SetModesAndRefframeForBlocks32 (mode, best8x8pdir[mode][mb_count], best8x8fwref[mode][mb_count], best8x8bwref[mode][mb_count]);
        //rate += writeMBLayer32 (1, &coeff_rate, 1, mb_count);


        /******************************************
         *save macroblock info for each 16x16 block
         *****************************************/
        i4p=cofAC32_RDCost[mb16_count]; cofAC32_RDCost[mb16_count]=img->cofAC; img->cofAC=i4p;  
        i3p=cofDC32_RDCost[mb16_count]; cofDC32_RDCost[mb16_count]=img->cofDC; img->cofDC=i3p;  
        i2p=cof32AC16x16_RDCost[mb16_count]; cof32AC16x16_RDCost[mb16_count]=img->cofAC16x16; img->cofAC16x16=i2p;
        cbp16_RDCost[mb16_count] = currMB->cbp;  
        cbp_blk16_RDCost[mb16_count] = currMB->cbp_blk;  
  
        if (((currMB->cbp & 15) == 0) && !(IS_OLDINTRA(currMB) || currMB->mb_type == I8MB))//originally done in set_stored_macroblock_parameters
          luma_transform_size_8x8_flag32_RDCost[mb16_count] = 0;
        else
          luma_transform_size_8x8_flag32_RDCost[mb16_count] = currMB->luma_transform_size_8x8_flag;
        
        mb16_count++;
        }//jj
        }//ii


        //prepare 32x32 block cbp     
        cbp_32_RDCost[mb32_count]=0;
        for(k=0; k<4; k++)
        {
          if(cbp16_RDCost[mb32_count*4+k] != 0)
          {
            cbp_32_RDCost[mb32_count] = 1;
            break;
          }
        }      

        mb32_count++;

      
   } // mbj             
   }//mbi

#ifdef MB32X32_MVC
nextinloop:
  if(loopbreak==1)
    continue; //for (iter= 0; iter < nbiter; iter++)
#endif

  //prepare 64x64 block cbp     
  MB64.cbp=0;
  for(k=0; k<4; k++)
  {
    if(cbp_32_RDCost[k] != 0)
    {
      MB64.cbp = 1;//will be used in writeMBLayer64
      break;
    }
  }      

  MB64.mb_type = mode;//will be used in writeMBLayer64

  mb16_count=0;
  mb32_count=0;
  for(mbi=CurrentMbAddr_y; mbi<=extend_y; mbi+=2)  
  {      
    for(mbj=CurrentMbAddr_x; mbj<=extend_x; mbj+=2)   
    {             
      MB32.cbp = cbp_32_RDCost[mb32_count];//will be used in writeMBLayer64

      for(ii=0; ii<2; ii++)
      {
        for(jj=0; jj<2; jj++)
        {
          encodeMbAddr = (mbi+ii)*img->PicWidthInMbs + (mbj+jj);
          
          start_macroblock32 (encodeMbAddr, FALSE, 0);//reset img->current_mb_nr and other global variables

          currMB   = &img->mb_data[img->current_mb_nr];   
     
            /******************************************
             *for the upper left 16x16 block, set b8mode
             *b8pdir for the whole 32x32 block in order
             *to write motion info in 1 shot
             * mv and ref_idx are set in SetModesAndRefframeForBlocks32
             *****************************************/
            if(mb32_count == 0 && mb16_count == 0)
            {
              int i;
              for(i=0;i<4;i++)
              {
                currMB->b8mode[i] = mode;
                currMB->b8pdir[i] = best16x16pdir[i];
              }    
            }

            prepare_data_for_write_RDCost(encodeMbAddr, mb16_count, mode);
            img->mb_data[encodeMbAddr].luma_transform_size_8x8_flag = saved_luma_transform_size_8x8_flag;

            rate += writeMBLayer64 (1, &coeff_rate, mb32_count, mb16_count&3);
            switchback_data_after_write(encodeMbAddr, mb16_count, mode);

            mb16_count++;
        }
      }
       
      mb32_count++;    
    } // j       
  }//i
      
      
      //=====   R E S T O R E   C O D I N G   S T A T E   =====
      //-------------------------------------------------------
      //restore some global variables
      set_MB_parameters(saved_mb32_nr);
      //restore some global variables

      reset_coding_state (cs_cm);

      
      rdcost = (double)distortion + lambda * max(0.5,(double)rate);



#ifndef SIMPLIFY_CODE
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
      if( (input->UseRDO_Q || img->slice_fractional_quant_flag) && IS_INTRA_MODE(mode) && IS_INTER_SLICE(img->type))
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
#endif    
    if (!(rdcost >= *min_rdcost ||
      ((img->qp + img->bitdepth_luma_qp_scale)==0 && img->lossless_qpprime_flag==1 && distortion!=0)))
    {               
#ifndef SIMPLIFY_CODE
      if (((img->MbaffFrameFlag) && (mode ? 0: ((img->type == B_SLICE) ? !currMB->cbp:1)))  // AFF and current is skip
        && (img->current_mb_nr & 0x01) //bottom
        && (prevMB->mb_type ? 0:((img->type == B_SLICE) ? !prevMB->cbp:1)) //top is skip
        && (!(field_flag_inference() == currMB->mb_field)))
        return_val = 0;
      else      
#endif
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
          
#ifdef MB32X32_MVC
          if ((input->mv_competition > 0) && (mode == 0) && (img->type==P_SLICE))
          {
            int by, bx;
            
            MB64.best_predictor_for_skip = mv_comp.current_predictor_for_skip;
            
            for (by = 0;by < 16;by++)
            {
              for (bx = 0;bx < 16;bx++)
              {
                all_mv [by][bx][0][0][0][0] = mv_comp.mv_pred_skip[MB64.best_predictor_for_skip][0];
                all_mv [by][bx][0][0][0][1] = mv_comp.mv_pred_skip[MB64.best_predictor_for_skip][1];
              }
            }
          }
#endif
          return_val = 1;
      }
      
      
    }
    else 
    {
#ifdef MB32X32_MVC
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
#ifdef MB32X32_MVC
      }
#endif
    }
#ifdef MB32X32_MVC
  }  //for (iter= 0; iter < nbiter; iter++)
#endif  
  
  return (return_val);
}

/*!
*************************************************************************************
* \brief
*    Store 32x32 block parameters
*************************************************************************************
*/
void store_macroblock_parameters32 (int mode, int mb_ext_level)
{
  int  j, ****i4p, ***i3p;
  Macroblock *currMB  = &img->mb_data[img->current_mb_nr];
  int        bframe   = (img->type==B_SLICE);
  int **i2p;
  //--- store best mode ---
  best_mode = mode;

  best_c_imode = currMB->c_ipred_mode;
  best_i16offset = img->i16offset;


  // If condition is not really necessary.
  bi_pred_me = (mode == 1) ? currMB->bi_pred_me : 0;  

  b8mode[0] = b8mode[1] = b8mode[2] = b8mode[3] = mode;
  memcpy(b8pdir, best16x16pdir, BLOCK_MULTIPLE * sizeof(int));

  
#ifndef SIMPLIFY_CODE
  // Residue Color Transform
  //for (k = 0, j=img->block_y; j<img->block_y+BLOCK_MULTIPLE; j++, k+=BLOCK_MULTIPLE)
  memcpy(b4_intra_pred_modes,currMB->intra_pred_modes, MB_BLOCK_PARTITIONS * sizeof(char));
  memcpy(b8_intra_pred_modes8x8,currMB->intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));
  for (j = 0 ; j<BLOCK_MULTIPLE; j++)
  {
    memcpy(&b4_ipredmode[j * BLOCK_MULTIPLE],&img->ipredmode[img->block_y + j][img->block_x],BLOCK_MULTIPLE * sizeof(char));
    memcpy(b8_ipredmode8x8[j],&img->ipredmode8x8[img->block_y + j][img->block_x],BLOCK_MULTIPLE * sizeof(char));
  }  
#endif
  //--- reconstructed blocks ----
  for (j = 0; j < (MB_BLOCK_SIZE<<mb_ext_level); j++)

  {
    memcpy(rec_mbY[j],&enc_picture->imgY[img->pix_y+j][img->pix_x], (MB_BLOCK_SIZE<<mb_ext_level) * sizeof(imgpel));
#ifndef SIMPLIFY_CODE
    if((img->type==SP_SLICE) && (si_frame_indicator==0 && sp2_frame_indicator==0))
      memcpy(lrec_rec[j],&lrec[img->pix_y+j][img->pix_x], MB_BLOCK_SIZE * sizeof(int));//store coefficients SP frame
#endif
  }
  
  if (img->AdaptiveRounding)
    store_adaptive_rounding_parameters (mode, currMB);
  
  if (img->yuv_format != YUV400)
  {
    for (j = 0; j<(img->mb_cr_size_y<<mb_ext_level); j++)
    {
      memcpy(rec_mbU[j],&enc_picture->imgUV[0][img->pix_c_y+j][img->pix_c_x], (img->mb_cr_size_x<<mb_ext_level) * sizeof(imgpel));
      memcpy(rec_mbV[j],&enc_picture->imgUV[1][img->pix_c_y+j][img->pix_c_x], (img->mb_cr_size_x<<mb_ext_level) * sizeof(imgpel));
#ifndef SIMPLIFY_CODE
      if((img->type==SP_SLICE) && (si_frame_indicator==0 && sp2_frame_indicator==0))
      {//store uv coefficients SP frame
        memcpy(lrec_rec_U[j],&lrec_uv[0][img->pix_c_y+j][img->pix_c_x], img->mb_cr_size_x * sizeof(int));
        memcpy(lrec_rec_V[j],&lrec_uv[1][img->pix_c_y+j][img->pix_c_x], img->mb_cr_size_x * sizeof(int));
      }
#endif
    }
  }
  
#ifndef SIMPLIFY_CODE
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
#endif
  //--- coeff, cbp, kac ---
  if (mode || bframe)
  {
    int i;
    for(i=0; i<(1<<(2*mb_ext_level)); i++)
    {
      i4p=cofAC32[i]; cofAC32[i]=cofAC32_RDCost[i]; cofAC32_RDCost[i]=i4p;
      i3p=cofDC32[i]; cofDC32[i]=cofDC32_RDCost[i]; cofDC32_RDCost[i]=i3p;    
      i2p=cof32AC16x16[i]; cof32AC16x16[i]=cof32AC16x16_RDCost[i]; cof32AC16x16_RDCost[i]=i2p;
    }


    if(input->UseExtMB == 2)//MB64X64
      memcpy(cbp_32, cbp_32_RDCost, (1<<(2*(mb_ext_level-1))) * sizeof(int)); 

    memcpy(cbp16, cbp16_RDCost, (1<<(2*mb_ext_level)) * sizeof(int)); 
    memcpy(cbp_blk16, cbp_blk16_RDCost, (1<<(2*mb_ext_level)) * sizeof(int64));
  }
  else
  {
    if(input->UseExtMB == 2)//MB64X64
      memset(cbp_32, 0, (1<<(2*(mb_ext_level-1))) * sizeof(int)); 

    memset(cbp16, 0, (1<<(2*mb_ext_level)) * sizeof(int)); 
    memset(cbp_blk16, 0, (1<<(2*mb_ext_level)) * sizeof(int64));
  }
  
  memcpy(luma_transform_size_8x8_flag32, luma_transform_size_8x8_flag32_RDCost, (1<<(2*mb_ext_level)) * sizeof(int)); 


  for (j = 0; j<(4<<mb_ext_level); j++)
    memcpy(frefframe32[j],&enc_picture->ref_idx[LIST_0][img->block_y+j][img->block_x], (BLOCK_MULTIPLE<<mb_ext_level) * sizeof(char));
  if (bframe)
  {
    for (j = 0; j<(4<<mb_ext_level); j++)
      memcpy(brefframe32[j],&enc_picture->ref_idx[LIST_1][img->block_y+j][img->block_x], (BLOCK_MULTIPLE<<mb_ext_level) * sizeof(char));
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


}

void compute_mode_RD_cost32(int mode, 
                          Macroblock *currMB, 
                          RD_PARAMS enc_mb, 
                          double *min_rdcost, 
                          double *min_rate, 
                          int i16mode, 
                          short bslice, 
                          short *inter_skip,
                          int   mb_ext_level)
{

  //--- transform size ---           
  currMB->luma_transform_size_8x8_flag = input->Transform8x8Mode==2
    ?  (mode >= 1 && mode <= 3)
    || (mode == 0 && bslice && active_sps->direct_8x8_inference_flag)
    || ((mode == P8x8) && (enc_mb.valid[4]))
    :  0; 

  //currMB->luma_transform_size_8x8_flag = 2;

  //store_coding_state (cs_cm); // RD  
  //SetModesAndRefframeForBlocks32 (mode);  //this func is called inside RDCost_for_macroblock32
  
  // Encode with coefficients
  img->NoResidueDirect = 0;  
#ifndef SIMPLIFY_CODE
  if (currMB->c_ipred_mode == DC_PRED_8 || (IS_INTRA(currMB) ))
#endif
  {
    while(1)
    {
#ifdef MV_COMPETITION             
#ifdef MB32X32_MVC
      if ((mode == 0) && (img->type==P_SLICE) && (input->mv_competition > 0))
        skip_mode = TRUE;
      else 
#endif      
        skip_mode = FALSE;
#endif      
      
      //printf("\n-----------------------RDCost_for_macroblock32--------------------\n\n");
      if (RDCost_for_macroblocks32 (enc_mb.lambda_md, mode, min_rdcost, min_rate, i16mode, mb_ext_level))
      {//printf("\n store macroblock\n\n\n\n");
#ifndef SIMPLIFY_CODE
        //Rate control
        if (input->RCEnable)
        {
          if(mode == P8x8)
            rc_store_diff(img->opix_x,img->opix_y,
            currMB->luma_transform_size_8x8_flag == 1 ? tr8x8.mpr8x8 : tr4x4.mpr8x8);
          else
            rc_store_diff(img->opix_x, img->opix_y, pred);
        }
#endif

        store_macroblock_parameters32 (mode, mb_ext_level);

#ifndef SIMPLIFY_CODE
        if(input->rdopt==2 && mode == 0 && input->EarlySkipEnable)
        {
          // check transform quantized coeff.
          if(currMB->cbp == 0)
            *inter_skip = 1;
        }
#endif        
      }
      
      // Go through transform modes. 
      // Note that if currMB->cbp is 0 one could choose to skip 8x8 mode
      // although this could be due to deadzoning decisions.
      //if (input->Transform8x8Mode==1 && currMB->cbp!=0) 
      if (input->Transform8x8Mode==1)
      {
        //=========== try mb_types 1,2,3 with 8x8 transform ===========
        if ((mode >= 1 && mode <= 3) && currMB->luma_transform_size_8x8_flag == 0)
        {
          //try with 8x8 transform size
          currMB->luma_transform_size_8x8_flag = 1;
          continue;
        }
        // try partition 16x16 with 16x16 transform
        else if ((mode == 1) && currMB->luma_transform_size_8x8_flag == 1)
        {
          //try with 16x16 transform size
          currMB->luma_transform_size_8x8_flag = 2;
          continue;
        }
        // try partition 16x16 with 16x16 transform
        else if ((mode == 2) && currMB->luma_transform_size_8x8_flag == 1)
        {
          //try with 16x16 transform size
          currMB->luma_transform_size_8x8_flag = 2;
          continue;
        }
        // try partition 16x16 with 16x16 transform
        else if ((mode == 3) && currMB->luma_transform_size_8x8_flag == 1)
        {
          //try with 16x16 transform size
          currMB->luma_transform_size_8x8_flag = 2;
          continue;
        }
        //=========== try DIRECT-MODE with 8x8 transform ===========
        else if (mode == 0 && bslice && active_sps->direct_8x8_inference_flag && currMB->luma_transform_size_8x8_flag == 0)
        {
          //try with 8x8 transform size
          currMB->luma_transform_size_8x8_flag = 1;
          continue;
        }
        else if (mode == 0 && bslice && active_sps->direct_8x8_inference_flag && currMB->luma_transform_size_8x8_flag == 1)
        {
          //try with 16x16 transform size
          currMB->luma_transform_size_8x8_flag = 2;
          continue;
        }
#ifndef SIMPLIFY_CODE
        //=========== try mb_type P8x8 for mode 4 with 4x4/8x8 transform ===========
        else if ((mode == P8x8) && (enc_mb.valid[4]) && (currMB->luma_transform_size_8x8_flag == 0))
        {
          currMB->luma_transform_size_8x8_flag = 1; //check 8x8 partition for transform size 8x8
          continue;
        }
#endif
        else
        {
          currMB->luma_transform_size_8x8_flag = 0;
          break;
        }
      }
      else
        break;
    }
    
    // Encode with no coefficients. Currently only for direct. This could be extended to all other modes as in example.
    //if (mode < P8x8 && (*inter_skip == 0) && enc_mb.valid[mode] && currMB->cbp && (currMB->cbp&15) != 15 && !input->nobskip)
    if ( bslice && mode == 0 && (*inter_skip == 0) && enc_mb.valid[mode] 
      && currMB->cbp && !input->nobskip) 
    {
      img->NoResidueDirect = 1;
      if (RDCost_for_macroblocks32 (enc_mb.lambda_md, mode, min_rdcost, min_rate, i16mode, mb_ext_level)) 
      {
#ifndef SIMPLIFY_CODE
        //Rate control
        if (input->RCEnable)
          rc_store_diff(img->opix_x,img->opix_y,pred);
#endif        

        store_macroblock_parameters32 (mode, mb_ext_level);
      }
    }
  }
};

/*!
*************************************************************************************
* \brief
*    computation of prediction list (including biprediction) cost
*************************************************************************************
*/
void list_prediction_cost32(int list, int block, int mode, RD_PARAMS enc_mb, int bmcost[5], char best_ref[2], int mb_ext_level)
{
  short ref;
  int mcost;
  int cur_list = list < BI_PRED ? enc_mb.list_offset[list] : enc_mb.list_offset[LIST_0];
  
  //--- get cost and reference frame for forward prediction ---
  
  if (list < BI_PRED)
  {
    for (ref=0; ref < listXsize[cur_list]; ref++)
    {
      if (!img->checkref || list || ref==0 || CheckReliabilityOfRef (block, list, ref, mode)) //img->checkref=0 
      {
        // limit the number of reference frames to 1 when switching SP frames are used
        if((!input->sp2_frame_indicator && !input->sp_output_indicator)||
          ((input->sp2_frame_indicator || input->sp_output_indicator) && (img->type!=P_SLICE && img->type!=SP_SLICE))||
          ((input->sp2_frame_indicator || input->sp_output_indicator) && ((img->type==P_SLICE || img->type==SP_SLICE) &&(ref==0))))
        {
          mcost  = (input->rdopt 
            ? REF_COST (enc_mb.lambda_mf, ref, cur_list) 
            : (int) (2 * enc_mb.lambda_me * min(ref, 1)));     
          
          mcost += motion_cost[mode][list][ref][block];
          if (mcost < bmcost[list])
          {
            bmcost[list]   = mcost;
            best_ref[list] = (char)ref;
          }
        }
      }
    }
  }
  else if (list == BI_PRED)
  {
    bmcost[list]  = (input->rdopt 
      ? (REF_COST  (enc_mb.lambda_mf, (short)best_ref[LIST_0], cur_list)
      +  REF_COST  (enc_mb.lambda_mf, (short)best_ref[LIST_1], cur_list + LIST_1)) 
      : (int) (2 * (enc_mb.lambda_me * (min((short)best_ref[LIST_0], 1) + min((short)best_ref[LIST_1], 1)))));    
    bmcost[list] += BIDPartitionCost32 (mode, block, (short)best_ref[LIST_0], (short)best_ref[LIST_1], enc_mb.lambda_mf, mb_ext_level);    
  }
  else
  {
    bmcost[list]  = (input->rdopt 
      ? (REF_COST (enc_mb.lambda_mf, 0, cur_list) 
      +  REF_COST (enc_mb.lambda_mf, 0, cur_list + LIST_1)) 
      : (int) (4 * enc_mb.lambda_me));
    bmcost[list] += BPredPartitionCost32(mode, block, 0, 0, enc_mb.lambda_mf, !(list&1), mb_ext_level);                
  }
}  

void assign_enc_picture_params32(int mode, int best_pdir, int block, int list_offset, int best_fw_ref, int best_bw_ref, int bframe, int mb_ext_level)
{
  int i,j;
  int block_x, block_y;
  short *cur_mv;
  
  if (mode==1)
  {
    assert(block==0);
    if (best_pdir==1)
    {
      for (j=img->block_y+(block&2); j<img->block_y+(block&2) + (BLOCK_MULTIPLE<<mb_ext_level); j++)   //block=0 b/c mode =1
      {
        block_x = img->block_x+(block&1)*2; //block=1 b/c mode =1
        
        memset(&enc_picture->ref_idx[LIST_0][j][block_x], -1 ,     (BLOCK_MULTIPLE<<mb_ext_level) * sizeof(char));
        memset(enc_picture->mv      [LIST_0][j][block_x],  0 , 2 * (BLOCK_MULTIPLE<<mb_ext_level) * sizeof(short));
        for (i=block_x; i<block_x + (BLOCK_MULTIPLE<<mb_ext_level); i++)
        {
          enc_picture->ref_pic_id [LIST_0][j][i]    = -1;
        }
      }
    }
    else if (img->bi_pred_me[mode])
    {
      for (j=0; j<(BLOCK_MULTIPLE<<mb_ext_level); j++)
      {
        block_y = img->block_y+(block&2)+j;  //block=0 b/c mode =1
        block_x = img->block_x+(block&1)*2;   //block=0 b/c mode =1
        memset(&enc_picture->ref_idx[LIST_0][block_y][block_x], 0, (BLOCK_MULTIPLE<<mb_ext_level) * sizeof(char));
        for (i=0; i<(BLOCK_MULTIPLE<<mb_ext_level); i++)
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
      for (j=0; j<(BLOCK_MULTIPLE<<mb_ext_level); j++)
      {
        block_y = img->block_y+(block&2)+j;  //block=0 b/c mode =1
        block_x = img->block_x+(block&1)*2;  //block=0 b/c mode =1
        memset(&enc_picture->ref_idx[LIST_0][block_y][block_x], best_fw_ref , (BLOCK_MULTIPLE<<mb_ext_level) * sizeof(char));
        for (i=0; i<(BLOCK_MULTIPLE<<mb_ext_level); i++)
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
        for (j=img->block_y+(block&2); j<img->block_y+(block&2) + (BLOCK_MULTIPLE<<mb_ext_level); j++)  //block=0 b/c mode =1 
        {
          block_x = img->block_x+(block&1)*2;  //block=0 b/c mode =1
          memset(&enc_picture->ref_idx[LIST_1][j][block_x], -1 , (BLOCK_MULTIPLE<<mb_ext_level) * sizeof(char));
          memset(enc_picture->mv[LIST_1][j][block_x], 0 , 2 * (BLOCK_MULTIPLE<<mb_ext_level) * sizeof(short));
          for (i=block_x; i<block_x + (BLOCK_MULTIPLE<<mb_ext_level); i++)
          {
            enc_picture->ref_pic_id [LIST_1][j][i] = -1;
          }
        }
      }
      else
      {
        if (img->bi_pred_me[mode])
        {
          for (j=0; j<(BLOCK_MULTIPLE<<mb_ext_level); j++)
          {
            block_y = img->block_y+(block&2)+j;  //block=0 b/c mode =1
            block_x = img->block_x+(block&1)*2;  //block=0 b/c mode =1
            memset(&enc_picture->ref_idx[LIST_1][block_y][block_x], 0, (BLOCK_MULTIPLE<<mb_ext_level) * sizeof(char));
            for (i=0; i<(BLOCK_MULTIPLE<<mb_ext_level); i++)
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
          for (j=0; j<(BLOCK_MULTIPLE<<mb_ext_level); j++)
          {
            block_y = img->block_y+(block&2)+j; //block=0 b/c mode =1
            block_x = img->block_x+(block&1)*2; // //block=0 b/c mode =1
            memset(&enc_picture->ref_idx[LIST_1][block_y][block_x], best_bw_ref, (BLOCK_MULTIPLE<<mb_ext_level) * sizeof(char));
            for (i=0; i<(BLOCK_MULTIPLE<<mb_ext_level); i++)
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
    for (j=0; j<(2<<mb_ext_level); j++)
    {
      block_y = img->block_y + block * (2<<mb_ext_level) + j;
      for (i=0; i<(BLOCK_MULTIPLE<<mb_ext_level); i++)
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
          cur_mv = img->all_mv[j+block*(2<<mb_ext_level)][i][LIST_0][best_fw_ref][mode];
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
              cur_mv = img->all_mv[j+ block*(2<<mb_ext_level)][i][LIST_1][best_bw_ref][mode];
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
    for (j=0; j<(BLOCK_MULTIPLE<<mb_ext_level); j++)
    {
      block_y = img->block_y+j;
      for (i=0; i<(2<<mb_ext_level); i++)
      {
        block_x = img->block_x + block*(2<<mb_ext_level) + i;
        if (best_pdir==1)
        {
          enc_picture->ref_idx    [LIST_0][block_y][block_x]    = -1;
          enc_picture->ref_pic_id [LIST_0][block_y][block_x]    = -1;
          enc_picture->mv         [LIST_0][block_y][block_x][0] = 0;
          enc_picture->mv         [LIST_0][block_y][block_x][1] = 0;
        }
        else
        {
          cur_mv = img->all_mv[j][block*(2<<mb_ext_level)+i][LIST_0][best_fw_ref][mode];
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
              cur_mv = img->all_mv[j][block*(2<<mb_ext_level)+i][LIST_1][best_bw_ref][mode];
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

/*!
*************************************************************************************
* \brief
*    Set reference frames and motion vectors
*************************************************************************************
*/
void SetRefAndMotionVectors32 (int block, int mode, int pdir, int fwref, int bwref, int mb_ext_level)
{
  int     i, j=0;
  int     bslice  = (img->type==B_SLICE);
  int     pmode   = (mode==1||mode==2||mode==3?mode:4);
  int     j0      = ((block >> 1)<<1);  //block == 0 
  int     i0      = ((block & 0x01)<<1);//block == 0 
  int     j1      = j0 + ((input->part_size[pmode][1])<<mb_ext_level);
  int     i1      = i0 + ((input->part_size[pmode][0])<<mb_ext_level);
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
      memset(&enc_picture->ref_idx[LIST_0][j][img->block_x + i0], -1, ((input->part_size[pmode][0]) << mb_ext_level) * sizeof(char));
      memset(&enc_picture->ref_idx[LIST_1][j][img->block_x + i0], -1, ((input->part_size[pmode][0]) <<mb_ext_level) * sizeof(char));
      memset(enc_picture->mv[LIST_0][j][img->block_x + i0], 0, 2*((input->part_size[pmode][0]) <<mb_ext_level) * sizeof(short));
      memset(enc_picture->mv[LIST_1][j][img->block_x + i0], 0, 2*((input->part_size[pmode][0]) <<mb_ext_level) * sizeof(short));
    }
    return;
  }
  
  if (!bslice)
  {
    for (j=j0; j<j1; j++)
    {
      block_y = img->block_y + j;
      memset(&enc_picture->ref_idx   [LIST_0][block_y][img->block_x + i0], fwref, ((input->part_size[pmode][0]) <<mb_ext_level) * sizeof(char));
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
*    Sets motion vectors for a macroblock
*************************************************************************************
*/
void SetMotionVectorsMB32 (Macroblock* currMB, int bframe, int mb_ext_level)
{
  int i, j, k, mode8, pdir8, ref, by, bx;
  short ******all_mv  = img->all_mv;
  //short ******pred_mv = img->pred_mv;
  int  bw_ref;
  int jdiv;

  if (!bframe)
  {
    for (j = 0; j<(4<<mb_ext_level); j++)
    {
#ifndef SIMPLIFY_CODE
      jmod = j & 0x01;
#endif

      jdiv = j >>   (1+mb_ext_level);
      by    = img->block_y+j;

      for (i = 0; i<(4<<mb_ext_level); i++)
      {
        mode8 = currMB->b8mode[k=2*jdiv+(i>>(1+mb_ext_level))];

#ifndef SIMPLIFY_CODE
        l     = 2*jmod + (i & 0x01);
#endif        
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
    for (j = 0; j<(4<<mb_ext_level); j++)
    {
#ifndef SIMPLIFY_CODE
      jmod = j & 0x01;
#endif
      jdiv = j >>   (1+mb_ext_level);
      by    = img->block_y+j;

      for (i = 0; i<(4<<mb_ext_level); i++)
      {
        mode8 = currMB->b8mode[k=2*jdiv+(i>>(1+mb_ext_level))];
#ifndef SIMPLIFY_CODE
        l     = 2*jmod + (i & 0x01);
#endif
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
 
#ifndef SIMPLIFY_CODE
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
#endif
}

/*!
*************************************************************************************
* \brief
*    Set stored 32x32 block parameters
*************************************************************************************
*/
void set_stored_macroblock_parameters32 (int mb_ext_level)
{
  int  i, j, k;
  Macroblock  *currMB  = &img->mb_data[img->current_mb_nr];
  int         mode     = best_mode; //best_mode is decided in store_macroblock_parameters32 
  int         bframe   = (img->type==B_SLICE);
  //char    **ipredmodes = img->ipredmode;
  
  imgpel        **imgY  = enc_picture->imgY;
  imgpel       ***imgUV = enc_picture->imgUV;
  int        block_x, block_y;  
  short   *cur_mv;

#ifndef SIMPLIFY_CODE
#ifdef ADAPTIVE_FD_SD_CODING
  currMB->SD_Coding_on_off=best_SD_Coding_on_off;
  memcpy(currMB->quantizer_indices,best_quantizer_indices,256*sizeof(int));
  memcpy(currMB->SD_or_FD,best_SD_or_FD,4*sizeof(int));
  currMB->SD_or_FD_t8x8=best_SD_or_FD_t8x8;
#endif
#ifdef ADAPTIVE_QUANTIZATION
  currMB->mb_iaqms_idx = best_mb_iaqms_idx;
#endif
#endif
  //===== reconstruction values =====
  for (j = 0; j < (MB_BLOCK_SIZE<<mb_ext_level); j++) 
  {
    memcpy(&imgY[img->pix_y+j][img->pix_x],rec_mbY[j], (MB_BLOCK_SIZE<<mb_ext_level) * sizeof(imgpel)); 
#ifndef SIMPLIFY_CODE
    if((img->type==SP_SLICE) &&(si_frame_indicator==0 && sp2_frame_indicator==0 ))
      memcpy(&lrec[img->pix_y+j][img->pix_x],lrec_rec[j], (MB_BLOCK_SIZE<<mb_ext_level) * sizeof(int)); //restore coeff SP frame  
    if(img->MbaffFrameFlag)
      memcpy(rdopt->rec_mbY[j],rec_mbY[j], (MB_BLOCK_SIZE<<mb_ext_level) * sizeof(imgpel)); //
#endif
  }
#ifndef SIMPLIFY_CODE
#ifdef ADAPTIVE_FD_SD_CODING
  if(img->MbaffFrameFlag)
  {
    rdopt->best_SD_Coding_on_off=best_SD_Coding_on_off;
    memcpy(rdopt->best_quantizer_indices,best_quantizer_indices,256*sizeof(int));
    memcpy(rdopt->best_SD_or_FD,best_SD_or_FD,4*sizeof(int));
    rdopt->best_SD_or_FD_t8x8=best_SD_or_FD_t8x8;
  }
#endif
#endif
  if (img->AdaptiveRounding)
  {
    update_offset_params(mode,luma_transform_size_8x8_flag);
  }
  
  if (img->yuv_format != YUV400)
  {
    for (j = 0; j<(img->mb_cr_size_y << mb_ext_level); j++)
    {
      memcpy(&imgUV[0][img->pix_c_y+j][img->pix_c_x],rec_mbU[j], (img->mb_cr_size_x <<mb_ext_level) * sizeof(imgpel));
      memcpy(&imgUV[1][img->pix_c_y+j][img->pix_c_x],rec_mbV[j], (img->mb_cr_size_x <<mb_ext_level) * sizeof(imgpel));
#ifndef SIMPLIFY_CODE
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
#endif
    }
  }
  
  //==== macroblock type ====
  currMB->mb_type = mode;
  if(mb_ext_level == 1)
  {
    MB32.mb_type    = mode;

    MB32.cbp=0;
    //prepare 32x32 block cbp      
    for(i=0; i<4; i++)
    {
      if(cbp16[i] != 0)
      {
        MB32.cbp = 1;
        break;
      }
    }      
  }
 
  if(mb_ext_level == 2)//MB64X64
  {
    int mbi;
    MB64.cbp=0;
    //prepare 32x32 block cbp      
    for(mbi=0; mbi<4; mbi++)
    {
      if(cbp_32[mbi] != 0)
      {
        MB64.cbp = 1;
        break;
      }
    }      
  }


#ifndef SIMPLIFY_CODE
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
#endif
  
  memcpy(currMB->b8mode,b8mode, BLOCK_MULTIPLE * sizeof(int));
  memcpy(currMB->b8pdir,b8pdir, BLOCK_MULTIPLE * sizeof(int));
  if(mb_ext_level == 1)
  {
    memcpy(MB32.b8mode,   b8mode, BLOCK_MULTIPLE * sizeof(int));
    memcpy(MB32.b8pdir,   b8pdir, BLOCK_MULTIPLE * sizeof(int));
  }
  if(mb_ext_level == 2)//MB64X64
  {
    memcpy(MB64.b8mode,   b8mode, BLOCK_MULTIPLE * sizeof(int));
    memcpy(MB64.b8pdir,   b8pdir, BLOCK_MULTIPLE * sizeof(int));
  }

#ifndef SIMPLIFY_CODE
  if(img->MbaffFrameFlag)
  {
    memcpy(rdopt->b8mode,b8mode, BLOCK_MULTIPLE * sizeof(int));
    memcpy(rdopt->b8pdir,b8pdir, BLOCK_MULTIPLE * sizeof(int));    
  }
#endif

  currMB->bi_pred_me = currMB->mb_type == 1 ? bi_pred_me : 0;  
  if(mb_ext_level == 1)
    MB32.bi_pred_me    = currMB->bi_pred_me;

  if(mb_ext_level == 2)//MB64X64
    MB64.bi_pred_me    = currMB->bi_pred_me;

#ifndef SIMPLIFY_CODE
  //if P8x8 mode and transform size 4x4 choosen, restore motion vector data for this transform size 
  if (mode == P8x8 && !luma_transform_size_8x8_flag && input->Transform8x8Mode)
    RestoreMV8x8(1);
#endif

  rdopt->luma_transform_size_8x8_flag  = currMB->luma_transform_size_8x8_flag;


#ifndef SIMPLIFY_CODE
  if (input->rdopt==3 && img->type!=B_SLICE)
  {
    //! save the MB Mode of every macroblock
    decs->dec_mb_mode[img->mb_x][img->mb_y] = mode;
  }
#endif
  //==== reference frames =====
  for (j = 0; j < (4<<mb_ext_level); j++)
  {
    block_y = img->block_y + j;
    for (i = 0; i < (4<<mb_ext_level); i++)
    {
      block_x = img->block_x + i;
      k = 2*(j >> (1+mb_ext_level))+(i >> (1+mb_ext_level));
      
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
          cur_mv = img->all_mv[j][i][LIST_0][(short)frefframe32[j][i]][currMB->b8mode[k]];
          
          enc_picture->ref_idx    [LIST_0][block_y][block_x]    = frefframe32[j][i];
          enc_picture->ref_pic_id [LIST_0][block_y][block_x]    = enc_picture->ref_pic_num[LIST_0 + currMB->list_offset][(short)frefframe32[j][i]];
          enc_picture->mv         [LIST_0][block_y][block_x][0] = cur_mv[0];
          enc_picture->mv         [LIST_0][block_y][block_x][1] = cur_mv[1];
          if(img->MbaffFrameFlag)
            rdopt->refar[LIST_0][j][i] = frefframe32[j][i];
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
    for (j=0; j<(4<<mb_ext_level); j++)
    {
      block_y = img->block_y + j;
      for (i=0; i<(4<<mb_ext_level); i++)
      {          
        block_x = img->block_x + i;
        k = 2*(j >> (1+mb_ext_level))+(i >> (1+mb_ext_level));
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
            cur_mv = img->all_mv[j][i][LIST_1][(short)brefframe32[j][i]][currMB->b8mode[k]];
            
            enc_picture->ref_idx    [LIST_1][block_y][block_x] = brefframe32[j][i];
            enc_picture->ref_pic_id [LIST_1][block_y][block_x] = enc_picture->ref_pic_num[LIST_1 + currMB->list_offset][(short)brefframe32[j][i]];
            enc_picture->mv         [LIST_1][block_y][block_x][0] = cur_mv[0];
            enc_picture->mv         [LIST_1][block_y][block_x][1] = cur_mv[1];
            if(img->MbaffFrameFlag)
              rdopt->refar[LIST_1][j][i] = brefframe32[j][i];
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
    for(j = img->block_y; j < img->block_y + (BLOCK_MULTIPLE<<mb_ext_level); j++)
      memset(&img->ipredmode[j][img->block_x], DC_PRED, (BLOCK_MULTIPLE<<mb_ext_level) * sizeof(char));
  }
  // Residue Color Transform
  else if (mode == I4MB)
  {
    memcpy(currMB->intra_pred_modes,b4_intra_pred_modes, MB_BLOCK_PARTITIONS * sizeof(char));
    for(j = 0; j < BLOCK_MULTIPLE; j++)
      memcpy(&img->ipredmode[img->block_y + j][img->block_x],&b4_ipredmode[BLOCK_MULTIPLE * j], BLOCK_MULTIPLE * sizeof(char));
  }
#ifndef SIMPLIFY_CODE    
  if(img->MbaffFrameFlag)
  {
    rdopt->c_ipred_mode = currMB->c_ipred_mode;
    rdopt->i16offset = img->i16offset;  
    memcpy(rdopt->intra_pred_modes,currMB->intra_pred_modes, MB_BLOCK_PARTITIONS * sizeof(char));
    memcpy(rdopt->intra_pred_modes8x8,currMB->intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));
    for(j = img->block_y; j < img->block_y +BLOCK_MULTIPLE; j++)
      memcpy(&rdopt->ipredmode[j][img->block_x],&ipredmodes[j][img->block_x], BLOCK_MULTIPLE * sizeof(char));
  }
#endif
  //==== motion vectors =====
  SetMotionVectorsMB32 (currMB, bframe, mb_ext_level);
  
#ifdef USE_INTRA_MDDT
  if(input->UseIntraMDDT)
  {
    memcpy(img->cofAC16[0], bestCofAC16[0], sizeof(int)*16*17);
    memcpy(img->cofAC16[1], bestCofAC16[1], sizeof(int)*16*17);
  }
#endif 
}





void encode_one_macroblock32 (int mb_ext_level)
{  
  int max_index;
  int         rerun, block, index=0, mode, i, j, ctr16x16;
  short       best_pdir;
  RD_PARAMS   enc_mb;
  double      min_rdcost, max_rdcost=1e30;
  char        best_ref[2] = {0, -1};
  int         bmcost[5] = {INT_MAX};
  int         cost=0;
  int         min_cost = INT_MAX, i16mode=0;
  int         intra1 = 0;
  //int         temp_cpb = 0;
  //int         best_transform_flag = 0;
  //int         cost8x8_direct = 0;  
  short       islice      = (img->type==I_SLICE);
  short       bslice      = (img->type==B_SLICE);
  short       pslice      = (img->type==P_SLICE) || (img->type==SP_SLICE);
  short       intra       = (islice || (pslice && img->mb_y==img->mb_y_upd && img->mb_y_upd!=img->mb_y_intra));
  
  short       runs        = (input->RestrictRef==1 && input->rdopt==3 && (pslice  || (bslice && img->nal_reference_idc>0)) ? 2 : 1);

  //int         pix_x, pix_y;
  Macroblock* currMB      = &img->mb_data[img->current_mb_nr];
  //int         prev_mb_nr  = FmoGetPreviousMBNr(img->current_mb_nr);
  //Macroblock* prevMB      = (prev_mb_nr >= 0) ? &img->mb_data[prev_mb_nr]:NULL ;
  
  //char   **ipredmodes = img->ipredmode;
  //short   *allmvs = img->all_mv[0][0][0][0][0];
#ifndef SIMPLIFY_CODE
  short   max_chroma_pred_mode;
#endif
  //int     ****i4p;  //for non-RD-opt. mode
  
  //int tmp_8x8_flag, tmp_no_mbpart;  
  // Residue Color Transform
  //int residue_R, residue_G, residue_B, temp;
  //int cr_cbp = 0;  
  // Fast Mode Decision
  short inter_skip = 0;
  //int cost16 = 0, mode16 = 0;
  double min_rate = 0;//, RDCost16 = DBL_MAX;
#ifdef ADAPTIVE_FD_SD_CODING
  memset(currMB->SD_or_FD,0,4*sizeof(int));
  currMB->SD_or_FD_t8x8=0;
#endif
  
#ifdef RDO_Q
  if(input->UseRDO_Q && active_pps->entropy_coding_mode_flag == CABAC)
  {
    estRunLevel_CABAC(LUMA_4x4);
    estRunLevel_CABAC(LUMA_8x8);
    estRunLevel_CABAC(LUMA_16AC);
#ifdef USE_INTRA_MDDT
    estRunLevel_CABAC(LUMA_16x16);
#endif
    estRunLevel_CABAC(LUMA_16x16P);
    estRunLevel_CABAC(LUMA_16x8P);
    estRunLevel_CABAC(LUMA_8x16P);
  }
#endif

#ifndef SIMPLIFY_CODE
  if(input->FMEnable == 1)
  {
    decide_intrabk_SAD();
  }
  else if (input->FMEnable ==2)
  {
    simplified_decide_intrabk_SAD();
  }
  intra |= RandomIntra (img->current_mb_nr);    // Forced Pseudo-Random Intra
#endif
 
  //===== Setup Macroblock encoding parameters =====
  init_enc_mb_params(currMB, &enc_mb, intra, bslice);

  //enc_mb.valid[0]     = 0;
  //enc_mb.valid[1]     = 0;
  //enc_mb.valid[2]     = 0;
  //enc_mb.valid[3]     = 0;



  // Perform multiple encodings if rdopt with losses is enabled
  for (rerun=0; rerun<runs; rerun++)
  {
#ifndef SIMPLIFY_CODE
    if (runs==2)
      input->rdopt= (rerun==0) ? 1 : 3;
#endif
    // reset chroma intra predictor to default
    currMB->c_ipred_mode = DC_PRED_8;
    
    //=====   S T O R E   C O D I N G   S T A T E   =====
    //---------------------------------------------------
    store_coding_state (cs_cm);
    
    if (!intra)
    {
      //===== set direct motion vectors =====
      best_mode = 1;
      if (bslice)
      {
        Get_Direct_Motion_Vectors32 (mb_ext_level);

#ifndef SIMPLIFY_CODE
        if (input->rdopt == 2 && enc_mb.valid[0])
        {
          best_mode = 0;
          currMB->c_ipred_mode=DC_PRED_8;
          min_rdcost = max_rdcost;
          compute_mode_RD_cost(0, currMB, enc_mb, &min_rdcost, &min_rate, i16mode, bslice, &inter_skip);
        }
#endif
      }
      
      //===== MOTION ESTIMATION FOR 32x32, 32x16, 16x32 BLOCKS =====      
      for (min_cost=INT_MAX, mode=1; mode<4; mode++)
      {
        bi_pred_me = 0;
        img->bi_pred_me[mode]=0;
        if (enc_mb.valid[mode] && !inter_skip)
        {
          for (cost=0, block=0; block<(mode==1?1:2); block++)
          {
#ifdef RDO_Q
            PartitionMotionSearch32 (mode, block, enc_mb.lambda_mf, 0, mb_ext_level);
#else
            PartitionMotionSearch32 (mode, block, enc_mb.lambda_mf, mb_ext_level);
#endif
            
            //--- set 4x4 block indizes (for getting MV) ---
            j = (block==1 && mode==2 ? (2<<mb_ext_level) : 0);
            i = (block==1 && mode==3 ? (2<<mb_ext_level) : 0);
            //--- get cost and reference frame for List 0 prediction ---
            bmcost[LIST_0] = INT_MAX;
            list_prediction_cost32(LIST_0, block, mode, enc_mb, bmcost, best_ref, mb_ext_level);
            if (bslice)
            {
              //--- get cost and reference frame for List 1 prediction ---
              bmcost[LIST_1] = INT_MAX;
              list_prediction_cost32(LIST_1, block, mode, enc_mb, bmcost, best_ref, mb_ext_level);
              // Compute bipredictive cost between best list 0 and best list 1 references
              list_prediction_cost32(BI_PRED, block, mode, enc_mb, bmcost, best_ref, mb_ext_level);
              // Finally, if mode 16x16, compute cost for bipredictive ME vectore
              if (input->BiPredMotionEstimation && mode == 1)
              {                
                list_prediction_cost32(BI_PRED_L0, block, mode, enc_mb, bmcost, 0, mb_ext_level);
                list_prediction_cost32(BI_PRED_L1, block, mode, enc_mb, bmcost, 0, mb_ext_level);
              }
              else
              {
                bmcost[BI_PRED_L0] = INT_MAX;
                bmcost[BI_PRED_L1] = INT_MAX;
              }
              
              // Determine prediction list based on mode cost
              determine_prediction_list(mode, bmcost, best_ref, &best_pdir, &cost, &bi_pred_me);
            }
            else // if (bslice)
            {
              best_pdir  = 0;
              cost      += bmcost[LIST_0];
            }
            assign_enc_picture_params32(mode, best_pdir, block, enc_mb.list_offset[LIST_0], best_ref[LIST_0], best_ref[LIST_1], bslice, mb_ext_level);
            //----- set reference frame and direction parameters -----
            if (mode==3)
            {
              best8x8fwref [3][block  ] = best8x8fwref [3][  block+2] = best_ref[LIST_0];
              best8x8pdir  [3][block  ] = best8x8pdir  [3][  block+2] = best_pdir;
              best8x8bwref [3][block  ] = best8x8bwref [3][  block+2] = best_ref[LIST_1];
            }
            else if (mode==2)
            {
              best8x8fwref [2][2*block] = best8x8fwref [2][2*block+1] = best_ref[LIST_0];
              best8x8pdir  [2][2*block] = best8x8pdir  [2][2*block+1] = best_pdir;
              best8x8bwref [2][2*block] = best8x8bwref [2][2*block+1] = best_ref[LIST_1];
            }
            else
            {
              best8x8fwref [1][0] = best8x8fwref [1][1] = best8x8fwref [1][2] = best8x8fwref [1][3] = best_ref[LIST_0];
              best8x8pdir  [1][0] = best8x8pdir  [1][1] = best8x8pdir  [1][2] = best8x8pdir  [1][3] = best_pdir;
              best8x8bwref [1][0] = best8x8bwref [1][1] = best8x8bwref [1][2] = best8x8bwref [1][3] = best_ref[LIST_1];
            }
            
            //--- set reference frames and motion vectors ---
            //generate context for block 1 
            if (mode>1 && block==0)
              SetRefAndMotionVectors32 (block, mode, best_pdir, best_ref[LIST_0], best_ref[LIST_1], mb_ext_level);            
          } // for (block=0; block<(mode==1?1:2); block++)
#ifndef SIMPLIFY_CODE          
          if(!input->rdopt)
          {
            currMB->luma_transform_size_8x8_flag = 0;
            if (input->Transform8x8Mode) //for inter rd-off, set 8x8 to do 8x8 transform
            {
              SetModesAndRefframeForBlocks32(mode);
              currMB->luma_transform_size_8x8_flag = TransformDecision(-1, &cost);
            }
          }          
          
          if(input->rdopt == 2 && mode == 1)
          {
            if(pslice)
              min_rdcost = max_rdcost;
            
            //=====   S T O R E   C O D I N G   S T A T E   =====
            //---------------------------------------------------
            //store_coding_state (cs_cm);
            
            for (ctr16x16=0, k=0; k<1; k++)
            {
              i16mode = 0; 
              
              //--- for INTER16x16 check all prediction directions ---
              if (bslice)
              {
                best8x8pdir[1][0] = best8x8pdir[1][1] = best8x8pdir[1][2] = best8x8pdir[1][3] = ctr16x16;
                if ( (bslice) && (input->BiPredMotionEstimation) 
                  && (ctr16x16 == 2 && img->bi_pred_me[mode] < 2 && mode == 1))
                  ctr16x16--;
                if (ctr16x16 < 2) 
                  index--;
                ctr16x16++;
              }
              
              currMB->c_ipred_mode=DC_PRED_8;
              compute_mode_RD_cost(mode, currMB, enc_mb, &min_rdcost, &min_rate, i16mode, bslice, &inter_skip);
              
              if ((input->BiPredMotionEstimation) && (bslice) && ctr16x16 == 2 
                && img->bi_pred_me[mode] < 2 && mode == 1 && best8x8pdir[1][0] == 2) 
                img->bi_pred_me[mode] = img->bi_pred_me[mode] + 1;
            } // for (ctr16x16=0, k=0; k<1; k++)
            
            if(pslice)
            {
              // Get SKIP motion vector and compare SKIP_MV with best motion vector of 16x16
              FindSkipModeMotionVector ();
              
              if(input->EarlySkipEnable)
              {
                //===== check for SKIP mode =====
                if ( currMB->cbp==0 && enc_picture->ref_idx[LIST_0][img->block_y][img->block_x]==0 &&
                  enc_picture->mv[LIST_0][img->block_y][img->block_x][0]==allmvs[0] &&
                  enc_picture->mv[LIST_0][img->block_y][img->block_x][1]==allmvs[1]               )
                {
                  inter_skip = 1;
                  best_mode = 0;
                }
              } // if(input->EarlySkipEnable)
            }
            
            // store variables.
            RDCost16 = min_rdcost;
            mode16 = best_mode;
            cost16 = cost;
          } // if(input->rdopt == 2 && mode == 1)
#endif
          if ((!inter_skip) && (cost < min_cost))
          {
            best_mode = mode;
            min_cost  = cost;
#ifndef SIMPLIFY_CODE
            best_transform_flag = currMB->luma_transform_size_8x8_flag;
#endif
          }
        } // if (enc_mb.valid[mode])
      } // for (mode=1; mode<4; mode++)
#ifndef SIMPLIFY_CODE
      if ((!inter_skip) && enc_mb.valid[P8x8])
      {
        giRDOpt_B8OnlyFlag = 1;
        
        tr8x8.cost8x8 = INT_MAX;
        tr4x4.cost8x8 = INT_MAX;
        //===== store coding state of macroblock =====
        store_coding_state (cs_mb);
        
        currMB->all_blk_8x8 = -1;
        
        if (input->Transform8x8Mode)
        {  
          tr8x8.cost8x8 = 0;
          //===========================================================
          // Check 8x8 partition with transform size 8x8 
          //===========================================================
          //=====  LOOP OVER 8x8 SUB-PARTITIONS  (Motion Estimation & Mode Decision) =====
          for (cost_direct=cbp8x8=cbp_blk8x8=cnt_nonz_8x8=0, block=0; block<4; block++)
          {
            submacroblock_mode_decision(enc_mb, &tr8x8, currMB, cofAC_8x8ts[block],
              &have_direct, bslice, block, &cost_direct, &cost, &cost8x8_direct, 1);
            best8x8mode       [block] = tr8x8.part8x8mode [block];
            best8x8pdir [P8x8][block] = tr8x8.part8x8pdir [block];
            best8x8fwref[P8x8][block] = tr8x8.part8x8fwref[block];
            best8x8bwref[P8x8][block] = tr8x8.part8x8bwref[block];
          }
          
          // following params could be added in RD_8x8DATA structure
          cbp8_8x8ts      = cbp8x8;
          cbp_blk8_8x8ts  = cbp_blk8x8;
          cnt_nonz8_8x8ts = cnt_nonz_8x8;
          currMB->luma_transform_size_8x8_flag = 0; //switch to 4x4 transform size
          
          //--- re-set coding state (as it was before 8x8 block coding) ---
          //reset_coding_state (cs_mb);        
        }// if (input->Transform8x8Mode)
        
        
        if (input->Transform8x8Mode != 2)  
        {
          tr4x4.cost8x8 = 0;
          //=================================================================
          // Check 8x8, 8x4, 4x8 and 4x4 partitions with transform size 4x4
          //=================================================================
          //=====  LOOP OVER 8x8 SUB-PARTITIONS  (Motion Estimation & Mode Decision) =====
          for (cost_direct=cbp8x8=cbp_blk8x8=cnt_nonz_8x8=0, block=0; block<4; block++)
          {
            submacroblock_mode_decision(enc_mb, &tr4x4, currMB, cofAC8x8[block],
              &have_direct, bslice, block, &cost_direct, &cost, &cost8x8_direct, 0);
            
            best8x8mode       [block] = tr4x4.part8x8mode [block];
            best8x8pdir [P8x8][block] = tr4x4.part8x8pdir [block];
            best8x8fwref[P8x8][block] = tr4x4.part8x8fwref[block];
            best8x8bwref[P8x8][block] = tr4x4.part8x8bwref[block];
          }          
          //--- re-set coding state (as it was before 8x8 block coding) ---
          // reset_coding_state (cs_mb);  
        }// if (input->Transform8x8Mode != 2)
        
        //--- re-set coding state (as it was before 8x8 block coding) ---
        reset_coding_state (cs_mb);
        
        
        // This is not enabled yet since mpr has reverse order.
        if (input->RCEnable)
          rc_store_diff(img->opix_x,img->opix_y,img->mpr);
        
        //check cost for P8x8 for non-rdopt mode
        if (!input->rdopt && (tr4x4.cost8x8 < min_cost || tr8x8.cost8x8 < min_cost))
        {
          best_mode = P8x8;
          if (input->Transform8x8Mode == 2)
          {
            min_cost = tr8x8.cost8x8;
            currMB->luma_transform_size_8x8_flag=1;
          }
          else if (input->Transform8x8Mode)
          {
            if (tr8x8.cost8x8 < tr4x4.cost8x8)
            {
              min_cost = tr8x8.cost8x8;
              currMB->luma_transform_size_8x8_flag=1;
            }
            else if(tr4x4.cost8x8 < tr8x8.cost8x8)
            {
              min_cost = tr4x4.cost8x8;
              currMB->luma_transform_size_8x8_flag=0;
            }
            else
            {
              if (GetBestTransformP8x8() == 0)
              {
                min_cost = tr4x4.cost8x8;
                currMB->luma_transform_size_8x8_flag=0;
              }
              else
              {
                min_cost = tr8x8.cost8x8;
                currMB->luma_transform_size_8x8_flag=1;
              }
            }
          }
          else
          {
            min_cost = tr4x4.cost8x8;
            currMB->luma_transform_size_8x8_flag=0;
          }
        }// if (!input->rdopt && (tr4x4.cost8x8 < min_cost || tr8x8.cost8x8 < min_cost))
        giRDOpt_B8OnlyFlag = 0;
      }
      else // if (enc_mb.valid[P8x8])
#endif
      {
        tr4x4.cost8x8 = INT_MAX;
      }
      
      // Find a motion vector for the Skip mode
      if(input->rdopt != 2 && pslice)
        FindSkipModeMotionVector32 (mb_ext_level);
    }   
    else // if (!intra)
    {
      min_cost = INT_MAX;
    }
    
    
    //========= C H O O S E   B E S T   M A C R O B L O C K   M O D E =========
    //-------------------------------------------------------------------------
    if (input->rdopt)
    {
      MB32.prev_qp = currMB->prev_qp;
      MB32.prev_delta_qp = currMB->prev_delta_qp; //for delta_qp encoding in cabac
      MB32.qp = currMB->qp;
      saved_mb32_nr = img->current_mb_nr;

      // store_coding_state (cs_cm);    
      if (!inter_skip)
      {
        //int mb_available_up;
        //int mb_available_left;
        //int mb_available_up_left;
#ifndef SIMPLIFY_CODE   
        if(input->rdopt == 2 && img->type!=I_SLICE)
        {
          min_rdcost = RDCost16;
          best_mode  = mode16;
        }
        else
#endif
          min_rdcost = max_rdcost;
        
        // if Fast High mode, compute  inter modes separate process for inter/intra
        max_index =  4;
        if (input->BiPredMotionEstimation)
          img->bi_pred_me[1] =0;  
        
#ifndef SIMPLIFY_CODE
        if (img->yuv_format != YUV400 && max_index != 5)
        {
          // precompute all new chroma intra prediction modes
          IntraChromaPrediction(&mb_available_up, &mb_available_left, &mb_available_up_left);
          max_chroma_pred_mode = PLANE_8;
        }
        else
          max_chroma_pred_mode = DC_PRED_8;


        for (currMB->c_ipred_mode=DC_PRED_8; currMB->c_ipred_mode<=max_chroma_pred_mode; currMB->c_ipred_mode++)
        {
          // bypass if c_ipred_mode is not allowed
          if ( (img->yuv_format != YUV400) &&
            (  ((!intra || !input->IntraDisableInterOnly) && input->ChromaIntraDisable == 1 && currMB->c_ipred_mode!=DC_PRED_8) 
            || (currMB->c_ipred_mode == VERT_PRED_8 && !mb_available_up) 
            || (currMB->c_ipred_mode == HOR_PRED_8  && !mb_available_left) 
            || (currMB->c_ipred_mode == PLANE_8     && (!mb_available_left || !mb_available_up || !mb_available_up_left))))
            continue;        
#endif
          //===== GET BEST MACROBLOCK MODE =====
          for (ctr16x16=0, index=0; index < max_index; index++)
          {
            mode = mb_mode_table[index];
            
            if (img->yuv_format != YUV400)
            {           
#ifndef SIMPLIFY_CODE
              if (input->rdopt == 2)
              {
                i16mode = 0;              
                // RDcost of mode 1 in P-slice and mode 0, 1 in B-slice are already available
                if(((bslice && mode == 0) || (!islice && mode == 1)))
                  continue;
              }
              else
#endif
              {
#ifndef SIMPLIFY_CODE
                // Residue Color Transform
                if(img->residue_transform_flag)
                {
                  mode = mb_mode_table_RCT[index];
                  if( mode == I16MB) 
                    i16mode = index -5;
                  // bypass if i16mode is not allowed
                  if (mode == I16MB &&
                    (  (i16mode==VERT_PRED_16 && !mb_available_up) 
                    || (i16mode==HOR_PRED_16  && !mb_available_left) 
                    || (i16mode==PLANE_16    && (!mb_available_left || !mb_available_up || !mb_available_up_left))))
                    continue;
                }
                else
#endif
                {
                  mode = mb_mode_table[index];
                  i16mode = 0; 
                }
              }
            }
            //--- for INTER16x16 check all prediction directions ---
            if (mode==1 && bslice)
            {
              best8x8pdir[1][0] = best8x8pdir[1][1] = best8x8pdir[1][2] = best8x8pdir[1][3] = ctr16x16;
              if ( (bslice) && (input->BiPredMotionEstimation) 
                && (ctr16x16 == 2 && img->bi_pred_me[mode] < 2 && mode == 1))
                ctr16x16--;
              if (ctr16x16 < 2) 
                index--;
              ctr16x16++;
            }
#ifndef SIMPLIFY_CODE
            // Skip intra modes in inter slices if best inter mode is 
            // a MB partition and cbp is 0.
            if (input->SkipIntraInInterSlices && !intra && mode >= I16MB 
              && best_mode <=3 && currMB->cbp == 0)
              continue;
#endif
            if (enc_mb.valid[mode])
              compute_mode_RD_cost32(mode, currMB, enc_mb, &min_rdcost, &min_rate, i16mode, bslice, &inter_skip, mb_ext_level);
            if ((input->BiPredMotionEstimation) && (bslice) && ctr16x16 == 2 
              && img->bi_pred_me[mode] < 2 && mode == 1 && best8x8pdir[1][0] == 2) 
              img->bi_pred_me[mode] = img->bi_pred_me[mode] + 1;

          }// for (ctr16x16=0, index=0; index<max_index; index++)
#ifndef SIMPLIFY_CODE
        }// for (currMB->c_ipred_mode=DC_PRED_8; currMB->c_ipred_mode<=max_chroma_pred_mode; currMB->c_ipred_mode++)

        // Selective Intra Coding
        if(img->type!=I_SLICE && input->rdopt == 2 && input->SelectiveIntraEnable && input->ProfileIDC<FREXT_HP)
        {
          fast_mode_intra_decision(&intra_skip, min_rate);
          
          if(!intra_skip)
          {
            // precompute all new chroma intra prediction modes
            if (img->yuv_format != YUV400)
            {
              // precompute all new chroma intra prediction modes
              IntraChromaPrediction(&mb_available_up, &mb_available_left, &mb_available_up_left);
              max_chroma_pred_mode = PLANE_8;
            }
            else
              max_chroma_pred_mode = DC_PRED_8;
            
            max_index = 9;
            
            for (currMB->c_ipred_mode=DC_PRED_8; currMB->c_ipred_mode<=max_chroma_pred_mode; currMB->c_ipred_mode++)
            {
              
              // bypass if c_ipred_mode is not allowed
              if ( (img->yuv_format != YUV400) &&
                (  ((!intra || !input->IntraDisableInterOnly) && input->ChromaIntraDisable == 1 && currMB->c_ipred_mode!=DC_PRED_8) 
                || (currMB->c_ipred_mode == VERT_PRED_8 && !mb_available_up) 
                || (currMB->c_ipred_mode == HOR_PRED_8  && !mb_available_left) 
                || (currMB->c_ipred_mode == PLANE_8     && (!mb_available_left || !mb_available_up || !mb_available_up_left))))
                continue;           
              
              //===== GET BEST MACROBLOCK MODE =====
              for (index = 5; index < max_index; index++)
              {
                mode = mb_mode_table[index];
                
                if (input->SkipIntraInInterSlices && !intra && mode >= I16MB 
                  && best_mode <=3 && currMB->cbp == 0)
                  continue;
                
                if (img->yuv_format != YUV400)
                {           
                  if (input->rdopt == 2)
                  {
                    i16mode = 0;              
                    // RDcost of mode 1 in P-slice and mode 0, 1 in B-slice are already available
                    if(((bslice && mode == 0) || (!islice && mode == 1)))
                      continue;
                  }
                  else
                  {
                    // Residue Color Transform
                    if(img->residue_transform_flag)
                    {
                      mode = mb_mode_table_RCT[index];
                      if( mode == I16MB) 
                        i16mode = index -5;
                      // bypass if i16mode is not allowed
                      if (mode == I16MB &&
                        (  (i16mode==VERT_PRED_16 && !mb_available_up) 
                        || (i16mode==HOR_PRED_16  && !mb_available_left) 
                        || (i16mode==PLANE_16    && (!mb_available_left || !mb_available_up || !mb_available_up_left))))
                        continue;
                    }
                    else
                    {
                      mode = mb_mode_table[index];
                      i16mode = 0; 
                    }
                  }
                }
                
                if (enc_mb.valid[mode])
                  compute_mode_RD_cost(mode, currMB, enc_mb, &min_rdcost, &min_rate, i16mode, bslice, &inter_skip);                            
              } // for (index = 5; index < max_index; index++)
            }
          }
        }          
#endif
      }
#ifdef BEST_NZ_COEFF
      for (j=0;j<4;j++)
        for (i=0; i<(4+img->num_blk8x8_uv); i++)
          img->nz_coeff[img->current_mb_nr][j][i] = gaaiMBAFF_NZCoeff[j][i]; 
#endif
   }
#ifndef SIMPLIFY_CODE
   else //rdopt off
   {
     tmp_8x8_flag = currMB->luma_transform_size_8x8_flag;  //save 8x8_flag
     tmp_no_mbpart = currMB->NoMbPartLessThan8x8Flag;      //save no-part-less
     
     if (img->yuv_format != YUV400)
       // precompute all chroma intra prediction modes
       IntraChromaPrediction(NULL, NULL, NULL);
     
     if (enc_mb.valid[0] && bslice) // check DIRECT MODE
     {
       if(have_direct)
       {
         switch(input->Transform8x8Mode)
         {
         case 1: // Mixture of 8x8 & 4x4 transform
           cost = ((cost8x8_direct < cost_direct) || !(enc_mb.valid[5] && enc_mb.valid[6] && enc_mb.valid[7])) 
             ? cost8x8_direct : cost_direct;
           break;
         case 2: // 8x8 Transform only
           cost = cost8x8_direct;
           break;
         default: // 4x4 Transform only
           cost = cost_direct;
           break;
         }
       }
       else
       { //!have_direct
         cost = Get_Direct_CostMB (enc_mb.lambda_mf);
       }
       if (cost!=INT_MAX)
       {
         cost -= (int)floor(16*enc_mb.lambda_me+0.4999);
       }
       
       if (cost <= min_cost)
       {
         if(active_sps->direct_8x8_inference_flag && input->Transform8x8Mode)
         {
           if(input->Transform8x8Mode==2)
             currMB->luma_transform_size_8x8_flag=1;
           else
           {
             if(cost8x8_direct < cost_direct)
               currMB->luma_transform_size_8x8_flag=1;
             else
               currMB->luma_transform_size_8x8_flag=0;
           }
         }
         else
           currMB->luma_transform_size_8x8_flag=0;
         
         //Rate control
         if (input->RCEnable)
           rc_store_diff(img->opix_x,img->opix_y,img->mpr);
         
         min_cost  = cost;
         best_mode = 0;
       }
       else
       {
         currMB->luma_transform_size_8x8_flag = tmp_8x8_flag; // restore if not best
         currMB->NoMbPartLessThan8x8Flag = tmp_no_mbpart; // restore if not best
       }        
     }
     
     if (enc_mb.valid[I8MB]) // check INTRA8x8
     {
       currMB->luma_transform_size_8x8_flag = 1; // at this point cost will ALWAYS be less than min_cost 
       
       currMB->mb_type = I8MB;
       temp_cpb = Mode_Decision_for_new_Intra8x8Macroblock (enc_mb.lambda_md, &cost);
       
       if (cost <= min_cost)
       {
         // Residue Color Transform
         if(img->residue_transform_flag)
         {            
           for(i=0; i<2; i++) 
           {
             for(j=0; j<4; j++) 
               for(k=0; k<4; k++)
                 if(cbp_chroma_block[i][j][k]) cr_cbp = 2;
           }            
           cr_cbp = dct_chroma_DC(0, cr_cbp);
           cr_cbp = dct_chroma_DC(1, cr_cbp);
           
           temp_cpb += (cr_cbp<<4);
           for(j=0; j<MB_BLOCK_SIZE; j++) 
           {
             pix_y = img->pix_y + j;
             for(i=0; i<MB_BLOCK_SIZE; i++) 
             {
               pix_x = img->pix_x + i;
               temp_imgU[j][i] = enc_picture->imgUV[0][pix_y][pix_x];
               temp_imgV[j][i] = enc_picture->imgUV[1][pix_y][pix_x];
             }
           }
         }          
         currMB->cbp = temp_cpb;
         
         //coeffs
         if (input->Transform8x8Mode != 2)
         {
           i4p=cofAC; cofAC=img->cofAC; img->cofAC=i4p;
         }
         
         for(j=0; j<MB_BLOCK_SIZE; j++) 
         {
           pix_y = img->pix_y + j;
           for(i=0; i<MB_BLOCK_SIZE; i++) 
           {
             pix_x = img->pix_x + i;
             temp_imgY[j][i] = enc_picture->imgY[pix_y][pix_x];
           }
         }
         
         //Rate control
         if (input->RCEnable)
           rc_store_diff(img->opix_x,img->opix_y,img->mpr);
         
         min_cost  = cost;
         best_mode = I8MB;
         tmp_8x8_flag = currMB->luma_transform_size_8x8_flag;
       } 
       else
         currMB->luma_transform_size_8x8_flag = tmp_8x8_flag; // restore if not best
     }
     
     if (enc_mb.valid[I4MB]) // check INTRA4x4
     {
       currMB->luma_transform_size_8x8_flag = 0;
       currMB->mb_type = I4MB;
       temp_cpb = Mode_Decision_for_Intra4x4Macroblock (enc_mb.lambda_md, &cost);
       
       if (cost <= min_cost)
       {
         // Residue Color Transform
         if(img->residue_transform_flag)
         {
           for(i=0; i<2; i++) 
           { 
             for(j=0; j<4; j++) 
               for(k=0; k<4; k++) 
                 if(cbp_chroma_block[i][j][k]) 
                   cr_cbp = 2;
           }
           cr_cbp = dct_chroma_DC(0, cr_cbp);
           cr_cbp = dct_chroma_DC(1, cr_cbp);
           
           temp_cpb += (cr_cbp<<4);
         }
         currMB->cbp = temp_cpb;
         
         //Rate control
         if (input->RCEnable)
           rc_store_diff(img->opix_x,img->opix_y,img->mpr);
         
         min_cost  = cost;
         best_mode = I4MB;
         tmp_8x8_flag = currMB->luma_transform_size_8x8_flag;
       } 
       else
       {
         currMB->luma_transform_size_8x8_flag = tmp_8x8_flag; // restore if not best
         //coeffs
         i4p=cofAC; cofAC=img->cofAC; img->cofAC=i4p;
       }
     }
     if (enc_mb.valid[I16MB]) // check INTRA16x16
     {
       currMB->luma_transform_size_8x8_flag = 0;
       intrapred_luma_16x16 ();
       cost = find_sad_16x16 (&i16mode);
       
       if (cost < min_cost)
       {
         //Rate control
         // should this access opix or pix?
         if (input->RCEnable)
           rc_store_diff(img->opix_x,img->opix_y,img->mprr_2[i16mode]);
         
         // Residue Color Transform
         if(img->residue_transform_flag)
         {
           for (j = 0; j < MB_BLOCK_SIZE; j++) 
           {
             pix_y = img->pix_y + j;
             for (i = 0; i < MB_BLOCK_SIZE; i++) 
             {
               pix_x = img->pix_x + i;
               residue_G = imgY_org    [pix_y][pix_x] - img->mprr_2   [i16mode]             [j][i];
               residue_B = imgUV_org[0][pix_y][pix_x] - img->mprr_c[0][currMB->c_ipred_mode][j][i];
               residue_R = imgUV_org[1][pix_y][pix_x] - img->mprr_c[1][currMB->c_ipred_mode][j][i];                
               /* Forward Residue Transform */
               resTrans_R[j][i] = residue_R - residue_B;
               temp             = residue_B + (resTrans_R[j][i] >> 1);
               resTrans_B[j][i] = residue_G-temp;
               resTrans_G[j][i] = temp+(resTrans_B[j][i] >> 1);                
               img->m7[j][i]    = resTrans_G[j][i];
             }
           }
         }
         
         best_mode   = I16MB;
         currMB->cbp = dct_luma_16x16 (i16mode);
         
         // Residue Color Transform
         if(img->residue_transform_flag)
         {
           for (j = 0; j < MB_BLOCK_SIZE; j++) 
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
           
           for (j=0; j < MB_BLOCK_SIZE; j++) 
           {
             for (i=0; i < MB_BLOCK_SIZE; i++) 
               rec_resR[j][i] = img->m7[j][i];
           } 
           currMB->cbp += (cr_cbp<<4);
           
           /* Inverse Residue Transform */
           for (j=0; j < MB_BLOCK_SIZE; j++) 
           {
             pix_y = img->pix_y + j;
             for (i=0; i < MB_BLOCK_SIZE; i++) 
             {
               pix_x = img->pix_x + i;
               temp      = rec_resG[j][i] - (rec_resB[j][i] >> 1);
               residue_G = rec_resB[j][i]+temp;
               residue_B = temp - (rec_resR[j][i]>>1);
               residue_R = residue_B+rec_resR[j][i];                
               enc_picture->imgY    [pix_y][pix_x] = 
                 min(img->max_imgpel_value   ,max(0, residue_G + (int) img->mprr_2[i16mode][j][i]));
               enc_picture->imgUV[0][pix_y][pix_x] = 
                 min(img->max_imgpel_value_uv,max(0, residue_B + (int) img->mprr_c[0][currMB->c_ipred_mode][j][i]));
               enc_picture->imgUV[1][pix_y][pix_x] = 
                 min(img->max_imgpel_value_uv,max(0, residue_R + (int) img->mprr_c[1][currMB->c_ipred_mode][j][i]));
             }
           }
         }
       }
       else
       {
         currMB->luma_transform_size_8x8_flag = tmp_8x8_flag; // restore
         currMB->NoMbPartLessThan8x8Flag = tmp_no_mbpart;     // restore
       }
     }       
   }
#endif
     if (rerun==0)
       intra1 = IS_INTRA(currMB);
  } // for (rerun=0; rerun<runs; rerun++) 
  
  //=====  S E T   F I N A L   M A C R O B L O C K   P A R A M E T E R S ======
  //---------------------------------------------------------------------------  
  if (input->rdopt)
  {   
    if (((cbp!=0 || best_mode==I16MB) && (best_mode!=IPCM) )) //cbp is set in store_macroblock_parameters 
      currMB->prev_cbp = 1;    
    else if ((cbp==0 && !input->RCEnable) || (best_mode==IPCM))
    {
      currMB->delta_qp = 0;
      currMB->qp       = currMB->prev_qp;
      set_chroma_qp(currMB);
      img->qp          = currMB->qp;
      currMB->prev_cbp = 0;
    }    

/*if(img->type == B_SLICE)
{
  printf("good\n");
}*/

    set_stored_macroblock_parameters32 (mb_ext_level);
  }
#ifndef SIMPLIFY_CODE
  else
  {
    //===== set parameters for chosen mode =====
    SetModesAndRefframeForBlocks32 (best_mode);
    
    if (best_mode==P8x8)
    {
      if (currMB->luma_transform_size_8x8_flag && (cbp8_8x8ts == 0) && input->Transform8x8Mode != 2)
        currMB->luma_transform_size_8x8_flag = 0;
      
      SetCoeffAndReconstruction8x8 (currMB);
      
      memset(currMB->intra_pred_modes, DC_PRED, MB_BLOCK_PARTITIONS * sizeof(char));
      for (k=0, j = img->block_y; j < img->block_y + BLOCK_MULTIPLE; j++)
        memset(&ipredmodes[j][img->block_x], DC_PRED, BLOCK_MULTIPLE * sizeof(char));
    }
    else
    {
      //===== set parameters for chosen mode =====
      if (best_mode == I8MB)
      {
        memcpy(currMB->intra_pred_modes,currMB->intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));
        for(j = img->block_y; j < img->block_y + BLOCK_MULTIPLE; j++)
          memcpy(&img->ipredmode[j][img->block_x],&img->ipredmode8x8[j][img->block_x], BLOCK_MULTIPLE * sizeof(char));
        
        //--- restore reconstruction for 8x8 transform ---
        for(j=0; j<MB_BLOCK_SIZE; j++) 
        {
          memcpy(&enc_picture->imgY[img->pix_y + j][img->pix_x],temp_imgY[j], MB_BLOCK_SIZE * sizeof(imgpel));
        }
        
        // Residue Color Transform
        if(img->residue_transform_flag)
        {          
          for(j=0; j<MB_BLOCK_SIZE; j++) 
          {
            pix_y = img->pix_c_y + j;
            for(i=0; i<MB_BLOCK_SIZE; i++) 
            {
              pix_x = img->pix_c_x + i;
              enc_picture->imgUV[0][pix_y][pix_x] = temp_imgU[j][i] ;
              enc_picture->imgUV[1][pix_y][pix_x] = temp_imgV[j][i] ;
            }
          }                
        }
      }
      
      if ((best_mode!=I4MB)&&(best_mode != I8MB))
      {
        memset(currMB->intra_pred_modes,DC_PRED, MB_BLOCK_PARTITIONS * sizeof(char));
        for(j = img->block_y; j < img->block_y + BLOCK_MULTIPLE; j++)
          memset(&ipredmodes[j][img->block_x],DC_PRED, BLOCK_MULTIPLE * sizeof(char));
        
        if (best_mode!=I16MB)
        {
          if((best_mode>=1) && (best_mode<=3))
            currMB->luma_transform_size_8x8_flag = best_transform_flag;
          LumaResidualCoding ();
          
          if((currMB->cbp==0)&&(best_mode==0))
            currMB->luma_transform_size_8x8_flag = 0;
          
          //Rate control
          if (input->RCEnable)
            rc_store_diff(img->opix_x,img->opix_y,img->mpr);
        }
      }
    }
    //check luma cbp for transform size flag
    if (((currMB->cbp&15) == 0) && !(IS_OLDINTRA(currMB) || currMB->mb_type == I8MB))
      currMB->luma_transform_size_8x8_flag = 0;
    
    // precompute all chroma intra prediction modes
    if (img->yuv_format != YUV400)
      IntraChromaPrediction(NULL, NULL, NULL);
    
    img->i16offset = 0;
    dummy = 0;
    
    // Residue Color Transform
    if ((!(img->residue_transform_flag && (best_mode==I4MB || best_mode==I16MB || best_mode==I8MB))) && img->yuv_format!=YUV400)
      ChromaResidualCoding (&dummy);
    
    if (best_mode==I16MB) 
    {
      img->i16offset = I16Offset  (currMB->cbp, i16mode);
    }
    
    SetMotionVectorsMB (currMB, bslice);
    
    //===== check for SKIP mode =====
    if ((pslice) && best_mode==1 && currMB->cbp==0 &&
      enc_picture->ref_idx[LIST_0][img->block_y][img->block_x]    == 0 &&
      enc_picture->mv     [LIST_0][img->block_y][img->block_x][0] == allmvs[0] &&
      enc_picture->mv     [LIST_0][img->block_y][img->block_x][1] == allmvs[1])
    {
      currMB->mb_type = currMB->b8mode[0] = currMB->b8mode[1] = currMB->b8mode[2] = currMB->b8mode[3] = 0;
      currMB->luma_transform_size_8x8_flag = 0;
    }
    
    if(img->MbaffFrameFlag)
      set_mbaff_parameters();
  }

  // Rate control
  if(input->RCEnable)
    update_rc(currMB, best_mode);
#endif

  rdopt->min_rdcost = input->rdopt ? min_rdcost : min_cost;
  
#ifndef SIMPLIFY_CODE   
  if ( (img->MbaffFrameFlag)
    && (img->current_mb_nr%2)
    && (currMB->mb_type ? 0:((bslice) ? !currMB->cbp:1))  // bottom is skip
    && (prevMB->mb_type ? 0:((bslice) ? !prevMB->cbp:1))
    && !(field_flag_inference() == enc_mb.curr_mb_field)) // top is skip
  {    
    rdopt->min_rdcost = 1e30;  // don't allow coding of a MB pair as skip if wrong inference
  }
  
  
  //===== Decide if this MB will restrict the reference frames =====
  if (input->RestrictRef)
    update_refresh_map(intra, intra1, currMB);  
  

  if(input->FMEnable == 1)
  {
    skip_intrabk_SAD(best_mode, listXsize[enc_mb.list_offset[LIST_0]]);
  }
  else if(input->FMEnable == 2)
  {
    simplified_skip_intrabk_SAD(best_mode, listXsize[enc_mb.list_offset[LIST_0]]);
  }
  

  //--- constrain intra prediction ---
  if(input->UseConstrainedIntraPred && (img->type==P_SLICE || img->type==B_SLICE))
  {
    if( !IS_INTRA(currMB) )
    {
      img->intra_block[img->current_mb_nr] = 0;
    }
    else
    {
      img->intra_block[img->current_mb_nr] = 1;
    }
    
  }
#endif
}
double encode_macroblock_cluster32(int CurrentMbAddr, int mb32count, int saved_prev_qp, int saved_prev_delta_qp, int saved_cluster_qp, int *delta_qp_sent)
{
    double rdcost32=INT_MAX, rdcost16=0;
    int i, j;
    int encodeMbAddr;
    int CurrentMbAddr_x, CurrentMbAddr_y, extend_x, extend_y, num_MB, nextMbAddr, mbcount=0;

//#ifdef MB32_DELTA_QP
    int saved_delta_qp_sent16, saved_delta_qp_sent32, saved_delta_qp_sent64;
    saved_delta_qp_sent16=saved_delta_qp_sent32 = saved_delta_qp_sent64 = *delta_qp_sent;
//#endif

    find_Mb_cluster(1, CurrentMbAddr, &CurrentMbAddr_x, &CurrentMbAddr_y, &extend_x, &extend_y, &num_MB, &nextMbAddr);

    /*****************
     *check 16x16 mode
     ****************/
    //MB32_16x16.mb_type = 4;
    //MB32_16x16.cbp = 0;
    
    MB32.mb_type = P8x8;
    for(i=CurrentMbAddr_y; i<=extend_y; i++)
    {
      for(j=CurrentMbAddr_x; j<=extend_x; j++)
      {
        encodeMbAddr = i*img->PicWidthInMbs + j;
        start_macroblock (encodeMbAddr, FALSE);
        /* Same qp is used for 32x32 block
           All inside 16x16 blocks share the same prev_qp, prev_delta_qp and delta_qp 
         */
        assign_qp_info(mbcount,  encodeMbAddr, saved_cluster_qp, saved_prev_qp, saved_prev_delta_qp);  //overwrite what have been done in start_macroblock for qp info
        if(mbcount == 0 && num_MB == 4)
          write_mb32_mbtype = 1;
        else
          write_mb32_mbtype=0;

#ifdef RDO_Q        
        curr_mbi = mb32count*4 + mbcount;
#endif

// SPEEDUP
        clusterof16_isbetterthan_one32_duringRDOQ = 0;    //default value

        encode_one_macroblock ();

        write_mb32_mbtype=0;

        if(img->type==P_SLICE && (img->mb_data[encodeMbAddr]).mb_type == 0 )
          assert(img->mb_data[encodeMbAddr].cbp ==0);

        //if(img->mb_data[encodeMbAddr].cbp != 0)
        //{          
        //  MB32_16x16.cbp = 1;
        //}
        //mb32_16x16_cbp[mbcount] = img->mb_data[encodeMbAddr].cbp;

        rdcost16 += rdopt->min_rdcost;
        //total_cost 
//#ifdef MB32_DELTA_QP
        //qp handling
        if(img->mb_data[encodeMbAddr].cbp > 0 || img->mb_data[encodeMbAddr].mb_type == I16MB)
        {
          *delta_qp_sent = saved_delta_qp_sent16 = 1;
        }
//#endif

// SPEEDUP
        saveMode.MBsize16[mb32count*4 + mbcount]   = img->mb_data[encodeMbAddr].mb_type;
        saveMode.MBsize32[mb32count] = P8x8;

        mbcount++;
      }//j
    }//i
    bestmode32=P8x8;


    /**************************
     *check 32x32, 32x16, 16x32
     *************************/    

    if(num_MB == 4)
    {
      //qp handling
      *delta_qp_sent = saved_delta_qp_sent32;//MB32_DELTA_QP
     
      start_macroblock32 (CurrentMbAddr, FALSE, 1);

      //assert(img->mb_data[encodeMbAddr].qp == saved_cluster_qp);   
      //assert(img->mb_data[encodeMbAddr].prev_qp == saved_prev_qp);
      //assert(img->mb_data[encodeMbAddr].prev_delta_qp == saved_prev_delta_qp); //in writeDquant_CABAC, last_dquant=currMB->prev_delta_qp


      assign_qp_info(mbcount,  CurrentMbAddr, saved_cluster_qp, saved_prev_qp, saved_prev_delta_qp);//overwrite what have been done in start_macroblock for qp info   

#ifdef RDO_Q        
        curr_mb32i = mb32count;
#endif

      encode_one_macroblock32 (1);
      rdcost32 = rdopt->min_rdcost;

      //if(CurrentMbAddr == 162)
      //  printf("mb32 mb_nr %3d, mode %4d rdcost32 %f\n\n", CurrentMbAddr, img->mb_data[img->current_mb_nr].mb_type, rdcost32);

//#ifdef MB32_DELTA_QP
       if(MB32.cbp !=0 )
       {
          saved_delta_qp_sent32 = 1;
       }
       if(rdcost32<rdcost16)
       {
         *delta_qp_sent = saved_delta_qp_sent32;
       }
       else
       {
         *delta_qp_sent = saved_delta_qp_sent16;
       }
//#endif

      if(rdcost32 < rdcost16)
      {
        bestmode32=img->mb_data[img->current_mb_nr].mb_type; //img->current_mb_nr point to the upper left 16x16 block of 32x32 block


        //if(bestmode32 == 2 || bestmode32 == 3)
        //  printf("good\n");
// SPEEDUP
        clusterof16_isbetterthan_one32_duringRDOQ = 0;
        saveMode.MBsize32[mb32count]        = bestmode32;
        saveMode.MBsize16[mb32count*4 + 0]  = -1;
        saveMode.MBsize16[mb32count*4 + 1]  = -1;
        saveMode.MBsize16[mb32count*4 + 2]  = -1;
        saveMode.MBsize16[mb32count*4 + 3]  = -1;
        
        return rdcost32;
      }
      else
      {
        if(input->UseExtMB == 2)//MB64X64
        {
          // we have to redo mode decision to rebuild the correct context for the remaining blocks inside the 64x64
          *delta_qp_sent = saved_delta_qp_sent64;
          mbcount = 0;

          for(i=CurrentMbAddr_y; i<=extend_y; i++)
          {
            for(j=CurrentMbAddr_x; j<=extend_x; j++)
            {
              encodeMbAddr = i*img->PicWidthInMbs + j;
              start_macroblock (encodeMbAddr, FALSE);
              /* Same qp is used for 32x32 block
                 All inside 16x16 blocks share the same prev_qp, prev_delta_qp and delta_qp 
               */
              assign_qp_info(mbcount,  encodeMbAddr, saved_cluster_qp, saved_prev_qp, saved_prev_delta_qp);  //overwrite what have been done in start_macroblock for qp info

              if(mbcount == 0 && num_MB == 4)
                write_mb32_mbtype = 1;
              else
                write_mb32_mbtype=0;

#ifdef RDO_Q        
              curr_mbi = mb32count*4 + mbcount;
#endif
// SPEEDUP
            clusterof16_isbetterthan_one32_duringRDOQ = 1;

              encode_one_macroblock ();
// SPEEDUP
            clusterof16_isbetterthan_one32_duringRDOQ = 0;

              write_mb32_mbtype=0;

//#ifdef MB32_DELTA_QP
              if(img->mb_data[encodeMbAddr].cbp > 0 || img->mb_data[encodeMbAddr].mb_type == I16MB)
              {
                *delta_qp_sent = saved_delta_qp_sent64 = 1;
              }
//#endif
// SPEEDUP
        saveMode.MBsize16[mb32count*4 + mbcount]   = img->mb_data[encodeMbAddr].mb_type;
        saveMode.MBsize32[mb32count] = P8x8;

              mbcount++;
            }//j
          }//i
        }

      }
    }
    return rdcost16;
}

double encode_macroblock_cluster64(int CurrentMbAddr, int saved_prev_qp, int saved_prev_delta_qp, int saved_cluster_qp, int *delta_qp_sent)//MB64X64
{
    double rdcost64=INT_MAX, rdcost32=0;

    int i, j;
    int encodeMbAddr;
    int CurrentMbAddr_x, CurrentMbAddr_y, extend_x, extend_y, num_MB, nextMbAddr, mb32count=0;

    find_Mb_cluster(2, CurrentMbAddr, &CurrentMbAddr_x, &CurrentMbAddr_y, &extend_x, &extend_y, &num_MB, &nextMbAddr);


    /*****************
     *check 32x32 mode
     ****************/

    //MB32_16x16.mb_type = 4;
    //MB32_16x16.cbp = 0;
    
    MB64.mb_type = P8x8;
    for(i=CurrentMbAddr_y; i<=extend_y; i+=2)
    {
      for(j=CurrentMbAddr_x; j<=extend_x; j+=2)
      {
        encodeMbAddr = i*img->PicWidthInMbs + j;
        start_macroblock (encodeMbAddr, FALSE);
        /* Same qp is used for 32x32 block
           All inside 16x16 blocks share the same prev_qp, prev_delta_qp and delta_qp 
         */
          //in RDOQ, currMB->qp and currMB->delta_qp may change if cbp=0;
        assign_qp_info(mb32count,  encodeMbAddr, saved_cluster_qp, saved_prev_qp, saved_prev_delta_qp); //overwrite what have been done in start_macroblock for qp info

        if(mb32count == 0 && num_MB == 16)
          write_mb64_mbtype = 1;
        else
          write_mb64_mbtype=0;

        rdcost32 += encode_macroblock_cluster32 (encodeMbAddr, mb32count, saved_prev_qp, saved_prev_delta_qp, saved_cluster_qp, delta_qp_sent);

        write_mb64_mbtype=0;

        if(i+1<=extend_y && j+1<=extend_x)
          mb32type[mb32count] = bestmode32;
        mb32count++;
      }//j
    }//i
    bestmode64=P8x8;

// SPEEDUP
    saveMode.MBsize64 = P8x8;


    /**************************
     *check 64x64, 64x32, 32x64
     *************************/    

    if(num_MB == 16)
    {
      *delta_qp_sent = 0;//MB32_DELTA_QP
     
      start_macroblock32 (CurrentMbAddr, FALSE, 2);
      assign_qp_info(mb32count,  CurrentMbAddr, saved_cluster_qp, saved_prev_qp, saved_prev_delta_qp);//overwrite what have been done in start_macroblock for qp info
      encode_one_macroblock32 (2);
      rdcost64 = rdopt->min_rdcost;
 
      if(rdcost64 < rdcost32)
      {
        bestmode64=img->mb_data[img->current_mb_nr].mb_type; //img->current_mb_nr point to the upper left 16x16 block of 32x32 block
// SPEEDUP
        saveMode.MBsize64 = bestmode64;
        memset(saveMode.MBsize32, -1 ,  4 * sizeof(int));
        memset(saveMode.MBsize16, -1 , 16 * sizeof(int));

        return rdcost64;
      }        
    }

    return rdcost32;
}

static void ComputeResidue (imgpel **curImg, int img_m7[16][16], int mb_y, int mb_x, int opix_y, int opix_x, int width, int height)
{
  static imgpel *imgOrg, *imgPred;
  static int    *m7;
  int i, j;

  for (j = mb_y; j < mb_y + height; j++)
  {
    imgOrg = &curImg[opix_y + j][opix_x];    
    imgPred = &(img->mpr[j][mb_x]);
    m7 = &img_m7[j][mb_x]; 
    for (i = 0; i < width; i++)
    {
      *m7++ = *imgOrg++ - *imgPred++;
    }
  }
}

void SetModesAndRefframe (int b8, short* p_dir, int* fw_mode, int* bw_mode, short* fw_ref, short* bw_ref);

void LumaResidualCodingBigBlock()
{
  int block8x8;
  int fw_mode, bw_mode;
  short p_dir, refframe, bw_ref;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int nonzero; 
  int pic_pix_y, pic_pix_x;
  int byy, block_y;
  int bxx, block_x;
  int mb_x, mb_y;

  currMB->cbp     = 0 ;
  currMB->cbp_blk = 0 ;

  assert(currMB->luma_transform_size_8x8_flag != 0);

  for (block8x8=0; block8x8<4; block8x8++)
  {    
    mb_y       = (block8x8 >> 1) << 3;
    mb_x       = (block8x8 & 0x01) << 3;

    SetModesAndRefframe (block8x8, &p_dir, &fw_mode, &bw_mode, &refframe, &bw_ref);

    if (((p_dir == 0 || p_dir == 2 )&& fw_mode < 5) || ((p_dir == 1 || p_dir == 2 ) && bw_mode < 5))
    {
      for (byy=0, block_y=mb_y; block_y<mb_y+8; byy+=4, block_y+=4)
      {
        pic_pix_y = img->opix_y + block_y;
        
        for (bxx=0, block_x=mb_x; block_x<mb_x+8; bxx+=4, block_x+=4)
        {
          pic_pix_x = img->opix_x + block_x;

          //===== prediction of 4x4 block =====
          LumaPrediction4x4 (block_x, block_y, p_dir, fw_mode, bw_mode, refframe, bw_ref);
        }
      }
      
      //===== compute prediction residual ======            
      ComputeResidue (imgY_org, img->m7, mb_y, mb_x, img->opix_y, img->opix_x + mb_x, 8, 8);
    }
  }
  if( (currMB->mb_type < 2 && currMB->luma_transform_size_8x8_flag == 2) || only_16x16_transform==1  )
  {
      nonzero = dct_luma_16x16_new();
      if (nonzero)
      {
        currMB->cbp_blk = 0xffff; // all sixteen 4x4 blocks contain coeff
        currMB->cbp = 0xf;        // alll 8x8 blocks contain coeff
      }
  }
  else if(currMB->mb_type == 2 && currMB->luma_transform_size_8x8_flag == 2)    // 16x8 transform
  {
    // top block
    nonzero = dct_luma_16x8(0);
    if (nonzero)
    {
      currMB->cbp_blk |= 0x00ff; // yye: double check
      currMB->cbp |= 0x3;        // yye: double check
    }
    // bottom block
    nonzero = dct_luma_16x8(1);
    if (nonzero)
    {
      currMB->cbp_blk |= 0xff00; // yye: double check
      currMB->cbp |= 0xc;        // yye: double check
    }
  }
  else if(currMB->mb_type == 3 && currMB->luma_transform_size_8x8_flag == 2)    // 8x16 transform
  {
    // top block
    nonzero = dct_luma_8x16(0);
    if (nonzero)
    {
      currMB->cbp_blk |= 0x3333; // yye: double check
      currMB->cbp |= 0x5;        // yye: double check
    }
    // bottom block
    nonzero = dct_luma_8x16(1);
    if (nonzero)
    {
      currMB->cbp_blk |= 0xcccc; // yye: double check
      currMB->cbp |= 0xa;        // yye: double check
    }
  }

}

#endif

