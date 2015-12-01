/*!
**************************************************************************************
* \file
*    slice.c
* \brief
*    generate the slice header, setup the bit buffer for slices,
*    and generates the slice NALU(s)

* \author
*    Main contributors (see contributors.h for copyright, address and affiliation details)
*      - Thomas Stockhammer            <stockhammer@ei.tum.de>
*      - Detlev Marpe                  <marpe@hhi.de>
*      - Stephan Wenger                <stewe@cs.tu-berlin.de>
*      - Alexis Michael Tourapis       <alexismt@ieee.org>
***************************************************************************************
*/

#include "contributors.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "global.h"
#include "header.h"
#include "rtp.h"
#include "fmo.h"
#include "vlc.h"
#include "image.h"
#include "cabac.h"
#include "elements.h"
#include "epzs.h"
#ifdef ADAPTIVE_FILTER
#include "adaptive_filter.h"
#endif
#ifdef ADAPTIVE_FD_SD_CODING
#include "rdopt_coding_state.h"
CSptr cs_slice_coding  = NULL;
CSptr cs_slice_coding1 = NULL;
#endif
#ifdef MV_COMPETITION
#include "mv_competition.h"
#endif
#ifdef ADAPTIVE_QUANTIZATION
#include "adaptive_quantization.h"
void InitAQMSStatParamForSlice();
void encode_one_macroblock_for_iaqms(int CurrentMbAddr, double *min_rdcost, int *best_iaqms_idx, int *saved_best_mode, int *best_QP, int *best_deltaQP, int *best_athr, int *best_cbp, int masterQP, int deltaQP);
int internal_param[10];
void init_internal_parameter(int CurrentMbAddr, int *best_iaqms_idx, int *saved_best_mode, int *best_QP, int *best_deltaQP, int *best_athr, int *best_cbp);
void set_internal_parameter();
void decide_last_coeff(int CurrentMbAddr, double *min_rdcost, int *best_iaqms_idx, int *saved_best_mode, int *best_QP, int *best_deltaQP, int *best_athr, int *best_cbp, int masterQP, int deltaQP);
#endif
#if (defined(ADAPTIVE_QUANTIZATION) || defined(RDO_Q))
#ifdef MB32X32
int mb32_16x16_cbp[4]; //each 16x16 block cbp of 32x32 block
int write_mb32_mbtype=0;
int mb32_mvd[2][BLOCK_MULTIPLE<<MAX_MB_EXT_LEVEL][BLOCK_MULTIPLE<<MAX_MB_EXT_LEVEL][2];          //!< indices correspond to [forw,backw][block_y][block_x][x,y]
extern int bestmode32;
extern Macroblock  MB32, MB32_16x16;
void SetFrameOffset32(int CurrentMbAddr);
int delta_qp_sent;
int eos_bit_32=1;
double encode_macroblock_cluster64(int CurrentMbAddr, int saved_cluster_qp, int saved_prev_qp, int saved_prev_delta_qp, int *delta_qp_sent);
int write_macroblock_cluster(int *NumberOfCodedMBs, int CurrentMbAddr);
double encode_macroblock_cluster32(int CurrentMbAddr, int mb32count, int saved_prev_qp, int saved_prev_delta_qp, int saved_cluster_qp, int *delta_qp_sent);
int write_macroblock_cluster32(int *NumberOfCodedMBs, int CurrentMbAddr, int mb32count, 
                               int saved_prev_qp, int saved_prev_delta_qp, int saved_cluster_qp, int *delta_qp_sent);
void SetAdaptiveFilter32(int CurrentMbAddr);
double encode_one_macroblock_each_quant(int CurrentMbAddr);
#else
void encode_one_macroblock_each_quant(int CurrentMbAddr);
#endif
#endif
#ifdef USE_INTRA_MDDT
#include <memory.h>
static void InitScanOrderForSlice();
extern void precompute_all_inner_product8x8();
extern void precompute_all_inner_product16x16();
#endif 
#ifdef SWITCHED_FILTERS
#include "switched_filters.h"
#endif  // SWITCHED_FILTERS

// Local declarations
static Slice *malloc_slice();
static void  free_slice(Slice *slice);
static void  init_slice(int start_mb_addr);
static void set_ref_pic_num();
extern ColocatedParams *Co_located;
extern StorablePicture **listX[6];
void poc_ref_pic_reorder(StorablePicture **list, unsigned num_ref_idx_lX_active, int *reordering_of_pic_nums_idc, int *abs_diff_pic_num_minus1, int *long_term_pic_idx, int weighted_prediction, int list_no);
void SetLagrangianMultipliers();
#ifdef RDO_Q
extern short best_mode;
#ifdef MB32X32
int curr_mbi, curr_mb32i;
#endif
#endif

/*!
************************************************************************
* \brief
*    init_ref_pic_list_reordering initializations should go here
************************************************************************
*/
void init_ref_pic_list_reordering()
{
  Slice* currSlice = img->currentSlice;

  currSlice->ref_pic_list_reordering_flag_l0 = 0;
  currSlice->ref_pic_list_reordering_flag_l1 = 0;
}


/*!
************************************************************************
*  \brief
*     This function generates the slice (and partition) header(s) 
*
*  \return number of bits used for the slice (and partition) header(s)
*
*  \par Side effects:
*      Adds slice/partition header symbols to the symbol buffer
*      increments Picture->no_slices, allocates memory for the
*      slice, sets img->currSlice
************************************************************************
*/
int start_slice()
{
  EncodingEnvironmentPtr eep;
  Slice *currSlice = img->currentSlice;
  Bitstream *currStream;
  int header_len = 0;
  int i;
  int NumberOfPartitions = (input->partition_mode == PAR_DP_1?1:3);

  //one  partition for IDR img
  if(img->currentPicture->idr_flag)
  {
    NumberOfPartitions = 1;
  }

  RTPUpdateTimestamp (img->tr);   // this has no side effects, just leave it for all NALs

  for (i=0; i<NumberOfPartitions; i++)
  {
    currStream = (currSlice->partArr[i]).bitstream;

    currStream->write_flag = 0;
    if (i==0)     // First partition
      header_len += SliceHeader (0);
    else          // Second/Third partition
      header_len += Partition_BC_Header(i);
    //! Initialize CABAC
    if (input->symbol_mode == CABAC)
    {
      eep = &((currSlice->partArr[i]).ee_cabac);
      if (currStream->bits_to_go != 8)
        header_len+=currStream->bits_to_go;
      writeVlcByteAlign(currStream);
      arienco_start_encoding(eep, currStream->streamBuffer, &(currStream->byte_pos));
      cabac_new_slice();
    } 
    else 
    {
      // Initialize CA-VLC
      CAVLC_init();
    }
  }
  if(input->symbol_mode == CABAC)
  {
    init_contexts();
  }
  return header_len;
}



/*!
************************************************************************
* \brief
*    This function terminates a slice (but doesn't write it out), 
*    the old terminate_slice (0)
* \return
*    0 if OK,                                                         \n
*    1 in case of error
*
************************************************************************
*/
int terminate_slice(int lastslice)
{
  static int MbWidthC  [4]= { 0, 8, 8,  16};
  static int MbHeightC [4]= { 0, 8, 16, 16};

  int bytes_written;
  Bitstream *currStream;
  Slice *currSlice = img->currentSlice;
  EncodingEnvironmentPtr eep;
  int i;
  int byte_pos_before_startcode_emu_prevention;
  int min_num_bytes=0;
  int stuffing_bytes=0;
  int RawMbBits;

  if (input->symbol_mode == CABAC)
    write_terminating_bit (1);      // only once, not for all partitions

#ifdef ADAPTIVE_LOOP_FILTER
  if(!input->UseAdaptiveLoopFilter)
  {
#endif

  for (i=0; i<currSlice->max_part_nr; i++)
  {
    currStream = (currSlice->partArr[i]).bitstream;
    if (input->symbol_mode == UVLC)
    {
      SODBtoRBSP(currStream);
      byte_pos_before_startcode_emu_prevention = currStream->byte_pos;
      currStream->byte_pos = RBSPtoEBSP(currStream->streamBuffer, 0 , currStream->byte_pos, 0);
      *(stats->em_prev_bits) += (currStream->byte_pos - byte_pos_before_startcode_emu_prevention) * 8;
    }
    else     // CABAC
    {
      eep = &((currSlice->partArr[i]).ee_cabac);
      // terminate the arithmetic code
      arienco_done_encoding(eep);
      currStream->bits_to_go = eep->Ebits_to_go;
      currStream->byte_buf = 0;
      bytes_written = currStream->byte_pos;
      img->bytes_in_picture += currStream->byte_pos;

      byte_pos_before_startcode_emu_prevention= currStream->byte_pos;
      if (lastslice && i==((currSlice->max_part_nr-1)))
      {
        RawMbBits = 256 * img->bitdepth_luma + 2 * MbWidthC[active_sps->chroma_format_idc] * MbHeightC[active_sps->chroma_format_idc] * img->bitdepth_chroma;
        min_num_bytes = ((96 * get_pic_bin_count()) - (RawMbBits * (int)img->PicSizeInMbs *3) + 1023) / 1024;
        if (min_num_bytes>img->bytes_in_picture)
        {
          stuffing_bytes = min_num_bytes - img->bytes_in_picture;
          printf ("CABAC stuffing words = %6d\n", stuffing_bytes/3);
        }
      }

      //      printf ("bytepos: %d\n", currStream->byte_pos);
      currStream->byte_pos = RBSPtoEBSP(currStream->streamBuffer, 0, currStream->byte_pos, currStream->byte_pos + stuffing_bytes);
      *(stats->em_prev_bits) += (currStream->byte_pos - byte_pos_before_startcode_emu_prevention) * 8;
    }           // CABAC
  }           // partition loop

#ifdef ADAPTIVE_LOOP_FILTER
  }
#endif

  if( input->symbol_mode == CABAC )
  {
    store_contexts();
  }

  if (img->type != I_SLICE || img->type != SI_SLICE)
    free_ref_pic_list_reordering_buffer (currSlice);
  return 0;   
}

/*!
************************************************************************
* \brief
*    Encodes one slice
* \par
*   returns the number of coded MBs in the SLice 
************************************************************************
*/
int encode_one_slice (int SliceGroupId, Picture *pic, int TotalCodedMBs)
{
  Boolean end_of_slice = FALSE;
  Boolean recode_macroblock;
  int len;
  int NumberOfCodedMBs = 0;
  int CurrentMbAddr;
  double FrameRDCost = DBL_MAX, FieldRDCost = DBL_MAX;
#ifdef RDO_Q
  int masterQP=0, deltaQP, bestQP;
  double min_rdcost=1e30; 
  int qp_left, qp_up;
  int best_deltaQP=0;
#ifdef MB32X32
  double rdcost_per_qp=0;
  int saved_CurrentMbAddr;
#endif
#endif
#ifdef ADAPTIVE_FD_SD_CODING
  int i,j;
  int slice_bits;
  int rate_only_FD;
  int rate_FD_and_SD;
  int SSD_only_FD;
  int SSD_FD_and_SD;
  double rd_cost_only_FD;
  double rd_cost_FD_and_SD;
  Macroblock*     currMB;
#endif

#ifdef ADAPTIVE_QUANTIZATION
  int best_iaqms_idx = 0;
  int best_cbp = 47;
  int best_athr = 0;
#endif

#ifdef USE_NEW_OFFSET
  if(input->UseNewOffset)
    ResetFrameOffset();
#endif

#ifdef ADAPTIVE_FILTER
  if(img->AdaptiveFilterFlag == 0)
    ResetAdaptiveFilter(); 
  else if (img->AdaptiveFilterFlag == 1)
  {
    if (input->UseAdaptiveFilter == 1)
      UnifiedOneForthPixWithNewFilter();
    else if (input->UseAdaptiveFilter == 2)
      UnifiedOneForthPixWithNewFilterVer(); // separable aif
#ifdef DIRECTIONAL_FILTER
    else if (input->UseAdaptiveFilter == 3)
    {
      if (input->ImpType == IMP_FLOAT32)
        UnifiedOneForthPixWith_1DAIF_float();
      else if (input->ImpType == IMP_INT16)
        UnifiedOneForthPixWith_1DAIF_int16();
    }
#ifdef E_DAIF
    else if (input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)
    {
#ifdef EDAIF2
      UnifiedOneForthPixWithNewFilter();
#else
      UnifiedOneForthPixWith_EDAIF();
#endif
    }
#endif  // E_DAIF
#ifdef EAIF
    else if (input->UseAdaptiveFilter == FILTER_TYPE_EAIF)
    {
      UnifiedOneForthPixWith_EAIF();
    }
#endif
#endif  // DIRECTIONAL_FILTER

    SwapUpsampledFrames();
  }
#endif


#ifdef ADAPTIVE_QUANTIZATION
  InitAQMSStatParamForSlice();
#endif

#ifdef ADAPTIVE_FD_SD_CODING
  rdopt_FDSD=&rddata_FDSD_coding;
  rdopt_FDSD_interlace=&rddata_FDSD_coding_interlace;
#endif

  img->cod_counter = 0;

#ifdef USE_INTRA_MDDT
  if(input->UseIntraMDDT)
  {
    InitScanOrderForSlice(); 
    precompute_all_inner_product8x8();
    precompute_all_inner_product16x16(); 
  }
#endif
  
#ifdef MB32X32
  if (input->UseExtMB!=0 && img->type != I_SLICE)
  {
    calcCoeff16x16I(img->dct_coeff16I);
    calcCoeff8x8I  (img->dct_coeff8I );
  }
#endif

  CurrentMbAddr = FmoGetFirstMacroblockInSlice (SliceGroupId);
  // printf ("\n\nEncode_one_slice: PictureID %d SliceGroupId %d  SliceID %d  FirstMB %d \n", img->tr, SliceGroupId, img->current_slice_nr, CurrentMbInScanOrder);

  init_slice (CurrentMbAddr);

  // <FTRD : Compatibility with hierarchical B slices
#ifdef MV_COMPETITION
  if(input->mv_competition > 0)
  {
    enc_picture->slice_type = img->type;
    if(img->type == B_SLICE /*&& input->HierarchicalCoding!=0*/) 
      init_mv_scale_hb();
  }
#endif
  // FTRD>

  Bytes_After_Header = img->currentSlice->partArr[0].bitstream->byte_pos;

  SetLagrangianMultipliers();

#ifdef SWITCHED_FILTERS
  if(img->type != I_SLICE)
  {    
    if((input->UseHPFilter == HPF_SIFO) || (input->UseHPFilter == HPF_SIFO_FPO))
    {
      if(img->filterParam == SIFO_FIRST_PASS_FPO)
      {
        if(img->type == B_SLICE)
        {
          setFirstPassSubpelOffset(LIST_0);
          setFirstPassSubpelOffset(LIST_1);
        }
        else
        {
          setFirstPassSubpelOffset(LIST_0);
        }
      }

      UnifiedOneForthPixFiltSel(img->filterParam);
      SwapFilteredFrames();
    }
  }
#endif  // SWITCHED_FILTERS

  if (input->symbol_mode==CABAC)
  {
    SetCtxModelNumber ();
#ifdef ADAPTIVE_LOOP_FILTER
    if (input->UseAdaptiveLoopFilter)
    {
      pic->slices[pic->no_slices-1]->model_number = img->model_number;
    }
#endif
  }

  img->checkref = (input->rdopt && input->RestrictRef && (img->type==P_SLICE || img->type==SP_SLICE));

  /*
  // Tian Dong: June 7, 2002 JVT-B042
  // When the pictures are put into different layers and subseq, not all the reference frames
  // in multi-frame buffer are valid for prediction. The acutual number of the valid reference
  // frames, fb->num_short_used, will be given by start_slice(sym).
  // Save the fb->short_used.
  if (input->NumFramesInELSubSeq)
  {
  short_used = fb->short_used;
  }
  */

  len = start_slice ();
  // Rate control
  img->NumberofHeaderBits +=len;

  // basic unit layer rate control
  if(img->BasicUnit<img->Frame_Total_Number_MB)
    img->NumberofBasicUnitHeaderBits +=len;

  //  printf("short size, used, num-used: (%d,%d,%d)\n", fb->short_size, fb->short_used, fb->num_short_used);

  /*
  // Tian Dong: June 7, 2002 JVT-B042
  if (input->NumFramesInELSubSeq)
  {
  fb->short_used = fb->num_short_used;
  }
  */
  // Update statistics
  stats->bit_slice += len;
  stats->bit_use_header[img->type] += len;
  // printf ("\n\n");

#ifdef USE_INTRA_MDDT
  if(input->UseIntraMDDT)
  {
    memset(img->quant_stat,    0, sizeof(long)*16);
    memset(img->quant_stat8x8, 0, sizeof(long)*64);
    memset(img->quant_stat16x16, 0, sizeof(long)*256);
  }
#endif 

  while (end_of_slice == FALSE) // loop over macroblocks
  {
#ifdef MB32X32
    if(input->UseExtMB!=0)
    {
      memset(saveBestMode.MBsize16, -1 , 16 * sizeof(int));
      memset(saveBestMode.MBsize32, -1 ,  4 * sizeof(int));
      saveBestMode.MBsize64 = -1;
      clusterof16_isbetterthan_one32_duringRDOQ   = 0;
      clusterof16_isbetterthan_one32_finalwriting = 0;
    }
#endif
    if (img->AdaptiveRounding && input->AdaptRndPeriod && (img->current_mb_nr % input->AdaptRndPeriod == 0))
    {
      CalculateOffsetParam();

      if(input->Transform8x8Mode)
      {
        CalculateOffset8Param();
      }
    }

    //sw paff
    if (!img->MbaffFrameFlag)
    {
#ifdef RDO_Q
      static int deltaQPTabB[] = {0, 1, -1, 2,  3};
      static int deltaQPTabP[] = {0, 1, -1, 2, -2};
      static int deltaQPTabSizeB = 5;
      static int deltaQPTabSizeP = 5; 
      int        deltaQPCnt; 
      int master_cr_qp=0, master_cr_cbp=0, cr_qp;

#ifdef ADAPTIVE_QUANTIZATION
      init_internal_parameter(CurrentMbAddr, &best_iaqms_idx, &saved_best_mode, &bestQP, &best_deltaQP, &best_athr, &best_cbp);
#endif
#endif
      recode_macroblock = FALSE;
      rdopt = &rddata_top_frame_mb;   // store data in top frame MB 

#ifndef RDO_Q
      start_macroblock (CurrentMbAddr, FALSE);
#else
#ifdef ADAPTIVE_QUANTIZATION
      if(input->UseRDO_Q || img->slice_fractional_quant_flag)
#else
      if(input->UseRDO_Q)
#endif
      {
        masterQP=img->masterQP=img->qp;

        Motion_Selected = 0;
        final_mb_encoding = 0;

        if(img->type == B_SLICE)
        {
          int qp_anchor; 
          Macroblock *currMB;
          start_macroblock (CurrentMbAddr, FALSE);

          currMB = &img->mb_data[img->current_mb_nr];
          qp_left = (currMB->mb_available_left) ? currMB->mb_available_left->qp:img->masterQP;
          qp_up = (currMB->mb_available_up) ? currMB->mb_available_up->qp:img->masterQP;
          qp_anchor = (qp_left+qp_up+1)>>1;
          min_rdcost = 1e30; 
          for (deltaQPCnt=0; deltaQPCnt<deltaQPTabSizeB; deltaQPCnt++)
          {
            deltaQP = deltaQPTabB[deltaQPCnt];
            img->qp=Clip3(-img->bitdepth_luma_qp_scale, 51, (masterQP+deltaQP));

            // skip this QP if deltaQP is too large based on local stats 
            if(!(img->qp-qp_anchor >= -1 && img->qp-qp_anchor <= 2) && currMB->mb_available_left && currMB->mb_available_up)
              continue;

#ifdef ADAPTIVE_QUANTIZATION
            if(!input->UseRDO_Q && deltaQP!=0) continue;
            if(img->slice_fractional_quant_flag)
            {
              encode_one_macroblock_for_iaqms(CurrentMbAddr, &min_rdcost, &best_iaqms_idx, &saved_best_mode, &bestQP, &best_deltaQP, &best_athr, &best_cbp, masterQP, deltaQP);
            }
            else
#endif
#ifdef MB32X32
            rdcost_per_qp = encode_one_macroblock_each_quant(CurrentMbAddr);
#else
            encode_one_macroblock_each_quant(CurrentMbAddr);
#endif       
            Motion_Selected = 1;

#ifdef MB32X32
            if(input->UseExtMB!=0)
            {
              if (rdcost_per_qp<min_rdcost)
              {
                min_rdcost=rdcost_per_qp;
                bestQP=img->qp;
                best_deltaQP = deltaQP;
                //saved_best_mode = bestmode32;//img->mb_data[img->current_mb_nr].mb_type;//best_mode;
  //SPEEDUP
                saveBestMode.MBsize64 = saveMode.MBsize64;
                memcpy(saveBestMode.MBsize32, saveMode.MBsize32,  4*sizeof(int));
                memcpy(saveBestMode.MBsize16, saveMode.MBsize16, 16*sizeof(int));
                saved_best_mode = saveBestMode.MBsize64;

              }
            }
            else
            {
#endif

              if (rdopt->min_rdcost<min_rdcost)
              {
                min_rdcost=rdopt->min_rdcost;
                bestQP=img->qp;
                best_deltaQP = deltaQP;
                saved_best_mode = img->mb_data[img->current_mb_nr].mb_type;//best_mode;
              }
#ifdef MB32X32
            }
#endif
          }
        }
        else
        {   // I and P slices
          int qp_anchor; 
          Macroblock *currMB; 

          start_macroblock (CurrentMbAddr, FALSE);

          currMB = &img->mb_data[img->current_mb_nr];
          qp_left = (currMB->mb_available_left) ? currMB->mb_available_left->qp:img->masterQP;
          qp_up = (currMB->mb_available_up) ? currMB->mb_available_up->qp:img->masterQP;
          qp_anchor = (qp_left+qp_up+1)>>1;

          min_rdcost = 1e30;
          for (deltaQPCnt=0; deltaQPCnt<deltaQPTabSizeP; deltaQPCnt++)
          {
            deltaQP = deltaQPTabP[deltaQPCnt];
            img->qp=Clip3(-img->bitdepth_luma_qp_scale, 51, (masterQP+deltaQP));

            if(deltaQP != 0 && !(img->qp-qp_anchor >= -2 && img->qp-qp_anchor <= 1) && currMB->mb_available_left && currMB->mb_available_up && img->type == P_SLICE)
              continue; 

#ifdef ADAPTIVE_QUANTIZATION
            if(!input->UseRDO_Q && deltaQP!=0) continue;
            if(img->slice_fractional_quant_flag)
            {
              encode_one_macroblock_for_iaqms(CurrentMbAddr, &min_rdcost, &best_iaqms_idx, &saved_best_mode, &bestQP, &best_deltaQP, &best_athr, &best_cbp, masterQP, deltaQP);
            }
            else
#endif
#ifdef MB32X32
              rdcost_per_qp = encode_one_macroblock_each_quant(CurrentMbAddr);
#else         
              encode_one_macroblock_each_quant(CurrentMbAddr);
#endif

            Motion_Selected = 1;
#ifdef ADAPTIVE_QUANTIZATION
            if(input->successive_Bframe==0 && input->UseRDO_Q)
#else
            if(input->successive_Bframe==0)
#endif
            {
              static byte QP_SCALE_CR[52]=
              {
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,
                12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,
                28,29,29,30,31,32,32,33,34,34,35,35,36,36,37,37,
                37,38,38,38,39,39,39,39
              };

              cr_qp = Clip3 ( -img->bitdepth_chroma_qp_scale, 51, img->qp + img->chroma_qp_offset[0] );
              cr_qp = (cr_qp < 0 ? cr_qp : QP_SCALE_CR[cr_qp]) + img->bitdepth_chroma_qp_scale;

              if(deltaQP == 0) // deltaQP =0 has to be test first!
              {
                master_cr_cbp = img->mb_data[img->current_mb_nr].cbp/16;
                master_cr_qp = cr_qp;
              }
              else
              { 
                if(deltaQP>0 && cr_qp > master_cr_qp && master_cr_cbp!=0 )
                  continue;
              }
            }

#ifdef MB32X32
            if(input->UseExtMB!=0)
            {
              if (rdcost_per_qp<min_rdcost)
              {
                min_rdcost=rdcost_per_qp;
                bestQP=img->qp;
                best_deltaQP = deltaQP;
                //saved_best_mode = bestmode32; //img->mb_data[img->current_mb_nr].mb_type;//best_mode;
  //SPEEDUP
                saveBestMode.MBsize64 = saveMode.MBsize64;
                memcpy(saveBestMode.MBsize32, saveMode.MBsize32,  4*sizeof(int));
                memcpy(saveBestMode.MBsize16, saveMode.MBsize16, 16*sizeof(int));
                saved_best_mode = saveBestMode.MBsize64;

              }
            }
            else
            {
#endif
              if (rdopt->min_rdcost<min_rdcost)
              {
                min_rdcost=rdopt->min_rdcost;
                bestQP=img->qp;
                best_deltaQP = deltaQP;
                saved_best_mode = img->mb_data[img->current_mb_nr].mb_type;//best_mode;
              }
#ifdef MB32X32
            }
#endif
          }
        }

        final_mb_encoding = 1;
        img->qp=Clip3(-img->bitdepth_luma_qp_scale, 51, (masterQP + best_deltaQP));
#ifdef ADAPTIVE_QUANTIZATION
        img->mb_iaqms_idx = best_iaqms_idx;
        input->disthres  = best_athr;
#endif
        
#ifdef MB32X32
        if(input->UseExtMB == 0)
#endif
          encode_one_macroblock_each_quant(CurrentMbAddr);

      }
      else
        start_macroblock (CurrentMbAddr, FALSE);
#endif

#ifdef ADAPTIVE_FD_SD_CODING
      currMB    = &img->mb_data[img->current_mb_nr];

      if (img->APEC_in_FD_and_SD==0)
      {
        currMB->SD_Coding_on_off=0;
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
        if(!input->UseRDO_Q && !(img->slice_fractional_quant_flag) )
#else
        if(!input->UseRDO_Q)
#endif
        {
#ifdef MB32X32
          if(input->UseExtMB > 0)
            encode_one_macroblock_each_quant(CurrentMbAddr);
          else
#endif  
            encode_one_macroblock ();
        }
#else
        encode_one_macroblock ();
#endif

#ifdef USE_NEW_OFFSET
#ifdef MB32X32
        if(input->UseNewOffset && input->UseExtMB==0)
#else
        if(input->UseNewOffset)
#endif
        {
          if(img->frameOffsetAvail == 0)
            SetFrameOffset(); 
        }
#endif

#ifdef ADAPTIVE_FILTER
#ifdef MB32X32
        if(input->UseExtMB==0  && UseAdaptiveFilterForCurrentFrame() && img->AdaptiveFilterFlag == 0)
#else
        if(UseAdaptiveFilterForCurrentFrame() && img->AdaptiveFilterFlag == 0)
#endif
        {
#ifdef EAIF
          if ((input->UseAdaptiveFilter == 1) || (input->UseAdaptiveFilter == 3) || (input->UseAdaptiveFilter == FILTER_TYPE_EDAIF) || (input->UseAdaptiveFilter == FILTER_TYPE_EAIF))
#else
          if ((input->UseAdaptiveFilter == 1) || (input->UseAdaptiveFilter == 3) || (input->UseAdaptiveFilter == FILTER_TYPE_EDAIF))
#endif
            SetAdaptiveFilter();
          else if (input->UseAdaptiveFilter == 2)
            SetAdaptiveFilterHor(); // separable aif
        }
#endif

#ifdef MB32X32
        if(input->UseExtMB>0)
        {
#ifdef RDO_Q
          if(input->UseRDO_Q)
          {
            Motion_Selected = 0;
            final_mb_encoding = 0;
          }
#endif

          clusterof16_isbetterthan_one32_duringRDOQ = 0;   //make it zero before Final encoding pass which will write the bitstream


          saved_CurrentMbAddr = CurrentMbAddr;

          CurrentMbAddr = write_macroblock_cluster (&NumberOfCodedMBs, CurrentMbAddr);//the return value is the next to-be-coded Mbaddr

          if(input->UseNewOffset && img->type != I_SLICE)
          {
            if(img->frameOffsetAvail == 0)
              SetFrameOffset32(saved_CurrentMbAddr);
          }
#ifdef ADAPTIVE_FILTER
          if(UseAdaptiveFilterForCurrentFrame() && img->AdaptiveFilterFlag == 0)
          {
              SetAdaptiveFilter32(saved_CurrentMbAddr);
          }
#endif

        }
        else
#endif
        write_one_macroblock (1);

#ifdef RDO_Q
        if(input->UseRDO_Q)
          img->qp=masterQP;
#endif
      }
      else
      {
        store_coding_state(cs_slice_coding1);
        slice_bits=stats->bit_slice;
        currMB->SD_Coding_on_off=1;
#ifdef ADAPTIVE_QUANTIZATION
        if(img->slice_fractional_quant_flag)
          encode_one_macroblock_each_quant(CurrentMbAddr);
        else
#endif
          encode_one_macroblock ();

        store_rdopt_data_for_FDSD_coding ();//copy all relevant data to rdopt_FDSD buffer

        write_one_macroblock (1);
        rate_only_FD=stats->bit_slice-slice_bits;
        SSD_only_FD=0;
        for (i=0;i<16;i++)
        {
          for (j=0;j<16;j++)
          {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
            SSD_only_FD+=SQR_DEPTH(enc_picture->imgY[img->pix_y+j][img->pix_x+i], imgY_org[img->opix_y+j][img->opix_x+i], input->BitDepthLuma, img->BitDepthIncrease);
#else
            SSD_only_FD+=img->quad[enc_picture->imgY[img->pix_y+j][img->pix_x+i]-imgY_org[img->opix_y+j][img->opix_x+i]];
#endif
          }
        }
        slice_bits=stats->bit_slice;
        currMB->SD_Coding_on_off=0;//Turn SD Coding off
        reset_coding_state(cs_slice_coding1);
#ifdef ADAPTIVE_QUANTIZATION
        if(img->slice_fractional_quant_flag)
          encode_one_macroblock_each_quant(CurrentMbAddr);
        else
#endif
          encode_one_macroblock ();

        write_one_macroblock (1);
        rate_FD_and_SD=stats->bit_slice-slice_bits;
        SSD_FD_and_SD=0;
        for (i=0;i<16;i++)
        {
          for (j=0;j<16;j++)
          {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
            SSD_FD_and_SD+=SQR_DEPTH(enc_picture->imgY[img->pix_y+j][img->pix_x+i], imgY_org[img->opix_y+j][img->opix_x+i], input->BitDepthLuma, img->BitDepthIncrease);
#else
            SSD_FD_and_SD+=img->quad[enc_picture->imgY[img->pix_y+j][img->pix_x+i]-imgY_org[img->opix_y+j][img->opix_x+i]];
#endif
          }
        }

        rd_cost_only_FD  =(double)SSD_only_FD  +img->lambda_md[img->type][img->qp] * (double)rate_only_FD;
        rd_cost_FD_and_SD=(double)SSD_FD_and_SD+img->lambda_md[img->type][img->qp] * (double)rate_FD_and_SD;

        if (rd_cost_only_FD<rd_cost_FD_and_SD)
        {
          reset_coding_state(cs_slice_coding1);
          slice_bits=stats->bit_slice;
          currMB->SD_Coding_on_off=1;

          restore_rdopt_data_for_FDSD_coding ();

          write_one_macroblock (1);
        }
#ifdef ADAPTIVE_FILTER
        if(UseAdaptiveFilterForCurrentFrame() && img->AdaptiveFilterFlag == 0)
        {
#ifdef EAIF
          if ((input->UseAdaptiveFilter == 1) || (input->UseAdaptiveFilter == 3) || (input->UseAdaptiveFilter == FILTER_TYPE_EDAIF) || (input->UseAdaptiveFilter == FILTER_TYPE_EAIF))
#else
          if ((input->UseAdaptiveFilter == 1) || ((input->UseAdaptiveFilter == 3) || (input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)))
#endif
            SetAdaptiveFilter();
          else
            SetAdaptiveFilterHor(); // separable aif
        }
#endif
      }
      if (currMB->written_SD_Coding_on_off==0)
      {
        currMB->SD_Coding_on_off=0;
      }
#else//ADAPTIVE_FD_SD_CODING
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
      if(!input->UseRDO_Q && !(img->slice_fractional_quant_flag) )
#else
      if(!input->UseRDO_Q)
#endif
        encode_one_macroblock ();
#else
      encode_one_macroblock ();
#endif
#ifdef ADAPTIVE_FILTER
      if(UseAdaptiveFilterForCurrentFrame() && img->AdaptiveFilterFlag == 0)
      {
        if ((input->UseAdaptiveFilter == 1) || ((input->UseAdaptiveFilter == 3) || (input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)))
          SetAdaptiveFilter();
        else
          SetAdaptiveFilterHor(); // separable aif
      }
#endif
      write_one_macroblock (1);

#ifdef RDO_Q
      if(input->UseRDO_Q)
        img->qp=masterQP;
#endif

#endif //ADAPTIVE_FD_SD_CODING
#ifdef ADAPTIVE_QUANTIZATION
      set_internal_parameter();
#endif

      terminate_macroblock (&end_of_slice, &recode_macroblock);

      //       printf ("encode_one_slice: mb %d,  slice %d,   bitbuf bytepos %d EOS %d\n", 
      //       img->current_mb_nr, img->current_slice_nr, 
      //       img->currentSlice->partArr[0].bitstream->byte_pos, end_of_slice);

      if (recode_macroblock == FALSE)       // The final processing of the macroblock has been done
      {
#ifdef MB32X32
        if(input->UseExtMB == 0)
#endif
          CurrentMbAddr = FmoGetNextMBNr (CurrentMbAddr); //taken care of by write_macroblock_cluster

        if (CurrentMbAddr == -1)   // end of slice
        {
          //          printf ("FMO End of Slice Group detected, current MBs %d, force end of slice\n", NumberOfCodedMBs+1);
          end_of_slice = TRUE;
#ifdef USE_NEW_OFFSET
          if(input->UseNewOffset)
          {
            if(img->frameOffsetAvail == 0)
            {
              CalcFrameOffset();
              img->frameOffsetAvail = 1; 
            }
          }
#endif
#ifdef ADAPTIVE_FILTER
          if(UseAdaptiveFilterForCurrentFrame() && img->AdaptiveFilterFlag == 0)
          {
#ifdef EAIF
            if((input->UseAdaptiveFilter == 1) || (input->UseAdaptiveFilter == 3) || (input->UseAdaptiveFilter == FILTER_TYPE_EAIF))
                FindFilterCoef();
#else
            if((input->UseAdaptiveFilter == 1) || (input->UseAdaptiveFilter == 3) || (input->UseAdaptiveFilter == FILTER_TYPE_EDAIF))
#endif
#ifdef EDAIF2
            else if(input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)
            {
              extern int g_nBitsThresh;
              g_nBitsThresh = stats->bit_slice;

              PreProcessFilterCoef();
              {
                SwapAIF(2,0);
                FindFilterCoef2();
                SwapAIF(2,1);
              }
            }
#endif
            else
              FindFilterCoefHor(); // separable aif
          }           
#endif    
        }
#ifdef MB32X32
        if(input->UseExtMB == 0)
#endif
        {
          NumberOfCodedMBs++;       // only here we are sure that the coded MB is actually included in the slice
          proceed2nextMacroblock ();
        }
      }
      else
      {
        //!Go back to the previous MB to recode it
        img->current_mb_nr = FmoGetPreviousMBNr(img->current_mb_nr);
        if(img->current_mb_nr == -1 )   // The first MB of the slice group  is too big,
          // which means it's impossible to encode picture using current slice bits restriction
        {
          snprintf (errortext, ET_SIZE, "Error encoding first MB with specified parameter, bits of current MB may be too big");
          error (errortext, 300);
        }
      }
    }
    else                      // TBD -- Addition of FMO
    {

      //! This following ugly code breaks slices, at least for a slice mode that accumulates a certain
      //! number of bits into one slice.  
      //! The suggested algorithm is as follows:
      //!
      //! SaveState (Bitstream, stats,  etc. etc.);
      //! BitsForThisMBPairInFrameMode = CodeMB (Upper, FRAME_MODE) + CodeMB (Lower, FRAME_MODE);
      //! DistortionForThisMBPairInFrameMode = CalculateDistortion(Upper) + CalculateDistortion (Lower);
      //! RestoreState();
      //! BitsForThisMBPairInFieldMode = CodeMB (Upper, FIELD_MODE) + CodeMB (Lower, FIELD_MODE);
      //! DistortionForThisMBPairInFrameMode = CalculateDistortion(Upper) + CalculateDistortion (Lower);
      //! FrameFieldMode = Decision (...)
      //! RestoreState()
      //! if (FrameFieldMode == FRAME) {
      //!   CodeMB (Upper, FRAME); CodeMB (Lower, FRAME);
      //! } else {
      //!   CodeMB (Upper FIELD); CodeMB (Lower, FIELD);
      //! }
      //!
      //! Open questions/issues:
      //!   1. CABAC/CA-VLC state:  It seems that the CABAC/CA_VLC states are changed during the
      //!      dummy encoding processes (for the R-D based selection), but that they are never
      //!      reset, once the selection is made.  I believe that this breaks the MB-adaptive
      //!      frame/field coding.  The necessary code for the state saves is readily available
      //!      in macroblock.c, start_macroblock() and terminate_macroblock() (this code needs
      //!      to be double checked that it works with CA-VLC as well
      //!   2. would it be an option to allocate Bitstreams with zero data in them (or copy the
      //!      already generated bitstream) for the "test coding"?  

      if (input->MbInterlace == ADAPTIVE_CODING)
      {
        //================ code MB pair as frame MB ================
        //----------------------------------------------------------
        recode_macroblock = FALSE;


        img->field_mode = 0;  // MB coded as frame
        img->top_field = 0;   // Set top field to 0

        //Rate control
        img->write_macroblock = 0;
        img->bot_MB = 0;   

        start_macroblock (CurrentMbAddr, FALSE);

        rdopt = &rddata_top_frame_mb; // store data in top frame MB
#ifdef ADAPTIVE_FD_SD_CODING
        currMB    = &img->mb_data[CurrentMbAddr];
        if (img->APEC_in_FD_and_SD==0)
        {
          currMB->SD_Coding_on_off=0;//Turn SD Coding off
          encode_one_macroblock ();
        }
        else
        {
          store_coding_state(cs_slice_coding1);
          slice_bits=stats->bit_slice;
          currMB->SD_Coding_on_off=1;//Turn SD Coding on
          encode_one_macroblock ();

          store_rdopt_data_for_FDSD_coding ();
          store_rdopt_data_for_FDSD_coding_interlace ();

          write_one_macroblock (1);     // write the Top MB data to the bitstream
          rate_only_FD=stats->bit_slice-slice_bits;
          SSD_only_FD=0;
          for (i=0;i<16;i++)
          {
            for (j=0;j<16;j++)
            {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
              SSD_only_FD+=SQR_DEPTH(enc_picture->imgY[img->pix_y+j][img->pix_x+i], imgY_org[img->opix_y+j][img->opix_x+i], input->BitDepthLuma, img->BitDepthIncrease);
#else
              SSD_only_FD+=img->quad[enc_picture->imgY[img->pix_y+j][img->pix_x+i]-imgY_org[img->opix_y+j][img->opix_x+i]];
#endif
            }
          }
          slice_bits=stats->bit_slice;
          reset_coding_state(cs_slice_coding1);
          currMB->SD_Coding_on_off=0;
          encode_one_macroblock ();
          write_one_macroblock (1);     // write the Top MB data to the bitstream
          rate_FD_and_SD=stats->bit_slice-slice_bits;
          reset_coding_state(cs_slice_coding1);
          SSD_FD_and_SD=0;
          for (i=0;i<16;i++)
          {
            for (j=0;j<16;j++)
            {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
              SSD_FD_and_SD+=SQR_DEPTH(enc_picture->imgY[img->pix_y+j][img->pix_x+i], imgY_org[img->opix_y+j][img->opix_x+i], input->BitDepthLuma, img->BitDepthIncrease);
#else
              SSD_FD_and_SD+=img->quad[enc_picture->imgY[img->pix_y+j][img->pix_x+i]-imgY_org[img->opix_y+j][img->opix_x+i]];
#endif
            }
          }
          rd_cost_only_FD  =(double)SSD_only_FD  +img->lambda_md[img->type][img->qp] * (double)rate_only_FD;
          rd_cost_FD_and_SD=(double)SSD_FD_and_SD+img->lambda_md[img->type][img->qp] * (double)rate_FD_and_SD;
          if (rd_cost_only_FD<rd_cost_FD_and_SD)
          {
            restore_rdopt_data_for_FDSD_coding();
            restore_rdopt_data_for_FDSD_coding_interlace ();
            slice_bits=stats->bit_slice;
            currMB->SD_Coding_on_off=1;
          }
        }
#else
        encode_one_macroblock ();     // code the MB as frame
#endif
        FrameRDCost = rdopt->min_rdcost;

        //***   Top MB coded as frame MB ***//

        //Rate control
        img->bot_MB = 1; //for Rate control

        // go to the bottom MB in the MB pair
        img->field_mode = 0;  // MB coded as frame  //GB

        start_macroblock (CurrentMbAddr+1, FALSE);
        rdopt = &rddata_bot_frame_mb; // store data in top frame MB
#ifdef ADAPTIVE_FD_SD_CODING
        currMB    = &img->mb_data[CurrentMbAddr+1];
        if (img->APEC_in_FD_and_SD==0)
        {
          currMB->SD_Coding_on_off=0;//Turn SD Coding off
          encode_one_macroblock ();
        }
        else
        {
          store_coding_state(cs_slice_coding1);
          slice_bits=stats->bit_slice;
          currMB->SD_Coding_on_off=1;//Turn SD Coding on
          encode_one_macroblock ();
          store_rdopt_data_for_FDSD_coding();
          store_rdopt_data_for_FDSD_coding_interlace ();
          write_one_macroblock (0);     // write the Bottom MB data to the bitstream
          rate_only_FD=stats->bit_slice-slice_bits;
          SSD_only_FD=0;
          for (i=0;i<16;i++)
          {
            for (j=0;j<16;j++)
            {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
              SSD_only_FD+=SQR_DEPTH(enc_picture->imgY[img->pix_y+j][img->pix_x+i], imgY_org[img->opix_y+j][img->opix_x+i], input->BitDepthLuma, img->BitDepthIncrease);
#else
              SSD_only_FD+=img->quad[enc_picture->imgY[img->pix_y+j][img->pix_x+i]-imgY_org[img->opix_y+j][img->opix_x+i]];
#endif
            }
          }
          slice_bits=stats->bit_slice;
          reset_coding_state(cs_slice_coding1);
          currMB->SD_Coding_on_off=0;
          encode_one_macroblock ();
          write_one_macroblock (0);     // write the Bottom MB data to the bitstream
          reset_coding_state(cs_slice_coding1);
          rate_FD_and_SD=stats->bit_slice-slice_bits;
          SSD_FD_and_SD=0;
          for (i=0;i<16;i++)
          {
            for (j=0;j<16;j++)
            {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
              SSD_FD_and_SD+=SQR_DEPTH(enc_picture->imgY[img->pix_y+j][img->pix_x+i], imgY_org[img->opix_y+j][img->opix_x+i], input->BitDepthLuma, img->BitDepthIncrease);
#else
              SSD_FD_and_SD+=img->quad[enc_picture->imgY[img->pix_y+j][img->pix_x+i]-imgY_org[img->opix_y+j][img->opix_x+i]];
#endif
            }
          }
          rd_cost_only_FD  =(double)SSD_only_FD  +img->lambda_md[img->type][img->qp] * (double)rate_only_FD;
          rd_cost_FD_and_SD=(double)SSD_FD_and_SD+img->lambda_md[img->type][img->qp] * (double)rate_FD_and_SD;
          if (rd_cost_only_FD<rd_cost_FD_and_SD)
          {
            restore_rdopt_data_for_FDSD_coding();
            restore_rdopt_data_for_FDSD_coding_interlace ();
            slice_bits=stats->bit_slice;
            currMB->SD_Coding_on_off=1;
          }
        }
#else
        encode_one_macroblock ();     // code the MB as frame
#endif
        FrameRDCost += rdopt->min_rdcost;

        //***   Bottom MB coded as frame MB ***//
      }

      if ((input->MbInterlace == ADAPTIVE_CODING) || (input->MbInterlace == FIELD_CODING))
      {
        //Rate control
        img->bot_MB = 0;

        //=========== start coding the MB pair as a field MB pair =============
        //---------------------------------------------------------------------
        img->field_mode = 1;  // MB coded as field
        img->top_field = 1;   // Set top field to 1
        img->buf_cycle <<= 1;
        input->num_ref_frames <<= 1;
        img->num_ref_idx_l0_active <<= 1;
        img->num_ref_idx_l0_active += 1;
        start_macroblock (CurrentMbAddr, TRUE);


        rdopt = &rddata_top_field_mb; // store data in top frame MB
        //        TopFieldIsSkipped = 0;        // set the top field MB skipped flag to 0
#ifdef ADAPTIVE_FD_SD_CODING
        currMB    = &img->mb_data[CurrentMbAddr];
        if (img->APEC_in_FD_and_SD==0)
        {
          currMB->SD_Coding_on_off=0;//Turn SD Coding off
          encode_one_macroblock ();
        }
        else
        {
          store_coding_state(cs_slice_coding1);
          slice_bits=stats->bit_slice;
          currMB->SD_Coding_on_off=1;//Turn SD Coding on
          encode_one_macroblock ();
          store_rdopt_data_for_FDSD_coding();
          store_rdopt_data_for_FDSD_coding_interlace ();

          write_one_macroblock (1);
          rate_only_FD=stats->bit_slice-slice_bits;
          reset_coding_state(cs_slice_coding1);
          SSD_only_FD=0;
          for (i=0;i<16;i++)
          {
            for (j=0;j<16;j++)
            {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
              SSD_only_FD+=SQR_DEPTH(enc_picture->imgY[img->pix_y+j][img->pix_x+i], imgY_org[img->opix_y+j][img->opix_x+i], input->BitDepthLuma, img->BitDepthIncrease);
#else
              SSD_only_FD+=img->quad[enc_picture->imgY[img->pix_y+j][img->pix_x+i]-imgY_org[img->opix_y+j][img->opix_x+i]];
#endif
            }
          }
          slice_bits=stats->bit_slice;
          currMB->SD_Coding_on_off=0;
          encode_one_macroblock ();
          write_one_macroblock (1);
          rate_FD_and_SD=stats->bit_slice-slice_bits;
          reset_coding_state(cs_slice_coding1);
          SSD_FD_and_SD=0;
          for (i=0;i<16;i++)
          {
            for (j=0;j<16;j++)
            {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
              SSD_FD_and_SD+=SQR_DEPTH(enc_picture->imgY[img->pix_y+j][img->pix_x+i], imgY_org[img->opix_y+j][img->opix_x+i], input->BitDepthLuma, img->BitDepthIncrease);
#else
              SSD_FD_and_SD+=img->quad[enc_picture->imgY[img->pix_y+j][img->pix_x+i]-imgY_org[img->opix_y+j][img->opix_x+i]];
#endif
            }
          }
          rd_cost_only_FD  =(double)SSD_only_FD  +img->lambda_md[img->type][img->qp] * (double)rate_only_FD;
          rd_cost_FD_and_SD=(double)SSD_FD_and_SD+img->lambda_md[img->type][img->qp] * (double)rate_FD_and_SD;
          if (rd_cost_only_FD<rd_cost_FD_and_SD)
          {
            restore_rdopt_data_for_FDSD_coding ();
            restore_rdopt_data_for_FDSD_coding_interlace ();
            slice_bits=stats->bit_slice;
            currMB->SD_Coding_on_off=1;
          }
        }
#else
        encode_one_macroblock ();     // code the MB as frame
#endif
        FieldRDCost = rdopt->min_rdcost;
        //***   Top MB coded as field MB ***//
        //Rate control
        img->bot_MB = 1;//for Rate control

        img->top_field = 0;   // Set top field to 0
        start_macroblock (CurrentMbAddr+1, TRUE);
        rdopt = &rddata_bot_field_mb; // store data in top frame MB
#ifdef ADAPTIVE_FD_SD_CODING
        currMB    = &img->mb_data[CurrentMbAddr+1];
        if (img->APEC_in_FD_and_SD==0)
        {
          currMB->SD_Coding_on_off=0;//Turn SD Coding off
          encode_one_macroblock ();
        }
        else
        {
          store_coding_state(cs_slice_coding1);
          slice_bits=stats->bit_slice;
          currMB->SD_Coding_on_off=1;//Turn SD Coding on
          encode_one_macroblock ();
          store_rdopt_data_for_FDSD_coding ();
          store_rdopt_data_for_FDSD_coding_interlace ();
          write_one_macroblock (0);
          rate_only_FD=stats->bit_slice-slice_bits;
          reset_coding_state(cs_slice_coding1);
          SSD_only_FD=0;
          for (i=0;i<16;i++)
          {
            for (j=0;j<16;j++)
            {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
              SSD_only_FD+=SQR_DEPTH(enc_picture->imgY[img->pix_y+j][img->pix_x+i], imgY_org[img->opix_y+j][img->opix_x+i], input->BitDepthLuma, img->BitDepthIncrease);
#else
              SSD_only_FD+=img->quad[enc_picture->imgY[img->pix_y+j][img->pix_x+i]-imgY_org[img->opix_y+j][img->opix_x+i]];
#endif
            }
          }
          slice_bits=stats->bit_slice;
          currMB->SD_Coding_on_off=0;
          write_one_macroblock (0);
          encode_one_macroblock ();
          rate_FD_and_SD=stats->bit_slice-slice_bits;
          reset_coding_state(cs_slice_coding1);
          SSD_FD_and_SD=0;
          for (i=0;i<16;i++)
          {
            for (j=0;j<16;j++)
            {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
              SSD_FD_and_SD+=SQR_DEPTH(enc_picture->imgY[img->pix_y+j][img->pix_x+i], imgY_org[img->opix_y+j][img->opix_x+i], input->BitDepthLuma, img->BitDepthIncrease);
#else
              SSD_FD_and_SD+=img->quad[enc_picture->imgY[img->pix_y+j][img->pix_x+i]-imgY_org[img->opix_y+j][img->opix_x+i]];
#endif
            }
          }
          rd_cost_only_FD  =(double)SSD_only_FD  +img->lambda_md[img->type][img->qp] * (double)rate_only_FD;
          rd_cost_FD_and_SD=(double)SSD_FD_and_SD+img->lambda_md[img->type][img->qp] * (double)rate_FD_and_SD;
          if (rd_cost_only_FD<rd_cost_FD_and_SD)
          {
            restore_rdopt_data_for_FDSD_coding();
            restore_rdopt_data_for_FDSD_coding_interlace ();
            slice_bits=stats->bit_slice;
            currMB->SD_Coding_on_off=1;
          }
        }
#else
        encode_one_macroblock ();     // code the MB as frame
#endif
        FieldRDCost += rdopt->min_rdcost;
        //***   Bottom MB coded as field MB ***//
      }

      //Rate control
      img->write_macroblock_frame = 0;  //Rate control

      //=========== decide between frame/field MB pair ============
      //-----------------------------------------------------------
      if ((input->MbInterlace == ADAPTIVE_CODING) && (FrameRDCost < FieldRDCost))
      {
        img->field_mode = 0;
        img->buf_cycle >>= 1;
        input->num_ref_frames >>= 1;
        MBPairIsField = 0;
        img->num_ref_idx_l0_active -= 1;
        img->num_ref_idx_l0_active >>= 1;

        //Rate control
        img->write_macroblock_frame = 1;  //for Rate control
      }
      else
      {
        img->field_mode = 1;
        MBPairIsField = 1;
      }

      //Rate control
      img->write_macroblock = 1;//Rate control 

      if (MBPairIsField)
        img->top_field = 1;
      else
        img->top_field = 0;

      //Rate control
      img->bot_MB = 0;// for Rate control

      // go back to the Top MB in the MB pair
      start_macroblock (CurrentMbAddr, img->field_mode);

      rdopt =  img->field_mode ? &rddata_top_field_mb : &rddata_top_frame_mb;
      copy_rdopt_data (0);  // copy the MB data for Top MB from the temp buffers

#ifdef ADAPTIVE_FD_SD_CODING
      write_one_macroblock (1);     // write the Top MB data to the bitstream
      currMB    = &img->mb_data[CurrentMbAddr];
      if (currMB->written_SD_Coding_on_off==0)
      {
        currMB->SD_Coding_on_off=0;
      }
#else
      write_one_macroblock (1);     // write the Top MB data to the bitstream
#endif

      terminate_macroblock (&end_of_slice, &recode_macroblock);     // done coding the Top MB 

      if (recode_macroblock == FALSE)       // The final processing of the macroblock has been done
      {
        CurrentMbAddr = FmoGetNextMBNr (CurrentMbAddr);
        if (CurrentMbAddr == -1)   // end of slice
        {
          end_of_slice = TRUE;
        }
        NumberOfCodedMBs++;       // only here we are sure that the coded MB is actually included in the slice
        proceed2nextMacroblock ();


        //Rate control
        img->bot_MB = 1;//for Rate control
        // go to the Bottom MB in the MB pair
        img->top_field = 0;
        start_macroblock (CurrentMbAddr, img->field_mode);

        rdopt = img->field_mode ? &rddata_bot_field_mb : &rddata_bot_frame_mb;
        copy_rdopt_data (1);  // copy the MB data for Bottom MB from the temp buffers

#ifdef ADAPTIVE_FD_SD_CODING
        write_one_macroblock (0);     // write the Bottom MB data to the bitstream
        currMB    = &img->mb_data[CurrentMbAddr];
        if (currMB->written_SD_Coding_on_off==0)
        {
          currMB->SD_Coding_on_off=0;
        }
#else
        write_one_macroblock (0);     // write the Bottom MB data to the bitstream
#endif

        terminate_macroblock (&end_of_slice, &recode_macroblock);     // done coding the Top MB 
        if (recode_macroblock == FALSE)       // The final processing of the macroblock has been done
        {
          CurrentMbAddr = FmoGetNextMBNr (CurrentMbAddr);
          if (CurrentMbAddr == -1)   // end of slice
          {
            end_of_slice = TRUE;
          }
          NumberOfCodedMBs++;       // only here we are sure that the coded MB is actually included in the slice
          proceed2nextMacroblock ();
        }
        else
        {
          //Go back to the beginning of the macroblock pair to recode it
          img->current_mb_nr = FmoGetPreviousMBNr(img->current_mb_nr);
          img->current_mb_nr = FmoGetPreviousMBNr(img->current_mb_nr);
          if(img->current_mb_nr == -1 )   // The first MB of the slice group  is too big,
            // which means it's impossible to encode picture using current slice bits restriction
          {
            snprintf (errortext, ET_SIZE, "Error encoding first MB with specified parameter, bits of current MB may be too big");
            error (errortext, 300);
          }
        }

      }

      else
      {
        //!Go back to the previous MB to recode it
        img->current_mb_nr = FmoGetPreviousMBNr(img->current_mb_nr);
        if(img->current_mb_nr == -1 )   // The first MB of the slice group  is too big,
          // which means it's impossible to encode picture using current slice bits restriction
        {
          snprintf (errortext, ET_SIZE, "Error encoding first MB with specified parameter, bits of current MB may be too big");
          error (errortext, 300);
        }
      }



      if (MBPairIsField)    // if MB Pair was coded as field the buffer size variables back to frame mode
      {
        img->buf_cycle >>= 1;
        input->num_ref_frames >>= 1;
        img->num_ref_idx_l0_active -= 1;
        img->num_ref_idx_l0_active >>= 1;
      }

      img->field_mode = img->top_field = 0; // reset to frame mode

      if (CurrentMbAddr == FmoGetLastCodedMBOfSliceGroup (FmoMB2SliceGroup (CurrentMbAddr)))
        end_of_slice = TRUE;        // just in case it doesn't get set in terminate_macroblock
    }
  }  
  /*
  // Tian Dong: June 7, 2002 JVT-B042
  // Restore the short_used
  if (input->NumFramesInELSubSeq)
  {
  fb->short_used = short_used;
  }
  */
  // separable aif (BEGIN)
#ifdef ADAPTIVE_FILTER      
  if(UseAdaptiveFilterForCurrentFrame() && img->AdaptiveFilterFlag == 0 && input->UseAdaptiveFilter == 2)
  {
    Macroblock*     currMB;
    UnifiedOneForthPixWithNewFilterHor ();
    CurrentMbAddr = FmoGetFirstMacroblockInSlice (SliceGroupId);
    end_of_slice = FALSE;
    while (end_of_slice == FALSE) // loop over macroblocks
    {
      set_MB_parameters (CurrentMbAddr);
      currMB    = &img->mb_data[img->current_mb_nr];
      SetAdaptiveFilterVer();
      terminate_macroblock (&end_of_slice, &recode_macroblock);
      CurrentMbAddr = FmoGetNextMBNr (CurrentMbAddr);
      if (CurrentMbAddr == -1)   // end of slice
      {
        end_of_slice = TRUE;
        FindFilterCoefVer();
      }
      proceed2nextMacroblock ();
    }
  }
  // separable aif (END)
  if (img->AdaptiveFilterFlag == 1)
    SwapUpsampledFrames();

#endif

#ifdef SWITCHED_FILTERS
  if((img->type != I_SLICE) && ((input->UseHPFilter == HPF_SIFO) || (input->UseHPFilter == HPF_SIFO_FPO)))
    SwapFilteredFrames();
#endif  // SWITCHED_FILTERS

  terminate_slice ( (NumberOfCodedMBs+TotalCodedMBs >= (int)img->PicSizeInMbs) );
  return NumberOfCodedMBs;
}

/*!
************************************************************************
* \brief
*    Initializes the parameters for a new slice and
*     allocates the memory for the coded slice in the Picture structure
*  \par Side effects:
*      Adds slice/partition header symbols to the symbol buffer
*      increments Picture->no_slices, allocates memory for the
*      slice, sets img->currSlice
************************************************************************
*/
static void init_slice (int start_mb_addr)
{
  int i;
  Picture *currPic = img->currentPicture;
  DataPartition *dataPart;
  Bitstream *currStream;
  Slice *currSlice;

  img->current_mb_nr = start_mb_addr;

  // Allocate new Slice in the current Picture, and set img->currentSlice
  assert (currPic != NULL);
  currPic->no_slices++;

  if (currPic->no_slices >= MAXSLICEPERPICTURE)
    error ("Too many slices per picture, increase MAXSLICEPERPICTURE in global.h.", -1);
#ifdef ADAPTIVE_FD_SD_CODING
  if (currPic->slices[currPic->no_slices-1]==NULL)
#endif
    currPic->slices[currPic->no_slices-1] = malloc_slice();

  currSlice = currPic->slices[currPic->no_slices-1];

  img->currentSlice = currSlice;

  currSlice->picture_id = img->tr % 256;
  currSlice->qp = img->qp;
  currSlice->start_mb_nr = start_mb_addr;
  currSlice->slice_too_big = dummy_slice_too_big;

  for (i = 0; i < currSlice->max_part_nr; i++)
  {
    dataPart = &(currSlice->partArr[i]);
    if (input->symbol_mode == UVLC)
      dataPart->writeSyntaxElement = writeSyntaxElement_UVLC;
    else
      dataPart->writeSyntaxElement = writeSyntaxElement_CABAC;

    currStream = dataPart->bitstream;
    currStream->bits_to_go = 8;
    currStream->byte_pos = 0;
    currStream->byte_buf = 0;
  }

  img->num_ref_idx_l0_active = active_pps->num_ref_idx_l0_active_minus1 + 1; 
  img->num_ref_idx_l1_active = active_pps->num_ref_idx_l1_active_minus1 + 1;

  // primary and redundant slices: number of references overriding.
  if(input->redundant_pic_flag)
  {
    if(!redundant_coding)
    {
      img->num_ref_idx_l0_active = min(img->number,input->NumRefPrimary);
    }
    else
    {
      // 1 reference picture for redundant slices 
      img->num_ref_idx_l0_active = 1;   
    }
  }

  // code now also considers fields. Issue whether we should account this within the appropriate input params directly
  if ((img->type == P_SLICE || img->type == SP_SLICE) && input->P_List0_refs)
  {
    img->num_ref_idx_l0_active = min(img->num_ref_idx_l0_active, input->P_List0_refs * ((img->structure !=0) + 1));
  }
  if (img->type == B_SLICE )
  {
    if (input->B_List0_refs)
    {
      img->num_ref_idx_l0_active = min(img->num_ref_idx_l0_active, input->B_List0_refs * ((img->structure !=0) + 1));
    }
    if (input->B_List1_refs)
    {
      img->num_ref_idx_l1_active = min(img->num_ref_idx_l1_active, input->B_List1_refs * ((img->structure !=0) + 1));
    }
  }
  // generate reference picture lists
  init_lists(img->type, (PictureStructure)img->structure);

  // assign list 0 size from list size
  img->num_ref_idx_l0_active = listXsize[0];
  img->num_ref_idx_l1_active = listXsize[1];

  //Perform memory management based on poc distances for HierarchicalCoding
  if (img->nal_reference_idc  && input->HierarchicalCoding && input->PocMemoryManagement && dpb.ref_frames_in_buffer==active_sps->num_ref_frames)
  {    
    poc_based_ref_management(img->frame_num);
  }

  if (input->EnableOpenGOP)
  {
    for (i = 0; i<listXsize[0]; i++)
    {    
      if (listX[0][i]->poc < img->last_valid_reference && img->ThisPOC > img->last_valid_reference)      
      {
        listXsize[0] = img->num_ref_idx_l0_active = max(1,i);
        break;
      }
    }

    for (i = 0; i<listXsize[1]; i++)
    {
      if (listX[1][i]->poc < img->last_valid_reference && img->ThisPOC > img->last_valid_reference)
      {
        listXsize[1] = img->num_ref_idx_l1_active = max(1,i);
        break;
      }
    }
  }

  init_ref_pic_list_reordering();

  //Perform reordering based on poc distances for HierarchicalCoding
  //if (img->type==P_SLICE && input->HierarchicalCoding && input->ReferenceReorder)
  if (img->type==P_SLICE && input->ReferenceReorder)
  {

    int i, num_ref;

    alloc_ref_pic_list_reordering_buffer(currSlice);

    if ((img->type != I_SLICE) && (img->type !=SI_SLICE))
    {
      for (i=0; i<img->num_ref_idx_l0_active + 1; i++)
      {
        currSlice->reordering_of_pic_nums_idc_l0[i] = 3;
        currSlice->abs_diff_pic_num_minus1_l0[i] = 0;
        currSlice->long_term_pic_idx_l0[i] = 0;
      }

      if (img->type == B_SLICE)
      {
        for (i=0; i<img->num_ref_idx_l1_active + 1; i++)
        {
          currSlice->reordering_of_pic_nums_idc_l1[i] = 3;
          currSlice->abs_diff_pic_num_minus1_l1[i] = 0;
          currSlice->long_term_pic_idx_l1[i] = 0;
        }
      }
    }

    if ((img->type != I_SLICE) && (img->type !=SI_SLICE))
    {
      num_ref = img->num_ref_idx_l0_active;
      poc_ref_pic_reorder(listX[LIST_0], 
        num_ref, 
        currSlice->reordering_of_pic_nums_idc_l0, 
        currSlice->abs_diff_pic_num_minus1_l0, 
        currSlice->long_term_pic_idx_l0, 0, LIST_0);

      //reference picture reordering
      reorder_ref_pic_list(listX[LIST_0], &listXsize[LIST_0], 
        img->num_ref_idx_l0_active - 1, 
        currSlice->reordering_of_pic_nums_idc_l0, 
        currSlice->abs_diff_pic_num_minus1_l0, 
        currSlice->long_term_pic_idx_l0);

      // This is not necessary since order is already poc based...  
      if (img->type == B_SLICE)
      {
        num_ref = img->num_ref_idx_l1_active;
        poc_ref_pic_reorder(listX[LIST_1], 
          num_ref, 
          currSlice->reordering_of_pic_nums_idc_l1, 
          currSlice->abs_diff_pic_num_minus1_l1, 
          currSlice->long_term_pic_idx_l1, 0, LIST_1);

        //reference picture reordering
        reorder_ref_pic_list(listX[LIST_1], &listXsize[LIST_1], 
          img->num_ref_idx_l1_active - 1, 
          currSlice->reordering_of_pic_nums_idc_l1, 
          currSlice->abs_diff_pic_num_minus1_l1, 
          currSlice->long_term_pic_idx_l1);
      }
    }
  }


  //if (img->MbaffFrameFlag)
  if (img->structure==FRAME)
    init_mbaff_lists();

  if (img->type != I_SLICE && (active_pps->weighted_pred_flag == 1 || (active_pps->weighted_bipred_idc > 0 && (img->type == B_SLICE))))
  {
    if (img->type==P_SLICE || img->type==SP_SLICE)
    {
      if (input->GenerateMultiplePPS && input->RDPictureDecision)
      {
        if (enc_picture==enc_frame_picture2)
          estimate_weighting_factor_P_slice (0);
        else
          estimate_weighting_factor_P_slice (1);
      }
      else
        estimate_weighting_factor_P_slice (0);

    }
    else
      estimate_weighting_factor_B_slice ();
  }

  set_ref_pic_num();

  if (img->type == B_SLICE)
    compute_colocated(Co_located, listX);
  if (img->type != I_SLICE && input->FMEnable == 3)
    EPZSSliceInit(EPZSCo_located, listX);
}


/*!
************************************************************************
* \brief
*    Allocates a slice structure along with its dependent data structures
* \return
*    Pointer to a Slice
************************************************************************
*/
static Slice *malloc_slice()
{
  int i;
  DataPartition *dataPart;
  Slice *slice;

  //  const int buffer_size = (img->width * img->height * 4); // AH 190202: There can be data expansion with 
  // low QP values. So, we make sure that buffer 
  // does not overflow. 4 is probably safe multiplier.
  const int buffer_size = 500 + img->FrameSizeInMbs * (128 + 256 * img->bitdepth_luma + 512 * img->bitdepth_chroma);
  // KS: this is approx. max. allowed code picture size

  if ((slice = (Slice *) calloc(1, sizeof(Slice))) == NULL) no_mem_exit ("malloc_slice: slice structure");

  if (input->symbol_mode == CABAC)
  {
    // create all context models
    slice->mot_ctx = create_contexts_MotionInfo();
    slice->tex_ctx = create_contexts_TextureInfo();
  }

  slice->max_part_nr = input->partition_mode==0?1:3;

  //for IDR img there should be only one partition
  if(img->currentPicture->idr_flag)
    slice->max_part_nr = 1;

  assignSE2partition[0] = assignSE2partition_NoDP;
  //ZL 
  //for IDR img all the syntax element should be mapped to one partition        
  if(!img->currentPicture->idr_flag&&input->partition_mode==1)
    assignSE2partition[1] =  assignSE2partition_DP;
  else
    assignSE2partition[1] =  assignSE2partition_NoDP;        



  slice->num_mb = 0;          // no coded MBs so far

  if ((slice->partArr = (DataPartition *) calloc(slice->max_part_nr, sizeof(DataPartition))) == NULL) no_mem_exit ("malloc_slice: partArr");
  for (i=0; i<slice->max_part_nr; i++) // loop over all data partitions
  {
    dataPart = &(slice->partArr[i]);
    if ((dataPart->bitstream = (Bitstream *) calloc(1, sizeof(Bitstream))) == NULL) no_mem_exit ("malloc_slice: Bitstream");
    if ((dataPart->bitstream->streamBuffer = (byte *) calloc(buffer_size, sizeof(byte))) == NULL) no_mem_exit ("malloc_slice: StreamBuffer");
    // Initialize storage of bitstream parameters
  }
  return slice;
}


/*!
************************************************************************
* \brief
*    Memory frees of all Slice structures and of its dependent
*    data structures
* \par Input:
*    Image Parameters struct struct img_par *img
************************************************************************
*/
void free_slice_list(Picture *currPic)
{
  int i;

  for (i=0; i<currPic->no_slices; i++)
  {
    free_slice (currPic->slices[i]);
    currPic->slices[i]=NULL;
  }
}


/*!
************************************************************************
* \brief
*    Memory frees of the Slice structure and of its dependent
*    data structures
* \param slice:
*    Slice to be freed
************************************************************************
*/
static void free_slice(Slice *slice)
{
  int i;
  DataPartition *dataPart;

  if (slice != NULL)
  {
    for (i=0; i<slice->max_part_nr; i++) // loop over all data partitions
    {
      dataPart = &(slice->partArr[i]);
      if (dataPart != NULL)
      {
        if (dataPart->bitstream->streamBuffer != NULL)
          free(dataPart->bitstream->streamBuffer);
        if (dataPart->bitstream != NULL)
          free(dataPart->bitstream);
      }
    }
    if (slice->partArr != NULL)
      free(slice->partArr);
    if (input->symbol_mode == CABAC)
    {
      delete_contexts_MotionInfo(slice->mot_ctx);
      delete_contexts_TextureInfo(slice->tex_ctx);
    }

    free(slice);
  }
}

void set_ref_pic_num()
{
  int i,j;
  StorablePicture *this_ref;

  //! need to add field ref_pic_num that handles field pair.

  for (i=0;i<listXsize[LIST_0];i++)
  {
    this_ref = listX[LIST_0][i];
    enc_picture->ref_pic_num        [LIST_0][i] = this_ref->poc * 2 + ((this_ref->structure==BOTTOM_FIELD)?1:0) ; 
    enc_picture->frm_ref_pic_num    [LIST_0][i] = this_ref->frame_poc * 2; 
    enc_picture->top_ref_pic_num    [LIST_0][i] = this_ref->top_poc * 2; 
    enc_picture->bottom_ref_pic_num [LIST_0][i] = this_ref->bottom_poc * 2 + 1; 
  }

  for (i=0;i<listXsize[LIST_1];i++)
  {
    this_ref = listX[LIST_1][i];
    enc_picture->ref_pic_num        [LIST_1][i] = this_ref->poc  *2 + ((this_ref->structure==BOTTOM_FIELD)?1:0);
    enc_picture->frm_ref_pic_num    [LIST_1][i] = this_ref->frame_poc * 2; 
    enc_picture->top_ref_pic_num    [LIST_1][i] = this_ref->top_poc * 2; 
    enc_picture->bottom_ref_pic_num [LIST_1][i] = this_ref->bottom_poc * 2 + 1; 
  }

  if (!active_sps->frame_mbs_only_flag && img->structure==FRAME)
  {
    for (j=2;j<6;j++)
      for (i=0;i<listXsize[j];i++)
      {    
        this_ref = listX[j][i];
        enc_picture->ref_pic_num[j][i] = this_ref->poc * 2 + ((this_ref->structure==BOTTOM_FIELD)?1:0);
        enc_picture->frm_ref_pic_num[j][i] = this_ref->frame_poc * 2 ;
        enc_picture->top_ref_pic_num[j][i] = this_ref->top_poc * 2 ;
        enc_picture->bottom_ref_pic_num[j][i] = this_ref->bottom_poc * 2 + 1;
      }
  }
}

/*!
************************************************************************
* \brief
*    decide reference picture reordering, Frame only
************************************************************************
*/
void poc_ref_pic_reorder(StorablePicture **list, unsigned num_ref_idx_lX_active, int *reordering_of_pic_nums_idc, int *abs_diff_pic_num_minus1, int *long_term_pic_idx, int weighted_prediction, int list_no)
{
  unsigned i,j,k;

  int currPicNum, picNumLXPred;

  int default_order[32];
  int re_order[32];
  int tmp_reorder[32];
  int list_sign[32];
  int reorder_stop, no_reorder;
  int poc_diff[32];
  int tmp_value, diff;  

  int abs_poc_dist;
  int maxPicNum, MaxFrameNum = 1 << (log2_max_frame_num_minus4 + 4);

  if (img->structure==FRAME)
  {
    maxPicNum  = MaxFrameNum;
    currPicNum = img->frame_num;
  }
  else
  {
    maxPicNum  = 2 * MaxFrameNum;
    currPicNum = 2 * img->frame_num + 1;
  }

  picNumLXPred = currPicNum;

  // First assign default list order. 
  for (i=0; i<num_ref_idx_lX_active; i++)
  {
    default_order[i] = list[i]->pic_num;
  }

  // Now access all references in buffer and assign them
  // to a pottential reordering list. For each one of these 
  // references compute the poc distance compared to current
  // frame.
  for (i=0; i<dpb.ref_frames_in_buffer; i++)
  {
    re_order[i] = dpb.fs_ref[i]->frame->pic_num;

    if (dpb.fs_ref[i]->is_used==3 && (dpb.fs_ref[i]->frame->used_for_reference)&&(!dpb.fs_ref[i]->frame->is_long_term))
    {
      abs_poc_dist = abs(dpb.fs_ref[i]->frame->poc - enc_picture->poc) ;
      poc_diff[i] = abs_poc_dist;
      if (list_no == LIST_0)
      {
        list_sign[i] = (enc_picture->poc < dpb.fs_ref[i]->frame->poc) ? +1 : -1;
      }
      else
      {
        list_sign[i] = (enc_picture->poc > dpb.fs_ref[i]->frame->poc) ? +1 : -1;
      }
    }
  }


  // now sort these references based on poc (temporal) distance
  for (i=0; i< dpb.ref_frames_in_buffer-1; i++)
  {
    for (j=i+1; j< dpb.ref_frames_in_buffer; j++)
    {      
      if (poc_diff[i]>poc_diff[j] || (poc_diff[i] == poc_diff[j] && list_sign[j] > list_sign[i]))
      {

        tmp_value = poc_diff[i];
        poc_diff[i] = poc_diff[j];
        poc_diff[j] = tmp_value;
        tmp_value  = re_order[i];
        re_order[i] = re_order[j];
        re_order[j] = tmp_value ;
        tmp_value  = list_sign[i];
        list_sign[i] = list_sign[j];        
        list_sign[j] = tmp_value ;
      }
    }
  }

  // Check versus default list to see if any
  // change has happened
  no_reorder = 1;
  for (i=0; i<num_ref_idx_lX_active; i++)
  {
    if (default_order[i] != re_order[i])
    {
      no_reorder = 0;
    }
  }

  // If different, then signal reordering
  if (no_reorder==0)
  {
    for (i=0; i<num_ref_idx_lX_active; i++)
    {
      diff = re_order[i]-picNumLXPred;
      if (diff <= 0)
      {
        reordering_of_pic_nums_idc[i] = 0;
        abs_diff_pic_num_minus1[i] = abs(diff)-1;
        if (abs_diff_pic_num_minus1[i] < 0)
          abs_diff_pic_num_minus1[i] = maxPicNum -1; 
      }
      else
      {
        reordering_of_pic_nums_idc[i] = 1;
        abs_diff_pic_num_minus1[i] = abs(diff)-1;
      }
      picNumLXPred = re_order[i];

      tmp_reorder[i] = re_order[i];

      k = i;
      for (j=i; j<num_ref_idx_lX_active; j++)
      {
        if (default_order[j] != re_order[i])
        {
          ++k;
          tmp_reorder[k] = default_order[j];
        }
      }
      reorder_stop = 1;
      for(j=i+1; j<num_ref_idx_lX_active; j++)
      {
        if (tmp_reorder[j] != re_order[j])
        {
          reorder_stop = 0;
          break;
        }
      }

      if (reorder_stop==1)
      {
        ++i;
        break;
      }


      for(j=0; j<num_ref_idx_lX_active; j++)
      {
        default_order[j] = tmp_reorder[j];
      }

    }
    reordering_of_pic_nums_idc[i] = 3;

    for(j=0; j<num_ref_idx_lX_active; j++)
    {
      default_order[j] = tmp_reorder[j];
    }

    if (list_no==0)
    {
      img->currentSlice->ref_pic_list_reordering_flag_l0=1;
    }
    else
    {
      img->currentSlice->ref_pic_list_reordering_flag_l1=1;
    }
  }
}

extern int QP2QUANT[40];

void SetLagrangianMultipliers()
{
  int qp, j;
  double qp_temp;
  double lambda_scale = 1.0 - Clip3(0.0,0.5,0.05 * (double) input->jumpd);;

  if (input->rdopt) // RDOPT on computation of Lagrangian multipliers
  {
    for (j = 0; j < 5; j++)
    {
      for (qp = -img->bitdepth_luma_qp_scale; qp < 52; qp++)
      {          
        qp_temp = (double)qp + img->bitdepth_luma_qp_scale - SHIFT_QP;

        if (input->UseExplicitLambdaParams) // consideration of explicit weights.
        {
          img->lambda_md[j][qp] = input->LambdaWeight[j] * pow (2, qp_temp/3.0);
          // Scale lambda due to hadamard qpel only consideration
          img->lambda_md[j][qp] = (input->hadamard == 2 ? 0.95 : 1.00) * img->lambda_md[j][qp];
          img->lambda_me[j][qp] = sqrt(img->lambda_md[j][qp]);
          img->lambda_mf[j][qp] = LAMBDA_FACTOR (img->lambda_me[j][qp]);
          if (j == B_SLICE)
          {
            img->lambda_md[5][qp] = input->LambdaWeight[5] * pow (2, qp_temp/3.0);
            img->lambda_md[5][qp] = (input->hadamard == 2 ? 0.95 : 1.00) * img->lambda_md[5][qp];
            img->lambda_me[5][qp] = sqrt(img->lambda_md[5][qp]);
            img->lambda_mf[5][qp] = LAMBDA_FACTOR (img->lambda_me[5][qp]);
          }
        }
        else
        {
#ifdef ADAPTIVE_QUANTIZATION
          if(input->UseAdaptiveQuantMatrix && img->type == I_SLICE && img->qp==qp)
            img->lambda_md[j][qp] = 0.57 * pow (2.0, qp_temp/3.0); 
          else
#endif
#ifdef RDO_Q
            if(input->UseRDO_Q && img->type == I_SLICE && img->qp==qp)
              img->lambda_md[j][qp] = 0.57 * pow (2.0, qp_temp/3.0); 
            else
#endif
#if defined(USE_INTRA_MDDT)
            if(input->UseIntraMDDT && img->type == I_SLICE && img->qp==qp)
              img->lambda_md[j][qp] = 0.57 * pow (2.0, qp_temp/3.0); 
            else
#endif
              if (input->successive_Bframe>0)
                img->lambda_md[j][qp] = 0.68 * pow (2, qp_temp/3.0) 
                * (j == B_SLICE ? Clip3(2.00,4.00,(qp_temp / 6.0)) : (j == SP_SLICE) ? Clip3(1.4,3.0,(qp_temp / 12.0)) : 1.0);
              else
                img->lambda_md[j][qp] = 0.85 * pow (2, qp_temp/3.0) 
                * ( (j == B_SLICE) ? 4.0 : (j == SP_SLICE) ? Clip3(1.4,3.0,(qp_temp / 12.0)) : 1.0);
          // Scale lambda due to hadamard qpel only consideration
          img->lambda_md[j][qp] = (input->hadamard == 2 ? 0.95 : 1.00) * img->lambda_md[j][qp];
          img->lambda_md[j][qp] = (j == B_SLICE && input->BRefPictures == 2 && img->b_frame_to_code == 0 ? 0.50 : 1.00) * img->lambda_md[j][qp];

          if (j == B_SLICE)
          {
            img->lambda_md[5][qp] = img->lambda_md[j][qp];

            if (input->HierarchicalCoding == 2)
              img->lambda_md[5][qp] *= (1.0 - min(0.4,0.2 * (double) gop_structure[img->b_frame_to_code-1].hierarchy_layer)) ;
            else
              img->lambda_md[5][qp] *= 0.80;
            img->lambda_md[5][qp] *= lambda_scale;
            img->lambda_me[5][qp] = sqrt(img->lambda_md[5][qp]);
            img->lambda_mf[5][qp] = LAMBDA_FACTOR (img->lambda_me[5][qp]);
          }
          else
            img->lambda_md[j][qp] *= lambda_scale;

          img->lambda_me[j][qp] = sqrt(img->lambda_md[j][qp]);  
          img->lambda_mf[j][qp] = LAMBDA_FACTOR (img->lambda_me[j][qp]);
        }

#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        qp_temp = (double)qp + img->bitdepth_luma_qp_scale - SHIFT_QP - 6*(img->BitDepthIncrease);

        if (input->UseExplicitLambdaParams) // consideration of explicit weights.
        {
          img->lambda_md[j][qp] = input->LambdaWeight[j] * pow (2, qp_temp/3.0);
          // Scale lambda due to hadamard qpel only consideration
          img->lambda_md[j][qp] = (input->hadamard == 2 ? 0.95 : 1.00) * img->lambda_md[j][qp];
          if (j == B_SLICE)
          {
            img->lambda_md[5][qp] = input->LambdaWeight[5] * pow (2, qp_temp/3.0);
            img->lambda_md[5][qp] = (input->hadamard == 2 ? 0.95 : 1.00) * img->lambda_md[5][qp];
          }
        }
        else
        {
#ifdef ADAPTIVE_QUANTIZATION
          if(input->UseAdaptiveQuantMatrix && img->type == I_SLICE && img->qp==qp)
            img->lambda_md[j][qp] = 0.57 * pow (2.0, qp_temp/3.0); 
          else
#endif
#ifdef RDO_Q
            if(input->UseRDO_Q && img->type == I_SLICE && img->qp==qp)
              img->lambda_md[j][qp] = 0.57 * pow (2.0, qp_temp/3.0); 
            else
#endif
#if defined(USE_INTRA_MDDT)
              if(input->UseIntraMDDT && img->type == I_SLICE && img->qp==qp)
                img->lambda_md[j][qp] = 0.57 * pow (2.0, qp_temp/3.0); 
              else
#endif
                if (input->successive_Bframe>0)
                  img->lambda_md[j][qp] = 0.68 * pow (2, qp_temp/3.0) 
                  * (j == B_SLICE ? Clip3(2.00,4.00,(qp_temp / 6.0)) : (j == SP_SLICE) ? Clip3(1.4,3.0,(qp_temp / 12.0)) : 1.0);
                else
                  img->lambda_md[j][qp] = 0.85 * pow (2, qp_temp/3.0) 
                  * ( (j == B_SLICE) ? 4.0 : (j == SP_SLICE) ? Clip3(1.4,3.0,(qp_temp / 12.0)) : 1.0);
          // Scale lambda due to hadamard qpel only consideration
          img->lambda_md[j][qp] = (input->hadamard == 2 ? 0.95 : 1.00) * img->lambda_md[j][qp];
          img->lambda_md[j][qp] = (j == B_SLICE && input->BRefPictures == 2 && img->b_frame_to_code == 0 ? 0.50 : 1.00) * img->lambda_md[j][qp];

          if (j == B_SLICE)
          {
            img->lambda_md[5][qp] = img->lambda_md[j][qp];

            if (input->HierarchicalCoding == 2)
              img->lambda_md[5][qp] *= (1.0 - min(0.4,0.2 * (double) gop_structure[img->b_frame_to_code-1].hierarchy_layer)) ;
            else
              img->lambda_md[5][qp] *= 0.80;
            img->lambda_md[5][qp] *= lambda_scale;
          }
          else
            img->lambda_md[j][qp] *= lambda_scale;

        }
#endif
      }
    }
  }
  else // RDOPT off computation of Lagrangian multipliers
  {
    for (j = 0; j < 6; j++)
    {
      for (qp = -img->bitdepth_luma_qp_scale; qp < 52; qp++)
      {
        img->lambda_md[j][qp] = img->lambda_me[j][qp] = QP2QUANT[max(0,qp-SHIFT_QP)];
        img->lambda_mf[j][qp] = LAMBDA_FACTOR (img->lambda_me[j][qp]);
      }
    }
  }
}

#ifdef USE_INTRA_MDDT

static const int SCANSTATS8x8[9][64] = 
{
  {85, 67, 46, 41, 32, 25, 15,  9, 44, 31, 26, 21, 15, 10,  5,  2, 27, 23, 18, 16, 10,  6,  3,  1, 18, 16, 15, 11,  8,  3,  2,  0, 13, 10, 10,  8,  6,  2,  1,  0,  7,  7,  5,  4,  3,  2,  1,  0,  4,  3,  3,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0, },
  {79, 32, 15, 10,  7,  4,  2,  1, 56, 23, 13,  9,  6,  3,  2,  0, 38, 19, 13,  8,  4,  3,  1,  0, 34, 18, 11,  7,  4,  2,  1,  0, 36, 15,  9,  6,  4,  2,  1,  0, 26, 10,  6,  5,  3,  2,  1,  0, 17,  6,  3,  2,  2,  1,  0,  0,  6,  1,  1,  0,  0,  0,  0,  0, },
  {88, 38, 21, 12,  8,  4,  2,  1, 53, 21, 15, 10,  6,  3,  2,  1, 26, 15, 12,  9,  4,  2,  1,  0, 20, 11, 10,  6,  5,  2,  1,  0, 15, 11,  8,  6,  4,  2,  1,  0, 12,  7,  6,  5,  3,  1,  1,  0,  7,  3,  2,  3,  1,  1,  0,  0,  2,  1,  0,  1,  1,  0,  0,  0, },
  {88, 67, 41, 27, 17, 13,  5,  3, 65, 45, 33, 25, 18, 11,  6,  1, 42, 35, 31, 24, 13, 12,  4,  1, 31, 24, 27, 21, 14,  9,  5,  2, 21, 24, 18, 15, 12,  7,  5,  2, 15, 11, 10, 10,  7,  7,  2,  1,  9,  6,  4,  4,  2,  2,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0, },
  {85, 69, 42, 29, 20, 12,  6,  2, 66, 53, 44, 29, 18, 10,  5,  1, 43, 46, 41, 29, 16,  7,  5,  1, 30, 30, 32, 26, 20,  9,  4,  1, 25, 21, 22, 23, 19,  9,  4,  1, 15, 14, 12, 15, 10,  7,  2,  1,  8,  6,  5,  5,  5,  3,  1,  0,  1,  1,  1,  0,  1,  0,  0,  0, },
  {89, 69, 42, 30, 22, 13,  6,  1, 63, 46, 38, 29, 22, 10,  7,  1, 38, 34, 28, 22, 15,  9,  7,  1, 25, 25, 22, 12, 11,  8,  4,  1, 21, 16, 15, 12,  7,  5,  2,  1, 13,  8,  9,  7,  3,  4,  1,  0,  4,  3,  3,  2,  2,  1,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0, },
  {86, 57, 34, 20, 15,  7,  2,  0, 69, 42, 29, 19, 13,  6,  3,  0, 49, 40, 30, 19, 10,  5,  2,  0, 36, 33, 27, 18,  9,  4,  3,  1, 24, 23, 22, 13, 10,  4,  3,  0, 16, 14, 13, 12,  5,  3,  1,  0,  7,  5,  4,  4,  5,  2,  0,  0,  1,  1,  0,  1,  0,  0,  0,  0, },
  {88, 70, 51, 37, 26, 15,  7,  2, 61, 46, 40, 31, 21, 12,  7,  2, 39, 33, 27, 24, 17, 10,  4,  1, 27, 21, 23, 17, 13,  8,  4,  1, 19, 16, 16, 12, 11,  5,  3,  1, 11, 12,  9,  8,  5,  4,  2,  0,  5,  4,  4,  4,  3,  2,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0, },
  {87, 57, 31, 20, 13,  7,  3,  1, 70, 44, 30, 18, 13,  6,  3,  1, 47, 40, 28, 19, 11,  5,  3,  1, 34, 32, 23, 16, 10,  5,  2,  1, 26, 26, 20, 14,  9,  3,  2,  1, 17, 15, 12, 10,  6,  2,  1,  0,  8,  8,  6,  4,  3,  1,  0,  0,  2,  1,  1,  1,  0,  0,  0,  0, },
};

static const int SCANSTATS4x4[9][16] = 
{
  {80, 66, 47, 27, 46, 35, 21,  8, 28, 23, 13,  4, 14, 11,  5,  2, },
  {75, 35, 19,  9, 66, 32, 17,  7, 60, 27, 14,  5, 42, 15,  7,  2, },
  {87, 43, 26, 12, 55, 33, 20,  8, 37, 27, 15,  6, 21, 12,  7,  2, },
  {86, 60, 36, 15, 63, 47, 31, 13, 39, 34, 24,  8, 18, 13,  9,  3, },
  {82, 54, 28,  9, 58, 53, 35,  9, 35, 41, 33, 10, 13, 15, 14,  5, },
  {85, 61, 36, 15, 56, 44, 30, 14, 33, 28, 18,  9, 15, 13,  8,  3, },
  {80, 47, 24, 11, 64, 43, 23,  9, 48, 41, 20,  9, 27, 23, 15,  4, },
  {85, 61, 39, 17, 56, 44, 31, 12, 36, 30, 19,  7, 17, 12,  8,  3, },
  {83, 44, 23,  8, 66, 41, 21,  9, 51, 40, 17,  7, 30, 24, 10,  3, },
};

static const int SCANSTATS16x16[4][256] = 
{
  {
    75,  40,  30,  26,  22,  19,  17,  14,  11,   9,   6,   4,   2,   2,   2,   1,  
      26,  16,  11,    8,   7,   6,   6,   4,   3,   2,   1,   1,   0,   0,   0,   0,
      13,   7,   5,   4,   3,   2,   2,   2,   1,   1,   0,   0,   0,   0,   0,   0,
      9,   5,   3,   2,   2,   2,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,
      7,   4,   2,   2,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      6,   3,   2,   1,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      5,   3,   2,   1,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      4,   3,   2,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      4,   3,   2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      3,   2,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      3,   2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      3,   2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  },
  {
    72,  24,  12,   9,   6,   5,   4,   3,   2,   2,   1,   0,   0,   0,   0,   0,
      42,  14,   7,   4,   3,   3,   3,   2,   1,   1,   0,   0,   0,   0,   0,   0,
      31,  11,   4,   2,   2,   1,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,
      25,   8,   3,   2,   1,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,
      21,   7,   3,   2,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      18,   6,   3,   1,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      16,   6,   3,   2,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      15,   6,   3,   2,   1,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,
      13,   6,   3,   2,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      11,   4,   2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      10,   4,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      8,   3,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      6,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      5,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      4,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      3,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    },
    {
      79,  42,  25,  18,  13,  11,   9,   7,   5,   4,   2,   1,   1,   0,   0,   0,
        43,  25,  16,  11,   9,   8,   7,   5,   4,   3,   1,   1,   0,   0,   0,   0,
        27,  16,   9,   6,   5,   4,   4,   3,   2,   1,   1,   0,   0,   0,   0,   0,
        19,  11,   6,   4,   3,   3,   2,   2,   1,   1,   0,   0,   0,   0,   0,   0,
        15,   9,   5,   4,   3,   2,   2,   1,   1,   0,   0,   0,   0,   0,   0,   0,
        13,   8,   5,   3,   2,   2,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,
        12,   8,   5,   3,   2,   2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,
        10,   8,   4,   3,   2,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,
        10,   7,   4,   3,   2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        8,   6,   3,   2,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        7,   4,   3,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        6,   4,   2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        4,   3,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        3,   2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        3,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    },
    {
      81,  42,  22,  14,   8,   6,   4,   3,   2,   2,   1,   0,   0,   0,   0,   0,
        43,  20,  12,   7,   5,   4,   3,   2,   2,   1,   0,   0,   0,   0,   0,   0,
        25,  12,   6,   4,   3,   2,   2,   1,   1,   0,   0,   0,   0,   0,   0,   0,
        15,   8,   4,   3,   2,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,
        11,   6,   3,   2,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        8,   5,   3,   2,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        7,   4,   2,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        6,   4,   2,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        6,   4,   2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        4,   3,   2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        4,   2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        3,   2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        2,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      },
};

/*! 
*************************************************************************************
* \brief
*   calculate the scanning order based on the input scanning stats
*
* \para calcScanOrder()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
static int  calcScanOrder(int *stats, int *orderX, int *orderY, int size, int width)
{
  int i, j, cOrder, cStats;
  int *stats1D = stats; 
  int order1D[256];
  int orderChanged = 0;

  for(i = 0; i < size; i++)
  {
    order1D[i] = orderX[i]+orderY[i]*width;
  }

  for (i=1; i < size; i++)
  {
    cStats = stats1D[i];
    cOrder = order1D[i];
    j = i;
    while ((j > 0) && (stats1D[j-1] < cStats) )
    {
      stats1D[j] = stats1D[j-1];
      order1D[j] = order1D[j-1];
      j = j - 1;
      orderChanged = 1;
    }
    stats1D[j] = cStats;
    order1D[j] = cOrder;
  }

  for(i = 0; i < size; i++)
  {
    orderX[i] = order1D[i]%width;
    orderY[i] = order1D[i]/width;
  }

  return orderChanged;
}

/*! 
*************************************************************************************
* \brief
*   At the beginning of each slice, initialize the scanning orders for all intra 
*   coding modes
*
* \para InitScanOrderForSlice()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
static void InitScanOrderForSlice()
{
  int ipredmode;
  int stats4x4[9][16];
  int stats8x8[9][64];
  int orderX4x4[16], orderY4x4[16];
  int orderX8x8[64], orderY8x8[64];
  int stats16x16[9][256];
  int orderX16x16[256], orderY16x16[256];

  for(ipredmode = 0; ipredmode < 9; ipredmode++)
  {
    int k;

    img->update4x4[ipredmode] = img->update8x8[ipredmode] = 1;
    for(k = 0; k < 16; k++)
    {
      img->scanStats4x4[ipredmode][k] = stats4x4[ipredmode][k] = SCANSTATS4x4[ipredmode][k]/2;
    }

    for(k = 0; k < 64; k++)
    {
      img->scanStats8x8[ipredmode][k] = stats8x8[ipredmode][k] = SCANSTATS8x8[ipredmode][k]/2;
    }
    img->update4x4Count[ipredmode] = img->update8x8Count[ipredmode] = 0;
    img->update4x4Thres[ipredmode] = 4;
    img->update8x8Thres[ipredmode] = 2;
  }

  for(ipredmode = 0; ipredmode < 9; ipredmode ++)
  {
    int i,dummy;

    {
      for(i = 0; i < 16; i++)
      {
        orderX4x4[i] = i%4;
        orderY4x4[i] = i/4;
      }
      dummy = calcScanOrder(stats4x4[ipredmode], orderX4x4, orderY4x4, 16, 4);
      for(i = 0; i < 16; i++)
      {
        img->scanOrder4x4[ipredmode][i][0] = orderX4x4[i];
        img->scanOrder4x4[ipredmode][i][1] = orderY4x4[i];
        img->scanStats4x4[ipredmode][i]    = stats4x4[ipredmode][i];
      }
    }

    {
      for(i = 0; i < 64; i++)
      {
        orderX8x8[i] = i%8;
        orderY8x8[i] = i/8;
      }
      dummy = calcScanOrder(stats8x8[ipredmode], orderX8x8, orderY8x8, 64, 8);    
      for(i = 0; i < 64; i++)
      {
        img->scanOrder8x8[ipredmode][i][0] = orderX8x8[i];
        img->scanOrder8x8[ipredmode][i][1] = orderY8x8[i];
        img->scanStats8x8[ipredmode][i]    = stats8x8[ipredmode][i];
      }
    }
  }

  for(ipredmode = 0; ipredmode < 4; ipredmode++)
  {
    int k, dummy;

    for(k = 0; k < 256; k++)
    {
      stats16x16[ipredmode][k] = SCANSTATS16x16[ipredmode][k]/2;
    }

    for(k = 0; k < 256; k++)
    {
      orderX16x16[k] = k%16;
      orderY16x16[k] = k/16;
    }
    dummy = calcScanOrder(stats16x16[ipredmode], orderX16x16, orderY16x16, 256, 16);
    for(k = 0;  k < 256; k++)
    {
      img->scanOrder16x16[ipredmode][k][0] = orderX16x16[k];
      img->scanOrder16x16[ipredmode][k][1] = orderY16x16[k];
    }
  }
}


/*! 
*************************************************************************************
* \brief
*   Update the scanning orders for all intra coding modes based on their coeff stats
*
* \para updateScanOrder()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
void updateScanOrder(int first)
{
  int ipredmode;
  int orderX[64], orderY[64];

  int orderChanged;

  for(ipredmode = 0; ipredmode < 9; ipredmode ++)
  {
    int i;

    if(img->update4x4[ipredmode])
    {
      for(i = 0; i < 16; i++)
      {
        orderX[i] = img->scanOrder4x4[ipredmode][i][0];
        orderY[i] = img->scanOrder4x4[ipredmode][i][1];
      }

      orderChanged = calcScanOrder(img->scanStats4x4[ipredmode], orderX, orderY, 16, 4);
      if(!orderChanged)
      {
        img->update4x4Thres[ipredmode] <<= 1;
      }
      else if(img->update4x4Thres[ipredmode] > 4) 
        img->update4x4Thres[ipredmode] >>= 1;
      if(orderChanged)
        for(i = 0; i < 16; i++)
        {
          img->scanOrder4x4[ipredmode][i][0] = orderX[i];
          img->scanOrder4x4[ipredmode][i][1] = orderY[i];
        }
    }
    if(img->update8x8[ipredmode])
    {
      for(i = 0; i < 64; i++)
      {
        orderX[i] = img->scanOrder8x8[ipredmode][i][0];
        orderY[i] = img->scanOrder8x8[ipredmode][i][1];
      }
      orderChanged = calcScanOrder(img->scanStats8x8[ipredmode], orderX, orderY, 64, 8);  
      if(!orderChanged)
      {
        img->update8x8Thres[ipredmode] <<= 1;
      }
      else if(img->update8x8Thres[ipredmode] > 2) 
        img->update8x8Thres[ipredmode] >>= 1;
      if(orderChanged)
        for(i = 0; i < 64; i++)
        {
          img->scanOrder8x8[ipredmode][i][0] = orderX[i];
          img->scanOrder8x8[ipredmode][i][1] = orderY[i];
        }
    }
  }
}

/*! 
*************************************************************************************
* \brief
*   Scale down the scanning stats during re-normalization
*
* \para scaleScanStats()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
static void scaleScanStats(int *stats, int size)
{
  int i;

  for(i = 0; i < size*size; i++)
    stats[i] >>= 1;
}

/*! 
*************************************************************************************
* \brief
*   Re-normalize the scanning stats as necessary
*
* \para normalizeScanStats()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
void normalizeScanStats()
{
  int ipredmode;
  int mb_type = img->mb_data[img->current_mb_nr].mb_type;

  if(mb_type == I4MB)
    for(ipredmode = 0; ipredmode < 9; ipredmode++)
    {
      if(img->scanStats4x4[ipredmode][0] >= 256)
      {
        scaleScanStats(img->scanStats4x4[ipredmode], 4);
      }
    }

  else if(mb_type == I8MB)
    for(ipredmode = 0; ipredmode < 9; ipredmode++)
    {
      if(img->scanStats8x8[ipredmode][0] >= 256)
      {
        scaleScanStats(img->scanStats8x8[ipredmode], 8);
      }
    }
}

#endif 
#if defined(ADAPTIVE_QUANTIZATION) || defined(RDO_Q)
/*!
************************************************************************
* \brief
*    encode one macroblock for RDO_Q or AQMS
*
************************************************************************
*/
#ifdef MB32X32
double encode_one_macroblock_each_quant(int CurrentMbAddr)
#else
void encode_one_macroblock_each_quant(int CurrentMbAddr)
#endif
{
#ifdef MB32X32
  int saved_cluster_qp, saved_prev_qp, saved_prev_delta_qp;
#endif
  CalculateQuantParam();
  CalculateOffsetParam();

  if(input->Transform8x8Mode)
  {
    CalculateQuant8Param();
    CalculateOffset8Param();
  }

#ifdef MB32X32

  if(input->UseExtMB!=0)
  {
    memset(saveMode.MBsize16, -1 , 16 * sizeof(int));
    memset(saveMode.MBsize32, -1 ,  4 * sizeof(int));
    saveMode.MBsize64 = -1;


    if(input->UseExtMB==0 || img->type == I_SLICE)
    {
      start_macroblock (CurrentMbAddr, FALSE);
      encode_one_macroblock ();
  //#if SPEEDUP
  //    saveMode.MBsize16[0] = img->mb_data[CurrentMbAddr].mb_type;
  //#endif
      return rdopt->min_rdcost;
    }
    else
    {
      double rd_cost;
  //#ifdef MB32_DELTA_QP
      //qp handling
      start_macroblock(CurrentMbAddr, FALSE); // the purpose of calling this function is to set img->mb_data[CurrentMbAddr].(prev_qp, prev_delta_qp, qp)
      saved_prev_qp = img->mb_data[CurrentMbAddr].prev_qp;
      saved_prev_delta_qp = img->mb_data[CurrentMbAddr].prev_delta_qp;
      saved_cluster_qp = img->mb_data[CurrentMbAddr].qp;
      delta_qp_sent = 0;
  //#endif
      if(input->UseExtMB == 2)
        rd_cost = encode_macroblock_cluster64(CurrentMbAddr, saved_prev_qp, saved_prev_delta_qp, saved_cluster_qp, &delta_qp_sent);//MB64X64
      else
        rd_cost = encode_macroblock_cluster32(CurrentMbAddr, 0, saved_prev_qp, saved_prev_delta_qp, saved_cluster_qp, &delta_qp_sent);

      return rd_cost;

    }
  }
  else
  {
#endif
  start_macroblock (CurrentMbAddr, FALSE);
  encode_one_macroblock ();
#ifdef MB32X32
  return 0;
  }
#endif
}
#endif
#ifdef ADAPTIVE_QUANTIZATION
/*! 
*************************************************************************************
* \brief
*    initialize skip stat
*************************************************************************************
*/
void InitAQMSStatParamForSlice()
{
  int i=0;
  for(i=0;i<6;i++)
  {
    dwSkipMbCountForSlice[i]=0;
    dwSkipEnableMbCountForSlice[i]=0;
  }
}
/*!
************************************************************************
* \brief
*    encode one macroblock for AQMS
*
************************************************************************
*/
void encode_one_macroblock_for_iaqms(int CurrentMbAddr, double *min_rdcost, int *best_iaqms_idx, int *saved_best_mode, int *best_QP, int *best_deltaQP, int *best_athr, int *best_cbp, int masterQP, int deltaQP)
{
  int iaqms_idx;
  static int saqms_idx_start[6] = { 0,  0,  0,  0,  0,  0}; 
  static int saqms_idx_last[6]  = { NUM_OF_IAQMS_MAT,  NUM_OF_IAQMS_MAT,  NUM_OF_IAQMS_MAT,  NUM_OF_IAQMS_MAT,  NUM_OF_IAQMS_MAT,  NUM_OF_IAQMS_MAT};

  start_macroblock (CurrentMbAddr, FALSE);

  if(img->type==B_SLICE)
    saqms_idx_last[B_SLICE] = 2;
  if(input->ProfileIDC < FREXT_HP && img->type == P_SLICE)    
    saqms_idx_last[P_SLICE] = 2;

  // modulated quantization matrix selection
  for (iaqms_idx=saqms_idx_start[img->type]; iaqms_idx<saqms_idx_last[img->type]; iaqms_idx++)
  {
    img->mb_iaqms_idx = iaqms_idx;
    img->qp=Clip3(-img->bitdepth_luma_qp_scale, 51, (masterQP+deltaQP));

    // This code is fast encoding techniques but efficiency is a little down for 720p60 sequences.
#ifdef RDO_Q
    if(input->UseRDO_Q && img->type==I_SLICE && img->qp > (int)gdwBaseQp && img->PicSizeInMbs > RES_SKIP_SIZE) continue;
#endif
    if(*best_cbp==0) continue;

    encode_one_macroblock_each_quant(CurrentMbAddr);

    Motion_Selected = 1;
    if (rdopt->min_rdcost<*min_rdcost)
    {
      *min_rdcost=rdopt->min_rdcost;
      *best_iaqms_idx = img->mb_data[img->current_mb_nr].mb_iaqms_idx;
      *saved_best_mode = img->mb_data[img->current_mb_nr].mb_type;//best_mode;
      *best_deltaQP = deltaQP;
      *best_QP = img->qp;
      *best_cbp = img->mb_data[img->current_mb_nr].cbp;
    }
  }
  img->mb_iaqms_idx = *best_iaqms_idx;
  decide_last_coeff(CurrentMbAddr, min_rdcost, best_iaqms_idx, saved_best_mode, best_QP, best_deltaQP, best_athr, best_cbp, masterQP, deltaQP);
}

/*!
************************************************************************
* \brief
*    store internal parameter
*
************************************************************************
*/
void init_internal_parameter(int CurrentMbAddr, int *best_iaqms_idx, int *saved_best_mode, int *best_QP, int *best_deltaQP, int *best_athr, int *best_cbp)
{
  // initialization
  img->mb_iaqms_idx = 0;
  *best_iaqms_idx = 0;
  *saved_best_mode = 0;
  *best_QP = img->qp;
  *best_deltaQP = 0;
  *best_athr = 0;
  *best_cbp = 47;

  internal_param[0] = input->disthres;
}

/*!
************************************************************************
* \brief
*    set internal parameter
*
************************************************************************
*/
void set_internal_parameter()
{
  input->disthres = internal_param[0];
}
/*!
************************************************************************
* \brief
*    encode one macroblock for AQMS
*
************************************************************************
*/
void decide_last_coeff(int CurrentMbAddr, double *min_rdcost, int *best_iaqms_idx, int *saved_best_mode, int *best_QP, int *best_deltaQP, int *best_athr, int *best_cbp, int masterQP, int deltaQP)
{

  static int sathr_length = 3;
  int iathr_idx;

#ifdef RDO_Q
  if(input->UseRDO_Q && img->PicSizeInMbs > RES_SKIP_SIZE) sathr_length = 2;
#endif

  // adaptive thresholding selection
  for (iathr_idx=1; iathr_idx<sathr_length; iathr_idx++)
  {
    input->disthres = iathr_idx;
    img->qp=Clip3(-img->bitdepth_luma_qp_scale, 51, (masterQP+deltaQP));

#ifdef RDO_Q
    if(input->UseRDO_Q && img->type==I_SLICE && img->qp > (int)gdwBaseQp && img->PicSizeInMbs > RES_SKIP_SIZE) continue;
#endif

    encode_one_macroblock_each_quant(CurrentMbAddr);

    Motion_Selected = 1;
    if (rdopt->min_rdcost<*min_rdcost)
    {
      *min_rdcost=rdopt->min_rdcost;
      *best_athr = iathr_idx;
      *saved_best_mode = img->mb_data[img->current_mb_nr].mb_type;
      *best_deltaQP = deltaQP;
      *best_QP = img->qp;
      *best_cbp = img->mb_data[img->current_mb_nr].cbp;
    }
  }
  input->disthres = *best_athr;
}

#endif
