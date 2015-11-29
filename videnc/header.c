
/*!
*************************************************************************************
* \file header.c
*
* \brief
*    H.264 Slice and Sequence headers
*
* \author
*    Main contributors (see contributors.h for copyright, address and affiliation details)
*      - Stephan Wenger                  <stewe@cs.tu-berlin.de>
*      - Karsten Suehring                <suehring@hhi.de>
*************************************************************************************
*/

#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "global.h"

#include "elements.h"
#include "header.h"
#include "rtp.h"
#include "mbuffer.h"
#include "defines.h"
#include "vlc.h"
#include "parset.h"
#ifdef ADAPTIVE_FILTER
#include "adaptive_filter.h"
#endif

#ifdef SWITCHED_FILTERS
#include "switched_filters.h"
#endif  // SWITCHED_FILTERS

// A little trick to avoid those horrible #if TRACE all over the source code
#if TRACE
#define SYMTRACESTRING(s) strncpy(sym->tracestring,s,TRACESTRING_SIZE)
#else
#define SYMTRACESTRING(s) // do nothing
#endif

int * assignSE2partition[2] ;
int assignSE2partition_NoDP[SE_MAX_ELEMENTS] =
{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
int assignSE2partition_DP[SE_MAX_ELEMENTS] =
{  0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 2, 2, 2, 2, 0, 0, 0, 0 } ;

static int ref_pic_list_reordering(Bitstream *bitstream);
static int dec_ref_pic_marking    (Bitstream *bitstream);
static int pred_weight_table      (Bitstream *bitstream);

/*!
********************************************************************************************
* \brief 
*    Write a slice header
*
* \return
*    number of bits used 
********************************************************************************************
*/
int SliceHeader()
{
  int dP_nr = assignSE2partition[input->partition_mode][SE_HEADER];
  Bitstream *bitstream = img->currentSlice->partArr[dP_nr].bitstream;
  Slice* currSlice = img->currentSlice;
  int len = 0;
  unsigned int field_pic_flag = 0, bottom_field_flag = 0;

  int num_bits_slice_group_change_cycle;
  float numtmp;  
#ifdef ADAPTIVE_FILTER
  int i, j, apply_adaptive_filter;
  extern int UseAllSubpelPositions;                  // 1 if FilterFlag for all independent positions is 1
  extern int SubpelPositionsPattern;
#ifndef EDAIF2
  extern int DiffQFilterCoef[15][SQR_FILTER];        // differences to be transmitted
  extern int POS_EQUATION_NUMBER[15];                // number of different coefficietns for each sub-pel position
  extern int FilterFlag[15];                        // Flags defining if Filter at the particular position calculated
  extern int FILTER_SIZE;
  extern int FILTER_OFFSET;
#else
  extern int DiffQFilterCoef[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];        // differences to be transmitted
  extern int POS_EQUATION_NUMBER[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];                // number of different coefficietns for each sub-pel position
  extern int FilterFlag[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];                        // Flags defining if Filter at the particular position calculated
#endif
#endif

  if (img->MbaffFrameFlag)
    len  = ue_v("SH: first_mb_in_slice", img->current_mb_nr >> 1,   bitstream);
  else
    len  = ue_v("SH: first_mb_in_slice", img->current_mb_nr,   bitstream);

  len += ue_v("SH: slice_type",        get_picture_type (),   bitstream);

  len += ue_v("SH: pic_parameter_set_id" , active_pps->pic_parameter_set_id ,bitstream);

  len += u_v (log2_max_frame_num_minus4 + 4,"SH: frame_num", img->frame_num, bitstream);

  if (!active_sps->frame_mbs_only_flag)
  {
    // field_pic_flag    u(1)
    field_pic_flag = (img->structure ==TOP_FIELD || img->structure ==BOTTOM_FIELD)?1:0;
    assert( field_pic_flag == img->fld_flag );
    len += u_1("SH: field_pic_flag", field_pic_flag, bitstream);

    if (field_pic_flag)
    {
      //bottom_field_flag     u(1)
      bottom_field_flag = (img->structure == BOTTOM_FIELD)?1:0;
      len += u_1("SH: bottom_field_flag" , bottom_field_flag ,bitstream);
    }
  }

  if (img->currentPicture->idr_flag)
  {
    // idr_pic_id
    len += ue_v ("SH: idr_pic_id", (img->number % 2), bitstream);
  }

  if (img->pic_order_cnt_type == 0)
  {
    if (active_sps->frame_mbs_only_flag)
    {
      img->pic_order_cnt_lsb = (img->toppoc & ~((((unsigned int)(-1)) << (log2_max_pic_order_cnt_lsb_minus4+4))) );
    }
    else
    {
      if (!field_pic_flag || img->structure == TOP_FIELD)
        img->pic_order_cnt_lsb = (img->toppoc & ~((((unsigned int)(-1)) << (log2_max_pic_order_cnt_lsb_minus4+4))) );
      else if ( img->structure == BOTTOM_FIELD )
        img->pic_order_cnt_lsb = (img->bottompoc & ~((((unsigned int)(-1)) << (log2_max_pic_order_cnt_lsb_minus4+4))) );
    }

    len += u_v (log2_max_pic_order_cnt_lsb_minus4+4, "SH: pic_order_cnt_lsb", img->pic_order_cnt_lsb, bitstream);

    if (img->pic_order_present_flag && !field_pic_flag)
    {
      len += se_v ("SH: delta_pic_order_cnt_bottom", img->delta_pic_order_cnt_bottom, bitstream);
    }
  }
  if (img->pic_order_cnt_type == 1 && !img->delta_pic_order_always_zero_flag)
  {
    len += se_v ("SH: delta_pic_order_cnt[0]", img->delta_pic_order_cnt[0], bitstream);

    if (img->pic_order_present_flag && !field_pic_flag)
    {
      len += se_v ("SH: delta_pic_order_cnt[1]", img->delta_pic_order_cnt[1], bitstream);
    }
  }

  if (active_pps->redundant_pic_cnt_present_flag)
  {
    len += ue_v ("SH: redundant_pic_cnt", img->redundant_pic_cnt, bitstream);
  }

  // Direct Mode Type selection for B pictures
  if (img->type==B_SLICE)
  {
    len +=  u_1 ("SH: direct_spatial_mv_pred_flag", img->direct_spatial_mv_pred_flag, bitstream);    
  }

  if ((img->type == P_SLICE) || (img->type == B_SLICE) || (img->type==SP_SLICE))
  {
    int override_flag;
    if ((img->type == P_SLICE) || (img->type==SP_SLICE))
    {
      override_flag = (img->num_ref_idx_l0_active != (active_pps->num_ref_idx_l0_active_minus1 +1)) ? 1 : 0;
    }
    else
    {
      override_flag = ((img->num_ref_idx_l0_active != (active_pps->num_ref_idx_l0_active_minus1 +1)) 
        || (img->num_ref_idx_l1_active != (active_pps->num_ref_idx_l1_active_minus1 +1))) ? 1 : 0;
    }

    len +=  u_1 ("SH: num_ref_idx_active_override_flag", override_flag, bitstream);

    if (override_flag) 
    {
      len += ue_v ("SH: num_ref_idx_l0_active_minus1", img->num_ref_idx_l0_active-1, bitstream);
      if (img->type==B_SLICE)
      {
        len += ue_v ("SH: num_ref_idx_l1_active_minus1", img->num_ref_idx_l1_active-1, bitstream);
      }
    }

  }
  len += ref_pic_list_reordering(bitstream);

  if (((img->type == P_SLICE || img->type == SP_SLICE) && active_pps->weighted_pred_flag) || 
    ((img->type == B_SLICE) && active_pps->weighted_bipred_idc == 1))  
  {
    len += pred_weight_table(bitstream);
  }

  if (img->nal_reference_idc)
    len += dec_ref_pic_marking(bitstream);

  if(input->symbol_mode==CABAC && img->type!=I_SLICE /*&& img->type!=SI_IMG*/)
  {
    len += ue_v("SH: cabac_init_idc", img->model_number, bitstream);
  }

  len += se_v("SH: slice_qp_delta", (currSlice->qp - 26 - active_pps->pic_init_qp_minus26), bitstream);  

  if (img->type==SP_SLICE /*|| img->type==SI_SLICE*/)
  {
    if (img->type==SP_SLICE) // Switch Flag only for SP pictures
    {
      len += u_1 ("SH: sp_for_switch_flag", (si_frame_indicator || sp2_frame_indicator), bitstream);   // 1 for switching SP, 0 for normal SP
    }
    len += se_v ("SH: slice_qs_delta", (img->qpsp - 26), bitstream );
  }

  if (active_pps->deblocking_filter_control_present_flag)
  {
    len += ue_v("SH: disable_deblocking_filter_idc",img->LFDisableIdc, bitstream);  // Turn loop filter on/off on slice basis 

    if (img->LFDisableIdc!=1)
    {
      len += se_v ("SH: slice_alpha_c0_offset_div2", img->LFAlphaC0Offset / 2, bitstream);

      len += se_v ("SH: slice_beta_offset_div2", img->LFBetaOffset / 2, bitstream);
    }
  }
#ifdef ADAPTIVE_QUANTIZATION
  len += u_1 ("SH: slice_fractional_quant_flag", (img->slice_fractional_quant_flag&1), bitstream);
  if(img->slice_fractional_quant_flag)
  {
    len += u_1 ("SH: slice_mqm_signaling_flag", (img->slice_mqm_signaling_flag&1), bitstream);
    if(img->slice_mqm_signaling_flag)
    {
      len+=se_v  ("SH: slice_scaling_model_param0",   (img->slice_modeling_qm_param0-32),    bitstream);
      len+=se_v  ("SH: slice_scaling_model_param1",   (img->slice_modeling_qm_param1-16),    bitstream);
    }
  }
#endif
#ifdef EIGHTH_PEL
  len += ue_v ("SH: motion vector resolution", input->mv_res==1?1:0, bitstream);
#endif
#ifdef ADAPTIVE_FD_SD_CODING
  len += u_1 ("SH: SD coding", img->APEC_in_FD_and_SD==1?1:0, bitstream);
  if (img->APEC_in_FD_and_SD)
  {
    len += u_1 ("SH: Quantizer"                , input->SD_Quantizer==1?1:0, bitstream);
  }
#endif
  if ( active_pps->num_slice_groups_minus1>0 &&
    active_pps->slice_group_map_type>=3 && active_pps->slice_group_map_type<=5)
  {
    numtmp=img->PicHeightInMapUnits*img->PicWidthInMbs/(float)(active_pps->slice_group_change_rate_minus1+1)+1;
    num_bits_slice_group_change_cycle = (int)ceil(log(numtmp)/log(2));

    //! img->slice_group_change_cycle can be changed before calling FmoInit()
    len += u_v (num_bits_slice_group_change_cycle, "SH: slice_group_change_cycle", img->slice_group_change_cycle, bitstream);
  }

  // NOTE: The following syntax element is actually part 
  //        Slice data bitstream A RBSP syntax

  if(input->partition_mode&&!img->currentPicture->idr_flag)
  {
    len += ue_v("DPA: slice_id", img->current_slice_nr, bitstream);
  }
#ifdef ADAPTIVE_FILTER
  if(input->UseAdaptiveFilter == 1 && img->AdaptiveFilterFlag == 1)
  {
#ifndef E_DAIF
    len += u_v(2, "SH: apply_adaptive_filter", input->UseAdaptiveFilter, bitstream);
#else
    len += u_v(3, "SH: apply_adaptive_filter", input->UseAdaptiveFilter, bitstream);
#endif

    len += u_1("SH: use_all_subpel_positions", UseAllSubpelPositions, bitstream);
    if(!UseAllSubpelPositions)
      len += ue_v("SH: positions_pattern", SubpelPositionsPattern, bitstream);

    if(FilterFlag[a_pos])
      for(i= 0; i < POS_EQUATION_NUMBER[a_pos]; i++)
        len += se_v("SH: a_pos", DiffQFilterCoef[a_pos][FILTER_SIZE*FILTER_OFFSET+i], bitstream);
    if(FilterFlag[b_pos])
      for(i= 0; i < POS_EQUATION_NUMBER[b_pos]; i++)
        len += se_v("SH: b_pos", DiffQFilterCoef[b_pos][FILTER_SIZE*FILTER_OFFSET+i], bitstream);
    if(FilterFlag[e_pos])
      for(i = 0; i < FILTER_SIZE; i++)
        for(j = i; j < FILTER_SIZE; j++)
          len += se_v("SH: e_pos", DiffQFilterCoef[e_pos][FILTER_SIZE*i+j], bitstream);
    if(FilterFlag[f_pos])
      for(i = 0; i < FILTER_SIZE; i++)
        for(j = 0; j < FILTER_SIZE/2; j++)
          len += se_v("SH: f_pos", DiffQFilterCoef[f_pos][FILTER_SIZE*i+j], bitstream);
    if(FilterFlag[j_pos])
      for(i = 0; i < FILTER_SIZE/2; i++)
        for(j = i; j < FILTER_SIZE/2; j++)
          len += se_v("SH: j_pos", DiffQFilterCoef[j_pos][FILTER_SIZE*i+j], bitstream);
  }
  else if(input->UseAdaptiveFilter == 2 && img->AdaptiveFilterFlag == 1) //separable aif
  {
#ifndef E_DAIF
    len += u_v(2, "SH: apply_adaptive_filter", input->UseAdaptiveFilter, bitstream);
#else
    len += u_v(3, "SH: apply_adaptive_filter", input->UseAdaptiveFilter, bitstream);
#endif

    len += u_1("SH: use_all_subpel_positions", UseAllSubpelPositions, bitstream);
    if(!UseAllSubpelPositions)
      len += ue_v("SH: positions_pattern", SubpelPositionsPattern, bitstream);

    // a_pos
    if(FilterFlag[a_pos])
      for(i= 0; i < FILTER_SIZE; i++)
        len += se_v("SH: a_pos", DiffQFilterCoef[a_pos][i], bitstream);

    // b_pos        
    if(FilterFlag[b_pos])
      for(i= 0; i < FILTER_SIZE>>1; i++)
        len += se_v("SH: b_pos", DiffQFilterCoef[b_pos][i], bitstream);

    // c_pos             
    if(FilterFlag[c_pos])
    {
      if (FilterFlag[a_pos])
        for(i= 0; i < FILTER_SIZE; i++)
          len += se_v("SH: c_pos", DiffQFilterCoef[c_pos][i] - DiffQFilterCoef[a_pos][FILTER_SIZE-i-1], bitstream);
      else 
        for(i= 0; i < FILTER_SIZE; i++)
          len += se_v("SH: c_pos", DiffQFilterCoef[c_pos][i], bitstream);
    }

    // d_pos             
    if(FilterFlag[d_pos])
      for(i= 0; i < FILTER_SIZE; i++)
        len += se_v("SH: d_pos", DiffQFilterCoef[d_pos][i], bitstream);

    // e_pos             
    if(FilterFlag[e_pos])
    {
      if (FilterFlag[d_pos])
        for(i= 0; i < FILTER_SIZE; i++)
          len += se_v("SH: e_pos", DiffQFilterCoef[e_pos][i] - DiffQFilterCoef[d_pos][i], bitstream);
      else
        for(i= 0; i < FILTER_SIZE; i++)
          len += se_v("SH: e_pos", DiffQFilterCoef[e_pos][i], bitstream);
    }

    // f_pos            
    if(FilterFlag[f_pos])
    {
      if (FilterFlag[e_pos])
        for(i= 0; i < FILTER_SIZE; i++)
          len += se_v("SH: f_pos", DiffQFilterCoef[f_pos][i] - DiffQFilterCoef[e_pos][i], bitstream);
      else
        for(i= 0; i < FILTER_SIZE; i++)
          len += se_v("SH: f_pos", DiffQFilterCoef[f_pos][i], bitstream);
    }

    // g_pos             
    if(FilterFlag[g_pos])
    {
      if (FilterFlag[f_pos])
        for(i= 0; i < FILTER_SIZE; i++)
          len += se_v("SH: g_pos", DiffQFilterCoef[g_pos][i] - DiffQFilterCoef[f_pos][i], bitstream);
      else
        for(i= 0; i < FILTER_SIZE; i++)
          len += se_v("SH: g_pos", DiffQFilterCoef[g_pos][i], bitstream);
    }

    // h_pos             
    if(FilterFlag[h_pos])
      for(i= 0; i < FILTER_SIZE>>1; i++)
        len += se_v("SH: h_pos", DiffQFilterCoef[h_pos][i], bitstream);

    // i_pos     
    if(FilterFlag[i_pos])
    {
      if (FilterFlag[h_pos])
        for(i= 0; i < FILTER_SIZE>>1; i++)
          len += se_v("SH: i_pos", DiffQFilterCoef[i_pos][i] - DiffQFilterCoef[h_pos][i], bitstream);
      else 
        for(i= 0; i < FILTER_SIZE>>1; i++)
          len += se_v("SH: i_pos", DiffQFilterCoef[i_pos][i], bitstream);
    }

    // j_pos            
    if(FilterFlag[j_pos])
    {
      if (FilterFlag[i_pos])
        for(i= 0; i < FILTER_SIZE>>1; i++)
          len += se_v("SH: j_pos", DiffQFilterCoef[j_pos][i] - DiffQFilterCoef[i_pos][i], bitstream);
      else
        for(i= 0; i < FILTER_SIZE>>1; i++)
          len += se_v("SH: j_pos", DiffQFilterCoef[j_pos][i], bitstream);
    }

    // k_pos             
    if(FilterFlag[k_pos])
    {
      if (FilterFlag[j_pos])
        for(i= 0; i < FILTER_SIZE>>1; i++)
          len += se_v("SH: k_pos", DiffQFilterCoef[k_pos][i] - DiffQFilterCoef[j_pos][i], bitstream);
      else 
        for(i= 0; i < FILTER_SIZE>>1; i++)
          len += se_v("SH: k_pos", DiffQFilterCoef[k_pos][i], bitstream);
    }
  }
#ifdef DIRECTIONAL_FILTER
  else if((input->UseAdaptiveFilter == 3) && img->AdaptiveFilterFlag == 1) //1D-AIF
  {
#ifndef E_DAIF
    len += u_v(2, "SH: apply_adaptive_filter", input->UseAdaptiveFilter, bitstream);
#else
    len += u_v(3, "SH: apply_adaptive_filter", input->UseAdaptiveFilter, bitstream);
#endif
    len += u_v(1, "SH: implementation_type", input->ImpType, bitstream);
    len += sendCoefsAIF(bitstream);
  }
#endif
#ifdef E_DAIF
  else if((input->UseAdaptiveFilter == FILTER_TYPE_EDAIF) && img->AdaptiveFilterFlag == 1) //1D-AIF
  {
    len += u_v(3, "SH: apply_adaptive_filter", input->UseAdaptiveFilter, bitstream);
    len += sendCoefs_EDAIF(bitstream);
  }
#endif 
#ifdef EAIF
  else if((input->UseAdaptiveFilter == FILTER_TYPE_EAIF) && img->AdaptiveFilterFlag == 1) //1D-AIF
  {
    len += u_v(3, "SH: apply_adaptive_filter", input->UseAdaptiveFilter, bitstream);
    len += sendCoefs_EAIF(bitstream);
  }
#endif
  else if(input->UseAdaptiveFilter && (img->type == P_SLICE || img->type == B_SLICE))
  {
    apply_adaptive_filter = 0;
#ifndef E_DAIF
    len += u_v(2, "SH: apply_adaptive_filter", apply_adaptive_filter, bitstream);
#else
    len += u_v(3, "SH: apply_adaptive_filter", apply_adaptive_filter, bitstream);
#endif
  }  
  //  printf("aif bits = %d\n",len - aif_temp);
#endif
#ifdef USE_INTRA_MDDT
  len += u_v(1, "SE: use_intra_MDDT", input->UseIntraMDDT, bitstream);
#endif

#ifdef USE_HP_FILTER  
#ifdef SWITCHED_FILTERS
  len += u_v(3, "SE: use_HP_filter", input->UseHPFilter, bitstream);
#else
  len += u_v(2, "SE: use_HP_filter", input->UseHPFilter, bitstream);
#endif
#endif

#ifdef SWITCHED_FILTERS
if(img->type != I_SLICE)
{	
  if((input->UseHPFilter == HPF_SIFO) || (input->UseHPFilter == HPF_SIFO_FPO))
  {
    // Save sequence or frame filters 
    if((img->filterParam == SIFO_FIRST_PASS) || 
      ((img->type == B_SLICE) && (img->filterParam == SIFO_SEQ_FILTER)))
    {     
      int sub_pos;
      len += u_v(1, "SH: use offsets", SIFO_FILTERS_ONLY, bitstream);  // 1 bit only
      len += u_v (2, "SH: bestFilter", NO_BEST_FILTER, bitstream);
      for(sub_pos = 1; sub_pos < 16; ++sub_pos)
      {
        if(img->filterSequence[sub_pos] == 0)
          len += u_1("SH: filterSequence", 0, bitstream);
        else 
        {
          len += u_1("SH: filterSequence", 1, bitstream);
          if(img->filterSequence[sub_pos] == 1)
            len += u_1("SH: filterSequence", 0, bitstream);
          else if(img->filterSequence[sub_pos] == 2)
            len += u_1("SH: filterSequence", 1, bitstream);
          else
            error("filterSequence should be less than 3", -1);
        }
      }
    }
    else if(img->filterParam == SIFO_FIRST_PASS_FPO)
    {    
      int sub_pos;

      len += u_v(1, "SH: use offsets", SIFO_USE_OFFSETS, bitstream);  // 1 bit only
      len += u_v (2, "SH: bestFilter", img->bestFilter, bitstream);
      if(img->bestFilter == NO_BEST_FILTER)
      {
        for(sub_pos = 1; sub_pos < 16; sub_pos++)
        {
          if(img->filterSequence[sub_pos] == 0)
            len += u_1("SH: filterSequence", img->filterSequence[sub_pos], bitstream);
          else
          {
            len += u_1("SH: filterSequence", 1, bitstream);
            if(img->filterSequence[sub_pos] == 1)
              len += u_1("SH: filterSequence", 0, bitstream);
            else if(img->filterSequence[sub_pos] == 2)
              len += u_1("SH: filterSequence", 1, bitstream);
            else
              error("filterSequence should be less than 3", -1);
          }
        }
      }
    } 
    else  // Frame filter
    {
      int sub_pos;

      if(img->filterParam == SIFO_FRAME_FILTER_WITH_OFFSET)
        len += u_v(1, "SH: use offsets", SIFO_USE_OFFSETS, bitstream);  // 1 bit only
      else
        len += u_v(1, "SH: use offsets", SIFO_FILTERS_ONLY, bitstream);   // 1 bit only

      len += u_v (2, "SH: bestFilter", NO_BEST_FILTER, bitstream);
      for(sub_pos = 1; sub_pos < 16; ++sub_pos)
      {
        if(img->filterFrame[sub_pos] == 0)
          len += u_1("SH: filterFrame", 0, bitstream);
        else 
        {
          len += u_1("SH: filterFrame", 1, bitstream);
          if(img->filterFrame[sub_pos] == 1)
            len += u_1("SH: filterFrame", 0, bitstream);
          else if(img->filterFrame[sub_pos] == 2)
            len += u_1("SH: filterFrame", 1, bitstream);
          else
            error("filterFrame should be less than 3", -1);
        }
      }
    }

    // Save offsets
    if((img->filterParam == SIFO_FRAME_FILTER_WITH_OFFSET) || 
      (img->filterParam == SIFO_FIRST_PASS_FPO))
    {
      int list;
      int listNo = (img->type == B_SLICE)? 2: 1;
      for(list = 0; list < listNo; ++list) 
      {
        int frame;
        int sub_pos;
        int nonzero = getNonzero(list);
        len += u_1("nonzero offset", nonzero, bitstream);
        if(nonzero)
        {
          for(frame = 0; frame < listXsize[list]; ++frame)
          {
            if(frame == 0)
            {     
              for(sub_pos = 0; sub_pos < 16; ++sub_pos)   
                len += se_v("SH: subpelOffset", img->subpelOffset[list][sub_pos], bitstream);
            }
            else
            {
              len += se_v("SH: imgOffset", img->imgOffset[list][frame], bitstream);
            }
          }
        }
      }
    }
  }
}
#endif  // SWITCHED_FILTERS

#ifdef ADAPTIVE_LOOP_FILTER
  currSlice->alf_flag_bit_to_go      = bitstream->bits_to_go;
  currSlice->alf_flag_stream_pointer = &(bitstream->streamBuffer[bitstream->byte_pos]);
  len += u_1("SH: adaptive_loop_filter_flag", 0, bitstream);

  // padding end of slice header for ALF
  if(bitstream->bits_to_go!=8)
  {
    while(bitstream->bits_to_go!=8)
    {
      len += u_1("SH: padding bits for ALF", 1, bitstream);
    }
  }
  currSlice->sh_byte_pos = bitstream->byte_pos;
#endif

  return len;
}

/*!
********************************************************************************************
* \brief 
*    Estimates the number of bits required for AIF coefficient transmitting 
*
* \return
*    Estimated number of bits used 
* \author
*    - Dmytro Rusanovskyy                   <dmytro.rusanovskyy@tut.fi>
********************************************************************************************
*/
#ifdef ADAPTIVE_FILTER
int nBitsAIF(void)
{
  int i, j;
  extern int UseAllSubpelPositions;                  // 1 if FilterFlag for all independent positions is 1
  extern int SubpelPositionsPattern;
#ifndef EDAIF2
  extern int FILTER_SIZE;
  extern int FILTER_OFFSET;
  extern int DiffQFilterCoef[15][SQR_FILTER];        // differences to be transmitted
  extern int POS_EQUATION_NUMBER[15];                // number of different coefficietns for each sub-pel position
  extern int FilterFlag[15];                        // Flags defining if Filter at the particular position calculated
#else
  extern int DiffQFilterCoef[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];        // differences to be transmitted
  extern int POS_EQUATION_NUMBER[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];                // number of different coefficietns for each sub-pel position
  extern int FilterFlag[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];                        // Flags defining if Filter at the particular position calculated
#endif
  
  Bitstream *bitstream;
  int len = 0;

  bitstream =  (Bitstream *) calloc(1, sizeof(Bitstream));
  bitstream->streamBuffer = (byte *) calloc(1000, sizeof(byte));


  len += u_1("SH: use_all_subpel_positions", UseAllSubpelPositions, bitstream);
  if(!UseAllSubpelPositions)
    len += ue_v("SH: positions_pattern", SubpelPositionsPattern, bitstream);

  if(FilterFlag[a_pos])
    for(i= 0; i < POS_EQUATION_NUMBER[a_pos]; i++)
      len += se_v("SH: a_pos", DiffQFilterCoef[a_pos][FILTER_SIZE*FILTER_OFFSET+i], bitstream);
  if(FilterFlag[b_pos])
    for(i= 0; i < POS_EQUATION_NUMBER[b_pos]; i++)
      len += se_v("SH: b_pos", DiffQFilterCoef[b_pos][FILTER_SIZE*FILTER_OFFSET+i], bitstream);
  if(FilterFlag[e_pos])
    for(i = 0; i < FILTER_SIZE; i++)
      for(j = i; j < FILTER_SIZE; j++)
        len += se_v("SH: e_pos", DiffQFilterCoef[e_pos][FILTER_SIZE*i+j], bitstream);
  if(FilterFlag[f_pos])
    for(i = 0; i < FILTER_SIZE; i++)
      for(j = 0; j < FILTER_SIZE/2; j++)
        len += se_v("SH: f_pos", DiffQFilterCoef[f_pos][FILTER_SIZE*i+j], bitstream);
  if(FilterFlag[j_pos])
    for(i = 0; i < FILTER_SIZE/2; i++)
      for(j = i; j < FILTER_SIZE/2; j++)
        len += se_v("SH: j_pos", DiffQFilterCoef[j_pos][FILTER_SIZE*i+j], bitstream);
  free(bitstream->streamBuffer);
  free(bitstream);
  return len;
}
// separable aif (BEGIN)
/*!
********************************************************************************************
* \brief 
*    Estimates the number of bits required for separable AIF coefficient transmitting 
*
* \return
*    Estimated number of bits used 
* \author
*    - Steffen Wittmann                   <steffen.wittmann-at-eu.panasonic.com>
********************************************************************************************
*/
int nBitsAIF_Sep(void)
{
  int i;
  extern int UseAllSubpelPositions;                  // 1 if FilterFlag for all independent positions is 1
  extern int SubpelPositionsPattern;
#ifndef EDAIF2
  extern int FILTER_SIZE;
  extern int DiffQFilterCoef[15][SQR_FILTER];        // differences to be transmitted
  extern int FilterFlag[15];                        // Flags defining if Filter at the particular position calculated
#else
  extern int DiffQFilterCoef[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];        // differences to be transmitted
  extern int FilterFlag[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];                        // Flags defining if Filter at the particular position calculated
#endif

  Bitstream *bitstream;
  int len = 0;

  bitstream =  (Bitstream *) calloc(1, sizeof(Bitstream));
  bitstream->streamBuffer = (byte *) calloc(1000, sizeof(byte));

  len += u_1("SH: use_all_subpel_positions", UseAllSubpelPositions, bitstream);
  if(!UseAllSubpelPositions)
    len += ue_v("SH: positions_pattern", SubpelPositionsPattern, bitstream);

  // a_pos
  if(FilterFlag[a_pos])
    for(i= 0; i < FILTER_SIZE; i++)
      len += se_v("SH: a_pos", DiffQFilterCoef[a_pos][i], bitstream);

  // b_pos        
  if(FilterFlag[b_pos])
    for(i= 0; i < FILTER_SIZE>>1; i++)
      len += se_v("SH: b_pos", DiffQFilterCoef[b_pos][i], bitstream);

  // c_pos             
  if(FilterFlag[c_pos])
  {
    if (FilterFlag[a_pos])
      for(i= 0; i < FILTER_SIZE; i++)
        len += se_v("SH: c_pos", DiffQFilterCoef[c_pos][i] - DiffQFilterCoef[a_pos][FILTER_SIZE-i-1], bitstream);
    else 
      for(i= 0; i < FILTER_SIZE; i++)
        len += se_v("SH: c_pos", DiffQFilterCoef[c_pos][i], bitstream);
  }

  // d_pos             
  if(FilterFlag[d_pos])
    for(i= 0; i < FILTER_SIZE; i++)
      len += se_v("SH: d_pos", DiffQFilterCoef[d_pos][i], bitstream);

  // e_pos             
  if(FilterFlag[e_pos])
  {
    if (FilterFlag[d_pos])
      for(i= 0; i < FILTER_SIZE; i++)
        len += se_v("SH: e_pos", DiffQFilterCoef[e_pos][i] - DiffQFilterCoef[d_pos][i], bitstream);
    else
      for(i= 0; i < FILTER_SIZE; i++)
        len += se_v("SH: e_pos", DiffQFilterCoef[e_pos][i], bitstream);
  }

  // f_pos            
  if(FilterFlag[f_pos])
  {
    if (FilterFlag[e_pos])
      for(i= 0; i < FILTER_SIZE; i++)
        len += se_v("SH: f_pos", DiffQFilterCoef[f_pos][i] - DiffQFilterCoef[e_pos][i], bitstream);
    else
      for(i= 0; i < FILTER_SIZE; i++)
        len += se_v("SH: f_pos", DiffQFilterCoef[f_pos][i], bitstream);
  }

  // g_pos             
  if(FilterFlag[g_pos])
  {
    if (FilterFlag[f_pos])
      for(i= 0; i < FILTER_SIZE; i++)
        len += se_v("SH: g_pos", DiffQFilterCoef[g_pos][i] - DiffQFilterCoef[f_pos][i], bitstream);
    else
      for(i= 0; i < FILTER_SIZE; i++)
        len += se_v("SH: g_pos", DiffQFilterCoef[g_pos][i], bitstream);
  }

  // h_pos             
  if(FilterFlag[h_pos])
    for(i= 0; i < FILTER_SIZE>>1; i++)
      len += se_v("SH: h_pos", DiffQFilterCoef[h_pos][i], bitstream);

  // i_pos     
  if(FilterFlag[i_pos])
  {
    if (FilterFlag[h_pos])
      for(i= 0; i < FILTER_SIZE>>1; i++)
        len += se_v("SH: i_pos", DiffQFilterCoef[i_pos][i] - DiffQFilterCoef[h_pos][i], bitstream);
    else 
      for(i= 0; i < FILTER_SIZE>>1; i++)
        len += se_v("SH: i_pos", DiffQFilterCoef[i_pos][i], bitstream);
  }

  // j_pos            
  if(FilterFlag[j_pos])
  {
    if (FilterFlag[i_pos])
      for(i= 0; i < FILTER_SIZE>>1; i++)
        len += se_v("SH: j_pos", DiffQFilterCoef[j_pos][i] - DiffQFilterCoef[i_pos][i], bitstream);
    else
      for(i= 0; i < FILTER_SIZE>>1; i++)
        len += se_v("SH: j_pos", DiffQFilterCoef[j_pos][i], bitstream);
  }

  // k_pos             
  if(FilterFlag[k_pos])
  {
    if (FilterFlag[j_pos])
      for(i= 0; i < FILTER_SIZE>>1; i++)
        len += se_v("SH: k_pos", DiffQFilterCoef[k_pos][i] - DiffQFilterCoef[j_pos][i], bitstream);
    else 
      for(i= 0; i < FILTER_SIZE>>1; i++)
        len += se_v("SH: k_pos", DiffQFilterCoef[k_pos][i], bitstream);
  }

  free(bitstream->streamBuffer);
  free(bitstream);
  return len;
}
// separable aif (END)
#endif
/*!
********************************************************************************************
* \brief 
*    writes the ref_pic_list_reordering syntax
*    based on content of according fields in img structure
*
* \return
*    number of bits used 
********************************************************************************************
*/
static int ref_pic_list_reordering(Bitstream *bitstream)
{
  Slice *currSlice = img->currentSlice;

  int i, len=0;

  // RPLR for redundant pictures
  if(input->redundant_pic_flag && redundant_coding)
  {
    currSlice->ref_pic_list_reordering_flag_l0 = 1;
    currSlice->reordering_of_pic_nums_idc_l0[0] = 0;
    currSlice->reordering_of_pic_nums_idc_l0[1] = 3;
    currSlice->abs_diff_pic_num_minus1_l0[0] = redundant_ref_idx - 1;
    currSlice->long_term_pic_idx_l0[0] = 0;
    reorder_ref_pic_list( listX[LIST_0], &listXsize[LIST_0], 
      img->num_ref_idx_l0_active-1, 
      currSlice->reordering_of_pic_nums_idc_l0, 
      currSlice->abs_diff_pic_num_minus1_l0, 
      currSlice->long_term_pic_idx_l0);
  }

  if ((img->type!=I_SLICE) /*&&(img->type!=SI_IMG)*/ )
  {
    len += u_1 ("SH: ref_pic_list_reordering_flag_l0", currSlice->ref_pic_list_reordering_flag_l0, bitstream);
    if (currSlice->ref_pic_list_reordering_flag_l0)
    {
      i=-1;
      do
      {
        i++;
        len += ue_v ("SH: reordering_of_pic_num_idc", currSlice->reordering_of_pic_nums_idc_l0[i], bitstream);
        if (currSlice->reordering_of_pic_nums_idc_l0[i]==0 ||
          currSlice->reordering_of_pic_nums_idc_l0[i]==1)
        {
          len += ue_v ("SH: abs_diff_pic_num_minus1_l0", currSlice->abs_diff_pic_num_minus1_l0[i], bitstream);
        }
        else
        {
          if (currSlice->reordering_of_pic_nums_idc_l0[i]==2)
          {
            len += ue_v ("SH: long_term_pic_idx_l0", currSlice->long_term_pic_idx_l0[i], bitstream);
          }
        }

      } while (currSlice->reordering_of_pic_nums_idc_l0[i] != 3);
    }
  }

  if (img->type==B_SLICE)
  {
    len += u_1 ("SH: ref_pic_list_reordering_flag_l1", currSlice->ref_pic_list_reordering_flag_l1, bitstream);
    if (currSlice->ref_pic_list_reordering_flag_l1)
    {
      i=-1;
      do
      {
        i++;
        len += ue_v ("SH: remapping_of_pic_num_idc", currSlice->reordering_of_pic_nums_idc_l1[i], bitstream);
        if (currSlice->reordering_of_pic_nums_idc_l1[i]==0 ||
          currSlice->reordering_of_pic_nums_idc_l1[i]==1)
        {
          len += ue_v ("SH: abs_diff_pic_num_minus1_l1", currSlice->abs_diff_pic_num_minus1_l1[i], bitstream);
        }
        else
        {
          if (currSlice->reordering_of_pic_nums_idc_l1[i]==2)
          {
            len += ue_v ("SH: long_term_pic_idx_l1", currSlice->long_term_pic_idx_l1[i], bitstream);
          }
        }
      } while (currSlice->reordering_of_pic_nums_idc_l1[i] != 3);
    }
  }

  return len;
}


/*!
************************************************************************
* \brief
*    write the memory management control operations
*
* \return
*    number of bits used 
************************************************************************
*/
static int dec_ref_pic_marking(Bitstream *bitstream)
{
  DecRefPicMarking_t *tmp_drpm;

  int val, len=0;

  if (img->currentPicture->idr_flag)
  {
    len += u_1("SH: no_output_of_prior_pics_flag", img->no_output_of_prior_pics_flag, bitstream);
    len += u_1("SH: long_term_reference_flag", img->long_term_reference_flag, bitstream);
  }
  else
  {
    img->adaptive_ref_pic_buffering_flag = (img->dec_ref_pic_marking_buffer!=NULL);

    len += u_1("SH: adaptive_ref_pic_buffering_flag", img->adaptive_ref_pic_buffering_flag, bitstream);

    if (img->adaptive_ref_pic_buffering_flag)
    {
      tmp_drpm = img->dec_ref_pic_marking_buffer;
      // write Memory Management Control Operation 
      do
      {
        if (tmp_drpm==NULL) error ("Error encoding MMCO commands", 500);

        val = tmp_drpm->memory_management_control_operation;
        len += ue_v("SH: memory_management_control_operation", val, bitstream);

        if ((val==1)||(val==3)) 
        {
          len += 1 + ue_v("SH: difference_of_pic_nums_minus1", tmp_drpm->difference_of_pic_nums_minus1, bitstream);
        }
        if (val==2)
        {
          len+= ue_v("SH: long_term_pic_num", tmp_drpm->long_term_pic_num, bitstream);
        }
        if ((val==3)||(val==6))
        {
          len+= ue_v("SH: long_term_frame_idx", tmp_drpm->long_term_frame_idx, bitstream);
        }
        if (val==4)
        {
          len += ue_v("SH: max_long_term_pic_idx_plus1", tmp_drpm->max_long_term_frame_idx_plus1, bitstream);
        }

        tmp_drpm=tmp_drpm->Next;

      } while (val != 0);

    }
  }
  return len;
}


/*!
************************************************************************
* \brief
*    write prediction weight table
*
* \return
*    number of bits used 
************************************************************************
*/
static int pred_weight_table(Bitstream *bitstream)
{
  int len = 0;
  int i,j;

  len += ue_v("SH: luma_log_weight_denom", luma_log_weight_denom, bitstream);

  if ( 0 != active_sps->chroma_format_idc)
  {
    len += ue_v("SH: chroma_log_weight_denom", chroma_log_weight_denom, bitstream);
  }

  for (i=0; i< img->num_ref_idx_l0_active; i++)
  {
    if ( (wp_weight[0][i][0] != 1<<luma_log_weight_denom) || (wp_offset[0][i][0] != 0) )
    {
      len += u_1 ("SH: luma_weight_flag_l0", 1, bitstream);

      len += se_v ("SH: luma_weight_l0", wp_weight[0][i][0], bitstream);

#ifdef  BUG_FIX_FOR_FREXT
      len += se_v ("SH: luma_offset_l0", (wp_offset[0][i][0]>>(img->bitdepth_luma-8)), bitstream);
#else
      len += se_v ("SH: luma_offset_l0", wp_offset[0][i][0], bitstream);
#endif
    }
    else
    {
      len += u_1 ("SH: luma_weight_flag_l0", 0, bitstream);
    }

    if (active_sps->chroma_format_idc!=0)
    {
      if ( (wp_weight[0][i][1] != 1<<chroma_log_weight_denom) || (wp_offset[0][i][1] != 0) || 
        (wp_weight[0][i][2] != 1<<chroma_log_weight_denom) || (wp_offset[0][i][2] != 0)  )
      {
        len += u_1 ("chroma_weight_flag_l0", 1, bitstream);
        for (j=1; j<3; j++)
        {
          len += se_v ("chroma_weight_l0", wp_weight[0][i][j] ,bitstream);

#ifdef  BUG_FIX_FOR_FREXT
          len += se_v ("chroma_offset_l0", (wp_offset[0][i][j]>>(img->bitdepth_chroma-8)) ,bitstream);
#else
          len += se_v ("chroma_offset_l0", wp_offset[0][i][j] ,bitstream);
#endif
        }
      }
      else
      {
        len += u_1 ("chroma_weight_flag_l0", 0, bitstream);
      }
    }
  }

  if (img->type == B_SLICE)
  {
    for (i=0; i< img->num_ref_idx_l1_active; i++)
    {
      if ( (wp_weight[1][i][0] != 1<<luma_log_weight_denom) || (wp_offset[1][i][0] != 0) )
      {
        len += u_1 ("SH: luma_weight_flag_l1", 1, bitstream);

        len += se_v ("SH: luma_weight_l1", wp_weight[1][i][0], bitstream);

#ifdef  BUG_FIX_FOR_FREXT
        len += se_v ("SH: luma_offset_l1", (wp_offset[1][i][0]>>(img->bitdepth_luma-8)), bitstream);
#else
        len += se_v ("SH: luma_offset_l1", wp_offset[1][i][0], bitstream);
#endif
      }
      else
      {
        len += u_1 ("SH: luma_weight_flag_l1", 0, bitstream);
      }

      if (active_sps->chroma_format_idc!=0)
      {
        if ( (wp_weight[1][i][1] != 1<<chroma_log_weight_denom) || (wp_offset[1][i][1] != 0) || 
          (wp_weight[1][i][2] != 1<<chroma_log_weight_denom) || (wp_offset[1][i][2] != 0) )
        {
          len += u_1 ("chroma_weight_flag_l1", 1, bitstream);
          for (j=1; j<3; j++)
          {
            len += se_v ("chroma_weight_l1", wp_weight[1][i][j] ,bitstream);

#ifdef  BUG_FIX_FOR_FREXT
            len += se_v ("chroma_offset_l1", (wp_offset[1][i][j]>>(img->bitdepth_chroma-8)) ,bitstream);
#else
            len += se_v ("chroma_offset_l1", wp_offset[1][i][j] ,bitstream);
#endif
          }
        }
        else
        {
          len += u_1 ("chroma_weight_flag_l1", 0, bitstream);
        }
      }
    }
  }
  return len;
}


/*!
************************************************************************
* \brief
*    Selects picture type and codes it to symbol
*
* \return
*    symbol value for picture type
************************************************************************
*/
int get_picture_type()
{
  // set this value to zero for transmission without signaling 
  // that the whole picture has the same slice type
  int same_slicetype_for_whole_frame = 5;

  switch (img->type)
  {
  case I_SLICE:
    return 2 + same_slicetype_for_whole_frame;
    break;
  case P_SLICE:
    return 0 + same_slicetype_for_whole_frame;
    break;
  case B_SLICE:
    return 1 + same_slicetype_for_whole_frame;
    break;
  case SP_SLICE:
    return 3 + same_slicetype_for_whole_frame;
    break;
  default:
    error("Picture Type not supported!",1);
    break;
  }

  return 0;
}



/*!
*****************************************************************************
*
* \brief 
*    int Partition_BC_Header () write the Partition type B, C header
*
* \return
*    Number of bits used by the partition header
*
* \par Parameters
*    PartNo: Partition Number to which the header should be written
*
* \par Side effects
*    Partition header as per VCEG-N72r2 is written into the appropriate 
*    partition bit buffer
*
* \par Limitations/Shortcomings/Tweaks
*    The current code does not support the change of picture parameters within
*    one coded sequence, hence there is only one parameter set necessary.  This
*    is hard coded to zero.
*
* \date
*    October 24, 2001
*
* \author
*    Stephan Wenger   stewe@cs.tu-berlin.de
*****************************************************************************/
int Partition_BC_Header(int PartNo)
{
  DataPartition *partition = &((img->currentSlice)->partArr[PartNo]);
  SyntaxElement symbol, *sym = &symbol;

  int len = 0;

  assert (PartNo > 0 && PartNo < img->currentSlice->max_part_nr);

  sym->type = SE_HEADER;         // This will be true for all symbols generated here
  sym->mapping = ue_linfo;       // Mapping rule: Simple code number to len/info
  sym->value2  = 0;

  //ZL 
  //changed according to the g050r1
  SYMTRACESTRING("RTP-PH: Slice ID");
  sym->value1 = img->current_slice_nr;
  len += writeSyntaxElement_UVLC (sym, partition);

  return len;
}

#ifdef E_DAIF
void generatePrefixCode_EDAIF(int *prefixCode, int *prefixCodeLen, int size)
{
  int i; 

  prefixCode[0] = 0;
  prefixCodeLen[0] = 1;

  for(i = 1; i < size-1; i++)
  {
    prefixCode[i] = (1<<(i+1))-2; 
    prefixCodeLen[i] = i+1;
  }

  prefixCode[size-1] = prefixCode[size-2]+1;
  prefixCodeLen[size-1] = prefixCodeLen[size-2];
}

void generateThres_EDAIF(int *thres, int max, int size)
{
  int i; 

  thres[0] = 0; 
  for(i = 0; i < size; i++)
    thres[i+1] = max >> (size-1-i);
}

int writeCoeff_EDAIF(int coeffQ, Bitstream *bitstream, int numQBits, int depth, char *tracestring)
{
  int prefixCode[100];
  int prefixCodeLen[100]; 
  int thres[100];

  int sign = (coeffQ >= 0) ? 0:1; 
  int mag = (sign==0) ? coeffQ:(-coeffQ); 

  int rangeQ = 1<<(numQBits-1); 

  int len = 0;
  int i; 
  int suffixLen; 
  int bin;

  generatePrefixCode_EDAIF(prefixCode, prefixCodeLen, depth);
  generateThres_EDAIF(thres, rangeQ, depth); 

  // find to which bin this value "mag" belongs
  for(i = 0; i < depth; i++)
  {
    if(mag >= thres[i] && mag < thres[i+1])
      break; 
  }
  bin = i;
  assert(bin < depth);

  suffixLen = numQBits-2;
  for(i = depth-1; i > bin; i--)
    suffixLen--;
  if(bin == 0)
    suffixLen++;

  len += u_v(prefixCodeLen[bin], tracestring, prefixCode[bin], bitstream); 
  len += u_v(suffixLen, tracestring, mag-thres[bin], bitstream); 
  // sign bit 
  if(mag) 
    len += u_1(tracestring, sign, bitstream); 

  return len; 
}

int sendAIFInteger(Bitstream *bitstream)
{
  int newStreamFlag; 
  int i; 
  int len = 0; 

  extern int DiffQFilterCoeffInt[SQR_FILTER_INT]; 
  extern int FilterFlagInt; 
  extern int numQBitsInt[SQR_FILTER_INT-1];

  if(!FilterFlagInt)
    return len; 

  newStreamFlag = 0;
  if (bitstream==NULL)
  {
    bitstream =  (Bitstream *) calloc(1, sizeof(Bitstream));
    bitstream->streamBuffer = (byte *) calloc(1000, sizeof(byte));
    newStreamFlag = 1;
  }

  for(i = 0; i < FILTER_SIZE_INT*FILTER_SIZE_INT; i++)
  {
    len += writeCoeff_EDAIF(DiffQFilterCoeffInt[i], bitstream, numQBitsInt[i], 6, "SH: full-pel filter");
  }

  if (newStreamFlag)
  {
    free(bitstream->streamBuffer);
    free(bitstream);
  }

  return len;
}

int sendAIFOffset(int sub_pos, Bitstream *bitstream)
{
  int newStreamFlag; 
  int len = 0; 
  int offsetI, offsetF; 
  int sign; 
  int offsetFracCodeLen[] = 
  {
    5, 4, 4, 2, 2, 2, 2, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
  };

  extern int FilterFlagInt; 
  extern double FilterCoef[15][SQR_FILTER];
  extern double FilterCoefInt [SQR_FILTER_INT];
  extern int DiffQFilterOffsetI[15], DiffQFilterOffsetF[15];
  extern int DiffQFilterOffsetIntI, DiffQFilterOffsetIntF; 
#ifndef EDAIF2
  extern int FILTER_SIZE; 
  extern int FilterFlag[15]; 
  extern int SymmetryPosition[15];
#else
  extern int FilterFlag[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA]; 
  extern int SymmetryPosition[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];
#endif

  newStreamFlag = 0;
  if (bitstream==NULL)
  {
    bitstream =  (Bitstream *) calloc(1, sizeof(Bitstream));
    bitstream->streamBuffer = (byte *) calloc(1000, sizeof(byte));
    newStreamFlag = 1;
  }

  if(sub_pos == -1) // full-pel offset
  {
    if(FilterFlagInt)
    {
      offsetI = DiffQFilterOffsetIntI; 
      offsetF = DiffQFilterOffsetIntF; 
      sign = FilterCoefInt[SQR(FILTER_SIZE_INT)] >= 0. ? 0:1;
      // send integer, using EXP-GOLOMB
      len += se_v("SH: offset int", offsetI, bitstream);

      // send fraction, using fixed-length code 
      if(offsetF != -1)
      {
        len += u_v(offsetFracCodeLen[offsetI], "SH: offset frac", offsetF, bitstream); 
      }
      if(offsetI || (offsetF > 0)) 
        len += u_1("SH: offset sign", sign, bitstream);
    }
  }
#ifndef EDAIF2
  else if(FilterFlag[sub_pos] && SymmetryPosition[sub_pos]) // subpel offset 
#else
  else if((FilterFlag[sub_pos]%2!=0) && SymmetryPosition[sub_pos]) // if filter is adaptive and primary 
#endif
  {
    offsetI = DiffQFilterOffsetI[sub_pos]; 
    offsetF = DiffQFilterOffsetF[sub_pos]; 
    sign = FilterCoef[sub_pos][SQR(FILTER_SIZE)] >= 0. ? 0:1;
    // send integer, using EXP-GOLOMB
    len += se_v("SH: offset int", offsetI, bitstream);

    // send fraction, using fixed-length code 
    if(offsetF != -1)
    {
      len += u_v(offsetFracCodeLen[offsetI], "SH: offset frac", offsetF, bitstream); 
    }
    if(offsetI || (offsetF > 0)) 
      len += u_1("SH: offset sign", sign, bitstream);
  }

  if (newStreamFlag)
  {
    free(bitstream->streamBuffer);
    free(bitstream);
  }

  return len; 
}
#endif

#ifdef DIRECTIONAL_FILTER
/*! 
*************************************************************************************
* \brief
*   Send 1D-AIF coefficient to the bitstream
*  
* \para <title>
*    <paragraph>
*
* \para
*    <another paragraph>
*
* \author
*    - Dmytro Rusanovskyy                   <dmytro.rusanovskyy@tut.fi>
*************************************************************************************
*/
int sendCoefsAIF(Bitstream *bitstream)
{
  int i, j;
  int newStreamFlag;
  extern int UseAllSubpelPositions;                  // 1 if FilterFlag for all independent positions is 1
  extern int SubpelPositionsPattern;
#ifndef EDAIF2
  extern int FILTER_SIZE;
  extern int FILTER_OFFSET;
  extern int FilterFlag[15]; 
  extern int SymmetryPosition[15];
  extern int DiffQFilterCoef[15][SQR_FILTER];        // differences to be transmitted
  extern int POS_EQUATION_NUMBER[15];                // number of different coefficietns for each sub-pel position
#else
  extern int FilterFlag[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA]; 
  extern int SymmetryPosition[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];
  extern int DiffQFilterCoef[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];        // differences to be transmitted
  extern int POS_EQUATION_NUMBER[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];                // number of different coefficietns for each sub-pel position
#endif

  //Bitstream *bitstream;
  int len = 0;

  newStreamFlag = 0;
  if (bitstream==NULL)
  {
    bitstream =  (Bitstream *) calloc(1, sizeof(Bitstream));
    bitstream->streamBuffer = (byte *) calloc(1000, sizeof(byte));
    newStreamFlag = 1;
  }

  len += u_1("SH: use_all_subpel_positions", UseAllSubpelPositions, bitstream);

  //  If simmetry is assumed FILTER_TYPE_Z17 transmits pattern otherwise  we transmit filterFlag[]
  if(!UseAllSubpelPositions)
    len += ue_v("SH: positions_pattern", SubpelPositionsPattern, bitstream);

  if(FilterFlag[a_pos])
    for(i= 0; i < POS_EQUATION_NUMBER[a_pos]; i++)
      len += se_v("SH: a_pos", DiffQFilterCoef[a_pos][FILTER_SIZE*FILTER_OFFSET+i], bitstream);
  if(FilterFlag[b_pos])
    for(i= 0; i < POS_EQUATION_NUMBER[b_pos]; i++)
      len += se_v("SH: b_pos", DiffQFilterCoef[b_pos][FILTER_SIZE*FILTER_OFFSET+i], bitstream);
  if(FilterFlag[e_pos])
    for(i = 0; i < FILTER_SIZE; i++)
    {
      j = i;//NW-SE
      len += se_v("SH: e_pos", DiffQFilterCoef[e_pos][FILTER_SIZE*i+j], bitstream);
    }
    if(FilterFlag[f_pos])
      for(i = 0; i < FILTER_SIZE; i++)
      {
        j = i;//NW-SE
        len += se_v("SH: f_pos", DiffQFilterCoef[f_pos][FILTER_SIZE*i+j], bitstream);
      }
      if(FilterFlag[j_pos])
        for(i = 0; i < FILTER_SIZE/2; i++)
        {
          j = i;//NW-SE
          len += se_v("SH: j_pos", DiffQFilterCoef[j_pos][FILTER_SIZE*i+j], bitstream);
        }

        if (newStreamFlag)
        {
          free(bitstream->streamBuffer);
          free(bitstream);
        }
        // printf("nAIF: %d\n",len);
        return len;
}


#ifdef E_DAIF
#ifdef EDAIF2
int sendEDAIF2(Bitstream *bitstream)
{
    extern int FilterFlag[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];													// Flags defining if Filter at the particular position calculated
    extern int SymmetryPosition[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];	// if 0, the position is copied from another one
    extern int DiffQFilterCoef[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];				// differences to be transmitted
    extern int POS_EQUATION_NUMBER[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];								// number of different coefficietns for each sub-pel position
    extern int Calc2Filt_Indexes[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];

    extern int RealTimeDecomposition;
    extern int  UseAllSubpelPositions;                  // 1 if FilterFlag for all independent positions is 1
    extern int nPairsSets[4];
    extern int realSymCommands[4][MAX_NUM_AIF];

    extern int numQBits1DH[], numQBits1DV[], numQBitsDIF[], numQBits2D[];

    int i,j, sub_pos, nQBits;
    unsigned len = 0;
    int newStreamFlag;

    len = 0;
    newStreamFlag = 0;
    if (bitstream==NULL)
    {
      bitstream =  (Bitstream *) calloc(1, sizeof(Bitstream));
      bitstream->streamBuffer = (byte *) calloc(1000, sizeof(byte));
      newStreamFlag = 1;
    }

    UseAllSubpelPositions = 0;// we use different filters here!
    for(sub_pos = 0; sub_pos <= o_pos; sub_pos++)
    {
      if((sub_pos == f_pos)||(sub_pos == i_pos)||(sub_pos == j_pos)||(sub_pos == k_pos)||(sub_pos == n_pos))
      {
        if((FilterFlag[sub_pos]==0)||(FilterFlag[sub_pos]==2))
          len += u_v(1, "SH: FilterFlag[sub_pos]", 0, bitstream);
        else
        {
          if(FilterFlag[sub_pos]==1)
            len += u_v(2, "SH: FilterFlag[sub_pos]", 2, bitstream); // Filter 1: codeword '10'
          else 
            len += u_v(2, "SH: FilterFlag[sub_pos]", 3, bitstream); // Filter 3: codeword '11'
        }
      }
      else
        len += u_v(1, "SH: FilterFlag[sub_pos]", FilterFlag[sub_pos], bitstream);
    }

    len += u_1("SH: RealTimeDecomposition", RealTimeDecomposition, bitstream); 
    if (RealTimeDecomposition == 1)
    {
      int counter = 0;
      for(i = 0;i<4;i++)
      {
        for (j = 0; j < nPairsSets[i]; j ++)
        {
          counter += realSymCommands[i][j];
        }
      }
      counter = (counter==20)?1:0;
      len += u_1("SH: RTD_Full", counter, bitstream); 
      if(!counter)
      {
        for(i = 0;i<4;i++)
        {
          for (j = 0; j < nPairsSets[i]; j ++)
          {
            len += u_1("SH: realSymCommands[i][j]", realSymCommands[i][j], bitstream); 
          }
        }
      }
    }

    for (sub_pos = 0; sub_pos < MAX_NUM_AIF; sub_pos++)
    {
      if((FilterFlag[sub_pos]%2!=0)&&(SymmetryPosition[sub_pos]))// if filter is adaptive and primary 
      {
        for(i= 0; i < POS_EQUATION_NUMBER[sub_pos]; i++)          
        {
          //if(DiffQFilterCoef[sub_pos][Calc2Filt_Indexes[sub_pos][i]]!=0)
          {
            if((sub_pos == a_pos)||(sub_pos == b_pos)||(sub_pos == c_pos))
              nQBits = numQBits1DH[Calc2Filt_Indexes[sub_pos][i]];  
            else if((sub_pos == d_pos)||(sub_pos == h_pos)||(sub_pos == l_pos))
              nQBits = numQBits1DV[Calc2Filt_Indexes[sub_pos][i]];  
            else if((sub_pos == e_pos)||(sub_pos == g_pos)||(sub_pos == m_pos)||(sub_pos == o_pos))
              nQBits = numQBitsDIF[Calc2Filt_Indexes[sub_pos][i]];  
            else if(FilterFlag[sub_pos]<=1)// DAIF
                nQBits = numQBitsDIF[Calc2Filt_Indexes[sub_pos][i]];  
            else 
              nQBits = numQBits2D[Calc2Filt_Indexes[sub_pos][i]];
            assert(nQBits != 0);
            len += writeCoeff_EDAIF(DiffQFilterCoef[sub_pos][Calc2Filt_Indexes[sub_pos][i]], 
              bitstream, nQBits, 6, "SH: DiffQFilterCoef");
          }
        }
      }
    }
    if (newStreamFlag)
    {
      free(bitstream->streamBuffer);
      free(bitstream);
    }

    //PrintFilterCoefInt(DiffQFilterCoef);
    //printf("EDAIF2: %d\n",len);
     return len;
}


#endif
int sendCoefs_EDAIF(Bitstream *bitstream)
{
  int len;
  int sub_pos;
  int newStreamFlag;
  extern int FilterFlagInt; 

  len = 0;
  newStreamFlag = 0;
  if (bitstream==NULL)
  {
    bitstream =  (Bitstream *) calloc(1, sizeof(Bitstream));
    bitstream->streamBuffer = (byte *) calloc(1000, sizeof(byte));
    newStreamFlag = 1;
  }

    len += sendEDAIF2(bitstream);


  // full-pel filter
  len += u_1("SH: FilterFlagInt", FilterFlagInt, bitstream); 
  len += sendAIFInteger(bitstream); 

  // sub-pel filter offsets
  for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
  {
    len += sendAIFOffset(sub_pos, bitstream); 
  }

  // full-pel filter offset
  len += sendAIFOffset(-1, bitstream); 

  if (newStreamFlag)
  {
    free(bitstream->streamBuffer);
    free(bitstream);
  }
  return len;
}

//
//  Estimates the cost in bit of sending the filter corresponding to a sub-pixel position
//
int estimateCostOfCoefs_EDAIF(int sub_pos)
{
#ifndef EDAIF2
  extern int FILTER_SIZE;
  extern int FILTER_OFFSET;
#endif
  extern int FilterFlag[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];													// Flags defining if Filter at the particular position calculated
  extern int SymmetryPosition[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];	// if 0, the position is copied from another one
  extern int DiffQFilterCoef[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];				// differences to be transmitted
  extern int Calc2Filt_Indexes[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];
  extern int POS_EQUATION_NUMBER[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];								// number of different coefficietns for each sub-pel position
  extern int numQBits1DH[], numQBits1DV[], numQBitsDIF[], numQBits2D[];
    int nQBits;

  int i,j;
  unsigned len = 0;
  int newStreamFlag;

  Bitstream *bitstream = (Bitstream *) calloc(1, sizeof(Bitstream));
  bitstream->streamBuffer = (byte *) calloc(1000, sizeof(byte));



  len = 0;
  newStreamFlag = 0;
  if (bitstream==NULL)
  {
    bitstream =  (Bitstream *) calloc(1, sizeof(Bitstream));
    bitstream->streamBuffer = (byte *) calloc(1000, sizeof(byte));
    newStreamFlag = 1;
  }

  j = sub_pos;

  if((FilterFlag[sub_pos]%2!=0)&&SymmetryPosition[sub_pos])// if filter is adaptive and primary 
  {
    for(i= 0; i < POS_EQUATION_NUMBER[sub_pos]; i++)          
    {
      {
        if((sub_pos == a_pos)||(sub_pos == b_pos)||(sub_pos == c_pos))
          nQBits = numQBits1DH[Calc2Filt_Indexes[sub_pos][i]];  
        else if((sub_pos == d_pos)||(sub_pos == h_pos)||(sub_pos == l_pos))
          nQBits = numQBits1DV[Calc2Filt_Indexes[sub_pos][i]];  
        else if((sub_pos == e_pos)||(sub_pos == g_pos)||(sub_pos == m_pos)||(sub_pos == o_pos))
          nQBits = numQBitsDIF[Calc2Filt_Indexes[sub_pos][i]];  
        else if(FilterFlag[sub_pos]<=1)// DAIF
            nQBits = numQBitsDIF[Calc2Filt_Indexes[sub_pos][i]];  
        else 
          nQBits = numQBits2D[Calc2Filt_Indexes[sub_pos][i]];
        assert(nQBits != 0);
        len += writeCoeff_EDAIF(DiffQFilterCoef[sub_pos][Calc2Filt_Indexes[sub_pos][i]], 
          bitstream, nQBits, 6, "SH: DiffQFilterCoef");
      }
    }
  }

  free(bitstream->streamBuffer);
  free(bitstream);

  // printf("nAIF: %d\n",len);
  return len;
}

#endif  //  E_DAIF

#ifdef EAIF
const int UniqueCoeff[15][SQR_FILTER] = 
{
  {
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    1, 1, 1, 1, 1, 1, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
  },  // a_pos
  {
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    1, 1, 1, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
  },  // b_pos
  {
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    1, 1, 1, 1, 1, 1, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
  },  // c_pos
  {
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    1, 1, 1, 1, 1, 1, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
  },  // d_pos
  {
    0, 0, 0, 0, 0, 0, 
    0, 0, 1, 1, 0, 0, 
    0, 1, 1, 1, 1, 0, 
    0, 1, 1, 1, 1, 0, 
    0, 0, 1, 1, 0, 0, 
    0, 0, 0, 0, 0, 0, 
  },  // e_pos
  {
    0, 0, 0, 0, 0, 0, 
    0, 0, 1, 0, 0, 0, 
    0, 1, 1, 0, 0, 0, 
    0, 1, 1, 0, 0, 0, 
    0, 0, 1, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
  },  // f_pos
  {
    0, 0, 0, 0, 0, 0, 
    0, 0, 1, 1, 0, 0, 
    0, 1, 1, 1, 1, 0, 
    0, 1, 1, 1, 1, 0, 
    0, 0, 1, 1, 0, 0, 
    0, 0, 0, 0, 0, 0, 
  },  // g_pos
  {
    0, 0, 1, 0, 0, 0, 
    0, 0, 1, 0, 0, 0, 
    0, 0, 1, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
  },  // h_pos
  {
    0, 0, 0, 0, 0, 0, 
    0, 0, 1, 1, 0, 0, 
    0, 1, 1, 1, 1, 0, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
  },  // i_pos
  {
    0, 0, 0, 0, 0, 0, 
    0, 0, 1, 0, 0, 0, 
    0, 1, 1, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
  },  // j_pos
  {
    0, 0, 0, 0, 0, 0, 
    0, 0, 1, 1, 0, 0, 
    0, 1, 1, 1, 1, 0, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
  },  // k_pos
  {
    0, 0, 1, 0, 0, 0, 
    0, 0, 1, 0, 0, 0, 
    0, 0, 1, 0, 0, 0, 
    0, 0, 1, 0, 0, 0, 
    0, 0, 1, 0, 0, 0, 
    0, 0, 1, 0, 0, 0, 
  },  // l_pos
  {
    0, 0, 0, 0, 0, 0, 
    0, 0, 1, 1, 0, 0, 
    0, 1, 1, 1, 1, 0, 
    0, 1, 1, 1, 1, 0, 
    0, 0, 1, 1, 0, 0, 
    0, 0, 0, 0, 0, 0, 
  },  // m_pos
  {
    0, 0, 0, 0, 0, 0, 
    0, 0, 1, 0, 0, 0, 
    0, 1, 1, 0, 0, 0, 
    0, 1, 1, 0, 0, 0, 
    0, 0, 1, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
  },  // n_pos
  {
    0, 0, 0, 0, 0, 0, 
    0, 0, 1, 1, 0, 0, 
    0, 1, 1, 1, 1, 0, 
    0, 1, 1, 1, 1, 0, 
    0, 0, 1, 1, 0, 0, 
    0, 0, 0, 0, 0, 0, 
  },  // o_pos
}; 

// estimate the cost of sending each EAIF sub-pel filter
int estimateCostOfCoefs_EAIF(int sub_pos)
{
#ifndef EDAIF2
  extern int FILTER_SIZE;
  extern int FILTER_OFFSET;
  extern int DiffQFilterCoef[15][SQR_FILTER];       // differences to be transmitted
  extern int POS_EQUATION_NUMBER[15];               // number of different coefficietns for each sub-pel position
#else
  extern int DiffQFilterCoef[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];       // differences to be transmitted
  extern int POS_EQUATION_NUMBER[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];               // number of different coefficietns for each sub-pel position
#endif

  extern int numQBits1DH[SQR_FILTER-1];
  extern int numQBits1DV[SQR_FILTER-1];
  extern int numQBits2D [SQR_FILTER-1];

  int len = 0;
  int i;
  Bitstream *bitstream = (Bitstream *) calloc(1, sizeof(Bitstream));

  bitstream->streamBuffer = (byte *) calloc(1000, sizeof(byte));

  if((sub_pos+1)/4==0)  // 1D horizontal
  {
    for(i = 0; i < POS_EQUATION_NUMBER[sub_pos]; i++)
    {
      len += writeCoeff_EDAIF(DiffQFilterCoef[sub_pos][FILTER_SIZE*FILTER_OFFSET+i], 
      bitstream, numQBits1DH[FILTER_SIZE*FILTER_OFFSET+i], 6, "SH: sub_pos 1D horiz");
    }
  }
  else if((sub_pos+1)%4 == 0) // 1D vertical
  {
    for(i = 0; i < POS_EQUATION_NUMBER[sub_pos]; i++)
      len += writeCoeff_EDAIF(DiffQFilterCoef[sub_pos][FILTER_SIZE*i+FILTER_OFFSET], 
      bitstream, numQBits1DV[FILTER_SIZE*i+FILTER_OFFSET], 6, "SH: sub_pos 1D vert");
  }
  else // 2D 
  {
    for(i = 0; i < SQR_FILTER-1; i++)
    {
      if(UniqueCoeff[sub_pos][i])
        len += writeCoeff_EDAIF(DiffQFilterCoef[sub_pos][i], bitstream, numQBits2D[i], 6, "SH: sub_pos 2D");
    }
  }

  free(bitstream->streamBuffer);
  free(bitstream);

  // printf("nAIF: %d\n",len);
  return len;
}

// send EAIF filter for each subpel position
int sendAIFSubpel(int sub_pos, Bitstream *bitstream)
{
  int newStreamFlag; 
  int i; 
  int len = 0; 

  extern int numQBits1DH[SQR_FILTER-1];
  extern int numQBits1DV[SQR_FILTER-1];
  extern int numQBits2D [SQR_FILTER-1];
  
#ifndef EDAIF2
  extern int FILTER_SIZE;
  extern int FILTER_OFFSET;
  extern int FilterFlag[15];
  extern int SymmetryPosition[15];
  extern int POS_EQUATION_NUMBER[15];
  extern int DiffQFilterCoef[15][SQR_FILTER];
#else
  extern int FilterFlag[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];
  extern int SymmetryPosition[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];
  extern int POS_EQUATION_NUMBER[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA];
  extern int DiffQFilterCoef[MAX_NUM_AIF+MAX_NUM_AIF_EXTRA][SQR_FILTER];
#endif

  newStreamFlag = 0;
  if (bitstream==NULL)
  {
    bitstream =  (Bitstream *) calloc(1, sizeof(Bitstream));
    bitstream->streamBuffer = (byte *) calloc(1000, sizeof(byte));
    newStreamFlag = 1;
  }
  if(FilterFlag[sub_pos] && SymmetryPosition[sub_pos])
  {
    if((sub_pos+1)/4==0)  // 1D horizontal
    {
      for(i = 0; i < POS_EQUATION_NUMBER[sub_pos]; i++)
      {
        len += writeCoeff_EDAIF(DiffQFilterCoef[sub_pos][FILTER_SIZE*FILTER_OFFSET+i], 
        bitstream, numQBits1DH[FILTER_SIZE*FILTER_OFFSET+i], 6, "SH: sub_pos 1D horiz");
      }
    }
    else if((sub_pos+1)%4 == 0) // 1D vertical
    {
      for(i = 0; i < POS_EQUATION_NUMBER[sub_pos]; i++)
        len += writeCoeff_EDAIF(DiffQFilterCoef[sub_pos][FILTER_SIZE*i+FILTER_OFFSET], 
        bitstream, numQBits1DV[FILTER_SIZE*i+FILTER_OFFSET], 6, "SH: sub_pos 1D vert");
    }
    else // 2D 
    {
      for(i = 0; i < SQR_FILTER-1; i++)
        if(UniqueCoeff[sub_pos][i])
          len += writeCoeff_EDAIF(DiffQFilterCoef[sub_pos][i], bitstream, numQBits2D[i], 6, "SH: sub_pos 2D");
    }
  }

  if (newStreamFlag)
	{
		free(bitstream->streamBuffer);
		free(bitstream);
	}

  return len; 
}

// send all EAIF filter coefficients (subpel, offset, and fullpel)
int sendCoefs_EAIF(Bitstream *bitstream)
{
  int len;
  int sub_pos;
  int newStreamFlag;

  extern int SymmetryPosition[]; 
  extern int FilterFlag[];
  extern int FilterFlagInt; 

  len = 0;
  newStreamFlag = 0;
  if (bitstream==NULL)
  {
    bitstream =  (Bitstream *) calloc(1, sizeof(Bitstream));
    bitstream->streamBuffer = (byte *) calloc(1000, sizeof(byte));
    newStreamFlag = 1;
  }

  // filter flags
  for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
  {
    if(SymmetryPosition[sub_pos])
    	len += u_1("SH: filter_flag", FilterFlag[sub_pos], bitstream);
  }
  // filter coeffs
  for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
  {
    len += sendAIFSubpel(sub_pos, bitstream); 
  }

  // full-pel filter
  len += u_1("SH: FilterFlagInt", FilterFlagInt, bitstream); 
  len += sendAIFInteger(bitstream); 

  // sub-pel filter offsets
  for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
  {
    len += sendAIFOffset(sub_pos, bitstream); 
  }

  // full-pel filter offset
  len += sendAIFOffset(-1, bitstream); 

  if (newStreamFlag)
  {
    free(bitstream->streamBuffer);
    free(bitstream);
  }
  // printf("nAIF: %d\n",len);
  return len;
}

#endif

#endif
