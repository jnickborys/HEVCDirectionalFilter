
/*!
*************************************************************************************
* \file header.c
*
* \brief
*    H.264 Slice headers
*
*************************************************************************************
*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "global.h"
#include "elements.h"
#include "defines.h"
#include "fmo.h"
#include "vlc.h"
#include "mbuffer.h"
#include "header.h"

#include "ctx_tables.h"

#ifdef ADAPTIVE_FILTER
#include "adaptive_filter.h"
#endif

#ifdef ADAPTIVE_QUANTIZATION
#include "adaptive_quantization.h"
#endif

#ifdef ADAPTIVE_LOOP_FILTER
#include "adaptive_loop_filter.h"
#endif

#ifdef SWITCHED_FILTERS
#include "switched_filters.h"
#endif

#ifdef EDAIF2
#include "aif_common.h" // we collect in a single file defines related to encoder and decoder
#endif

extern StorablePicture *dec_picture;

#if TRACE
#define SYMTRACESTRING(s) strncpy(sym.tracestring,s,TRACESTRING_SIZE)
#else
#define SYMTRACESTRING(s) // to nothing
#endif

extern int UsedBits;

static void ref_pic_list_reordering();
static void pred_weight_table();


/*!
************************************************************************
* \brief
*    calculate Ceil(Log2(uiVal))
************************************************************************
*/
unsigned CeilLog2( unsigned uiVal)
{
  unsigned uiTmp = uiVal-1;
  unsigned uiRet = 0;

  while( uiTmp != 0 )
  {
    uiTmp >>= 1;
    uiRet++;
  }
  return uiRet;
}


/*!
************************************************************************
* \brief
*    read the first part of the header (only the pic_parameter_set_id)
* \return
*    Length of the first part of the slice header (in bits)
************************************************************************
*/
int FirstPartOfSliceHeader()
{
  Slice *currSlice = img->currentSlice;
  int dP_nr = assignSE2partition[currSlice->dp_mode][SE_HEADER];
  DataPartition *partition = &(currSlice->partArr[dP_nr]);
  Bitstream *currStream = partition->bitstream;
  int tmp;

  UsedBits= partition->bitstream->frame_bitoffset; // was hardcoded to 31 for previous start-code. This is better.

  // Get first_mb_in_slice
  currSlice->start_mb_nr = ue_v ("SH: first_mb_in_slice", currStream);

  tmp = ue_v ("SH: slice_type", currStream);

  if (tmp>4) tmp -=5;

  img->type = currSlice->picture_type = (SliceType) tmp;

  currSlice->pic_parameter_set_id = ue_v ("SH: pic_parameter_set_id", currStream);

  return UsedBits;
}

/*!
************************************************************************
* \brief
*    read the scond part of the header (without the pic_parameter_set_id 
* \return
*    Length of the second part of the Slice header in bits
************************************************************************
*/
int RestOfSliceHeader()
{
  Slice *currSlice = img->currentSlice;
  int dP_nr = assignSE2partition[currSlice->dp_mode][SE_HEADER];
  DataPartition *partition = &(currSlice->partArr[dP_nr]);
  Bitstream *currStream = partition->bitstream;

  int val, len;
#ifdef ADAPTIVE_FILTER
  int sub_pos, i, j;
#ifndef EDAIF2
  extern int FILTER_SIZE;
  extern int FILTER_OFFSET;
#endif
  extern int UseAllSubpelPositions;
  extern int SubpelPositionsPattern;
  extern int DiffQFilterCoef[15][SQR_FILTER];            // differences to be received
  extern int FilterFlag[15];
  extern int POS_EQUATION_NUMBER[15];
#endif

  img->frame_num = u_v (active_sps->log2_max_frame_num_minus4 + 4, "SH: frame_num", currStream);

  /* Tian Dong: frame_num gap processing, if found */
  if (img->idr_flag)
  {
    img->pre_frame_num = img->frame_num;
    // picture error concealment
    img->last_ref_pic_poc = 0;
    assert(img->frame_num == 0);
  }

  if (active_sps->frame_mbs_only_flag)
  {
    img->structure = FRAME;
    img->field_pic_flag=0;
  }
  else
  {
    // field_pic_flag   u(1)
    img->field_pic_flag = u_1("SH: field_pic_flag", currStream);
    if (img->field_pic_flag)
    {
      // bottom_field_flag  u(1)
      img->bottom_field_flag = u_1("SH: bottom_field_flag", currStream);

      img->structure = img->bottom_field_flag ? BOTTOM_FIELD : TOP_FIELD;
    }
    else
    {
      img->structure = FRAME;
      img->bottom_field_flag=0;
    }
  }

  currSlice->structure = (PictureStructure)img->structure;

  img->MbaffFrameFlag=(active_sps->mb_adaptive_frame_field_flag && (img->field_pic_flag==0));

  if (img->structure == FRAME       ) assert (img->field_pic_flag == 0);
  if (img->structure == TOP_FIELD   ) assert (img->field_pic_flag == 1 && img->bottom_field_flag == 0);
  if (img->structure == BOTTOM_FIELD) assert (img->field_pic_flag == 1 && img->bottom_field_flag == 1);

  if (img->idr_flag)
  {
    img->idr_pic_id = ue_v("SH: idr_pic_id", currStream);
  }

  if (active_sps->pic_order_cnt_type == 0)
  {
    img->pic_order_cnt_lsb = u_v(active_sps->log2_max_pic_order_cnt_lsb_minus4 + 4, "SH: pic_order_cnt_lsb", currStream);
    if( active_pps->pic_order_present_flag  ==  1 &&  !img->field_pic_flag )
      img->delta_pic_order_cnt_bottom = se_v("SH: delta_pic_order_cnt_bottom", currStream);
    else
      img->delta_pic_order_cnt_bottom = 0;  
  }
  if( active_sps->pic_order_cnt_type == 1 && !active_sps->delta_pic_order_always_zero_flag ) 
  {
    img->delta_pic_order_cnt[ 0 ] = se_v("SH: delta_pic_order_cnt[0]", currStream);
    if( active_pps->pic_order_present_flag  ==  1  &&  !img->field_pic_flag )
      img->delta_pic_order_cnt[ 1 ] = se_v("SH: delta_pic_order_cnt[1]", currStream);
  }else
  {
    if (active_sps->pic_order_cnt_type == 1)
    {
      img->delta_pic_order_cnt[ 0 ] = 0;
      img->delta_pic_order_cnt[ 1 ] = 0;
    }
  }

  //! redundant_pic_cnt is missing here
  if (active_pps->redundant_pic_cnt_present_flag)
  {
    img->redundant_pic_cnt = ue_v ("SH: redundant_pic_cnt", currStream);
  }

  if(img->type==B_SLICE)
  {
    img->direct_spatial_mv_pred_flag = u_1 ("SH: direct_spatial_mv_pred_flag", currStream);
  }

  img->num_ref_idx_l0_active = active_pps->num_ref_idx_l0_active_minus1 + 1;
  img->num_ref_idx_l1_active = active_pps->num_ref_idx_l1_active_minus1 + 1;

  if(img->type==P_SLICE || img->type == SP_SLICE || img->type==B_SLICE)
  {
    val = u_1 ("SH: num_ref_idx_override_flag", currStream);
    if (val)
    {
      img->num_ref_idx_l0_active = 1 + ue_v ("SH: num_ref_idx_l0_active_minus1", currStream);

      if(img->type==B_SLICE)
      {
        img->num_ref_idx_l1_active = 1 + ue_v ("SH: num_ref_idx_l1_active_minus1", currStream);
      }
    }
  }
  if (img->type!=B_SLICE)
  {
    img->num_ref_idx_l1_active = 0;
  }

  ref_pic_list_reordering();

  img->apply_weights = ((active_pps->weighted_pred_flag && (currSlice->picture_type == P_SLICE || currSlice->picture_type == SP_SLICE) )
    || ((active_pps->weighted_bipred_idc > 0 ) && (currSlice->picture_type == B_SLICE)));

  if ((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
    (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE)))
  {
    pred_weight_table();
  }

  if (img->nal_reference_idc)
    dec_ref_pic_marking(currStream);

  if (active_pps->entropy_coding_mode_flag && img->type!=I_SLICE && img->type!=SI_SLICE)
  {
    img->model_number = ue_v("SH: cabac_init_idc", currStream);
  }
  else 
  {
    img->model_number = 0;
  }

  val = se_v("SH: slice_qp_delta", currStream);
  currSlice->qp = img->qp = 26 + active_pps->pic_init_qp_minus26 + val;
  if ((img->qp < -img->bitdepth_luma_qp_scale) || (img->qp > 51))
    error ("slice_qp_delta makes slice_qp_y out of range", 500);


  currSlice->slice_qp_delta = val;  

  if(img->type==SP_SLICE || img->type == SI_SLICE) 
  {
    if(img->type==SP_SLICE)
    {
      img->sp_switch = u_1 ("SH: sp_for_switch_flag", currStream);
    }
    val = se_v("SH: slice_qs_delta", currStream);
    img->qpsp = 26 + active_pps->pic_init_qs_minus26 + val;
    if ((img->qpsp < 0) || (img->qpsp > 51))
      error ("slice_qs_delta makes slice_qs_y out of range", 500);
  }

  if (active_pps->deblocking_filter_control_present_flag)
  {
    currSlice->LFDisableIdc = ue_v ("SH: disable_deblocking_filter_idc", currStream);

    if (currSlice->LFDisableIdc!=1)
    {
      currSlice->LFAlphaC0Offset = 2 * se_v("SH: slice_alpha_c0_offset_div2", currStream);
      currSlice->LFBetaOffset = 2 * se_v("SH: slice_beta_offset_div2", currStream);
    }
    else
    {
      currSlice->LFAlphaC0Offset = currSlice->LFBetaOffset = 0;
    }
  }
  else 
  {
    currSlice->LFDisableIdc = currSlice->LFAlphaC0Offset = currSlice->LFBetaOffset = 0;
  }
#ifdef ADAPTIVE_QUANTIZATION
  img->slice_fractional_quant_flag = u_1 ("SH: slice_fractional_quant_flag", currStream);
  if(img->slice_fractional_quant_flag) 
  {
    img->slice_mqm_signaling_flag = u_1 ("SH: slice_mqm_signaling_flag", currStream);
    if(img->slice_mqm_signaling_flag)
    {
      img->slice_modeling_qm_param0 = se_v  ("SH: slice_modeling_qm_param0",  currStream) + 32;
      img->slice_modeling_qm_param1 = se_v  ("SH: slice_modeling_qm_param1",  currStream) + 16;
    }
    else
    {
      img->slice_modeling_qm_param0 = 32;
      img->slice_modeling_qm_param1 = 16;
    }

    SetAutoScalingListDataForMQM(img->slice_modeling_qm_param0, img->slice_modeling_qm_param1);
    CreateModelingMatrixForIAQMS();
  }
#endif
#ifdef EIGHTH_PEL
  // 1/8-pel motion vector resolution (munderl)
  img->mv_res = ue_v( "SH: motion vector resolution" , currStream);
#endif
#ifdef ADAPTIVE_FD_SD_CODING
  img->Allow_SD_Coding = u_1( "SH: SD coding" , currStream);
  if (img->Allow_SD_Coding)
  {
    img->SD_Quantizer  = u_1( "SH: Quantizer" , currStream);
  }
#endif

  if (active_pps->num_slice_groups_minus1>0 && active_pps->slice_group_map_type>=3 &&
    active_pps->slice_group_map_type<=5)
  {
    len = (active_sps->pic_height_in_map_units_minus1+1)*(active_sps->pic_width_in_mbs_minus1+1)/ 
      (active_pps->slice_group_change_rate_minus1+1);
    if (((active_sps->pic_height_in_map_units_minus1+1)*(active_sps->pic_width_in_mbs_minus1+1))% 
      (active_pps->slice_group_change_rate_minus1+1))
      len +=1;

    len = CeilLog2(len+1);

    img->slice_group_change_cycle = u_v (len, "SH: slice_group_change_cycle", currStream);
  }
  img->PicHeightInMbs = img->FrameHeightInMbs / ( 1 + img->field_pic_flag );
  img->PicSizeInMbs   = img->PicWidthInMbs * img->PicHeightInMbs;
  img->FrameSizeInMbs = img->PicWidthInMbs * img->FrameHeightInMbs;
#ifdef ADAPTIVE_FILTER
  ResetAdaptiveFilter();

  if(UseAdaptiveFilterForCurrentFrame())
  {
#ifdef E_DAIF
    img->AdaptiveFilterFlag = u_v(3,"apply_adaptive_filter", currStream);
#else
    img->AdaptiveFilterFlag = u_v(2,"apply_adaptive_filter", currStream);
#endif
    if(img->AdaptiveFilterFlag == FILTER_TYPE_2D_NS) // 2D aif
    {
      UseAllSubpelPositions = u_1("use_all_subpel_positions", currStream);
      if(UseAllSubpelPositions)
      {
        for(sub_pos=a_pos; sub_pos <= o_pos; sub_pos++)
          FilterFlag[sub_pos] = 1;
        SubpelPositionsPattern =   (FilterFlag[j_pos]<<4) + (FilterFlag[f_pos]<<3) + (FilterFlag[e_pos]<<2) + 
          (FilterFlag[b_pos]<<1) +  FilterFlag[a_pos];
      }
      else
      {
        SubpelPositionsPattern = ue_v("positions_pattern", currStream);
        FilterFlag[a_pos] = 
          FilterFlag[c_pos] = 
          FilterFlag[d_pos] =
          FilterFlag[l_pos] =  SubpelPositionsPattern       & 1;
        FilterFlag[b_pos] = 
          FilterFlag[h_pos] = (SubpelPositionsPattern >> 1) & 1;
        FilterFlag[e_pos] = 
          FilterFlag[g_pos] = 
          FilterFlag[m_pos] =
          FilterFlag[o_pos] = (SubpelPositionsPattern >> 2) & 1;
        FilterFlag[f_pos] = 
          FilterFlag[i_pos] = 
          FilterFlag[k_pos] =
          FilterFlag[n_pos] = (SubpelPositionsPattern >> 3) & 1;
        FilterFlag[j_pos] = (SubpelPositionsPattern >> 4) & 1;
      }
      if(FilterFlag[a_pos])
        for(i= 0; i < POS_EQUATION_NUMBER[a_pos]; i++)
          DiffQFilterCoef[a_pos][FILTER_SIZE*FILTER_OFFSET+i] = se_v("a_pos", currStream);
      if(FilterFlag[b_pos])
        for(i= 0; i < POS_EQUATION_NUMBER[b_pos]; i++)
          DiffQFilterCoef[b_pos][FILTER_SIZE*FILTER_OFFSET+i] = se_v("b_pos", currStream);
      if(FilterFlag[e_pos])
        for(i = 0; i < FILTER_SIZE; i++)
          for(j = i; j < FILTER_SIZE; j++)
            DiffQFilterCoef[e_pos][FILTER_SIZE*i+j] = se_v("e_pos", currStream);
      if(FilterFlag[f_pos])
        for(i = 0; i < FILTER_SIZE; i++)
          for(j = 0; j < FILTER_SIZE/2; j++)
            DiffQFilterCoef[f_pos][FILTER_SIZE*i+j] = se_v("f_pos", currStream);
      if(FilterFlag[j_pos])
        for(i = 0; i < FILTER_SIZE/2; i++)
          for(j = i; j < FILTER_SIZE/2; j++)
            DiffQFilterCoef[j_pos][FILTER_SIZE*i+j] = se_v("j_pos", currStream);
      for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
        ExtendFilterCoefficientsInt(sub_pos, DiffQFilterCoef);
    }
    else if (img->AdaptiveFilterFlag == FILTER_TYPE_2D_S) // separable aif
    {
      UseAllSubpelPositions = u_1("use_all_subpel_positions", currStream);
      if(UseAllSubpelPositions)
      {
        for(sub_pos=a_pos; sub_pos <= o_pos; sub_pos++)
          FilterFlag[sub_pos] = 1;
        SubpelPositionsPattern = 2047;
      }
      else
      {
        SubpelPositionsPattern = ue_v ("positions_pattern", currStream);
        FilterFlag[a_pos] = SubpelPositionsPattern       & 1;
        FilterFlag[b_pos] = (SubpelPositionsPattern >> 1)    & 1; 
        FilterFlag[c_pos] = (SubpelPositionsPattern >> 2)    & 1; 
        FilterFlag[d_pos] = (SubpelPositionsPattern >> 3)    & 1; 
        FilterFlag[e_pos] = (SubpelPositionsPattern >> 4)    & 1; 
        FilterFlag[f_pos] = (SubpelPositionsPattern >> 5)    & 1; 
        FilterFlag[g_pos] = (SubpelPositionsPattern >> 6)    & 1; 
        FilterFlag[h_pos] = (SubpelPositionsPattern >> 7)    & 1;  
        FilterFlag[i_pos] = (SubpelPositionsPattern >> 8)    & 1; 
        FilterFlag[j_pos] = (SubpelPositionsPattern >> 9)    & 1; 
        FilterFlag[k_pos] = (SubpelPositionsPattern >> 10)  & 1;  
        FilterFlag[l_pos] = FilterFlag[d_pos];  
        FilterFlag[m_pos] = FilterFlag[e_pos]; 
        FilterFlag[n_pos] = FilterFlag[f_pos]; 
        FilterFlag[o_pos] = FilterFlag[g_pos]; 
      }

      // a_pos
      if(FilterFlag[a_pos])
        for(i= 0; i < FILTER_SIZE; i++)
          DiffQFilterCoef[a_pos][i] = se_v("a_pos", currStream);

      // b_pos
      if(FilterFlag[b_pos])
      {
        for(i= 0; i < FILTER_SIZE >> 1; i++)
          DiffQFilterCoef[b_pos][i] = se_v("b_pos", currStream);
        DiffQFilterCoef[b_pos][3] = DiffQFilterCoef[b_pos][2];
        DiffQFilterCoef[b_pos][4] = DiffQFilterCoef[b_pos][1];                 
        DiffQFilterCoef[b_pos][5] = DiffQFilterCoef[b_pos][0]; 
      }

      // c_pos
      if(FilterFlag[c_pos])
      {
        if (FilterFlag[a_pos])
          for(i= 0; i < FILTER_SIZE; i++)
            DiffQFilterCoef[c_pos][i] = se_v("c_pos", currStream) + DiffQFilterCoef[a_pos][FILTER_SIZE-i-1];
        else
          for(i= 0; i < FILTER_SIZE; i++)
            DiffQFilterCoef[c_pos][i] = se_v("c_pos", currStream);
      }

      // d_pos
      if(FilterFlag[d_pos])
        for(i= 0; i < FILTER_SIZE; i++)
          DiffQFilterCoef[d_pos][i] = se_v("d_pos", currStream);

      // e_pos
      if(FilterFlag[e_pos])
      {
        if (FilterFlag[d_pos])
          for(i= 0; i < FILTER_SIZE; i++)
            DiffQFilterCoef[e_pos][i] = se_v("e_pos", currStream) + DiffQFilterCoef[d_pos][i];
        else
          for(i= 0; i < FILTER_SIZE; i++)
            DiffQFilterCoef[e_pos][i] = se_v("e_pos", currStream);
      }

      // f_pos
      if(FilterFlag[f_pos])
      {
        if (FilterFlag[e_pos])
          for(i= 0; i < FILTER_SIZE; i++)
            DiffQFilterCoef[f_pos][i] = se_v("f_pos", currStream) + DiffQFilterCoef[e_pos][i];
        else
          for(i= 0; i < FILTER_SIZE; i++)
            DiffQFilterCoef[f_pos][i] = se_v("f_pos", currStream);
      }

      // g_pos
      if(FilterFlag[g_pos])
      {
        if (FilterFlag[f_pos])
          for(i= 0; i < FILTER_SIZE; i++)
            DiffQFilterCoef[g_pos][i] = se_v("g_pos", currStream) + DiffQFilterCoef[f_pos][i];
        else
          for(i= 0; i < FILTER_SIZE; i++)
            DiffQFilterCoef[g_pos][i] = se_v("g_pos", currStream);
      }

      // h_pos
      if(FilterFlag[h_pos])
      {
        for(i= 0; i < FILTER_SIZE>>1; i++)
          DiffQFilterCoef[h_pos][i] = se_v("h_pos", currStream);
        DiffQFilterCoef[h_pos][3] = DiffQFilterCoef[h_pos][2];
        DiffQFilterCoef[h_pos][4] = DiffQFilterCoef[h_pos][1];                 
        DiffQFilterCoef[h_pos][5] = DiffQFilterCoef[h_pos][0]; 
      }

      // i_pos
      if(FilterFlag[i_pos])
      {
        if(FilterFlag[h_pos])
          for(i= 0; i < FILTER_SIZE>>1; i++)
            DiffQFilterCoef[i_pos][i] = se_v("i_pos", currStream) + DiffQFilterCoef[h_pos][i];
        else
          for(i= 0; i < FILTER_SIZE>>1; i++)
            DiffQFilterCoef[i_pos][i] = se_v("i_pos", currStream);
        DiffQFilterCoef[i_pos][3] = DiffQFilterCoef[i_pos][2];
        DiffQFilterCoef[i_pos][4] = DiffQFilterCoef[i_pos][1];                 
        DiffQFilterCoef[i_pos][5] = DiffQFilterCoef[i_pos][0]; 
      }

      // j_pos
      if(FilterFlag[j_pos])
      {
        if(FilterFlag[i_pos])
          for(i= 0; i < FILTER_SIZE>>1; i++)
            DiffQFilterCoef[j_pos][i] = se_v("j_pos", currStream) + DiffQFilterCoef[i_pos][i];
        else
          for(i= 0; i < FILTER_SIZE>>1; i++)
            DiffQFilterCoef[j_pos][i] = se_v("j_pos", currStream);
        DiffQFilterCoef[j_pos][3] = DiffQFilterCoef[j_pos][2];
        DiffQFilterCoef[j_pos][4] = DiffQFilterCoef[j_pos][1];                 
        DiffQFilterCoef[j_pos][5] = DiffQFilterCoef[j_pos][0]; 
      }

      // k_pos
      if(FilterFlag[k_pos])
      {
        if(FilterFlag[j_pos])
          for(i= 0; i < FILTER_SIZE>>1; i++)
            DiffQFilterCoef[k_pos][i] = se_v("k_pos", currStream) + DiffQFilterCoef[j_pos][i];
        else
          for(i= 0; i < FILTER_SIZE>>1; i++)
            DiffQFilterCoef[k_pos][i] = se_v("k_pos", currStream);
        DiffQFilterCoef[k_pos][3] = DiffQFilterCoef[k_pos][2];
        DiffQFilterCoef[k_pos][4] = DiffQFilterCoef[k_pos][1];                 
        DiffQFilterCoef[k_pos][5] = DiffQFilterCoef[k_pos][0]; 
      }

      // l_pos
      if(FilterFlag[l_pos])
        for(i= 0; i < FILTER_SIZE; i++)
          DiffQFilterCoef[l_pos][i] = DiffQFilterCoef[d_pos][FILTER_SIZE-i-1];

      // m_pos
      if(FilterFlag[m_pos])
        for(i= 0; i < FILTER_SIZE; i++)
          DiffQFilterCoef[m_pos][i] = DiffQFilterCoef[e_pos][FILTER_SIZE-i-1];

      // n_pos
      if(FilterFlag[n_pos])
        for(i= 0; i < FILTER_SIZE; i++)
          DiffQFilterCoef[n_pos][i] = DiffQFilterCoef[f_pos][FILTER_SIZE-i-1];

      // o_pos
      if(FilterFlag[o_pos])
        for(i= 0; i < FILTER_SIZE; i++)
          DiffQFilterCoef[o_pos][i] = DiffQFilterCoef[g_pos][FILTER_SIZE-i-1];      
    }
#ifdef DIRECTIONAL_FILTER
    else if (img->AdaptiveFilterFlag == FILTER_TYPE_1D) // 1D-AIF
    {
      img->ImpType = u_v(1,"implem_type", currStream);
      initFilterCustom(img->AdaptiveFilterFlag);
      readFilterCoefs(img->AdaptiveFilterFlag, currStream);
    }
#endif  // DIRECTIONAL_FILTER
#ifdef E_DAIF
    else if ((img->AdaptiveFilterFlag == FILTER_TYPE_EDAIF)) 
    {
      initFilterCustom(img->AdaptiveFilterFlag);
      readFilterCoefs_EDAIF(img->AdaptiveFilterFlag, currStream);
    }
#endif  // E_DAIF
#ifdef EAIF
    else if ((img->AdaptiveFilterFlag == FILTER_TYPE_EAIF)) 
    {
      initFilterCustom(img->AdaptiveFilterFlag);
      readFilterCoefs_EAIF(img->AdaptiveFilterFlag, currStream);
    }
#endif
  }
  else
  {
    img->AdaptiveFilterFlag = -1;
  }
#endif
#ifdef USE_INTRA_MDDT
  img->UseIntraMDDT = u_v(1, "use_intra_MDDT", currStream);
#endif

#ifdef USE_HP_FILTER    
#ifdef SWITCHED_FILTERS
  img->use_high_precision_flag = u_v(3, "use_high_precision_flag" , currStream);
#else
  img->use_high_precision_flag = u_v(2,"use_high_precision_flag" , currStream);
#endif
#endif

#ifdef ADAPTIVE_LOOP_FILTER
#ifdef SWITCHED_FILTERS
  if((img->use_high_precision_flag != HPF_SIFO) && (img->use_high_precision_flag != HPF_SIFO_FPO))
  {
#endif
  active_alfps->alf_flag = u_1("SH: adaptive_loop_filter_flag", currStream);

  // padding end of slice header for ALF
  if ( currStream->frame_bitoffset%8 != 0 )
  {
    while(currStream->frame_bitoffset%8 != 0)
    {
      u_1("SH: padding bits for ALF", currStream);
    }
  }

  // adaptive loopfilter parameters
  if(active_alfps->alf_flag)
  {
    read_ALFParameterSet(active_alfps, partition);
  }
#ifdef SWITCHED_FILTERS
  }
#endif
#endif

  return UsedBits;
}


/*!
************************************************************************
* \brief
*    read the reference picture reordering information
************************************************************************
*/
static void ref_pic_list_reordering()
{
  Slice *currSlice = img->currentSlice;
  int dP_nr = assignSE2partition[currSlice->dp_mode][SE_HEADER];
  DataPartition *partition = &(currSlice->partArr[dP_nr]);
  Bitstream *currStream = partition->bitstream;
  int i, val;

  alloc_ref_pic_list_reordering_buffer(currSlice);

  if (img->type!=I_SLICE && img->type!=SI_SLICE)
  {
    val = currSlice->ref_pic_list_reordering_flag_l0 = u_1 ("SH: ref_pic_list_reordering_flag_l0", currStream);

    if (val)
    {
      i=0;
      do
      {
        val = currSlice->reordering_of_pic_nums_idc_l0[i] = ue_v("SH: reordering_of_pic_nums_idc_l0", currStream);
        if (val==0 || val==1)
        {
          currSlice->abs_diff_pic_num_minus1_l0[i] = ue_v("SH: abs_diff_pic_num_minus1_l0", currStream);
        }
        else
        {
          if (val==2)
          {
            currSlice->long_term_pic_idx_l0[i] = ue_v("SH: long_term_pic_idx_l0", currStream);
          }
        }
        i++;
        // assert (i>img->num_ref_idx_l0_active);
      } while (val != 3);
    }
  }

  if (img->type==B_SLICE)
  {
    val = currSlice->ref_pic_list_reordering_flag_l1 = u_1 ("SH: ref_pic_list_reordering_flag_l1", currStream);

    if (val)
    {
      i=0;
      do
      {
        val = currSlice->reordering_of_pic_nums_idc_l1[i] = ue_v("SH: reordering_of_pic_nums_idc_l1", currStream);
        if (val==0 || val==1)
        {
          currSlice->abs_diff_pic_num_minus1_l1[i] = ue_v("SH: abs_diff_pic_num_minus1_l1", currStream);
        }
        else
        {
          if (val==2)
          {
            currSlice->long_term_pic_idx_l1[i] = ue_v("SH: long_term_pic_idx_l1", currStream);
          }
        }
        i++;
        // assert (i>img->num_ref_idx_l1_active);
      } while (val != 3);
    }
  }

  // set reference index of redundant slices.
  if(img->redundant_pic_cnt)
  {
    redundant_slice_ref_idx = currSlice->abs_diff_pic_num_minus1_l0[0] + 1;
  }
}

/*!
************************************************************************
* \brief
*    read the weighted prediction tables
************************************************************************
*/
static void pred_weight_table()
{
  Slice *currSlice = img->currentSlice;
  int dP_nr = assignSE2partition[currSlice->dp_mode][SE_HEADER];
  DataPartition *partition = &(currSlice->partArr[dP_nr]);
  Bitstream *currStream = partition->bitstream;
  int luma_weight_flag_l0, luma_weight_flag_l1, chroma_weight_flag_l0, chroma_weight_flag_l1;
  int i,j;

  img->luma_log2_weight_denom = ue_v ("SH: luma_log2_weight_denom", currStream);
  img->wp_round_luma = img->luma_log2_weight_denom ? 1<<(img->luma_log2_weight_denom - 1): 0;

  if ( 0 != active_sps->chroma_format_idc)
  {
    img->chroma_log2_weight_denom = ue_v ("SH: chroma_log2_weight_denom", currStream);
    img->wp_round_chroma = img->chroma_log2_weight_denom ? 1<<(img->chroma_log2_weight_denom - 1): 0;
  }

  reset_wp_params(img);

  for (i=0; i<img->num_ref_idx_l0_active; i++)
  {
    luma_weight_flag_l0 = u_1("SH: luma_weight_flag_l0", currStream);

    if (luma_weight_flag_l0)
    {
      img->wp_weight[0][i][0] = se_v ("SH: luma_weight_l0", currStream);
      img->wp_offset[0][i][0] = se_v ("SH: luma_offset_l0", currStream);
#ifdef  BUG_FIX_FOR_FREXT
      img->wp_offset[0][i][0] = img->wp_offset[0][i][0]<<(img->bitdepth_luma-8);
#endif
    }
    else
    {
      img->wp_weight[0][i][0] = 1<<img->luma_log2_weight_denom;
      img->wp_offset[0][i][0] = 0;
    }

    if (active_sps->chroma_format_idc != 0)
    {
      chroma_weight_flag_l0 = u_1 ("SH: chroma_weight_flag_l0", currStream);

      for (j=1; j<3; j++)
      {
        if (chroma_weight_flag_l0)
        {
          img->wp_weight[0][i][j] = se_v("SH: chroma_weight_l0", currStream);
          img->wp_offset[0][i][j] = se_v("SH: chroma_offset_l0", currStream);
#ifdef  BUG_FIX_FOR_FREXT
          img->wp_offset[0][i][j] = img->wp_offset[0][i][j]<<(img->bitdepth_chroma-8);
#endif
        }
        else
        {
          img->wp_weight[0][i][j] = 1<<img->chroma_log2_weight_denom;
          img->wp_offset[0][i][j] = 0;
        }
      }
    }
  }
  if ((img->type == B_SLICE) && active_pps->weighted_bipred_idc == 1)
  {
    for (i=0; i<img->num_ref_idx_l1_active; i++)
    {
      luma_weight_flag_l1 = u_1("SH: luma_weight_flag_l1", currStream);

      if (luma_weight_flag_l1)
      {
        img->wp_weight[1][i][0] = se_v ("SH: luma_weight_l1", currStream);
        img->wp_offset[1][i][0] = se_v ("SH: luma_offset_l1", currStream);
#ifdef  BUG_FIX_FOR_FREXT
        img->wp_offset[1][i][0] = img->wp_offset[1][i][0]<<(img->bitdepth_luma-8);
#endif
      }
      else
      {
        img->wp_weight[1][i][0] = 1<<img->luma_log2_weight_denom;
        img->wp_offset[1][i][0] = 0;
      }

      if (active_sps->chroma_format_idc != 0)
      {
        chroma_weight_flag_l1 = u_1 ("SH: chroma_weight_flag_l1", currStream);

        for (j=1; j<3; j++)
        {
          if (chroma_weight_flag_l1)
          {
            img->wp_weight[1][i][j] = se_v("SH: chroma_weight_l1", currStream);
            img->wp_offset[1][i][j] = se_v("SH: chroma_offset_l1", currStream);
#ifdef  BUG_FIX_FOR_FREXT
            img->wp_offset[1][i][j] = img->wp_offset[1][i][j]<<(img->bitdepth_chroma-8);
#endif
          }
          else
          {
            img->wp_weight[1][i][j] = 1<<img->chroma_log2_weight_denom;
            img->wp_offset[1][i][j] = 0;
          }
        }
      }
    }  
  }
}


/*!
************************************************************************
* \brief
*    read the memory control operations
************************************************************************
*/
void dec_ref_pic_marking(Bitstream *currStream)
{
  int val;

  DecRefPicMarking_t *tmp_drpm,*tmp_drpm2;

  // free old buffer content
  while (img->dec_ref_pic_marking_buffer)
  { 
    tmp_drpm=img->dec_ref_pic_marking_buffer;

    img->dec_ref_pic_marking_buffer=tmp_drpm->Next;
    free (tmp_drpm);
  } 

  if (img->idr_flag)
  {
    img->no_output_of_prior_pics_flag = u_1("SH: no_output_of_prior_pics_flag", currStream);
    img->long_term_reference_flag = u_1("SH: long_term_reference_flag", currStream);
  }
  else
  {
    img->adaptive_ref_pic_buffering_flag = u_1("SH: adaptive_ref_pic_buffering_flag", currStream);
    if (img->adaptive_ref_pic_buffering_flag)
    {
      // read Memory Management Control Operation 
      do
      {
        tmp_drpm=(DecRefPicMarking_t*)calloc (1,sizeof (DecRefPicMarking_t));
        tmp_drpm->Next=NULL;

        val = tmp_drpm->memory_management_control_operation = ue_v("SH: memory_management_control_operation", currStream);

        if ((val==1)||(val==3)) 
        {
          tmp_drpm->difference_of_pic_nums_minus1 = ue_v("SH: difference_of_pic_nums_minus1", currStream);
        }
        if (val==2)
        {
          tmp_drpm->long_term_pic_num = ue_v("SH: long_term_pic_num", currStream);
        }

        if ((val==3)||(val==6))
        {
          tmp_drpm->long_term_frame_idx = ue_v("SH: long_term_frame_idx", currStream);
        }
        if (val==4)
        {
          tmp_drpm->max_long_term_frame_idx_plus1 = ue_v("SH: max_long_term_pic_idx_plus1", currStream);
        }

        // add command
        if (img->dec_ref_pic_marking_buffer==NULL) 
        {
          img->dec_ref_pic_marking_buffer=tmp_drpm;
        }
        else
        {
          tmp_drpm2=img->dec_ref_pic_marking_buffer;
          while (tmp_drpm2->Next!=NULL) tmp_drpm2=tmp_drpm2->Next;
          tmp_drpm2->Next=tmp_drpm;
        }

      }while (val != 0);

    }
  }
}

/*!
************************************************************************
* \brief
*    To calculate the poc values
*        based upon JVT-F100d2
*  POC200301: Until Jan 2003, this function will calculate the correct POC
*    values, but the management of POCs in buffered pictures may need more work.
* \return
*    none
************************************************************************
*/
void decode_poc(struct img_par *img)
{
  int i;
  // for POC mode 0:
  unsigned int MaxPicOrderCntLsb = (1<<(active_sps->log2_max_pic_order_cnt_lsb_minus4+4));

  switch ( active_sps->pic_order_cnt_type )
  {
  case 0: // POC MODE 0
    // 1st
    if(img->idr_flag)
    {
      img->PrevPicOrderCntMsb = 0;
      img->PrevPicOrderCntLsb = 0;
    }
    else
    {
      if (img->last_has_mmco_5) 
      {
        if (img->last_pic_bottom_field)
        {
          img->PrevPicOrderCntMsb = 0;
          img->PrevPicOrderCntLsb = 0;
        }
        else
        {
          img->PrevPicOrderCntMsb = 0;
          img->PrevPicOrderCntLsb = img->toppoc;
        }
      }
    }
    // Calculate the MSBs of current picture
    if( img->pic_order_cnt_lsb  <  img->PrevPicOrderCntLsb  &&  
      ( img->PrevPicOrderCntLsb - img->pic_order_cnt_lsb )  >=  ( MaxPicOrderCntLsb / 2 ) )
      img->PicOrderCntMsb = img->PrevPicOrderCntMsb + MaxPicOrderCntLsb;
    else if ( img->pic_order_cnt_lsb  >  img->PrevPicOrderCntLsb  &&
      ( img->pic_order_cnt_lsb - img->PrevPicOrderCntLsb )  >  ( MaxPicOrderCntLsb / 2 ) )
      img->PicOrderCntMsb = img->PrevPicOrderCntMsb - MaxPicOrderCntLsb;
    else
      img->PicOrderCntMsb = img->PrevPicOrderCntMsb;

    // 2nd

    if(img->field_pic_flag==0)
    {           //frame pix
      img->toppoc = img->PicOrderCntMsb + img->pic_order_cnt_lsb;
      img->bottompoc = img->toppoc + img->delta_pic_order_cnt_bottom;
      img->ThisPOC = img->framepoc = (img->toppoc < img->bottompoc)? img->toppoc : img->bottompoc; // POC200301
    }
    else if (img->bottom_field_flag==0)
    {  //top field 
      img->ThisPOC= img->toppoc = img->PicOrderCntMsb + img->pic_order_cnt_lsb;
    }
    else
    {  //bottom field
      img->ThisPOC= img->bottompoc = img->PicOrderCntMsb + img->pic_order_cnt_lsb;
    }
    img->framepoc=img->ThisPOC;

    if ( img->frame_num!=img->PreviousFrameNum)
      img->PreviousFrameNum=img->frame_num;

    if(img->nal_reference_idc)
    {
      img->PrevPicOrderCntLsb = img->pic_order_cnt_lsb;
      img->PrevPicOrderCntMsb = img->PicOrderCntMsb;
    }

    break;

  case 1: // POC MODE 1
    // 1st
    if(img->idr_flag)
    {
      img->FrameNumOffset=0;     //  first pix of IDRGOP, 
      img->delta_pic_order_cnt[0]=0;                        //ignore first delta
      if(img->frame_num)
        error("frame_num not equal to zero in IDR picture", -1020);
    }
    else 
    {
      if (img->last_has_mmco_5)
      {
        img->PreviousFrameNumOffset = 0;
        img->PreviousFrameNum = 0;
      }
      if (img->frame_num<img->PreviousFrameNum)
      {             //not first pix of IDRGOP
        img->FrameNumOffset = img->PreviousFrameNumOffset + img->MaxFrameNum;
      }
      else 
      {
        img->FrameNumOffset = img->PreviousFrameNumOffset;
      }
    }

    // 2nd
    if(active_sps->num_ref_frames_in_pic_order_cnt_cycle) 
      img->AbsFrameNum = img->FrameNumOffset+img->frame_num;
    else 
      img->AbsFrameNum=0;
    if( (!img->nal_reference_idc) && img->AbsFrameNum>0)
      img->AbsFrameNum--;

    // 3rd
    img->ExpectedDeltaPerPicOrderCntCycle=0;

    if(active_sps->num_ref_frames_in_pic_order_cnt_cycle)
      for(i=0;i<(int) active_sps->num_ref_frames_in_pic_order_cnt_cycle;i++)
        img->ExpectedDeltaPerPicOrderCntCycle += active_sps->offset_for_ref_frame[i];

    if(img->AbsFrameNum)
    {
      img->PicOrderCntCycleCnt = (img->AbsFrameNum-1)/active_sps->num_ref_frames_in_pic_order_cnt_cycle;
      img->FrameNumInPicOrderCntCycle = (img->AbsFrameNum-1)%active_sps->num_ref_frames_in_pic_order_cnt_cycle;
      img->ExpectedPicOrderCnt = img->PicOrderCntCycleCnt*img->ExpectedDeltaPerPicOrderCntCycle;
      for(i=0;i<=(int)img->FrameNumInPicOrderCntCycle;i++)
        img->ExpectedPicOrderCnt += active_sps->offset_for_ref_frame[i];
    }
    else 
      img->ExpectedPicOrderCnt=0;

    if(!img->nal_reference_idc)
      img->ExpectedPicOrderCnt += active_sps->offset_for_non_ref_pic;

    if(img->field_pic_flag==0)
    {           //frame pix
      img->toppoc = img->ExpectedPicOrderCnt + img->delta_pic_order_cnt[0];
      img->bottompoc = img->toppoc + active_sps->offset_for_top_to_bottom_field + img->delta_pic_order_cnt[1];
      img->ThisPOC = img->framepoc = (img->toppoc < img->bottompoc)? img->toppoc : img->bottompoc; // POC200301
    }
    else if (img->bottom_field_flag==0)
    {  //top field 
      img->ThisPOC = img->toppoc = img->ExpectedPicOrderCnt + img->delta_pic_order_cnt[0];
    } 
    else
    {  //bottom field
      img->ThisPOC = img->bottompoc = img->ExpectedPicOrderCnt + active_sps->offset_for_top_to_bottom_field + img->delta_pic_order_cnt[0];
    }
    img->framepoc=img->ThisPOC;

    img->PreviousFrameNum=img->frame_num;
    img->PreviousFrameNumOffset=img->FrameNumOffset;

    break;


  case 2: // POC MODE 2
    if(img->idr_flag) // IDR picture
    {
      img->FrameNumOffset=0;     //  first pix of IDRGOP, 
      img->ThisPOC = img->framepoc = img->toppoc = img->bottompoc = 0;
      if(img->frame_num) 
        error("frame_num not equal to zero in IDR picture", -1020);
    }
    else
    {
      if (img->last_has_mmco_5)
      {
        img->PreviousFrameNum = 0;
        img->PreviousFrameNumOffset = 0;
      }
      if (img->frame_num<img->PreviousFrameNum)
        img->FrameNumOffset = img->PreviousFrameNumOffset + img->MaxFrameNum;
      else 
        img->FrameNumOffset = img->PreviousFrameNumOffset;


      img->AbsFrameNum = img->FrameNumOffset+img->frame_num;
      if(!img->nal_reference_idc)
        img->ThisPOC = (2*img->AbsFrameNum - 1);
      else
        img->ThisPOC = (2*img->AbsFrameNum);

      if (img->field_pic_flag==0)
        img->toppoc = img->bottompoc = img->framepoc = img->ThisPOC;
      else if (img->bottom_field_flag==0)
        img->toppoc = img->framepoc = img->ThisPOC;
      else img->bottompoc = img->framepoc = img->ThisPOC;
    }

    img->PreviousFrameNum=img->frame_num;
    img->PreviousFrameNumOffset=img->FrameNumOffset;
    break;


  default:
    //error must occurs
    assert( 1==0 );
    break;
  }
}

/*!
************************************************************************
* \brief
*    A little helper for the debugging of POC code
* \return
*    none
************************************************************************
*/
int dumppoc(struct img_par *img) 
{
  printf ("\nPOC locals...\n");
  printf ("toppoc                                %d\n", img->toppoc);
  printf ("bottompoc                             %d\n", img->bottompoc);
  printf ("frame_num                             %d\n", img->frame_num);
  printf ("field_pic_flag                        %d\n", img->field_pic_flag);
  printf ("bottom_field_flag                     %d\n", img->bottom_field_flag);
  printf ("POC SPS\n");
  printf ("log2_max_frame_num_minus4             %d\n", active_sps->log2_max_frame_num_minus4);         // POC200301
  printf ("log2_max_pic_order_cnt_lsb_minus4     %d\n", active_sps->log2_max_pic_order_cnt_lsb_minus4);
  printf ("pic_order_cnt_type                    %d\n", active_sps->pic_order_cnt_type);
  printf ("num_ref_frames_in_pic_order_cnt_cycle %d\n", active_sps->num_ref_frames_in_pic_order_cnt_cycle);
  printf ("delta_pic_order_always_zero_flag      %d\n", active_sps->delta_pic_order_always_zero_flag);
  printf ("offset_for_non_ref_pic                %d\n", active_sps->offset_for_non_ref_pic);
  printf ("offset_for_top_to_bottom_field        %d\n", active_sps->offset_for_top_to_bottom_field);
  printf ("offset_for_ref_frame[0]               %d\n", active_sps->offset_for_ref_frame[0]);
  printf ("offset_for_ref_frame[1]               %d\n", active_sps->offset_for_ref_frame[1]);
  printf ("POC in SLice Header\n");
  printf ("pic_order_present_flag                %d\n", active_pps->pic_order_present_flag);
  printf ("delta_pic_order_cnt[0]                %d\n", img->delta_pic_order_cnt[0]);
  printf ("delta_pic_order_cnt[1]                %d\n", img->delta_pic_order_cnt[1]);
  printf ("delta_pic_order_cnt[2]                %d\n", img->delta_pic_order_cnt[2]);
  printf ("idr_flag                              %d\n", img->idr_flag);
  printf ("MaxFrameNum                           %d\n", img->MaxFrameNum);

  return 0;
}

/*!
************************************************************************
* \brief
*    return the poc of img as per (8-1) JVT-F100d2
*  POC200301
************************************************************************
*/
int picture_order(struct img_par *img)
{
  if (img->field_pic_flag==0) // is a frame
    return img->framepoc;
  else if (img->bottom_field_flag==0) // top field
    return img->toppoc;
  else // bottom field
    return img->bottompoc;
}

#ifdef E_DAIF
int readCoeff1(char *tracestring, int numQBits, Bitstream *bitstream)
{
  int prefix, suffix, suffixLen; 
  int coeffQ; 
  int sign = 0;
  int mag; 
  int rangeQ = 1<<(numQBits-1); 

  int i, bit; 

  // read prefix code
  for(i = 0; i < 5; i++)
  {
    bit = u_1(tracestring, bitstream);
    if(!bit)
      break; 
  }
  prefix = i; 

  // read suffix code
  suffixLen = numQBits-2-(5-prefix);
  if(prefix == 0)
    suffixLen++;
  suffix = u_v(suffixLen, tracestring, bitstream); 

  // get the magnitude 
  if(prefix == 0)
    mag = 0;
  else
    mag = rangeQ>>(6-prefix); 
  mag += suffix; 

  // read sign bit 
  if(mag) 
    sign = u_1(tracestring, bitstream); 

  coeffQ = sign ? -mag:mag; 

  return coeffQ; 
}

void readAIFInteger(Bitstream *currStream)
{
  int i;

  extern int DiffQFilterCoeffInt[SQR_FILTER_INT]; 
  extern int numQBitsInt[SQR_FILTER_INT-1];

  for(i = 0; i < FILTER_SIZE_INT*FILTER_SIZE_INT; i++)
  {
    DiffQFilterCoeffInt[i] = readCoeff1("SH: full-pel filter", numQBitsInt[i], currStream);
  }
}

void readAIFOffset(int sub_pos, Bitstream *bitstream)
{
  int offsetI, offsetF; 
  int sign = 0; 
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

  extern int DiffQFilterOffsetI[15], DiffQFilterOffsetF[15], DiffQFilterOffsetSign[15];
  extern int DiffQFilterOffsetIntI, DiffQFilterOffsetIntF, DiffQFilterOffsetIntSign; 

  if(sub_pos == -1) // full-pel offset
  {
    // integer part
    offsetI = se_v("full-pel offset int", bitstream);

    if(offsetFracCodeLen[offsetI])
      offsetF = u_v(offsetFracCodeLen[offsetI], "full-pel offset frac", bitstream); 
    else
      offsetF = 0; 
    if(offsetI || (offsetF > 0)) 
      sign = u_1("full-pel offset sign", bitstream);
    DiffQFilterOffsetIntI = offsetI; 
    DiffQFilterOffsetIntF = offsetF; 
    DiffQFilterOffsetIntSign = sign ? -1:1; 
  }
  else // subpel offset 
  {
    // send integer, using EXP-GOLOMB
    offsetI = se_v("sub-pel offset int", bitstream);

    // send fraction, using fixed-length code 
    if(offsetFracCodeLen[offsetI])
    {
      offsetF = u_v(offsetFracCodeLen[offsetI], "sub-pel offset frac", bitstream); 
    }
    else 
      offsetF = 0;
    if(offsetI || (offsetF > 0)) 
      sign = u_1("sub-pel offset sign", bitstream);
    DiffQFilterOffsetI[sub_pos] = offsetI; 
    DiffQFilterOffsetF[sub_pos] = offsetF; 
    DiffQFilterOffsetSign[sub_pos] = sign ? -1:1; 
  }
}
#endif  // E_DAIF

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

#ifndef EDAIF2 
void fillAIFSubpel(int sub_pos, FILTER_SIZE, FILTER_OFFSET, int *QFCoeff)
#else
void fillAIFSubpel(int sub_pos, int *QFCoeff)
#endif 
{
  int i, j; 
  if(sub_pos == b_pos)
  {
    for(i = FILTER_SIZE/2; i < FILTER_SIZE; i++) 
      QFCoeff[FILTER_SIZE*FILTER_OFFSET+i] = QFCoeff[FILTER_SIZE*FILTER_OFFSET+(FILTER_SIZE-1-i)];
  }
  else if(sub_pos == h_pos)
  {
    for(i = FILTER_SIZE/2; i < FILTER_SIZE; i++) 
      QFCoeff[FILTER_SIZE*i+FILTER_OFFSET] = QFCoeff[FILTER_SIZE*(FILTER_SIZE-1-i)+FILTER_OFFSET];
  }
  else if(sub_pos == f_pos)
  {
    for(j = 0; j < FILTER_SIZE; j++) 
    for(i = FILTER_SIZE/2; i < FILTER_SIZE; i++) 
      QFCoeff[FILTER_SIZE*j+i] = QFCoeff[FILTER_SIZE*j+(FILTER_SIZE-1-i)];
  }
  else if(sub_pos == i_pos)
  {
    for(j = FILTER_SIZE/2; j < FILTER_SIZE; j++) 
    for(i = 0; i < FILTER_SIZE; i++) 
      QFCoeff[FILTER_SIZE*j+i] = QFCoeff[FILTER_SIZE*(FILTER_SIZE-1-j)+i];
  }
  else if(sub_pos == j_pos)
  {
    for(j = 0; j < FILTER_SIZE/2; j++) 
    for(i = FILTER_SIZE/2; i < FILTER_SIZE; i++) 
      QFCoeff[FILTER_SIZE*j+i] = QFCoeff[FILTER_SIZE*j+(FILTER_SIZE-1-i)];
    for(j = FILTER_SIZE/2; j < FILTER_SIZE; j++) 
    for(i = 0; i < FILTER_SIZE; i++) 
      QFCoeff[FILTER_SIZE*j+i] = QFCoeff[FILTER_SIZE*(FILTER_SIZE-1-j)+i];
  }
  else
    return; 
}

void readAIFSubpel(int sub_pos, Bitstream *bitstream)
{
  int i; 

  extern int numQBits1DH[SQR_FILTER-1];
  extern int numQBits1DV[SQR_FILTER-1];
  extern int numQBits2D [SQR_FILTER-1];
  extern int POS_EQUATION_NUMBER[15];
  extern int DiffQFilterCoef[15][SQR_FILTER];
#ifndef EDAIF2
  extern int FILTER_SIZE;
  extern int FILTER_OFFSET;
#endif

  if((sub_pos+1)/4==0)  // 1D horizontal
  {
    for(i = 0; i < POS_EQUATION_NUMBER[sub_pos]; i++)
      DiffQFilterCoef[sub_pos][FILTER_SIZE*FILTER_OFFSET+i] = 
      readCoeff1("sub_pos 1D horiz", numQBits1DH[FILTER_SIZE*FILTER_OFFSET+i], bitstream);
  }
  else if((sub_pos+1)%4 == 0) // 1D vertical
  {
    for(i = 0; i < POS_EQUATION_NUMBER[sub_pos]; i++)
      DiffQFilterCoef[sub_pos][FILTER_SIZE*i+FILTER_OFFSET] = 
      readCoeff1("sub_pos 1D vert", numQBits1DV[FILTER_SIZE*i+FILTER_OFFSET], bitstream);
  }
  else // 2D 
  {
    for(i = 0; i < SQR_FILTER-1; i++)
      if(UniqueCoeff[sub_pos][i])  
        DiffQFilterCoef[sub_pos][i] = 
        readCoeff1("sub_pos 2D", numQBits2D[i], bitstream);
  }
#ifndef EDAIF2
  fillAIFSubpel(sub_pos, FILTER_SIZE, FILTER_OFFSET, DiffQFilterCoef[sub_pos]);
#else
  fillAIFSubpel(sub_pos, DiffQFilterCoef[sub_pos]);
#endif

}

#endif


#ifdef EDAIF2
void readEDAIF2(Bitstream *currStream)
{
  extern int POS_EQUATION_NUMBER[15];
  extern int DiffQFilterCoef[15][SQR_FILTER];
  extern int FilterFlag[15];            

  extern int UseAllSubpelPositions;
  extern int SubpelPositionsPattern;
  extern int  RealTimeDecomposition;
  extern int nPairsSets[4];
  extern int realSymCommands[4][MAX_NUM_AIF];
  extern int SymmetryPosition[MAX_NUM_AIF];
  extern int numQBits1DH[], numQBits1DV[], numQBits2D[], numQBitsDIF[];
  extern int Calc2Filt_Indexes[MAX_NUM_AIF][SQR_FILTER];
  int nQBits; 
  int DecompFlag;

  int i,j,sub_pos;
  int bitFlag;

  
  for(sub_pos = 0; sub_pos <= o_pos; sub_pos++)
  {
    if((sub_pos == f_pos)||(sub_pos == i_pos)||(sub_pos == j_pos)||(sub_pos == k_pos)||(sub_pos == n_pos))
    {
      bitFlag = u_v(1,"FilterFlag-1", currStream);;
      if(bitFlag == 0)
        FilterFlag[sub_pos] = 0;
      else
      { // use AIF
        bitFlag = u_v(1,"FilterFlag-2", currStream);
        if (bitFlag == 1) 
          FilterFlag[sub_pos] = 3; // Filter 3: codeword '11'
        else 
          FilterFlag[sub_pos] = 1; // Filter 1: codeword '10'
      }
    }
    else
      FilterFlag[sub_pos] = u_v(1,"FilterFlag", currStream);;
  }

  SetCalc2FilterIdx();
  RealTimeDecomposition = u_1("RealTimeDecomposition", currStream);

  if (RealTimeDecomposition == 1)
  {
    DecompFlag = u_1("FullDecomposition", currStream);
    if(DecompFlag)
    {
      for(i = 0;i<4;i++)
        for (j = 0; j < nPairsSets[i]; j ++)
          realSymCommands[i][j] = 1;
    }else
    {
      for(i = 0;i<4;i++)
      {
        for (j = 0; j < nPairsSets[i]; j ++)
        {
          realSymCommands[i][j] = u_1("realSymCommands[i][j]", currStream);
        }
      }
    }
    ArrangeSymmetryEDAIF2();
  }// end FilterDecomposition



  for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
  {
     
    if((FilterFlag[sub_pos]%2!=0)&&(SymmetryPosition[sub_pos]))
    {
      for(i= 0; i < POS_EQUATION_NUMBER[sub_pos]; i++)
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

        DiffQFilterCoef[sub_pos][i] = readCoeff1("SH: DiffQFilterCoef", nQBits, currStream);
      }
    }
  }

//  PrintFilterCoefInt(DiffQFilterCoef);
}

#endif
#ifdef DIRECTIONAL_FILTER
void readFilterCoefs(int filterID,Bitstream *currStream)
{
  extern int POS_EQUATION_NUMBER[15];
  extern int DiffQFilterCoef[15][SQR_FILTER];
  extern int FilterFlag[15];            

#ifndef EDAIF2
  extern int FILTER_SIZE;
  extern int FILTER_OFFSET;
#endif
  extern int UseAllSubpelPositions;
  extern int SubpelPositionsPattern;

  int i,j,sub_pos;
  if (filterID==FILTER_TYPE_1D)
  {
    UseAllSubpelPositions = u_1("use_all_subpel_positions", currStream);
    if(UseAllSubpelPositions)
    {
      for(sub_pos=a_pos; sub_pos <= o_pos; sub_pos++)
        FilterFlag[sub_pos] = 1;
      SubpelPositionsPattern =   (FilterFlag[j_pos]<<4) + (FilterFlag[f_pos]<<3) + (FilterFlag[e_pos]<<2) + 
        (FilterFlag[b_pos]<<1) +  FilterFlag[a_pos];
    }
    else
    {
      SubpelPositionsPattern = ue_v("positions_pattern", currStream);
      FilterFlag[a_pos] = 
        FilterFlag[c_pos] = 
        FilterFlag[d_pos] =
        FilterFlag[l_pos] =  SubpelPositionsPattern       & 1;
      FilterFlag[b_pos] = 
        FilterFlag[h_pos] = (SubpelPositionsPattern >> 1) & 1;
      FilterFlag[e_pos] = 
        FilterFlag[g_pos] = 
        FilterFlag[m_pos] =
        FilterFlag[o_pos] = (SubpelPositionsPattern >> 2) & 1;
      FilterFlag[f_pos] = 
        FilterFlag[i_pos] = 
        FilterFlag[k_pos] =
        FilterFlag[n_pos] = (SubpelPositionsPattern >> 3) & 1;
      FilterFlag[j_pos] = (SubpelPositionsPattern >> 4) & 1;
    }
    if(FilterFlag[a_pos])
      for(i= 0; i < POS_EQUATION_NUMBER[a_pos]; i++)
        DiffQFilterCoef[a_pos][FILTER_SIZE*FILTER_OFFSET+i] = se_v("a_pos", currStream);
    if(FilterFlag[b_pos])
      for(i= 0; i < POS_EQUATION_NUMBER[b_pos]; i++)
        DiffQFilterCoef[b_pos][FILTER_SIZE*FILTER_OFFSET+i] = se_v("b_pos", currStream);
    if(FilterFlag[e_pos])
      for(i = 0; i < FILTER_SIZE; i++)
      {
        j = i;//NW-SE
        DiffQFilterCoef[e_pos][FILTER_SIZE*i+j] = se_v("e_pos", currStream);
      }
      if(FilterFlag[f_pos])
      {
        for(i = 0; i < FILTER_SIZE/*/2*/; i++)
        {
          j = i;//NW-SE
          DiffQFilterCoef[f_pos][FILTER_SIZE*i+j] = se_v("f_pos", currStream);
          /*            }
          for(i = FILTER_SIZE/2; i< FILTER_SIZE; i++)
          {
          j = FILTER_SIZE-i-1;//NW-SE*
          DiffQFilterCoef[f_pos][FILTER_SIZE*i+j] = se_v("f_pos", currStream);*/
        }
      }

      if(FilterFlag[j_pos])
        for(i = 0; i < FILTER_SIZE/2; i++)
        {
          j = i;//NW-SE
          DiffQFilterCoef[j_pos][FILTER_SIZE*i+j] = se_v("j_pos", currStream);
        }


  }
  //PrintFilterCoefInt(DiffQFilterCoef);
  for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
    ExtendFilterCoefficientsInt(sub_pos, DiffQFilterCoef);
  //PrintFilterCoefInt(DiffQFilterCoef);

}
#endif

#ifdef E_DAIF
void readFilterCoefs_EDAIF(int filterID,Bitstream *currStream)
{
  int sub_pos;
  extern int SymmetryPosition[15];
  extern int FilterFlag[15];   
  extern  int DiffQFilterOffsetI[15];
  extern int DiffQFilterOffsetF[15];
  extern int DiffQFilterOffsetSign[15];

  extern int FilterFlagInt;

  if (!(filterID == FILTER_TYPE_EDAIF))
  {
    printf("ERROR: Wrong filter type @readFilterCoefs_EDAIF!\n");
    return;
  }

  readEDAIF2(currStream);

  // full-pel filter flag
  FilterFlagInt = u_1("FilterFlagInt", currStream); 
  if(FilterFlagInt)
    readAIFInteger(currStream);

  // sub-pel filter offsets
  for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
  {
    if((FilterFlag[sub_pos]%2!=0) && (SymmetryPosition[sub_pos]))
      readAIFOffset(sub_pos, currStream); 
    else
    {
      DiffQFilterOffsetI[sub_pos] = 0; 
      DiffQFilterOffsetF[sub_pos] = 0; 
      DiffQFilterOffsetSign[sub_pos] = 0; 

    }
  }
  // full-pel filter flag
  if(FilterFlagInt)
    readAIFOffset(-1, currStream); 
}
#endif  // E_DAIF

#ifdef EAIF

void readFilterCoefs_EAIF(int filterID,Bitstream *currStream)
{
  extern int POS_EQUATION_NUMBER[15];
  extern int DiffQFilterCoef[15][SQR_FILTER];
  extern int FilterFlag[15];			  		
#ifndef EDAIF2
  extern int FILTER_SIZE;
  extern int FILTER_OFFSET;
#endif
  extern int UseAllSubpelPositions;
  extern int SubpelPositionsPattern;
  extern int FilterFlagInt;
  extern int SymmetryPosition[15];

  int sub_pos;

  // sub-pel filter flags
  for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
  {
    if(SymmetryPosition[sub_pos])
    	FilterFlag[sub_pos] = u_1("sub-pel filter flag", currStream);
  }

  // filter coeffs
  for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
  {
    if(FilterFlag[sub_pos] && SymmetryPosition[sub_pos])
      readAIFSubpel(sub_pos, currStream); 
  }

  // full-pel filter flag
  FilterFlagInt = u_1("FilterFlagInt", currStream); 
  if(FilterFlagInt)
    readAIFInteger(currStream);

  // sub-pel filter offsets
  for(sub_pos = a_pos; sub_pos <= o_pos; sub_pos++)
  {
    if(FilterFlag[sub_pos] && SymmetryPosition[sub_pos])
      readAIFOffset(sub_pos, currStream); 
  }
  // full-pel filter flag
  if(FilterFlagInt)
    readAIFOffset(-1, currStream); 
}
#endif

