/*!
************************************************************************
* \file mv_competition.c
*
* \brief
*    competition based motion vector coding
*
* \author
*  Joel Jung              <joelb.jung@orange-ftgroup.com>           \n
*  Guillaume Laroche      <guillaume.laroche@orange-ftgroup.com>    \n
*
* This software may be used for research purposes only.
* Orange FTR&D would also appreciate that any technical report, conference 
* or journal paper which uses the software includes one of the following references:
*    - J. Jung, G. Laroche, "Competition-Based Scheme for Motion Vector Selection and Coding: new results",
*      VCEG contribution VCEG-???, Geneva, November 2006.
*    - J. Jung, G. Laroche, "Competition-Based Scheme for Motion Vector Selection and Coding",
*      VCEG contribution VCEG-AC06, Klagenfurt, July 2006.
*    - G. Laroche, J. Jung, B. Pesquet-Popescu, "A Spatio Temporal Competing Scheme for the 
*      Rate Distortion Optimized Selection and Coding of Motion Vectors",
*      EUSIPCO'06, Firenza, September 2006.
************************************************************************
*/


// CONFIGURATION
////////////////

// Configuration of the .Cfg file:
// In section KTA STUFF:
//
// MV_Competition        = 0 disabled
//                       = 1 enabled with default parameters for motion vectors of p and b frames and skip
//                       = 2 enabled with user parameters
//
// Predictors_skip, Predictors_MVp, Predictors_MVb are user parameters, describes the predictors that are enabled.
//
// Order of predictors for Predictors_skip: H.264 Median - ExtendedSpatial - a - b - c - 0 - Collocated - Empty
// Order of predictors for Predictors_MVp and MVb: H.264 Median - Empty - a - b - c - 0 - Collocated - Empty 
// Examples:
//          Predictors_skip     = '10100000' : 2 predictors enabled : 'H.264 Median' and 'a'
//          Predictors_skip     = '10000010' : 2 predictors enabled : 'H.264 Median' and 'Colocated'
//          Predictors_skip     = '10000000' : 1 predictor enabled : 'H.264 Median'
//
//          Predictors_MVp      = '10000000' : 1 predictor enabled : 'H.264 Median'
//          Predictors_MVp      = '10000010' : 2 predictors enabled : 'H.264 Median' and 'Colocated'
//

// KNOWN LIMITATIONS
////////////////////
// Enabling 'Empty' predictors or '00000000' configuration will lead to errors.
// Not compatible with:                
//    - interlace coding
//    - RdOpt != 1     

// VERSIONS
//////////
//////////
// 1.0  31 October 06 - first implementation for P and B frames 
// 1.1  9 February 07 - compatibility added with UseFME=3 and other tools
// 1.2  20 June 08    - compatibility added with RDO_Q
//                      compatibility added with AIF=1, 2, 3
//                      compatibility added with 1/8pel
//                      compatibility added with APEC
//                      compatibility added with hierarchical B slices
//                      MVCompetition chroma bug fix (slightly modifies the results)
//                      cleaning - improvements to increase compliance with JM coding rules


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "defines.h"
#include "global.h"
#include "image.h"
#include "mb_access.h"
#include "biaridecod.h"
#include "vlc.h"



#ifdef MV_COMPETITION
#include "mv_competition.h"

MV_Competition mv_comp;
extern int assignSE2partition[][SE_MAX_ELEMENTS];

int *****mv_predictors;        // Stores the predictor for each sub-partition and for each predictor  
// predMVXY[block_x][block_y][mode][Ref_Frame][list][prediction number][xy];


/*!
************************************************************************
* \brief
*    Read MVP for Skip Mode
*
************************************************************************
*/
int readPredictorForSkip(struct img_par  *img, struct inp_par *inp)
{
  SyntaxElement currSE;
  Slice *currSlice    = img->currentSlice;
  DataPartition *dP;
  int *partMap        = assignSE2partition[currSlice->dp_mode];
  
  currSE.type = SE_MV_PREDICTOR;
  dP = &(currSlice->partArr[partMap[currSE.type]]);
  currSE.value2 = (mv_comp.nb_mode_for_skip /2) + (mv_comp.nb_mode_for_skip % 2);
  
  if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
  {
    currSE.value1 = u_v(currSE.value2, "predictor for the SKIP mode", dP->bitstream);  
  }
  else
  {  
    currSE.reading = read_MV_predictor_CABAC;
    dP->readSyntaxElement(&currSE,img,inp,dP);
  }
  
  return currSE.value1;
}

/*!
************************************************************************
* \brief
*    Read MVP for Skip Mode
*
************************************************************************
*/
int readPredictorMV(struct img_par  *img, struct inp_par *inp)
{
  SyntaxElement currSE;
  Slice *currSlice    = img->currentSlice;
  DataPartition *dP;
  int *partMap        = assignSE2partition[currSlice->dp_mode];
  int maxmode=0;
  
  currSE.type = SE_MV_PREDICTOR;
  dP = &(currSlice->partArr[partMap[currSE.type]]);
  if(img->type == P_SLICE)
  {
    maxmode = mv_comp.nb_mode_for_mvp;
    currSE.value2 = (mv_comp.nb_mode_for_mvp / 2) + (mv_comp.nb_mode_for_mvp % 2);
  }
  else if (img->type == B_SLICE)
  {
    maxmode = mv_comp.nb_mode_for_mvb;
    currSE.value2 = (mv_comp.nb_mode_for_mvb / 2) + (mv_comp.nb_mode_for_mvb % 2);
  }
  
  if(maxmode > 1)
  {
    if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
    {
      currSE.value1 = u_v(currSE.value2, "predictor MV", dP->bitstream);  
    }
    else
    {  
      currSE.reading = read_MV_predictor_CABAC;
      dP->readSyntaxElement(&currSE,img,inp,dP);
    }
  }else currSE.value1 = 0;
  
  //  if(img->type == B_SLICE) printf("MB%i\t%i\n",img->current_mb_nr,currSE.value1);
  return currSE.value1;
}
/*!
************************************************************************
* \brief
*    Read Motion vector predictor in H.264 bitstream
*
************************************************************************
*/
void read_MV_predictor_CABAC( SyntaxElement *se,
                             struct inp_par *inp,
                             struct img_par *img,
                             DecodingEnvironmentPtr dep_dp)
{
  Macroblock          *currMB = &img->mb_data[img->current_mb_nr];
  int act_ctx = 0;
  MotionInfoContexts *ctx = (img->currentSlice)->mot_ctx;
  int nb_symbol;
  int two_power = 2;
  
  if(se->value2 != 0)
  {
    if(currMB->mb_type == 0)
      se->value1 = biari_decode_symbol(dep_dp, &ctx->mv_predictor_skip_contexts[act_ctx]);
    else 
      if(img->type == P_SLICE)
        se->value1 = biari_decode_symbol(dep_dp, &ctx->mv_predictor_mvp_contexts[act_ctx]);
      else //B_SLICE
        se->value1 = biari_decode_symbol(dep_dp, &ctx->mv_predictor_mvb_contexts[act_ctx]);
      
      nb_symbol = se->value2-1;
      
      while(nb_symbol != 0)
      {
        if(currMB->mb_type == 0)
          se->value1 = (se->value1) + (two_power * biari_decode_symbol(dep_dp, &ctx->mv_predictor_skip_contexts[act_ctx]));
        else 
          if(img->type == P_SLICE)
            se->value1 = (se->value1) + (two_power * biari_decode_symbol(dep_dp, &ctx->mv_predictor_mvp_contexts[act_ctx]));
          else //B_SLICE
            se->value1 = (se->value1) + (two_power * biari_decode_symbol(dep_dp, &ctx->mv_predictor_mvb_contexts[act_ctx]));
          
          nb_symbol--;
          two_power = two_power *2;
      }  
  }
  
}


void reinit_MV_Competition()
{
  int i;
  
  mv_comp.nb_mode_for_skip = 0;
  mv_comp.nb_mode_for_mvp = 0;
  mv_comp.nb_mode_for_mvb = 0;
  
  for (i=0; i<MAX_MV_PREDICTOR; i++)
  {
    if (mv_comp.predictors_skip[i] == 1)
    {
      mv_comp.predictor_for_skip[mv_comp.nb_mode_for_skip++]=i;
    }
    if (mv_comp.predictors_mvp[i] == 1)
    {
      mv_comp.predictor_for_mvp[mv_comp.nb_mode_for_mvp++]=i;
    }
    if (mv_comp.predictors_mvb[i] == 1)
    {
      mv_comp.predictor_for_mvb[mv_comp.nb_mode_for_mvb++]=i;
    }
  }
  if(img->type ==  B_SLICE)
    if( type_of_the_previous_frame != B_SLICE)  
      successive_Bframe = (listX[LIST_1][0]->poc/2 - listX[LIST_0][0]->poc/2)/(dec_picture->frame_poc/2-listX[LIST_0][0]->poc/2)-1;
    
    type_of_the_previous_frame = img->type;
}

void init_MV_Competition()
{
  int i,j, k, l;
  
  mv_comp.predictor_for_skip = (int *) malloc(sizeof(int) * MAX_MV_PREDICTOR);
  mv_comp.predictor_for_mvp  = (int *) malloc(sizeof(int) * MAX_MV_PREDICTOR);
  mv_comp.predictor_for_mvb  = (int *) malloc(sizeof(int) * MAX_MV_PREDICTOR);
  
#ifdef MB32X32_MVC
  mv_predictors = (int*****)malloc((4<<MAX_MB_EXT_LEVEL)*sizeof(int****)) ;
  for(l=0;l<(4<<MAX_MB_EXT_LEVEL);l++)
  {
    (mv_predictors)[l] = (int****)malloc((4<<MAX_MB_EXT_LEVEL)*sizeof(int***)) ;
    for(k=0;k<(4<<MAX_MB_EXT_LEVEL);k++)
#else
  mv_predictors = (int*****)malloc(4*sizeof(int****)) ;
  for(l=0;l<4;l++)
  {
    (mv_predictors)[l] = (int****)malloc(4*sizeof(int***)) ;
    for(k=0;k<4;k++)
#endif
    {
      (mv_predictors)[l][k] = (int***)malloc(MAXMODE*sizeof(int**)) ;
      for(j=0;j<(MAXMODE);j++)
      {
        mv_predictors[l][k][j] = (int**)malloc(15*sizeof(int*)) ;                     //int predMVXY[block_x][block_y][MAXMODE][15][2];  
        for (i=0; i<15; i++)
        {
          (mv_predictors)[l][k][j][i] = (int*)malloc(2*sizeof(int));
        } 
      }
    }
  }
}

void init_MV_Competition_mv_previous_tab()
{
  int i,j, k;
  
  if(successive_Bframe > 1)
  { 
    if(mv_previous_B_frame == NULL)
    {
      mv_previous_B_frame = (short****)malloc(2*sizeof(short***));
      for (k=0; k<2; k++)
      {
        mv_previous_B_frame[k] = (short***)malloc((img->height/4)*sizeof(short**));
        for(i=0;i<img->height/4;i+=1)
        {
          mv_previous_B_frame[k][i] = (short**)malloc((img->width/4)*sizeof(short*));
          for(j=0;j<img->width/4;j+=1)
          {
            mv_previous_B_frame[k][i][j] = (short*)malloc(2*sizeof(short));
          }
        }
      }
    }
    if(ref_idx_previous_B_frame == NULL)
    {
      ref_idx_previous_B_frame = (char***)malloc(2*sizeof(char**));
      for (k=0; k<2; k++)
      {
        ref_idx_previous_B_frame[k] = (char**)malloc((img->height/4)*sizeof(char*));
        for(i=0;i<img->height/4;i+=1)
        {
          ref_idx_previous_B_frame[k][i] = (char*)malloc((img->width/4)*sizeof(char));
        }
      }
    }
  }
}
void close_MV_Competition()
{
  int i, j, k, l;
  free (mv_comp.predictor_for_skip);
  free (mv_comp.predictor_for_mvp);
  free (mv_comp.predictor_for_mvb);
  
#ifdef MB32X32_MVC
  for(l=0;l<(4<<MAX_MB_EXT_LEVEL);l++)
  {
    for(k=0;k<(4<<MAX_MB_EXT_LEVEL);k++)
#else
  for(l=0;l<4;l++)
  {
    for(k=0;k<4;k++)
#endif
    {
      for(j=0;j<(MAXMODE);j++)
      {
        for (i=0; i<15; i++)
        {
          free(mv_predictors[l][k][j][i]);
        }
        free(mv_predictors[l][k][j]);
      }
      free(mv_predictors[l][k]);
    }
    free(mv_predictors[l]);
  }
  free(mv_predictors);
  
  if(successive_Bframe > 1)
  { 
    for (k=0; k<2; k++)
    {
      for(i=0;i<img->height/4;i+=1)
      {
        for(j=0;j<img->width/4;j+=1)
        {
          free(mv_previous_B_frame[k][i][j]);
        }
        free(mv_previous_B_frame[k][i]);
      }
      free(mv_previous_B_frame[k]);
    }
    free(mv_previous_B_frame);
    
    
    for (k=0; k<2; k++)
    {
      for(i=0;i<img->height/4;i+=1)
      {
        free(ref_idx_previous_B_frame[k][i]);
      }
      free(ref_idx_previous_B_frame[k]);
    }
    free(ref_idx_previous_B_frame);
    
  }
}

short collocated_mv_available(int pos_y, int pos_x, int list)
{
  short return_val = TRUE;
  
  if (listX[list][0]->ref_idx[0][pos_y*4][pos_x*4] == -1)  // Collocated block in INTRA
  {
    
    return_val = FALSE;
  }
  
  return (return_val);
}

short collocated_mv_available_previous_B_frame(int pos_y, int pos_x, int list)
{
  short return_val = TRUE;
  
  // Macrobloc belongs to the picture ?
  //if ((pos_v < 0) || (pos_h < 0) || (pos_v>=(img->height/16)) || (pos_h>=(img->width/16)))
  //  return_val = 0;
  
  // Macrobloc has a motion vector ?
  //else 
  //  if (coding_modes_previous_frame[img->current_mb_nr * 4] >= 9)
  if (ref_idx_previous_B_frame[list][pos_y*4][pos_x*4] == -1)  // Collocated block in INTRA
  {
    return_val = FALSE;
  }
  
  return (return_val);
}

short read_index_for_skip_mode()
{
  short mv_x_ref, mv_y_ref;
  short i;
  short return_val = FALSE;      // Assumes it can be guessed
  short nb_of_available = 0;    // Number of modes available
  
  mv_x_ref = mv_comp.mv_pred_skip[0][0];
  mv_y_ref = mv_comp.mv_pred_skip[0][1];
  
  // If all the mv_skip_jj (for all modes except 'best_mode_skip' are equal to NOT_AVAILABLE
  // then the mode can be guessed
  
  
  for (i=0; i<mv_comp.nb_mode_for_skip; i++)
  {
    //Check if among all modes enabled, more than one is AVAILABLE 
    if (mv_comp.mv_pred_skip[i][0] != NOT_AVAILABLE) 
      nb_of_available++;
  }
  
  if (nb_of_available == 1)
    return_val = FALSE;
  else
  {
    for (i=1; i<mv_comp.nb_mode_for_skip; i++)
    {
      if ((mv_comp.mv_pred_skip[i][0] != mv_x_ref) 
        || (mv_comp.mv_pred_skip[i][1] != mv_y_ref))
        
        return_val = TRUE;    // One of the component for one prediction is different
      // We can return 0, the mode cannot be guessed
    }
  }
  
  return (return_val);
}

void SetMotionVectorPredictor_Competition (struct img_par  *img,
                                           struct inp_par *inp,
                                           short           *pmv_x,
                                           short           *pmv_y,
                                           char            ref_frame,
                                           byte            list,
                                           char            ***refPic,
                                           short           ****tmp_mv,
                                           int             block_x,
                                           int             block_y,
                                           int             blockshape_x,
                                           int             blockshape_y)
{
  int mb_x                 = 4*block_x;
  int mb_y                 = 4*block_y;
  int mb_nr                = img->current_mb_nr;
  
  //  int  pred_vec=0;
  
  short  pmv2[2] = {0,0};
  short predictor;
  
  int maxmode=0;
  int *predictors=NULL;  
  
  int modeblock;
  int prediction_mode;
  
  
  if (img->type == I_SLICE) 
    maxmode = 1;
  else if (img->type == P_SLICE)
  {
    maxmode = mv_comp.nb_mode_for_mvp;
    predictors = mv_comp.predictor_for_mvp;
  }
  else if (img->type == B_SLICE)
  {
    maxmode = mv_comp.nb_mode_for_mvb;
    predictors = mv_comp.predictor_for_mvb;
  }
  
  modeblock= determine_mode_block(blockshape_x,blockshape_y);
  
  predictor =0;
  
  for (prediction_mode = 0; prediction_mode<maxmode; prediction_mode++)
  {
    if (predictors[prediction_mode] == PRED_H264_MEDIAN)
    {
      SetMotionVectorPredictorMedian(img,pmv_x, pmv_y, ref_frame, list, refPic, tmp_mv, block_x, block_y, blockshape_x, blockshape_y);
      pmv2[0] = *pmv_x;
      pmv2[1] = *pmv_y;
      
    }
    
    else if (predictors[prediction_mode] == PRED_ZERO)
    {
      pmv2[0] = 0;
      pmv2[1] = 0;
    }
    
    else if (predictors[prediction_mode] == PRED_A)
    {
      PixelPos block_a;
      
      getLuma4x4Neighbour(img->current_mb_nr, block_x, block_y,           -1,  0, &block_a);
      pmv2[0] = block_a.available  ? tmp_mv[list][block_a.pos_y][block_a.pos_x][0] : 0;
      pmv2[1] = block_a.available  ? tmp_mv[list][block_a.pos_y][block_a.pos_x][1] : 0;
    }
    
    else if (predictors[prediction_mode] == PRED_B)
    {
      PixelPos  block_b;
      
      getLuma4x4Neighbour(mb_nr, block_x, block_y,            0, -1, &block_b);
      pmv2[0] = block_b.available  ? tmp_mv[list][block_b.pos_y][block_b.pos_x][0] : 0;
      pmv2[1] = block_b.available  ? tmp_mv[list][block_b.pos_y][block_b.pos_x][1] : 0;  
    }
    
    else if (predictors[prediction_mode] == PRED_C)
    {
      PixelPos block_c, block_d;
      
      getLuma4x4Neighbour(mb_nr, block_x, block_y, blockshape_x, -1, &block_c);
      getLuma4x4Neighbour(mb_nr, block_x, block_y,           -1, -1, &block_d);
      
      if (mb_y > 0)
      {
        if (mb_x < 8) 
        {
          if (mb_y==8)
          {
            if (blockshape_x == 16)      block_c.available  = 0;
          }
          else
          {
            if (mb_x+blockshape_x == 8)  block_c.available  = 0;
          }
        }
        else
        {
          if (mb_x+blockshape_x == 16)   block_c.available  = 0;
        }
      }
      
      if (!block_c.available)
        block_c=block_d;
      
      pmv2[0] = block_c.available  ? tmp_mv[list][block_c.pos_y][block_c.pos_x][0] : 0;
      pmv2[1] = block_c.available  ? tmp_mv[list][block_c.pos_y][block_c.pos_x][1] : 0;
    }
    else if (predictors[prediction_mode] == PRED_COLOCATED)
    {  
      int y=(int)img->current_mb_nr/(img->width/16);  
      int x=(int)img->current_mb_nr%(img->width/16);
      
      if (img->type == P_SLICE)
      {
        if (collocated_mv_available(y, x, LIST_0) == TRUE)
        {
          
          pmv2[0] = listX[0][0]->mv[LIST_0][y*4+block_y][x*4+block_x][0] * (ref_frame+1) / (listX[0][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
          pmv2[1] = listX[0][0]->mv[LIST_0][y*4+block_y][x*4+block_x][1] * (ref_frame+1) / (listX[0][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
        }
        else
        {
          pmv2[0] = 0; 
          pmv2[1] = 0; 
        }
      } 
      else if (img->type == B_SLICE)
      {  
        // <FTRD : Compatibility with hierarchical B frames
        //SetMotionVectorPredictor_collocated_B_SLICE(pmv_x, pmv_y, y, x, list, ref_frame, block_y, block_x);
        
        SetMotionVectorPredictor_collocated_HB_SLICE(pmv_x, pmv_y, y, x, list, ref_frame, block_y, block_x);
        
        // FTRD>
        
        pmv2[0] = *pmv_x;
        pmv2[1] = *pmv_y;
      }
    }
    
    
    mv_predictors[block_x][block_y][modeblock][prediction_mode][0]=pmv2[0];
    mv_predictors[block_x][block_y][modeblock][prediction_mode][1]=pmv2[1];
  }
  
  
{
  short tab_for_skip_mode = 0;
  short current_predictor_index;
  
  int pred0_x, pred0_y, current_pred_x, current_pred_y; 
  
  pred0_x = mv_predictors[block_x][block_y][modeblock][0][0];
  pred0_y = mv_predictors[block_x][block_y][modeblock][0][1];
  
  for(current_predictor_index=1;current_predictor_index<maxmode ;current_predictor_index++)
  {
    current_pred_x = mv_predictors[block_x][block_y][modeblock][current_predictor_index][0];
    current_pred_y = mv_predictors[block_x][block_y][modeblock][current_predictor_index][1];
    
    if ((current_pred_x == pred0_x) && (current_pred_y == pred0_y))
      tab_for_skip_mode++;  
  }
  
  if(tab_for_skip_mode == maxmode-1) 
    predictor = 0;
  else
    predictor = readPredictorMV(img,inp);
}
*pmv_x = mv_predictors[block_x][block_y][modeblock][predictor][0]; 
*pmv_y = mv_predictors[block_x][block_y][modeblock][predictor][1];
}

#ifdef MB32X32_MVC
extern int ResetAvailabilityOfNeighbors32(int mb_x, int mb_y, int *reset_block_x, int *reset_block_y, int blockshape_x, int blockshape_y);
extern void SetMotionVectorPredictor32 (struct img_par  *img,
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
                               int             mb_ext_level);

void SetMotionVectorPredictor_Competition32 (struct img_par  *img,
                                           struct inp_par *inp,
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
  int mb_x                 = 4*block_x;
  int mb_y                 = 4*block_y;
  int mb_nr                = img->current_mb_nr;
  
  //  int  pred_vec=0;
  
  short  pmv2[2] = {0,0};
  short predictor;
  
  int maxmode=0;
  int *predictors=NULL;  
  
  int modeblock;
  int prediction_mode;
  
//===Rahul
    int reset_mb_nr, reset_block_x=block_x, reset_block_y=block_y;
    reset_mb_nr = ResetAvailabilityOfNeighbors32(mb_x, mb_y, &reset_block_x, &reset_block_y, blockshape_x, blockshape_y);
    
//    block_x = reset_block_x;
//    block_y = reset_block_y;

//  img->current_mb_nr will give the top-left MB number of a current 32x32 block
//  reset_mb_nr will give the actual MB number 
//===

  if (img->type == I_SLICE) 
    maxmode = 1;
  else if (img->type == P_SLICE)
  {
    maxmode = mv_comp.nb_mode_for_mvp;
    predictors = mv_comp.predictor_for_mvp;
  }
  else if (img->type == B_SLICE)
  {
    maxmode = mv_comp.nb_mode_for_mvb;
    predictors = mv_comp.predictor_for_mvb;
  }
  
  modeblock = determine_mode_block(blockshape_x >> mb_ext_level,blockshape_y >> mb_ext_level);
  
  predictor =0;
  
  for (prediction_mode = 0; prediction_mode<maxmode; prediction_mode++)
  {
    if (predictors[prediction_mode] == PRED_H264_MEDIAN)
    {
//      SetMotionVectorPredictorMedian  (img,pmv_x, pmv_y, ref_frame, list, refPic, tmp_mv, block_x, block_y, blockshape_x, blockshape_y);

//      SetMotionVectorPredictorMedian()  is same as SetMotionVectorPredictor()..dont know why to have separate functions????
      SetMotionVectorPredictor32 (img,pmv_x, pmv_y, ref_frame, list, refPic, tmp_mv, reset_block_x, reset_block_y, blockshape_x, blockshape_y, mb_ext_level);
      pmv2[0] = *pmv_x;
      pmv2[1] = *pmv_y;
    }
    
    else if (predictors[prediction_mode] == PRED_ZERO)
    {
      pmv2[0] = 0;
      pmv2[1] = 0;
    }
    
    else if (predictors[prediction_mode] == PRED_A)
    {
      PixelPos block_a;
      
      getLuma4x4Neighbour(reset_mb_nr, reset_block_x, reset_block_y,           -1,  0, &block_a);
      pmv2[0] = block_a.available  ? tmp_mv[list][block_a.pos_y][block_a.pos_x][0] : 0;
      pmv2[1] = block_a.available  ? tmp_mv[list][block_a.pos_y][block_a.pos_x][1] : 0;
    }
    
    else if (predictors[prediction_mode] == PRED_B)
    {
      PixelPos  block_b;
      
      getLuma4x4Neighbour(reset_mb_nr, reset_block_x, reset_block_y,            0, -1, &block_b);
      pmv2[0] = block_b.available  ? tmp_mv[list][block_b.pos_y][block_b.pos_x][0] : 0;
      pmv2[1] = block_b.available  ? tmp_mv[list][block_b.pos_y][block_b.pos_x][1] : 0;  
    }
    
    else if (predictors[prediction_mode] == PRED_C)
    {
      PixelPos block_c, block_d;
      
      getLuma4x4Neighbour(reset_mb_nr, reset_block_x, reset_block_y, blockshape_x, -1, &block_c);
      getLuma4x4Neighbour(reset_mb_nr, reset_block_x, reset_block_y,           -1, -1, &block_d);
      
      if (mb_y > 0)
      {
        if (mb_x < (8<<mb_ext_level)) 
        {
          if (mb_y==(8<<mb_ext_level))
          {
            if (blockshape_x == (16<<mb_ext_level))      block_c.available  = 0;
          }
          else
          {
            if (mb_x+blockshape_x == (8<<mb_ext_level))  block_c.available  = 0;
          }
        }
        else
        {
          if (mb_x+blockshape_x == (16<<mb_ext_level))   block_c.available  = 0;
        }
      }
      
      if (!block_c.available)
        block_c=block_d;
      
      pmv2[0] = block_c.available  ? tmp_mv[list][block_c.pos_y][block_c.pos_x][0] : 0;
      pmv2[1] = block_c.available  ? tmp_mv[list][block_c.pos_y][block_c.pos_x][1] : 0;
    }
    else if (predictors[prediction_mode] == PRED_COLOCATED)
    {  
      int y=(int)reset_mb_nr/(img->width/16);  // Vertical
      int x=(int)reset_mb_nr%(img->width/16);  // Horizontal
      //Rahul---(x,y) is the 16x16MB location in the frame
      int ytop=(int)mb_nr/(img->width/16);  // Vertical
      int xtop=(int)mb_nr%(img->width/16);  // Horizontal

      if (img->type == P_SLICE)
      {
        if (collocated_mv_available(y, x, LIST_0) == TRUE)
        {
//        pmv2[0] = listX[0][0]->mv[LIST_0][y*4+block_y][x*4+block_x][0] * (ref_frame+1) / (listX[0][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
//        pmv2[1] = listX[0][0]->mv[LIST_0][y*4+block_y][x*4+block_x][1] * (ref_frame+1) / (listX[0][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
//        pmv2[0] = listX[0][0]->mv[LIST_0][y*4][x*4][0] * (ref_frame+1) / (listX[0][0]->ref_idx[0][y*4][x*4] + 1);
//        pmv2[1] = listX[0][0]->mv[LIST_0][y*4][x*4][1] * (ref_frame+1) / (listX[0][0]->ref_idx[0][y*4][x*4] + 1);
          pmv2[0] = listX[0][0]->mv[LIST_0][ytop*4+block_y][xtop*4+block_x][0] * (ref_frame+1) / (listX[0][0]->ref_idx[0][ytop*4+block_y][xtop*4+block_x] + 1);
          pmv2[1] = listX[0][0]->mv[LIST_0][ytop*4+block_y][xtop*4+block_x][1] * (ref_frame+1) / (listX[0][0]->ref_idx[0][ytop*4+block_y][xtop*4+block_x] + 1);

        }
        else
        {
          pmv2[0] = 0; 
          pmv2[1] = 0; 
        }
      } 
      else if (img->type == B_SLICE)
      {  
        // <FTRD : Compatibility with hierarchical B frames
        //SetMotionVectorPredictor_collocated_B_SLICE(pmv_x, pmv_y, y, x, list, ref_frame, block_y, block_x);

//      SetMotionVectorPredictor_collocated_HB_SLICE(pmv_x, pmv_y, y   , x   , list, ref_frame, block_y, block_x);
        SetMotionVectorPredictor_collocated_HB_SLICE(pmv_x, pmv_y, ytop, xtop, list, ref_frame, block_y, block_x);
        
        // FTRD>
        
        pmv2[0] = *pmv_x;
        pmv2[1] = *pmv_y;
      }
    }
    
    
    mv_predictors[block_x][block_y][modeblock][prediction_mode][0]=pmv2[0];
    mv_predictors[block_x][block_y][modeblock][prediction_mode][1]=pmv2[1];
  }
  
  
{
  short tab_for_skip_mode = 0;
  short current_predictor_index;
  
  int pred0_x, pred0_y, current_pred_x, current_pred_y; 
  
  pred0_x = mv_predictors[block_x][block_y][modeblock][0][0];
  pred0_y = mv_predictors[block_x][block_y][modeblock][0][1];
  
  for(current_predictor_index=1;current_predictor_index<maxmode ;current_predictor_index++)
  {
    current_pred_x = mv_predictors[block_x][block_y][modeblock][current_predictor_index][0];
    current_pred_y = mv_predictors[block_x][block_y][modeblock][current_predictor_index][1];
    
    if ((current_pred_x == pred0_x) && (current_pred_y == pred0_y))
      tab_for_skip_mode++;  
  }
  
  if(tab_for_skip_mode == maxmode-1) 
    predictor = 0;
  else
    predictor = readPredictorMV(img,inp);
}
*pmv_x = mv_predictors[block_x][block_y][modeblock][predictor][0]; 
*pmv_y = mv_predictors[block_x][block_y][modeblock][predictor][1];
}


#endif
void SetMotionVectorPredictorMedian (struct img_par  *img,
                                     short           *pmv_x,
                                     short           *pmv_y,
                                     char            ref_frame,
                                     byte            list,
                                     char            ***refPic,
                                     short           ****tmp_mv,
                                     int             block_x,
                                     int             block_y,
                                     int             blockshape_x,
                                     int             blockshape_y)
{
  int mb_x                 = BLOCK_SIZE*block_x;
  int mb_y                 = BLOCK_SIZE*block_y;
  int mb_nr                = img->current_mb_nr;
  
  int mv_a, mv_b, mv_c, pred_vec=0;
  int mvPredType, rFrameL=0, rFrameU=0, rFrameUR=0;
  int hv;
  
  
  PixelPos block_a, block_b, block_c, block_d;
  
  
  getLuma4x4Neighbour(mb_nr, block_x, block_y,           -1,  0, &block_a);
  getLuma4x4Neighbour(mb_nr, block_x, block_y,            0, -1, &block_b);
  getLuma4x4Neighbour(mb_nr, block_x, block_y, blockshape_x, -1, &block_c);
  getLuma4x4Neighbour(mb_nr, block_x, block_y,           -1, -1, &block_d);
  
  if (mb_y > 0)
  {
    if (mb_x < 8)  // first column of 8x8 blocks
    {
      if (mb_y==8)
      {
        if (blockshape_x == 16)      block_c.available  = 0;
        else                         block_c.available &= 1;
      }
      else
      {
        if (mb_x+blockshape_x != 8)  block_c.available &= 1;
        else                         block_c.available  = 0;
      }
    }
    else
    {
      if (mb_x+blockshape_x != 16)   block_c.available &= 1;
      else                           block_c.available  = 0;
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
  /*else
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
                                     }*/
  
  
  // Prediction if only one of the neighbors uses the reference frame
  /// we are checking
  //
  if(rFrameL == ref_frame && rFrameU != ref_frame && rFrameUR != ref_frame)       mvPredType = MVPRED_L;
  else if(rFrameL != ref_frame && rFrameU == ref_frame && rFrameUR != ref_frame)  mvPredType = MVPRED_U;
  else if(rFrameL != ref_frame && rFrameU != ref_frame && rFrameUR == ref_frame)  mvPredType = MVPRED_UR;
  // Directional predictions 
  if(blockshape_x == 8 && blockshape_y == 16)
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
  else if(blockshape_x == 16 && blockshape_y == 8)
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
    /*  if (img->mb_data[img->current_mb_nr].mb_field)
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
      else*/
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

int determine_mode_block(int blockshape_x,int blockshape_y)
{
  if((blockshape_x==16) && ( blockshape_y==16))
    return 1;
  if((blockshape_x==16) && ( blockshape_y==8))
    return 2;
  if((blockshape_x==8) && ( blockshape_y==16))
    return 3;
  if((blockshape_x==8) && ( blockshape_y==8))
    return 4;
  if((blockshape_x==8) && ( blockshape_y==4))
    return 5;
  if((blockshape_x==4) && ( blockshape_y==8))
    return 6;
  if((blockshape_x==4) && ( blockshape_y==4))
    return 7;
  return 0;
}


void SetMotionVectorPredictor_collocated_B_SLICE (short           *pmv_x, 
                                                  short           *pmv_y,
                                                  int y,
                                                  int x,
                                                  int list,
                                                  int ref_frame,
                                                  int block_y, 
                                                  int block_x)
{
  int colocated_MV_x, colocated_MV_y;
  int b_frame_to_code = Bframe_ctr%successive_Bframe;
  
  *pmv_x = 0;
  *pmv_y = 0;
  
  if(list==LIST_0)
  {
    if((b_frame_to_code)==0)//if this B_frame is the first B_frame of the GOP
    {
      if (collocated_mv_available(y, x, LIST_1) == TRUE)
      { 
        colocated_MV_x = listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][0] * (ref_frame+1) / (listX[LIST_1][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
        colocated_MV_y = listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][1] * (ref_frame+1) / (listX[LIST_1][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
        
        *pmv_x = colocated_MV_x/(successive_Bframe+1)+(ref_frame)*colocated_MV_x;
        *pmv_y = colocated_MV_y/(successive_Bframe+1)+(ref_frame)*colocated_MV_y;
      }
      else
      {
        *pmv_x = 0;
        *pmv_y = 0;
      }
    }
    else
    {
      if (collocated_mv_available_previous_B_frame(y, x, LIST_1) == TRUE)
      {
        colocated_MV_x = mv_previous_B_frame[LIST_1][y*4+block_y][x*4+block_x][0];
        colocated_MV_y = mv_previous_B_frame[LIST_1][y*4+block_y][x*4+block_x][1];
        
        *pmv_x = -((b_frame_to_code + 1))*colocated_MV_x / (successive_Bframe+1-((b_frame_to_code + 1)-1));
        *pmv_y = -((b_frame_to_code + 1))*colocated_MV_y / (successive_Bframe+1-((b_frame_to_code + 1)-1));
        
        if ((collocated_mv_available(y, x, LIST_1) == TRUE) && (ref_frame != 0))
        {
          *pmv_x = *pmv_x  + (ref_frame*listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][0] / (listX[LIST_1][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1));
          *pmv_y = *pmv_y  + (ref_frame*listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][1] / (listX[LIST_1][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1));
        }
      }
      else 
      {
        if (collocated_mv_available_previous_B_frame(y, x, LIST_0) == TRUE)
        {
          colocated_MV_x = ((b_frame_to_code + 1)-1) * mv_previous_B_frame[LIST_0][y*4+block_y][x*4+block_x][0] / (ref_idx_previous_B_frame[LIST_0][y*4+block_y][x*4+block_x] * (successive_Bframe+1) + ((b_frame_to_code + 1)-1));
          colocated_MV_y = ((b_frame_to_code + 1)-1) * mv_previous_B_frame[LIST_0][y*4+block_y][x*4+block_x][1] / (ref_idx_previous_B_frame[LIST_0][y*4+block_y][x*4+block_x] * (successive_Bframe+1) + ((b_frame_to_code + 1)-1));
          
          *pmv_x = ((b_frame_to_code + 1))*colocated_MV_x / ((b_frame_to_code + 1)-1);  
          *pmv_y = ((b_frame_to_code + 1))*colocated_MV_y / ((b_frame_to_code + 1)-1);
          
          if ((collocated_mv_available(y, x, LIST_1) == TRUE) && (ref_frame != 0))
          {
            *pmv_x = *pmv_x + (ref_frame*listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][0] / (listX[LIST_1][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1));
            *pmv_y = *pmv_y + (ref_frame*listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][1] / (listX[LIST_1][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1));
          }
        }
        else
        {
          *pmv_x = 0;
          *pmv_y = 0;
        }
      }
    }
    
  }
  else
  {//LIST_1
    if((b_frame_to_code)==0)//if this B_frame is the first B_frame of the GOP
    {
      if (collocated_mv_available(y, x, LIST_1) == TRUE)
      { 
        colocated_MV_x = listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][0] * (ref_frame+1) / (listX[LIST_1][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
        colocated_MV_y = listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][1] * (ref_frame+1) / (listX[LIST_1][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
        
        *pmv_x = -successive_Bframe*colocated_MV_x/(successive_Bframe+1);
        *pmv_y = -successive_Bframe*colocated_MV_y/(successive_Bframe+1);
      }
      else
      {
        *pmv_x = 0;
        *pmv_y = 0;
      }
    }
    else
    {
      if (collocated_mv_available_previous_B_frame(y, x, LIST_1) == TRUE)
      {
        colocated_MV_x = mv_previous_B_frame[LIST_1][y*4+block_y][x*4+block_x][0];
        colocated_MV_y = mv_previous_B_frame[LIST_1][y*4+block_y][x*4+block_x][1];
        
        *pmv_x = (successive_Bframe+1-(b_frame_to_code + 1))*colocated_MV_x / (successive_Bframe+1-((b_frame_to_code + 1)-1));
        *pmv_y = (successive_Bframe+1-(b_frame_to_code + 1))*colocated_MV_y / (successive_Bframe+1-((b_frame_to_code + 1)-1));
        
      }
      else
      {
        if (collocated_mv_available_previous_B_frame(y, x, LIST_0) == TRUE)
        {
          colocated_MV_x = ((b_frame_to_code + 1)-1) * mv_previous_B_frame[LIST_0][y*4+block_y][x*4+block_x][0] / (ref_idx_previous_B_frame[LIST_0][y*4+block_y][x*4+block_x] * (successive_Bframe+1) + ((b_frame_to_code + 1)-1));
          colocated_MV_y = ((b_frame_to_code + 1)-1) * mv_previous_B_frame[LIST_0][y*4+block_y][x*4+block_x][1] / (ref_idx_previous_B_frame[LIST_0][y*4+block_y][x*4+block_x] * (successive_Bframe+1) + ((b_frame_to_code + 1)-1));
          
          *pmv_x = -(successive_Bframe+1-(b_frame_to_code + 1))*colocated_MV_x / ((b_frame_to_code + 1)-1);
          *pmv_y = -(successive_Bframe+1-(b_frame_to_code + 1))*colocated_MV_y / ((b_frame_to_code + 1)-1);            
        }
        else
        {
          *pmv_x = 0;
          *pmv_y = 0;
        }
      }
    }
  }
  
}

void Copy_MV_B_frame(short ****  mv, char  ***   ref_idx)
{ 
  int i, j;
  
  for(i=0;i<img->height/4;i+=1)
  {
    for(j=0;j<img->width/4;j+=1)
    {
      mv_previous_B_frame[LIST_0][i][j][0] = mv[LIST_0][i][j][0];
      mv_previous_B_frame[LIST_0][i][j][1] = mv[LIST_0][i][j][1];
      ref_idx_previous_B_frame[LIST_0][i][j] = ref_idx[LIST_0][i][j];
      
      mv_previous_B_frame[LIST_1][i][j][0] = mv[LIST_1][i][j][0];
      mv_previous_B_frame[LIST_1][i][j][1] = mv[LIST_1][i][j][1];
      ref_idx_previous_B_frame[LIST_1][i][j] = ref_idx[LIST_1][i][j];
    }
  }
  
}

/*!
************************************************************************
* \brief
*    Set motion vector predictor
************************************************************************
*/

void SetMotionVectorPredictor_Skip (short  *pmv_x, short *pmv_y, char   ***refPic, short  ****tmp_mv, short  ref_frame,
                                    int    list, int    block_x, int    block_y, int    blockshape_x, int    blockshape_y,
                                    short mode_skip)
{ 
  int zeroMotionAbove;
  int zeroMotionLeft;
  //PixelPos mb_a, mb_b;
  int      a_mv_y = 0;
  int      a_ref_idx = 0;
  int      b_mv_y = 0;
  int      b_ref_idx = 0;
  //  int i=(int)img->current_mb_nr/(img->width/16);
  //  int j=(int)img->current_mb_nr%(img->width/16);
  
  
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  
  PixelPos block_a, block_b, block_c, block_d;
  int mb_nr = img->current_mb_nr;
  
  getLuma4x4Neighbour(mb_nr, block_x, block_y,           -1,  0, &block_a);
  getLuma4x4Neighbour(mb_nr, block_x, block_y,            0, -1, &block_b);
  getLuma4x4Neighbour(mb_nr, block_x, block_y, blockshape_x, -1, &block_c);
  getLuma4x4Neighbour(mb_nr, block_x, block_y,           -1, -1, &block_d);
  
  //  Same as in the original JM... might not be optimal
  if (mv_comp.predictor_for_skip[mode_skip] == PRED_H264_MEDIAN)
  {
    if (block_a.available)
    {
      a_mv_y    = dec_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][1];
      a_ref_idx = dec_picture->ref_idx[LIST_0][block_a.pos_y][block_a.pos_x];
      
      if (currMB->mb_field && !img->mb_data[block_a.mb_addr].mb_field)
      {
        a_mv_y    /=2;
        a_ref_idx *=2;
      }
      if (!currMB->mb_field && img->mb_data[block_a.mb_addr].mb_field)
      {
        a_mv_y    *=2;
        a_ref_idx >>=1;
      }
    }
    
    if (block_b.available)
    {
      b_mv_y    = dec_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][1];
      b_ref_idx = dec_picture->ref_idx[LIST_0][block_b.pos_y][block_b.pos_x];
      
      if (currMB->mb_field && !img->mb_data[block_b.mb_addr].mb_field)
      {
        b_mv_y    /=2;
        b_ref_idx *=2;
      }
      if (!currMB->mb_field && img->mb_data[block_b.mb_addr].mb_field)
      {
        b_mv_y    *=2;
        b_ref_idx >>=1;
      }
    }
    
    zeroMotionLeft  = !block_a.available ? 1 : a_ref_idx==0 && dec_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][0]==0 && a_mv_y==0 ? 1 : 0;
    zeroMotionAbove = !block_b.available ? 1 : b_ref_idx==0 && dec_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][0]==0 && b_mv_y==0 ? 1 : 0;
    
    if (zeroMotionAbove || zeroMotionLeft)
    {
      
      *pmv_x = 0;
      *pmv_y = 0;
    }
    else
    {
      SetMotionVectorPredictor (img, pmv_x, pmv_y, 0, LIST_0, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);
    }    
  }
  
  else if (mv_comp.predictor_for_skip[mode_skip] == PRED_A)
  {  
    if (block_a.available) /*&& (img->mb_data[img->current_mb_nr-1].mb_type<9))*/ 
    {
      *pmv_x = dec_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][0];
      *pmv_y = dec_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][1];
    }
    else 
    {
      *pmv_x = NOT_AVAILABLE;
      *pmv_y = NOT_AVAILABLE;
    }
  }
  else if (mv_comp.predictor_for_skip[mode_skip] == PRED_B)
  {  
    if (block_b.available) /*&& (img->mb_data[img->current_mb_nr-1].mb_type<9))*/ 
    {
      *pmv_x = dec_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][0];
      *pmv_y = dec_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][1];
    }
    else 
    {
      *pmv_x = NOT_AVAILABLE;
      *pmv_y = NOT_AVAILABLE;
    }
  }
  else if (mv_comp.predictor_for_skip[mode_skip] == PRED_C)
  {  
    if (block_c.available) /*&& (img->mb_data[img->current_mb_nr-1].mb_type<9))*/ 
    {
      *pmv_x = dec_picture->mv[LIST_0][block_c.pos_y][block_c.pos_x][0];
      *pmv_y = dec_picture->mv[LIST_0][block_c.pos_y][block_c.pos_x][1];
    }
    else 
    {
      *pmv_x = NOT_AVAILABLE;
      *pmv_y = NOT_AVAILABLE;
    }
  }
  else if (mv_comp.predictor_for_skip[mode_skip] == PRED_ZERO)
  {  
    *pmv_x = 0;
    *pmv_y = 0;
  }
  else if (mv_comp.predictor_for_skip[mode_skip] == PRED_COLOCATED)
  {
    int y=(int)img->current_mb_nr/(img->width/16);  // Vertical
    int x=(int)img->current_mb_nr%(img->width/16);  // Horizontal
    
    if (collocated_mv_available(y, x, LIST_0) == TRUE)
    {
      
      *pmv_x = listX[0][0]->mv[LIST_0][y*4][x*4][0] / (listX[0][0]->ref_idx[0][y*4][x*4] + 1);
      *pmv_y = listX[0][0]->mv[LIST_0][y*4][x*4][1] / (listX[0][0]->ref_idx[0][y*4][x*4] + 1);
    }
    else
    {
      *pmv_x = NOT_AVAILABLE;
      *pmv_y = NOT_AVAILABLE;
    }
  }
  else if (mv_comp.predictor_for_skip[mode_skip] == PRED_EXTENDEDSPATIAL)  
  {
    int mv_a, mv_b, mv_c;
    
    if ((block_a.available) && (block_b.available) && (block_c.available))  
    {
      mv_a = dec_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][0];
      mv_b = dec_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][0];
      mv_c = dec_picture->mv[LIST_0][block_c.pos_y][block_c.pos_x][0];      
      *pmv_x = mv_a+mv_b+mv_c-min(mv_a,min(mv_b,mv_c))-max(mv_a,max(mv_b,mv_c));
      
      mv_a = dec_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][1];
      mv_b = dec_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][1];
      mv_c = dec_picture->mv[LIST_0][block_c.pos_y][block_c.pos_x][1];      
      *pmv_y = mv_a+mv_b+mv_c-min(mv_a,min(mv_b,mv_c))-max(mv_a,max(mv_b,mv_c));
    }
    else if (block_a.available) 
    {
      *pmv_x = dec_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][0];        
      *pmv_y = dec_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][1];        
      
    }
    else if (block_b.available)  
    {
      *pmv_x = dec_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][0];
      *pmv_y = dec_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][1];
      
    }
    else if (block_c.available)  
    {
      *pmv_x = dec_picture->mv[LIST_0][block_c.pos_y][block_c.pos_x][0];
      *pmv_y = dec_picture->mv[LIST_0][block_c.pos_y][block_c.pos_x][1];
      
    }
    else
    {
      *pmv_x = 0;
      *pmv_y = 0;
    }
  }
  
  
}

#ifdef MB32X32_MVC
void SetMotionVectorPredictor_Skip32 (short  *pmv_x, short *pmv_y, char   ***refPic, short  ****tmp_mv, short  ref_frame,
                                    int    list, int    block_x, int    block_y, int    blockshape_x, int    blockshape_y,
                                    short mode_skip, int mb_ext_level)
{ 
  int zeroMotionAbove;
  int zeroMotionLeft;
  //PixelPos mb_a, mb_b;
  int      a_mv_y = 0;
  int      a_ref_idx = 0;
  int      b_mv_y = 0;
  int      b_ref_idx = 0;
  //  int i=(int)img->current_mb_nr/(img->width/16);
  //  int j=(int)img->current_mb_nr%(img->width/16);
  
  int mb_x                 = BLOCK_SIZE*block_x;
  int mb_y                 = BLOCK_SIZE*block_y;
  
  
  PixelPos block_a, block_b, block_c, block_d;
  
  {
    int reset_mb_nr, reset_block_x=block_x, reset_block_y=block_y;

    reset_mb_nr = ResetAvailabilityOfNeighbors32(mb_x, mb_y, &reset_block_x, &reset_block_y, blockshape_x, blockshape_y);

    getLuma4x4Neighbour(reset_mb_nr, reset_block_x, reset_block_y,           -1,  0, &block_a);
    getLuma4x4Neighbour(reset_mb_nr, reset_block_x, reset_block_y,            0, -1, &block_b);
    getLuma4x4Neighbour(reset_mb_nr, reset_block_x, reset_block_y, blockshape_x, -1, &block_c);
    getLuma4x4Neighbour(reset_mb_nr, reset_block_x, reset_block_y,           -1, -1, &block_d);
  }

/*
  getLuma4x4Neighbour(mb_nr, block_x, block_y,           -1,  0, &block_a);
  getLuma4x4Neighbour(mb_nr, block_x, block_y,            0, -1, &block_b);
  getLuma4x4Neighbour(mb_nr, block_x, block_y, blockshape_x, -1, &block_c);
  getLuma4x4Neighbour(mb_nr, block_x, block_y,           -1, -1, &block_d);
*/  
  //  Same as in the original JM... might not be optimal
  if (mv_comp.predictor_for_skip[mode_skip] == PRED_H264_MEDIAN)
  {
    if (block_a.available)
    {
      a_mv_y    = dec_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][1];
      a_ref_idx = dec_picture->ref_idx[LIST_0][block_a.pos_y][block_a.pos_x];
    }
    
    if (block_b.available)
    {
      b_mv_y    = dec_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][1];
      b_ref_idx = dec_picture->ref_idx[LIST_0][block_b.pos_y][block_b.pos_x];
    }
    
    zeroMotionLeft  = !block_a.available ? 1 : a_ref_idx==0 && dec_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][0]==0 && a_mv_y==0 ? 1 : 0;
    zeroMotionAbove = !block_b.available ? 1 : b_ref_idx==0 && dec_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][0]==0 && b_mv_y==0 ? 1 : 0;
    
    if (zeroMotionAbove || zeroMotionLeft)
    {
      *pmv_x = 0;
      *pmv_y = 0;
    }
    else
    {
      SetMotionVectorPredictor32 (img, pmv_x, pmv_y, 0, LIST_0, dec_picture->ref_idx, dec_picture->mv, 0, 0, blockshape_x, blockshape_y, mb_ext_level);
    }    
  }
  
  else if (mv_comp.predictor_for_skip[mode_skip] == PRED_A)
  {  
    if (block_a.available) /*&& (img->mb_data[img->current_mb_nr-1].mb_type<9))*/ 
    {
      *pmv_x = dec_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][0];
      *pmv_y = dec_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][1];
    }
    else 
    {
      *pmv_x = NOT_AVAILABLE;
      *pmv_y = NOT_AVAILABLE;
    }
  }
  else if (mv_comp.predictor_for_skip[mode_skip] == PRED_B)
  {  
    if (block_b.available) /*&& (img->mb_data[img->current_mb_nr-1].mb_type<9))*/ 
    {
      *pmv_x = dec_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][0];
      *pmv_y = dec_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][1];
    }
    else 
    {
      *pmv_x = NOT_AVAILABLE;
      *pmv_y = NOT_AVAILABLE;
    }
  }
  else if (mv_comp.predictor_for_skip[mode_skip] == PRED_C)
  {  
    if (block_c.available) /*&& (img->mb_data[img->current_mb_nr-1].mb_type<9))*/ 
    {
      *pmv_x = dec_picture->mv[LIST_0][block_c.pos_y][block_c.pos_x][0];
      *pmv_y = dec_picture->mv[LIST_0][block_c.pos_y][block_c.pos_x][1];
    }
    else 
    {
      *pmv_x = NOT_AVAILABLE;
      *pmv_y = NOT_AVAILABLE;
    }
  }
  else if (mv_comp.predictor_for_skip[mode_skip] == PRED_ZERO)
  {  
    *pmv_x = 0;
    *pmv_y = 0;
  }
  else if (mv_comp.predictor_for_skip[mode_skip] == PRED_COLOCATED)
  {
    int y=(int)img->current_mb_nr/(img->width/16);  // Vertical
    int x=(int)img->current_mb_nr%(img->width/16);  // Horizontal
    
    if (collocated_mv_available(y, x, LIST_0) == TRUE)
    {
      
      *pmv_x = listX[0][0]->mv[LIST_0][y*4][x*4][0] / (listX[0][0]->ref_idx[0][y*4][x*4] + 1);
      *pmv_y = listX[0][0]->mv[LIST_0][y*4][x*4][1] / (listX[0][0]->ref_idx[0][y*4][x*4] + 1);
    }
    else
    {
      *pmv_x = NOT_AVAILABLE;
      *pmv_y = NOT_AVAILABLE;
    }
  }
  else if (mv_comp.predictor_for_skip[mode_skip] == PRED_EXTENDEDSPATIAL)  
  {
    int mv_a, mv_b, mv_c;
    
    if ((block_a.available) && (block_b.available) && (block_c.available))  
    {
      mv_a = dec_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][0];
      mv_b = dec_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][0];
      mv_c = dec_picture->mv[LIST_0][block_c.pos_y][block_c.pos_x][0];      
      *pmv_x = mv_a+mv_b+mv_c-min(mv_a,min(mv_b,mv_c))-max(mv_a,max(mv_b,mv_c));
      
      mv_a = dec_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][1];
      mv_b = dec_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][1];
      mv_c = dec_picture->mv[LIST_0][block_c.pos_y][block_c.pos_x][1];      
      *pmv_y = mv_a+mv_b+mv_c-min(mv_a,min(mv_b,mv_c))-max(mv_a,max(mv_b,mv_c));
    }
    else if (block_a.available) 
    {
      *pmv_x = dec_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][0];        
      *pmv_y = dec_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][1];        
    }
    else if (block_b.available)  
    {
      *pmv_x = dec_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][0];
      *pmv_y = dec_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][1];
    }
    else if (block_c.available)  
    {
      *pmv_x = dec_picture->mv[LIST_0][block_c.pos_y][block_c.pos_x][0];
      *pmv_y = dec_picture->mv[LIST_0][block_c.pos_y][block_c.pos_x][1];
    }
    else
    {
      *pmv_x = 0;
      *pmv_y = 0;
    }
  }
}
#endif



// <FTRD : Compatibility with hierarchical B frames
void init_mvscale_hb()
{
  int i, j, k, l;
  int iTRb = 0 ; int iTRp = 0;
  int prescale = 0;
  int num_ref_max = 0;
  
  //printf("Encoding picture : %d\n",dec_picture->poc/2);
  //printf("[point_list_cur_B] [ref_num_B] [point_list_coloc] [ref_num_coloc] [mv_scale]\n");
  
  // Mv in current B frame points to LIST_0 or LIST_1
  for(i = 0 ; i < 2 ; i++)
  {
    for(k = 0 ; k <listXsize[i] ; k++)
    {
      iTRb = Clip3(-128, 127, dec_picture->poc - listX[i][k]->poc);
      
      // The colocated vector points to LIST_0 or LIST_1
      for(j=0;j<2;j++)
      {
        if(listX[LIST_1][0]->slice_type == B_SLICE) 
        {
#ifdef MB32X32_MVC
//          num_ref_max = 2; 
          if(j==0) num_ref_max = img->num_ref_idx_l0_active + 1;  
          if(j==1) num_ref_max = img->num_ref_idx_l1_active + 1;  
#else
          if(j==0) num_ref_max = img->num_ref_idx_l0_active;  
          if(j==1) num_ref_max = img->num_ref_idx_l0_active;  
#endif
        }
        if(listX[LIST_1][0]->slice_type == P_SLICE) 
        {
          num_ref_max = dpb.num_ref_frames;
        }
        
        for(l = 0 ; l < num_ref_max ; l++)
        {
          iTRp = (int)Clip3(-128, 127, listX[LIST_1][0]->poc - (listX[LIST_1][0]->ref_pic_num[0][j][l]/2));
          
          if(iTRp==0) img->mvscale_hb[i][j][k][l]=9999;
          else 
          {
            prescale=(16384 + absm(iTRp/2))/iTRp;
            //printf("B = %d , ref_B = %d , colc = %d , ref_col = %d\n",enc_picture->poc/2,
            img->mvscale_hb[i][j][k][l] = Clip3(-2048, 2047,(iTRb*prescale + 32) >> 6);
            //printf("LIST_%d   %d   LIST_%d   %d   %f\n",i,k,j,l,(float)iTRb/(float)iTRp);
          }
        }
      }
    }
  }
}
// FTRD>

// <FTRD : Compatibility with hierarchical B frames
void SetMotionVectorPredictor_collocated_HB_SLICE (short *pmv_x,short *pmv_y,int y,int x,int list,int ref_frame,int block_y,int block_x)
{
  int colocated_MV_x = 0, colocated_MV_y = 0;
  int ref_frame_coloc = 0;
  *pmv_x = 0;
  *pmv_y = 0;
  
  // printf("%d, %d\n",listX[LIST_1][0]->ref_idx[LIST_0][y*4+block_y][x*4+block_x],listX[LIST_1][0]->ref_idx[LIST_1][y*4+block_y][x*4+block_x]);
  
  // B frame mv points to LIST_0
  if(list == LIST_0)
  {
    // If colocated mv points to LIST_0
    if(listX[LIST_1][0]->ref_idx[LIST_0][y*4+block_y][x*4+block_x] != -1)
    {
      colocated_MV_x = listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][0];
      colocated_MV_y = listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][1];
      
      ref_frame_coloc = listX[LIST_1][0]->ref_idx[LIST_0][y*4+block_y][x*4+block_x];
      
      // Scaling
      *pmv_x = (img->mvscale_hb[LIST_0][LIST_0][ref_frame][ref_frame_coloc] * colocated_MV_x + 128) >> 8;
      *pmv_y = (img->mvscale_hb[LIST_0][LIST_0][ref_frame][ref_frame_coloc] * colocated_MV_y + 128) >> 8;
    }
    
    else
    {
      // If colocated mv points to LIST_1
      if(listX[LIST_1][0]->ref_idx[LIST_1][y*4+block_y][x*4+block_x] != -1)
      {
        colocated_MV_x = listX[LIST_1][0]->mv[LIST_1][y*4+block_y][x*4+block_x][0];
        colocated_MV_y = listX[LIST_1][0]->mv[LIST_1][y*4+block_y][x*4+block_x][1];
        
        ref_frame_coloc = listX[LIST_1][0]->ref_idx[LIST_1][y*4+block_y][x*4+block_x];
        
        // Scaling
        *pmv_x = (img->mvscale_hb[LIST_0][LIST_1][ref_frame][ref_frame_coloc] * colocated_MV_x + 128) >> 8;
        *pmv_y = (img->mvscale_hb[LIST_0][LIST_1][ref_frame][ref_frame_coloc] * colocated_MV_y + 128) >> 8;
      }
      
      // No colocated mv is available.
      else
      {
        *pmv_x = 0;
        *pmv_y = 0;
      }
    }
  }
  
  // B frame mv points to LIST_1
  else
  {
    // If colocated mv points to LIST_0
    if(listX[LIST_1][0]->ref_idx[LIST_0][y*4+block_y][x*4+block_x] != -1)
    {
      colocated_MV_x = listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][0];
      colocated_MV_y = listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][1];
      
      ref_frame_coloc = listX[LIST_1][0]->ref_idx[LIST_0][y*4+block_y][x*4+block_x];
      
      // Scaling
      *pmv_x = (img->mvscale_hb[LIST_1][LIST_0][ref_frame][ref_frame_coloc] * colocated_MV_x + 128) >> 8;
      *pmv_y = (img->mvscale_hb[LIST_1][LIST_0][ref_frame][ref_frame_coloc] * colocated_MV_y + 128) >> 8;
    }
    
    else
    {
      // If colocated mv points to LIST_1
      if(listX[LIST_1][0]->ref_idx[LIST_1][y*4+block_y][x*4+block_x] != -1)
      {
        colocated_MV_x = listX[LIST_1][0]->mv[LIST_1][y*4+block_y][x*4+block_x][0];
        colocated_MV_y = listX[LIST_1][0]->mv[LIST_1][y*4+block_y][x*4+block_x][1];
        
        ref_frame_coloc = listX[LIST_1][0]->ref_idx[LIST_1][y*4+block_y][x*4+block_x];
        
        // Scaling
        *pmv_x = (img->mvscale_hb[LIST_1][LIST_1][ref_frame][ref_frame_coloc] * colocated_MV_x + 128) >> 8;
        *pmv_y = (img->mvscale_hb[LIST_1][LIST_1][ref_frame][ref_frame_coloc] * colocated_MV_y + 128) >> 8;
      }
      
      // No colocated mv is available.
      else
      {
        *pmv_x = 0;
        *pmv_y = 0;
      }
    }
  } 
}
// FTRD>
#endif
