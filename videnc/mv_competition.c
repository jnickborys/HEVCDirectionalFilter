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
*                              VCEG contribution VCEG-???, Geneva, November 2006.
*    - J. Jung, G. Laroche, "Competition-Based Scheme for Motion Vector Selection and Coding",
*                              VCEG contribution VCEG-AC06, Klagenfurt, July 2006.
*      - G. Laroche, J. Jung, B. Pesquet-Popescu, "A Spatio Temporal Competing Scheme for the 
*                Rate Distortion Optimized Selection and Coding of Motion Vectors",
*                EUSIPCO'06, Firenza, September 2006.
************************************************************************
*/


// CONFIGURATION
////////////////

// Configuration of the .Cfg file:
// In section KTA STUFF:
//
// MV_Competition     = 0 disabled
//                       = 1 enabled with default parameters for motion vectors of p and b frames and skip
//                       = 2 enabled with user parameters
//
// Predictors_skip, Predictors_MVp, Predictors_MVb are user parameters, describes the predictors that are enabled.
//
// Order of predictors for Predictors_skip: H.264 Median - ExtendedSpatial - a - b - c - 0 - Collocated - Empty
// Order of predictors for Predictors_MVp and MVb: H.264 Median - Empty - a - b - c - 0 - Collocated - Empty 
// Examples:
//      Predictors_skip     = '10100000' : 2 predictors enabled : 'H.264 Median' and 'a'
//      Predictors_skip     = '10000010' : 2 predictors enabled : 'H.264 Median' and 'Colocated'
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
#include "global.h"
#include "mbuffer.h"
#include "mb_access.h"
#include "cabac.h"
#include "string.h"
#include "vlc.h"
#include "fmo.h"
#include "fast_me.h"
#include "image.h"

#ifdef MV_COMPETITION
#include "mv_competition.h"

int skip_mode;
int *******mv_predictors;          // Stores the predictor for each sub-partition and for each predictor  
// predMVXY[block_x][block_y][mode][Ref_Frame][list][prediction number][xy];
int ******send_index_mv_prediction;        // int predModeMV2[mode][10][list][4][4][MAX_MODEPREDMV+1]; 
int *****predModeMV;          // contient les modes de prediction des predicteur MV [mode][ref][list][i][j]; int predModeMV[10][18*22];
int  *ranking;            //int  ranking[MAX_VAL];
int MAX_VAL ;    

extern int *mvbits;
extern int *assignSE2partition[2];
extern void SetMotionVectorPredictor (short  pmv[2],
                                      char   **refPic,
                                      short  ***tmp_mv,
                                      short  ref_frame,
                                      int    list,
                                      int    block_x,
                                      int    block_y,
                                      int    blockshape_x,
                                      int    blockshape_y);

MV_Competition mv_comp;



/*! 
*************************************************************************************
* \brief
*    initialization of motion vector competition arrays
*
* \param 
*    
*
* \return
*    
*************************************************************************************
*/

void init_MV_Competition()
{
  int l,n,m,i;
  int  j, k=0;
  int nLength;
  
  nLength = strlen(input->predictors_skip)-2;
  
  mv_comp.predictor_for_skip = (int *) malloc(sizeof(int) * MAX_MV_PREDICTOR);
  mv_comp.predictor_for_mvp  = (int *) malloc(sizeof(int) * MAX_MV_PREDICTOR);
  mv_comp.predictor_for_mvb  = (int *) malloc(sizeof(int) * MAX_MV_PREDICTOR);
  
  if (input->mv_competition == 0)
  {
    mv_comp.nb_mode_for_skip = 1;
    mv_comp.nb_mode_for_mvp  = 1;
    mv_comp.nb_mode_for_mvb  = 1;
    
    mv_comp.predictor_for_skip[0] = PRED_H264_MEDIAN;
    mv_comp.predictor_for_mvp[0]  = PRED_H264_MEDIAN;
    mv_comp.predictor_for_mvb[0]  = PRED_H264_MEDIAN;
  }
  else if (input->mv_competition == 1)
  {
    mv_comp.nb_mode_for_skip = 2;
    mv_comp.nb_mode_for_mvp  = 2;
    mv_comp.nb_mode_for_mvb  = 2;
    
    for (i=1; i<8; i++)
    {
      mv_comp.predictor_for_skip[i] = 0;
      mv_comp.predictor_for_mvp[i]  = 0;
      mv_comp.predictor_for_mvb[i]  = 0;
    }
    
    mv_comp.predictor_for_skip[0] = PRED_EXTENDEDSPATIAL;
    mv_comp.predictor_for_mvp[0]  = PRED_H264_MEDIAN;
    mv_comp.predictor_for_mvb[0]  = PRED_H264_MEDIAN;
    mv_comp.predictor_for_skip[1] = PRED_A;
    mv_comp.predictor_for_mvp[1]  = PRED_COLOCATED;
    mv_comp.predictor_for_mvb[1]  = PRED_COLOCATED;
  }
  else
  {
    mv_comp.nb_mode_for_skip = 0;
    mv_comp.nb_mode_for_mvp  = 0;
    mv_comp.nb_mode_for_mvb  = 0;
    
    // Check the number of predictors + enabled predictors for skip, motion vectors for P and B frames
    for(j=1;j<nLength+1;j++)
    {
      switch (input->predictors_skip[j])
      {
      case '1':
        mv_comp.predictor_for_skip[mv_comp.nb_mode_for_skip++] = j-1;
        break;
      case '0':
        break;
      default:
        snprintf(errortext, ET_SIZE, "Invalid Predictors_skip parameter.");
        error (errortext, 400);
        break;
      }
      
      switch (input->predictors_mvp[j])
      {
      case '1':
        mv_comp.predictor_for_mvp[mv_comp.nb_mode_for_mvp++] = j-1;
        break;
      case '0':
        break;
      default:
        snprintf(errortext, ET_SIZE, "Invalid Predictors_mvp parameter.");
        error (errortext, 400);
        break;
      }
      
      switch (input->predictors_mvb[j])
      {
      case '1':
        mv_comp.predictor_for_mvb[mv_comp.nb_mode_for_mvb++] = j-1;
        break;
      case '0':
        break;
      default:
        snprintf(errortext, ET_SIZE, "Invalid Predictors_mvb parameter.");
        error (errortext, 400);
        break;
      }
      
    }
  }
  
  mv_comp.mv_pred_skip = (int**) malloc(mv_comp.nb_mode_for_skip * sizeof(int*));  
  for(j=0;j<(mv_comp.nb_mode_for_skip);j++)
  {
    mv_comp.mv_pred_skip[j] = (int*) malloc(2 * sizeof(int)) ;             
  }
  
#ifdef MB32X32_MVC
  mv_predictors = (int*******)malloc((4<<MAX_MB_EXT_LEVEL)*sizeof(int******)) ;              //block_x
  for(l=0; l<(4<<MAX_MB_EXT_LEVEL); l++)
  {
    (mv_predictors)[l] = (int******)malloc((4<<MAX_MB_EXT_LEVEL)*sizeof(int*****)) ;            //block_y
    for(k=0; k<(4<<MAX_MB_EXT_LEVEL); k++)
#else
  mv_predictors = (int*******)malloc(4*sizeof(int******)) ;              //block_x
  for(l=0; l<4; l++)
  {
    (mv_predictors)[l] = (int******)malloc(4*sizeof(int*****)) ;            //block_y
    for(k=0; k<4; k++)
#endif
    {
      (mv_predictors)[l][k] = (int*****)malloc(MAXMODE*sizeof(int****)) ;      //mode
      for(j=0; j<(MAXMODE); j++)
      {
        (mv_predictors)[l][k][j]=(int****)malloc((input->num_ref_frames)*sizeof(int***)); //REF FRAME
        for(n=0; n<input->num_ref_frames; n++)
        {
          (mv_predictors)[l][k][j][n] = (int***)malloc(2*sizeof(int**)) ;      //LIST_0, LIST_1
          for(m=0; m<(2); m++)
          {
            (mv_predictors)[l][k][j][n][m] = (int**)malloc(15*sizeof(int*)) ;
            for (i=0; i<15; i++)
            {
              (mv_predictors)[l][k][j][n][m][i] = (int*)malloc(2*sizeof(int));  //x,y
            }
          }
        }
      }
    }
  }
  
  send_index_mv_prediction = (int******)malloc(MAXMODE*sizeof(int*****));
  
  for (k=0; k<MAXMODE; k++)
  {
    (send_index_mv_prediction)[k] = (int*****)malloc(img->num_ref_frames*sizeof(int****)) ;   // int predModeMV2[MAXMODE][10][list][18][22][MAX_MODEPREDMV+1]; 
    for (i=0; i<img->num_ref_frames; i++)
    {
      (send_index_mv_prediction)[k][i] = (int****)malloc((2)*sizeof(int***)); 
      for (m=0; m<2; m++)
      {
#ifdef MB32X32_MVC
        (send_index_mv_prediction)[k][i][m] = (int***)malloc((4<<MAX_MB_EXT_LEVEL)*sizeof(int**)); 
        for(l=0; l<(4<<MAX_MB_EXT_LEVEL); l++)
        {
          (send_index_mv_prediction)[k][i][m][l] = (int**)malloc((4<<MAX_MB_EXT_LEVEL)*sizeof(int*)); 
          for(j=0; j<(4<<MAX_MB_EXT_LEVEL); j++)
            (send_index_mv_prediction)[k][i][m][l][j] = (int*)malloc((max(mv_comp.nb_mode_for_mvp,mv_comp.nb_mode_for_mvb)+1)*sizeof(int));
#else
        (send_index_mv_prediction)[k][i][m] = (int***)malloc(4*sizeof(int**)); 
        for(l=0; l<4; l++)
        {
          (send_index_mv_prediction)[k][i][m][l] = (int**)malloc(4*sizeof(int*)); 
          for(j=0; j<4; j++)
            (send_index_mv_prediction)[k][i][m][l][j] = (int*)malloc((max(mv_comp.nb_mode_for_mvp,mv_comp.nb_mode_for_mvb)+1)*sizeof(int));
#endif
        }
      }
    }
  }
  
  predModeMV= (int*****)malloc(MAXMODE*sizeof(int****)) ;
  for(j=0; j<(MAXMODE); j++)
  {
    (predModeMV)[j]= (int****)malloc(img->num_ref_frames*sizeof(int***)) ;       //int predModeMV[MAXMODE][10][list][18][22];
    for (i=0; i<img->num_ref_frames; i++)
    {
      (predModeMV)[j][i] = (int***)malloc(2*sizeof(int**));
      for (m=0; m<2; m++)
      {
#ifdef MB32X32_MVC
        (predModeMV)[j][i][m] = (int**)malloc((4<<MAX_MB_EXT_LEVEL)*sizeof(int*));
        for (k=0; k<(4<<MAX_MB_EXT_LEVEL); k++)
          (predModeMV)[j][i][m][k] = (int*)malloc((4<<MAX_MB_EXT_LEVEL)*sizeof(int));
#else
        (predModeMV)[j][i][m] = (int**)malloc(4*sizeof(int*));
        for (k=0; k<4; k++)
          (predModeMV)[j][i][m][k] = (int*)malloc(4*sizeof(int));
#endif
      }
    }
  }
  
  if(input->successive_Bframe > 1)
  {
    mv_previous_B_frame = (short****)malloc(2*sizeof(short***));
    for (m=0; m<2; m++)
    {
      mv_previous_B_frame[m] = (short***)malloc((img->height/4)*sizeof(short**));
      for(i=0;i<img->height/4;i+=1)
      {
        mv_previous_B_frame[m][i] = (short**)malloc((img->width/4)*sizeof(short*));
        for(j=0;j<img->width/4;j+=1)
          mv_previous_B_frame[m][i][j] = (short*)malloc(2*sizeof(short));
      }
    }
    
    ref_idx_previous_B_frame = (char***)malloc(2*sizeof(char**));
    for (m=0; m<2; m++)
    {
      ref_idx_previous_B_frame[m] = (char**)malloc((img->height/4)*sizeof(char*));
      for(i=0;i<img->height/4;i+=1)
      {
        ref_idx_previous_B_frame[m][i] = (char*)malloc((img->width/4)*sizeof(char));
      }
    }
  }
  
  ranking =(int*)malloc (MAX_VAL*sizeof(int));                       //int  ranking[MAX_VAL];
  MAX_VAL=  max(mv_comp.nb_mode_for_mvp,mv_comp.nb_mode_for_mvb)*4;
  }
  
  
  /*! 
  *************************************************************************************
  * \brief
  *    free motion vector competition arrays
  *
  * \param 
  *    
  *
  * \return
  *    
  *************************************************************************************
  */
  
  void close_MV_Competition ()
  {
    int i, j, k, l, m, n;
    for(j=0;j<(mv_comp.nb_mode_for_skip);j++)
      free(mv_comp.mv_pred_skip[j]);
    
    free(mv_comp.predictor_for_skip);
    free(mv_comp.predictor_for_mvp );
    free(mv_comp.predictor_for_mvb);
    
#ifdef MB32X32_MVC
    for(l=0; l<(4<<MAX_MB_EXT_LEVEL); l++)
    {
      for(k=0; k<(4<<MAX_MB_EXT_LEVEL); k++)
#else
    for(l=0; l<4; l++)
    {
      for(k=0; k<4; k++)
#endif
      {
        for(j=0; j<(MAXMODE); j++)
        {
          for(n=0; n<input->num_ref_frames; n++)
          {
            for(m=0; m<(2); m++)
            {
              for (i=0; i<15; i++)
                free(mv_predictors[l][k][j][n][m][i]);
              
              free(mv_predictors[l][k][j][n][m]);
            }
            free(mv_predictors[l][k][j][n]);
          }
          free(mv_predictors[l][k][j]);
        }
        free(mv_predictors[l][k]);
      }
      free(mv_predictors[l]);
    }
    free(mv_predictors);
    
    for (k=0; k<MAXMODE; k++)
    {
      for (i=0; i<img->num_ref_frames; i++)
      {
        for (m=0; m<2; m++)
        {
#ifdef MB32X32_MVC
          for(l=0; l<(4<<MAX_MB_EXT_LEVEL); l++)
          {
            for(j=0; j<(4<<MAX_MB_EXT_LEVEL); j++)
#else
          for(l=0; l<4; l++)
          {
            for(j=0; j<4; j++)
#endif
              free(send_index_mv_prediction[k][i][m][l][j]);
            free(send_index_mv_prediction[k][i][m][l]);
          }
          free(send_index_mv_prediction[k][i][m]);
        }
        free(send_index_mv_prediction[k][i]);
      }
      free(send_index_mv_prediction[k]);
    }
    free(send_index_mv_prediction);
    
    for(j=0; j<(MAXMODE); j++)
    {
      for (i=0; i<img->num_ref_frames; i++)
      {
        for (m=0; m<2; m++)
        {
#ifdef MB32X32_MVC
          for (k=0; k<(4<<MAX_MB_EXT_LEVEL); k++)
#else
          for (k=0; k<4; k++)
#endif
            free(predModeMV[j][i][m][k]);
          free(predModeMV[j][i][m]);
        }
        free(predModeMV[j][i]);
      }
      free(predModeMV[j]);
    }
    free(predModeMV);
    
    if(input->successive_Bframe > 1)
    {
      for (m=0; m<2; m++)
      {
        for(i=0;i<img->height/4;i+=1)
        {
          for(j=0;j<img->width/4;j+=1)
          {
            free(mv_previous_B_frame[m][i][j]);
          }
          free(mv_previous_B_frame[m][i]);
        }
        free(mv_previous_B_frame[m]);
      }
      free(mv_previous_B_frame);
      
      for (m=0; m<2; m++)
      {
        for(i=0;i<img->height/4;i+=1)
        {
          free(ref_idx_previous_B_frame[m][i]);
        }
        free(ref_idx_previous_B_frame[m]);
      }
      free(ref_idx_previous_B_frame);
    }
    free(ranking);
  }
  
  /*! 
  *************************************************************************************
  * \brief
  *    writes the index of the selected motion vector predictor for the skip mode
  *
  * \param 
  *    
  *
  * \return int
  *    
  *************************************************************************************
  */
  
  int write_predictor_index_for_skip_mode()
  {
    int            rate       = 0;
    Macroblock*     currMB     = &img->mb_data[img->current_mb_nr];
    
    if (input->symbol_mode == CABAC)
    { currMB->send_index_predictor_skip = send_index_for_skip_prediction(currMB->best_predictor_for_skip);
    if (currMB->send_index_predictor_skip == TRUE)
      rate += writeMotionVectorPredictorSKIP(currMB->best_predictor_for_skip);
    }
    else
    {
      int i=0;
      Macroblock*     prevMB;
      int cod_counter;
      
      if((FmoGetNextMBNr(img->current_mb_nr) == -1)&&(currMB->mb_type == 0))
        cod_counter = img->cod_counter -1;
      else
        cod_counter = img->cod_counter;
      for( i = cod_counter ; i > 0 ; i--)
      {
        prevMB = &img->mb_data[img->current_mb_nr-i];
        
        if(prevMB->send_index_predictor_skip)
          rate += writeMotionVectorPredictorSKIP(prevMB->best_predictor_for_skip);
      }
      if((FmoGetNextMBNr(img->current_mb_nr) == -1)&&(currMB->mb_type == 0))
        if(currMB->send_index_predictor_skip)
          rate += writeMotionVectorPredictorSKIP(currMB->best_predictor_for_skip);
    }
    
    return rate;
  }
  
#ifdef MB32X32_MVC
  extern Macroblock  MB32, MB64;

  int write_predictor_index_for_skip_mode32()
  {
    int            rate       = 0;
    if (input->symbol_mode == CABAC)
    { 
	    MB32.send_index_predictor_skip = send_index_for_skip_prediction(MB32.best_predictor_for_skip);
      if (MB32.send_index_predictor_skip == TRUE)
        rate += writeMotionVectorPredictorSKIP(MB32.best_predictor_for_skip);
    }
    return rate;
  }
  int write_predictor_index_for_skip_mode64()
  {
    int            rate       = 0;
    if (input->symbol_mode == CABAC)
    { 
	    MB64.send_index_predictor_skip = send_index_for_skip_prediction(MB64.best_predictor_for_skip);
      if (MB64.send_index_predictor_skip == TRUE)
        rate += writeMotionVectorPredictorSKIP(MB64.best_predictor_for_skip);
    }
    return rate;
  }
#endif
  /*! 
  *************************************************************************************
  * \brief
  *    
  *
  * \param 
  *    
  *
  * \return
  *    
  *************************************************************************************
  */
  
  int writeMotionVectorPredictorSKIP (int predictor)
  {
    int            rate       = 0;
    Macroblock*    currMB     = &img->mb_data[img->current_mb_nr];
    SyntaxElement* currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
    int*           bitCount   = currMB->bitcounter;
    Slice*         currSlice  = img->currentSlice;
    const int*     partMap    = assignSE2partition[input->partition_mode];
    DataPartition* dataPart;
    
    currSE->value1 = predictor;
    currSE->type = SE_MV_PREDICTOR;
    currSE->value2 = (mv_comp.nb_mode_for_skip / 2) + (mv_comp.nb_mode_for_skip % 2);     
    
    if (input->symbol_mode == UVLC)
    {
      dataPart = &(currSlice->partArr[partMap[currSE->type]]);
      currSE->len = u_v(currSE->value2, "predictor for the SKIP mode", currSE->value1,dataPart->bitstream);
    }
    else
    {
      currSE->writing = writeMB_Motion_predictor_CABAC;
      dataPart = &(currSlice->partArr[partMap[currSE->type]]);
      dataPart->writeSyntaxElement (currSE, dataPart);
    }
    
    bitCount[BITS_SKIP_PRED] += currSE->len;
    rate                     += currSE->len;
    currSE++;  
    currMB->currSEnr++;
    
    return rate;
  }
  
  
  /*! 
  *************************************************************************************
  * \brief
  *    
  *
  * \param 
  *    
  *
  * \return
  *    
  *************************************************************************************
  */
  int writeIndexForMotionVectorPrediction (int  refframe,  int  list_idx,
    int  mv_mode,
    int i,
    int j)
  {
    int            rate       = 0;
    
    Macroblock*    currMB     = &img->mb_data[img->current_mb_nr];
    int*           bitCount   = currMB->bitcounter;
    SyntaxElement* currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
    Slice*         currSlice  = img->currentSlice;
    const int*     partMap    = assignSE2partition[input->partition_mode];
    DataPartition* dataPart;
    
    int maxmode = 0;
    
    if (img->type == I_SLICE) 
      maxmode = 1;
    else if (img->type == P_SLICE)
      maxmode = mv_comp.nb_mode_for_mvp;
    else if (img->type == B_SLICE)
      maxmode = mv_comp.nb_mode_for_mvb;
    
    currSE->value1 = predModeMV[mv_mode][refframe][list_idx][j][i];  
    currSE->type = SE_MV_PREDICTOR;
    currSE->value2 = (maxmode / 2) + (maxmode % 2);     
    
    if (maxmode > 1)
    {
      if (input->symbol_mode == UVLC)
      {
        dataPart = &(currSlice->partArr[partMap[currSE->type]]);
        currSE->len = u_v(currSE->value2, "predictor of MV", currSE->value1,dataPart->bitstream);
      }
      else
      {
        currSE->writing = writeMB_Motion_predictor_CABAC;
        dataPart = &(currSlice->partArr[partMap[currSE->type]]);
        dataPart->writeSyntaxElement (currSE, dataPart);
      }
    }
    
    bitCount[BITS_MV_PRED] += currSE->len;
    rate                   += currSE->len;
    currSE++;  
    currMB->currSEnr++;
    
    return rate;
  }
  
  /*! 
  *************************************************************************************
  * \brief
  *    This function is used to arithmetically encode the Motion vector 
  *   predictor
  *
  * \param 
  *    
  *
  * \return
  *    
  *************************************************************************************
  */
  
  void writeMB_Motion_predictor_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
  {
    Macroblock          *currMB = &img->mb_data[img->current_mb_nr];
    int act_ctx = 0;
    MotionInfoContexts *ctx  = (img->currentSlice)->mot_ctx;
    int nb_symbol;
    int curr_symbol;
    
    if(se->value2 != 0)
    {
      nb_symbol = se->value2;
      curr_symbol = se->value1;
      
      while(nb_symbol!=0)
      {
        if(currMB->mb_type == 0)
          biari_encode_symbol(eep_dp, (signed short) (curr_symbol % 2), &ctx->mv_predictor_skip_contexts[act_ctx]);
        else 
          if(img->type == P_SLICE)
            biari_encode_symbol(eep_dp, (signed short) (curr_symbol % 2), &ctx->mv_predictor_mvp_contexts[act_ctx]);
          else if(img->type == B_SLICE)
            biari_encode_symbol(eep_dp, (signed short) (curr_symbol % 2), &ctx->mv_predictor_mvb_contexts[act_ctx]);
          
          nb_symbol--;
          curr_symbol = curr_symbol / 2;
      }
    }
  }
  
  
  /*! 
  *************************************************************************************
  * \brief
  *    Returns 0 if no motion vector is available for the previous frame
  * at macrobloc position (pos_j, pos_i), 1 otherwise
  * because of INTRA_MODE or out of frame (checked here)Writes motion vectors predictor for SKIP
  *
  * \param 
  *    
  *
  * \return
  *    
  *************************************************************************************
  */
  short collocated_mv_available(int pos_y, int pos_x, int list)
  {
    short return_val = TRUE;
    
    if (listX[list][0]->ref_idx[0][pos_y*4][pos_x*4] == -1)  // Collocated block in INTRA
    {
      return_val = FALSE;
    }
    
    return (return_val);
  }
  
  /*! 
  *************************************************************************************
  * \brief
  *    Writes motion vectors predictor for SKIP
  *
  * \param 
  *    
  *
  * \return
  *    
  *************************************************************************************
  */
  short collocated_mv_available_previous_B_frame(int pos_y, int pos_x, int list)
  {
    short return_val = TRUE;
    
    if (ref_idx_previous_B_frame[list][pos_y*4][pos_x*4] == -1)  // Collocated block in INTRA
      return_val = FALSE;
    
    return (return_val);
  }
  
  /*! 
  *************************************************************************************
  * \brief
  *    
  *
  * \param 
  *    
  *
  * \return
  *    
  *************************************************************************************
  */
  short send_index_for_skip_prediction(int best_mode_skip)
  {
    short mv_x_ref, mv_y_ref;
    short i;
    short return_val = FALSE;      // Assumes it can be guessed
    short nb_of_available = 0;        // Number of modes available
    
    mv_x_ref = mv_comp.mv_pred_skip[best_mode_skip][0];
    mv_y_ref = mv_comp.mv_pred_skip[best_mode_skip][1];
    
    for (i=0; i<mv_comp.nb_mode_for_skip; i++)
    {
      //Check if among all modes enabled, one single is AVAILABLE (all others are NOT_AVAILABLE)
      if (mv_comp.mv_pred_skip[i][0] != NOT_AVAILABLE) 
        nb_of_available++;
    }
    
    if (nb_of_available > 1)  //Otherwise does not need to be sent: only one is available
    {
      for (i=0; i<mv_comp.nb_mode_for_skip; i++)
      {
        if ((mv_comp.mv_pred_skip[i][0] != mv_x_ref) 
          || (mv_comp.mv_pred_skip[i][1] != mv_y_ref))
          
          return_val = TRUE;    // One of the component for one prediction is different
        // We can return 0, the mode cannot be guessed
      }
    }
    
    return (return_val);
  }
  
  /*! 
  *************************************************************************************
  * \brief
  *    Set motion vector predictor
  *
  * \param 
  *    
  *
  * \return
  *    
  *************************************************************************************
  */
  void SetMotionVectorPredictor_Competition (short  pmv[2],
    char   **refPic,
    short  ***tmp_mv,
    short  ref_frame,
    int    list,
    int    block_x,
    int    block_y,
    int    blockshape_x,
    int    blockshape_y)
  {
    
    int mb_x                 = 4*block_x;
    int mb_y                 = 4*block_y;
    short  pmv2[2];
    short predictor;
    
    int maxmode=0;
    int modeblock;
    int prediction_mode;
    int *predictors=NULL;
    
    modeblock = determine_mode_block(blockshape_x,blockshape_y);
    
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
    
    predictor = 0;
    
    for (prediction_mode = 0; prediction_mode<maxmode; prediction_mode++)
    {
      if (predictors[prediction_mode] == PRED_H264_MEDIAN)  
      {
        SetMotionVectorPredictor(pmv2,refPic,tmp_mv,ref_frame,list,block_x,block_y,blockshape_x,blockshape_y);
        
        if (predictor == 0 )
        { 
          pmv[0]=pmv2[0];
          pmv[1]=pmv2[1];
        }
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
        pmv2[0] = block_a.available  ? tmp_mv[block_a.pos_y][block_a.pos_x][0] : 0;
        pmv2[1] = block_a.available  ? tmp_mv[block_a.pos_y][block_a.pos_x][1] : 0;  
      }
      else if (predictors[prediction_mode] == PRED_B)  
      {
        PixelPos  block_b;
        
        getLuma4x4Neighbour(img->current_mb_nr, block_x, block_y,            0, -1, &block_b);
        pmv2[0] = block_b.available  ? tmp_mv[block_b.pos_y][block_b.pos_x][0] : 0;
        pmv2[1] = block_b.available  ? tmp_mv[block_b.pos_y][block_b.pos_x][1] : 0;
      }
      else if (predictors[prediction_mode] == PRED_C)
      {
        PixelPos block_c, block_d;
        
        getLuma4x4Neighbour(img->current_mb_nr, block_x, block_y, blockshape_x, -1, &block_c);
        getLuma4x4Neighbour(img->current_mb_nr, block_x, block_y,           -1, -1, &block_d);
        
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
        
        pmv2[0] = block_c.available  ? tmp_mv[block_c.pos_y][block_c.pos_x][0] : 0;
        pmv2[1] = block_c.available  ? tmp_mv[block_c.pos_y][block_c.pos_x][1] : 0;
      }
      else if (predictors[prediction_mode] == PRED_COLOCATED)
      {    
        int y=(int)img->current_mb_nr/(img->width/16);  // Vertical
        int x=(int)img->current_mb_nr%(img->width/16);  // Horizontal
        
        if (img->type == P_SLICE)
        {
          
          if (collocated_mv_available(y, x, LIST_0) == TRUE)
          {
            // The collocated is the top left MV of the macroblock, whatever the subpartition...
            // To be improved.
            pmv2[0] = listX[0][0]->mv[LIST_0][y*4+block_y][x*4+block_x][0] * (ref_frame+1) / (listX[0][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
            pmv2[1] = listX[0][0]->mv[LIST_0][y*4+block_y][x*4+block_x][1] * (ref_frame+1) / (listX[0][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
          }
          else
          {
            pmv2[0] = 0;//NOT_AVAILABLE;
            pmv2[1] = 0;//NOT_AVAILABLE;
          }
        } 
        else if (img->type == B_SLICE)
        {
          // <FTRD : Compatibility with hierarchical B slices
          //  if(input->HierarchicalCoding == 0)
          //    {
          //    SetMotionVectorPredictor_collocated_B_SLICE(pmv2, y, x, list, ref_frame, block_y, block_x);
          //    }
          //    else
          //    {
          SetMotionVectorPredictor_collocated_HB_SLICE(pmv2, y, x, list, ref_frame, block_y, block_x);
          //    }
          // FTRD>
        }
      }
      mv_predictors[block_x][block_y][modeblock][ref_frame][list][prediction_mode][0] = pmv2[0];
      mv_predictors[block_x][block_y][modeblock][ref_frame][list][prediction_mode][1] = pmv2[1];
  }
  
  // Check if the index for the motion vector prediction needs to be sent or not
  // Stored in predModeMV2
  send_index_for_mv_prediction(block_x, block_y, modeblock, ref_frame, list, maxmode);
}

#ifdef MB32X32_MVC
int ResetAvailabilityOfNeighbors32(int mb_x, int mb_y, int *reset_block_x, int *reset_block_y, int blockshape_x, int blockshape_y);
void SetMotionVectorPredictor32 (short  pmv[2],
                               char   **refPic,
                               short  ***tmp_mv,
                               short  ref_frame,
                               int    list,
                               int    block_x,
                               int    block_y,
                               int    blockshape_x,
                               int    blockshape_y,
                               int    mb_ext_level);

void SetMotionVectorPredictor_Competition32 (short  pmv[2],
    char   **refPic,
    short  ***tmp_mv,
    short  ref_frame,
    int    list,
    int    block_x,
    int    block_y,
    int    blockshape_x,
    int    blockshape_y,
    short  mb_ext_level) 
{
    
    int mb_x                 = 4*block_x;
    int mb_y                 = 4*block_y;
    short  pmv2[2];
    short predictor;
    
    int maxmode=0;
    int modeblock;
    int prediction_mode;
    int *predictors=NULL;
    
//===Rahul
    int reset_mb_nr, reset_block_x=block_x, reset_block_y=block_y;
    reset_mb_nr = ResetAvailabilityOfNeighbors32(mb_x, mb_y, &reset_block_x, &reset_block_y, blockshape_x, blockshape_y);

//    block_x = reset_block_x;
//    block_y = reset_block_y;
//  img->current_mb_nr will give the top-left MB number of a current 32x32 block
//  reset_mb_nr will give the actual MB number 
// block_x and block_y are 4x4 index inside a Big Block....
// for eg, in 32x32 MB block_x and block_y range from 0 to 7...
//===

    modeblock = determine_mode_block(blockshape_x >> mb_ext_level,blockshape_y >> mb_ext_level);

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

    predictor = 0;
    
    for (prediction_mode = 0; prediction_mode<maxmode; prediction_mode++)
    {
      if (predictors[prediction_mode] == PRED_H264_MEDIAN)  
      {
//        SetMotionVectorPredictor(pmv2,refPic,tmp_mv,ref_frame,list,block_x,block_y,blockshape_x,blockshape_y);
        SetMotionVectorPredictor32(pmv2,refPic,tmp_mv,ref_frame,list,reset_block_x,reset_block_y,blockshape_x,blockshape_y,mb_ext_level);

        if (predictor == 0 )
        { 
          pmv[0]=pmv2[0];
          pmv[1]=pmv2[1];
        }
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
        pmv2[0] = block_a.available  ? tmp_mv[block_a.pos_y][block_a.pos_x][0] : 0;
        pmv2[1] = block_a.available  ? tmp_mv[block_a.pos_y][block_a.pos_x][1] : 0;  
      }
      else if (predictors[prediction_mode] == PRED_B)  
      {
        PixelPos  block_b;
        
        getLuma4x4Neighbour(reset_mb_nr, reset_block_x, reset_block_y,            0, -1, &block_b);
        pmv2[0] = block_b.available  ? tmp_mv[block_b.pos_y][block_b.pos_x][0] : 0;
        pmv2[1] = block_b.available  ? tmp_mv[block_b.pos_y][block_b.pos_x][1] : 0;
      }
      else if (predictors[prediction_mode] == PRED_C)
      {
        PixelPos block_c, block_d;
        
        getLuma4x4Neighbour(reset_mb_nr, reset_block_x, reset_block_y, blockshape_x, -1, &block_c);
        getLuma4x4Neighbour(reset_mb_nr, reset_block_x, reset_block_y,           -1, -1, &block_d);
        
        if (mb_y > 0)
        {
          if (mb_x < (8<<mb_ext_level))  // first column of 8x8 blocks
          {
            if (mb_y==(8<<mb_ext_level))
            {
              if (blockshape_x == 16<<mb_ext_level)      block_c.available  = 0;
              else                                       block_c.available &= 1;
            }
            else
            {
              if (mb_x+blockshape_x != (8<<mb_ext_level))  block_c.available &= 1;
              else                                         block_c.available  = 0;
            }
          }
          else
          {
            if (mb_x+blockshape_x != 16<<mb_ext_level)   block_c.available &= 1;
            else                                         block_c.available  = 0;
          }
        }
        
        if (!block_c.available)
        {
          block_c=block_d;
        }
        
        pmv2[0] = block_c.available  ? tmp_mv[block_c.pos_y][block_c.pos_x][0] : 0;
        pmv2[1] = block_c.available  ? tmp_mv[block_c.pos_y][block_c.pos_x][1] : 0;
      }
      else if (predictors[prediction_mode] == PRED_COLOCATED)
      {    
        int y=(int)reset_mb_nr/(img->width/16);  // Vertical
        int x=(int)reset_mb_nr%(img->width/16);  // Horizontal
        int ytop=(int)img->current_mb_nr/(img->width/16);  // Vertical
        int xtop=(int)img->current_mb_nr%(img->width/16);  // Horizontal

        //Rahul---(x,y) is the 16x16MB location (even within Bigblocks)
        if (img->type == P_SLICE)
        {
          
          if (collocated_mv_available(y, x, LIST_0) == TRUE)
          {
            // The collocated is the top left MV of the macroblock, whatever the subpartition...
            // To be improved.
//          pmv2[0] = listX[0][0]->mv[LIST_0][y*4+block_y][x*4+block_x][0] * (ref_frame+1) / (listX[0][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
//          pmv2[1] = listX[0][0]->mv[LIST_0][y*4+block_y][x*4+block_x][1] * (ref_frame+1) / (listX[0][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
          pmv2[0] = listX[0][0]->mv[LIST_0][ytop*4+block_y][xtop*4+block_x][0] * (ref_frame+1) / (listX[0][0]->ref_idx[0][ytop*4+block_y][xtop*4+block_x] + 1);
          pmv2[1] = listX[0][0]->mv[LIST_0][ytop*4+block_y][xtop*4+block_x][1] * (ref_frame+1) / (listX[0][0]->ref_idx[0][ytop*4+block_y][xtop*4+block_x] + 1);
//            pmv2[0] = listX[0][0]->mv[LIST_0][y*4][x*4][0] * (ref_frame+1) / (listX[0][0]->ref_idx[0][y*4][x*4] + 1);
//            pmv2[1] = listX[0][0]->mv[LIST_0][y*4][x*4][1] * (ref_frame+1) / (listX[0][0]->ref_idx[0][y*4][x*4] + 1);
          }
          else
          {
            pmv2[0] = 0;//NOT_AVAILABLE;
            pmv2[1] = 0;//NOT_AVAILABLE;
          }
        } 
        else if (img->type == B_SLICE)
        {

//        SetMotionVectorPredictor_collocated_HB_SLICE(pmv2, y   , x   , list, ref_frame, block_y, block_x);
          SetMotionVectorPredictor_collocated_HB_SLICE(pmv2, ytop, xtop, list, ref_frame, block_y, block_x);
        }
      }
      mv_predictors[block_x][block_y][modeblock][ref_frame][list][prediction_mode][0] = pmv2[0];
      mv_predictors[block_x][block_y][modeblock][ref_frame][list][prediction_mode][1] = pmv2[1];

   }
  
  // Check if the index for the motion vector prediction needs to be sent or not
  // Stored in predModeMV2
  send_index_for_mv_prediction(block_x, block_y, modeblock, ref_frame, list, maxmode);
}

#if 0// original...not working
  void SetMotionVectorPredictor_Competition32 (short  pmv[2],
    char   **refPic,
    short  ***tmp_mv,
    short  ref_frame,
    int    list,
    int    block_x,
    int    block_y,
    int    blockshape_x,
    int    blockshape_y,
    short  mb_ext_level)
  {
    
    int mb_x                 = 4*block_x;
    int mb_y                 = 4*block_y;
    short  pmv2[2];
    short predictor;
    
    int maxmode=0;
    int modeblock;
    int prediction_mode;
    int *predictors=NULL;
    
//===Rahul
    int reset_mb_nr, reset_block_x=block_x, reset_block_y=block_y;
    reset_mb_nr = ResetAvailabilityOfNeighbors32(mb_x, mb_y, &reset_block_x, &reset_block_y, blockshape_x, blockshape_y);

    block_x = reset_block_x;
    block_y = reset_block_y;
//  img->current_mb_nr will give the top-left MB number of a current 32x32 block
//  reset_mb_nr will give the actual MB number 
//===

    modeblock = determine_mode_block(blockshape_x >> mb_ext_level,blockshape_y >> mb_ext_level);

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

    predictor = 0;
    
    for (prediction_mode = 0; prediction_mode<maxmode; prediction_mode++)
    {
      if (predictors[prediction_mode] == PRED_H264_MEDIAN)  
      {
//        SetMotionVectorPredictor(pmv2,refPic,tmp_mv,ref_frame,list,block_x,block_y,blockshape_x,blockshape_y);
        SetMotionVectorPredictor32(pmv2,refPic,tmp_mv,ref_frame,list,block_x,block_y,blockshape_x,blockshape_y,mb_ext_level);

        if (predictor == 0 )
        { 
          pmv[0]=pmv2[0];
          pmv[1]=pmv2[1];
        }
      }
      else if (predictors[prediction_mode] == PRED_ZERO)
      {
        pmv2[0] = 0;
        pmv2[1] = 0;
      }
      else if (predictors[prediction_mode] == PRED_A)  
      {  
        PixelPos block_a;
        getLuma4x4Neighbour(reset_mb_nr, block_x, block_y,           -1,  0, &block_a);
        pmv2[0] = block_a.available  ? tmp_mv[block_a.pos_y][block_a.pos_x][0] : 0;
        pmv2[1] = block_a.available  ? tmp_mv[block_a.pos_y][block_a.pos_x][1] : 0;  
      }
      else if (predictors[prediction_mode] == PRED_B)  
      {
        PixelPos  block_b;
        
        getLuma4x4Neighbour(reset_mb_nr, block_x, block_y,            0, -1, &block_b);
        pmv2[0] = block_b.available  ? tmp_mv[block_b.pos_y][block_b.pos_x][0] : 0;
        pmv2[1] = block_b.available  ? tmp_mv[block_b.pos_y][block_b.pos_x][1] : 0;
      }
      else if (predictors[prediction_mode] == PRED_C)
      {
        PixelPos block_c, block_d;
        
        getLuma4x4Neighbour(reset_mb_nr, block_x, block_y, blockshape_x, -1, &block_c);
        getLuma4x4Neighbour(reset_mb_nr, block_x, block_y,           -1, -1, &block_d);
        
        if (mb_y > 0)
        {
          if (mb_x < (8<<mb_ext_level))  // first column of 8x8 blocks
          {
            if (mb_y==(8<<mb_ext_level))
            {
              if (blockshape_x == 16<<mb_ext_level)      block_c.available  = 0;
              else                                       block_c.available &= 1;
            }
            else
            {
              if (mb_x+blockshape_x != (8<<mb_ext_level))  block_c.available &= 1;
              else                                         block_c.available  = 0;
            }
          }
          else
          {
            if (mb_x+blockshape_x != 16<<mb_ext_level)   block_c.available &= 1;
            else                                         block_c.available  = 0;
          }
        }
        
        if (!block_c.available)
        {
          block_c=block_d;
        }
        
        pmv2[0] = block_c.available  ? tmp_mv[block_c.pos_y][block_c.pos_x][0] : 0;
        pmv2[1] = block_c.available  ? tmp_mv[block_c.pos_y][block_c.pos_x][1] : 0;
      }
      else if (predictors[prediction_mode] == PRED_COLOCATED)
      {    
        int y=(int)reset_mb_nr/(img->width/16);  // Vertical
        int x=(int)reset_mb_nr%(img->width/16);  // Horizontal
        //Rahul---(x,y) is the 16x16MB location (even within Bigblocks)
        if (img->type == P_SLICE)
        {
          
          if (collocated_mv_available(y, x, LIST_0) == TRUE)
          {
            // The collocated is the top left MV of the macroblock, whatever the subpartition...
            // To be improved.
//          pmv2[0] = listX[0][0]->mv[LIST_0][y*4+block_y][x*4+block_x][0] * (ref_frame+1) / (listX[0][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
//          pmv2[1] = listX[0][0]->mv[LIST_0][y*4+block_y][x*4+block_x][1] * (ref_frame+1) / (listX[0][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
            pmv2[0] = listX[0][0]->mv[LIST_0][y*4][x*4][0] * (ref_frame+1) / (listX[0][0]->ref_idx[0][y*4][x*4] + 1);
            pmv2[1] = listX[0][0]->mv[LIST_0][y*4][x*4][1] * (ref_frame+1) / (listX[0][0]->ref_idx[0][y*4][x*4] + 1);
          }
          else
          {
            pmv2[0] = 0;//NOT_AVAILABLE;
            pmv2[1] = 0;//NOT_AVAILABLE;
          }
        } 
        else if (img->type == B_SLICE)
        {
          SetMotionVectorPredictor_collocated_HB_SLICE(pmv2, y, x, list, ref_frame, block_y, block_x);
        }
      }
      mv_predictors[block_x][block_y][modeblock][ref_frame][list][prediction_mode][0] = pmv2[0];
      mv_predictors[block_x][block_y][modeblock][ref_frame][list][prediction_mode][1] = pmv2[1];

   }
  
  // Check if the index for the motion vector prediction needs to be sent or not
  // Stored in predModeMV2
  send_index_for_mv_prediction(block_x, block_y, modeblock, ref_frame, list, maxmode);
}

#endif
#endif

/*! 
*************************************************************************************
* \brief
*    
*
* \param 
*    
*
* \return
*    
*************************************************************************************
*/

void send_index_for_mv_prediction(int block_x, int block_y, int modeblock, short ref_frame, int list, int maxmode) 
{
  short egality_of_predictors = 0;
  short current_predictor_index;
  short nb_of_available = 0;
  short i;
  
  
  int pred0_x, pred0_y, current_pred_x, current_pred_y; 
  
  pred0_x = mv_predictors[block_x][block_y][modeblock][ref_frame][list][0][0];
  pred0_y = mv_predictors[block_x][block_y][modeblock][ref_frame][list][0][1];
  
  for (i=0; i<maxmode; i++)
  {
    //Check if among all modes enabled, one single is AVAILABLE (all others are NOT_AVAILABLE)
    if (mv_predictors[block_x][block_y][modeblock][ref_frame][list][i][0] != NOT_AVAILABLE) 
      nb_of_available++;
  }
  
  //  if (nb_of_available > 1)  //Otherwise does not need to be sent: only one is available
  {
    for(current_predictor_index=1;current_predictor_index<maxmode ;current_predictor_index++)
    {
      current_pred_x = mv_predictors[block_x][block_y][modeblock][ref_frame][list][current_predictor_index][0];
      current_pred_y = mv_predictors[block_x][block_y][modeblock][ref_frame][list][current_predictor_index][1];
      
      if ((current_pred_x == pred0_x) && (current_pred_y == pred0_y))
        egality_of_predictors++;  
    }
  }
  
  if(egality_of_predictors == maxmode-1)
    send_index_mv_prediction[modeblock][ref_frame][list][block_y][block_x][maxmode]=1;
  else 
    send_index_mv_prediction[modeblock][ref_frame][list][block_y][block_x][maxmode]=0;
}

/*! 
*************************************************************************************
* \brief
*    return le mode associé au blockshape x et y
*
* \param 
*    
*
* \return
*    
*************************************************************************************
*/

int determine_mode_block(int blockshape_x,int blockshape_y)
{
  short return_val=0;
  
  if((blockshape_x==16) && ( blockshape_y==16))
    return_val = 1;
  else if((blockshape_x==16) && ( blockshape_y==8))
    return_val = 2;
  else if((blockshape_x==8) && ( blockshape_y==16))
    return_val = 3;
  else if((blockshape_x==8) && ( blockshape_y==8))
    return_val = 4;
  else if((blockshape_x==8) && ( blockshape_y==4))
    return_val = 5;
  else if((blockshape_x==4) && ( blockshape_y==8))
    return_val =  6;
  else if((blockshape_x==4) && ( blockshape_y==4))
    return_val = 7;
  
  return (return_val);
}

/*! 
*************************************************************************************
* \brief
*    return  min cost for all prediction of MV 
*
* \param 
*    
*
* \return
*    
*************************************************************************************
*/

int compute_mv_cost_for_each_predictor(int f, int s, int cx, int cy, int px, int py, int ref, int list, int modeblock, int mb_x, int mb_y)  
{
  int min=MAX_VALUE;
  int i;
  
  
  int val, p;  
  
  int block_x = mb_x/4;
  int block_y = mb_y/4;
  int maxmode = 0;
  
  if (img->type == P_SLICE)
    maxmode = mv_comp.nb_mode_for_mvp;
  else if (img->type == B_SLICE)
    maxmode = mv_comp.nb_mode_for_mvb;
  
  for(i=0;i<maxmode;i++)
  {
    p = send_index_mv_prediction[modeblock][ref][list][block_y][block_x][maxmode];
    val = MV_COST(f,s,cx,cy,mv_predictors[block_x][block_y][modeblock][ref][list][i][0],mv_predictors[block_x][block_y][modeblock][ref][list][i][1])
      +    MVMODE_COST(f,i,p);
    
    if(val<min) 
      min = val;
    
  }
  return min ;
}

/*! 
*************************************************************************************
* \brief
*    
*
* \param 
*    
*
* \return
*    
*************************************************************************************
*/

int  MVMODE_COST(int f,int mode,int skip) 
{ 
  int nbbits=1;
  
  if(skip==1) nbbits=0;
  
  return nbbits*MVMODE_COST_FACTOR;
  
}

/*! 
*************************************************************************************
* \brief
*    get best prediction for motion vector 
*
* \param 
*    
*
* \return
*    
*************************************************************************************
*/

void GetBestPredVector (short  *pmv,short mv_x,short mv_y,int /*double*/ lambda,int ref, int list,int modeblock,int mb_x,int mb_y)
{
  //  Macroblock*    currMB     = &img->mb_data[img->current_mb_nr];
  //  int stop = 0;
  int current_predictor;
  
  int maxmode = 0;
  int current_cost;
  int nb_best_mode = 0;
  
  int best_cost = compute_mv_cost_for_each_predictor(lambda,0,mv_x,mv_y,0,0,ref,list,modeblock,mb_x,mb_y);
  
  int block_x = mb_x/4;
  int block_y = mb_y/4;
  
  if (img->type == I_SLICE) 
    maxmode = 1;
  else if (img->type == P_SLICE)
    maxmode = mv_comp.nb_mode_for_mvp;
  else if (img->type == B_SLICE)
    maxmode = mv_comp.nb_mode_for_mvb;
  
  for(current_predictor=0; current_predictor<maxmode; current_predictor++)
  {       
    current_cost = (MV_COST(lambda,0,mv_x,mv_y,mv_predictors[block_x][block_y][modeblock][ref][list][current_predictor][0],mv_predictors[block_x][block_y][modeblock][ref][list][current_predictor][1]))+
      MVMODE_COST(lambda,current_predictor,send_index_mv_prediction[modeblock][ref][list][block_y][block_x][maxmode]);
    
    if (current_cost == best_cost)
    {
      if(nb_best_mode == 0)
      {
        predModeMV[modeblock][ref][list][block_y][block_x]=current_predictor;      
        pmv[0]=(short)mv_predictors[block_x][block_y][modeblock][ref][list][current_predictor][0];         
        pmv[1]=(short)mv_predictors[block_x][block_y][modeblock][ref][list][current_predictor][1];         
      }
      send_index_mv_prediction[modeblock][ref][list][block_y][block_x][current_predictor]=1; 
      nb_best_mode++;
    }
    else 
    {
      send_index_mv_prediction[modeblock][ref][list][block_y][block_x][current_predictor]=0;
    }
  }
  
}

/*! 
*************************************************************************************
* \brief
*    
*
* \param 
*    
*
* \return
*    
*************************************************************************************
*/

void SetMotionVectorPredictor_collocated_B_SLICE (short  pmv[2], int y, int x, int list, int ref_frame, int block_y, int block_x)
{
  int colocated_MV_x, colocated_MV_y;
  
  pmv[0] = 0;
  pmv[1] = 0;
  
  if(list==LIST_0)
  {
    if(img->b_frame_to_code==1)//if this B frame is the first B frame of the GOP
    {
      if (collocated_mv_available(y, x, LIST_1) == TRUE)
      {
        colocated_MV_x = listX[LIST_1][0]->mv[0][y*4+block_y][x*4+block_x][0] * (ref_frame+1) / (listX[LIST_1][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
        colocated_MV_y = listX[LIST_1][0]->mv[0][y*4+block_y][x*4+block_x][1] * (ref_frame+1) / (listX[LIST_1][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
        pmv[0] = (colocated_MV_x/(input->successive_Bframe+1))+(ref_frame)*colocated_MV_x;
        pmv[1] = (colocated_MV_y/(input->successive_Bframe+1))+(ref_frame)*colocated_MV_y;
      }
      else 
      {
        pmv[0] = 0;//NOT_AVAILABLE;
        pmv[1] = 0;//NOT_AVAILABLE;
      }
    }
    else// not the first B frame: use the  motion vector field of previous B frame
    {
      if (collocated_mv_available_previous_B_frame(y, x, LIST_1) == TRUE)
      {
        colocated_MV_x = mv_previous_B_frame[LIST_1][y*4+block_y][x*4+block_x][0];// / (ref_idx_previous_B_frame[LIST_1][y*4+block_y][x*4+block_x] + 1);
        colocated_MV_y = mv_previous_B_frame[LIST_1][y*4+block_y][x*4+block_x][1];// / (ref_idx_previous_B_frame[LIST_1][y*4+block_y][x*4+block_x] + 1);
        
        pmv[0] = -(img->b_frame_to_code)*colocated_MV_x / (input->successive_Bframe+1-(img->b_frame_to_code-1));
        pmv[1] = -(img->b_frame_to_code)*colocated_MV_y / (input->successive_Bframe+1-(img->b_frame_to_code-1));
        
        if ((collocated_mv_available(y, x, LIST_1) == TRUE) && (ref_frame != 0))
        {
          pmv[0] = pmv[0] + (ref_frame*listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][0] / (listX[LIST_1][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1));
          pmv[1] = pmv[1] + (ref_frame*listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][1] / (listX[LIST_1][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1));
        }
      }
      else 
      {
        if (collocated_mv_available_previous_B_frame(y, x, LIST_0) == TRUE)
        {
          colocated_MV_x = (img->b_frame_to_code-1) * mv_previous_B_frame[LIST_0][y*4+block_y][x*4+block_x][0] / (ref_idx_previous_B_frame[LIST_0][y*4+block_y][x*4+block_x] * (input->successive_Bframe+1) + (img->b_frame_to_code-1));
          colocated_MV_y = (img->b_frame_to_code-1) * mv_previous_B_frame[LIST_0][y*4+block_y][x*4+block_x][1] / (ref_idx_previous_B_frame[LIST_0][y*4+block_y][x*4+block_x] * (input->successive_Bframe+1) + (img->b_frame_to_code-1));
          
          pmv[0]=(img->b_frame_to_code)*colocated_MV_x / (img->b_frame_to_code-1);  
          pmv[1]=(img->b_frame_to_code)*colocated_MV_y / (img->b_frame_to_code-1);
          
          if ((collocated_mv_available(y, x, LIST_1) == TRUE) && (ref_frame != 0))
          {
            pmv[0] = pmv[0] + (ref_frame*listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][0] / (listX[LIST_1][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1));
            pmv[1] = pmv[1] + (ref_frame*listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][1] / (listX[LIST_1][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1));
          }
        }
        else
        {
          pmv[0] = 0;//NOT_AVAILABLE;
          pmv[1] = 0;//NOT_AVAILABLE;
        }
      }
    }
    
  }
  else
  {//LIST_1 Refference frame multiple are not implemented for LIST_1
    if(img->b_frame_to_code==1)//if this B frame is the first B frame of the GOP
    {
      if (collocated_mv_available(y, x, LIST_1) == TRUE)
      {
        colocated_MV_x = listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][0] * (ref_frame+1) / (listX[LIST_1][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
        colocated_MV_y = listX[LIST_1][0]->mv[LIST_0][y*4+block_y][x*4+block_x][1] * (ref_frame+1) / (listX[LIST_1][0]->ref_idx[0][y*4+block_y][x*4+block_x] + 1);
        pmv[0] =  -input->successive_Bframe*colocated_MV_x/(input->successive_Bframe+1);
        pmv[1] =  -input->successive_Bframe*colocated_MV_y/(input->successive_Bframe+1);
      }
      else 
      {
        pmv[0] = 0;//NOT_AVAILABLE;
        pmv[1] = 0;//NOT_AVAILABLE;
      }
    }
    else// not the first B frame: use the  motion vector field of previous B frame
    {
      if (collocated_mv_available_previous_B_frame(y, x, LIST_1) == TRUE)
      {
        colocated_MV_x = mv_previous_B_frame[LIST_1][y*4+block_y][x*4+block_x][0];//  / (ref_idx_previous_B_frame[LIST_1][y*4+block_y][x*4+block_x] + 1);
        colocated_MV_y = mv_previous_B_frame[LIST_1][y*4+block_y][x*4+block_x][1];//  / (ref_idx_previous_B_frame[LIST_1][y*4+block_y][x*4+block_x] + 1);
        
        pmv[0]=(input->successive_Bframe+1-img->b_frame_to_code)*colocated_MV_x / (input->successive_Bframe+1-(img->b_frame_to_code-1));
        pmv[1]=(input->successive_Bframe+1-img->b_frame_to_code)*colocated_MV_y / (input->successive_Bframe+1-(img->b_frame_to_code-1));
        
      }
      else
      {
        if (collocated_mv_available_previous_B_frame(y, x, LIST_0) == TRUE)
        {
          colocated_MV_x = (img->b_frame_to_code-1) * mv_previous_B_frame[LIST_0][y*4+block_y][x*4+block_x][0] / (ref_idx_previous_B_frame[LIST_0][y*4+block_y][x*4+block_x] * (input->successive_Bframe+1) + (img->b_frame_to_code-1));
          colocated_MV_y = (img->b_frame_to_code-1) * mv_previous_B_frame[LIST_0][y*4+block_y][x*4+block_x][1] / (ref_idx_previous_B_frame[LIST_0][y*4+block_y][x*4+block_x] * (input->successive_Bframe+1) + (img->b_frame_to_code-1));
          
          pmv[0]=-(input->successive_Bframe+1-img->b_frame_to_code)*colocated_MV_x / (img->b_frame_to_code-1);
          pmv[1]=-(input->successive_Bframe+1-img->b_frame_to_code)*colocated_MV_y / (img->b_frame_to_code-1);            
        }
        else
        {
          pmv[0] = 0;//NOT_AVAILABLE;
          pmv[1] = 0;//NOT_AVAILABLE;
        }
      }
    }
    
    
  }
}

/*! 
*************************************************************************************
* \brief
*    
*
* \param 
*    
*
* \return
*    
*************************************************************************************
*/

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
*************************************************************************************
* \brief
*    
*
* \param 
*    
*
* \return
*    
*************************************************************************************
*/

// <FTRD Compatibility with APEC
void get_mem_pred_mv_indice (int ****** send_index_mv_prediction_RDDATA, int ****** predModeMV_RDDATA)
{
  int i, j, k, l, m;
  
  if ((*send_index_mv_prediction_RDDATA = (int*****)calloc(MAXMODE,sizeof(int****)))== NULL)
    no_mem_exit ("get_mem_pred_mv_indice: send_index_mv_prediction_RDDATA");
  for (k=0; k<MAXMODE; k++)
  {
    if (((*send_index_mv_prediction_RDDATA)[k] = (int****)calloc(img->num_ref_frames,sizeof(int***)) )== NULL)
      no_mem_exit ("get_mem_pred_mv_indice: send_index_mv_prediction_RDDATA");
    for (i=0; i<img->num_ref_frames; i++)
    { 
      if (((*send_index_mv_prediction_RDDATA)[k][i] = (int***)calloc((2),sizeof(int**)))== NULL)
        no_mem_exit ("get_mem_pred_mv_indice: send_index_mv_prediction_RDDATA");
      for (m=0; m<2; m++)
      { 
        if (((*send_index_mv_prediction_RDDATA)[k][i][m] = (int**)calloc(4,sizeof(int*)))== NULL)
          no_mem_exit ("get_mem_pred_mv_indice: send_index_mv_prediction_RDDATA");
        for(l=0; l<4; l++)
        {
          if (((*send_index_mv_prediction_RDDATA)[k][i][m][l] = (int*)calloc(4,sizeof(int)))== NULL)
            no_mem_exit ("get_mem_pred_mv_indice: send_index_mv_prediction_RDDATA");
            /*     for(j=0; j<4; j++)
            if (((*send_index_mv_prediction_RDDATA)[k][i][m][l][j] = (int*)calloc((max(mv_comp.nb_mode_for_mvp,mv_comp.nb_mode_for_mvb)+1),sizeof(int)))== NULL)
          no_mem_exit ("get_mem_pred_mv_indice: send_index_mv_prediction_RDDATA");*/
        }
      } 
    }
  }
  
  
  if ((*predModeMV_RDDATA= (int*****)calloc(MAXMODE,sizeof(int****)))== NULL)
    no_mem_exit ("get_mem_pred_mv_indice: predModeMV_RDDATA");
  for(j=0; j<(MAXMODE); j++)
  {
    if (((*predModeMV_RDDATA)[j]= (int****)calloc(img->num_ref_frames,sizeof(int***)))== NULL)
      no_mem_exit ("get_mem_pred_mv_indice: predModeMV_RDDATA");
    for (i=0; i<img->num_ref_frames; i++)
    {
      if (((*predModeMV_RDDATA)[j][i] = (int***)calloc(2,sizeof(int**)))== NULL)
        no_mem_exit ("get_mem_pred_mv_indice: predModeMV_RDDATA");
      for (m=0; m<2; m++)
      {
        if (((*predModeMV_RDDATA)[j][i][m] = (int**)calloc(4,sizeof(int*)))== NULL)
          no_mem_exit ("get_mem_pred_mv_indice: predModeMV_RDDATA");
        for (k=0; k<4; k++)
        {
          if (((*predModeMV_RDDATA)[j][i][m][k] = (int*)calloc(4,sizeof(int)))== NULL)
            no_mem_exit ("get_mem_pred_mv_indice: predModeMV_RDDATA");
        }
      }
    }
  }
  
}

/*! 
*************************************************************************************
* \brief
*    
*
* \param 
*    
*
* \return
*    
*************************************************************************************
*/

void free_mem_pred_mv_indice (int ***** send_index_mv_prediction_RDDATA, int ***** predModeMV_RDDATA)
{
  int i, j, k, l, m;
  
  for (k=0; k<MAXMODE; k++)
  {
    for (i=0; i<img->num_ref_frames; i++)
    { 
      for (m=0; m<2; m++)
      { 
        for(l=0; l<4; l++)
        {
        /* for(j=0; j<4; j++)
        {     
        free(send_index_mv_prediction_RDDATA[k][i][m][l][j]);
        }*/
          free(send_index_mv_prediction_RDDATA[k][i][m][l]);
        }
        free(send_index_mv_prediction_RDDATA[k][i][m]);
      }
      free(send_index_mv_prediction_RDDATA[k][i]);
    }
    free(send_index_mv_prediction_RDDATA[k]);
  }
  free(send_index_mv_prediction_RDDATA);
  
  
  for(j=0; j<(MAXMODE); j++)
  {
    for (i=0; i<img->num_ref_frames; i++)
    {
      for (m=0; m<2; m++)
      {
        for (k=0; k<4; k++)
        {
          free(predModeMV_RDDATA[j][i][m][k]);
        }
        free(predModeMV_RDDATA[j][i][m]);
      }
      free(predModeMV_RDDATA[j][i]);
    }
    free(predModeMV_RDDATA[j]);
  }
  free(predModeMV_RDDATA);
  
}

// FTRD>

/*! 
*************************************************************************************
* \brief
*    
*
* \param 
*    
*
* \return
*    
*************************************************************************************
*/

// <FTRD : Compatibility with hierarchical B slices
void init_mv_scale_hb()
{
  int i, j, k, l;
  int iTRb = 0 ; int iTRp = 0;
  int prescale = 0;
  int num_ref_max = 0;
  
  // Mv in current B frame points to LIST_0 or LIST_1
  for(i = 0 ; i < 2 ; i++)
  {
    for(k = 0 ; k <listXsize[i] ; k++)
    {
      iTRb = Clip3(-128, 127, enc_picture->poc - listX[i][k]->poc);
      
      // The colocated vector points to LIST_0 or LIST_1
      for(j=0;j<2;j++)
      {
        if(listX[LIST_1][0]->slice_type == B_SLICE) 
        {
          if(j==0) num_ref_max = input->B_List0_refs;
          if(j==1) num_ref_max = input->B_List1_refs;
        }
        if(listX[LIST_1][0]->slice_type == P_SLICE) 
        {
          num_ref_max = input->num_ref_frames;
        }
        
        for(l = 0 ; l < num_ref_max ; l++)
        {
          iTRp = (int) Clip3(-128, 127, listX[LIST_1][0]->poc - (listX[LIST_1][0]->ref_pic_num[j][l]/2));
          
          if(iTRp==0) img->mvscale_hb[i][j][k][l]=9999;
          else 
          {
            prescale=(16384 + absm(iTRp/2))/iTRp;
            img->mvscale_hb[i][j][k][l] = Clip3(-2048, 2047,(iTRb*prescale + 32) >> 6);
          }
        }
      }
    }
  }
  
}
// FTRD>

/*! 
*************************************************************************************
* \brief
*    Set colocated motion vector for hirarchical B frames
*
* \param 
*    
*
* \return
*    
*************************************************************************************
*/

void SetMotionVectorPredictor_collocated_HB_SLICE (short  pmv[2], int y, int x, int list, int ref_frame, int block_y, int block_x)
{
  int colocated_MV_x = 0, colocated_MV_y = 0;
  int ref_frame_coloc = 0;
  pmv[0] = 0; pmv[1] = 0; 
  
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
      pmv[0] = (img->mvscale_hb[LIST_0][LIST_0][ref_frame][ref_frame_coloc] * colocated_MV_x + 128) >> 8;
      pmv[1] = (img->mvscale_hb[LIST_0][LIST_0][ref_frame][ref_frame_coloc] * colocated_MV_y + 128) >> 8;
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
        pmv[0] = (img->mvscale_hb[LIST_0][LIST_1][ref_frame][ref_frame_coloc] * colocated_MV_x + 128) >> 8;
        pmv[1] = (img->mvscale_hb[LIST_0][LIST_1][ref_frame][ref_frame_coloc] * colocated_MV_y + 128) >> 8;
      }
      
      // No colocated mv is available.
      else
      {
        pmv[0] = 0;
        pmv[1] = 0;
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
      pmv[0] = (img->mvscale_hb[LIST_1][LIST_0][ref_frame][ref_frame_coloc] * colocated_MV_x + 128) >> 8;
      pmv[1] = (img->mvscale_hb[LIST_1][LIST_0][ref_frame][ref_frame_coloc] * colocated_MV_y + 128) >> 8;
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
        pmv[0] = (img->mvscale_hb[LIST_1][LIST_1][ref_frame][ref_frame_coloc] * colocated_MV_x + 128) >> 8;
        pmv[1] = (img->mvscale_hb[LIST_1][LIST_1][ref_frame][ref_frame_coloc] * colocated_MV_y + 128) >> 8;
      }
      
      // No colocated mv is available.
      else
      {
        pmv[0] = 0;
        pmv[1] = 0;
      }
    }
  }
}

#endif
