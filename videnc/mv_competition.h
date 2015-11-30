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

#ifndef _MV_COMPETITION_H
#define _MV_COMPETITION_H

#include "global.h"
#include "mbuffer.h"

#define DISABLED      -1
#define ENABLED        0 
#define NOT_AVAILABLE    999

#define MVMODE_COST_FACTOR  1        

#define PRED_H264_MEDIAN       0
#define PRED_EXTENDEDSPATIAL   1
#define PRED_A                 2
#define PRED_B                 3
#define PRED_C                 4  
#define PRED_ZERO              5  
#define PRED_COLOCATED         6

typedef struct
{
  int nb_mode_for_skip;            // Number of active predictors for the skip mode 
  //this variable is in macroblock structure  int best_predictor_for_skip;  
  // Index of the best predictor for the skip mode
  int current_predictor_for_skip;  // Currently tested predictor for the skip mode
  int *predictor_for_skip;         // Contains the name of the active predictors (PRED_H264_MEDIAN, PRED_A, etc) 
  int **mv_pred_skip;              // [Mode_Skip][hv] Contains the motion vectors for each skip predictors
  int nb_mode_for_mvp;             // Number of active predictors for the motion vectors of P frames 
  int nb_mode_for_mvb;             // Number of active predictors for the motion vectors of B frames
  int *predictor_for_mvp;          // Contains the name of the active predictors for P frames 
  int *predictor_for_mvb;          // Contains the name of the active predictors for B frames 
} MV_Competition;


void  init_MV_Competition();
void  close_MV_Competition();
int   write_predictor_index_for_skip_mode();
int   writeMotionVectorPredictorSKIP (int);
void  writeMB_Motion_predictor_CABAC(SyntaxElement *, EncodingEnvironmentPtr);
short collocated_mv_available(int, int, int);
short send_index_for_skip_prediction(int);
void  send_index_for_mv_prediction(int, int, int, short, int, int); 
int   determine_mode_block(int,int);
int   compute_mv_cost_for_each_predictor(int, int, int, int, int, int, int, int, int, int, int); 
int   MVMODE_COST(int,int,int);
void  GetBestPredVector (short*, short, short, int, int, int, int, int, int);
int   writeIndexForMotionVectorPrediction (int, int, int, int, int);
void  SetMotionVectorPredictor_Competition (short  pmv[2], char   **, short  ***,short, int, int, int, int, int);
void  SetMotionVectorPredictor_collocated_B_SLICE (short  pmv[2], int, int, int, int, int, int);

#ifdef MB32X32_MVC
void  SetMotionVectorPredictor_Competition32 (short  pmv[2], char   **, short  ***,short, int, int, int, int, int, short  mb_ext_level);
int write_predictor_index_for_skip_mode32();
int write_predictor_index_for_skip_mode64();
#endif
// <FTRD : Compatibility with hierarchical B slices
void    init_mv_scale_hb();
void  SetMotionVectorPredictor_collocated_HB_SLICE (short  pmv[2], int, int, int, int, int, int);
// FTRD>

void  Copy_MV_B_frame(short ****  mv, char  ***   ref_idx);
short collocated_mv_available_previous_B_frame(int pos_y, int pos_x, int list);
short ****    mv_previous_B_frame;      //!< motion vector       [list][subblock_x][subblock_y][component]
char  ***     ref_idx_previous_B_frame; //!< reference picture   [list][subblock_y][subblock_x]
int   pass_with_writing; 

// <FTRD Compatibility with APEC
void get_mem_pred_mv_indice (int ****** send_index_mv_prediction_RDDATA, int ****** predModeMV_RDDATA);
void free_mem_pred_mv_indice (int ***** send_index_mv_prediction_RDDATA, int ***** predModeMV_RDDATA);
// FTRD>
#endif
