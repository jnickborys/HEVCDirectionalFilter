#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

// #define NDEBUG
#include <assert.h>
#include <memory.h>

#include "global.h"
#include "memalloc.h"
#include "malloc.h"
#include "defines.h"
#include "mbuffer.h"
#include "image.h"
#include "adaptive_filter.h"
#include "switched_filters.h"

#ifdef SWITCHED_FILTERS
static double AccErrorP[NUM_SIFO][16];                              // [Filter][Sppos]
static double SequenceAccErrorP[NUM_SIFO][16];                      // [Filter][Sppos]
static int SamplesSubPelOffsetP[16];                                // [Sppos]
static int SamplesFrameOffsetP[MAX_REFERENCE_PICTURES];             // [Frame]
static double AccFrameOffsetP[MAX_REFERENCE_PICTURES][16][NUM_SIFO];// [Frame][Sppos][Filter]

// Prediction error used in frame filter selection
static double AccErrorB[NUM_SIFO][NUM_SIFO][16][16];                // [FilterF][FilterB][SpposF][SpposB]
static int SamplesB[NUM_SIFO][NUM_SIFO][16][16];                    // [FilterF][FilterB][SpposF][SpposB]
static int BestCombFilterB[16];                                     // [Sppos]

// Prediction error used in Sequence filter selection
static double SequenceAccErrorB[NUM_SIFO][NUM_SIFO][16][16];        // [FilterF][FilterB][SpposF][SpposB]
static int SequenceBestCombFilterB[16];                             // [Sppos]

// Frame offsets
static int FrameOffset[2][MAX_REFERENCE_PICTURES];                  // [Lists][Frame]
static double AccFrameOffset[2][MAX_REFERENCE_PICTURES];            // [Lists][Frame]
static int SamplesFrameOffset[2][MAX_REFERENCE_PICTURES];           // [Lists][Frame]

// SubPel offsets
static int SubPelOffset[2][16];                                     // [Lists][Sppos]
static double AccSubPelOffset[2][16];                               // [Lists][Sppos]
static int SamplesSubPelOffset[2][16];                              // [Lists][Sppos]

static double SIFO_FILTER[NUM_SIFO][15][SQR_FILTER];

extern pic_parameter_set_rbsp_t *PicParSet[MAXPPS];
extern void GenerateFullPelRepresentation(pel_t ** Fourthpel, pel_t * Fullpel, int xsize, int ysize);
extern void frame_picture (Picture *frame, int method);

extern double STANDARD_2D_FILTER[15][SQR_FILTER];

//
//
//
#ifdef HPF_COMBO
double LCDIF[15][SQR_FILTER] = {
{
   0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0, 0.0/128.0,  
   0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0, 0.0/128.0,  
   3.0/128.0, -15.0/128.0, 111.0/128.0, 37.0/128.0, -10.0/128.0, 2.0/128.0, 
   0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0, 0.0/128.0,  
   0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0, 0.0/128.0,  
   0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0, 0.0/128.0
},  // a_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0, 0.0/128.0,  
  3.0/128.0, -17.0/128.0, 78.0/128.0, 78.0/128.0, -17.0/128.0, 3.0/128.0, 
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0, 0.0/128.0
},  // b_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,   0.0/128.0, 0.0/128.0,  
  2.0/128.0, -10.0/128.0, 37.0/128.0, 111.0/128.0, -15.0/128.0, 3.0/128.0, 
  0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,   0.0/128.0, 0.0/128.0
},  // c_pos
{
  0.0/128.0, 0.0/128.0,  -3.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0,  
  0.0/128.0, 0.0/128.0, -15.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0,  
  0.0/128.0, 0.0/128.0, 111.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0,  
  0.0/128.0, 0.0/128.0,  37.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0,  
  0.0/128.0, 0.0/128.0, -10.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0,  
  0.0/128.0, 0.0/128.0,   2.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0
},  // d_pos
{
  3.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0, -15.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0, 111.0/128.0,  0.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0,   0.0/128.0, 37.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0, -10.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0, 2.0/128.0
},  // e_pos
{
  3.0/256.0,   0.0/256.0,   0.0/256.0,   0.0/256.0,   0.0/256.0, 3.0/256.0,  
  0.0/256.0, -15.0/256.0,   0.0/256.0,   0.0/256.0, -15.0/256.0, 0.0/256.0,  
  0.0/256.0,   0.0/256.0, 111.0/256.0, 111.0/256.0,   0.0/256.0, 0.0/256.0,  
  0.0/256.0,   0.0/256.0,  37.0/256.0,  37.0/256.0,   0.0/256.0, 0.0/256.0,  
  0.0/256.0, -10.0/256.0,   0.0/256.0,   0.0/256.0, -10.0/256.0, 0.0/256.0,  
  2.0/256.0,   0.0/256.0,   0.0/256.0,   0.0/256.0,   0.0/256.0, 2.0/256.0
},  // f_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,   0.0/128.0, 3.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0, -15.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0, 111.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0, 37.0/128.0,   0.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0, -10.0/128.0,  0.0/128.0,   0.0/128.0,   0.0/128.0, 0.0/128.0,  
  2.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,   0.0/128.0, 0.0/128.0,  
},  // g_pos
{
  0.0/128.0, 0.0/128.0,   3.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0,  
  0.0/128.0, 0.0/128.0, -17.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0,  
  0.0/128.0, 0.0/128.0,  78.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0,  
  0.0/128.0, 0.0/128.0,  78.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0,  
  0.0/128.0, 0.0/128.0, -17.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0,  
  0.0/128.0, 0.0/128.0,   3.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0
},  // h_pos
{
  3.0/256.0,   0.0/256.0,   0.0/256.0,  0.0/256.0,   0.0/256.0, 2.0/256.0,  
  0.0/256.0, -15.0/256.0,   0.0/256.0,  0.0/256.0, -10.0/256.0, 0.0/256.0,  
  0.0/256.0,   0.0/256.0, 111.0/256.0, 37.0/256.0,   0.0/256.0, 0.0/256.0,  
  0.0/256.0,   0.0/256.0, 111.0/256.0, 37.0/256.0,   0.0/256.0, 0.0/256.0,  
  0.0/256.0, -15.0/256.0,   0.0/256.0,  0.0/256.0, -10.0/256.0, 0.0/256.0,  
  3.0/256.0,   0.0/256.0,   0.0/256.0,  0.0/256.0,   0.0/256.0, 2.0/256.0,  
},  // i_pos
{
  3.0/256.0,   0.0/256.0,  0.0/256.0,  0.0/256.0,   0.0/256.0, 3.0/256.0,  
  0.0/256.0, -17.0/256.0,  0.0/256.0,  0.0/256.0, -17.0/256.0, 0.0/256.0,  
  0.0/256.0,   0.0/256.0, 78.0/256.0, 78.0/256.0,   0.0/256.0, 0.0/256.0,  
  0.0/256.0,   0.0/256.0, 78.0/256.0, 78.0/256.0,   0.0/256.0, 0.0/256.0,  
  0.0/256.0, -17.0/256.0,  0.0/256.0,  0.0/256.0, -17.0/256.0, 0.0/256.0,  
  3.0/256.0,   0.0/256.0,  0.0/256.0,  0.0/256.0,   0.0/256.0, 3.0/256.0,  
},  // j_pos
{
  2.0/256.0,   0.0/256.0,  0.0/256.0,   0.0/256.0,   0.0/256.0, 3.0/256.0,  
  0.0/256.0, -10.0/256.0,  0.0/256.0,   0.0/256.0, -15.0/256.0, 0.0/256.0,  
  0.0/256.0,   0.0/256.0, 37.0/256.0, 111.0/256.0,   0.0/256.0, 0.0/256.0,  
  0.0/256.0,   0.0/256.0, 37.0/256.0, 111.0/256.0,   0.0/256.0, 0.0/256.0,  
  0.0/256.0, -10.0/256.0,  0.0/256.0,   0.0/256.0, -15.0/256.0, 0.0/256.0,  
  2.0/256.0,   0.0/256.0,  0.0/256.0,   0.0/256.0,   0.0/256.0, 3.0/256.0,  
},  // k_pos
{
  0.0/128.0, 0.0/128.0,   2.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0,  
  0.0/128.0, 0.0/128.0, -10.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0,  
  0.0/128.0, 0.0/128.0,  37.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0,  
  0.0/128.0, 0.0/128.0, 111.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0,  
  0.0/128.0, 0.0/128.0, -15.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0,  
  0.0/128.0, 0.0/128.0,   3.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0
},  // l_pos
{
  0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0, 2.0/128.0,  
  0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0, -10.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0,   0.0/128.0, 37.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0, 111.0/128.0,  0.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0, -15.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0, 0.0/128.0,  
  3.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0, 0.0/128.0,  
},  // m_pos
{
  2.0/256.0,   0.0/256.0,   0.0/256.0,   0.0/256.0,   0.0/256.0, 2.0/256.0,  
  0.0/256.0, -10.0/256.0,   0.0/256.0,   0.0/256.0, -10.0/256.0, 0.0/256.0,  
  0.0/256.0,   0.0/256.0,  37.0/256.0,  37.0/256.0,   0.0/256.0, 0.0/256.0,  
  0.0/256.0,   0.0/256.0, 111.0/256.0, 111.0/256.0,   0.0/256.0, 0.0/256.0,  
  0.0/256.0, -15.0/256.0,   0.0/256.0,   0.0/256.0, -15.0/256.0, 0.0/256.0,  
  3.0/256.0,   0.0/256.0,   0.0/256.0,   0.0/256.0,   0.0/256.0, 3.0/256.0,  
},  // n_pos
{
  2.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0, -10.0/128.0,  0.0/128.0,   0.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0, 37.0/128.0,   0.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0, 111.0/128.0,   0.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0, -15.0/128.0, 0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,   0.0/128.0, 3.0/128.0,  
}  // o_pos
};
#endif

//
//
//
double SYMMETRIC_1[15][SQR_FILTER] = {
{
   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
  -1.0/128.0, -2.0/128.0, 97.0/128.0, 42.0/128.0, -11.0/128.0,  3.0/128.0, 
   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0
},  // a_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
  3.0/128.0, -16.0/128.0, 77.0/128.0, 77.0/128.0, -16.0/128.0,  3.0/128.0, 
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0
},  // b_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  3.0/128.0, -11.0/128.0, 42.0/128.0, 97.0/128.0, -2.0/128.0, -1.0/128.0, 
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0
},  // c_pos
{
  0.0/128.0,  0.0/128.0,   -1.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,   -2.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,   97.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,   42.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,  -11.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,    3.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0
},  // d_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  -7.0/128.0,  4.0/128.0, -2.0/128.0, -3.0/128.0,  0.0/128.0,  
  0.0/128.0,   4.0/128.0, 66.0/128.0, 36.0/128.0,  1.0/128.0,  0.0/128.0,  
  0.0/128.0,  -2.0/128.0, 36.0/128.0,  7.0/128.0, -6.0/128.0,  0.0/128.0,  
  0.0/128.0,  -3.0/128.0,  1.0/128.0, -6.0/128.0,  2.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0
},  // e_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  -4.0/128.0, -1.0/128.0, -1.0/128.0, -4.0/128.0,  0.0/128.0,  
  0.0/128.0,  -2.0/128.0, 56.0/128.0, 56.0/128.0, -2.0/128.0,  0.0/128.0,  
  0.0/128.0,  -7.0/128.0, 25.0/128.0, 25.0/128.0, -7.0/128.0,  0.0/128.0,  
  0.0/128.0,   1.0/128.0, -4.0/128.0, -4.0/128.0,  1.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0
},  // f_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  -3.0/128.0, -2.0/128.0,  4.0/128.0, -7.0/128.0,  0.0/128.0,  
  0.0/128.0,   1.0/128.0, 36.0/128.0, 66.0/128.0,  4.0/128.0,  0.0/128.0,  
  0.0/128.0,  -6.0/128.0,  7.0/128.0, 36.0/128.0, -2.0/128.0,  0.0/128.0,  
  0.0/128.0,   2.0/128.0, -6.0/128.0,  1.0/128.0, -3.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0
},  // g_pos
{
  0.0/128.0,  0.0/128.0,    3.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,  -16.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,   77.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,   77.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,  -16.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,    3.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0
},  // h_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  -4.0/128.0, -2.0/128.0, -7.0/128.0,  1.0/128.0,  0.0/128.0,  
  0.0/128.0,  -1.0/128.0, 56.0/128.0, 25.0/128.0, -4.0/128.0,  0.0/128.0,  
  0.0/128.0,  -1.0/128.0, 56.0/128.0, 25.0/128.0, -4.0/128.0,  0.0/128.0,  
  0.0/128.0,  -4.0/128.0, -2.0/128.0, -7.0/128.0,  1.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0
},  // i_pos
{
  0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  -6.0/128.0, -6.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  -6.0/128.0,  44.0/128.0, 44.0/128.0, -6.0/128.0,  0.0/128.0,  
  0.0/128.0,  -6.0/128.0,  44.0/128.0, 44.0/128.0, -6.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  -6.0/128.0, -6.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0
},  // j_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,   1.0/128.0, -7.0/128.0, -2.0/128.0, -4.0/128.0,  0.0/128.0,  
  0.0/128.0,  -4.0/128.0, 25.0/128.0, 56.0/128.0, -1.0/128.0,  0.0/128.0,  
  0.0/128.0,  -4.0/128.0, 25.0/128.0, 56.0/128.0, -1.0/128.0,  0.0/128.0,  
  0.0/128.0,   1.0/128.0, -7.0/128.0, -2.0/128.0, -4.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0
},  // k_pos
{
  0.0/128.0,  0.0/128.0,    3.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,  -11.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,   42.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,   97.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,   -2.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,   -1.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0
},  // l_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  -3.0/128.0,  1.0/128.0, -6.0/128.0,  2.0/128.0,  0.0/128.0,  
  0.0/128.0,  -2.0/128.0, 36.0/128.0,  7.0/128.0, -6.0/128.0,  0.0/128.0,  
  0.0/128.0,   4.0/128.0, 66.0/128.0, 36.0/128.0,  1.0/128.0,  0.0/128.0,  
  0.0/128.0,  -7.0/128.0,  4.0/128.0, -2.0/128.0, -3.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0
},  // m_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,   1.0/128.0, -4.0/128.0, -4.0/128.0,  1.0/128.0,  0.0/128.0,  
  0.0/128.0,  -7.0/128.0, 25.0/128.0, 25.0/128.0, -7.0/128.0,  0.0/128.0,  
  0.0/128.0,  -2.0/128.0, 56.0/128.0, 56.0/128.0, -2.0/128.0,  0.0/128.0,  
  0.0/128.0,  -4.0/128.0, -1.0/128.0, -1.0/128.0, -4.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0
},  // n_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,   2.0/128.0, -6.0/128.0,  1.0/128.0, -3.0/128.0,  0.0/128.0,  
  0.0/128.0,  -6.0/128.0,  7.0/128.0, 36.0/128.0, -2.0/128.0,  0.0/128.0,  
  0.0/128.0,   1.0/128.0, 36.0/128.0, 66.0/128.0,  4.0/128.0,  0.0/128.0,  
  0.0/128.0,  -3.0/128.0, -2.0/128.0,  4.0/128.0, -7.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0
}  // o_pos
};


//
//
//
double SYMMETRIC_2[15][SQR_FILTER] = {
{
  0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
  4.0/128.0, -14.0/128.0, 108.0/128.0, 38.0/128.0, -12.0/128.0,  4.0/128.0, 
  0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0
},  // a_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
  6.0/128.0, -23.0/128.0, 81.0/128.0, 81.0/128.0, -23.0/128.0,  6.0/128.0, 
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0
},  // b_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  
  4.0/128.0, -12.0/128.0, 38.0/128.0, 108.0/128.0, -14.0/128.0, 4.0/128.0, 
  0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0
},  // c_pos
{
  0.0/128.0,  0.0/128.0,    4.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,  -14.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,  108.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,   38.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,  -12.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,    4.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0
},  // d_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  -6.0/128.0,  2.0/128.0, -4.0/128.0, -3.0/128.0,  0.0/128.0,  
  0.0/128.0,   2.0/128.0, 71.0/128.0, 38.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  -4.0/128.0, 38.0/128.0,  7.0/128.0, -7.0/128.0,  0.0/128.0,  
  0.0/128.0,  -3.0/128.0,  0.0/128.0, -7.0/128.0,  4.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0
},  // e_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  -3.0/128.0, -2.0/128.0, -2.0/128.0, -3.0/128.0,  0.0/128.0,  
  0.0/128.0,  -3.0/128.0, 58.0/128.0, 58.0/128.0, -3.0/128.0,  0.0/128.0,  
  0.0/128.0,  -9.0/128.0, 26.0/128.0, 26.0/128.0, -9.0/128.0,  0.0/128.0,  
  0.0/128.0,   2.0/128.0, -5.0/128.0, -5.0/128.0,  2.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0
},  // f_pos
{
  0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
  0.0/128.0,  -3.0/128.0,  -4.0/128.0,  2.0/128.0,  -6.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  38.0/128.0, 71.0/128.0,   2.0/128.0,  0.0/128.0,  
  0.0/128.0,  -7.0/128.0,   7.0/128.0, 38.0/128.0,  -4.0/128.0,  0.0/128.0,  
  0.0/128.0,   4.0/128.0,  -7.0/128.0,  0.0/128.0,  -3.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0
},  // g_pos
{
  0.0/128.0,  0.0/128.0,    6.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,  -23.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,   81.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,   81.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,  -23.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,    6.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0
},  // h_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  -3.0/128.0, -3.0/128.0, -9.0/128.0,  2.0/128.0,  0.0/128.0,  
  0.0/128.0,  -2.0/128.0, 58.0/128.0, 26.0/128.0, -5.0/128.0,  0.0/128.0,  
  0.0/128.0,  -2.0/128.0, 58.0/128.0, 26.0/128.0, -5.0/128.0,  0.0/128.0,  
  0.0/128.0,  -3.0/128.0, -3.0/128.0, -9.0/128.0,  2.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0
},  // i_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,   1.0/128.0, -7.0/128.0, -7.0/128.0,  1.0/128.0,  0.0/128.0,  
  0.0/128.0,  -7.0/128.0, 45.0/128.0, 45.0/128.0, -7.0/128.0,  0.0/128.0,  
  0.0/128.0,  -7.0/128.0, 45.0/128.0, 45.0/128.0, -7.0/128.0,  0.0/128.0,  
  0.0/128.0,   1.0/128.0, -7.0/128.0, -7.0/128.0,  1.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0
},  // j_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,   2.0/128.0, -9.0/128.0, -3.0/128.0, -3.0/128.0,  0.0/128.0,  
  0.0/128.0,  -5.0/128.0, 26.0/128.0, 58.0/128.0, -2.0/128.0,  0.0/128.0,  
  0.0/128.0,  -5.0/128.0, 26.0/128.0, 58.0/128.0, -2.0/128.0,  0.0/128.0,  
  0.0/128.0,   2.0/128.0, -9.0/128.0, -3.0/128.0, -3.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0
},  // k_pos
{
  0.0/128.0,  0.0/128.0,    4.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,  -12.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,   38.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,  108.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,  -14.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  0.0/128.0,    4.0/128.0, 0.0/128.0,  0.0/128.0,  0.0/128.0
},  // l_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  -3.0/128.0,  0.0/128.0,  -7.0/128.0,  4.0/128.0,  0.0/128.0,  
  0.0/128.0,  -4.0/128.0, 38.0/128.0,   7.0/128.0, -7.0/128.0,  0.0/128.0,  
  0.0/128.0,   2.0/128.0, 71.0/128.0,  38.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,  -6.0/128.0,  2.0/128.0,  -4.0/128.0, -3.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0
},  // m_pos
{
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  
  0.0/128.0,   2.0/128.0, -5.0/128.0, -5.0/128.0,  2.0/128.0,  0.0/128.0,  
  0.0/128.0,  -9.0/128.0, 26.0/128.0, 26.0/128.0, -9.0/128.0,  0.0/128.0,  
  0.0/128.0,  -3.0/128.0, 58.0/128.0, 58.0/128.0, -3.0/128.0,  0.0/128.0,  
  0.0/128.0,  -3.0/128.0, -2.0/128.0, -2.0/128.0, -3.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0,  0.0/128.0
},  // n_pos
{
  0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0,  
  0.0/128.0,   4.0/128.0,  -7.0/128.0,  0.0/128.0,  -3.0/128.0,  0.0/128.0,  
  0.0/128.0,  -7.0/128.0,   7.0/128.0, 38.0/128.0,  -4.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,  38.0/128.0, 71.0/128.0,   2.0/128.0,  0.0/128.0,  
  0.0/128.0,  -3.0/128.0,  -4.0/128.0,  2.0/128.0,  -6.0/128.0,  0.0/128.0,  
  0.0/128.0,   0.0/128.0,   0.0/128.0,  0.0/128.0,   0.0/128.0,  0.0/128.0
}  // o_pos
};


//
//  Clip a value between low and high
//
int Clip3Fun(int low, int high, int val)
{
  return ((val < low)? low: (val > high)? high: val);
}


//
//  Determine noSamples
//
int get_thDC(int noMB)
{
  int thDC[4] = NO_SAMPLES;
  int mb[4]   = {   99,  396, 3600, 8100};
  int i;
  int last = 1;

  for(i = 0; i < 4; ++i)
    if(noMB <= mb[i])
      return thDC[i];
    else
      last = thDC[i];
  return last; 
}

//
//  Offset between two frames
//
double getDCdiff(int list, int ref)
{
  int i, j;
  double DCOrg = 0.0;
  double DCRef = 0.0;

  for(i = 0; i < input->img_height; ++i) //y
  {
    for(j = 0; j < input->img_width; ++j) //x
    {
      DCOrg += imgY_org[i][j];
      DCRef += listX[list][ref]->imgY[i][j];
    }
  }

  return ((DCOrg - DCRef) / (double)(input->img_height * input->img_width));
}

//
//  Block-based offset between two frames
//
void blockDC(int list, int *DCmin, int *DCmax, int noSamples)
{
  int i, j, ii, jj, iiLimit, jjLimit;
  int ref_frame, noPixels;
  double DCOrg = 0.0;
  double DCRef = 0.0;
  double DCdiff, temp, err0, err1;
  int DCInt;
  int DCinterval[2][NO_DC_VAL] = {{0}};

  memset(DCinterval, 0, 2 * NO_DC_VAL * sizeof(int));
  ref_frame = 0;

  for(i = 0; i < input->img_height; i+=BL_SIZE)
  {
    for(j = 0; j < input->img_width; j+=BL_SIZE)
    {
      DCOrg = 0.0;
      DCRef = 0.0;
      noPixels = 0;

      for(ii = i; ii < i + BL_SIZE; ii++)
      {
        iiLimit = min(ii, input->img_height - 1);
        for(jj = j; jj < j + BL_SIZE; jj++)
        {
          jjLimit = min(jj, input->img_width - 1);
          DCOrg += imgY_org[iiLimit][jjLimit];
          DCRef += listX[list][ref_frame]->imgY[iiLimit][jjLimit];
          noPixels++;
        }
      }
      DCdiff = (DCOrg - DCRef) / (double)noPixels;
      DCOrg = DCOrg / (double)noPixels;
      DCRef = DCRef / (double)noPixels;

      err0 = 0.0;
      for(ii = i; ii < i + BL_SIZE; ii++)
      {
        iiLimit = min(ii, input->img_height - 1);
        for(jj = j; jj < j + BL_SIZE; jj++)
        {
          jjLimit = min(jj, input->img_width - 1);
          temp = imgY_org[iiLimit][jjLimit] - listX[list][ref_frame]->imgY[iiLimit][jjLimit];
          err0 += temp * temp;
        }
      }

      err1 = 0.0;
      for(ii = i; ii < i + BL_SIZE; ii++)
      {
        iiLimit = min(ii, input->img_height - 1);
        for(jj = j; jj < j + BL_SIZE; jj++)
        {
          jjLimit = min(jj, input->img_width - 1);
          temp = (imgY_org[iiLimit][jjLimit] - DCOrg) - (listX[list][ref_frame]->imgY[iiLimit][jjLimit] - DCRef);
          err1 += temp * temp;
        }
      }
 
      DCInt = min((int)(fabs(DCdiff) + 0.5), NO_DC_VAL - 1);

      if(err1 < (err0 - err1))
      {    
        if(DCdiff < 0)
        {
          DCinterval[0][DCInt]++;
          DCInt = -DCInt;
        }
        else
        {
          DCinterval[1][DCInt]++;
        }
      }
    }
  }

  (*DCmin) = 0;
  (*DCmax) = 0;

#ifdef USE_DCMINMAX
  {
    int flag;
    for(i=0, flag=0; i<NO_DC_VAL; i++)
    {
      if(DCinterval[1][i]>=noSamples)
      {
        (*DCmax)=(*DCmax)+1;
        flag = 1;
      }
      else if(flag)
      {
        i = NO_DC_VAL;
      }
    }
    for(i=0, flag=0; i<NO_DC_VAL; i++)
    {
      if(DCinterval[0][i]>=noSamples)
      {
        (*DCmin)=(*DCmin)+1;
        flag = 1;
      }
      else if(flag)
      {
        i = NO_DC_VAL;
      }
    }
  }
#else
  for(i = 0; i < NO_DC_VAL; i++)
  {
    (*DCmax) = (DCinterval[1][i] >= noSamples)? (*DCmax) + 1: (*DCmax);
    (*DCmin) = (DCinterval[0][i] >= noSamples)? (*DCmin) + 1: (*DCmin);
  }
#endif

  (*DCmin) = -(*DCmin);
}

//
//  Reset subpel offsets for the first pass
//
void resetFirstPassSubpelOffset(int list)
{
  int i;

  img->noOffsets[list] = img->noMEOffsets[list] = 0;
  img->DCmin[list] = img->DCmax[list] = 0;
  img->roundFact[list] = 0;
  img->firstVal[list] = img->secondVal[list] = 0;
  img->noOffsetsSecond[list] = 0;
  img->currMEOffset[list] = 0;

  for(i = 0; i < 16; ++i) 
    img->offsetMETab[list][i] = 0;
}

//
//  Set subpel offsets for the first pass
//
void setFirstPassSubpelOffset(int list)
{
  int subPelPosOrder[15] = {5, 15, 13, 7, 9, 6, 11, 14, 1, 4, 3, 12, 10, 2, 8};
  int i, sign, DCmin, DCmax, thDC;
  int ref_frame, noOffsets, noOffsetsSecond, DCint, firstVal, secondVal, addOffset;
  double DCdiff[5];
                                      //  e   o   m   g   i   f    k    n   a   d   c   l   j   b   h
  static const int search_points_x[15] = {1, -1,  1, -1,  1,  2,  -1,   2,  1,  0, -1,  0,  2,  2,  0};
  static const int search_points_y[15] = {1, -1, -1,  1,  2,  1,   2,  -1,  0,  1,  0, -1,  2,  0,  2};

  memset(DCdiff, 0, 5 * sizeof(double));
  memset(img->offsetMETab[list], 0, 16 * sizeof(int));

  memcpy(img->search_point_qp_x_offset, search_points_x, 15 * sizeof(int));
  memcpy(img->search_point_qp_y_offset, search_points_y, 15 * sizeof(int));

  for(ref_frame = 0; ref_frame < listXsize[list]; ref_frame++)
  {    
    DCdiff[ref_frame] = getDCdiff(list, ref_frame);
    // printf("DCDIFF[%d] = %lf\n", ref_frame,  DCdiff[ref_frame]);

    if(ref_frame==0)
    {
      sign = (DCdiff[ref_frame] >= 0)? 1: -1;
      noOffsets = (int)(fabs(DCdiff[ref_frame]) + 0.5);
      DCint = noOffsets * sign;
 
      if(noOffsets < 2)
      {
        DCmax = noOffsets;
        DCmin = -noOffsets;
      }
      else
      {
        thDC = get_thDC(img->FrameSizeInMbs);
        blockDC(list, &DCmin, &DCmax, thDC); 
        DCmax = max(DCint, DCmax);
        DCmin = min(DCint, DCmin);
      }

      img->roundFact[list] = (int)ceil((double)(DCmax - DCmin) / 15.0);
      img->roundFact[list] = max(img->roundFact[list], 1);      
      img->DCmin[list] = (int)(fabs(DCmin) / img->roundFact[list] + 0.5);
      img->DCmax[list] = (int)(fabs(DCmax) / img->roundFact[list] + 0.5);
      DCmin = -img->roundFact[list] * img->DCmin[list];
      DCmax = img->roundFact[list] * img->DCmax[list];

      noOffsetsSecond = 0;
      if(noOffsets < 2)
        noOffsetsSecond = min((int)(fabs(DCdiff[ref_frame]) / 0.1 + 0.5), 15);

      if(noOffsetsSecond > 0)
      {
        firstVal = 0;
        secondVal = sign;
        if(fabs(DCdiff[ref_frame]) > 0.5)
        {
          firstVal = sign; 
          secondVal = 0;
          noOffsetsSecond = 16 - noOffsetsSecond;
        }
      }
      else
      {
        firstVal = 0; 
        secondVal = 0;
      }
      addOffset = -firstVal;

      if(((firstVal == DCmin) && (secondVal == DCmax)) || ((firstVal == DCmax) && (secondVal == DCmin)))
      {
        DCmin = 0; 
        DCmax = 0;
        img->DCmin[list] = 0; 
        img->DCmax[list] = 0;
      }

      for(i = 0; i < 16; i++)
        img->subpelOffset[list][i] = firstVal;

      img->noOffsets[list] = 0;
      for(i = DCmin; i <= DCmax; i += img->roundFact[list])
      {
        if((i != firstVal) && (i != secondVal))
        {
          img->subpelOffset[list][subPelPosOrder[img->noOffsets[list]]] = i;
          img->offsetMETab[list][img->noOffsets[list]] = i + addOffset;
          (img->noOffsets[list])++;
          if(img->noOffsets[list] == 15)
            break;
        }
      }

      img->noMEOffsets[list] = img->noOffsets[list];
      if((img->noOffsets[list] < 15) && (noOffsetsSecond > 0))
      {
        img->offsetMETab[list][img->noMEOffsets[list]] = secondVal + addOffset;
        (img->noMEOffsets[list])++;
      }

      for(i = 0; i < noOffsetsSecond; i++)
      {
        if(img->noOffsets[list] == 15)
          break;
        img->subpelOffset[list][subPelPosOrder[img->noOffsets[list]]] = secondVal;
        (img->noOffsets[list])++;
      }

      img->firstVal[list] = firstVal;
      img->secondVal[list] = secondVal;
      img->noOffsetsSecond[list] = noOffsetsSecond;
    }
    else
    {
      sign = (DCdiff[ref_frame] >= 0)? 1: -1;

#ifdef THRESHOLD_FRAME_OFFSETS
#define TFO 2.0
      img->imgOffset[list][ref_frame] = (fabs(DCdiff[ref_frame]) < TFO)? 0: 1;
#else
      img->imgOffset[list][ref_frame] = min((int)(fabs(DCdiff[ref_frame]) + 0.5), 1);
#endif  // THRESHOLD_FRAME_OFFSETS

      if(sign < 0)
        img->imgOffset[list][ref_frame] = -img->imgOffset[list][ref_frame];
    }
  }
}

//
//  Find a decision vector different in one position and with smallest SAD 
//
static double get_min_d(double err[NUM_SIFO][NUM_SIFO][16][16], 
                        int old_d[16], int new_d[16], double currSAD)
{
  int spp, flt, sppF, sppB;
  double minSAD = currSAD;
  int tmp_d[16];
  double tmpSAD;

  for(spp = 1; spp < 16; ++spp)
  {
    // Copy a fresh old_d into tmp_d        
    memcpy(tmp_d, old_d, sizeof(int) * 16);

    // Try the remaining (NUM_SIFO - 1) filters on the position spp
    for(flt = 0; flt < NUM_SIFO - 1; ++flt)
    {
      tmp_d[spp] = (tmp_d[spp] + 1) % NUM_SIFO;
      tmpSAD = 0;
      for(sppF = 1; sppF < 16; ++sppF)
        for(sppB = 1; sppB < 16; ++sppB)
          tmpSAD += err[tmp_d[sppF]][tmp_d[sppB]][sppF][sppB];

      if(tmpSAD < minSAD)
      {
        minSAD = tmpSAD;
        memcpy(new_d, tmp_d, sizeof(int) * 16);
      }
    }
  }
  return minSAD;
}


//
//  Compute filter combination for frame or sequence
//
static void ComputeFilterCombination_B_gd(double err[NUM_SIFO][NUM_SIFO][16][16], int out[16])
{
  int min_d[16] = {0};
  double oldSAD = 0.0;
  double minSAD;
  int sppF, sppB;
  int tmp_d[16] = {0};

  // Compute error for min_d = 0^n
  for(sppF = 1; sppF < 16; ++sppF)
    for(sppB = 1; sppB < 16; ++sppB)
      oldSAD += err[0][0][sppF][sppB];

  // While SAD can be improved
  while((minSAD = get_min_d(err, min_d, tmp_d, oldSAD)) < oldSAD)
  {
    oldSAD = minSAD;        
    memcpy(min_d, tmp_d, sizeof(int) * 16);
  }

  // Copy result
  memcpy(out, min_d, sizeof(int) * 16);
}


//
//  Check whether a unique filter has a RD advantage over the filter combination
//
static void ComputeBestFilter_B(double err[NUM_SIFO][NUM_SIFO][16][16], int combination[16])
{
  int bits = 0;
  int i;
  double errComb = 0.0;
  double errBest = 0.0;
  int sppF, sppB;
  int best = -1;

  // Lowest error with a single filter
  for(i = 0; i < NUM_SIFO; ++i)
  {
    double errTmp = 0.0;
    for(sppF = 1; sppF < 16; ++sppF)
      for(sppB = 1; sppB < 16; ++sppB)
        errTmp += err[i][i][sppF][sppB];
    if((best == -1) || (errTmp < errBest))
    {
      best = i;
      errBest = errTmp;
    }
  }

  // Compute bit cost of sending the filters
  for(i = 1; i < 16; ++i)
    bits += (combination[i] == 0)? 1: 2;

  // Error with filter combination
  for(sppF = 1; sppF < 16; ++sppF)
    for(sppB = 1; sppB < 16; ++sppB)
      errComb += err[combination[sppF]][combination[sppB]][sppF][sppB];


  if(errBest < (errComb + img->lambda_md[img->type][img->qp] * (double)bits))
  {
    for(i = 0; i < 16; ++i){
      combination[i] = best;
    }
    img->bestFilter = best + 1;
  }
  else
  {
    img->bestFilter = 0;
  }
}

//
//  Compute filter combination for frame and sequence
//
static void ComputeFilterCombination_B()
{
  ComputeFilterCombination_B_gd(AccErrorB, BestCombFilterB);
  ComputeFilterCombination_B_gd(SequenceAccErrorB, SequenceBestCombFilterB);
  // Compute best single sequence filter
  ComputeBestFilter_B(SequenceAccErrorB, SequenceBestCombFilterB);
}

//
//
//
static void ResetAll(int slice)
{
  int a, b, c, d;
  static int firstP = 1;
  static int firstB = 1;

  if(slice == P_SLICE)
  {
    memset(SamplesSubPelOffsetP, 0, 16 * sizeof(int));
    memset(SamplesFrameOffsetP, 0, MAX_REFERENCE_PICTURES * sizeof(int));

    resetFirstPassSubpelOffset(LIST_0);

    for(a = 0; a < NUM_SIFO; ++a)
      memset(AccErrorP[a], 0, 16 * sizeof(double));

    for(a = 0; a < MAX_REFERENCE_PICTURES; ++a)
      for(b = 0; b < 16; ++b)
        for(c = 0; c < NUM_SIFO; ++c)
          AccFrameOffsetP[a][b][c] = 0;

    if(firstP == 1)
    {
      firstP = 0;

      for(a = 0; a < NUM_SIFO; ++a)
        memset(SequenceAccErrorP[a], 0, 16 * sizeof(double));
    }
    else
    {
      for(a = 0; a < 16; ++a)
        for(b = 0; b < NUM_SIFO; ++b)
          SequenceAccErrorP[b][a] *= NORM_FACTOR;
    }
  }
  else if(slice == B_SLICE)
  {
    // Error
    for(a = 0; a < NUM_SIFO; ++a)
      for(b = 0; b < NUM_SIFO; ++b)
        for(c = 0; c < 16; ++c)
          for(d = 0; d < 16; ++d)
          {
            AccErrorB[a][b][c][d] = 0.0;
            SamplesB[a][b][c][d] = 0;
          }

    // Frame and SubPel offsets
    for(a = 0; a < 2; ++a)
    {
      memset(AccFrameOffset[a], 0, MAX_REFERENCE_PICTURES * sizeof(double));
      memset(FrameOffset[a], 0, MAX_REFERENCE_PICTURES * sizeof(int));
      memset(SamplesFrameOffset[a], 0, MAX_REFERENCE_PICTURES * sizeof(int));
      memset(AccSubPelOffset[a], 0, 16 * sizeof(double));
      memset(SubPelOffset[a], 0, 16 * sizeof(int));
      memset(SamplesSubPelOffset[a], 0, 16 * sizeof(int));
    }

    memset(BestCombFilterB, 0, 16 * sizeof(int));

    resetFirstPassSubpelOffset(LIST_0);
    resetFirstPassSubpelOffset(LIST_1);

    if(firstB == 1)
    {    
      firstB = 0;
      memset(SequenceBestCombFilterB, 0, 16 * sizeof(int));

      for(a = 0; a < NUM_SIFO; ++a)
        for(b = 0; b < NUM_SIFO; ++b)
          for(c = 0; c < 16; ++c)
            for(d = 0; d < 16; ++d)
              SequenceAccErrorB[a][b][c][d] = 0;
    }
    else
    {
      for(a = 0; a < NUM_SIFO; ++a)
        for(b = 0; b < NUM_SIFO; ++b)
          for(c = 0; c < 16; ++c)
            for(d = 0; d < 16; ++d)
              SequenceAccErrorB[a][b][c][d] *= NORM_FACTOR;
    }
  }  
}

//
//  Return the value of a pixel interpolated with filter filterNo 
//
static int GetInterpolatedPixel(unsigned short ** img, int icoord, int jcoord, 
                         int img_width, int img_height, int sub_pos, int filterNo)
{
  int pos_x, pos_y;
  double ipVal;
  int out;
  int ii, jj;
  int i, j;
#ifndef EDAIF2
#define FILTER_SIZE       6                     // Undefined at the end of the function
#define FILTER_OFFSET     (6 / 2 - 1)           // Undefined at the end of the function
#endif

  assert(filterNo < NUM_SIFO);

  ipVal = 0.0;
  i = icoord - 4 * IMG_PAD_SIZE;
  j = jcoord - 4 * IMG_PAD_SIZE;

  if(sub_pos)
  {
#ifdef HPF_COMBO
    if((sub_pos == 7)  && (filterNo == 0))		// g
    {	
      // Strong filter position
      static const int f0[4][4] = {{  0,  5,  5,  0},
                                   {  5, 22, 22,  5},
                                   {  5, 22, 22,  5},
                                   {  0,  5,  5,  0}};
      int sfp = 0;
      int l, k;
      for(l = 0; l < 4; ++l)
      {
        for(k = 0; k < 4; ++k)
        {
          pos_y = Clip3Fun(0, img_height - 1, i / 4 - 1 + l);
          if((i < 0) && (pos_y > 0))
            pos_y -= 1;
          pos_x = Clip3Fun(0, img_width - 1, j / 4 - 1 + k);
          if((j < 0) && (pos_x > 0))
            pos_x -=1;
          sfp += f0[l][k] * img[pos_y][pos_x];
        }
      }
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
      out = Clip3Fun(0, (1 << input->BitDepthLuma) - 1, (sfp + 64) >> 7);
#else
      out = Clip3Fun(0, 255, (sfp + 64) >> 7);
#endif
    }
    else
    {
#endif  // HPF_COMBO
      for(ii = 0; ii < FILTER_SIZE; ++ii)
      {
        for(jj = 0; jj < FILTER_SIZE; ++jj)
        {
          pos_y = Clip3Fun(0, img_height - 1, i / 4 - FILTER_OFFSET + ii);
          if((i < 0) && (pos_y > 0))
            pos_y -= 1;
          pos_x = Clip3Fun(0, img_width - 1, j / 4 - FILTER_OFFSET + jj);
          if((j < 0) && (pos_x > 0))
            pos_x -=1;
          ipVal += (SIFO_FILTER[filterNo][sub_pos - 1][FILTER_SIZE * ii + jj] * img[pos_y][pos_x]);
        }
      }
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
      out = Clip3Fun(0, (1 << input->BitDepthLuma) - 1, (int)(ipVal + 0.5));
#else
      out = Clip3Fun(0, 255, (int)(ipVal + 0.5));
#endif
#ifdef HPF_COMBO
    }
#endif  // HPF_COMBO
  }
  else
  {
    out = img[Clip3Fun(0, img_height - 1, i / 4)][Clip3Fun(0, img_width - 1, j / 4)];
  }

#undef FILTER_OFFSET
#undef FILTER_SIZE

  return out;
}

//
//
//
static void AccumulateError_P(int x_orig, int y_orig)
{
  int x, y, xj, yi, temp, valOrg;
  int mvx = 0,  mvy = 0;
  int sub_pos = 0, mvx_sub, mvy_sub;
  int ref_frame = 0;
  int subblock = 0;
  int out4Y_width  = (img->width + 2 * IMG_PAD_SIZE) * 4 - 1;  
  int out4Y_height = (img->height + 2 * IMG_PAD_SIZE) * 4 - 1;
  unsigned short ** imgY;
  int f;

  for(subblock = 0; subblock < 16; ++subblock)
  {
    x = x_orig + 4 * (subblock % 4);
    y = y_orig + 4 * (subblock / 4);
    mvx = enc_picture->mv[LIST_0][y / 4][x / 4][0];
    mvy = enc_picture->mv[LIST_0][y / 4][x / 4][1];

    mvx_sub = (mvx >= 0)? mvx % 4: (4 - abs(mvx) % 4) % 4;
    mvy_sub = (mvy >= 0)? mvy % 4: (4 - abs(mvy) % 4) % 4;
    sub_pos = mvx_sub + 4 * mvy_sub;    // pos 0..15 in a 4x4 block 

    ref_frame = enc_picture->ref_idx[LIST_0][y / 4][x / 4];

    if(ref_frame != -1)
    {
      imgY = listX[LIST_0][ref_frame]->imgY;

      for(yi = 0; yi < 4; ++yi)
      {
        for(xj = 0; xj < 4; ++xj)
        {
          int xcoord = Clip3Fun(0, out4Y_width, 4 * (x + xj) + 4 * IMG_PAD_SIZE + mvx);
          int ycoord = Clip3Fun(0, out4Y_height, 4 * (y + yi) + 4 * IMG_PAD_SIZE + mvy);

          valOrg = imgY_org[y+yi][x+xj];

          for(f = 0; f < NUM_SIFO; ++f)
          {
            temp = valOrg - GetInterpolatedPixel(imgY, ycoord, xcoord, img->width, img->height, sub_pos, f);
            AccFrameOffsetP[ref_frame][sub_pos][f] += temp;
            AccErrorP[f][sub_pos] += temp * temp;
          }

          SamplesFrameOffsetP[ref_frame]++;

          if(ref_frame == 0)
            SamplesSubPelOffsetP[sub_pos]++;
        }
      }
    }
  }
}


//
//  Update table of offsets used in ME
//
void updateMEOffsets(int list)
{
                                // 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
  const int search_points[2][15]={{1,  2, -1,  0,  1,  2, -1,  0,  1,  2, -1,  0,  1,  2, -1},
                                  {0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2, -1, -1, -1, -1}};
  int i, j;

  memset(img->search_point_qp_x_offset, 0, 15 * sizeof(int));
  memset(img->search_point_qp_y_offset, 0, 15 * sizeof(int));
  memset(img->offsetMETab, 0, 16 * sizeof(int));

  img->noOffsets[list] = 0; 
  for(i = 1; i < 16; ++i)
  {
    if(img->subpelOffset[list][i] != 0)
    {
      img->search_point_qp_x_offset[img->noOffsets[list]] = search_points[0][i-1];
      img->search_point_qp_y_offset[img->noOffsets[list]] = search_points[1][i-1];
      (img->noOffsets[list])++;   
    }
  }

  img->noMEOffsets[list] = 0; 
  for(i = 1; i < 16; ++i)
  {
    if(img->subpelOffset[list][i] != 0)
    {
      if(img->noMEOffsets[list] == 0)
      {
        img->offsetMETab[list][img->noMEOffsets[list]] = img->subpelOffset[list][i] - img->subpelOffset[list][0];
        (img->noMEOffsets[list])++;
      }
      else
      {
        for(j = 0; j < img->noMEOffsets[list]; ++j)
        {
          if(img->offsetMETab[list][j] != img->subpelOffset[list][i])
          {
            img->offsetMETab[list][img->noMEOffsets[list]] = img->subpelOffset[list][i] - img->subpelOffset[list][0];
            (img->noMEOffsets[list])++;
            break;
          } 
        }
      }
    }
  }
}

//
//
//
static void ComputeFiltersAndOffsets_P()
{
  int i, frame, ind, offset;
  double combAccSamplesFrameOffsetP[MAX_REFERENCE_PICTURES], combAccFrameOffsetP[16], minVal;
    
  // Filter selection based on the current frame
  for(i = 0; i < 16; ++i)
  {
    img->filterFrame[i] = 0;
    minVal = AccErrorP[0][i];
    for(ind = 1; ind < NUM_SIFO; ++ind)
    {
      if(AccErrorP[ind][i] < minVal)
      {
        img->filterFrame[i] = ind;
        minVal = AccErrorP[ind][i];
      }
    }
  }

  // Frame offsets corresponding to the frame filters
  for(frame = 0; frame < listXsize[LIST_0]; ++frame)
  {
    combAccFrameOffsetP[frame] = 0;
    for(i = 0; i < 16; ++i)
      combAccFrameOffsetP[frame] += AccFrameOffsetP[frame][i][img->filterFrame[i]];

    img->imgOffset[LIST_0][frame] = 0;
    if(SamplesFrameOffsetP[frame] > 0)
    {
      offset = (int)(fabs(combAccFrameOffsetP[frame]) / (double)SamplesFrameOffsetP[frame] + 0.5);
      img->imgOffset[LIST_0][frame] = (combAccFrameOffsetP[frame] >= 0)? offset: -offset;
    }
  }

  // Sub-pel offsets
  for(i = 0; i < 16; ++i)
  {
    combAccSamplesFrameOffsetP[i] = AccFrameOffsetP[0][i][img->filterFrame[i]];

    img->subpelOffset[LIST_0][i] = 0;
    if(SamplesSubPelOffsetP[i] > 0)
    {
      offset = (int)((double)fabs(combAccSamplesFrameOffsetP[i]) / (double)SamplesSubPelOffsetP[i] + 0.5);
      img->subpelOffset[LIST_0][i] = (combAccSamplesFrameOffsetP[i] >= 0)? offset: -offset;
    }
  }
  
  // Update ME offsets
  // updateMEOffsets(LIST_0);
}

//
//  Check whether a unique filter has a RD advantage over the filter combination
//
static void ComputeBestFilter_P(double err[NUM_SIFO][16], int combination[16])
{
  int bits = 0;
  int i;
  double errComb = 0.0;
  double errBest = 0.0;
  int spp;
  int best = -1;

  // Lowest error with a single filter
  for(i = 0; i < NUM_SIFO; ++i)
  {
    double errTmp = 0.0;
    for(spp = 1; spp < 16; ++spp)
      errTmp += err[i][spp];
    if((best == -1) || (errTmp < errBest))
    {
      best = i;
      errBest = errTmp;
    }
  }

  // Compute bit cost of sending the filters
  for(i = 1; i < 16; ++i)
    bits += (combination[i] == 0)? 1: 2;

  // Error with filter combination
  for(spp = 1; spp < 16; ++spp)
      errComb += err[combination[spp]][spp];


  if(errBest < (errComb + img->lambda_md[img->type][img->qp] * (double)bits))
  {
    for(i = 0; i < 16; ++i){
      combination[i] = best;
    }
    img->bestFilter = best + 1;
  }
  else
  {
    img->bestFilter = 0;
  }
}


//
//
//
void UpdateSequenceFilters_P(void)
{
  int i, f, ind;
  double minVal;

  // Filter selection based on the sequence 
  for(i = 0; i < 16; ++i)
  {
    for(f = 0; f < NUM_SIFO; ++f)
      SequenceAccErrorP[f][i] += AccErrorP[f][i];

    img->filterSequence[i] = 0;
    minVal = SequenceAccErrorP[0][i];
    for(ind = 1; ind < NUM_SIFO; ++ind)
    {
      if(SequenceAccErrorP[ind][i] < minVal)
      {
        img->filterSequence[i] = ind;
        minVal = SequenceAccErrorP[ind][i];
      }
    }
  }

  ComputeBestFilter_P(SequenceAccErrorP, img->filterSequence);
}


//
//  Accumulate squared error to compute best filter combination
//
static void AccumulateError_B(int x_orig, int y_orig)
{
  int x, y, xj, yi, valOrg, i, j;
  int mvx[2],  mvy[2];
  int sub_pos[2], mvx_sub, mvy_sub;
  int ref_frame[2];
  int subblock=0;
  int out4Y_width  = (img->width + 2 * IMG_PAD_SIZE) * 4 - 1;  
  int out4Y_height = (img->height + 2 * IMG_PAD_SIZE) * 4 - 1;
  double tmp;
  int apply_weights = ((img->filterParam == SIFO_SEQ_FILTER) && (active_pps->pic_parameter_set_id == 2));
  unsigned short ** imgY0;
  unsigned short ** imgY1;
  unsigned short ** imgY;
  int f;

  for(subblock = 0; subblock < 16; ++subblock)
  {
    x = x_orig + 4 * (subblock % 4);
    y = y_orig + 4 * (subblock / 4);

    //get motion information
    for (i = 0; i < 2; ++i)
    {
      mvx[i] = enc_picture->mv[LIST_0+i][y/4][x/4][0];
      mvy[i] = enc_picture->mv[LIST_0+i][y/4][x/4][1];

      mvx_sub = (mvx[i] >= 0)? mvx[i] % 4: (4 - abs(mvx[i]) % 4) % 4;
      mvy_sub = (mvy[i] >= 0)? mvy[i] % 4: (4 - abs(mvy[i]) % 4) % 4;

      sub_pos[i] = mvx_sub + 4 * mvy_sub;    // Pos 0..15 in a 4x4 block
      ref_frame[i] = enc_picture->ref_idx[LIST_0+i][y/4][x/4];
    }

    if ((ref_frame[0] == -1) && (ref_frame[1] == -1))
      continue;     // Intra macroblocks are not used for calculation

    if ((ref_frame[0] != -1) && (ref_frame[1] != -1)) // Bi-predicted
    {
      int tempF[NUM_SIFO], tempB[NUM_SIFO];
      int fw_ref_idx = ref_frame[0];
      int bw_ref_idx = ref_frame[1];
      
      imgY0 = listX[LIST_0][fw_ref_idx]->imgY;
      imgY1 = listX[LIST_1][bw_ref_idx]->imgY;

      for(yi = 0; yi < 4; ++yi)
      {
        for(xj = 0; xj < 4; ++xj)
        {
          int xcoord0 = Clip3Fun(0 , out4Y_width,  4 * (x + xj) + 4 * IMG_PAD_SIZE + mvx[0]);
          int ycoord0 = Clip3Fun(0 , out4Y_height, 4 * (y + yi) + 4 * IMG_PAD_SIZE + mvy[0]);
          int xcoord1 = Clip3Fun(0 , out4Y_width,  4 * (x + xj) + 4 * IMG_PAD_SIZE + mvx[1]);
          int ycoord1 = Clip3Fun(0 , out4Y_height, 4 * (y + yi) + 4 * IMG_PAD_SIZE + mvy[1]);

          valOrg   = imgY_org[y + yi][x + xj];

          for(f = 0; f < NUM_SIFO; ++f)
          {
            tempF[f] = GetInterpolatedPixel(imgY0, ycoord0, xcoord0, img->width, img->height, sub_pos[0], f);          
            tempB[f] = GetInterpolatedPixel(imgY1, ycoord1, xcoord1, img->width, img->height, sub_pos[1], f);
          }

          for(i = 0; i < NUM_SIFO; ++i)
          {
            for(j = 0; j < NUM_SIFO; ++j)
            {
              if(apply_weights)
              {            
                tmp = valOrg - clip1a(((wbp_weight[0][fw_ref_idx][bw_ref_idx][0] * tempF[i] + 
                  wbp_weight[1][fw_ref_idx][bw_ref_idx][0] * tempB[j] + 
                  (wp_offset[0][fw_ref_idx][0] + wp_offset[1][bw_ref_idx][0]) * 2 * wp_luma_round +
                  2 * wp_luma_round * (1 - img->bipred_rounding_control)) >> (luma_log_weight_denom + 1))); 
              }
              else
              {
                tmp = valOrg - ((tempF[i] + tempB[j] + 1 - img->bipred_rounding_control) >> 1); 
              }

              AccErrorB[i][j][sub_pos[0]][sub_pos[1]] += (tmp * tmp);
              SamplesB[i][j][sub_pos[0]][sub_pos[1]]++;
              SequenceAccErrorB[i][j][sub_pos[0]][sub_pos[1]] += (tmp * tmp);
            }
          }
        } // x
      } // y
    }
    else //uni-directional
    {
      int ref_index = (ref_frame[0]!= -1)? LIST_0: LIST_1;
      int ind = (ref_frame[0]!= -1)? 0: 1;

      imgY = listX[ref_index][ref_frame[ind]]->imgY;

      for(yi = 0; yi < 4; ++yi)
      {    
        for(xj = 0; xj < 4; ++xj)
        {
          int xcoordind = Clip3Fun(0 , out4Y_width,  4 * (x + xj) + 4 * IMG_PAD_SIZE + mvx[ind]);
          int ycoordind = Clip3Fun(0 , out4Y_height, 4 * (y + yi) + 4 * IMG_PAD_SIZE + mvy[ind]);

          valOrg = imgY_org[y+yi][x+xj];

          for(i = 0; i < NUM_SIFO; ++i)
          {
            // Unidirectional error is accumulated along the diagonal
            tmp = valOrg - GetInterpolatedPixel(imgY, ycoordind, xcoordind, img->width, img->height, sub_pos[ind], i);
            AccErrorB[i][i][sub_pos[ind]][sub_pos[ind]] += (tmp * tmp);
            SamplesB[i][i][sub_pos[ind]][sub_pos[ind]]++;
            SequenceAccErrorB[i][i][sub_pos[ind]][sub_pos[ind]] += (tmp * tmp);
          }
        } // x
      } // y
    } // if bidirectional
  } // for(subblock = 0; subblock < 16; subblock++)
}

//
//  Accumulate values to compute offset for F reference frames
//
static void AccumulateOffset0_B(int x_orig, int y_orig)
{
  int x, y, xj, yi, valOrg, i;
  int mvx[2], mvy[2];
  int sub_pos[2], mvx_sub, mvy_sub;
  int ref_frame[2];
  int subblock = 0;
  int out4Y_width  = (img->width + 2 * IMG_PAD_SIZE) * 4 - 1;  
  int out4Y_height = (img->height + 2 * IMG_PAD_SIZE) * 4 - 1;
  double tmp;
  int filterF;
  unsigned short ** imgY;

  if(img->type == B_SLICE)
  {
    for(subblock = 0; subblock < 16; ++subblock)
    {
      x = x_orig + 4 * (subblock % 4);
      y = y_orig + 4 * (subblock / 4);

      //get motion information
      for (i = 0; i < 2; ++i)
      {
        mvx[i] = enc_picture->mv[LIST_0+i][y/4][x/4][0];
        mvy[i] = enc_picture->mv[LIST_0+i][y/4][x/4][1];

        mvx_sub = (mvx[i] >= 0)? mvx[i] % 4: (4 - abs(mvx[i]) % 4) % 4;
        mvy_sub = (mvy[i] >= 0)? mvy[i] % 4: (4 - abs(mvy[i]) % 4) % 4;

        sub_pos[i] = mvx_sub + 4 * mvy_sub;    // pos 0..15 in a 4x4 block
        ref_frame[i] = enc_picture->ref_idx[LIST_0+i][y/4][x/4];
      }

      if(ref_frame[0] != -1)
      {
        imgY = listX[LIST_0][ref_frame[0]]->imgY;

        for(yi = 0; yi < 4; ++yi)
        {
          for(xj = 0; xj < 4;++xj)
          {
            int xcoord0 = Clip3Fun(0 , out4Y_width,  4 * (x + xj) + 4 * IMG_PAD_SIZE + mvx[0]);
            int ycoord0 = Clip3Fun(0 , out4Y_height, 4 * (y + yi) + 4 * IMG_PAD_SIZE + mvy[0]);

            valOrg = imgY_org[y+yi][x+xj];
            filterF = BestCombFilterB[sub_pos[0]];      // Use best filter

            tmp = valOrg - GetInterpolatedPixel(imgY, ycoord0, xcoord0, img->width, img->height, sub_pos[0], filterF);

            AccFrameOffset[LIST_0][ref_frame[0]] += tmp;
            SamplesFrameOffset[LIST_0][ref_frame[0]]++;
            if(ref_frame[0] == 0)
            {
              AccSubPelOffset[LIST_0][sub_pos[0]] += tmp;
              SamplesSubPelOffset[LIST_0][sub_pos[0]]++;
            }
          } // x
        } // y
      }
    } // for(subblock = 0; subblock < 16; subblock++)
  } // if (img->type==B_SLICE)
}

//
//  Computes the offsets for the F reference frames
//
static void ComputeOffset0(void)
{
  int i;
  int offset;

  // Frame offsets
  for (i = 0; i < listXsize[LIST_0]; ++i)
  {
    if(SamplesFrameOffset[LIST_0][i] > 0)
    {
      offset = (int)((double)fabs(AccFrameOffset[LIST_0][i]) / (double)SamplesFrameOffset[LIST_0][i] + 0.5);
      if (AccFrameOffset[LIST_0][i] < 0)
        offset = -offset;
    }
    else
    {
      offset = 0;
    }
    FrameOffset[LIST_0][i] = offset; // Frame offset
  }

  // SubPel offsets
  for(i = 0; i < 16; ++i)
  {
    if(SamplesSubPelOffset[LIST_0][i] > 0)
    {
      offset = (int)((double)fabs(AccSubPelOffset[LIST_0][i]) / (double)SamplesSubPelOffset[LIST_0][i] + 0.5);
      if (AccSubPelOffset[LIST_0][i] < 0)
        offset = -offset;
    }
    else
    {
      offset = 0;
    }
    SubPelOffset[LIST_0][i] = offset;
  }
}


//
//  Accumulate values to compute offset for B reference frames
//
static void AccumulateOffset1_B(int x_orig, int y_orig)
{
  int x, y, xj, yi, valOrg,i;
  int mvx[2],  mvy[2];
  int sub_pos[2], mvx_sub, mvy_sub;
  int ref_frame[2];
  int subblock = 0;
  int out4Y_width  = (img->width + 2 * IMG_PAD_SIZE) * 4 - 1;  
  int out4Y_height = (img->height + 2 * IMG_PAD_SIZE) * 4 - 1;
  double tmp, tmp1, tmp2;
  double w0, w1, denom, round;
  unsigned short ** imgY0;
  unsigned short ** imgY1;
  unsigned short ** imgY;
  int tempF, tempB, tempFsp;
  int filterF, filterB;

  int apply_weights = ((img->filterParam == SIFO_SEQ_FILTER) && (active_pps->pic_parameter_set_id == 2));

  for(subblock = 0; subblock < 16; ++subblock)
  {
    x = x_orig+4*(subblock%4);
    y = y_orig+4*(subblock/4);

    //get motion information
    for (i = 0; i < 2; ++i)
    {
      mvx[i] = enc_picture->mv[LIST_0+i][y/4][x/4][0];
      mvy[i] = enc_picture->mv[LIST_0+i][y/4][x/4][1];

      mvx_sub = (mvx[i] >= 0)? mvx[i] % 4: (4 - abs(mvx[i]) % 4) % 4;
      mvy_sub = (mvy[i] >= 0)? mvy[i] % 4: (4 - abs(mvy[i]) % 4) % 4;

      sub_pos[i] = mvx_sub + 4 * mvy_sub;    // pos 0..15 in a 4x4 block
      ref_frame[i] = enc_picture->ref_idx[LIST_0+i][y/4][x/4];
    }

    if ((ref_frame[0] != -1) && (ref_frame[1] != -1)) // bi-directional
    {
      int fw_ref_idx = ref_frame[0];
      int bw_ref_idx = ref_frame[1];

      imgY0 = listX[LIST_0][fw_ref_idx]->imgY;
      imgY1 = listX[LIST_1][bw_ref_idx]->imgY;

      if(apply_weights)
      {
        w0 = wbp_weight[0][fw_ref_idx][bw_ref_idx][0];
        w1 = wbp_weight[1][fw_ref_idx][bw_ref_idx][0];
        denom = (double)(1 << (luma_log_weight_denom + 1));
        round = 2.0 * wp_luma_round * (1 - img->bipred_rounding_control);
      }
      else
      {
        w0 = 1.0;
        w1 = 1.0;
        denom = 2.0;
        round = 1.0 - img->bipred_rounding_control;
      }

      for(yi = 0; yi < 4; ++yi)
      {
        for(xj = 0; xj < 4; ++xj)
        {
          int xcoord0 = Clip3Fun(0 , out4Y_width,  4 * (x + xj) + 4 * IMG_PAD_SIZE + mvx[0]);
          int ycoord0 = Clip3Fun(0 , out4Y_height, 4 * (y + yi) + 4 * IMG_PAD_SIZE + mvy[0]);
          int xcoord1 = Clip3Fun(0 , out4Y_width,  4 * (x + xj) + 4 * IMG_PAD_SIZE + mvx[1]);
          int ycoord1 = Clip3Fun(0 , out4Y_height, 4 * (y + yi) + 4 * IMG_PAD_SIZE + mvy[1]);

          valOrg   = imgY_org[y+yi][x+xj];

          filterF = BestCombFilterB[sub_pos[0]];
          filterB = BestCombFilterB[sub_pos[1]];

          tempF = GetInterpolatedPixel(imgY0, ycoord0, xcoord0, img->width, img->height, sub_pos[0], filterF) 
            + FrameOffset[LIST_0][fw_ref_idx];
          tempB = GetInterpolatedPixel(imgY1, ycoord1, xcoord1, img->width, img->height, sub_pos[1], filterB);
       
          tmp1 = valOrg - clip1a((int)(((w0 * tempF) + (w1 * tempB) + round) / denom)); 
          AccFrameOffset[LIST_1][bw_ref_idx] += (denom / w1) * tmp1;
          SamplesFrameOffset[LIST_1][bw_ref_idx]++;

          // Sub pel calculation is for B, so it is for sub_pos[1]
          if(ref_frame[1] == 0)
          {
            tempFsp = GetInterpolatedPixel(imgY0, ycoord0, xcoord0, img->width, img->height, sub_pos[0], filterF) 
              + SubPelOffset[LIST_0][sub_pos[0]];
            tmp2 = valOrg - clip1a((int)(((w0 * tempFsp) + (w1 * tempB) + round) / denom)); 
            AccSubPelOffset[LIST_1][sub_pos[1]] += (denom / w1) * tmp2;
            SamplesSubPelOffset[LIST_1][sub_pos[1]]++;
          }
        } // x
      } // y
    }
    else if(ref_frame[1] != -1)   // Unidirectional, consider only B
    {
      imgY = listX[LIST_1][ref_frame[1]]->imgY;

      for(yi = 0; yi < 4; ++yi)
      {
        for(xj = 0; xj < 4; ++xj)
        {
          int xcoord1 = Clip3Fun(0 , out4Y_width,  4 * (x + xj) + 4 * IMG_PAD_SIZE + mvx[1]);
          int ycoord1 = Clip3Fun(0 , out4Y_height, 4 * (y + yi) + 4 * IMG_PAD_SIZE + mvy[1]);

          valOrg = imgY_org[y+yi][x+xj];

          filterB = BestCombFilterB[sub_pos[1]];
          tempB = GetInterpolatedPixel(imgY, ycoord1, xcoord1, img->width, img->height, sub_pos[1], filterB);
          tmp = valOrg - tempB;    // Use best filter

          AccFrameOffset[LIST_1][ref_frame[1]] += tmp;
          SamplesFrameOffset[LIST_1][ref_frame[1]]++;
          if(ref_frame[1] == 0)
          {
            AccSubPelOffset[LIST_1][sub_pos[1]] += tmp;
            SamplesSubPelOffset[LIST_1][sub_pos[1]]++;
          }
        } // x
      } // y
    }
  }
}

//
//  Compute the offsets for the B reference frames
//
static void ComputeOffset1(void)
{
  int i;
  int offset;

  // Frame offsets
  for(i = 0; i < listXsize[LIST_1]; ++i)
  {
    if(SamplesFrameOffset[LIST_1][i] > 0)
    {
      offset = (int)((double)fabs(AccFrameOffset[LIST_1][i]) / (double)SamplesFrameOffset[LIST_1][i] + 0.5);
      if (AccFrameOffset[LIST_1][i] < 0)
        offset = -offset;
    }
    else
    {
      offset = 0;
    }
    FrameOffset[LIST_1][i] = offset; // Frame offset
  }

  // SubPel offsets
  for(i = 0; i < 16; ++i)
  {
    if(SamplesSubPelOffset[LIST_1][i] > 0)
    {
      offset = (int)((double)fabs(AccSubPelOffset[LIST_1][i]) / (double)SamplesSubPelOffset[LIST_1][i] + 0.5);
      if(AccSubPelOffset[LIST_1][i] < 0)
        offset = -offset;
    }
    else
    {
      offset = 0;
    }
    SubPelOffset[LIST_1][i] = offset;
  }
}


//
//
//
void SwapFilteredFrames(void)
{
  unsigned frame;

  for(frame = 0; frame < dpb.ref_frames_in_buffer; ++frame)
  {
    dpb.fs_ref[frame]->frame->p_imgY_filt = dpb.fs_ref[frame]->frame->imgY_filt;
    dpb.fs_ref[frame]->frame->imgY_filt = dpb.fs_ref[frame]->frame->imgY;
    dpb.fs_ref[frame]->frame->imgY = dpb.fs_ref[frame]->frame->p_imgY_filt;

    dpb.fs_ref[frame]->frame->p_imgY_11_filt = dpb.fs_ref[frame]->frame->imgY_11_filt;
    dpb.fs_ref[frame]->frame->imgY_11_filt = dpb.fs_ref[frame]->frame->imgY_11;
    dpb.fs_ref[frame]->frame->imgY_11 = dpb.fs_ref[frame]->frame->p_imgY_11_filt;
  }
}


//
//
//
static void ScanMacroblocksAndApplyFunction(void fun(int, int))
{
  int i, j;

  if((img->type == P_SLICE) || (img->type == B_SLICE))
  {
    // Loop over all macroblocks and call fun()
    for (i = 0; i < img->height; i = i + MB_BLOCK_SIZE)
      for (j = 0; j < img->width; j = j + MB_BLOCK_SIZE)
        fun(j, i);
  }
}

//
//  Update the filter choice for the sequence
//
void UpdateSequenceFilters_B(void)
{
  int i;
  for(i = 0; i < 16; ++i)
    img->filterSequence[i] = SequenceBestCombFilterB[i];
}

// 
//  Main entry point
//
int ComputeFiltersAndOffsets(void)
{
  int i;

  if(img->type == P_SLICE)
  {
    ResetAll(img->type);
    ScanMacroblocksAndApplyFunction(AccumulateError_P);       // Accumulate squared error for filter decision and offsets
    ComputeFiltersAndOffsets_P();
  }
  else if(img->type == B_SLICE)
  {
    ResetAll(img->type);

    ScanMacroblocksAndApplyFunction(AccumulateError_B);       // Accumulate squared error for frame and sequence
    ComputeFilterCombination_B();                             // Compute filter combination (unique filter per encoded frame)

    ScanMacroblocksAndApplyFunction(AccumulateOffset0_B);     // Accumulate offset data for F reference frame
    ComputeOffset0();                                         // Compute offsets for F (one offset per REFERENCE frame)

    ScanMacroblocksAndApplyFunction(AccumulateOffset1_B);     // Accumulate offset data for B reference frame
    ComputeOffset1();                                         // Compute offsets for B (one offset per REFERENCE frame)

    // Copy filter settings
    for(i = 0; i < 16; ++i)
      img->filterFrame[i] = BestCombFilterB[i];

    // Copy offsets
    for(i = 0; i < 16; ++i)
    {
      img->subpelOffset[LIST_0][i] = SubPelOffset[LIST_0][i];
      img->subpelOffset[LIST_1][i] = SubPelOffset[LIST_1][i];
    }

    // Update ME offsets
    // updateMEOffsets(LIST_0);
    // updateMEOffsets(LIST_1);

    for(i = 0; i < listXsize[LIST_0]; ++i)
      img->imgOffset[LIST_0][i] = FrameOffset[LIST_0][i];

    for(i = 0; i < listXsize[LIST_1]; ++i)
      img->imgOffset[LIST_1][i] = FrameOffset[LIST_1][i];
  }

  return 1; 
}

//
//
//
void initSIFOFilters()
{
  static int first = 1;

  if(first == 1)
  {  
    int i, j;
    first = 0;

    memset(img->filterSequence, 0, 16 * sizeof(int));
    memset(img->filterFrame, 0, 16 * sizeof(int));

    for(i = 0; i < 2; ++i)
    {
      memset(img->subpelOffset[i], 0, 16 * sizeof(int));
      memset(img->imgOffset[i], 0, MAX_REFERENCE_PICTURES * sizeof(int));
    }

    // Assign filters
    for(i = 0; i < 15; ++i)
      for(j = 0; j < SQR_FILTER; ++j)
      {
        SIFO_FILTER[0][i][j] = STANDARD_2D_FILTER[i][j];
#ifdef HPF_COMBO
        SIFO_FILTER[1][i][j] = LCDIF[i][j];
#else
        SIFO_FILTER[1][i][j] = SYMMETRIC_1[i][j];
#endif  // HPF_COMBO
        SIFO_FILTER[2][i][j] = SYMMETRIC_2[i][j];
      }
  }
}


//
//
//
static void UnifiedOneForthPixSwitchedFilters(int *filter, int *subpelOffset, int *imgOffset, int list)
{
  int i, j, offset;
  int frame;
  imgpel **out4Y;
  int img_width, img_height;
  int sub_pos, x_sub, y_sub;
  int icoord, jcoord;
  unsigned short ** imgY;

  for(frame = 0; frame < listXsize[list]; ++frame)
  {
    img_width  = listX[list][frame]->size_x;
    img_height = listX[list][frame]->size_y;

    out4Y  = listX[list][frame]->imgY_ups;
    imgY = listX[list][frame]->imgY;

    for(i = -4 * IMG_PAD_SIZE; i < (IMG_PAD_SIZE + img_height) * 4; ++i)
    {
      for(j = -4 * IMG_PAD_SIZE; j < (IMG_PAD_SIZE + img_width) * 4; ++j)
      {
        icoord = i + (IMG_PAD_SIZE << 2);
        jcoord = j + (IMG_PAD_SIZE << 2);

        x_sub = jcoord % 4;             // x-sub-coordinate in a 4x4block
        y_sub = icoord % 4;             // y-sub-coordinate in a 4x4block

        sub_pos = x_sub + (y_sub << 2);    // pos 0..15 in a 4x4 block

        offset = (frame == 0)? subpelOffset[sub_pos]: imgOffset[frame];

        assert(filter[sub_pos] < NUM_SIFO);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        out4Y[icoord][jcoord] = (int) Clip3Fun(0, (1 << img->bitdepth_luma) - 1, offset + 
          GetInterpolatedPixel(imgY, icoord, jcoord, img_width, img_height, sub_pos, filter[sub_pos]));
#else
        out4Y[icoord][jcoord] = (int) Clip3Fun(0, 255, offset + 
          GetInterpolatedPixel(imgY, icoord, jcoord, img_width, img_height, sub_pos, filter[sub_pos]));
#endif
      }
    }

    offset = (frame == 0)? subpelOffset[0]: imgOffset[frame];

    for(i = 0; i < img_height; ++i)
      for(j = 0; j < img_width; ++j)
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        listX[list][frame]->imgY_filt[i][j] = (int)Clip3(0, (1 << img->bitdepth_luma) - 1, listX[list][frame]->imgY[i][j] + offset);
#else
        listX[list][frame]->imgY_filt[i][j] = (int)Clip3(0, 255, listX[list][frame]->imgY[i][j] + offset);
#endif

    GenerateFullPelRepresentation(out4Y, listX[list][frame]->imgY_11_filt, img_width, img_height);  
  }
}


//
//
//
void UnifiedOneForthPixFiltSel(int filterParam)
{
  int zeroSubpelOffsets[16] = {0};
  int zeroImgOffset[MAX_REFERENCE_PICTURES] = {0};

  if(img->type == P_SLICE)
  {
    if(filterParam == SIFO_FIRST_PASS)
    {
      // SIFO_FIRST_PASS -> Sequence filters and zero offsets
      UnifiedOneForthPixSwitchedFilters(img->filterSequence, zeroSubpelOffsets, zeroImgOffset, LIST_0);
    }
    else if(filterParam == SIFO_FIRST_PASS_FPO)
    {
      // SIFO_FIRST_PASS_FPO -> Sequence filters and first pass offsets
      UnifiedOneForthPixSwitchedFilters(img->filterSequence, img->subpelOffset[LIST_0], img->imgOffset[LIST_0], LIST_0);
    }
    else if(filterParam == SIFO_FRAME_FILTER)
    {
      // SIFO_FRAME_FILTER -> Frame filters and zero offsets
      UnifiedOneForthPixSwitchedFilters(img->filterFrame, zeroSubpelOffsets, zeroImgOffset, LIST_0);
    }
    else
    {      
      // SIFO_FRAME_FILTER_WITH_OFFSET -> Frame filters and offsets
      UnifiedOneForthPixSwitchedFilters(img->filterFrame, img->subpelOffset[LIST_0], img->imgOffset[LIST_0], LIST_0);
    }
  }
  else if(img->type == B_SLICE)
  {
    if((filterParam == SIFO_FIRST_PASS) || (filterParam == SIFO_SEQ_FILTER))
    {
      // SIFO_FIRST_PASS || SIFO_SEQ_FILTER -> Sequence filters and zero offsets
      UnifiedOneForthPixSwitchedFilters(img->filterSequence, zeroSubpelOffsets, zeroImgOffset, LIST_0);
      UnifiedOneForthPixSwitchedFilters(img->filterSequence, zeroSubpelOffsets, zeroImgOffset, LIST_1);
    }
    else if(filterParam == SIFO_FIRST_PASS_FPO)
    {
      // SIFO_FIRST_PASS_FPO -> Sequence filters and first pass offsets
      UnifiedOneForthPixSwitchedFilters(img->filterSequence, img->subpelOffset[LIST_0], img->imgOffset[LIST_0], LIST_0);
      UnifiedOneForthPixSwitchedFilters(img->filterSequence, img->subpelOffset[LIST_1], img->imgOffset[LIST_1], LIST_1);
    }
    else if(filterParam == SIFO_FRAME_FILTER)
    {
      // SIFO_FRAME_FILTER -> Frame filters and zero offsets
      UnifiedOneForthPixSwitchedFilters(img->filterFrame, zeroSubpelOffsets, zeroImgOffset, LIST_0);
      UnifiedOneForthPixSwitchedFilters(img->filterFrame, zeroSubpelOffsets, zeroImgOffset, LIST_1);
    }
    else
    {
      // SIFO_FRAME_FILTER_WITH_OFFSET -> Frame filters and offsets
      UnifiedOneForthPixSwitchedFilters(img->filterFrame, img->subpelOffset[LIST_0], img->imgOffset[LIST_0], LIST_0);
      UnifiedOneForthPixSwitchedFilters(img->filterFrame, img->subpelOffset[LIST_1], img->imgOffset[LIST_1], LIST_1);
    }
  }
  else
  {
    assert(0);    // Called on a I-slice?
  }
}


//
//
//
static int TryImplicitWP(void)
{
  int td, tb;
  float b_scaling = 0.5;

  if(listXsize[LIST_0] && listXsize[LIST_1])
  {
    td = Clip3(-128, 127, (listX[LIST_1][0]->poc - listX[LIST_0][0]->poc));
    tb = Clip3(-128, 127, (enc_picture->poc - listX[LIST_0][0]->poc));
    if (td != 0)
      b_scaling = (float) tb / (float) td;

    if(b_scaling == 0.5)
      return 0;
    else
      return 1;
  }
  else
  {
    return 0;
  }
}


//
//  Used in HPF4
//
void rdPictureCodingFilterP_IPass(void)
{
  int previntras = intras;
  int rd_qp = img->qp;
  int prevtype = img->type;

  // printf("I Macroblocks %d (%2.0lf%%)\n", intras, (intras * 100.0 ) / img->FrameSizeInMbs);
  if(input->GenerateMultiplePPS && ((intras * 100 ) / img->FrameSizeInMbs >= 75))
  {
    img->qp = rd_qp - 1;
    active_pps = PicParSet[0];
    img->type = I_SLICE;
    img->filterParam = SIFO_INTRA;    // It shouldn't matter, img->filterParam is not used in I slices 

    frame_picture(frame_pic_2, 1);
    img->rd_pass = picture_coding_decision(frame_pic_1, frame_pic_2, rd_qp);

    // if(img->rd_pass)
    //   printf("I Macroblocks (%2.0lf%%), so Frame %d forced to Intra...\n", (previntras * 100.0 ) / img->FrameSizeInMbs, frame_no);
  }
  else
  {
    img->rd_pass = 0;
    enc_frame_picture2 = NULL;
  }
  
  if(img->rd_pass == 0)
  {
    enc_picture = enc_frame_picture;
    intras = previntras;
    frame_pic = frame_pic_1;
    img->qp = rd_qp;
    img->type = prevtype;
  }
  else
  {
    enc_picture = enc_frame_picture2;
    frame_pic = frame_pic_2;
  }
} 


//
//  Used in HPF3
//
void rdPictureCodingFilterP(void)
{
  int rd_qp = img->qp;
  int previntras = intras;
  int prevtype = img->type;  
  int IMblocks = 0;

  active_pps = PicParSet[0];

  // Compute frame filters
  ComputeFiltersAndOffsets();

  img->filterParam = SIFO_FRAME_FILTER_WITH_OFFSET;
  frame_picture (frame_pic_2, 1);
  img->rd_pass = picture_coding_decision(frame_pic_1, frame_pic_2, rd_qp);
  previntras = (img->rd_pass == 1)? intras: previntras;

  IMblocks = previntras;
  // printf("I Macroblocks %d (%2.0lf%%)\n", IMblocks, (IMblocks * 100.0 ) / img->FrameSizeInMbs);

  if(input->GenerateMultiplePPS && ((IMblocks * 100.0 ) / img->FrameSizeInMbs >= 75))
  {
    img->qp = rd_qp - 1;
    img->type = I_SLICE;
    img->filterParam = SIFO_INTRA;    // It shouldn't matter, img->filterParam is not used in I slices
    
    frame_picture(frame_pic_3, 2);

    if(img->rd_pass == 0)
      img->rd_pass = 2 * picture_coding_decision(frame_pic_1, frame_pic_3, rd_qp);
    else
      img->rd_pass += picture_coding_decision(frame_pic_2, frame_pic_3, rd_qp);

    // if(img->rd_pass == 2)
    //   printf("I Macroblocks (%2.0lf%%), so Frame %d forced to Intra...\n", (IMblocks * 100.0 ) / img->FrameSizeInMbs, frame_no);

  }

  if(img->rd_pass == 0)
  {
    enc_picture = enc_frame_picture;
    frame_pic = frame_pic_1;
    intras = previntras;
    img->qp = rd_qp;
    img->type = prevtype;
  }
  else if(img->rd_pass == 1)
  {
    enc_picture = enc_frame_picture2;
    frame_pic = frame_pic_2;
    intras = previntras;
    img->qp = rd_qp;
    img->type = prevtype;
  }       
  else
  {
    enc_picture = enc_frame_picture3;
    frame_pic = frame_pic_3;
  }

  // Update sequence filters
  UpdateSequenceFilters_P();
}


//
//
//
void rdPictureCodingFilterB(void)
{
  int rd_qp = img->qp;
  int previntras = intras;

  img->write_macroblock = 0;
  if(input->GenerateMultiplePPS && TryImplicitWP())
  {    
    // Implicit WP (POC) with sequence filters
    active_pps = PicParSet[2];
    img->filterParam = SIFO_SEQ_FILTER;
    frame_picture(frame_pic_2, 1);
    img->rd_pass = picture_coding_decision(frame_pic_1, frame_pic_2, rd_qp);
  }
  else
  {
    active_pps = PicParSet[0];
    img->rd_pass = 0;
    enc_frame_picture2 = NULL;
  }

  if(img->rd_pass == 0)
  {
    enc_picture = enc_frame_picture;
    intras = previntras;
    active_pps = PicParSet[0];
  }
  else
  {
    enc_picture = enc_frame_picture2;
    previntras = intras;
  }

  // Compute filters on Pass 0 or Pass 1
  ComputeFiltersAndOffsets();

  // Frame filters and sub-pel offsets
  img->write_macroblock = 0;
  img->filterParam = SIFO_FRAME_FILTER_WITH_OFFSET;
  frame_picture(frame_pic_3, 2);

  if(img->rd_pass == 0)
    img->rd_pass = 2 * picture_coding_decision(frame_pic_1, frame_pic_3, rd_qp);  // 0 or 2
  else
    img->rd_pass += picture_coding_decision(frame_pic_2, frame_pic_3, rd_qp);     // 1 or 2

  if(img->rd_pass == 0)
  {
    enc_picture = enc_frame_picture;
  }
  else if(img->rd_pass == 1)
  {
    enc_picture = enc_frame_picture2;
  }   
  else
  {
    enc_picture = enc_frame_picture3;
    previntras = intras;
  }

  // Filters and sub-pel offsets with QP+1
  img->write_macroblock = 0;
  img->qp = rd_qp + 1;
  img->filterParam = SIFO_FRAME_FILTER_WITH_OFFSET;

  if((img->rd_pass == 0) || (img->rd_pass == 1))
  {
    img->filterParam = SIFO_SEQ_FILTER;                 // Try sequence filters with QP+1
    frame_picture(frame_pic_3, 2);                      // Reuse frame_pic_3
  }
  else if(img->rd_pass == 2)
  {
    img->filterParam = SIFO_FRAME_FILTER_WITH_OFFSET;   // Try frame filters and offsets with QP+1
    frame_picture(frame_pic_2, 1);                      // Reuse frame_pic_2
  }   

  if(img->rd_pass == 0)
  {
    img->rd_pass = picture_coding_decision(frame_pic_1, frame_pic_3, rd_qp);
    if(img->rd_pass == 0)
    {
      enc_picture = enc_frame_picture;
      intras = previntras;
      img->qp = rd_qp;
      img->rd_pass = 0;
    }
    else
    {
      enc_picture = enc_frame_picture3;
      img->qp = rd_qp + 1;
      img->rd_pass = 2;
    }   
  }
  else if(img->rd_pass == 1)
  {
    img->rd_pass = picture_coding_decision(frame_pic_2, frame_pic_3, rd_qp);
    if(img->rd_pass == 0)
    {
      enc_picture = enc_frame_picture2;
      intras = previntras;
      img->qp = rd_qp;
      img->rd_pass = 1;
    }
    else
    {
      enc_picture = enc_frame_picture3;
      img->qp = rd_qp + 1;
      img->rd_pass = 2;
    }   
  }
  else
  {
    img->rd_pass = picture_coding_decision(frame_pic_3, frame_pic_2, rd_qp);
    if(img->rd_pass == 0)
    {
      enc_picture = enc_frame_picture3;
      intras = previntras;
      img->qp = rd_qp;
      img->rd_pass = 2;
    }
    else
    {
      enc_picture = enc_frame_picture2;
      img->qp = rd_qp + 1;      
      img->rd_pass = 1;
    }
  } 

  // Update sequence filters
  UpdateSequenceFilters_B();
}



//
//  Return whether any offset is non zero
//
int getNonzero(int list)
{
  int sub_pos;
  int frame;
  for(frame = 0; frame < listXsize[list]; ++frame)
  {       
    if(frame == 0)
    {     
      for(sub_pos = 0; sub_pos < 16; ++sub_pos)   
        if(img->subpelOffset[list][sub_pos] != 0)
          return 1;
    }
    else
    {
      if(img->imgOffset[list][frame] != 0)
        return 1;
    }
  }
  return 0;
}

#endif




