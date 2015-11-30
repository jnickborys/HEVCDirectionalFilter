#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <assert.h>

#include <time.h>
#include <sys/timeb.h>

#include "global.h"
#include "image.h"

#include "refbuf.h"
#include "memalloc.h"
#include "mb_access.h"
#include "fast_me.h"
#include "simplified_fast_me.h"

#ifdef SWITCHED_FILTERS
#include "switched_filters.h"
#endif

#ifdef MV_COMPETITION
#include "mv_competition.h"
extern MV_Competition mv_comp;
#endif

#include "epzs.h"
#define EPZSREF 1 //also defined in epzs.c

#define SIMPLIFY_CODE

#ifdef MB32X32

void get_mb_pos (int mb_addr, int *x, int*y);
int computeBiPredSad2(pel_t** cur_pic,
                             int blocksize_y,
                             int blocksize_x, 
                             int blockshape_x,
                             int mcost,
                             int min_mcost,
                             int cand_x1, int cand_y1, 
                             int cand_x2, int cand_y2);
int computeBiPredSad1(pel_t** cur_pic,
                             int blocksize_y,
                             int blocksize_x, 
                             int blockshape_x,
                             int mcost,
                             int min_mcost,
                             int cand_x1, int cand_y1, 
                             int cand_x2, int cand_y2);
short EPZSSpatialPredictors (PixelPos block_a, 
                                    PixelPos block_b,
                                    PixelPos block_c, 
                                    PixelPos block_d,
                                    int list, 
                                    int list_offset, 
                                    short ref,
                                    char **refPic, 
                                    short ***tmp_mv,
                                    EPZSStructure * predictor);
void EPZSSpatialMemPredictors (int list, 
                                      short ref,
                                      int blocktype,
                                      int pic_x, 
                                      int bs_x, 
                                      int bs_y, 
                                      int by,
                                      int *prednum,
                                      int img_width,
                                      EPZSStructure * predictor);
void
EPZSTemporalPredictors (int list,         // <--  current list
                        int list_offset,  // <--  list offset for MBAFF
                        short ref,        // <--  current reference frame
                        int o_block_x,  // <--  absolute x-coordinate of regarded AxB block
                        int o_block_y,  // <--  absolute y-coordinate of regarded AxB block
                        EPZSStructure * predictor, 
                        int *prednum,
                        int block_available_left, 
                        int block_available_up,
                        int block_available_right, 
                        int block_available_below,
                        int blockshape_x, 
                        int blockshape_y,
                        int stopCriterion, 
                        int min_mcost);
void EPZSWindowPredictors (int mv_x, int mv_y, EPZSStructure *predictor, int *prednum, int extended);
void EPZSBlockTypePredictors (int block_x, int block_y, int blocktype, int ref, int list,
                                     EPZSStructure * predictor, int *prednum);
int                                            //  ==> minimum motion cost after search
EPZSPelBlockMotionSearch32 (pel_t ** cur_pic,    // <--  original pixel values for the AxB block
                          short ref,          // <--  reference picture 
                          int list,           // <--  reference list
                          int list_offset,    // <--  offset for Mbaff
                          char ***refPic,    // <--  reference array
                          short ****tmp_mv,   // <--  mv array
                          int pic_pix_x,      // <--  absolute x-coordinate of regarded AxB block
                          int pic_pix_y,      // <--  absolute y-coordinate of regarded AxB block
                          int blocktype,      // <--  block type (1-16x16 ... 7-4x4)
                          short pred_mv_x,    // <--  motion vector predictor (x) in sub-pel units
                          short pred_mv_y,    // <--  motion vector predictor (y) in sub-pel units
                          short *mv_x,        // <--> in: search center (x) / out: motion vector (x) - in pel units
                          short *mv_y,        // <--> in: search center (y) / out: motion vector (y) - in pel units
                          int search_range,    // <--  1-d search range in pel units
                          int min_mcost,      // <--  minimum motion cost (cost for center or huge value)
                          int lambda_factor,      // <--  lagrangian parameter for determining motion cost
                          int max_ext_level);
int                                               //  ==> minimum motion cost after search
EPZSSubPelBlockMotionSearch32 (pel_t**   orig_pic,      // <--  original pixel values for the AxB block
                             short     ref,           // <--  reference frame (0... or -1 (backward))
                             int       list,          // <--  reference picture list 
                             int       pic_pix_x,     // <--  absolute x-coordinate of regarded AxB block
                             int       pic_pix_y,     // <--  absolute y-coordinate of regarded AxB block
                             int       blocktype,     // <--  block type (1-16x16 ... 7-4x4)
                             int       pred_mv_x,     // <--  motion vector predictor (x) in sub-pel units
                             int       pred_mv_y,     // <--  motion vector predictor (y) in sub-pel units
                             short*    mv_x,          // <--> in: search center (x) / out: motion vector (x) - in pel units
                             short*    mv_y,          // <--> in: search center (y) / out: motion vector (y) - in pel units
                             int       search_pos2,   // <--  search positions for    half-pel search  (default: 9)
                             int       search_pos4,   // <--  search positions for quarter-pel search  (default: 9)
#ifdef EIGHTH_PEL
                             int       search_pos8,   // <--  search positions for eighth-pel search  (default: 9)
#endif
                             int       min_mcost,     // <--  minimum motion cost (cost for center or huge value)
                             int       lambda_factor,  // <--  lagrangian parameter for determining motion cost
                             int       max_ext_level);
int                                                //  ==> minimum motion cost after search
EPZSBiPredBlockMotionSearch32 (pel_t ** cur_pic,    // <--  original pixel values for the AxB block
                             short  ref,          // <--  reference picture 
                             int    list,         // <--  reference list
                             int    list_offset,  // <--  offset for Mbaff
                             char  ***refPic,    // <--  reference array
                             short  ****tmp_mv,   // <--  mv array
                             int    pic_pix_x,    // <--  absolute x-coordinate of regarded AxB block
                             int    pic_pix_y,    // <--  absolute y-coordinate of regarded AxB block
                             int    blocktype,    //<--  block type (1-16x16 ... 7-4x4)
                             short  pred_mv_x1,   // <--  motion vector predictor (x) in sub-pel units
                             short  pred_mv_y1,   // <--  motion vector predictor (y) in sub-pel units
                             short  pred_mv_x2,   // <--  motion vector predictor (x) in sub-pel units
                             short  pred_mv_y2,   // <--  motion vector predictor (y) in sub-pel units
                             short  *mv_x,        // <--> in: search center (x) / out: motion vector (x) - in pel units
                             short  *mv_y,        // <--> in: search center (y) / out: motion vector (y) - in pel units
                             short  *s_mv_x,      // <--> in: search center (x) / out: motion vector (x) - in pel units
                             short  *s_mv_y,      // <--> in: search center (y) / out: motion vector (y) - in pel units
                             int    search_range,  // <--  1-d search range in pel units
                             int    min_mcost,    // <--  minimum motion cost (cost for center or huge value)
                             int    lambda_factor,        // <--  lagrangian parameter for determining motion cost
                             int    max_ext_level);


int computeSad(pel_t** cur_pic,
                      int blocksize_y,
                      int blocksize_x, 
                      int blockshape_x,
                      int mcost,
                      int min_mcost,
                      int cand_x,
                      int cand_y);
void CheckAvailabilityOfNeighbors32(int mb32);

static pel_t *(*get_line_p1) (pel_t**, int, int, int, int);
static pel_t *(*get_line_p2) (pel_t**, int, int, int, int);


extern StorablePicture **listX[6];
extern int**** motion_cost;
#if defined MB32X32 && defined RDO_Q
extern int motion32_cost4[5][8][2][2][4];
extern int motion32_cost8[5][8][2][2][4];

#endif
extern int curr_mb32i;
extern ColocatedParams *Co_located;
extern short*  spiral_search_x;
extern short*  spiral_search_y;
extern short*  spiral_hpel_search_x;
extern short*  spiral_hpel_search_y;
extern int*    mvbits;
extern int medthres[8];
extern EPZSStructure *searchPattern, *searchPatternD, *predictor;
extern pel_t *ref_pic;
extern pel_t *ref_pic1;
extern pel_t *ref_pic2;

#if EPZSREF
extern short ******EPZSMotion;  //!< Array for storing Motion Vectors
#else
extern short *****EPZSMotion;  //!< Array for storing Motion Vectors
#endif
extern const  int LEVELMVLIMIT[17][6];


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
void SetMotionVectorPredictor32 (short  pmv[2],
                               char   **refPic,
                               short  ***tmp_mv,
                               short  ref_frame,
                               int    list,
                               int    block_x,
                               int    block_y,
                               int    blockshape_x,
                               int    blockshape_y,
                               int    mb_ext_level)
{
  int mb_x                 = 4*block_x;
  int mb_y                 = 4*block_y;
  //int mb_nr                = img->current_mb_nr;
  
  int mv_a=0, mv_b=0, mv_c=0, pred_vec=0;
  int mvPredType, rFrameL=0, rFrameU=0, rFrameUR=0;
  int hv;
  
  PixelPos block_a, block_b, block_c, block_d;
    
#ifndef SIMPLIFY_CODE
  int dsr_temp_search_range[2];
  int dsr_mv_avail, dsr_mv_max, dsr_mv_sum, dsr_small_search_range;
  // neighborhood SAD init
  if(input->FMEnable == 1) 
  {
    SAD_a=0;
    SAD_b=0; 
    SAD_c=0; 
    SAD_d=0;
  }
#endif




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
    if (mb_x < (8<<mb_ext_level))  // first column of 16x16 blocks
    {
      if (mb_y==(8<<mb_ext_level))
      {//2nd partition in 32x16 mode
        if (blockshape_x == (16<<mb_ext_level))      block_c.available  = 0;
      }
      else
      {
        error("never reach here!\n", -1);
        if (mb_x+blockshape_x == 16)  block_c.available = 0;
      }
    }
    else
    {
      error("never reach here!\n", -1);
      if (mb_x+blockshape_x == 32)   block_c.available = 0;
    }
  }
  
  if (!block_c.available)
  {
    block_c=block_d;
  }
  
  mvPredType = MVPRED_MEDIAN;
  
  if (!img->MbaffFrameFlag)
  {
    rFrameL    = block_a.available    ? refPic[block_a.pos_y][block_a.pos_x] : -1;
    rFrameU    = block_b.available    ? refPic[block_b.pos_y][block_b.pos_x] : -1;
    rFrameUR   = block_c.available    ? refPic[block_c.pos_y][block_c.pos_x] : -1;
  }
#ifndef SIMPLIFY_CODE
  else
  {
    if (img->mb_data[img->current_mb_nr].mb_field)
    {
      rFrameL    = block_a.available    ? 
        img->mb_data[block_a.mb_addr].mb_field ? 
        refPic[block_a.pos_y][block_a.pos_x]:
      refPic[block_a.pos_y][block_a.pos_x] * 2: 
      -1;
      rFrameU    = block_b.available    ? 
        img->mb_data[block_b.mb_addr].mb_field ? 
        refPic[block_b.pos_y][block_b.pos_x]:
      refPic[block_b.pos_y][block_b.pos_x] * 2: 
      -1;
      rFrameUR    = block_c.available    ? 
        img->mb_data[block_c.mb_addr].mb_field ? 
        refPic[block_c.pos_y][block_c.pos_x]:
      refPic[block_c.pos_y][block_c.pos_x] * 2: 
      -1;
    }
    else
    {
      rFrameL    = block_a.available    ? 
        img->mb_data[block_a.mb_addr].mb_field ? 
        refPic[block_a.pos_y][block_a.pos_x] >>1:
      refPic[block_a.pos_y][block_a.pos_x] : 
      -1;
      rFrameU    = block_b.available    ? 
        img->mb_data[block_b.mb_addr].mb_field ? 
        refPic[block_b.pos_y][block_b.pos_x] >>1:
      refPic[block_b.pos_y][block_b.pos_x] : 
      -1;
      rFrameUR    = block_c.available    ? 
        img->mb_data[block_c.mb_addr].mb_field ? 
        refPic[block_c.pos_y][block_c.pos_x] >>1:
      refPic[block_c.pos_y][block_c.pos_x] : 
      -1;
    }
  }
#endif
  /* Prediction if only one of the neighbors uses the reference frame
  * we are checking
  */
  if(rFrameL == ref_frame && rFrameU != ref_frame && rFrameUR != ref_frame)       mvPredType = MVPRED_L;
  else if(rFrameL != ref_frame && rFrameU == ref_frame && rFrameUR != ref_frame)  mvPredType = MVPRED_U;
  else if(rFrameL != ref_frame && rFrameU != ref_frame && rFrameUR == ref_frame)  mvPredType = MVPRED_UR;
  // Directional predictions 

  if(blockshape_x == (8<<mb_ext_level) && blockshape_y == (16<<mb_ext_level)) 
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

  else if(blockshape_x == (16<<mb_ext_level) && blockshape_y == (8<<mb_ext_level)) 
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
#ifndef SIMPLIFY_CODE  
  // neighborhood SAD prediction
  if(input->FMEnable == 1 &&(input->DynamicSearchRange == 1 || input->BiPredMotionEstimation == 1))
  {
    SAD_a = block_a.available ? ((list==1) ? (fastme_l1_cost_flag[FME_blocktype][block_a.pos_y][block_a.pos_x]) : (fastme_l0_cost_flag[FME_blocktype][block_a.pos_y][block_a.pos_x])) : 0;
    SAD_b = block_b.available ? ((list==1) ? (fastme_l1_cost_flag[FME_blocktype][block_b.pos_y][block_b.pos_x]) : (fastme_l0_cost_flag[FME_blocktype][block_b.pos_y][block_b.pos_x])) : 0;
    SAD_d = block_d.available ? ((list==1) ? (fastme_l1_cost_flag[FME_blocktype][block_d.pos_y][block_d.pos_x]) : (fastme_l0_cost_flag[FME_blocktype][block_d.pos_y][block_d.pos_x])) : 0;
    SAD_c = block_c.available ? ((list==1) ? (fastme_l1_cost_flag[FME_blocktype][block_c.pos_y][block_c.pos_x]) : (fastme_l0_cost_flag[FME_blocktype][block_c.pos_y][block_c.pos_x])) : SAD_d;
  }
#endif
  for (hv=0; hv < 2; hv++)
  {
    if (!img->MbaffFrameFlag || hv==0)
    {
      mv_a = block_a.available  ? tmp_mv[block_a.pos_y][block_a.pos_x][hv] : 0;
      mv_b = block_b.available  ? tmp_mv[block_b.pos_y][block_b.pos_x][hv] : 0;
      mv_c = block_c.available  ? tmp_mv[block_c.pos_y][block_c.pos_x][hv] : 0;
    }
#ifndef SIMPLIFY_CODE
    else
    {
      if (img->mb_data[img->current_mb_nr].mb_field)
      {
        mv_a = block_a.available  ? img->mb_data[block_a.mb_addr].mb_field
          ? tmp_mv[block_a.pos_y][block_a.pos_x][hv]
          : tmp_mv[block_a.pos_y][block_a.pos_x][hv] / 2
          : 0;
        mv_b = block_b.available  ? img->mb_data[block_b.mb_addr].mb_field
          ? tmp_mv[block_b.pos_y][block_b.pos_x][hv]
          : tmp_mv[block_b.pos_y][block_b.pos_x][hv] / 2
          : 0;
        mv_c = block_c.available  ? img->mb_data[block_c.mb_addr].mb_field
          ? tmp_mv[block_c.pos_y][block_c.pos_x][hv]
          : tmp_mv[block_c.pos_y][block_c.pos_x][hv] / 2
          : 0;
      }
      else
      {
        mv_a = block_a.available  ? img->mb_data[block_a.mb_addr].mb_field
          ? tmp_mv[block_a.pos_y][block_a.pos_x][hv] * 2
          : tmp_mv[block_a.pos_y][block_a.pos_x][hv]
          : 0;
        mv_b = block_b.available  ? img->mb_data[block_b.mb_addr].mb_field
          ? tmp_mv[block_b.pos_y][block_b.pos_x][hv] * 2
          : tmp_mv[block_b.pos_y][block_b.pos_x][hv]
          : 0;
        mv_c = block_c.available  ? img->mb_data[block_c.mb_addr].mb_field
          ? tmp_mv[block_c.pos_y][block_c.pos_x][hv] * 2
          : tmp_mv[block_c.pos_y][block_c.pos_x][hv]
          : 0;
      }
    }
#endif
    switch (mvPredType)
    {
    case MVPRED_MEDIAN:
      if(!(block_b.available || block_c.available))
      {
        pred_vec = mv_a;
      }
      else
      {
        pred_vec = mv_a+mv_b+mv_c-min(mv_a,min(mv_b,mv_c))-max(mv_a,max(mv_b,mv_c));
      }
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
    
    pmv[hv] = pred_vec;


#ifndef SIMPLIFY_CODE
    //Dynamic Search Range
    if (input->FMEnable == 1 && input->DynamicSearchRange)
    {
      dsr_mv_avail=block_a.available+block_b.available+block_c.available;
      if(dsr_mv_avail < 2)
      {
        dsr_temp_search_range[hv] = input->search_range;
      }
      else
      {
        dsr_mv_max = max(abs(mv_a),max(abs(mv_b),abs(mv_c)));
        dsr_mv_sum = (abs(mv_a)+abs(mv_b)+abs(mv_c));
        if(dsr_mv_sum == 0) dsr_small_search_range = (input->search_range + 4) >> 3;
        else if(dsr_mv_sum > 3 ) dsr_small_search_range = (input->search_range + 2) >>2;
        else dsr_small_search_range = (3*input->search_range + 8) >> 4;
        dsr_temp_search_range[hv]=min(input->search_range,max(dsr_small_search_range,dsr_mv_max<<1));    
        if(max(SAD_a,max(SAD_b,SAD_c)) > Threshold_DSR_MB[FME_blocktype])
          dsr_temp_search_range[hv] = input->search_range;
      }       
    }
#endif
  }
#ifndef SIMPLIFY_CODE  
  //Dynamic Search Range
  if (input->FMEnable == 1 && input->DynamicSearchRange)
    dsr_new_search_range = max(dsr_temp_search_range[0],dsr_temp_search_range[1]);
#endif
}

/*!
***********************************************************************
* \brief
*    Calculate SA(T)D for 8x8
***********************************************************************
*/
int
find_SATD32 (int c_diff[MB_PIXELS], int blocktype, int mb_ext_level)
{
  int i, sad=0;
  
  switch(blocktype)
  {
    //32x32
  case 1: 

    for(i=0; i<((2<<mb_ext_level)*(2<<mb_ext_level)); i++)
      sad += SATD8X8 (&c_diff[i*64],       input->hadamard);
    break;
    //32x16 16x32
  case 2:
  case 3: 

    for(i=0; i<((2<<mb_ext_level)*(1<<mb_ext_level)); i++)
      sad += SATD8X8 (&c_diff[i*64],       input->hadamard);
    break;
#ifndef SIMPLIFY_CODE
    //8x8
  case 4: 
    sad  = SATD8X8 (c_diff, input->hadamard);
    break;
#endif
    //8x4 4x8
  default:
    sad=-1;
    break;
  }
  
  return sad;
}

/*!
***********************************************************************
* \brief
*    Bipred Sub pixel block motion search
***********************************************************************
*/
int                                               //  ==> minimum motion cost after search
SubPelBlockSearchBiPred32 (pel_t**   orig_pic,      // <--  original pixel values for the AxB block
                         short     ref,           // <--  reference frame (0... or -1 (backward))
                         int       list,          // <--  reference picture list 
                         int       pic_pix_x,     // <--  absolute x-coordinate of regarded AxB block
                         int       pic_pix_y,     // <--  absolute y-coordinate of regarded AxB block
                         int       blocktype,     // <--  block type (1-16x16 ... 7-4x4)
                         short     pred_mv_x,     // <--  motion vector predictor (x) in sub-pel units
                         short     pred_mv_y,     // <--  motion vector predictor (y) in sub-pel units
                         short*    mv_x,          // <--> in: search center (x) / out: motion vector (x) - in pel units
                         short*    mv_y,          // <--> in: search center (y) / out: motion vector (y) - in pel units
                         short*    s_mv_x,          // <--> in: search center (x) / out: motion vector (x) - in pel units
                         short*    s_mv_y,          // <--> in: search center (y) / out: motion vector (y) - in pel units
                         int       search_pos2,   // <--  search positions for    half-pel search  (default: 9)
                         int       search_pos4,   // <--  search positions for quarter-pel search  (default: 9)
#ifdef EIGHTH_PEL
                         int       search_pos8,   // <--  search positions for eighth-pel search   (default: 9)
#endif
                         int       min_mcost,     // <--  minimum motion cost (cost for center or huge value)
                         int       lambda_factor,   // <--  lagrangian parameter for determining motion cost
                         int       mb_ext_level)
{
  int   j, i, k;
  int   c_diff[MB_PIXELS<<(2*MAX_MB_EXT_LEVEL)];
  int   diff[16], *d;  
  int   pos, best_pos, mcost, abort_search;
  int   y0, y1, y2, y3;
  int   x0;
  int   ry0, ry4, ry8, ry12, rx0;
  int   sy0, sy4, sy8, sy12, sx0;
  
  int   cand_mv_x, cand_mv_y;
#ifdef SWITCHED_FILTERS
  int   best_pos_sp16 = 0;
  int   min_mcost_sp16 = min_mcost;
  short mv_x_sp16 = 0;
  short mv_y_sp16 = 0;
#endif

  int   blocksize_x     = input->blc_size[blocktype][0]<<mb_ext_level;
  int   blocksize_y     = input->blc_size[blocktype][1]<<mb_ext_level;
  int   pic4_pix_x      = ((pic_pix_x + IMG_PAD_SIZE)<< 2);
  int   pic4_pix_y      = ((pic_pix_y + IMG_PAD_SIZE)<< 2);
  int   min_pos2        = (input->hadamard ? 0 : 1);
  int   max_pos2        = (input->hadamard ? max(1,search_pos2) : search_pos2);
  int   list_offset   = img->mb_data[img->current_mb_nr].list_offset; 
  int   apply_weights =  (active_pps->weighted_bipred_idc );  
  short weightSpic = (apply_weights ? (list == 0? wbp_weight[list_offset    ][ref][0  ][0]: wbp_weight[list_offset + 1][0  ][ref][0]) : 1);
  short weightRpic = (apply_weights ? (list == 0? wbp_weight[list_offset + 1][ref][0  ][0]: wbp_weight[list_offset    ][0  ][ref][0]) : 1);
  short offsetSpic = (apply_weights ? (list == 0?  wp_offset[list_offset    ][ref]     [0]:  wp_offset[list_offset + 1][0  ]     [0]) : 0);
  short offsetRpic = (apply_weights ? (list == 0?  wp_offset[list_offset + 1][ref]     [0]:  wp_offset[list_offset    ][0  ]     [0]) : 0);

#ifdef USE_HP_FILTER//BROUND
  short offsetBi;
#else
  short offsetBi=(offsetRpic + offsetSpic + 1)>>1;
#endif

  pel_t weightedpel;
  int   test8x8transform = input->Transform8x8Mode && blocktype <= 4 && input->hadamard;
  int   cmv_x, cmv_y;
  int   smv_x = *s_mv_x + pic4_pix_x;
  int   smv_y = *s_mv_y + pic4_pix_y;
  
  StorablePicture *ref_picture = listX[list+list_offset][ref];
  pel_t **ref1_pic = ref_picture->imgY_ups;      
  pel_t **ref2_pic = listX[list==0? 1 +list_offset: list_offset][0  ]->imgY_ups;  
  pel_t *ref_line_p1,*ref_line_p2;
  pel_t *orig_line;
  int img_width  = ((ref_picture->size_x + 2*IMG_PAD_SIZE - 1)<<2);
  int img_height = ((ref_picture->size_y + 2*IMG_PAD_SIZE - 1)<<2);
  int max_pos_x4 = ((ref_picture->size_x - blocksize_x + 2*IMG_PAD_SIZE)<<2);
  int max_pos_y4 = ((ref_picture->size_y - blocksize_y + 2*IMG_PAD_SIZE)<<2);
#ifdef  USE_HP_FILTER//BROUND
  int   round;
  img->bipred_rounding_control = (img->nal_reference_idc!=0);
  round = 1 - img->bipred_rounding_control;
  if(input->UseHPFilter != 0)
  {
    offsetBi=offsetRpic + offsetSpic;
  }
  else
  {
    offsetBi=(offsetRpic + offsetSpic + 1)>>1;
  }
#endif
  
  /*********************************
  *****                       *****
  *****  HALF-PEL REFINEMENT  *****
  *****                       *****
  *********************************/
  
  //===== set function for getting pixel values =====
  if ((pic4_pix_x + *mv_x > 2) && (pic4_pix_x + *mv_x < max_pos_x4 - 1) &&
    (pic4_pix_y + *mv_y > 2) && (pic4_pix_y + *mv_y < max_pos_y4 - 1))
  {
    get_line_p2 = FastLine4X;
  } 
  else
  {
    get_line_p2 = UMVLine4X2;
  }
  
  if ((pic4_pix_x + *s_mv_x > 2) && (pic4_pix_x + *s_mv_x < max_pos_x4 - 1) &&
    (pic4_pix_y + *s_mv_y > 2) && (pic4_pix_y + *s_mv_y < max_pos_y4 - 1))
  {
    get_line_p1 = FastLine4X;
  }
  else
  {
    get_line_p1 = UMVLine4X;    
  }
  
  
  //===== loop over search positions =====
#ifdef SWITCHED_FILTERS
  for (best_pos = 0, pos = min_pos2; 
    pos < max_pos2 + (((img->filterParam == SIFO_FIRST_PASS_FPO) && (ref == 0)) ? img->noOffsets[list]: 0); 
    pos++)
  {
    if(pos < max_pos2)
    {
      cand_mv_x = *mv_x + (spiral_hpel_search_x[pos]);    // quarter-pel units
      cand_mv_y = *mv_y + (spiral_hpel_search_y[pos]);    // quarter-pel units
    }
    else
    {
      cand_mv_x = *mv_x + img->search_point_qp_x_offset[pos-max_pos2];    // quarter-pel units
      cand_mv_y = *mv_y + img->search_point_qp_y_offset[pos-max_pos2];    // quarter-pel units
    }
#else
  for (best_pos = 0, pos = min_pos2; pos < max_pos2; pos++)
  {
    cand_mv_x = *mv_x + (spiral_hpel_search_x[pos]);    // quarter-pel units
    cand_mv_y = *mv_y + (spiral_hpel_search_y[pos]);    // quarter-pel units
#endif
    
    //----- set motion vector cost -----
#ifdef MB32X32_MVC
    if (input->mv_competition > 0)
      mcost = compute_mv_cost_for_each_predictor (lambda_factor, 0, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y, ref, list ^ 1, blocktype, pic_pix_x-img->opix_x, pic_pix_y-img->opix_y);
    else
#endif
      mcost = MV_COST (lambda_factor, 0, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);
    
    if (mcost >= min_mcost) continue;
    
    cmv_x = cand_mv_x + pic4_pix_x;
    cmv_y = cand_mv_y + pic4_pix_y;
    
    //----- add up SATD -----
    for (y0=0, abort_search=0; y0<blocksize_y && !abort_search; y0+=BLOCK_SIZE)
    {
      ry0 = (y0<<2) + cmv_y;
      ry4  = ry0 + 4;
      ry8  = ry4 + 4;
      ry12 = ry8 + 4;      
      sy0 = (y0<<2) + smv_y;
      sy4  = sy0 + 4;
      sy8  = sy4 + 4;
      sy12 = sy8 + 4;
      y1 = y0 + 1;
      y2 = y1 + 1;
      y3 = y2 + 1;
      
      if (apply_weights)
      {
        for (x0=0; x0<blocksize_x; x0+=BLOCK_SIZE)
        {
          rx0 = (x0<<2) + cmv_x;            
          sx0 = (x0<<2) + smv_x;
          d   = diff;
          

#ifdef  USE_HP_FILTER//BROUND
          if(input->UseHPFilter != 0)
          {
          orig_line = &orig_pic [y0][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy0, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry0, rx0,  img_height, img_width);
          weightedpel = clip1a (((( weightSpic * (*(ref_line_p1    )) + weightRpic * (*(ref_line_p2     ))) 
            +(offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
             + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)));
          *d++ = *orig_line++ - weightedpel;          
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line - weightedpel;
          
          orig_line = &orig_pic [y1][x0];
          ref_line_p1 = get_line_p1 (ref1_pic, sy4, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry4, rx0,  img_height, img_width);
          weightedpel = clip1a (((( weightSpic * (*(ref_line_p1    )) + weightRpic * (*(ref_line_p2     ))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;          
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line - weightedpel;
          
          
          orig_line = &orig_pic [y2][x0];
          ref_line_p1 = get_line_p1 (ref1_pic, sy8, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry8, rx0,  img_height, img_width);
          weightedpel = clip1a (((( weightSpic * (*(ref_line_p1    )) + weightRpic * (*(ref_line_p2     ))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;          
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line - weightedpel;
          
          orig_line = &orig_pic [y3][x0];
          ref_line_p1 = get_line_p1 (ref1_pic, sy12, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry12, rx0,  img_height, img_width);
          weightedpel = clip1a (((( weightSpic * (*(ref_line_p1    )) + weightRpic * (*(ref_line_p2     ))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;          
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d   = *orig_line - weightedpel;
          }
          else
          {
#endif

          orig_line = &orig_pic [y0][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy0, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry0, rx0,  img_height, img_width);
          weightedpel = clip1a (((( weightSpic * (*(ref_line_p1    )) + weightRpic * (*(ref_line_p2     ))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;          
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line - weightedpel;
          
          orig_line = &orig_pic [y1][x0];
          ref_line_p1 = get_line_p1 (ref1_pic, sy4, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry4, rx0,  img_height, img_width);
          weightedpel = clip1a (((( weightSpic * (*(ref_line_p1    )) + weightRpic * (*(ref_line_p2     ))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;          
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line - weightedpel;
          
          
          orig_line = &orig_pic [y2][x0];
          ref_line_p1 = get_line_p1 (ref1_pic, sy8, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry8, rx0,  img_height, img_width);
          weightedpel = clip1a (((( weightSpic * (*(ref_line_p1    )) + weightRpic * (*(ref_line_p2     ))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;          
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line - weightedpel;
          
          orig_line = &orig_pic [y3][x0];
          ref_line_p1 = get_line_p1 (ref1_pic, sy12, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry12, rx0,  img_height, img_width);
          weightedpel = clip1a (((( weightSpic * (*(ref_line_p1    )) + weightRpic * (*(ref_line_p2     ))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;          
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d   = *orig_line - weightedpel;
#ifdef USE_HP_FILTER//BROUND
          }
#endif
          
          if (!test8x8transform)
          {
            if ((mcost += SATD (diff, input->hadamard)) > min_mcost)
            {
              abort_search = 1;
              break;
            }
          }
          else
          {
            int x8=x0/8, y8=y0/8;
            i = (x0&0x7) +  x8*64 + (y8 * blocksize_x/8)*64 ;
            for(k=0, j=y0; j<4 + y0; j++, k+=4)
              memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], 4*sizeof(int));
          }
        }
      }
      else
      {   
#ifdef SWITCHED_FILTERS
        img->bipred_rounding_control = (img->nal_reference_idc != 0);
        round = img->currMEOffset[0] + img->currMEOffset[1] + (1 - img->bipred_rounding_control);
#endif
        for (x0=0; x0<blocksize_x; x0+=BLOCK_SIZE)
        {         
          rx0 = (x0<<2) + cmv_x;            
          sx0 = (x0<<2) + smv_x;
          d   = diff;
          
#ifdef  USE_HP_FILTER//BROUND
          if(input->UseHPFilter != 0)
          {
          orig_line = &orig_pic [y0][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy0, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry0, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  ((*(ref_line_p1     ) + *(ref_line_p2     ) + round) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          
          orig_line = &orig_pic [y1][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy4, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry4, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  ((*(ref_line_p1     ) + *(ref_line_p2     ) + round) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          
          orig_line = &orig_pic [y2][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy8, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry8, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  ((*(ref_line_p1     ) + *(ref_line_p2     ) + round) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          
          orig_line = &orig_pic [y3][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy12, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry12, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  ((*(ref_line_p1     ) + *(ref_line_p2     ) + round) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d   = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          }
          else
          {
#endif
          orig_line = &orig_pic [y0][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy0, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry0, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  ((*(ref_line_p1     ) + *(ref_line_p2     ) + 1) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          *d++ = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          
          orig_line = &orig_pic [y1][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy4, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry4, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  ((*(ref_line_p1     ) + *(ref_line_p2     ) + 1) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          *d++ = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          
          orig_line = &orig_pic [y2][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy8, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry8, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  ((*(ref_line_p1     ) + *(ref_line_p2     ) + 1) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          *d++ = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          
          orig_line = &orig_pic [y3][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy12, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry12, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  ((*(ref_line_p1     ) + *(ref_line_p2     ) + 1) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          *d   = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
#ifdef  USE_HP_FILTER//BROUND
          }
#endif
          
          if (!test8x8transform)
          {
            if ((mcost += SATD (diff, input->hadamard)) > min_mcost)
            {
              abort_search = 1;
              break;
            }
          }
          else
          {
            int x8=x0/8, y8=y0/8;
            i = (x0&0x7) +  x8*64 + (y8 * blocksize_x/8)*64 ;
            for(k=0, j=y0; j<4 + y0; j++, k+=4)
              memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], 4*sizeof(int));
          }
        }
      }
    }  
    
    if(test8x8transform)
      mcost += find_SATD32 (c_diff, blocktype, mb_ext_level);   
    
#ifdef SWITCHED_FILTERS
    if((img->filterParam == SIFO_FIRST_PASS_FPO) || (pos < max_pos2))
    {
      // Update cost for regular search
      if (mcost < min_mcost)
      {
        min_mcost = mcost;
        best_pos  = pos;
      }      
    }
    else
    {
      // Update results of offset search
      if(mcost < min_mcost_sp16)
      {
        min_mcost_sp16 = mcost;
        best_pos_sp16  = pos;
        mv_x_sp16 = *mv_x + img->search_point_qp_x_offset[pos - max_pos2];
        mv_y_sp16 = *mv_y + img->search_point_qp_y_offset[pos - max_pos2];
      }  
    }
#else
    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      best_pos  = pos;
    }      
#endif
  }
  
  if (best_pos)
  {
    *mv_x += (spiral_hpel_search_x [best_pos]);
    *mv_y += (spiral_hpel_search_y [best_pos]);
  }
  
  test8x8transform = input->Transform8x8Mode && blocktype <= 4 && input->hadamard;
  
  /************************************
  *****                          *****
  *****  QUARTER-PEL REFINEMENT  *****
  *****                          *****
  ************************************/
  //===== set function for getting pixel values =====
  if ((pic4_pix_x + *mv_x > 0) && (pic4_pix_x + *mv_x < max_pos_x4) &&
    (pic4_pix_y + *mv_y > 0) && (pic4_pix_y + *mv_y < max_pos_y4))
  {
    get_line_p2 = FastLine4X;
  }
  else
  {
    get_line_p2 = UMVLine4X2;    
  }
  
  if ((pic4_pix_x + *s_mv_x > 0) && (pic4_pix_x + *s_mv_x < max_pos_x4) &&
    (pic4_pix_y + *s_mv_y > 0) && (pic4_pix_y + *s_mv_y < max_pos_y4))
  {
    get_line_p1 = FastLine4X;
  }
  else
  {
    get_line_p1 = UMVLine4X;    
  }
  
  
  //===== loop over search positions =====
  for (best_pos = 0, pos = 1; pos < search_pos4; pos++)
  {
    cand_mv_x = *mv_x + spiral_search_x[pos];    // quarter-pel units
    cand_mv_y = *mv_y + spiral_search_y[pos];    // quarter-pel units
    
    //----- set motion vector cost -----
#ifdef MB32X32_MVC
    if (input->mv_competition > 0)
      mcost = compute_mv_cost_for_each_predictor (lambda_factor, 0, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y, ref, list ^ 1, blocktype, pic_pix_x-img->opix_x, pic_pix_y-img->opix_y);
    else
#endif
      mcost = MV_COST (lambda_factor, 0, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);
    
    if (mcost >= min_mcost) continue;
    cmv_x = cand_mv_x + pic4_pix_x;
    cmv_y = cand_mv_y + pic4_pix_y;
    
    //----- add up SATD -----
    for (y0=0, abort_search=0; y0<blocksize_y && !abort_search; y0+=BLOCK_SIZE)
    {
      ry0 = (y0<<2) + cmv_y;
      ry4  = ry0 + 4;
      ry8  = ry4 + 4;
      ry12 = ry8 + 4;      
      sy0 = (y0<<2) + smv_y;
      sy4  = sy0 + 4;
      sy8  = sy4 + 4;
      sy12 = sy8 + 4;
      y1 = y0 + 1;
      y2 = y1 + 1;
      y3 = y2 + 1;
      
      if (apply_weights)
      {
        
        for (x0=0; x0<blocksize_x; x0+=BLOCK_SIZE)
        {
          rx0 = (x0<<2) + cmv_x;            
          sx0 = (x0<<2) + smv_x;
          d   = diff;
          

#ifdef  USE_HP_FILTER//BROUND
          if(input->UseHPFilter != 0)
          {
          orig_line = &orig_pic [y0][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy0, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry0, rx0,  img_height, img_width);
          weightedpel = clip1a (((( weightSpic * (*(ref_line_p1    )) + weightRpic * (*(ref_line_p2     ))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;          
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line - weightedpel;
          
          orig_line = &orig_pic [y1][x0];
          ref_line_p1 = get_line_p1 (ref1_pic, sy4, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry4, rx0,  img_height, img_width);
          weightedpel = clip1a (((( weightSpic * (*(ref_line_p1    )) + weightRpic * (*(ref_line_p2     ))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;          
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line - weightedpel;
          
          
          orig_line = &orig_pic [y2][x0];
          ref_line_p1 = get_line_p1 (ref1_pic, sy8, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry8, rx0,  img_height, img_width);
          weightedpel = clip1a (((( weightSpic * (*(ref_line_p1    )) + weightRpic * (*(ref_line_p2     ))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;          
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line - weightedpel;
          
          orig_line = &orig_pic [y3][x0];
          ref_line_p1 = get_line_p1 (ref1_pic, sy12, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry12, rx0,  img_height, img_width);
          weightedpel = clip1a (((( weightSpic * (*(ref_line_p1    )) + weightRpic * (*(ref_line_p2     ))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;          
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + (offsetBi<<luma_log_weight_denom)
            + 2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)) );
          *d   = *orig_line - weightedpel;
          }
          else
          {
#endif          
          orig_line = &orig_pic [y0][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy0, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry0, rx0,  img_height, img_width);
          weightedpel = clip1a (((( weightSpic * (*(ref_line_p1    )) + weightRpic * (*(ref_line_p2     ))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;          
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line - weightedpel;
          
          orig_line = &orig_pic [y1][x0];
          ref_line_p1 = get_line_p1 (ref1_pic, sy4, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry4, rx0,  img_height, img_width);
          weightedpel = clip1a (((( weightSpic * (*(ref_line_p1    )) + weightRpic * (*(ref_line_p2     ))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;          
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line - weightedpel;
          
          
          orig_line = &orig_pic [y2][x0];
          ref_line_p1 = get_line_p1 (ref1_pic, sy8, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry8, rx0,  img_height, img_width);
          weightedpel = clip1a (((( weightSpic * (*(ref_line_p1    )) + weightRpic * (*(ref_line_p2     ))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;          
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line - weightedpel;
          
          orig_line = &orig_pic [y3][x0];
          ref_line_p1 = get_line_p1 (ref1_pic, sy12, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry12, rx0,  img_height, img_width);
          weightedpel = clip1a (((( weightSpic * (*(ref_line_p1    )) + weightRpic * (*(ref_line_p2     ))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;          
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d++ = *orig_line++ - weightedpel;
          weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4))) 
            + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          *d   = *orig_line - weightedpel;
#ifdef  USE_HP_FILTER//BROUND
          }
#endif
          
          if (!test8x8transform)
          {
            if ((mcost += SATD (diff, input->hadamard)) > min_mcost)
            {
              abort_search = 1;
              break;
            }
          }
          else
          {
            int x8=x0/8, y8=y0/8;
            i = (x0&0x7) +  x8*64 + (y8 * blocksize_x/8)*64 ;
            for(k=0, j=y0; j<4 + y0; j++, k+=4)
              memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], 4*sizeof(int));
          }
        }
      }
      else
      {      
#ifdef SWITCHED_FILTERS
        img->bipred_rounding_control = (img->nal_reference_idc != 0);
        round = img->currMEOffset[0] + img->currMEOffset[1] + (1 - img->bipred_rounding_control);
#endif
        for (x0=0; x0<blocksize_x; x0+=BLOCK_SIZE)
        {
          rx0 = (x0<<2) + cmv_x;            
          sx0 = (x0<<2) + smv_x;
          d   = diff;
          
#ifdef  USE_HP_FILTER//BROUND
          if(input->UseHPFilter != 0)
          {
          orig_line = &orig_pic [y0][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy0, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry0, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  ((*(ref_line_p1     ) + *(ref_line_p2     ) + round) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          
          orig_line = &orig_pic [y1][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy4, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry4, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  ((*(ref_line_p1     ) + *(ref_line_p2     ) + round) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          
          orig_line = &orig_pic [y2][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy8, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry8, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  ((*(ref_line_p1     ) + *(ref_line_p2     ) + round) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          
          orig_line = &orig_pic [y3][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy12, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry12, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  ((*(ref_line_p1     ) + *(ref_line_p2     ) + round) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d   = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          }
          else
          {
#endif
          orig_line = &orig_pic [y0][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy0, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry0, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  ((*(ref_line_p1     ) + *(ref_line_p2     ) + 1) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          *d++ = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          
          orig_line = &orig_pic [y1][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy4, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry4, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  ((*(ref_line_p1     ) + *(ref_line_p2     ) + 1) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          *d++ = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          
          orig_line = &orig_pic [y2][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy8, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry8, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  ((*(ref_line_p1     ) + *(ref_line_p2     ) + 1) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          *d++ = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          
          orig_line = &orig_pic [y3][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy12, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry12, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  ((*(ref_line_p1     ) + *(ref_line_p2     ) + 1) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          *d++ = *orig_line++  -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
          *d   = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
#ifdef  USE_HP_FILTER//BROUND
          }
#endif
          
          if (!test8x8transform)
          {
            if ((mcost += SATD (diff, input->hadamard)) > min_mcost)
            {
              abort_search = 1;
              break;
            }
          }
          else
          {
            int x8=x0/8, y8=y0/8;
            i = (x0&0x7) +  x8*64 + (y8 * blocksize_x/8)*64 ;
            for(k=0, j=y0; j<4 + y0; j++, k+=4)
              memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], 4*sizeof(int));
          }
        }
      }
    }
    if(test8x8transform)
      mcost += find_SATD32 (c_diff, blocktype, mb_ext_level);
    
    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      best_pos  = pos;
    }
    
  }
  
#ifdef SWITCHED_FILTERS
  if(best_pos || best_pos_sp16)
  {
    if(best_pos_sp16 && (min_mcost_sp16 < min_mcost))
    {
      // Replace 1/4 pixel motion vectors with the result of the offset search
      min_mcost = min_mcost_sp16;
      *mv_x = mv_x_sp16;
      *mv_y = mv_y_sp16;
    }
    else
    {
      *mv_x += spiral_search_x [best_pos];
      *mv_y += spiral_search_y [best_pos];
    }
  }
#else
  if (best_pos)
  {
    *mv_x += spiral_search_x [best_pos];
    *mv_y += spiral_search_y [best_pos];
  }
#endif
  //===== return minimum motion cost =====
  return min_mcost;
}


/*!
***********************************************************************
* \brief
*    Block motion search
***********************************************************************
*/
int                                         //!< minimum motion cost after search
BlockMotionSearch32 (short     ref,           //!< reference idx
                   int       list,          //!< reference pciture list
                   int       mb_x,          //!< x-coordinate inside macroblock   (in the unit of pel)
                   int       mb_y,          //!< y-coordinate inside macroblock   (in the unit of pel)
                   int       blocktype,     //!< block type (1-16x16 ... 7-4x4)
                   int       search_range,  //!< 1-d search range for integer-position search
                   int       lambda_factor, //!< lagrangian parameter for determining motion cost
                   int       mb_ext_level) 
{
  short     pred_mv_x, pred_mv_y, mv_x, mv_y;
  int       i, j;
#ifdef EIGHTH_PEL
  //int       mvshift   = (img->mv_res ? 3 : 2);
#endif
  int       max_value = INT_MAX;
  int       min_mcost = max_value;
  
  int       block_x   = (mb_x>>2);   
  int       block_y   = (mb_y>>2);
  
  int       bsx       = input->blc_size[blocktype][0]<<mb_ext_level; //blc_size = part_sizes * BLOCK_SIZE 
  int       bsy       = input->blc_size[blocktype][1]<<mb_ext_level;

  int       pic_pix_x = img->opix_x + mb_x;
  int       pic_pix_y = img->opix_y + mb_y;
  
  
#ifdef WIN32
  struct _timeb tstruct1;
  struct _timeb tstruct2;
#else
  struct timeb tstruct1;
  struct timeb tstruct2;
#endif

#ifdef SWITCHED_FILTERS
  static int min_mcost_offset;
#endif

  int me_tmp_time;
  short*    pred_mv = img->pred_mv[block_y][block_x][list][ref][blocktype];
  short****** all_mv    = img->all_mv;  
  int list_offset = ((img->MbaffFrameFlag) && (img->mb_data[img->current_mb_nr].mb_field)) ? img->current_mb_nr % 2 ? 4 : 2 : 0;
  int *prevSad = (input->FMEnable == 3)? EPZSDistortion[list + list_offset][blocktype - 1]: NULL;

  static pel_t   orig_val [(MB_BLOCK_SIZE<<MAX_MB_EXT_LEVEL)*(MB_BLOCK_SIZE<<MAX_MB_EXT_LEVEL)];
  static pel_t  *orig_pic  [MB_BLOCK_SIZE<<MAX_MB_EXT_LEVEL];
  for(i=0; i<(MB_BLOCK_SIZE<<mb_ext_level); i++)
  {
    orig_pic[i] = orig_val+(MB_BLOCK_SIZE<<mb_ext_level)*i;
  }

#ifdef WIN32
  _ftime( &tstruct1 );    // start time ms
#else
  ftime(&tstruct1);
#endif
  


  //==================================
  //=====   GET ORIGINAL BLOCK   =====
  //==================================
  for (j = 0; j < bsy; j++)  
    memcpy(orig_pic[j],&imgY_org[pic_pix_y+j][pic_pix_x], bsx *sizeof(imgpel));

#ifndef SIMPLIFY_CODE
  if(input->FMEnable == 1)
  {
    FME_blocktype = blocktype;
    bipred_flag = 0;
  }
  else if (input->FMEnable == 2)
  {
    simplified_setup_FME(ref, list, block_y, block_x, blocktype, all_mv );
  }
#endif

  //===========================================
  //=====   GET MOTION VECTOR PREDICTOR   =====
  //===========================================
#ifdef MB32X32_MVC
  if (input->mv_competition > 0)
    SetMotionVectorPredictor_Competition32 (pred_mv, enc_picture->ref_idx[list], enc_picture->mv[list], ref, list, block_x, block_y, bsx, bsy, mb_ext_level);
  else
#endif
  SetMotionVectorPredictor32 (pred_mv, enc_picture->ref_idx[list], enc_picture->mv[list], ref, list, block_x, block_y, bsx, bsy, mb_ext_level);  
  
#ifndef SIMPLIFY_CODE
  //Dynamic Search Range
  if (input->FMEnable == 1 && input->DynamicSearchRange)
  {
#ifdef _FULL_SEARCH_RANGE_
    if      (input->full_search == 2) search_range = dsr_new_search_range;
    else if (input->full_search == 1) search_range = dsr_new_search_range /  (min(ref,1)+1);
    else                              search_range = dsr_new_search_range / ((min(ref,1)+1) * min(2,blocktype));
#else
    search_range = dsr_new_search_range / ((min(ref,1)+1) * min(2,blocktype));
#endif
  }
#endif

  pred_mv_x = pred_mv[0];
  pred_mv_y = pred_mv[1];
  
  //==================================
  //=====   INTEGER-PEL SEARCH   =====
  //==================================
  
#ifndef SIMPLIFY_CODE
  if (input->FMEnable == 1)
  {
#ifdef EIGHTH_PEL
    mv_x = pred_mv_x / (1 << mvshift);
    mv_y = pred_mv_y / (1 << mvshift);
#else
    mv_x = pred_mv_x / 4;
    mv_y = pred_mv_y / 4;
#endif

    if (!input->rdopt)
    {
      //--- adjust search center so that the (0,0)-vector is inside ---
      mv_x = max (-search_range, min (search_range, mv_x));
      mv_y = max (-search_range, min (search_range, mv_y));
    }
    
    mv_x = Clip3(-2047 + search_range, 2047 - search_range, mv_x);
    mv_y = Clip3(LEVELMVLIMIT[img->LevelIndex][0] + search_range, LEVELMVLIMIT[img->LevelIndex][1]  - search_range, mv_y);
    
    
    min_mcost = FastIntegerPelBlockMotionSearch(orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
      pred_mv_x, pred_mv_y, &mv_x, &mv_y, search_range,
      min_mcost, lambda_factor);
  }
  else if (input->FMEnable == 2)
  {
    mv_x = pred_mv_x / 4;
    mv_y = pred_mv_y / 4;
    
    if (!input->rdopt)
    {
      //--- adjust search center so that the (0,0)-vector is inside ---
      mv_x = max (-search_range, min (search_range, mv_x));
      mv_y = max (-search_range, min (search_range, mv_y));
    }
    
    mv_x = Clip3(-2047 + search_range, 2047 - search_range, mv_x);
    mv_y = Clip3(LEVELMVLIMIT[img->LevelIndex][0] + search_range, LEVELMVLIMIT[img->LevelIndex][1]  - search_range, mv_y);
    
    
    min_mcost = simplified_FastIntegerPelBlockMotionSearch (orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
      pred_mv_x, pred_mv_y, &mv_x, &mv_y, search_range,
      min_mcost, lambda_factor);
    for (i=0; i < (bsx>>2); i++)
    {
      for (j=0; j < (bsy>>2); j++)
      {
        if(list == 0) 
        {
          simplified_fastme_l0_cost[blocktype][(img->pix_y>>2)+block_y+j][(img->pix_x>>2)+block_x+i] = min_mcost;
        }
        else
        {
          simplified_fastme_l1_cost[blocktype][(img->pix_y>>2)+block_y+j][(img->pix_x>>2)+block_x+i] = min_mcost;
        }
      }
    }
  }
  //--- perform motion search using EPZS schemes---
  else if (input->FMEnable == 3)
#endif
  {    
    //--- set search center ---
    mv_x = pred_mv_x / 4;
    mv_y = pred_mv_y / 4;
    //mv_x = (pred_mv_x + 2)>> 2;
    //mv_y = (pred_mv_y + 2)>> 2;    
    if (!input->rdopt)
    {
      //--- adjust search center so that the (0,0)-vector is inside ---
      mv_x = max (-search_range, min (search_range, mv_x));
      mv_y = max (-search_range, min (search_range, mv_y));
    }
    
    
    mv_x = Clip3(-2047 + search_range, 2047 - search_range, mv_x);
    mv_y = Clip3(LEVELMVLIMIT[img->LevelIndex][0] + search_range, LEVELMVLIMIT[img->LevelIndex][1]  - search_range, mv_y);
    

    min_mcost = EPZSPelBlockMotionSearch32 (orig_pic, ref, list, list_offset, 
      enc_picture->ref_idx, enc_picture->mv,pic_pix_x, pic_pix_y, blocktype,
      pred_mv_x, pred_mv_y, &mv_x, &mv_y, search_range, min_mcost, lambda_factor, mb_ext_level); 
#ifdef SWITCHED_FILTERS
    if(img->filterParam == SIFO_FIRST_PASS_FPO)
    {
      if(img->currMEOffset[list] == 0)
      {
        min_mcost_offset = min_mcost;
      }
      else
      {
        if(min_mcost < min_mcost_offset)
        {
          min_mcost_offset = min_mcost;
          if(input->EPZSSpatialMem)  
          {
            EPZSMotion[list + list_offset][ref][blocktype - 1][block_y][pic_pix_x >> 2][0] = (short) mv_x;
            EPZSMotion[list + list_offset][ref][blocktype - 1][block_y][pic_pix_x >> 2][1] = (short) mv_y;
          }
        }
      }
    }
#endif  
    
  }
#ifndef SIMPLIFY_CODE
  else
  {
#ifndef _FAST_FULL_ME_
    
    //--- set search center ---
    mv_x = pred_mv_x / 4;
    mv_y = pred_mv_y / 4;
    if (!input->rdopt)
    {
      //--- adjust search center so that the (0,0)-vector is inside ---
      mv_x = max (-search_range, min (search_range, mv_x));
      mv_y = max (-search_range, min (search_range, mv_y));
    }
    
    mv_x = Clip3(-2047 + search_range, 2047 - search_range, mv_x);
    mv_y = Clip3(LEVELMVLIMIT[img->LevelIndex][0] + search_range, LEVELMVLIMIT[img->LevelIndex][1]  - search_range, mv_y);
    
    //--- perform motion search ---
    min_mcost = FullPelBlockMotionSearch     (orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
      pred_mv_x, pred_mv_y, &mv_x, &mv_y, search_range,
      min_mcost, lambda_factor);
    
#else
    // comments:   - orig_pic is not used  -> be careful
    //             - search center is automatically determined
    min_mcost = FastFullPelBlockMotionSearch (orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
      pred_mv_x, pred_mv_y, &mv_x, &mv_y, search_range,
      min_mcost, lambda_factor);
    
#endif // #ifndef _FAST_FULL_ME_
  }
#endif
  //===== convert search center to quarter-pel units =====
  mv_x <<= 2;
  mv_y <<= 2;
  
  //==============================
  //=====   SUB-PEL SEARCH   =====
  //==============================
  if (!input->DisableSubpelME)
  {
    if (input->FMEnable != 3 || (ref == 0 || img->structure != FRAME || (ref > 0 && min_mcost < 3.5 * prevSad[pic_pix_x >> 2])))
    {
      
      if (input->hadamard == 1)
      {
        min_mcost = max_value;
      }
#ifndef SIMPLIFY_CODE      
      if (input->FMEnable == 1)
      {
        if(blocktype >3)
        {
#ifdef EIGHTH_PEL
          min_mcost =  FastSubPelBlockMotionSearch (orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
            pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9, 9, min_mcost, lambda_factor, (input->Transform8x8Mode && blocktype <= 4 && input->hadamard));
#else
          min_mcost =  FastSubPelBlockMotionSearch (orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
            pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9, min_mcost, lambda_factor, (input->Transform8x8Mode && blocktype <= 4 && input->hadamard));
#endif
        }
        else
        {
#ifdef EIGHTH_PEL
          min_mcost =  SubPelBlockMotionSearch (orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
            pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9, 9, min_mcost, lambda_factor);
#else
          min_mcost =  SubPelBlockMotionSearch (orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
            pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9,  min_mcost, lambda_factor);
#endif
        }
      }
      else if (input->FMEnable == 2)
      {
        if(blocktype > 1)
        {
#ifdef EIGHTH_PEL
          min_mcost =  simplified_FastSubPelBlockMotionSearch (orig_pic, ref, list, pic_pix_x, pic_pix_y,
            blocktype, pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9, 9, min_mcost, lambda_factor,
            (input->Transform8x8Mode && blocktype <= 4 && input->hadamard));
#else
          min_mcost =  simplified_FastSubPelBlockMotionSearch (orig_pic, ref, list, pic_pix_x, pic_pix_y,
            blocktype, pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9,  min_mcost, lambda_factor,
            (input->Transform8x8Mode && blocktype <= 4 && input->hadamard));
#endif
        }
        else
        {
#ifdef EIGHTH_PEL
          min_mcost =  simplified_FastFullSubPelBlockMotionSearch (orig_pic, ref, list, pic_pix_x, pic_pix_y,
            blocktype, pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9, 9, min_mcost, lambda_factor);
#else
          min_mcost =  simplified_FastFullSubPelBlockMotionSearch (orig_pic, ref, list, pic_pix_x, pic_pix_y,
            blocktype, pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9,  min_mcost, lambda_factor);
#endif
        }
      }
      else
#endif
      {
        if (input->EPZSSubPelME)
#ifdef EIGHTH_PEL

          min_mcost =  EPZSSubPelBlockMotionSearch32 (orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
          pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9, 9, min_mcost, lambda_factor, mb_ext_level);            
#else
          min_mcost =  EPZSSubPelBlockMotionSearch32 (orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
          pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9, min_mcost, lambda_factor, mb_ext_level);            
#endif
        else
        {
          error("should not get here !\n", -1);
#ifdef EIGHTH_PEL
          min_mcost =  SubPelBlockMotionSearch (orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
          pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9, 9, min_mcost, lambda_factor);
#else
        min_mcost =  SubPelBlockMotionSearch (orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
          pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9,  min_mcost, lambda_factor);
#endif
        }
      }
    }
  }
  
#ifndef SIMPLIFY_CODE
  if (!input->rdopt)
  {
    // Get the skip mode cost
    if (blocktype == 1 && (img->type == P_SLICE||img->type == SP_SLICE))
    {
      int cost;
      
      FindSkipModeMotionVector ();
      
      cost  = GetSkipCostMB (lambda_factor);
      cost -= ((lambda_factor + 4096) >> 13);
      
      if (cost < min_mcost)
      {
        min_mcost = cost;
        mv_x      = img->all_mv [0][0][0][0][0][0];
        mv_y      = img->all_mv [0][0][0][0][0][1];
      }
    }
  }
#endif
  //===============================================
  //=====   SET MV'S AND RETURN MOTION COST   =====
  //===============================================
  
#ifdef MB32X32_MVC
  if (input->mv_competition > 0)
    GetBestPredVector( pred_mv,mv_x,mv_y,lambda_factor,ref,list,blocktype,mb_x,mb_y);  
#endif
 
#ifdef SWITCHED_FILTERS
  if(img->filterParam == SIFO_FIRST_PASS_FPO)
  {
    if(img->currMEOffset[list] == 0)
    {
      for(j = block_y; j < block_y + (bsy>>2); j++)    
      {
        for(i = block_x; i < block_x + (bsx>>2); i++)
        {
          all_mv[j][i][list][ref][blocktype][0] = mv_x;
          all_mv[j][i][list][ref][blocktype][1] = mv_y;
        }
      }
    }
    else
    {
      if(min_mcost < img->min_mcost_offset[list])
      {
        for(j = block_y; j < block_y + (bsy>>2); j++)    
        {
          for(i = block_x; i < block_x + (bsx>>2); i++)
          {
            all_mv[j][i][list][ref][blocktype][0] = mv_x;
            all_mv[j][i][list][ref][blocktype][1] = mv_y;
          }
        }
      }
    }
  }
  else
  {
    for (j=block_y; j < block_y + (bsy>>2); j++)    
    {
      for (i=block_x; i < block_x + (bsx>>2); i++)
      {
        all_mv[j][i][list][ref][blocktype][0] = mv_x;
        all_mv[j][i][list][ref][blocktype][1] = mv_y;
      }
    }
  }
#else
  for (j=block_y; j < block_y + (bsy>>2); j++)      
  {
    for (i=block_x; i < block_x + (bsx>>2); i++)
    {
      all_mv[j][i][list][ref][blocktype][0] = mv_x;
      all_mv[j][i][list][ref][blocktype][1] = mv_y;
    }
  }
#endif
  
  if (img->type==B_SLICE && input->BiPredMotionEstimation!=0 && (blocktype == 1) && (ref==0))
  {
    
    short   ******bipred_mv = list ? img->bipred_mv1 : img->bipred_mv2;
    int     min_mcostbi = max_value;
    short   bimv_x, bimv_y, tempmv_x ,tempmv_y;
    short   pred_mv_x1, pred_mv_y1;
    short   pred_mv_x2 = 0, pred_mv_y2 = 0;
    short   iterlist=list;
    short   pred_mv_bi[2];
    
#ifdef SWITCHED_FILTERS
    short old_mv_x = mv_x;
    short old_mv_y = mv_y;
    short mv_x_b16, mv_y_b16, bimv_x_b16, bimv_y_b16;
    int min_mcostbi_b16 = max_value;
    int k;
    int other_list = list ^ 1;
#endif

    if(input->FMEnable == 1)
      bipred_flag = 1;

#ifdef MB32X32_MVC
    if (input->mv_competition > 0)
      SetMotionVectorPredictor_Competition32 (pred_mv_bi, enc_picture->ref_idx[list ^ 1], enc_picture->mv[(list == LIST_0? LIST_1: LIST_0)], 0, (list == LIST_0? LIST_1: LIST_0), block_x, block_y, bsx, bsy, mb_ext_level);
    else
#endif
      SetMotionVectorPredictor32             (pred_mv_bi, enc_picture->ref_idx[list ^ 1], enc_picture->mv[(list == LIST_0? LIST_1: LIST_0)], 0, (list == LIST_0? LIST_1: LIST_0), block_x, block_y, bsx, bsy, mb_ext_level);
    
    mv_x=(mv_x + 2)>>2;
    mv_y=(mv_y + 2)>>2;
    
#ifdef SWITCHED_FILTERS
    img->currMEOffset[other_list] = 0;
#endif

    for (i=0;i<=input->BiPredMERefinements;i++)
    {
      if (i%2)
      {
        pred_mv_x2=pred_mv[0];
        pred_mv_y2=pred_mv[1]; 
        pred_mv_x1=pred_mv_bi[0];
        pred_mv_y1=pred_mv_bi[1]; 
        tempmv_x=bimv_x;
        tempmv_y=bimv_y;        
        bimv_x=mv_x;
        bimv_y=mv_y;
        iterlist= list ^ 1;
        
      }
      else
      {
        pred_mv_x1=pred_mv[0];
        pred_mv_y1=pred_mv[1]; 
        pred_mv_x2=pred_mv_bi[0];
        pred_mv_y2=pred_mv_bi[1]; 
        
        if (i!=0)
        {
          tempmv_x=bimv_x;
          tempmv_y=bimv_y;        
          bimv_x=mv_x;
          bimv_y=mv_y;
        }
        else
        {
          tempmv_x=mv_x;
          tempmv_y=mv_y;        
          bimv_x = (pred_mv_x2 + 2)>>2;
          bimv_y = (pred_mv_y2 + 2)>>2;
        }
        
        iterlist=list;
      }
      if (input->FMEnable == 3)        
      {
        min_mcostbi = EPZSBiPredBlockMotionSearch32 (orig_pic, ref, iterlist, list_offset, 
          enc_picture->ref_idx, enc_picture->mv, 
          pic_pix_x, pic_pix_y, blocktype, 
          pred_mv_x1, pred_mv_y1, pred_mv_x2, pred_mv_y2, 
          &bimv_x, &bimv_y, &tempmv_x, &tempmv_y, 
          input->BiPredMESearchRange, min_mcostbi, lambda_factor, mb_ext_level);
      }
#ifndef SIMPLIFY_CODE
      //bipred mode
      else  if(input->FMEnable == 1)
      {
        min_mcostbi = FastBipredIntegerPelBlockMotionSearch (orig_pic, ref, iterlist, 
          pic_pix_x, pic_pix_y, blocktype, 
          pred_mv_x1, pred_mv_y1, pred_mv_x2, pred_mv_y2, 
          &bimv_x, &bimv_y, &tempmv_x, &tempmv_y, 
          input->BiPredMESearchRange>>i, min_mcostbi, lambda_factor);
      }
      else
      {
        min_mcostbi = FullPelBlockMotionBiPred (orig_pic, ref, iterlist, 
          pic_pix_x, pic_pix_y, blocktype, 
          pred_mv_x1, pred_mv_y1, pred_mv_x2, pred_mv_y2, 
          &bimv_x, &bimv_y, &tempmv_x, &tempmv_y, 
          input->BiPredMESearchRange>>i, min_mcostbi, lambda_factor);
      }
#endif
      mv_x=tempmv_x;
      mv_y=tempmv_y;        
    }
    
    mv_x=tempmv_x << 2;
    mv_y=tempmv_y << 2;
    bimv_x = bimv_x << 2;
    bimv_y = bimv_y << 2;
    
    if (input->BiPredMESubPel && !input->DisableSubpelME)
    {
      if (input->hadamard)
      {
        min_mcostbi = max_value;
      }
      
#ifdef EIGHTH_PEL

      min_mcostbi =  SubPelBlockSearchBiPred32 (orig_pic, 0, iterlist, pic_pix_x, pic_pix_y, blocktype,
        pred_mv_x2, pred_mv_y2, &bimv_x, &bimv_y, &mv_x, &mv_y, 9, 9, 9,
        min_mcostbi, lambda_factor, mb_ext_level);
#else
      min_mcostbi =  SubPelBlockSearchBiPred32 (orig_pic, 0, iterlist, pic_pix_x, pic_pix_y, blocktype,
        pred_mv_x2, pred_mv_y2, &bimv_x, &bimv_y, &mv_x, &mv_y, 9, 9,
        min_mcostbi, lambda_factor, mb_ext_level);
#endif
    }
    
    if (input->BiPredMESubPel==2 && !input->DisableSubpelME)
    {
      if (input->hadamard)
      {
        min_mcostbi = max_value;
      }
      
#ifdef EIGHTH_PEL
      min_mcostbi =  SubPelBlockSearchBiPred32 (orig_pic, 0, (iterlist == LIST_0? LIST_1: LIST_0), pic_pix_x, pic_pix_y, blocktype,
        pred_mv_x, pred_mv_y, &mv_x, &mv_y, &bimv_x, &bimv_y, 9, 9, 9,
        min_mcostbi, lambda_factor, mb_ext_level);
#else
      min_mcostbi =  SubPelBlockSearchBiPred32 (orig_pic, 0, (iterlist == LIST_0? LIST_1: LIST_0), pic_pix_x, pic_pix_y, blocktype,
        pred_mv_x, pred_mv_y, &mv_x, &mv_y, &bimv_x, &bimv_y, 9, 9,
        min_mcostbi, lambda_factor, mb_ext_level);
#endif
    }
    
#ifdef SWITCHED_FILTERS
    if(img->filterParam == SIFO_FIRST_PASS_FPO)
    {
      // Store results of the default search
      mv_x_b16 = mv_x;
      mv_y_b16 = mv_y;
      bimv_x_b16 = bimv_x;
      bimv_y_b16 = bimv_y;
      min_mcostbi_b16 = min_mcostbi;

      // Try offsets
      if((ref == 0) && (img->noMEOffsets[other_list] > 0))
      {
        for(k = 0; k < img->noMEOffsets[other_list]; ++k)
        {
          // Restore original values
          mv_x = old_mv_x;
          mv_y = old_mv_y;
          min_mcostbi = max_value;

          img->currMEOffset[other_list] = img->offsetMETab[other_list][k];
          for (i=0; i<=input->BiPredMERefinements; i++)
          {
            if (i%2)
            {
              pred_mv_x2=pred_mv[0];
              pred_mv_y2=pred_mv[1]; 
              pred_mv_x1=pred_mv_bi[0];
              pred_mv_y1=pred_mv_bi[1]; 
              tempmv_x=bimv_x;
              tempmv_y=bimv_y;        
              bimv_x=mv_x;
              bimv_y=mv_y;
              iterlist= list ^ 1;
            }
            else
            {
              pred_mv_x1=pred_mv[0];
              pred_mv_y1=pred_mv[1]; 
              pred_mv_x2=pred_mv_bi[0];
              pred_mv_y2=pred_mv_bi[1]; 

              if (i!=0)
              {
                tempmv_x=bimv_x;
                tempmv_y=bimv_y;        
                bimv_x=mv_x;
                bimv_y=mv_y;
              }
              else
              {
                tempmv_x=mv_x;
                tempmv_y=mv_y;        
                bimv_x = (pred_mv_x2 + 2)>>2;
                bimv_y = (pred_mv_y2 + 2)>>2;
              }

              iterlist=list;
            }
            if (input->FMEnable == 3)        
            {
              min_mcostbi = EPZSBiPredBlockMotionSearch32 (orig_pic, ref, iterlist, list_offset, 
                enc_picture->ref_idx, enc_picture->mv, 
                pic_pix_x, pic_pix_y, blocktype, 
                pred_mv_x1, pred_mv_y1, pred_mv_x2, pred_mv_y2, 
                &bimv_x, &bimv_y, &tempmv_x, &tempmv_y, 
                input->BiPredMESearchRange, min_mcostbi, lambda_factor, mb_ext_level);
            }

            mv_x=tempmv_x;
            mv_y=tempmv_y;        
          }

          mv_x=tempmv_x << 2;
          mv_y=tempmv_y << 2;
          bimv_x = bimv_x << 2;
          bimv_y = bimv_y << 2;

          if (input->BiPredMESubPel && !input->DisableSubpelME)
          {
            if (input->hadamard)
            {
              min_mcostbi = max_value;
            }

            min_mcostbi =  SubPelBlockSearchBiPred32 (orig_pic, 0, iterlist, pic_pix_x, pic_pix_y, blocktype,
              pred_mv_x2, pred_mv_y2, &bimv_x, &bimv_y, &mv_x, &mv_y, 9, 9, 9,
              min_mcostbi, lambda_factor, mb_ext_level);
          }

          if (input->BiPredMESubPel==2 && !input->DisableSubpelME)
          {
            if (input->hadamard)
            {
              min_mcostbi = max_value;
            }
            min_mcostbi =  SubPelBlockSearchBiPred32 (orig_pic, 0, (iterlist == LIST_0? LIST_1: LIST_0), pic_pix_x, pic_pix_y, blocktype,
              pred_mv_x, pred_mv_y, &mv_x, &mv_y, &bimv_x, &bimv_y, 9, 9, 9,
              min_mcostbi, lambda_factor, mb_ext_level);

          }

          // Check cost
          if(min_mcostbi < min_mcostbi_b16)
          {
            min_mcostbi_b16 = min_mcostbi;
            mv_x_b16 = mv_x;
            mv_y_b16 = mv_y;
            bimv_x_b16 = bimv_x;
            bimv_y_b16 = bimv_y;
          }
        }
      }

      // Returns results
      mv_x = mv_x_b16;
      mv_y = mv_y_b16;
      bimv_x = bimv_x_b16;
      bimv_y = bimv_y_b16;
      min_mcostbi = min_mcostbi_b16;
    }
#endif // SWITCHED_FILTERS

    
#ifdef SWITCHED_FILTERS
    if(img->filterParam == SIFO_FIRST_PASS_FPO)
    {
      if(img->currMEOffset[list] == 0)
      {
        img->min_mcost_offset_bi[list] = min_mcostbi;
        for(j=block_y; j < block_y + (bsy>>2); j++)
        {
          for(i=block_x ; i < block_x + (bsx>>2); i++)      
          {
            bipred_mv[j][i][iterlist    ][0][blocktype][0] = mv_x;
            bipred_mv[j][i][iterlist    ][0][blocktype][1] = mv_y;
            bipred_mv[j][i][iterlist ^ 1][0][blocktype][0] = bimv_x;
            bipred_mv[j][i][iterlist ^ 1][0][blocktype][1] = bimv_y;        
          }
        }
      }
      else
      {
        if(min_mcostbi < img->min_mcost_offset_bi[list])
        {
          img->min_mcost_offset_bi[list] = min_mcostbi;
          for (j=block_y; j < block_y + (bsy>>2); j++)
          {
            for (i=block_x ; i < block_x + (bsx>>2); i++)      
            {
              bipred_mv[j][i][iterlist    ][0][blocktype][0] = mv_x;
              bipred_mv[j][i][iterlist    ][0][blocktype][1] = mv_y;
              bipred_mv[j][i][iterlist ^ 1][0][blocktype][0] = bimv_x;
              bipred_mv[j][i][iterlist ^ 1][0][blocktype][1] = bimv_y;        
            }
          }
        }
      }
    }
    else
    {
      for (j=block_y; j < block_y + (bsy>>2); j++)
      {
        for (i=block_x ; i < block_x + (bsx>>2); i++)      
        {
          bipred_mv[j][i][iterlist    ][0][blocktype][0] = mv_x;
          bipred_mv[j][i][iterlist    ][0][blocktype][1] = mv_y;
          bipred_mv[j][i][iterlist ^ 1][0][blocktype][0] = bimv_x;
          bipred_mv[j][i][iterlist ^ 1][0][blocktype][1] = bimv_y;        
        }
      }
    }
#else
    for (j=block_y; j < block_y + (bsy>>2); j++) 
    {
      for (i=block_x ; i < block_x + (bsx>>2); i++)      
      {
        bipred_mv[j][i][iterlist    ][0][blocktype][0] = mv_x;
        bipred_mv[j][i][iterlist    ][0][blocktype][1] = mv_y;
        bipred_mv[j][i][iterlist ^ 1][0][blocktype][0] = bimv_x;
        bipred_mv[j][i][iterlist ^ 1][0][blocktype][1] = bimv_y;        
      }
    }
#endif
  }
  
  
#ifdef WIN32
  _ftime(&tstruct2);   // end time ms
#else
  ftime(&tstruct2);    // end time ms
#endif
  
  me_tmp_time=(int)((tstruct2.time*1000+tstruct2.millitm) - (tstruct1.time*1000+tstruct1.millitm)); 
  me_tot_time += me_tmp_time;
  me_time += me_tmp_time;
  
  return min_mcost;
}


/*!
************************************************************************
* \brief
*    Motion search for a partition
*    applicable to inter 32x32, 32x16 and 16x32 only
*    blocktype 0-32x32, 1-32x32, 2-32x16, 3-16x32
************************************************************************
*/
void
PartitionMotionSearch32 (int    blocktype,
                       int    block8x8, //MB32x32: block8x8=32x32, 32x16, 16x32, 16x16; MB16x16: block8x8=16x16, 16x8, 8x16, 8x8  
                       int    lambda_factor,
#ifdef RDO_Q
                       int   transform8x8,
#endif
                       int  mb_ext_level
                       )
{
  //upper left corner co-or in the unit of 4x4 pel per block 
  static int  bx0_ext[3][5][4] = { {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,2,0,0}, {0,2,0,2}},
                               {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,4,0,0}, {0,4,0,4}},
                               {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,8,0,0}, {0,8,0,8}}
                             };
  static int  by0_ext[3][5][4] = { {{0,0,0,0}, {0,0,0,0}, {0,2,0,0}, {0,0,0,0}, {0,0,2,2}},
                               {{0,0,0,0}, {0,0,0,0}, {0,4,0,0}, {0,0,0,0}, {0,0,4,4}},
                               {{0,0,0,0}, {0,0,0,0}, {0,8,0,0}, {0,0,0,0}, {0,0,8,8}}
                             };
  
  char  **ref_array;
  short ***mv_array;
  short *all_mv;
  short ref;
  int   v, h, mcost, search_range, i, j;
  int   pic_block_x, pic_block_y;
  int   bslice    = (img->type==B_SLICE);
  int   parttype  = (blocktype<4?blocktype:4);

  int   step_h0   = (input->part_size[ parttype][0])<<mb_ext_level; // in the unit of 4x4 pel per block (at or above 8x8 level)
  int   step_v0   = (input->part_size[ parttype][1])<<mb_ext_level; 
  int   step_h    = (input->part_size[blocktype][0])<<mb_ext_level; // all levels
  int   step_v    = (input->part_size[blocktype][1])<<mb_ext_level;
  int   list;
  int   numlists  = bslice ? 2 : 1;
  int   list_offset = img->mb_data[img->current_mb_nr].list_offset; 
  int   block_y;
  int   *m_cost = NULL;

  int   by = by0_ext[mb_ext_level][parttype][block8x8]; //scalable to MB32X32 because by0 and bx0 have been adapted to MB32X32
  int   bx = bx0_ext[mb_ext_level][parttype][block8x8];

#ifdef SWITCHED_FILTERS
  int   mcost_temp;
#endif
  int   mb32i = (mb_ext_level==1) ? curr_mb32i : 4;//mb_ext_level=1: 0~3; mb_ext_level=2: 4

  //===== LOOP OVER REFERENCE FRAMES =====
  for (list=0; list<numlists;list++)
  {
    for (ref=0; ref < listXsize[list+list_offset]; ref++)
    {
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
#ifdef TRELLIS_ONLY
      if(0)
#else
      if(input->UseRDO_Q || img->slice_fractional_quant_flag)
#endif
#else
      if(input->UseRDO_Q)
#endif
      {
        if (Motion_Selected == 0)  
          m_cost = &motion_cost[blocktype][list][ref][block8x8];
        else
        {
          if (transform8x8) motion_cost[blocktype][list][ref][block8x8] = motion32_cost8[mb32i][blocktype][list][ref][block8x8];
          else motion_cost[blocktype][list][ref][block8x8] = motion32_cost4[mb32i][blocktype][list][ref][block8x8];
        }
      }
      else
#endif
        m_cost = &motion_cost[blocktype][list][ref][block8x8];
      
      //----- set search range ---
#ifdef _FULL_SEARCH_RANGE_
      if      (input->full_search == 2) 
        search_range = input->search_range;
      else if (input->full_search == 1) 
        search_range = input->search_range /  (min(ref,1)+1);
      else                              
        search_range = input->search_range / ((min(ref,1)+1) * min(2,blocktype));
#else
      search_range = input->search_range / ((min(ref,1)+1) * min(2,blocktype));
#endif
      
      //----- set arrays -----
      ref_array = enc_picture->ref_idx[list];
      mv_array  = enc_picture->mv[list];
      
      //----- init motion cost -----
      //motion_cost[blocktype][list][ref][block8x8] = 0;
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
#ifdef TRELLIS_ONLY
      if(1)
#else
      if (!(input->UseRDO_Q  || img->slice_fractional_quant_flag) || Motion_Selected == 0)
#endif
#else
      if (!(input->UseRDO_Q) || Motion_Selected == 0)
#endif
#endif
        *m_cost = 0;
      
      //===== LOOP OVER SUB MACRO BLOCK partitions
      //for 32x32, 32x16 and 16x32 step_h(v)0=step_h(v)  
      for (v=by; v<by + step_v0; v += step_v)// by: the coordinator of 8x8 block in the unit of 4x4 in a 16x16 MB 
      {
        pic_block_y = img->block_y + v;
        
        for (h=bx; h<bx+step_h0; h+=step_h)
        {
          all_mv = img->all_mv[v][h][list][ref][blocktype];
          pic_block_x = img->block_x + h;
          
          //--- motion search for block ---
          
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
#ifdef TRELLIS_ONLY
          if(0)
#else
          if(input->UseRDO_Q || img->slice_fractional_quant_flag) 
#endif
#else
          if(input->UseRDO_Q)
#endif
          {
            int       bsx       = (input->blc_size[blocktype][0])<<mb_ext_level; //blc_size = part_sizes * BLOCK_SIZE 
            int       bsy       = (input->blc_size[blocktype][1])<<mb_ext_level;
            if (Motion_Selected == 0)
            { 
#ifdef SWITCHED_FILTERS
              img->currMEOffset[list] = 0;
              mcost = BlockMotionSearch32   (ref, list, h<<2, v<<2, blocktype, search_range, lambda_factor, mb_ext_level); //scalable to MB32X32
              if(img->filterParam == SIFO_FIRST_PASS_FPO)
              {
                img->min_mcost_offset[list] = mcost;
                if (img->noMEOffsets[list] > 0 && (ref == 0))
                {
                  for (i=0; i<img->noMEOffsets[list]; i++)
                  {
                    img->currMEOffset[list] = img->offsetMETab[list][i];
                    mcost_temp = BlockMotionSearch32   (ref, list, h<<2, v<<2, blocktype, search_range, lambda_factor, mb_ext_level); //scalable to MB32X32
                    mcost = (mcost_temp < mcost)? mcost_temp: mcost;
                    img->min_mcost_offset[list] = mcost;
                  }
                }
              }
#else
              mcost = BlockMotionSearch32   (ref, list, h<<2, v<<2, blocktype, search_range, lambda_factor, mb_ext_level); //scalable to MB32X32
#endif
              *m_cost += mcost;
              
              if (transform8x8)
              {
                for (j=v; j < v + (bsy>>2); j++)    
                {
                  for (i=h; i < h + (bsx>>2); i++)
                  {
                    img->all32_mv8[mb32i][j][i][list][ref][blocktype][0] = img->all_mv[j][i][list][ref][blocktype][0];
                    img->all32_mv8[mb32i][j][i][list][ref][blocktype][1] = img->all_mv[j][i][list][ref][blocktype][1];

                    img->bipred32_mv1_8[mb32i][j][i][list][ref][blocktype][0] = img->bipred_mv1[j][i][list][ref][blocktype][0];
                    img->bipred32_mv1_8[mb32i][j][i][list][ref][blocktype][1] = img->bipred_mv1[j][i][list][ref][blocktype][1];

                    img->bipred32_mv2_8[mb32i][j][i][list][ref][blocktype][0] = img->bipred_mv2[j][i][list][ref][blocktype][0];
                    img->bipred32_mv2_8[mb32i][j][i][list][ref][blocktype][1] = img->bipred_mv2[j][i][list][ref][blocktype][1];
                  }
                }
              }
              else
              {
                for (j=v; j < v + (bsy>>2); j++)    
                {
                  for (i=h; i < h + (bsx>>2); i++)
                  {
                    img->all32_mv4[mb32i][j][i][list][ref][blocktype][0] = img->all_mv[j][i][list][ref][blocktype][0];
                    img->all32_mv4[mb32i][j][i][list][ref][blocktype][1] = img->all_mv[j][i][list][ref][blocktype][1];

                    img->bipred32_mv1_4[mb32i][j][i][list][ref][blocktype][0] = img->bipred_mv1[j][i][list][ref][blocktype][0];
                    img->bipred32_mv1_4[mb32i][j][i][list][ref][blocktype][1] = img->bipred_mv1[j][i][list][ref][blocktype][1];

                    img->bipred32_mv2_4[mb32i][j][i][list][ref][blocktype][0] = img->bipred_mv2[j][i][list][ref][blocktype][0];
                    img->bipred32_mv2_4[mb32i][j][i][list][ref][blocktype][1] = img->bipred_mv2[j][i][list][ref][blocktype][1];
                  }
                }
              }
            }
            else
            {
              if (transform8x8)
              {
                for (j=v; j < v + (bsy>>2); j++)    
                {
                  for (i=h; i < h + (bsx>>2); i++)
                 {
                    img->all_mv[j][i][list][ref][blocktype][0] = img->all32_mv8[mb32i][j][i][list][ref][blocktype][0];
                    img->all_mv[j][i][list][ref][blocktype][1] = img->all32_mv8[mb32i][j][i][list][ref][blocktype][1];

                    img->bipred_mv1[j][i][list][ref][blocktype][0] = img->bipred32_mv1_8[mb32i][j][i][list][ref][blocktype][0];
                    img->bipred_mv1[j][i][list][ref][blocktype][1] = img->bipred32_mv1_8[mb32i][j][i][list][ref][blocktype][1];

                    img->bipred_mv2[j][i][list][ref][blocktype][0] = img->bipred32_mv2_8[mb32i][j][i][list][ref][blocktype][0];
                    img->bipred_mv2[j][i][list][ref][blocktype][1] = img->bipred32_mv2_8[mb32i][j][i][list][ref][blocktype][1];
                  }
                }
              }
              else
              {
                for (j=v; j < v + (bsy>>2); j++)    
                {
                  for (i=h; i < h + (bsx>>2); i++)
                  {
                    img->all_mv[j][i][list][ref][blocktype][0] = img->all32_mv4[mb32i][j][i][list][ref][blocktype][0];
                    img->all_mv[j][i][list][ref][blocktype][1] = img->all32_mv4[mb32i][j][i][list][ref][blocktype][1];

                    img->bipred_mv1[j][i][list][ref][blocktype][0] = img->bipred32_mv1_4[mb32i][j][i][list][ref][blocktype][0];
                    img->bipred_mv1[j][i][list][ref][blocktype][1] = img->bipred32_mv1_4[mb32i][j][i][list][ref][blocktype][1];

                    img->bipred_mv2[j][i][list][ref][blocktype][0] = img->bipred32_mv2_4[mb32i][j][i][list][ref][blocktype][0];
                    img->bipred_mv2[j][i][list][ref][blocktype][1] = img->bipred32_mv2_4[mb32i][j][i][list][ref][blocktype][1];
                  }
                }
              }
            }
          }
          else
          {
#ifdef SWITCHED_FILTERS
            img->currMEOffset[list] = 0;
            mcost = BlockMotionSearch32     (ref, list, h<<2, v<<2, blocktype, search_range, lambda_factor, mb_ext_level); //scalable to MB32X32
            if(img->filterParam == SIFO_FIRST_PASS_FPO)
            {
              img->min_mcost_offset[list] = mcost;
              if(img->noMEOffsets[list] > 0 && ref == 0)
              {
                for(i=0; i<img->noMEOffsets[list]; i++)
                {
                  img->currMEOffset[list] = img->offsetMETab[list][i];
                  mcost_temp = BlockMotionSearch32     (ref, list, h<<2, v<<2, blocktype, search_range, lambda_factor, mb_ext_level); //scalable to MB32X32
                  mcost = (mcost_temp < mcost)? mcost_temp: mcost;
                  img->min_mcost_offset[list] = mcost;
                }
              }    
            }
#else
            mcost = BlockMotionSearch32     (ref, list, h<<2, v<<2, blocktype, search_range, lambda_factor, mb_ext_level); //scalable to MB32X32
#endif

            *m_cost += mcost;
          }
          
#else
           mcost = BlockMotionSearch32     (ref, list, h<<2, v<<2, blocktype, search_range, lambda_factor, mb_ext_level); 
          *m_cost += mcost;
#endif
          //--- set motion vectors and reference frame (for motion vector prediction) ---
          for (j=0; j<step_v; j++) 
          {
            block_y = pic_block_y+j;  
            memset(&ref_array [block_y][pic_block_x], ref, step_h * sizeof(char));
            
            for (i=0; i<step_h; i++)
            {
              memcpy(mv_array  [block_y][pic_block_x+i], all_mv, 2* sizeof(short));
              //mv_array  [block_y][pic_block_x+i][0] = all_mv[0];
              //mv_array  [block_y][pic_block_x+i][1] = all_mv[1];
            }
          }
        }
      }
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
#ifdef TRELLIS_ONLY         
      if(0)
#else
      if(input->UseRDO_Q || img->slice_fractional_quant_flag)

#endif
#else
      if(input->UseRDO_Q)
#endif
      {
        if (Motion_Selected == 0)
        {
          if (transform8x8)

            motion32_cost8[mb32i][blocktype][list][ref][block8x8]= motion_cost[blocktype][list][ref][block8x8];
          else
            motion32_cost4[mb32i][blocktype][list][ref][block8x8]= motion_cost[blocktype][list][ref][block8x8];
          
          
        }
        else
        {
          if (transform8x8)
            motion_cost[blocktype][list][ref][block8x8]= motion32_cost8[mb32i][blocktype][list][ref][block8x8];
          else
            motion_cost[blocktype][list][ref][block8x8]= motion32_cost4[mb32i][blocktype][list][ref][block8x8];
        }
      }
#endif
    }
  }
}


/*!
***********************************************************************
* \brief
*    Motion Cost for Bidirectional modes
***********************************************************************
*/
int BIDPartitionCost32 (int   blocktype,
                      int   block8x8,
                      short fw_ref,
                      short bw_ref,
                      int   lambda_factor, int mb_ext_level)
{
  static int  bx0_ext[3][5][4] = { {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,2,0,0}, {0,2,0,2}},
                               {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,4,0,0}, {0,4,0,4}},
                               {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,8,0,0}, {0,8,0,8}}
                             };
  static int  by0_ext[3][5][4] = { {{0,0,0,0}, {0,0,0,0}, {0,2,0,0}, {0,0,0,0}, {0,0,2,2}},
                               {{0,0,0,0}, {0,0,0,0}, {0,4,0,0}, {0,0,0,0}, {0,0,4,4}},
                               {{0,0,0,0}, {0,0,0,0}, {0,8,0,0}, {0,0,0,0}, {0,0,8,8}}
                             };
  int   diff[64];
  int   curr_blk[MB_BLOCK_SIZE<<MAX_MB_EXT_LEVEL][MB_BLOCK_SIZE<<MAX_MB_EXT_LEVEL]; // ABT pred.error buffer
  int   bsx       = min(input->blc_size[blocktype][0],8); 
  int   bsy       = min(input->blc_size[blocktype][1],8);
  
  int   pic_pix_x, pic_pix_y, block_x, block_y;
  int   v, h, mcost, i, j, k;
  int   mvd_bits  = 0;
  int   parttype  = (blocktype<4?blocktype:4);
  int   step_h0   = (input->part_size[ parttype][0])<<mb_ext_level;
  int   step_v0   = (input->part_size[ parttype][1])<<mb_ext_level;
  int   step_h    = (input->part_size[blocktype][0])<<mb_ext_level;
  int   step_v    = (input->part_size[blocktype][1])<<mb_ext_level;
  int   bxx, byy;                               // indexing curr_blk
  int   bx = bx0_ext[mb_ext_level][parttype][block8x8];
  int   by = by0_ext[mb_ext_level][parttype][block8x8];
  short   ******all_mv = img->all_mv;
  short   ******  p_mv = img->pred_mv;
  
  //----- cost for motion vector bits -----
  for (v=by; v<by + step_v0; v+=step_v)  
  {
    for (h=bx; h<bx + step_h0; h+=step_h)
    {
      mvd_bits += mvbits[ all_mv [v][h][LIST_0][fw_ref][blocktype][0] - p_mv[v][h][LIST_0][fw_ref][blocktype][0] ];
      mvd_bits += mvbits[ all_mv [v][h][LIST_0][fw_ref][blocktype][1] - p_mv[v][h][LIST_0][fw_ref][blocktype][1] ];
      
      mvd_bits += mvbits[ all_mv [v][h][LIST_1][bw_ref][blocktype][0] - p_mv[v][h][LIST_1][bw_ref][blocktype][0] ];
      mvd_bits += mvbits[ all_mv [v][h][LIST_1][bw_ref][blocktype][1] - p_mv[v][h][LIST_1][bw_ref][blocktype][1] ];
    }
  }
  
  mcost = WEIGHTED_COST (lambda_factor, mvd_bits);
  
  //----- cost of residual signal -----
  for (byy=0, v=by; v<by + step_v0; byy+=4, v++)  
  {
    pic_pix_y = img->opix_y + (block_y = (v<<2));  
    for (bxx=0, h=bx; h<bx + step_h0; bxx+=4, h++)
    {
      pic_pix_x = img->opix_x + (block_x = (h<<2));
      LumaPrediction4x4 (block_x, block_y, 2, blocktype, blocktype, fw_ref, bw_ref);
      
      for (k=j=0; j<4; j++)
      {
        for (  i=0; i<4; i++)
          diff[k++] = curr_blk[byy+j][bxx+i] = 
          imgY_org[pic_pix_y+j][pic_pix_x+i] - img->mpr[j+block_y][i+block_x];
      }
      if ((!input->Transform8x8Mode) || (blocktype>4))
        mcost += SATD (diff, input->hadamard);
    }
  }
  if (input->Transform8x8Mode && (blocktype<=4))  // tchen 4-29-04
  {
    for (byy=0; byy < (input->blc_size[parttype][1]<<mb_ext_level); byy+=bsy)
      for (bxx=0; bxx< (input->blc_size[parttype][0]<<mb_ext_level); bxx+=bsx)
      {        
        for (k=0, j=byy;j<byy + 8;j++, k += 8)
          memcpy(&diff[k], &(curr_blk[j][bxx]), 8 * sizeof(int));          
        
        mcost += SATD8X8(diff, input->hadamard);
      }
  }
  return mcost;
}

/*!
***********************************************************************
* \brief
*    Motion Cost for Bidirectional modes
***********************************************************************
*/
int BPredPartitionCost32 (int   blocktype,
                        int   block8x8,
                        short fw_ref,
                        short bw_ref,
                        int   lambda_factor,
                        int   list,
                        int mb_ext_level)
{
  static int  bx0_ext[3][5][4] = { {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,2,0,0}, {0,2,0,2}},
                                   {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,4,0,0}, {0,4,0,4}},
                                   {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,8,0,0}, {0,8,0,8}}
                                 };
  static int  by0_ext[3][5][4] = { {{0,0,0,0}, {0,0,0,0}, {0,2,0,0}, {0,0,0,0}, {0,0,2,2}},
                                   {{0,0,0,0}, {0,0,0,0}, {0,4,0,0}, {0,0,0,0}, {0,0,4,4}},
                                   {{0,0,0,0}, {0,0,0,0}, {0,8,0,0}, {0,0,0,0}, {0,0,8,8}}
                                 };
  int (*bx0)[4]=bx0_ext[mb_ext_level];
  int (*by0)[4]=by0_ext[mb_ext_level];
  int   diff[64];
  int   curr_blk[MB_BLOCK_SIZE<<MAX_MB_EXT_LEVEL][MB_BLOCK_SIZE<<MAX_MB_EXT_LEVEL]; // ABT pred.error buffer
  int   bsx       = min(input->blc_size[blocktype][0],8); 
  int   bsy       = min(input->blc_size[blocktype][1],8);
  
  int   pic_pix_x, pic_pix_y, block_x, block_y;
  int   v, h, mcost, i, j, k;
  int   mvd_bits  = 0;
  int   parttype  = (blocktype<4?blocktype:4);
  int   step_h0   = (input->part_size[ parttype][0])<<mb_ext_level;
  int   step_v0   = (input->part_size[ parttype][1])<<mb_ext_level;
  int   step_h    = (input->part_size[blocktype][0])<<mb_ext_level;
  int   step_v    = (input->part_size[blocktype][1])<<mb_ext_level;
  int   bxx, byy;                               // indexing curr_blk
  
  short   ******all_mv = list ? img->bipred_mv1 : img->bipred_mv2;
  short   ******  p_mv = img->pred_mv;
  
  for (v=by0[parttype][block8x8]; v<by0[parttype][block8x8]+step_v0; v+=step_v)
  {
    for (h=bx0[parttype][block8x8]; h<bx0[parttype][block8x8]+step_h0; h+=step_h)
    {
      mvd_bits += mvbits[ all_mv [v][h][LIST_0][fw_ref][blocktype][0] - p_mv[v][h][LIST_0][fw_ref][blocktype][0] ];
      mvd_bits += mvbits[ all_mv [v][h][LIST_0][fw_ref][blocktype][1] - p_mv[v][h][LIST_0][fw_ref][blocktype][1] ];
      
      mvd_bits += mvbits[ all_mv [v][h][LIST_1][bw_ref][blocktype][0] - p_mv[v][h][LIST_1][bw_ref][blocktype][0] ];
      mvd_bits += mvbits[ all_mv [v][h][LIST_1][bw_ref][blocktype][1] - p_mv[v][h][LIST_1][bw_ref][blocktype][1] ];
    }
  }
  mcost = WEIGHTED_COST (lambda_factor, mvd_bits);
  
  //----- cost of residual signal -----
  for (byy=0, v=by0[parttype][block8x8]; v<by0[parttype][block8x8]+step_v0; byy+=4, v++) 
  {
    
    pic_pix_y = img->opix_y + (block_y = (v<<2));      

    for (bxx=0, h=bx0[parttype][block8x8]; h<bx0[parttype][block8x8]+step_h0; bxx+=4, h++)
    {
      pic_pix_x = img->opix_x + (block_x = (h<<2));        
      LumaPrediction4x4Bi (block_x, block_y, 2, blocktype, blocktype, fw_ref, bw_ref, list);
      
      for (k=j=0; j<4; j++)
      {
        for (  i=0; i<4; i++)
          diff[k++] = curr_blk[byy+j][bxx+i] = 
          imgY_org[pic_pix_y+j][pic_pix_x+i] - img->mpr[j+block_y][i+block_x];            
      }

#ifndef SIMPLIFY_CODE
      if ((!input->Transform8x8Mode) || (blocktype>4))
      {
        mcost += SATD (diff, input->hadamard);          
      }
#endif
    }
  }

  if (input->Transform8x8Mode && (blocktype<=4))  // tchen 4-29-04
  {
    for (byy=0; byy < (input->blc_size[parttype][1]<<mb_ext_level); byy+=bsy)
      for (bxx=0; bxx< (input->blc_size[parttype][0]<<mb_ext_level); bxx+=bsx)
      {
        for (k=0, j=byy;j<byy + 8;j++, k += 8)
          memcpy(&diff[k], &(curr_blk[j][bxx]), 8 * sizeof(int));          
        
        mcost += SATD8X8(diff, input->hadamard);
      }
  }
  return mcost;
}


/*!
************************************************************************
* \brief
*    Find motion vector for the Skip mode
************************************************************************
*/
void FindSkipModeMotionVector32 (int mb_ext_level)
{
  int   bx, by;
  short ******all_mv = img->all_mv;
  
  short pmv[2];
  
  int zeroMotionAbove;
  int zeroMotionLeft;
  PixelPos mb_a, mb_b;
  int      a_mv_y = 0;
  int      a_ref_idx = 0;
  int      b_mv_y = 0;
  int      b_ref_idx = 0;
  
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  
  
  
  // Set all the possible values for the skipped motion vector (in the case of the skip mode) according
  // to each possible predictor
  
#ifdef MB32X32_MVC
  int prediction_mode_for_skip;
  
  for (prediction_mode_for_skip = 0; prediction_mode_for_skip<mv_comp.nb_mode_for_skip; prediction_mode_for_skip++)
  {
    if (mv_comp.predictor_for_skip[prediction_mode_for_skip] == PRED_H264_MEDIAN)  
    {
#endif
      getLuma4x4Neighbour(img->current_mb_nr,0,0,-1, 0,&mb_a);
      getLuma4x4Neighbour(img->current_mb_nr,0,0, 0,-1,&mb_b);
      
      if (mb_a.available)
      {
        a_mv_y    = enc_picture->mv[LIST_0][mb_a.pos_y][mb_a.pos_x][1];
        a_ref_idx = enc_picture->ref_idx[LIST_0][mb_a.pos_y][mb_a.pos_x];
        
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
        b_mv_y    = enc_picture->mv[LIST_0][mb_b.pos_y][mb_b.pos_x][1];
        b_ref_idx = enc_picture->ref_idx[LIST_0][mb_b.pos_y][mb_b.pos_x];
        
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
      
      zeroMotionLeft  = !mb_a.available ? 1 : a_ref_idx==0 && enc_picture->mv[LIST_0][mb_a.pos_y][mb_a.pos_x][0]==0 && a_mv_y==0 ? 1 : 0;
      zeroMotionAbove = !mb_b.available ? 1 : b_ref_idx==0 && enc_picture->mv[LIST_0][mb_b.pos_y][mb_b.pos_x][0]==0 && b_mv_y==0 ? 1 : 0;
      
      if (zeroMotionAbove || zeroMotionLeft)
      {
        for (by = 0;by < (4<<mb_ext_level);by++)
          for (bx = 0;bx < (4<<mb_ext_level);bx++)
          {
            memset(all_mv [by][bx][0][0][0], 0, 2* sizeof(short));
          }
          
#ifdef MB32X32_MVC
          if (input->mv_competition > 0)
          {
            mv_comp.mv_pred_skip[prediction_mode_for_skip][0] = 0;
            mv_comp.mv_pred_skip[prediction_mode_for_skip][1] = 0;
          }
#endif
      }
      else
      {
        SetMotionVectorPredictor32 (pmv, enc_picture->ref_idx[LIST_0], enc_picture->mv[LIST_0], 0, LIST_0, 0, 0, (16<<mb_ext_level), (16<<mb_ext_level), mb_ext_level);
        
        for (by = 0;by < (4<<mb_ext_level);by++)
          for (bx = 0;bx < (4<<mb_ext_level);bx++)
          {
            memcpy(all_mv [by][bx][0][0][0], pmv, 2* sizeof(short));
            //all_mv [by][bx][0][0][0][0] = pmv[0];
            //all_mv [by][bx][0][0][0][1] = pmv[1];
          }
#ifdef MB32X32_MVC
          if (input->mv_competition > 0)
          {
            mv_comp.mv_pred_skip[prediction_mode_for_skip][0] = pmv[0];
            mv_comp.mv_pred_skip[prediction_mode_for_skip][1] = pmv[1];
          }
#endif
      }
#ifdef MB32X32_MVC
    }
    else   
    {
      PixelPos block_a, block_b, block_c, block_d;
      int block_x = 0, block_y = 0;
      int blockshape_x = 16 << mb_ext_level;  //rahul  (<< mb_ext_level)//, blockshape_y = 16;
      int blockshape_y = 16 << mb_ext_level;  //rahul  (<< mb_ext_level)//, blockshape_y = 16;
      int hv, mv_a, mv_b, mv_c;
      int y=(int)img->current_mb_nr/(img->width/16);  // Vertical
      int x=(int)img->current_mb_nr%(img->width/16);  // Horizontal
      //      int ref_frame = 0;
      int mb_x                 = 4*block_x;
      int mb_y                 = 4*block_y;      
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
      if (mv_comp.predictor_for_skip[prediction_mode_for_skip] == PRED_A)
      {
        for (hv=0; hv < 2; hv++)
        {
          if (block_a.available)  
          {
            mv_a = enc_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][hv];
            pmv[hv] = mv_a;
          }
          else 
            pmv[hv] = NOT_AVAILABLE;
        }
        
        mv_comp.mv_pred_skip[prediction_mode_for_skip][0] = pmv[0];
        mv_comp.mv_pred_skip[prediction_mode_for_skip][1] = pmv[1];
        
        if (mv_comp.nb_mode_for_skip == 1)
        {
          for (by = 0;by < 4<<mb_ext_level;by++)
            for (bx = 0;bx < 4<<mb_ext_level;bx++)
              memcpy(all_mv [by][bx][0][0][0], pmv, 2* sizeof(short));
        }
      }
      else if (mv_comp.predictor_for_skip[prediction_mode_for_skip] == PRED_B)
      {
        for (hv=0; hv < 2; hv++)
        {
          if (block_b.available)  
          {
            mv_b = enc_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][hv];
            pmv[hv] = mv_b;
            
          }
          else 
            pmv[hv] = NOT_AVAILABLE;
        }
        
        mv_comp.mv_pred_skip[prediction_mode_for_skip][0] = pmv[0];
        mv_comp.mv_pred_skip[prediction_mode_for_skip][1] = pmv[1];
        
        if (mv_comp.nb_mode_for_skip == 1)
        {
          for (by = 0;by < 4<<mb_ext_level;by++)
            for (bx = 0;bx < 4<<mb_ext_level;bx++)
              memcpy(all_mv [by][bx][0][0][0], pmv, 2* sizeof(short));
        }
      }
      else if (mv_comp.predictor_for_skip[prediction_mode_for_skip] == PRED_C)
      {
        for (hv=0; hv < 2; hv++)
        {
          if (block_c.available)  
          {
            mv_c = enc_picture->mv[LIST_0][block_c.pos_y][block_c.pos_x][hv];
            pmv[hv] = mv_c;
            
          }
          else 
            pmv[hv] = NOT_AVAILABLE;
        }
        
        mv_comp.mv_pred_skip[prediction_mode_for_skip][0] = pmv[0];
        mv_comp.mv_pred_skip[prediction_mode_for_skip][1] = pmv[1];
        
        if (mv_comp.nb_mode_for_skip == 1)
        {
          for (by = 0;by < 4<<mb_ext_level;by++)
            for (bx = 0;bx < 4<<mb_ext_level;bx++)
              memcpy(all_mv [by][bx][0][0][0], pmv, 2* sizeof(short));
        }
      }
      else if (mv_comp.predictor_for_skip[prediction_mode_for_skip] == PRED_ZERO)
      {
        pmv[0] = 0;
        pmv[1] = 0;
        
        mv_comp.mv_pred_skip[prediction_mode_for_skip][0] = pmv[0];
        mv_comp.mv_pred_skip[prediction_mode_for_skip][1] = pmv[1];
        
        if (mv_comp.nb_mode_for_skip == 1)
        {
          for (by = 0;by < 4<<mb_ext_level;by++)
            for (bx = 0;bx < 4<<mb_ext_level;bx++)
              memcpy(all_mv [by][bx][0][0][0], pmv, 2* sizeof(short));
        }
      }
      else if (mv_comp.predictor_for_skip[prediction_mode_for_skip] == PRED_COLOCATED)
      {
        if (collocated_mv_available(y, x, LIST_0) == TRUE)
        {
          // The collocated is the top left MV of the macroblock, whatever the subpartition...
          // To be improved.
          pmv[0] = listX[0][0]->mv[LIST_0][y*4][x*4][0] / (listX[0][0]->ref_idx[0][y*4][x*4] + 1);
          pmv[1] = listX[0][0]->mv[LIST_0][y*4][x*4][1] / (listX[0][0]->ref_idx[0][y*4][x*4] + 1);
        }
        else
        {
          pmv[0] = NOT_AVAILABLE;
          pmv[1] = NOT_AVAILABLE;
        }
        
        mv_comp.mv_pred_skip[prediction_mode_for_skip][0] = pmv[0];
        mv_comp.mv_pred_skip[prediction_mode_for_skip][1] = pmv[1];
        
        if (mv_comp.nb_mode_for_skip == 1)
        {
          for (by = 0;by < 4<<mb_ext_level;by++)
            for (bx = 0;bx < 4<<mb_ext_level;bx++)
              memcpy(all_mv [by][bx][0][0][0], pmv, 2* sizeof(short));
        }
      }
      else if (mv_comp.predictor_for_skip[prediction_mode_for_skip] == PRED_EXTENDEDSPATIAL)
      {
        for (hv=0; hv < 2; hv++)
        {
          if ((block_a.available) && (block_b.available) && (block_c.available))  
          {
            mv_a = enc_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][hv];
            mv_b = enc_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][hv];
            mv_c = enc_picture->mv[LIST_0][block_c.pos_y][block_c.pos_x][hv];  
            pmv[hv] = mv_a+mv_b+mv_c-min(mv_a,min(mv_b,mv_c))-max(mv_a,max(mv_b,mv_c));
          }
          else if (block_a.available) 
          {
            mv_a = enc_picture->mv[LIST_0][block_a.pos_y][block_a.pos_x][hv];
            pmv[hv] = mv_a;
          }
          else if (block_b.available)  
          {
            mv_b = enc_picture->mv[LIST_0][block_b.pos_y][block_b.pos_x][hv];    
            pmv[hv] = mv_b;
          }
          else if (block_c.available)  
          {
            mv_c = enc_picture->mv[LIST_0][block_c.pos_y][block_c.pos_x][hv];
            pmv[hv] = mv_c;
          }
          else 
            pmv[hv] = 0;
          
          mv_comp.mv_pred_skip[prediction_mode_for_skip][0] = pmv[0];
          mv_comp.mv_pred_skip[prediction_mode_for_skip][1] = pmv[1];
          
          if (mv_comp.nb_mode_for_skip == 1)
          {
            for (by = 0;by < 4<<mb_ext_level;by++)
              for (bx = 0;bx < 4<<mb_ext_level;bx++)
                memcpy(all_mv [by][bx][0][0][0], pmv, 2 * sizeof(short));
          }
        }    
      }
    }
  }
#endif
}

/*!
************************************************************************
* \brief
*    Calculate Direct Motion Vectors  *****
************************************************************************
*/
void Get_Direct_Motion_Vectors32 (int mb_ext_level)
{
  
  int   block_x, block_y, pic_block_x, pic_block_y, opic_block_x, opic_block_y;
  short ****all_mvs;
#ifndef SIMPLIFY_CODE
  int   mv_scale;
  int refList; 
  int ref_idx;   
#endif
  byte  **   moving_block;
  short ****   co_located_mv;
  char  ***    co_located_ref_idx;
  int64 ***    co_located_ref_id;
  char  **     ref_pic_l0 = enc_picture->ref_idx[LIST_0];
  char  **     ref_pic_l1 = enc_picture->ref_idx[LIST_1];
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  
  if (currMB->list_offset)
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
  
  if (img->direct_spatial_mv_pred_flag)  //spatial direct mode copy from decoder
  {
    
    short l0_refA=0, l0_refB=0, l0_refD, l0_refC=0;
    short l1_refA=0, l1_refB=0, l1_refD, l1_refC=0; 
    short l0_refX,l1_refX;
    short pmvfw[2]={0,0},pmvbw[2]={0,0};
    
    PixelPos mb_a, mb_b, mb_d, mb_c;              
    
    CheckAvailabilityOfNeighbors32(mb_ext_level);

    getLuma4x4Neighbour(img->current_mb_nr, 0, 0, -1,  0,&mb_a);
    getLuma4x4Neighbour(img->current_mb_nr, 0, 0,  0, -1,&mb_b);
    getLuma4x4Neighbour(img->current_mb_nr, 0, 0, (16<<mb_ext_level), -1,&mb_c);
    getLuma4x4Neighbour(img->current_mb_nr, 0, 0, -1, -1,&mb_d);
    
    if (!img->MbaffFrameFlag)
    {
      l0_refA = mb_a.available ? ref_pic_l0[mb_a.pos_y][mb_a.pos_x] : -1;
      l0_refB = mb_b.available ? ref_pic_l0[mb_b.pos_y][mb_b.pos_x] : -1;
      l0_refD = mb_d.available ? ref_pic_l0[mb_d.pos_y][mb_d.pos_x] : -1;
      l0_refC = mb_c.available ? ref_pic_l0[mb_c.pos_y][mb_c.pos_x] : l0_refD;      
      
      l1_refA = mb_a.available ? ref_pic_l1[mb_a.pos_y][mb_a.pos_x] : -1;
      l1_refB = mb_b.available ? ref_pic_l1[mb_b.pos_y][mb_b.pos_x] : -1;
      l1_refD = mb_d.available ? ref_pic_l1[mb_d.pos_y][mb_d.pos_x] : -1;
      l1_refC = mb_c.available ? ref_pic_l1[mb_c.pos_y][mb_c.pos_x] : l1_refD;      
    }
#ifndef SIMPLIFY_CODE
    else
    {
      if (currMB->mb_field)
      {
        l0_refA = mb_a.available 
          ? (img->mb_data[mb_a.mb_addr].mb_field  || ref_pic_l0[mb_a.pos_y][mb_a.pos_x] < 0
          ?  ref_pic_l0[mb_a.pos_y][mb_a.pos_x]
          :  ref_pic_l0[mb_a.pos_y][mb_a.pos_x] * 2) : -1;
        
        l0_refB = mb_b.available 
          ? (img->mb_data[mb_b.mb_addr].mb_field || ref_pic_l0[mb_b.pos_y][mb_b.pos_x] < 0
          ?  ref_pic_l0[mb_b.pos_y][mb_b.pos_x] 
          :  ref_pic_l0[mb_b.pos_y][mb_b.pos_x] * 2) : -1;
        
        l0_refD = mb_d.available 
          ? (img->mb_data[mb_d.mb_addr].mb_field || ref_pic_l0[mb_d.pos_y][mb_d.pos_x] < 0
          ?  ref_pic_l0[mb_d.pos_y][mb_d.pos_x] 
          :  ref_pic_l0[mb_d.pos_y][mb_d.pos_x] * 2) : -1;
        
        l0_refC = mb_c.available 
          ? (img->mb_data[mb_c.mb_addr].mb_field || ref_pic_l0[mb_c.pos_y][mb_c.pos_x] < 0
          ?  ref_pic_l0[mb_c.pos_y][mb_c.pos_x] 
          :  ref_pic_l0[mb_c.pos_y][mb_c.pos_x] * 2) : l0_refD;
        
        l1_refA = mb_a.available 
          ? (img->mb_data[mb_a.mb_addr].mb_field || ref_pic_l1[mb_a.pos_y][mb_a.pos_x] < 0
          ?  ref_pic_l1[mb_a.pos_y][mb_a.pos_x] 
          :  ref_pic_l1[mb_a.pos_y][mb_a.pos_x] * 2) : -1;
        
        l1_refB = mb_b.available 
          ? (img->mb_data[mb_b.mb_addr].mb_field || ref_pic_l1[mb_b.pos_y][mb_b.pos_x] < 0
          ?  ref_pic_l1[mb_b.pos_y][mb_b.pos_x] 
          :  ref_pic_l1[mb_b.pos_y][mb_b.pos_x] * 2) : -1;
        
        l1_refD = mb_d.available 
          ? (img->mb_data[mb_d.mb_addr].mb_field || ref_pic_l1[mb_d.pos_y][mb_d.pos_x] < 0
          ?  ref_pic_l1[mb_d.pos_y][mb_d.pos_x] 
          :  ref_pic_l1[mb_d.pos_y][mb_d.pos_x] * 2) : -1;
        
        l1_refC = mb_c.available 
          ? (img->mb_data[mb_c.mb_addr].mb_field || ref_pic_l1[mb_c.pos_y][mb_c.pos_x] < 0
          ?  ref_pic_l1[mb_c.pos_y][mb_c.pos_x] 
          :  ref_pic_l1[mb_c.pos_y][mb_c.pos_x] * 2) : l1_refD;
      }
      else
      {
        l0_refA = mb_a.available 
          ? (img->mb_data[mb_a.mb_addr].mb_field || ref_pic_l0[mb_a.pos_y][mb_a.pos_x]  < 0 
          ?  ref_pic_l0[mb_a.pos_y][mb_a.pos_x] >> 1 
          :  ref_pic_l0[mb_a.pos_y][mb_a.pos_x]) : -1;
        
        l0_refB = mb_b.available 
          ? (img->mb_data[mb_b.mb_addr].mb_field || ref_pic_l0[mb_b.pos_y][mb_b.pos_x] < 0 
          ?  ref_pic_l0[mb_b.pos_y][mb_b.pos_x] >> 1 
          :  ref_pic_l0[mb_b.pos_y][mb_b.pos_x]) : -1;
        
        l0_refD = mb_d.available 
          ? (img->mb_data[mb_d.mb_addr].mb_field || ref_pic_l0[mb_d.pos_y][mb_d.pos_x] < 0 
          ?  ref_pic_l0[mb_d.pos_y][mb_d.pos_x] >> 1 
          :  ref_pic_l0[mb_d.pos_y][mb_d.pos_x]) : -1;      
        
        l0_refC = mb_c.available 
          ? (img->mb_data[mb_c.mb_addr].mb_field || ref_pic_l0[mb_c.pos_y][mb_c.pos_x] < 0 
          ?  ref_pic_l0[mb_c.pos_y][mb_c.pos_x] >> 1 
          :  ref_pic_l0[mb_c.pos_y][mb_c.pos_x]) : l0_refD;      
        
        l1_refA = mb_a.available 
          ? (img->mb_data[mb_a.mb_addr].mb_field || ref_pic_l1[mb_a.pos_y][mb_a.pos_x] < 0 
          ?  ref_pic_l1[mb_a.pos_y][mb_a.pos_x] >> 1 
          :  ref_pic_l1[mb_a.pos_y][mb_a.pos_x]) : -1;
        
        l1_refB = mb_b.available 
          ? (img->mb_data[mb_b.mb_addr].mb_field || ref_pic_l1[mb_b.pos_y][mb_b.pos_x] < 0 
          ?  ref_pic_l1[mb_b.pos_y][mb_b.pos_x] >> 1 
          :  ref_pic_l1[mb_b.pos_y][mb_b.pos_x]) : -1;
        
        l1_refD = mb_d.available 
          ? (img->mb_data[mb_d.mb_addr].mb_field || ref_pic_l1[mb_d.pos_y][mb_d.pos_x] < 0 
          ?  ref_pic_l1[mb_d.pos_y][mb_d.pos_x] >> 1 
          :  ref_pic_l1[mb_d.pos_y][mb_d.pos_x]) : -1;
        
        l1_refC = mb_c.available 
          ? (img->mb_data[mb_c.mb_addr].mb_field || ref_pic_l1[mb_c.pos_y][mb_c.pos_x] < 0 
          ?  ref_pic_l1[mb_c.pos_y][mb_c.pos_x] >> 1
          :  ref_pic_l1[mb_c.pos_y][mb_c.pos_x]) : l1_refD;
      }
    }
#endif

    l0_refX = (l0_refA >= 0 && l0_refB >= 0) ? min(l0_refA,l0_refB): max(l0_refA,l0_refB);
    l0_refX = (l0_refX >= 0 && l0_refC >= 0) ? min(l0_refX,l0_refC): max(l0_refX,l0_refC);
    
    l1_refX = (l1_refA >= 0 && l1_refB >= 0) ? min(l1_refA,l1_refB): max(l1_refA,l1_refB);
    l1_refX = (l1_refX >= 0 && l1_refC >= 0) ? min(l1_refX,l1_refC): max(l1_refX,l1_refC);        
    
    if (l0_refX >=0)
      SetMotionVectorPredictor32 (pmvfw, enc_picture->ref_idx[LIST_0], enc_picture->mv[LIST_0], l0_refX, LIST_0, 0, 0, (16<<mb_ext_level), (16<<mb_ext_level), mb_ext_level);
    if (l1_refX >=0)
      SetMotionVectorPredictor32 (pmvbw, enc_picture->ref_idx[LIST_1], enc_picture->mv[LIST_1], l1_refX, LIST_1, 0, 0, (16<<mb_ext_level), (16<<mb_ext_level), mb_ext_level);

    for (block_y=0; block_y<(4<<mb_ext_level); block_y++)
    {
      pic_block_y  = (img->pix_y  >> 2) + block_y;
      opic_block_y = (img->opix_y >> 2) + block_y;
      
      for (block_x=0; block_x<(4<<mb_ext_level); block_x++)
      {
        pic_block_x  = (img->pix_x  >> 2) + block_x;
        opic_block_x = (img->opix_x >> 2) + block_x;
        
        all_mvs = img->all_mv[block_y][block_x];
        
        if (l0_refX >=0)
        {
          if (!l0_refX  && !moving_block[opic_block_y][opic_block_x])
          {
            
            memset(all_mvs[LIST_0][0][0], 0, 2* sizeof(short));
            //all_mvs[LIST_0][0][0][0] = 0;
            //all_mvs[LIST_0][0][0][1] = 0;            
            direct_ref_idx[LIST_0][pic_block_y][pic_block_x]=0;       
          }
          else
          {
            all_mvs[LIST_0][l0_refX][0][0] = pmvfw[0];
            all_mvs[LIST_0][l0_refX][0][1] = pmvfw[1];
            direct_ref_idx[LIST_0][pic_block_y][pic_block_x]= (char)l0_refX;              
          }
        }
        else
        {
          all_mvs[LIST_0][0][0][0] = 0;
          all_mvs[LIST_0][0][0][1] = 0;
          direct_ref_idx[LIST_0][pic_block_y][pic_block_x]=-1;          
        }
        
        if (l1_refX >=0)
        {
          if(l1_refX==0 && !moving_block[opic_block_y][opic_block_x])
          {                  
            all_mvs[LIST_1][0][0][0] = 0;
            all_mvs[LIST_1][0][0][1] = 0;
            direct_ref_idx[LIST_1][pic_block_y][pic_block_x]= (char)l1_refX;     
          }
          else
          {
            all_mvs[LIST_1][l1_refX][0][0] = pmvbw[0];
            all_mvs[LIST_1][l1_refX][0][1] = pmvbw[1];
            direct_ref_idx[LIST_1][pic_block_y][pic_block_x]= (char)l1_refX;
          }               
        }
        else
        {      
          direct_ref_idx[LIST_1][pic_block_y][pic_block_x]=-1;
          
          all_mvs[LIST_1][0][0][0] = 0;
          all_mvs[LIST_1][0][0][1] = 0;
        }
#ifndef SIMPLIFY_CODE        
        // Test Level Limits if satisfied.
        if (img->MbaffFrameFlag 
          && (all_mvs[LIST_0][l0_refX < 0? 0 : l0_refX][0][0] < -8192 
          ||  all_mvs[LIST_0][l0_refX < 0? 0 : l0_refX][0][0] >  8191 
          ||  all_mvs[LIST_0][l0_refX < 0? 0 : l0_refX][0][1] < LEVELMVLIMIT[img->LevelIndex][4] 
          ||  all_mvs[LIST_0][l0_refX < 0? 0 : l0_refX][0][1] > LEVELMVLIMIT[img->LevelIndex][5] 
          ||  all_mvs[LIST_1][l1_refX < 0? 0 : l1_refX][0][0] < -8192 
          ||  all_mvs[LIST_1][l1_refX < 0? 0 : l1_refX][0][0] > 8191 
          ||  all_mvs[LIST_1][l1_refX < 0? 0 : l1_refX][0][1] < LEVELMVLIMIT[img->LevelIndex][4] 
          ||  all_mvs[LIST_1][l1_refX < 0? 0 : l1_refX][0][1] > LEVELMVLIMIT[img->LevelIndex][5])) 
        { 
          direct_ref_idx[LIST_0][pic_block_y][pic_block_x] = -1; 
          direct_ref_idx[LIST_1][pic_block_y][pic_block_x] = -1; 
          direct_pdir           [pic_block_y][pic_block_x] = -1; 
        } 
        else 
#endif
        { 
          if (l0_refX < 0 && l1_refX < 0)
          {
            direct_ref_idx[LIST_0][pic_block_y][pic_block_x] = 
              direct_ref_idx[LIST_1][pic_block_y][pic_block_x] = 0;
          }
          if      (direct_ref_idx[LIST_1][pic_block_y][pic_block_x] == -1) 
            direct_pdir[pic_block_y][pic_block_x] = 0;
          else if (direct_ref_idx[LIST_0][pic_block_y][pic_block_x] == -1) 
            direct_pdir[pic_block_y][pic_block_x] = 1;
          else                                                           
            direct_pdir[pic_block_y][pic_block_x] = 2;
        }        
      }
    }
  }
#ifndef SIMPLIFY_CODE   
  else
  {
    int64 *refpic = enc_picture->ref_pic_num[LIST_0 +currMB->list_offset];
    
    //temporal direct mode copy from decoder
    for (block_y = 0; block_y < 4; block_y++)
    {
      pic_block_y  = (img->pix_y  >> 2) + block_y;
      opic_block_y = (img->opix_y >> 2) + block_y;
      
      for (block_x = 0; block_x < 4; block_x++)
      {
        pic_block_x  = (img->pix_x>>2) + block_x;
        opic_block_x = (img->opix_x>>2) + block_x;
        all_mvs = img->all_mv[block_y][block_x];
        
        refList = (co_located_ref_idx[LIST_0][opic_block_y][opic_block_x]== -1 ? LIST_1 : LIST_0);
        ref_idx = co_located_ref_idx[refList][opic_block_y][opic_block_x];
        
        // next P is intra mode
        if (ref_idx==-1)
        {
          all_mvs[LIST_0][0][0][0] = 0;
          all_mvs[LIST_0][0][0][1] = 0;
          all_mvs[LIST_1][0][0][0] = 0;
          all_mvs[LIST_1][0][0][1] = 0;
          direct_ref_idx[LIST_0][pic_block_y][pic_block_x] = 0;
          direct_ref_idx[LIST_1][pic_block_y][pic_block_x] = 0;
          direct_pdir[pic_block_y][pic_block_x] = 2;
        }
        // next P is skip or inter mode
        else 
        {
          int mapped_idx=INVALIDINDEX;
          int iref; 
          
          for (iref=0;iref<min(img->num_ref_idx_l0_active,listXsize[LIST_0+currMB->list_offset]);iref++)
          {
            if (refpic[iref]==co_located_ref_id[refList ][opic_block_y][opic_block_x])
            {
              mapped_idx=iref;
              break;
            }
            else //! invalid index. Default to zero even though this case should not happen
            {                        
              mapped_idx=INVALIDINDEX;
            }
          }
          
          if (mapped_idx !=INVALIDINDEX)
          {
            mv_scale = img->mvscale[LIST_0+currMB->list_offset][mapped_idx];
            
            if (mv_scale==9999)
            {
              // forward
              all_mvs[LIST_0][0][0][0] = co_located_mv[refList][opic_block_y][opic_block_x][0];
              all_mvs[LIST_0][0][0][1] = co_located_mv[refList][opic_block_y][opic_block_x][1];
              // backward
              all_mvs[LIST_1][0][0][0] = 0;
              all_mvs[LIST_1][0][0][1] = 0;
            }
            else
            {
              // forward
              all_mvs[LIST_0][mapped_idx][0][0] = (mv_scale * co_located_mv[refList][opic_block_y][opic_block_x][0] + 128) >> 8;
              all_mvs[LIST_0][mapped_idx][0][1] = (mv_scale * co_located_mv[refList][opic_block_y][opic_block_x][1] + 128) >> 8;
              // backward
              all_mvs[LIST_1][         0][0][0] = ((mv_scale - 256)* co_located_mv[refList][opic_block_y][opic_block_x][0] + 128) >> 8;
              all_mvs[LIST_1][         0][0][1] = ((mv_scale - 256)* co_located_mv[refList][opic_block_y][opic_block_x][1] + 128) >> 8;
            }
            
            // Test Level Limits if satisfied.
            if ( all_mvs[LIST_0][mapped_idx][0][0] < -8192 
              || all_mvs[LIST_0][mapped_idx][0][0] >  8191 
              || all_mvs[LIST_0][mapped_idx][0][1] < LEVELMVLIMIT[img->LevelIndex][4] 
              || all_mvs[LIST_0][mapped_idx][0][1] > LEVELMVLIMIT[img->LevelIndex][5] 
              || all_mvs[LIST_1][0][0][0] < -8192 
              || all_mvs[LIST_1][0][0][0] > 8191 
              || all_mvs[LIST_1][0][0][1] < LEVELMVLIMIT[img->LevelIndex][4] 
              || all_mvs[LIST_1][0][0][1] > LEVELMVLIMIT[img->LevelIndex][5]) 
            { 
              direct_ref_idx[LIST_0][pic_block_y][pic_block_x] = -1; 
              direct_ref_idx[LIST_1][pic_block_y][pic_block_x] = -1; 
              direct_pdir[pic_block_y][pic_block_x] = -1; 
            } 
            else 
            { 
              direct_ref_idx[LIST_0][pic_block_y][pic_block_x] = mapped_idx; 
              direct_ref_idx[LIST_1][pic_block_y][pic_block_x] = 0; 
              direct_pdir[pic_block_y][pic_block_x] = 2; 
            }
          }
          else
          {
            direct_ref_idx[LIST_0][pic_block_y][pic_block_x] = -1;
            direct_ref_idx[LIST_1][pic_block_y][pic_block_x] = -1;
            direct_pdir[pic_block_y][pic_block_x] = -1;
          }
        }
      }
    }
  }
#endif
}
#endif
