
/*!
*************************************************************************************
* \file mv-search.c
*
* \brief
*    Motion Vector Search, unified for B and P Pictures
*
* \author
*    Main contributors (see contributors.h for copyright, address and affiliation details)
*      - Stephan Wenger                  <stewe@cs.tu-berlin.de>
*      - Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
*      - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
*      - Stephan Wenger                  <stewe@cs.tu-berlin.de>
*      - Jani Lainema                    <jani.lainema@nokia.com>
*      - Detlev Marpe                    <marpe@hhi.de>
*      - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
*      - Heiko Schwarz                   <hschwarz@hhi.de>
*      - Alexis Michael Tourapis         <alexismt@ieee.org>
*
*************************************************************************************
*/

#include "contributors.h"

#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#include "global.h"

#include "image.h"
#include "mv-search.h"
#include "refbuf.h"
#include "memalloc.h"
#include "mb_access.h"
#include "fast_me.h"
#include "simplified_fast_me.h"

#include "epzs.h"

#ifdef ADAPTIVE_FILTER
#include "adaptive_filter.h"
#endif

#ifdef MV_COMPETITION
#include "mv_competition.h"
extern MV_Competition mv_comp;
#endif

#include <time.h>
#include <sys/timeb.h>

#ifdef SWITCHED_FILTERS
#include "switched_filters.h"
#endif

// These procedure pointers are used by motion_search() and one_eigthpel()
static pel_t *(*PelYline_11) (pel_t *, int, int, int, int);
static pel_t *(*get_line) (pel_t**, int, int, int, int);
static pel_t *(*get_line_p1) (pel_t**, int, int, int, int);
static pel_t *(*get_line_p2) (pel_t**, int, int, int, int);

// Statistics, temporary
int     max_mvd;
short*  spiral_search_x;
short*  spiral_search_y;
short*  spiral_hpel_search_x;
short*  spiral_hpel_search_y;

int*    mvbits;
int*    refbits;
unsigned int*  byte_abs;
int**** motion_cost;
int     byte_abs_range;
static const short block_type_shift_factor[8] = {0, 0, 1, 1, 2, 3, 3, 1}; // last one relaxed to 1 instead 4

#ifdef RDO_Q
#ifdef MB32X32
extern int curr_mbi;
int motion16_cost4[1<<(2*MAX_MB_EXT_LEVEL)][8][2][2][4];
int motion16_cost8[1<<(2*MAX_MB_EXT_LEVEL)][8][2][2][4];
int motion32_cost4[5][8][2][2][4];
int motion32_cost8[5][8][2][2][4];
#endif
int motion_cost4[8][2][2][4];
int motion_cost8[8][2][2][4];
#endif

#ifdef SWITCHED_FILTERS
extern short ******EPZSMotion;  //!< Array for storing Motion Vectors
#endif

void SetMotionVectorPredictor (short  pmv[2],
                               char   **refPic,
                               short  ***tmp_mv,
                               short  ref_frame,
                               int    list,
                               int    block_x,
                               int    block_y,
                               int    blockshape_x,
                               int    blockshape_y);

extern ColocatedParams *Co_located;

#ifdef _FAST_FULL_ME_

/*****
*****  static variables for fast integer motion estimation
*****
*/
static int  **search_setup_done;  //!< flag if all block SAD's have been calculated yet
static int  **search_center_x;    //!< absolute search center for fast full motion search
static int  **search_center_y;    //!< absolute search center for fast full motion search
static int  **pos_00;             //!< position of (0,0) vector
static distpel  *****BlockSAD;        //!< SAD for all blocksize, ref. frames and motion vectors
static int  **max_search_range;


/*!
***********************************************************************
* \brief
*    function creating arrays for fast integer motion estimation
***********************************************************************
*/
void
InitializeFastFullIntegerSearch ()
{
  int  i, j, k, list;
  int  search_range = input->search_range;
  int  max_pos      = (2*search_range+1) * (2*search_range+1);
  
  if ((BlockSAD = (distpel*****)malloc (2 * sizeof(distpel****))) == NULL)
    no_mem_exit ("InitializeFastFullIntegerSearch: BlockSAD");
  
  for (list=0; list<2;list++)
  {
    if ((BlockSAD[list] = (distpel****)malloc ((img->max_num_references) * sizeof(distpel***))) == NULL)
      no_mem_exit ("InitializeFastFullIntegerSearch: BlockSAD");
    for (i = 0; i < img->max_num_references; i++)
    {
      if ((BlockSAD[list][i] = (distpel***)malloc (8 * sizeof(distpel**))) == NULL)
        no_mem_exit ("InitializeFastFullIntegerSearch: BlockSAD");
      for (j = 1; j < 8; j++)
      {
        if ((BlockSAD[list][i][j] = (distpel**)malloc (16 * sizeof(distpel*))) == NULL)
          no_mem_exit ("InitializeFastFullIntegerSearch: BlockSAD");
        for (k = 0; k < 16; k++)
        {
          if ((BlockSAD[list][i][j][k] = (distpel*)malloc (max_pos * sizeof(distpel))) == NULL)
            no_mem_exit ("InitializeFastFullIntegerSearch: BlockSAD");
        }        
      }
    }
  }
  
  if ((search_setup_done = (int**)malloc (2*sizeof(int*)))==NULL)
    no_mem_exit ("InitializeFastFullIntegerSearch: search_setup_done");
  if ((search_center_x = (int**)malloc (2*sizeof(int*)))==NULL)
    no_mem_exit ("InitializeFastFullIntegerSearch: search_center_x");
  if ((search_center_y = (int**)malloc (2*sizeof(int*)))==NULL)
    no_mem_exit ("InitializeFastFullIntegerSearch: search_center_y");
  if ((pos_00 = (int**)malloc (2*sizeof(int*)))==NULL)
    no_mem_exit ("InitializeFastFullIntegerSearch: pos_00");
  if ((max_search_range = (int**)malloc (2*sizeof(int*)))==NULL)
    no_mem_exit ("InitializeFastFullIntegerSearch: max_search_range");
  
  for (list=0; list<2; list++)
  {
    if ((search_setup_done[list] = (int*)malloc ((img->max_num_references)*sizeof(int)))==NULL)
      no_mem_exit ("InitializeFastFullIntegerSearch: search_setup_done");
    if ((search_center_x[list] = (int*)malloc ((img->max_num_references)*sizeof(int)))==NULL)
      no_mem_exit ("InitializeFastFullIntegerSearch: search_center_x");
    if ((search_center_y[list] = (int*)malloc ((img->max_num_references)*sizeof(int)))==NULL)
      no_mem_exit ("InitializeFastFullIntegerSearch: search_center_y");
    if ((pos_00[list] = (int*)malloc ((img->max_num_references)*sizeof(int)))==NULL)
      no_mem_exit ("InitializeFastFullIntegerSearch: pos_00");
    if ((max_search_range[list] = (int*)malloc ((img->max_num_references)*sizeof(int)))==NULL)
      no_mem_exit ("InitializeFastFullIntegerSearch: max_search_range");
  }
  
  // assign max search ranges for reference frames
  if (input->full_search == 2)
  {
    for (list=0;list<2;list++)
      for (i=0; i<img->max_num_references; i++)  
        max_search_range[list][i] = search_range;
  }
  else
  {
    for (list=0;list<2;list++)
    {
      max_search_range[list][0] = search_range;
      for (i=1; i< img->max_num_references; i++)  max_search_range[list][i] = search_range / 2;
    }
  }
  
}



/*!
***********************************************************************
* \brief
*    function for deleting the arrays for fast integer motion estimation
***********************************************************************
*/
void
ClearFastFullIntegerSearch ()
{
  int  i, j, k, list;
  
  for (list=0; list<2; list++)
  {
    for (i = 0; i < img->max_num_references; i++)
    {
      for (j = 1; j < 8; j++)
      {
        for (k = 0; k < 16; k++)
        {
          free (BlockSAD[list][i][j][k]);
        }
        free (BlockSAD[list][i][j]);
      }
      free (BlockSAD[list][i]);
    }
    free (BlockSAD[list]);
  }
  free (BlockSAD);
  
  for (list=0; list<2; list++)
  {
    free (search_setup_done[list]);
    free (search_center_x[list]);
    free (search_center_y[list]);
    free (pos_00[list]);
    free (max_search_range[list]);
  }
  free (search_setup_done);
  free (search_center_x);
  free (search_center_y);
  free (pos_00);
  free (max_search_range);
  
}


/*!
***********************************************************************
* \brief
*    function resetting flags for fast integer motion estimation
*    (have to be called in start_macroblock())
***********************************************************************
*/
void
ResetFastFullIntegerSearch ()
{
  int i,list;
  
  for (list=0; list<2; list++)
    for (i = 0; i < img->max_num_references; i++)
      search_setup_done [list][i] = 0;
}

/*!
***********************************************************************
* \brief
*    calculation of SAD for larger blocks on the basis of 4x4 blocks
***********************************************************************
*/
void
SetupLargerBlocks (int list, int refindex, int max_pos)
{
#define ADD_UP_BLOCKS()   _o=*_bo; _i=*_bi; _j=*_bj; for(pos=0;pos<max_pos;pos++) _o[pos] = _i[pos] + _j[pos];
#define INCREMENT(inc)    _bo+=inc; _bi+=inc; _bj+=inc;
  
  distpel   pos, **_bo, **_bi, **_bj;
  register distpel *_o,   *_i,   *_j;
  
  //--- blocktype 6 ---
  _bo = BlockSAD[list][refindex][6];
  _bi = BlockSAD[list][refindex][7];
  _bj = _bi + 4;
  ADD_UP_BLOCKS(); INCREMENT(1);
  ADD_UP_BLOCKS(); INCREMENT(1);
  ADD_UP_BLOCKS(); INCREMENT(1);
  ADD_UP_BLOCKS(); INCREMENT(5);
  ADD_UP_BLOCKS(); INCREMENT(1);
  ADD_UP_BLOCKS(); INCREMENT(1);
  ADD_UP_BLOCKS(); INCREMENT(1);
  ADD_UP_BLOCKS();
  
  //--- blocktype 5 ---
  _bo = BlockSAD[list][refindex][5];
  _bi = BlockSAD[list][refindex][7];
  _bj = _bi + 1;
  ADD_UP_BLOCKS(); INCREMENT(2);
  ADD_UP_BLOCKS(); INCREMENT(2);
  ADD_UP_BLOCKS(); INCREMENT(2);
  ADD_UP_BLOCKS(); INCREMENT(2);
  ADD_UP_BLOCKS(); INCREMENT(2);
  ADD_UP_BLOCKS(); INCREMENT(2);
  ADD_UP_BLOCKS(); INCREMENT(2);
  ADD_UP_BLOCKS();
  
  //--- blocktype 4 ---
  _bo = BlockSAD[list][refindex][4];
  _bi = BlockSAD[list][refindex][6];
  _bj = _bi + 1;
  ADD_UP_BLOCKS(); INCREMENT(2);
  ADD_UP_BLOCKS(); INCREMENT(6);
  ADD_UP_BLOCKS(); INCREMENT(2);
  ADD_UP_BLOCKS();
  
  //--- blocktype 3 ---
  _bo = BlockSAD[list][refindex][3];
  _bi = BlockSAD[list][refindex][4];
  _bj = _bi + 8;
  ADD_UP_BLOCKS(); INCREMENT(2);
  ADD_UP_BLOCKS();
  
  //--- blocktype 2 ---
  _bo = BlockSAD[list][refindex][2];
  _bi = BlockSAD[list][refindex][4];
  _bj = _bi + 2;
  ADD_UP_BLOCKS(); INCREMENT(8);
  ADD_UP_BLOCKS();
  
  //--- blocktype 1 ---
  _bo = BlockSAD[list][refindex][1];
  _bi = BlockSAD[list][refindex][3];
  _bj = _bi + 2;
  ADD_UP_BLOCKS();
}


/*!
***********************************************************************
* \brief
*    Setup the fast search for an macroblock
***********************************************************************
*/
void SetupFastFullPelSearch (short ref, int list)  // <--  reference frame parameter, list0 or 1
{
  short   pmv[2];
  pel_t   orig_blocks[256], *orgptr=orig_blocks, *refptr;
  int     offset_x, offset_y, x, y, range_partly_outside, ref_x, ref_y, pos, abs_x, abs_y, bindex, blky;
  int     LineSadBlk0, LineSadBlk1, LineSadBlk2, LineSadBlk3;
  int     max_width, max_height;
  int     img_width, img_height;
#ifdef EIGHTH_PEL
  int     mvshift       = (img->mv_res ? 3 : 2);
#endif
  StorablePicture *ref_picture;
  pel_t   *ref_pic;
  
  distpel**   block_sad = BlockSAD[list][ref][7];
  int     search_range  = max_search_range[list][ref];
  int     max_pos       = (2*search_range+1) * (2*search_range+1);
  
  int     list_offset   = img->mb_data[img->current_mb_nr].list_offset; 
  
  int     apply_weights = ( (active_pps->weighted_pred_flag && (img->type == P_SLICE || img->type == SP_SLICE)) ||
    (active_pps->weighted_bipred_idc && (img->type == B_SLICE)));
  
  
  ref_picture     = listX[list+list_offset][ref];
  
  //===== Use weighted Reference for ME ====
  if (apply_weights && input->UseWeightedReferenceME)
    ref_pic       = ref_picture->imgY_11_w;
  else
    ref_pic       = ref_picture->imgY_11;
  
  max_width     = ref_picture->size_x - 17;
  max_height    = ref_picture->size_y - 17;
  
  img_width     = ref_picture->size_x;
  img_height    = ref_picture->size_y;
  
  //===== get search center: predictor of 16x16 block =====
  SetMotionVectorPredictor (pmv, enc_picture->ref_idx[list], enc_picture->mv[list], ref, list, 0, 0, 16, 16);
  
#ifdef EIGHTH_PEL
  search_center_x[list][ref] = pmv[0] / (1 << mvshift);
  search_center_y[list][ref] = pmv[1] / (1 << mvshift);
#else
  search_center_x[list][ref] = pmv[0] / 4;
  search_center_y[list][ref] = pmv[1] / 4;
#endif
  
  if (!input->rdopt)
  {
    //--- correct center so that (0,0) vector is inside ---
    search_center_x[list][ref] = max(-search_range, min(search_range, search_center_x[list][ref]));
    search_center_y[list][ref] = max(-search_range, min(search_range, search_center_y[list][ref]));
  }
  search_center_x[list][ref] = Clip3(-2047 + search_range, 2047 - search_range, search_center_x[list][ref]);
  search_center_y[list][ref] = Clip3(LEVELMVLIMIT[img->LevelIndex][0] + search_range, LEVELMVLIMIT[img->LevelIndex][1]  - search_range, search_center_y[list][ref]);
  
  search_center_x[list][ref] += img->opix_x;
  search_center_y[list][ref] += img->opix_y;
  
  offset_x = search_center_x[list][ref];
  offset_y = search_center_y[list][ref];
  
  //===== copy original block for fast access =====
  for   (y = img->opix_y; y < img->opix_y+16; y++)
    for (x = img->opix_x; x < img->opix_x+16; x++)
      *orgptr++ = imgY_org [y][x];
    
    
    //===== check if whole search range is inside image =====
    if (offset_x >= search_range && offset_x <= max_width  - search_range &&
      offset_y >= search_range && offset_y <= max_height - search_range   )
    {
      range_partly_outside = 0; PelYline_11 = FastLine16Y_11;
    }
    else
    {
      range_partly_outside = 1;
    }
    
    //===== determine position of (0,0)-vector =====
    if (!input->rdopt)
    {
      ref_x = img->opix_x - offset_x;
      ref_y = img->opix_y - offset_y;
      
      for (pos = 0; pos < max_pos; pos++)
      {
        if (ref_x == spiral_search_x[pos] &&
          ref_y == spiral_search_y[pos])
        {
          pos_00[list][ref] = pos;
          break;
        }
      }
    }
    
    //===== loop over search range (spiral search): get blockwise SAD =====
    for (pos = 0; pos < max_pos; pos++)
    {
      abs_y = offset_y + spiral_search_y[pos];
      abs_x = offset_x + spiral_search_x[pos];
      
      if (range_partly_outside)
      {
        if (abs_y >= 0 && abs_y <= max_height &&
          abs_x >= 0 && abs_x <= max_width    )
        {
          PelYline_11 = FastLine16Y_11;
        }
        else
        {
          PelYline_11 = UMVLine16Y_11;
        }
      }
      
      orgptr = orig_blocks;
      bindex = 0;
      for (blky = 0; blky < 4; blky++)
      {
        LineSadBlk0 = LineSadBlk1 = LineSadBlk2 = LineSadBlk3 = 0;
        for (y = 0; y < 4; y++)
        {
          refptr = PelYline_11 (ref_pic, abs_y++, abs_x, img_height, img_width);
          
          LineSadBlk0 += byte_abs [*refptr++ - *orgptr++];
          LineSadBlk0 += byte_abs [*refptr++ - *orgptr++];
          LineSadBlk0 += byte_abs [*refptr++ - *orgptr++];
          LineSadBlk0 += byte_abs [*refptr++ - *orgptr++];
          LineSadBlk1 += byte_abs [*refptr++ - *orgptr++];
          LineSadBlk1 += byte_abs [*refptr++ - *orgptr++];
          LineSadBlk1 += byte_abs [*refptr++ - *orgptr++];
          LineSadBlk1 += byte_abs [*refptr++ - *orgptr++];
          LineSadBlk2 += byte_abs [*refptr++ - *orgptr++];
          LineSadBlk2 += byte_abs [*refptr++ - *orgptr++];
          LineSadBlk2 += byte_abs [*refptr++ - *orgptr++];
          LineSadBlk2 += byte_abs [*refptr++ - *orgptr++];
          LineSadBlk3 += byte_abs [*refptr++ - *orgptr++];
          LineSadBlk3 += byte_abs [*refptr++ - *orgptr++];
          LineSadBlk3 += byte_abs [*refptr++ - *orgptr++];
          LineSadBlk3 += byte_abs [*refptr++ - *orgptr++];
        }
        block_sad[bindex++][pos] = LineSadBlk0;
        block_sad[bindex++][pos] = LineSadBlk1;
        block_sad[bindex++][pos] = LineSadBlk2;
        block_sad[bindex++][pos] = LineSadBlk3;
      }
    }
    
    
    //===== combine SAD's for larger block types =====
    SetupLargerBlocks (list, ref, max_pos);
    
    //===== set flag marking that search setup have been done =====
    search_setup_done[list][ref] = 1;
}
#endif // _FAST_FULL_ME_

/*!
************************************************************************
* \brief
*    Set motion vector predictor
************************************************************************
*/
void SetMotionVectorPredictor (short  pmv[2],
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
  int mb_nr                = img->current_mb_nr;
  
  int mv_a, mv_b, mv_c, pred_vec=0;
  int mvPredType, rFrameL, rFrameU, rFrameUR;
  int hv;
  
  PixelPos block_a, block_b, block_c, block_d;
  
  // added for bipred mode
  int *** fastme_l0_cost_flag = (input->FMEnable == 1) ? (bipred_flag ? fastme_l0_cost_bipred:fastme_l0_cost) : NULL;
  int *** fastme_l1_cost_flag = (input->FMEnable == 1) ? (bipred_flag ? fastme_l1_cost_bipred:fastme_l1_cost) : NULL;    
  
  //Dynamic Search Range
  
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
      }
      else
      {
        if (mb_x+blockshape_x == 8)  block_c.available = 0;
      }
    }
    else
    {
      if (mb_x+blockshape_x == 16)   block_c.available = 0;
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
  
  /* Prediction if only one of the neighbors uses the reference frame
  * we are checking
  */
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
  
  // neighborhood SAD prediction
  if(input->FMEnable == 1 &&(input->DynamicSearchRange == 1 || input->BiPredMotionEstimation == 1))
  {
    SAD_a = block_a.available ? ((list==1) ? (fastme_l1_cost_flag[FME_blocktype][block_a.pos_y][block_a.pos_x]) : (fastme_l0_cost_flag[FME_blocktype][block_a.pos_y][block_a.pos_x])) : 0;
    SAD_b = block_b.available ? ((list==1) ? (fastme_l1_cost_flag[FME_blocktype][block_b.pos_y][block_b.pos_x]) : (fastme_l0_cost_flag[FME_blocktype][block_b.pos_y][block_b.pos_x])) : 0;
    SAD_d = block_d.available ? ((list==1) ? (fastme_l1_cost_flag[FME_blocktype][block_d.pos_y][block_d.pos_x]) : (fastme_l0_cost_flag[FME_blocktype][block_d.pos_y][block_d.pos_x])) : 0;
    SAD_c = block_c.available ? ((list==1) ? (fastme_l1_cost_flag[FME_blocktype][block_c.pos_y][block_c.pos_x]) : (fastme_l0_cost_flag[FME_blocktype][block_c.pos_y][block_c.pos_x])) : SAD_d;
  }
  for (hv=0; hv < 2; hv++)
  {
    if (!img->MbaffFrameFlag || hv==0)
    {
      mv_a = block_a.available  ? tmp_mv[block_a.pos_y][block_a.pos_x][hv] : 0;
      mv_b = block_b.available  ? tmp_mv[block_b.pos_y][block_b.pos_x][hv] : 0;
      mv_c = block_c.available  ? tmp_mv[block_c.pos_y][block_c.pos_x][hv] : 0;
    }
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
  }
  
  //Dynamic Search Range
  if (input->FMEnable == 1 && input->DynamicSearchRange)
    dsr_new_search_range = max(dsr_temp_search_range[0],dsr_temp_search_range[1]);
  
}

/*!
************************************************************************
* \brief
*    Initialize the motion search
************************************************************************
*/
void
Init_Motion_Search_Module ()
{
  int bits, i, imin, imax, k, l;
  
  int search_range               = input->search_range;
  int max_search_points          = max(9, (2*search_range+1)*(2*search_range+1));
  int max_ref_bits               = 1 + 2 * (int)floor(log(max(16,img->max_num_references+1)) / log(2) + 1e-10);
  int max_ref                    = (1<<((max_ref_bits>>1)+1))-1;
#ifdef EIGHTH_PEL
  int number_of_subpel_positions = (img->mv_res?8:4) * (2*search_range+3);
#else
  int number_of_subpel_positions = 4 * (2*search_range+3);
#endif
  int max_mv_bits                = 3 + 2 * (int)ceil (log(number_of_subpel_positions+1) / log(2) + 1e-10);
  max_mvd                        = (1<<( max_mv_bits >>1)   )-1;
  byte_abs_range                 = (img->max_imgpel_value > img->max_imgpel_value_uv) ? (img->max_imgpel_value + 1) * 64 : (img->max_imgpel_value_uv + 1) * 64;
  
  
  //=====   CREATE ARRAYS   =====
  //-----------------------------
  if ((spiral_search_x = (short*)calloc(max_search_points, sizeof(short))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: spiral_search_x");
  if ((spiral_search_y = (short*)calloc(max_search_points, sizeof(short))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: spiral_search_y");
  if ((spiral_hpel_search_x = (short*)calloc(max_search_points, sizeof(short))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: spiral_hpel_search_x");
  if ((spiral_hpel_search_y = (short*)calloc(max_search_points, sizeof(short))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: spiral_hpel_search_y");
  if ((mvbits = (int*)calloc(2*max_mvd+1, sizeof(int))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: mvbits");
  if ((refbits = (int*)calloc(max_ref, sizeof(int))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: refbits");
  if ((byte_abs = (unsigned int*)calloc(byte_abs_range, sizeof(unsigned int))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: byte_abs");
  
  get_mem4Dint (&motion_cost, 8, 2, img->max_num_references, 4);
  
  //--- set array offsets ---
  mvbits   += max_mvd;
  byte_abs += byte_abs_range/2;
  
  
  //=====   INIT ARRAYS   =====
  //---------------------------
  //--- init array: motion vector bits ---
  mvbits[0] = 1;
  for (bits=3; bits<=max_mv_bits; bits+=2)
  {
    imax = 1    << (bits >> 1);
    imin = imax >> 1;
    
    for (i = imin; i < imax; i++)   mvbits[-i] = mvbits[i] = bits;
  }
  //--- init array: reference frame bits ---
  refbits[0] = 1;
  for (bits=3; bits<=max_ref_bits; bits+=2)
  {
    imax = (1   << ((bits >> 1) + 1)) - 1;
    imin = imax >> 1;
    
    for (i = imin; i < imax; i++)   refbits[i] = bits;
  }
  //--- init array: absolute value ---
  byte_abs[0] = 0;
  for (i=1; i<byte_abs_range/2; i++)   byte_abs[i] = byte_abs[-i] = i;
  //--- init array: search pattern ---
  spiral_search_x[0] = spiral_search_y[0] = 0;
  spiral_hpel_search_x[0] = spiral_hpel_search_y[0] = 0;
  for (k=1, l=1; l<=max(1,search_range); l++)
  {
    for (i=-l+1; i< l; i++)
    {
      spiral_search_x[k] =  i;  spiral_search_y[k] = -l;
      spiral_hpel_search_x[k] =  i<<1;  spiral_hpel_search_y[k++] = -l<<1;
      spiral_search_x[k] =  i;  spiral_search_y[k] =  l;      
      spiral_hpel_search_x[k] =  i<<1;  spiral_hpel_search_y[k++] =  l<<1;
    }
    for (i=-l;   i<=l; i++)
    {
      spiral_search_x[k] = -l;  spiral_search_y[k] =  i;
      spiral_hpel_search_x[k] = -l<<1;  spiral_hpel_search_y[k++] = i<<1;
      spiral_search_x[k] =  l;  spiral_search_y[k] =  i;      
      spiral_hpel_search_x[k] =  l<<1;  spiral_hpel_search_y[k++] = i<<1;
    }
  }
  
#ifdef _FAST_FULL_ME_
  if(!input->FMEnable)
    InitializeFastFullIntegerSearch ();
#endif
}


/*!
************************************************************************
* \brief
*    Free memory used by motion search
************************************************************************
*/
void
Clear_Motion_Search_Module ()
{
  //--- correct array offset ---
  mvbits   -= max_mvd;
  byte_abs -= byte_abs_range/2;
  
  //--- delete arrays ---
  free (spiral_search_x);
  free (spiral_search_y);
  free (spiral_hpel_search_x);
  free (spiral_hpel_search_y);
  free (mvbits);
  free (refbits);
  free (byte_abs);
  free_mem4Dint (motion_cost, 8, 2);
  
#ifdef _FAST_FULL_ME_
  if(!input->FMEnable)
    ClearFastFullIntegerSearch ();
#endif
}



/*!
***********************************************************************
* \brief
*    Full pixel block motion search
***********************************************************************
*/
int                                               //  ==> minimum motion cost after search
FullPelBlockMotionSearch (pel_t**   orig_pic,     // <--  original pixel values for the AxB block
                          short     ref,          // <--  reference frame (0... or -1 (backward))
                          int       list,         // <--  current list
                          int       pic_pix_x,    // <--  absolute x-coordinate of regarded AxB block
                          int       pic_pix_y,    // <--  absolute y-coordinate of regarded AxB block
                          int       blocktype,    // <--  block type (1-16x16 ... 7-4x4)
                          short     pred_mv_x,    // <--  motion vector predictor (x) in sub-pel units
                          short     pred_mv_y,    // <--  motion vector predictor (y) in sub-pel units
                          short*    mv_x,         // <--> in: search center (x) / out: motion vector (x) - in pel units
                          short*    mv_y,         // <--> in: search center (y) / out: motion vector (y) - in pel units
                          int       search_range, // <--  1-d search range in pel units
                          int       min_mcost,    // <--  minimum motion cost (cost for center or huge value)
                          int       lambda_factor)       // <--  lagrangian parameter for determining motion cost
{
  int   pos, cand_x, cand_y, y, x4, mcost;
  
  pel_t *orig_line, *ref_line;
  pel_t *(*get_ref_line)(int, pel_t*, int, int, int, int);
  
  int   list_offset   = img->mb_data[img->current_mb_nr].list_offset; 
  pel_t *ref_pic      = listX[list+list_offset][ref]->imgY_11;
  int   img_width     = listX[list+list_offset][ref]->size_x;
  int   img_height    = listX[list+list_offset][ref]->size_y;
  
  int   best_pos      = 0;                                        // position with minimum motion cost
  int   max_pos       = (2*search_range+1)*(2*search_range+1);    // number of search positions
#ifdef EIGHTH_PEL
  int   mvshift       = (img->mv_res ? 3 : 2);                  // motion vector shift for getting sub-pel units
#endif
  int   blocksize_y   = input->blc_size[blocktype][1];            // vertical block size
  int   blocksize_x   = input->blc_size[blocktype][0];            // horizontal block size
  int   blocksize_x4  = blocksize_x >> 2;                         // horizontal block size in 4-pel units
#ifdef EIGHTH_PEL
  int   pred_x        = (pic_pix_x << mvshift) + pred_mv_x;       // predicted position x (in sub-pel units)
  int   pred_y        = (pic_pix_y << mvshift) + pred_mv_y;       // predicted position y (in sub-pel units)
#else
  int   pred_x        = (pic_pix_x << 2) + pred_mv_x;       // predicted position x (in sub-pel units)
  int   pred_y        = (pic_pix_y << 2) + pred_mv_y;       // predicted position y (in sub-pel units)
#endif
  int   center_x      = pic_pix_x + *mv_x;                        // center position x (in pel units)
  int   center_y      = pic_pix_y + *mv_y;                        // center position y (in pel units)
  int   check_for_00  = (blocktype==1 && !input->rdopt && img->type!=B_SLICE && ref==0);
  
  int  apply_weights = ( (active_pps->weighted_pred_flag && (img->type == P_SLICE || img->type == SP_SLICE)) ||
    (active_pps->weighted_bipred_idc && (img->type == B_SLICE)));  
  
  //===== Use weighted Reference for ME ====
  if (apply_weights && input->UseWeightedReferenceME)
    ref_pic       = listX[list+list_offset][ref]->imgY_11_w;
  else
    ref_pic       = listX[list+list_offset][ref]->imgY_11;
  
  //===== set function for getting reference picture lines =====
  if ((center_x > search_range) && (center_x < img_width -1-search_range-blocksize_x) &&
    (center_y > search_range) && (center_y < img_height-1-search_range-blocksize_y)   )
  {
    get_ref_line = FastLineX;
  }
  else
  {
    get_ref_line = UMVLineX;
  }
  
  
  //===== loop over all search positions =====
  for (pos=0; pos<max_pos; pos++)
  {
    //--- set candidate position (absolute position in pel units) ---
    cand_x = center_x + spiral_search_x[pos];
    cand_y = center_y + spiral_search_y[pos];
    
    //--- initialize motion cost (cost for motion vector) and check ---
#ifdef EIGHTH_PEL
    mcost = MV_COST (lambda_factor, mvshift, cand_x, cand_y, pred_x, pred_y);
#else
    mcost = MV_COST (lambda_factor, 2, cand_x, cand_y, pred_x, pred_y);
#endif
    
    if (check_for_00 && cand_x==pic_pix_x && cand_y==pic_pix_y)
    {
      mcost -= WEIGHTED_COST (lambda_factor, 16);
    }
    if (mcost >= min_mcost)   continue;
    
    //--- add residual cost to motion cost ---
    for (y=0; y<blocksize_y; y++)
    {
      ref_line  = get_ref_line (blocksize_x, ref_pic, cand_y+y, cand_x, img_height, img_width);
      orig_line = orig_pic [y];
      
      for (x4=0; x4<blocksize_x4; x4++)
      {
        mcost += byte_abs[ *orig_line++ - *ref_line++ ];
        mcost += byte_abs[ *orig_line++ - *ref_line++ ];
        mcost += byte_abs[ *orig_line++ - *ref_line++ ];
        mcost += byte_abs[ *orig_line++ - *ref_line++ ];
      }
      
      if (mcost >= min_mcost)
      {
        break;
      }
    }
    
    //--- check if motion cost is less than minimum cost ---
    if (mcost < min_mcost)
    {
      best_pos  = pos;
      min_mcost = mcost;
    }
  }
  
  
  //===== set best motion vector and return minimum motion cost =====
  if (best_pos)
  {
    *mv_x += spiral_search_x[best_pos];
    *mv_y += spiral_search_y[best_pos];
  }
  return min_mcost;
}


#ifdef _FAST_FULL_ME_
/*!
***********************************************************************
* \brief
*    Fast Full pixel block motion search
***********************************************************************
*/
int                                                   //  ==> minimum motion cost after search
FastFullPelBlockMotionSearch (pel_t**   orig_pic,     // <--  not used
                              short     ref,          // <--  reference frame (0... or -1 (backward))
                              int       list,
                              int       pic_pix_x,    // <--  absolute x-coordinate of regarded AxB block
                              int       pic_pix_y,    // <--  absolute y-coordinate of regarded AxB block
                              int       blocktype,    // <--  block type (1-16x16 ... 7-4x4)
                              short     pred_mv_x,    // <--  motion vector predictor (x) in sub-pel units
                              short     pred_mv_y,    // <--  motion vector predictor (y) in sub-pel units
                              short*    mv_x,         //  --> motion vector (x) - in pel units
                              short*    mv_y,         //  --> motion vector (y) - in pel units
                              int       search_range, // <--  1-d search range in pel units
                              int       min_mcost,    // <--  minimum motion cost (cost for center or huge value)
                              int       lambda_factor)       // <--  lagrangian parameter for determining motion cost
{
  int   pos, offset_x, offset_y, cand_x, cand_y, mcost;
  
#ifdef EIGHTH_PEL
  int   mvshift       = (img->mv_res ? 3 : 2);                            // motion vector shift for getting sub-pel units
#endif
  int   max_pos       = (2*search_range+1)*(2*search_range+1);              // number of search positions
  int   best_pos      = 0;                                                  // position with minimum motion cost
  int   block_index;                                                        // block index for indexing SAD array
  distpel*  block_sad;                                                          // pointer to SAD array
  
  block_index   = (pic_pix_y-img->opix_y)+((pic_pix_x-img->opix_x)>>2); // block index for indexing SAD array
  block_sad     = BlockSAD[list][ref][blocktype][block_index];         // pointer to SAD array
  
  //===== set up fast full integer search if needed / set search center =====
  if (!search_setup_done[list][ref])
  {
    SetupFastFullPelSearch (ref, list);
  }
  
  offset_x = search_center_x[list][ref] - img->opix_x;
  offset_y = search_center_y[list][ref] - img->opix_y;
  
  //===== cost for (0,0)-vector: it is done before, because MVCost can be negative =====
  if (!input->rdopt)
  {
#ifdef EIGHTH_PEL
    mcost = block_sad[pos_00[list][ref]] + MV_COST (lambda_factor, mvshift, 0, 0, pred_mv_x, pred_mv_y);
#else
    mcost = block_sad[pos_00[list][ref]] + MV_COST (lambda_factor, 2, 0, 0, pred_mv_x, pred_mv_y);
#endif
    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      best_pos  = pos_00[list][ref];
    }
  }
  
  //===== loop over all search positions =====
  for (pos=0; pos<max_pos; pos++, block_sad++)
  {
    //--- check residual cost ---
    if (*block_sad < min_mcost)
    {
      //--- get motion vector cost ---
      cand_x = offset_x + spiral_search_x[pos];
      cand_y = offset_y + spiral_search_y[pos];
      mcost  = *block_sad;
#ifdef EIGHTH_PEL
#ifdef MV_COMPETITION
      if (input->mv_competition > 0)
        mcost += compute_mv_cost_for_each_predictor (lambda_factor, mvshift, cand_x, cand_y, pred_mv_x, pred_mv_y, ref, list, blocktype, pic_pix_x-img->opix_x, pic_pix_y-img->opix_y);
      else
#endif
        mcost += MV_COST (lambda_factor, mvshift, cand_x, cand_y, pred_mv_x, pred_mv_y);
      
      
#else
      
      mcost += MV_COST (lambda_factor, 2, cand_x, cand_y, pred_mv_x, pred_mv_y);
      
#endif
      
      //--- check motion cost ---
      if (mcost < min_mcost)
      {
        min_mcost = mcost;
        best_pos  = pos;
      }
    }
  }
  
  //===== set best motion vector and return minimum motion cost =====
  *mv_x = offset_x + spiral_search_x[best_pos];
  *mv_y = offset_y + spiral_search_y[best_pos];
  return min_mcost;
}
#endif


/*!
***********************************************************************
* \brief
*    Calculate SA(T)D
***********************************************************************
*/
int
SATD (int* diff, int use_hadamard)
{
  int k, satd = 0, m[16], d[16];
#ifdef ADAPTIVE_FD_SD_CODING
  int sad=0;
  int i;
#endif
  
  if (use_hadamard)
  {
#ifdef ADAPTIVE_FD_SD_CODING
    if (img->type!=I_SLICE)
    {
      for (i = 0; i < 16; i++)
      {
        sad += 2*byte_abs [diff [i]];
      }
    }
#endif
    /*===== hadamard transform =====*/
    m[ 0] = diff[ 0] + diff[12];
    m[ 1] = diff[ 1] + diff[13];
    m[ 2] = diff[ 2] + diff[14];
    m[ 3] = diff[ 3] + diff[15];
    m[ 4] = diff[ 4] + diff[ 8];
    m[ 5] = diff[ 5] + diff[ 9];
    m[ 6] = diff[ 6] + diff[10];
    m[ 7] = diff[ 7] + diff[11];
    m[ 8] = diff[ 4] - diff[ 8];
    m[ 9] = diff[ 5] - diff[ 9];
    m[10] = diff[ 6] - diff[10];
    m[11] = diff[ 7] - diff[11];
    m[12] = diff[ 0] - diff[12];
    m[13] = diff[ 1] - diff[13];
    m[14] = diff[ 2] - diff[14];
    m[15] = diff[ 3] - diff[15];
    
    d[ 0] = m[ 0] + m[ 4];
    d[ 1] = m[ 1] + m[ 5];
    d[ 2] = m[ 2] + m[ 6];
    d[ 3] = m[ 3] + m[ 7];
    d[ 4] = m[ 8] + m[12];
    d[ 5] = m[ 9] + m[13];
    d[ 6] = m[10] + m[14];
    d[ 7] = m[11] + m[15];
    d[ 8] = m[ 0] - m[ 4];
    d[ 9] = m[ 1] - m[ 5];
    d[10] = m[ 2] - m[ 6];
    d[11] = m[ 3] - m[ 7];
    d[12] = m[12] - m[ 8];
    d[13] = m[13] - m[ 9];
    d[14] = m[14] - m[10];
    d[15] = m[15] - m[11];
    
    m[ 0] = d[ 0] + d[ 3];
    m[ 1] = d[ 1] + d[ 2];
    m[ 2] = d[ 1] - d[ 2];
    m[ 3] = d[ 0] - d[ 3];
    m[ 4] = d[ 4] + d[ 7];
    m[ 5] = d[ 5] + d[ 6];
    m[ 6] = d[ 5] - d[ 6];
    m[ 7] = d[ 4] - d[ 7];
    m[ 8] = d[ 8] + d[11];
    m[ 9] = d[ 9] + d[10];
    m[10] = d[ 9] - d[10];
    m[11] = d[ 8] - d[11];
    m[12] = d[12] + d[15];
    m[13] = d[13] + d[14];
    m[14] = d[13] - d[14];
    m[15] = d[12] - d[15];
    
    d[ 0] = m[ 0] + m[ 1];
    d[ 1] = m[ 0] - m[ 1];
    d[ 2] = m[ 2] + m[ 3];
    d[ 3] = m[ 3] - m[ 2];
    d[ 4] = m[ 4] + m[ 5];
    d[ 5] = m[ 4] - m[ 5];
    d[ 6] = m[ 6] + m[ 7];
    d[ 7] = m[ 7] - m[ 6];
    d[ 8] = m[ 8] + m[ 9];
    d[ 9] = m[ 8] - m[ 9];
    d[10] = m[10] + m[11];
    d[11] = m[11] - m[10];
    d[12] = m[12] + m[13];
    d[13] = m[12] - m[13];
    d[14] = m[14] + m[15];
    d[15] = m[15] - m[14];
    
    //===== sum up =====
    // Table lookup is faster than abs macro
    if (input->FMEnable == 2)
    {
      for (k=0; k<16; ++k)
      {
        satd += byte_abs [d [k]];
      }
    }
    else
    {
      for (k=0; k<16; ++k)
      {
        //satd += absm(d[k]);
        satd += byte_abs [d [k]];
      }
    }
    satd = ((satd+1)>>1);
    
#ifdef ADAPTIVE_FD_SD_CODING
    if (img->type!=I_SLICE && img->APEC_in_FD_and_SD && sad<satd)
    {
      satd=sad;
    }
#endif
  }
  else
  {
    /*===== sum up =====*/
    for (k = 0; k < 16; k++)
    {
      satd += byte_abs [diff [k]];
    }
  }
  
  return satd;
}

/*!
***********************************************************************
* \brief
*    Calculate SA(T)D for 8x8
***********************************************************************
*/
int
SATD8X8 (int* diff, int use_hadamard)
{
  int i, j, jj, sad=0;
  int m1[8][8], m2[8][8], m3[8][8];
#ifdef ADAPTIVE_FD_SD_CODING
  int satd=0;
  if (img->type!=I_SLICE)
  {
    for (i = 0; i < 64; i++)
    {
      satd += 2*byte_abs [diff [i]];
    }
  }
#endif
  
  if(use_hadamard)
  {
    //horizontal
    for (j=0; j < 8; j++)
    {
      jj = j << 3;
      m2[j][0] = diff[jj  ] + diff[jj+4];
      m2[j][1] = diff[jj+1] + diff[jj+5];
      m2[j][2] = diff[jj+2] + diff[jj+6];
      m2[j][3] = diff[jj+3] + diff[jj+7];
      m2[j][4] = diff[jj  ] - diff[jj+4];
      m2[j][5] = diff[jj+1] - diff[jj+5];
      m2[j][6] = diff[jj+2] - diff[jj+6];
      m2[j][7] = diff[jj+3] - diff[jj+7];
      
      m1[j][0] = m2[j][0] + m2[j][2];
      m1[j][1] = m2[j][1] + m2[j][3];
      m1[j][2] = m2[j][0] - m2[j][2];
      m1[j][3] = m2[j][1] - m2[j][3];
      m1[j][4] = m2[j][4] + m2[j][6];
      m1[j][5] = m2[j][5] + m2[j][7];
      m1[j][6] = m2[j][4] - m2[j][6];
      m1[j][7] = m2[j][5] - m2[j][7];
      
      m2[j][0] = m1[j][0] + m1[j][1];
      m2[j][1] = m1[j][0] - m1[j][1];
      m2[j][2] = m1[j][2] + m1[j][3];
      m2[j][3] = m1[j][2] - m1[j][3];
      m2[j][4] = m1[j][4] + m1[j][5];
      m2[j][5] = m1[j][4] - m1[j][5];
      m2[j][6] = m1[j][6] + m1[j][7];
      m2[j][7] = m1[j][6] - m1[j][7];
    }
    
    //vertical
    for (i=0; i < 8; i++)
    {
      m3[0][i] = m2[0][i] + m2[4][i];
      m3[1][i] = m2[1][i] + m2[5][i];
      m3[2][i] = m2[2][i] + m2[6][i];
      m3[3][i] = m2[3][i] + m2[7][i];
      m3[4][i] = m2[0][i] - m2[4][i];
      m3[5][i] = m2[1][i] - m2[5][i];
      m3[6][i] = m2[2][i] - m2[6][i];
      m3[7][i] = m2[3][i] - m2[7][i];
      
      m1[0][i] = m3[0][i] + m3[2][i];
      m1[1][i] = m3[1][i] + m3[3][i];
      m1[2][i] = m3[0][i] - m3[2][i];
      m1[3][i] = m3[1][i] - m3[3][i];
      m1[4][i] = m3[4][i] + m3[6][i];
      m1[5][i] = m3[5][i] + m3[7][i];
      m1[6][i] = m3[4][i] - m3[6][i];
      m1[7][i] = m3[5][i] - m3[7][i];
      
      m2[0][i] = m1[0][i] + m1[1][i];
      m2[1][i] = m1[0][i] - m1[1][i];
      m2[2][i] = m1[2][i] + m1[3][i];
      m2[3][i] = m1[2][i] - m1[3][i];
      m2[4][i] = m1[4][i] + m1[5][i];
      m2[5][i] = m1[4][i] - m1[5][i];
      m2[6][i] = m1[6][i] + m1[7][i];
      m2[7][i] = m1[6][i] - m1[7][i];
    }
    for (j=0; j < 8; j++)
      for (i=0; i < 8; i++)
        sad += (absm(m2[j][i]));
      
      sad=((sad+2)>>2);
      
#ifdef ADAPTIVE_FD_SD_CODING
      if (img->type!=I_SLICE && img->APEC_in_FD_and_SD && satd<sad)
      {
        sad=satd;
      }
#endif
  }
  else
  {
    for(i=0; i<64; i++)
      sad += byte_abs [diff [i]];
  }
  
  return sad;
}

/*!
***********************************************************************
* \brief
*    Calculate SA(T)D for 8x8
***********************************************************************
*/
int
find_SATD (int c_diff[MB_PIXELS], int blocktype)
{
  int sad=0;
  
  switch(blocktype)
  {
    //16x16
  case 1: 
    sad  = SATD8X8 (c_diff,       input->hadamard);
    sad += SATD8X8 (&c_diff[ 64], input->hadamard);
    sad += SATD8X8 (&c_diff[128], input->hadamard);
    sad += SATD8X8 (&c_diff[192], input->hadamard);
    break;
    //16x8 8x16
  case 2:
  case 3: 
    sad  = SATD8X8 (c_diff,      input->hadamard);
    sad += SATD8X8 (&c_diff[64], input->hadamard);
    break;
    //8x8
  case 4: 
    sad  = SATD8X8 (c_diff, input->hadamard);
    break;
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
*    Sub pixel block motion search
***********************************************************************
*/
int                                               //  ==> minimum motion cost after search
SubPelBlockMotionSearch (pel_t**   orig_pic,      // <--  original pixel values for the AxB block
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
                         int       search_pos8,   // <--  search positions for  eighth-pel search  (default: 9)
#endif
                         int       min_mcost,     // <--  minimum motion cost (cost for center or huge value)
                         int       lambda_factor  // <--  lagrangian parameter for determining motion cost
                         )
{
  int   j, i, k;
  int   c_diff[MB_PIXELS];
  int   diff[16], *d;
  int   pos, best_pos, mcost, abort_search;
  int   y0, y1, y2, y3;
  int   x0;
  int   ry0, ry4, ry8, ry12, rx0;
  
  int   cand_mv_x, cand_mv_y;
  int   y_offset, ypels =(128 - 64 * (blocktype == 3));
#ifdef EIGHTH_PEL
  int   mv_shift        = (img->mv_res ? 1 : 0);
#endif
  
  int   check_position0 = (!input->rdopt && img->type!=B_SLICE && ref==0 && blocktype==1 && *mv_x==0 && *mv_y==0 && input->hadamard);
  int   blocksize_x     = input->blc_size[blocktype][0];
  int   blocksize_y     = input->blc_size[blocktype][1];
  int   pic4_pix_x      = ((pic_pix_x + IMG_PAD_SIZE)<< 2);
  int   pic4_pix_y      = ((pic_pix_y + IMG_PAD_SIZE)<< 2);
  int   min_pos2        = (input->hadamard == 1 ? 0 : 1);
  int   max_pos2        = (input->hadamard ? max(1,search_pos2) : search_pos2);
  int   list_offset     = img->mb_data[img->current_mb_nr].list_offset; 
  int   apply_weights   = ((active_pps->weighted_pred_flag && (img->type == P_SLICE || img->type == SP_SLICE)) ||
    (active_pps->weighted_bipred_idc && (img->type == B_SLICE)));  
  int   halfpelhadamard  = input->hadamard == 2 ? 0 : input->hadamard;
  int   qpelstart        = input->hadamard == 2 ? 0 : 1;
  int   test8x8transform = input->Transform8x8Mode && blocktype <= 4 && halfpelhadamard;
  int   cmv_x, cmv_y;
  
  StorablePicture *ref_picture = listX[list+list_offset][ref];
  pel_t **ref_pic = (apply_weights && input->UseWeightedReferenceME)? ref_picture->imgY_ups_w : ref_picture->imgY_ups;      
  pel_t *ref_line;
  pel_t *orig_line;  
  int img_width  = ((ref_picture->size_x + 2*IMG_PAD_SIZE - 1)<<2);
  int img_height = ((ref_picture->size_y + 2*IMG_PAD_SIZE - 1)<<2);
  int max_pos_x4 = ((ref_picture->size_x - blocksize_x + 2*IMG_PAD_SIZE)<<2);
  int max_pos_y4 = ((ref_picture->size_y - blocksize_y + 2*IMG_PAD_SIZE)<<2);
  
  /*********************************
  *****                       *****
  *****  HALF-PEL REFINEMENT  *****
  *****                       *****
  *********************************/
  
  //===== set function for getting pixel values =====
  if ((pic4_pix_x + *mv_x > 1) && (pic4_pix_x + *mv_x < max_pos_x4 - 1) &&
    (pic4_pix_y + *mv_y > 1) && (pic4_pix_y + *mv_y < max_pos_y4 - 1)   )
  {
    get_line = FastLine4X;
  }
  else
  {
    get_line = UMVLine4X;    
  }
  
  //===== loop over search positions =====
  for (best_pos = 0, pos = min_pos2; pos < max_pos2; pos++)
  {
    cand_mv_x = *mv_x + (spiral_hpel_search_x[pos]);    // quarter-pel units
    cand_mv_y = *mv_y + (spiral_hpel_search_y[pos]);    // quarter-pel units
    
    //----- set motion vector cost -----
    
    
    
#ifdef EIGHTH_PEL
    mcost = MV_COST (lambda_factor, mv_shift, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);
#else
    mcost = MV_COST (lambda_factor, 0, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);
#endif
    
    
    if (check_position0 && pos==0)
    {
      mcost -= WEIGHTED_COST (lambda_factor, 16);
    }
    
    if (mcost >= min_mcost) continue;
    
    cmv_x = cand_mv_x + pic4_pix_x;
    cmv_y = cand_mv_y + pic4_pix_y;
    
    //----- add up SATD -----
    for (y0=0, abort_search=0; y0<blocksize_y && !abort_search; y0+=4)
    {
      y_offset = (y0>7)*ypels;
      ry0  = (y0<<2) + cmv_y;
      ry4  = ry0 + 4;
      ry8  = ry4 + 4;
      ry12 = ry8 + 4;
      y1 = y0 + 1;
      y2 = y1 + 1;
      y3 = y2 + 1;
      
      for (x0=0; x0<blocksize_x; x0+=BLOCK_SIZE)
      {
        rx0 = (x0<<2) + cmv_x;
        d   = diff;
        
        orig_line = &orig_pic [y0][x0];    
        ref_line  = get_line (ref_pic, ry0, rx0, img_height, img_width);
        *d++      = *orig_line++ - *(ref_line     );
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line   - *(ref_line += 4);
        
        orig_line = &orig_pic [y1][x0];    
        ref_line  = get_line (ref_pic, ry4, rx0, img_height, img_width);
        *d++      = *orig_line++ - *(ref_line     );
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line   - *(ref_line += 4);
        
        orig_line = &orig_pic [y2][x0];
        ref_line  = get_line (ref_pic, ry8, rx0, img_height, img_width);
        *d++      = *orig_line++ - *(ref_line     );
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line   - *(ref_line += 4);
        
        orig_line = &orig_pic [y3][x0];    
        ref_line  = get_line (ref_pic, ry12, rx0, img_height, img_width);
        *d++      = *orig_line++ - *(ref_line     );
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line   - *(ref_line += 4);
        
        if (!test8x8transform)
        {
          if ((mcost += SATD (diff, halfpelhadamard)) > min_mcost)
          {
            abort_search = 1;
            break;
          }
        }
        else
        {
          i = (x0&0x7) +  (x0>7) * 64 + y_offset;
          for(k=0, j=y0; j<BLOCK_SIZE + y0; j++, k+=BLOCK_SIZE)
            memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], BLOCK_SIZE*sizeof(int));
        }
      }
    }
    
    
    if(test8x8transform)
      mcost += find_SATD (c_diff, blocktype);
    
    
    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      best_pos  = pos;
    }
  }
  if (best_pos)
  {
    *mv_x += (spiral_hpel_search_x [best_pos]);
    *mv_y += (spiral_hpel_search_y [best_pos]);
  }
  
  if (input->hadamard == 2)
    min_mcost = INT_MAX;
  
  test8x8transform = input->Transform8x8Mode && blocktype <= 4 && input->hadamard;
  
  /************************************
  *****                          *****
  *****  QUARTER-PEL REFINEMENT  *****
  *****                          *****
  ************************************/
  //===== set function for getting pixel values =====
  if ((pic4_pix_x + *mv_x > 0) && (pic4_pix_x + *mv_x < max_pos_x4) &&
    (pic4_pix_y + *mv_y > 0) && (pic4_pix_y + *mv_y < max_pos_y4)   )
  {
    get_line = FastLine4X;
  }
  else
  {
    get_line = UMVLine4X;    
  }
  
  //===== loop over search positions =====
  for (best_pos = 0, pos = qpelstart; pos < search_pos4; pos++)
  {
    cand_mv_x = *mv_x + spiral_search_x[pos];    // quarter-pel units
    cand_mv_y = *mv_y + spiral_search_y[pos];    // quarter-pel units
    
    //----- set motion vector cost -----
#ifdef EIGHTH_PEL
    mcost = MV_COST (lambda_factor, mv_shift, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);
#else
    mcost = MV_COST (lambda_factor, 0, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);
#endif
    
    
    if (mcost >= min_mcost) continue;
    cmv_x = cand_mv_x + pic4_pix_x;
    cmv_y = cand_mv_y + pic4_pix_y;
    
    //----- add up SATD -----
    for (y0=0, abort_search=0; y0<blocksize_y && !abort_search; y0+=4)
    {
      y_offset = (y0>7)*ypels;
      ry0 = (y0<<2) + cmv_y;
      ry4  = ry0 + 4;
      ry8  = ry4 + 4;
      ry12 = ry8 + 4;
      y1 = y0 + 1;
      y2 = y1 + 1;
      y3 = y2 + 1;
      
      for (x0=0; x0<blocksize_x; x0+=BLOCK_SIZE)
      {
        rx0  = (x0<<2) + cmv_x;
        d    = diff;
        
        orig_line = &orig_pic [y0][x0];    
        ref_line  = get_line (ref_pic, ry0, rx0, img_height, img_width);
        *d++      = *orig_line++ - *(ref_line     );
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line   - *(ref_line += 4);
        
        orig_line = &orig_pic [y1][x0];    
        ref_line  = get_line (ref_pic, ry4, rx0, img_height, img_width);
        *d++      = *orig_line++ - *(ref_line     );
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line   - *(ref_line += 4);
        
        orig_line = &orig_pic [y2][x0];
        ref_line  = get_line (ref_pic, ry8, rx0, img_height, img_width);
        *d++      = *orig_line++ - *(ref_line     );
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line   - *(ref_line += 4);
        
        orig_line = &orig_pic [y3][x0];    
        ref_line  = get_line (ref_pic, ry12, rx0, img_height, img_width);
        *d++      = *orig_line++ - *(ref_line     );
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line++ - *(ref_line += 4);
        *d        = *orig_line   - *(ref_line += 4);
        
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
          i = (x0&0x7) + (x0>7) * 64 + y_offset;
          for(k=0, j=y0; j<y0 + BLOCK_SIZE; j++, k+=BLOCK_SIZE)
            memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], BLOCK_SIZE*sizeof(int));
        }
      }
    }
    
    if(test8x8transform)
      mcost += find_SATD (c_diff, blocktype);
    
    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      best_pos  = pos;
    }
  }
  if (best_pos)
  {
    *mv_x += spiral_search_x [best_pos];
    *mv_y += spiral_search_y [best_pos];
  }
  
  
#ifdef EIGHTH_PEL
  /***********************************
  *****                         *****
  *****  EIGHTH-PEL REFINEMENT  *****
  *****                         *****
  ***********************************/
  if (img->mv_res)
  {
    //===== convert search center to eighth-pel units (double resolution) =====
    *mv_x <<= 1;
    *mv_y <<= 1;
    pic4_pix_x <<= 1;
    pic4_pix_y <<= 1;
    
    //===== set function for getting pixel values =====
    if ((pic4_pix_x + *mv_x > 1) && (pic4_pix_x + *mv_x < 2*max_pos_x4 - 2) &&
      (pic4_pix_y + *mv_y > 1) && (pic4_pix_y + *mv_y < 2*max_pos_y4 - 2))
      get_line = FastLine4X8;
    else
      get_line = UMVLine4X8;
    
    //===== loop over search positions =====
    for (best_pos = 0, pos = 1; pos < search_pos8; pos++)
    {
      cand_mv_x = *mv_x + spiral_search_x[pos];    // eighth-pel units
      cand_mv_y = *mv_y + spiral_search_y[pos];    // eighth-pel units
      
      //----- set motion vector cost -----
      
      mcost = MV_COST (lambda_factor, 0, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);
      
      if (mcost >= min_mcost) continue;
      cmv_x = cand_mv_x + pic4_pix_x;
      cmv_y = cand_mv_y + pic4_pix_y;
      
      //----- add up SATD -----
      for (y0=0, abort_search=0; y0<blocksize_y && !abort_search; y0+=BLOCK_SIZE)
      {
        y_offset = (y0>7)*ypels;
        ry0 = (y0<<3) + cmv_y;
        ry4  = ry0 + 8;
        ry8  = ry4 + 8;
        ry12 = ry8 + 8;
        y1 = y0 + 1;
        y2 = y1 + 1;
        y3 = y2 + 1;
        
        for (x0=0; x0<blocksize_x; x0+=BLOCK_SIZE)
        {
          rx0 = (x0<<3) + cmv_x;
          d   = diff;
          
          orig_line = &orig_pic [y0][x0];
          ref_line  = get_line (ref_pic, ry0, rx0, img_height, img_width);
          *d++      = *orig_line++ - *(ref_line      );
          *d++      = *orig_line++ - *(ref_line += 4 );
          *d++      = *orig_line++ - *(ref_line += 4 );
          *d++      = *orig_line   - *(ref_line += 4 );
          
          orig_line = &orig_pic [y1][x0];
          ref_line  = get_line (ref_pic, ry4, rx0, img_height, img_width);
          *d++      = *orig_line++ - *(ref_line      );
          *d++      = *orig_line++ - *(ref_line += 4 );
          *d++      = *orig_line++ - *(ref_line += 4 );
          *d++      = *orig_line   - *(ref_line += 4 );
          
          orig_line = &orig_pic [y2][x0];
          ref_line  = get_line (ref_pic, ry8, rx0, img_height, img_width);
          *d++      = *orig_line++ - *(ref_line      );
          *d++      = *orig_line++ - *(ref_line += 4 );
          *d++      = *orig_line++ - *(ref_line += 4 );
          *d++      = *orig_line   - *(ref_line += 4 );
          
          orig_line = &orig_pic [y3][x0];
          ref_line  = get_line (ref_pic, ry12, rx0, img_height, img_width);
          *d++      = *orig_line++ - *(ref_line      );
          *d++      = *orig_line++ - *(ref_line += 4 );
          *d++      = *orig_line++ - *(ref_line += 4 );
          *d        = *orig_line   - *(ref_line += 4 );
          
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
            i = (x0&0x7) +  (x0>7) * 64 + y_offset;
            for(k=0, j=y0; j<BLOCK_SIZE + y0; j++, k+=BLOCK_SIZE)
              memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], BLOCK_SIZE*sizeof(int));
          }
        }
      }
      if (test8x8transform)
        mcost += find_SATD (c_diff, blocktype);
      
      if (mcost < min_mcost)
      {
        min_mcost = mcost;
        best_pos  = pos;
      }
    }
    if (best_pos)
    {
      *mv_x += spiral_search_x [best_pos];
      *mv_y += spiral_search_y [best_pos];
    }
  }
#endif
  
  //===== return minimum motion cost =====
  return min_mcost;
}


/*!
***********************************************************************
* \brief
*    Sub pixel block motion search enhanced
***********************************************************************
*/
int                                               //  ==> minimum motion cost after search
simplified_FastFullSubPelBlockMotionSearch (pel_t**   orig_pic,      // <--  original pixel values for the AxB block
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
                                            int       lambda_factor  // <--  lagrangian parameter for determining motion cost
                                            )
{
  int   j, i, k;
  int   c_diff[MB_PIXELS];
  int   diff[16], *d;
  int   pos, best_pos, mcost, abort_search;
  int   y0, y1, y2, y3;
  int   x0;
  int   ry0, ry4, ry8, ry12, rx0;
#ifdef EIGHTH_PEL
  int   mv_shift        = (img->mv_res ? 1 : 0);
#endif
  int   cand_mv_x, cand_mv_y;
  int   y_offset, ypels =(128 - 64 * (blocktype == 3));
  
  int   check_position0 = (!input->rdopt && img->type!=B_SLICE && ref==0 && blocktype==1 && *mv_x==0 && *mv_y==0 && input->hadamard);
  int   blocksize_x     = input->blc_size[blocktype][0];
  int   blocksize_y     = input->blc_size[blocktype][1];
  int   pic4_pix_x      = ((pic_pix_x + IMG_PAD_SIZE)<< 2);
  int   pic4_pix_y      = ((pic_pix_y + IMG_PAD_SIZE)<< 2);
  int   min_pos2        = (input->hadamard == 1 ? 0 : 1);
  int   max_pos2        = (input->hadamard ? max(1,search_pos2) : search_pos2);
  int   list_offset     = img->mb_data[img->current_mb_nr].list_offset; 
  int   apply_weights   = ((active_pps->weighted_pred_flag && (img->type == P_SLICE || img->type == SP_SLICE)) ||
    (active_pps->weighted_bipred_idc && (img->type == B_SLICE)));  
  int   halfpelhadamard  = input->hadamard == 2 ? 0 : input->hadamard;
  int   qpelstart        = input->hadamard == 2 ? 0 : 1;
  int   test8x8transform = input->Transform8x8Mode && blocktype <= 4 && halfpelhadamard;
  int   cmv_x, cmv_y;
  
  StorablePicture *ref_picture = listX[list+list_offset][ref];
  pel_t **ref_pic = (apply_weights && input->UseWeightedReferenceME)? ref_picture->imgY_ups_w : ref_picture->imgY_ups;      
  pel_t *ref_line;
  pel_t *orig_line;  
  int img_width  = ((ref_picture->size_x + 2*IMG_PAD_SIZE - 1)<<2);
  int img_height = ((ref_picture->size_y + 2*IMG_PAD_SIZE - 1)<<2);
  int max_pos_x4 = ((ref_picture->size_x - blocksize_x + 2*IMG_PAD_SIZE)<<2);
  int max_pos_y4 = ((ref_picture->size_y - blocksize_y + 2*IMG_PAD_SIZE)<<2);
  
  /*********************************
  *****                       *****
  *****  HALF-PEL REFINEMENT  *****
  *****                       *****
  *********************************/
  
  //===== set function for getting pixel values =====
  if ((pic4_pix_x + *mv_x > 1) && (pic4_pix_x + *mv_x < max_pos_x4 - 1) &&
    (pic4_pix_y + *mv_y > 1) && (pic4_pix_y + *mv_y < max_pos_y4 - 1)   )
  {
    get_line = FastLine4X;
  }
  else
  {
    get_line = UMVLine4X;    
  }
  
  //===== loop over search positions =====
  for (best_pos = 0, pos = min_pos2; pos < max_pos2; pos++)
  {
    cand_mv_x = *mv_x + (spiral_hpel_search_x[pos]);    // quarter-pel units
    cand_mv_y = *mv_y + (spiral_hpel_search_y[pos]);    // quarter-pel units
    
    //----- set motion vector cost -----
#ifdef EIGHTH_PEL
    mcost = MV_COST (lambda_factor, mv_shift, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);
#else
    mcost = MV_COST (lambda_factor, 0, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);
#endif
    
    if (check_position0 && pos==0)
    {
      mcost -= WEIGHTED_COST (lambda_factor, 16);
    }
    
    if (mcost >= min_mcost) continue;
    
    cmv_x = cand_mv_x + pic4_pix_x;
    cmv_y = cand_mv_y + pic4_pix_y;
    
    //----- add up SATD -----
    for (y0=0, abort_search=0; y0<blocksize_y && !abort_search; y0+=4)
    {
      y_offset = (y0>7)*ypels;
      ry0  = (y0<<2) + cmv_y;
      ry4  = ry0 + 4;
      ry8  = ry4 + 4;
      ry12 = ry8 + 4;
      y1 = y0 + 1;
      y2 = y1 + 1;
      y3 = y2 + 1;
      
      for (x0=0; x0<blocksize_x; x0+=BLOCK_SIZE)
      {
        rx0 = (x0<<2) + cmv_x;
        d   = diff;
        
        orig_line = &orig_pic [y0][x0];    
        ref_line  = get_line (ref_pic, ry0, rx0, img_height, img_width);
        *d++      = *orig_line++ - *(ref_line     );
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line   - *(ref_line += 4);
        
        orig_line = &orig_pic [y1][x0];    
        ref_line  = get_line (ref_pic, ry4, rx0, img_height, img_width);
        *d++      = *orig_line++ - *(ref_line     );
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line   - *(ref_line += 4);
        
        orig_line = &orig_pic [y2][x0];
        ref_line  = get_line (ref_pic, ry8, rx0, img_height, img_width);
        *d++      = *orig_line++ - *(ref_line     );
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line   - *(ref_line += 4);
        
        orig_line = &orig_pic [y3][x0];    
        ref_line  = get_line (ref_pic, ry12, rx0, img_height, img_width);
        *d++      = *orig_line++ - *(ref_line     );
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line   - *(ref_line += 4);
        
        if (!test8x8transform)
        {
          if ((mcost += SATD (diff, halfpelhadamard)) >= min_mcost)
          {
            abort_search = 1;
            break;
          }
        }
        else
        {
          i = (x0&0x7) +  (x0>7) * 64 + y_offset;
          for(k=0, j=y0; j<BLOCK_SIZE + y0; j++, k+=BLOCK_SIZE)
            memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], BLOCK_SIZE*sizeof(int));
        }
      }
    }
    
    if(test8x8transform)
      mcost += find_SATD (c_diff, blocktype);
    
    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      best_pos  = pos;
    }
    if (min_mcost < (SubPelThreshold3>>block_type_shift_factor[blocktype])) 
    {
      break;
    }
  }
  
  if (best_pos)
  {
    *mv_x += (spiral_hpel_search_x [best_pos]);
    *mv_y += (spiral_hpel_search_y [best_pos]);
  }
  
  if ((*mv_x == 0) && (*mv_y == 0) && (pred_mv_x == 0 && pred_mv_y == 0) &&
    (min_mcost < (SubPelThreshold1>>block_type_shift_factor[blocktype])) ) 
  {
    best_pos = 0;
    return min_mcost;
  }
  
  if (input->hadamard == 2)
    min_mcost = INT_MAX;
  
  test8x8transform = input->Transform8x8Mode && blocktype <= 4 && input->hadamard;
  
  /************************************
  *****                          *****
  *****  QUARTER-PEL REFINEMENT  *****
  *****                          *****
  ************************************/
  //===== set function for getting pixel values =====
  if ((pic4_pix_x + *mv_x > 0) && (pic4_pix_x + *mv_x < max_pos_x4) &&
    (pic4_pix_y + *mv_y > 0) && (pic4_pix_y + *mv_y < max_pos_y4)   )
  {
    get_line = FastLine4X;
  }
  else
  {
    get_line = UMVLine4X;    
  }
  
  //===== loop over search positions =====
  for (best_pos = 0, pos = qpelstart; pos < search_pos4; pos++)
  {
    cand_mv_x = *mv_x + spiral_search_x[pos];    // quarter-pel units
    cand_mv_y = *mv_y + spiral_search_y[pos];    // quarter-pel units
    
    //----- set motion vector cost -----
#ifdef EIGHTH_PEL
    mcost = MV_COST (lambda_factor, mv_shift, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);
#else
    mcost = MV_COST (lambda_factor, 0, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);
#endif
    
    if (mcost >= min_mcost) continue;
    cmv_x = cand_mv_x + pic4_pix_x;
    cmv_y = cand_mv_y + pic4_pix_y;
    
    //----- add up SATD -----
    for (y0=0, abort_search=0; y0<blocksize_y && !abort_search; y0+=4)
    {
      y_offset = (y0>7)*ypels;
      ry0 = (y0<<2) + cmv_y;
      ry4  = ry0 + 4;
      ry8  = ry4 + 4;
      ry12 = ry8 + 4;
      y1 = y0 + 1;
      y2 = y1 + 1;
      y3 = y2 + 1;
      
      for (x0=0; x0<blocksize_x; x0+=BLOCK_SIZE)
      {
        rx0  = (x0<<2) + cmv_x;
        d    = diff;
        
        orig_line = &orig_pic [y0][x0];    
        ref_line  = get_line (ref_pic, ry0, rx0, img_height, img_width);
        *d++      = *orig_line++ - *(ref_line     );
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line   - *(ref_line += 4);
        
        orig_line = &orig_pic [y1][x0];    
        ref_line  = get_line (ref_pic, ry4, rx0, img_height, img_width);
        *d++      = *orig_line++ - *(ref_line     );
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line   - *(ref_line += 4);
        
        orig_line = &orig_pic [y2][x0];
        ref_line  = get_line (ref_pic, ry8, rx0, img_height, img_width);
        *d++      = *orig_line++ - *(ref_line     );
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line   - *(ref_line += 4);
        
        orig_line = &orig_pic [y3][x0];    
        ref_line  = get_line (ref_pic, ry12, rx0, img_height, img_width);
        *d++      = *orig_line++ - *(ref_line     );
        *d++      = *orig_line++ - *(ref_line += 4);
        *d++      = *orig_line++ - *(ref_line += 4);
        *d        = *orig_line   - *(ref_line += 4);
        
        if (!test8x8transform)
        {
          if ((mcost += SATD (diff, input->hadamard)) >= min_mcost)
          {
            abort_search = 1;
            break;
          }
        }
        else
        {
          i = (x0&0x7) + (x0>7) * 64 + y_offset;
          for(k=0, j=y0; j<y0 + BLOCK_SIZE; j++, k+=BLOCK_SIZE)
            memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], BLOCK_SIZE*sizeof(int));
        }
      }
    }
    
    if(test8x8transform)
      mcost += find_SATD (c_diff, blocktype);
    
    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      best_pos  = pos;
    }
    if (min_mcost < (SubPelThreshold3>>block_type_shift_factor[blocktype])) 
    {
      break;
    }
  }
  
  if (best_pos) 
  {
    *mv_x += spiral_search_x [best_pos];
    *mv_y += spiral_search_y [best_pos];
  }
  
  
  
  
  
#ifdef EIGHTH_PEL
  /***********************************
  *****                         *****
  *****  EIGHTH-PEL REFINEMENT  *****
  *****                         *****
  ***********************************/
  if (img->mv_res)
  {
    //===== convert search center to eighth-pel units (double resolution) =====
    *mv_x <<= 1;
    *mv_y <<= 1;
    pic4_pix_x <<= 1;
    pic4_pix_y <<= 1;
    
    //===== set function for getting pixel values =====
    if ((pic4_pix_x + *mv_x > 1) && (pic4_pix_x + *mv_x < 2*max_pos_x4 - 2) &&
      (pic4_pix_y + *mv_y > 1) && (pic4_pix_y + *mv_y < 2*max_pos_y4 - 2))
      get_line = FastLine4X8;
    else
      get_line = UMVLine4X8;
    
    //===== loop over search positions =====
    for (best_pos = 0, pos = 1; pos < search_pos8; pos++)
    {
      cand_mv_x = *mv_x + spiral_search_x[pos];    // eighth-pel units
      cand_mv_y = *mv_y + spiral_search_y[pos];    // eighth-pel units
      
      //----- set motion vector cost -----
      
      mcost = MV_COST (lambda_factor, 0, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);
      
      if (mcost >= min_mcost) continue;
      cmv_x = cand_mv_x + pic4_pix_x;
      cmv_y = cand_mv_y + pic4_pix_y;
      
      //----- add up SATD -----
      for (y0=0, abort_search=0; y0<blocksize_y && !abort_search; y0+=BLOCK_SIZE)
      {
        y_offset = (y0>7)*ypels;
        ry0 = (y0<<3) + cmv_y;
        ry4  = ry0 + 8;
        ry8  = ry4 + 8;
        ry12 = ry8 + 8;
        y1 = y0 + 1;
        y2 = y1 + 1;
        y3 = y2 + 1;
        
        for (x0=0; x0<blocksize_x; x0+=BLOCK_SIZE)
        {
          rx0 = (x0<<3) + cmv_x;
          d   = diff;
          
          orig_line = &orig_pic [y0][x0];
          ref_line  = get_line (ref_pic, ry0, rx0, img_height, img_width);
          *d++      = *orig_line++ - *(ref_line      );
          *d++      = *orig_line++ - *(ref_line += 4 );
          *d++      = *orig_line++ - *(ref_line += 4 );
          *d++      = *orig_line   - *(ref_line += 4 );
          
          orig_line = &orig_pic [y1][x0];
          ref_line  = get_line (ref_pic, ry4, rx0, img_height, img_width);
          *d++      = *orig_line++ - *(ref_line      );
          *d++      = *orig_line++ - *(ref_line += 4 );
          *d++      = *orig_line++ - *(ref_line += 4 );
          *d++      = *orig_line   - *(ref_line += 4 );
          
          orig_line = &orig_pic [y2][x0];
          ref_line  = get_line (ref_pic, ry8, rx0, img_height, img_width);
          *d++      = *orig_line++ - *(ref_line      );
          *d++      = *orig_line++ - *(ref_line += 4 );
          *d++      = *orig_line++ - *(ref_line += 4 );
          *d++      = *orig_line   - *(ref_line += 4 );
          
          orig_line = &orig_pic [y3][x0];
          ref_line  = get_line (ref_pic, ry12, rx0, img_height, img_width);
          *d++      = *orig_line++ - *(ref_line      );
          *d++      = *orig_line++ - *(ref_line += 4 );
          *d++      = *orig_line++ - *(ref_line += 4 );
          *d        = *orig_line   - *(ref_line += 4 );
          
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
            i = (x0&0x7) +  (x0>7) * 64 + y_offset;
            for(k=0, j=y0; j<BLOCK_SIZE + y0; j++, k+=BLOCK_SIZE)
              memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], BLOCK_SIZE*sizeof(int));
          }
        }
      }
      if (test8x8transform)
        mcost += find_SATD (c_diff, blocktype);
      
      if (mcost < min_mcost)
      {
        min_mcost = mcost;
        best_pos  = pos;
      }
    }
    if (best_pos)
    {
      *mv_x += spiral_search_x [best_pos];
      *mv_y += spiral_search_y [best_pos];
    }
  }
  
#endif
  
  
  
  
  
  
  
  //===== return minimum motion cost =====
  return min_mcost;
}


/*!
***********************************************************************
* \brief
*    Full pixel block motion search
***********************************************************************
*/
int                                                //  ==> minimum motion cost after search
FullPelBlockMotionBiPred (pel_t**   orig_pic,      // <--  original pixel values for the AxB block
                          short       ref,         // <--  reference frame (0... or -1 (backward))
                          int       list,
                          int       pic_pix_x,     // <--  absolute x-coordinate of regarded AxB block
                          int       pic_pix_y,     // <--  absolute y-coordinate of regarded AxB block
                          int       blocktype,     // <--  block type (1-16x16 ... 7-4x4)
                          short     pred_mv_x1,    // <--  motion vector predictor (x) in sub-pel units
                          short     pred_mv_y1,    // <--  motion vector predictor (y) in sub-pel units
                          short     pred_mv_x2,    // <--  motion vector predictor (x) in sub-pel units
                          short     pred_mv_y2,    // <--  motion vector predictor (y) in sub-pel units
                          short*    mv_x,          // <--> in: search center (x) / out: motion vector (x) - in pel units
                          short*    mv_y,          // <--> in: search center (y) / out: motion vector (y) - in pel units
                          short*    s_mv_x,        // <--> in: search center (x) / out: motion vector (x) - in pel units
                          short*    s_mv_y,        // <--> in: search center (y) / out: motion vector (y) - in pel units
                          int       search_range,  // <--  1-d search range in pel units
                          int       min_mcost,     // <--  minimum motion cost (cost for center or huge value)
                          int       lambda_factor) // <--  lagrangian parameter for determining motion cost
{
  int   pos, cand_x, cand_y, y, x4, mcost;
  
  pel_t *orig_line, *ref2_line, *ref1_line;
  pel_t *(*get_ref_line1)(int, pel_t*, int, int, int, int);
  pel_t *(*get_ref_line2)(int, pel_t*, int, int, int, int);
  
  int   list_offset   = img->mb_data[img->current_mb_nr].list_offset; 
  pel_t *ref1_pic     = listX[list + list_offset          ][ref]->imgY_11;
  pel_t *ref2_pic     = listX[list ==0? 1 + list_offset: list_offset][ 0 ]->imgY_11;  
  int   img_width     = listX[list+list_offset            ][ref]->size_x;
  int   img_height    = listX[list+list_offset            ][ref]->size_y;
  
  int   best_pos      = 0;                                        // position with minimum motion cost
  int   max_pos       = (2*search_range+1)*(2*search_range+1);    // number of search positions
  int   blocksize_y   = input->blc_size[blocktype][1];            // vertical block size
  int   blocksize_x   = input->blc_size[blocktype][0];            // horizontal block size
  int   blocksize_x4  = blocksize_x >> 2;                         // horizontal block size in 4-pel units
  int   pred_x1        = (pic_pix_x << 2) + pred_mv_x1;       // predicted position x (in sub-pel units)
  int   pred_y1        = (pic_pix_y << 2) + pred_mv_y1;       // predicted position y (in sub-pel units)
  int   pred_x2        = (pic_pix_x << 2) + pred_mv_x2;       // predicted position x (in sub-pel units)
  int   pred_y2        = (pic_pix_y << 2) + pred_mv_y2;       // predicted position y (in sub-pel units)
  short center_x      = pic_pix_x + *mv_x;                      // center position x (in pel units)
  short center_y      = pic_pix_y + *mv_y;                      // center position y (in pel units)
  short ref1_center_x = pic_pix_x + *s_mv_x;                      // mvx of second pred (in pel units)
  short ref1_center_y = pic_pix_y + *s_mv_y;                      // mvy of second pred (in pel units)
  
  
  int   bi_diff; 
  short apply_weights   = (active_pps->weighted_bipred_idc>0);  
  short weightSpic = (apply_weights ? (list == 0? wbp_weight[list_offset    ][ref][0  ][0]: wbp_weight[list_offset + 1][0  ][ref][0]) : 1<<luma_log_weight_denom);
  short weightRpic = (apply_weights ? (list == 0? wbp_weight[list_offset + 1][ref][0  ][0]: wbp_weight[list_offset    ][0  ][ref][0]) : 1<<luma_log_weight_denom);
  short offsetSpic = (apply_weights ? (list == 0?  wp_offset[list_offset    ][ref]     [0]:  wp_offset[list_offset + 1][0  ]     [0]) : 0);
  short offsetRpic = (apply_weights ? (list == 0?  wp_offset[list_offset + 1][ref]     [0]:  wp_offset[list_offset    ][0  ]     [0]) : 0);
  short weightedpel,pixel1,pixel2;
#ifdef  USE_HP_FILTER//BROUND
  short offsetBi;
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
#else
  short offsetBi=(offsetRpic + offsetSpic + 1)>>1;
#endif
  //===== set function for getting reference picture lines =====
  if ((center_x > search_range) && (center_x < img_width -1-search_range-blocksize_x) &&
    (center_y > search_range) && (center_y < img_height-1-search_range-blocksize_y)   )
  {
    get_ref_line2 = FastLineX;
  }
  else
  {
    get_ref_line2 = UMVLineX2;
  }
  
  //===== set function for getting reference picture lines =====
  if ((ref1_center_x > search_range) && (ref1_center_x < img_width -1-search_range-blocksize_x) &&
    (ref1_center_y > search_range) && (ref1_center_y < img_height-1-search_range-blocksize_y)   )
  {
    get_ref_line1 = FastLineX;
  }
  else
  {
    get_ref_line1 = UMVLineX;
  }
  
  //===== loop over all search positions =====
  for (pos=0; pos<max_pos; pos++)
  {
    //--- set candidate position (absolute position in pel units) ---
    cand_x = center_x + spiral_search_x[pos];
    cand_y = center_y + spiral_search_y[pos];
    
    //--- initialize motion cost (cost for motion vector) and check ---
    
    mcost =  MV_COST (lambda_factor, 2, ref1_center_x, ref1_center_y, pred_x1, pred_y1);
    mcost += MV_COST (lambda_factor, 2,        cand_x,        cand_y, pred_x2, pred_y2);
    
    if (mcost >= min_mcost)   continue;
    
    //--- add residual cost to motion cost ---
    if (apply_weights)
    {
      for (y=0; y<blocksize_y; y++)
      {
        ref2_line  = get_ref_line2 (blocksize_x, ref2_pic,        cand_y+y,        cand_x, img_height, img_width);
        ref1_line  = get_ref_line1 (blocksize_x, ref1_pic, ref1_center_y+y, ref1_center_x, img_height, img_width);
        orig_line = orig_pic [y];
        
        for (x4=0; x4<blocksize_x4; x4++)
        { 
          pixel1=weightSpic * (*ref1_line++);
          pixel2=weightRpic * (*ref2_line++);
#ifdef  USE_HP_FILTER//BROUND
          if(input->UseHPFilter != 0)
          {
            weightedpel =  Clip3 (0, img->max_imgpel_value ,((pixel1 + pixel2 +
                                (offsetBi * 2 * wp_luma_round) +
                                2 * wp_luma_round * round) >> (luma_log_weight_denom + 1)));
          }
          else
          {
            weightedpel =  Clip3 (0, img->max_imgpel_value ,((pixel1 + pixel2 +
                                2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          }
#else
          weightedpel =  Clip3 (0, img->max_imgpel_value ,((pixel1 + pixel2 +
                                2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
#endif
          bi_diff = (*orig_line++)  - weightedpel;
          mcost += byte_abs[bi_diff];
          
          pixel1=weightSpic * (*ref1_line++);
          pixel2=weightRpic * (*ref2_line++);
#ifdef  USE_HP_FILTER//BROUND
          if(input->UseHPFilter != 0)
          {
            weightedpel =  Clip3 (0, img->max_imgpel_value ,((pixel1 + pixel2 +
                                (offsetBi * 2 * wp_luma_round) +
                                2 * wp_luma_round  * round) >> (luma_log_weight_denom + 1)));
          }
          else
          {
            weightedpel =  Clip3 (0, img->max_imgpel_value ,((pixel1 + pixel2 +
                                2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          }
#else
          weightedpel =  Clip3 (0, img->max_imgpel_value ,((pixel1 + pixel2 +
                                2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
#endif
          bi_diff = (*orig_line++)  - weightedpel;
          mcost += byte_abs[bi_diff];
          
          pixel1=weightSpic * (*ref1_line++);
          pixel2=weightRpic * (*ref2_line++);
#ifdef  USE_HP_FILTER//BROUND
          if(input->UseHPFilter != 0)
          {
            weightedpel =  Clip3 (0, img->max_imgpel_value ,((pixel1 + pixel2 +
                                (offsetBi * 2 * wp_luma_round) +
                                2 * wp_luma_round  * round) >> (luma_log_weight_denom + 1)));
          }
          else
          {
            weightedpel =  Clip3 (0, img->max_imgpel_value ,((pixel1 + pixel2 +
                                2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          }
#else
          weightedpel =  Clip3 (0, img->max_imgpel_value ,((pixel1 + pixel2 +
                                2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
#endif
          bi_diff = (*orig_line++)  - weightedpel;
          mcost += byte_abs[bi_diff];
          
          pixel1=weightSpic * (*ref1_line++);
          pixel2=weightRpic * (*ref2_line++);
#ifdef  USE_HP_FILTER//BROUND
          if(input->UseHPFilter != 0)
          {
            weightedpel =  Clip3 (0, img->max_imgpel_value ,((pixel1 + pixel2 +
                                (offsetBi * 2 * wp_luma_round) +
                                2 * wp_luma_round  * round) >> (luma_log_weight_denom + 1)));
          }
          else
          {
            weightedpel =  Clip3 (0, img->max_imgpel_value ,((pixel1 + pixel2 +
                                2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
          }

#else
          weightedpel =  Clip3 (0, img->max_imgpel_value ,((pixel1 + pixel2 +
                                2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
#endif
          bi_diff = (*orig_line++)  - weightedpel;
          mcost += byte_abs[bi_diff];
        }    
        
        if (mcost >= min_mcost)
        {
          break;
        }        
      }
    }
    else
    {
      for (y=0; y<blocksize_y; y++)
      {
        ref2_line = get_ref_line2 (blocksize_x, ref2_pic,        cand_y+y,        cand_x, img_height, img_width);
        ref1_line = get_ref_line1 (blocksize_x, ref1_pic, ref1_center_y+y, ref1_center_x, img_height, img_width);
        orig_line = orig_pic [y];
        
        for (x4=0; x4<blocksize_x4; x4++)
        {         
#ifdef  USE_HP_FILTER//BROUND
          if(input->UseHPFilter != 0)
            round = 1 - img->bipred_rounding_control;
          else
            round =0;

          bi_diff = (*orig_line++) - (((*ref1_line++) + (*ref2_line++) + round)>>1) ;
          mcost += byte_abs[bi_diff];
          bi_diff = (*orig_line++) - (((*ref1_line++) + (*ref2_line++) + round)>>1) ;
          mcost += byte_abs[bi_diff];
          bi_diff = (*orig_line++) - (((*ref1_line++) + (*ref2_line++) + round)>>1) ;
          mcost += byte_abs[bi_diff];
          bi_diff = (*orig_line++) - (((*ref1_line++) + (*ref2_line++) + round)>>1) ;
          mcost += byte_abs[bi_diff];
#else
          bi_diff = (*orig_line++) - (((*ref1_line++) + *ref2_line++)>>1) ;
          mcost += byte_abs[bi_diff];
          bi_diff = (*orig_line++) - (((*ref1_line++) + *ref2_line++)>>1) ;
          mcost += byte_abs[bi_diff];
          bi_diff = (*orig_line++) - (((*ref1_line++) + *ref2_line++)>>1) ;
          mcost += byte_abs[bi_diff];
          bi_diff = (*orig_line++) - (((*ref1_line++) + *ref2_line++)>>1) ;
          mcost += byte_abs[bi_diff];
#endif
        }    
        
        if (mcost >= min_mcost)
        {
          break;
        }
      }
    }
    
    //--- check if motion cost is less than minimum cost ---
    if (mcost < min_mcost)
    {
      best_pos  = pos;
      min_mcost = mcost;
    }
  }
  
  
  //===== set best motion vector and return minimum motion cost =====
  if (best_pos)
  {
    *mv_x += spiral_search_x[best_pos];
    *mv_y += spiral_search_y[best_pos];
  }
  return min_mcost;
}

/*!
***********************************************************************
* \brief
*    Bipred Sub pixel block motion search
***********************************************************************
*/
int                                               //  ==> minimum motion cost after search
SubPelBlockSearchBiPred (pel_t**   orig_pic,      // <--  original pixel values for the AxB block
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
                         int       lambda_factor   // <--  lagrangian parameter for determining motion cost
                         )
{
  int   j, i, k;
  int   c_diff[MB_PIXELS];
  int   diff[16], *d;  
  int   pos, best_pos, mcost, abort_search;
  int   y0, y1, y2, y3;
  int   x0;
  int   ry0, ry4, ry8, ry12, rx0;
  int   sy0, sy4, sy8, sy12, sx0;
  
  int   cand_mv_x, cand_mv_y;
  int   y_offset, ypels =(128 - 64 * (blocktype == 3));
  
#ifdef SWITCHED_FILTERS
  int   best_pos_sp16 = 0;
  int   min_mcost_sp16 = min_mcost;
  short mv_x_sp16 = 0;
  short mv_y_sp16 = 0;
#endif

  int   blocksize_x     = input->blc_size[blocktype][0];
  int   blocksize_y     = input->blc_size[blocktype][1];
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
#ifdef MV_COMPETITION
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
      y_offset = (y0>7)*ypels;
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
            i = (x0&0x7) + (x0>7) * 64 + y_offset;
            for(k=0, j=y0; j<y0 + BLOCK_SIZE; j++, k+=BLOCK_SIZE)
              memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], BLOCK_SIZE*sizeof(int));
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
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1     ) + *(ref_line_p2     ) + round) >> 1);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line    -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          
          orig_line = &orig_pic [y1][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy4, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry4, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1     ) + *(ref_line_p2     ) + round) >> 1);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line    -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          
          orig_line = &orig_pic [y2][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy8, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry8, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1     ) + *(ref_line_p2     ) + round) >> 1);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line    -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          
          orig_line = &orig_pic [y3][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy12, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry12, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1     ) + *(ref_line_p2     ) + round) >> 1);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d   = *orig_line    -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
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
            i = (x0&0x7) + (x0>7) * 64 + y_offset;
            for(k=0, j=y0; j<BLOCK_SIZE + y0; j++, k+=BLOCK_SIZE)
              memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], BLOCK_SIZE*sizeof(int));
          }
        }
      }
    }  
    
    if(test8x8transform)
      mcost += find_SATD (c_diff, blocktype);
    
#ifdef SWITCHED_FILTERS
    if((img->filterParam != SIFO_FIRST_PASS_FPO) || (pos < max_pos2))
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
#ifdef MV_COMPETITION
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
      y_offset = (y0>7)*ypels;
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
            i = (x0&0x7) + (x0>7) * 64 + y_offset;
            for(k=0, j=y0; j<y0 + BLOCK_SIZE; j++, k+=BLOCK_SIZE)
              memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], BLOCK_SIZE*sizeof(int));
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
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1     ) + *(ref_line_p2     ) + round) >> 1);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line    -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          
          orig_line = &orig_pic [y1][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy4, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry4, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1     ) + *(ref_line_p2     ) + round) >> 1);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line    -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          
          orig_line = &orig_pic [y2][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy8, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry8, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1     ) + *(ref_line_p2     ) + round) >> 1);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line    -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          
          orig_line = &orig_pic [y3][x0];    
          ref_line_p1 = get_line_p1 (ref1_pic, sy12, sx0,  img_height, img_width);
          ref_line_p2 = get_line_p2 (ref2_pic, ry12, rx0,  img_height, img_width);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1     ) + *(ref_line_p2     ) + round) >> 1);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d++ = *orig_line++  -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
          *d   = *orig_line    -  clip1a((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + round) >> 1);
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
            i = (x0&0x7) + (x0>7) * 64 + y_offset;
            for(k=0, j=y0; j<y0 + BLOCK_SIZE; j++, k+=BLOCK_SIZE)
              memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], BLOCK_SIZE*sizeof(int));
          }
        }
      }
    }
    if(test8x8transform)
      mcost += find_SATD (c_diff, blocktype);
    
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

#ifdef EIGHTH_PEL
  /************************************
  *****                          *****
  *****  EIGHTH-PEL REFINEMENT  *****
  *****                          *****
  ************************************/
  if (img->mv_res)
  {
    //===== convert search center to eighth-pel units (double resolution) =====
    *mv_x <<= 1;
    *mv_y <<= 1;
    smv_x <<= 1;
    smv_y <<= 1;
    *s_mv_x <<= 1;
    *s_mv_y <<= 1;
    pic4_pix_x <<= 1;
    pic4_pix_y <<= 1;
    
    //===== set function for getting pixel values =====
    if ((pic4_pix_x + *mv_x > 1) && (pic4_pix_x + *mv_x < 2*max_pos_x4 - 2) &&
      (pic4_pix_y + *mv_y > 1) && (pic4_pix_y + *mv_y < 2*max_pos_y4 - 2))
      get_line_p2 = FastLine4X8;
    else
      get_line_p2 = UMVLine4X8;
    
    if ((pic4_pix_x + *s_mv_x > 1) && (pic4_pix_x + *s_mv_x < 2*max_pos_x4 - 2) &&
      (pic4_pix_y + *s_mv_y > 1) && (pic4_pix_y + *s_mv_y < 2*max_pos_y4 - 2))
      get_line_p1 = FastLine4X8;
    else
      get_line_p1 = UMVLine4X8;
    
    //===== loop over search positions =====
    for (best_pos = 0, pos = 1; pos < search_pos8; pos++)
    {
      cand_mv_x = *mv_x + spiral_search_x[pos];    // quarter-pel units
      cand_mv_y = *mv_y + spiral_search_y[pos];    // quarter-pel units
      
      //----- set motion vector cost -----
#ifdef MV_COMPETITION
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
        y_offset = (y0>7)*ypels;
        ry0 = (y0<<2) + cmv_y;
        ry4  = ry0 + 8;
        ry8  = ry4 + 8;
        ry12 = ry8 + 8;
        sy0 = (y0<<2) + smv_y;
        sy4  = sy0 + 8;
        sy8  = sy4 + 8;
        sy12 = sy8 + 8;
        y1 = y0 + 1;
        y2 = y1 + 1;
        y3 = y2 + 1;
        
        if (apply_weights)
        {
          
          for (x0=0; x0<blocksize_x; x0+=BLOCK_SIZE)
          {
            rx0 = (x0<<3) + cmv_x;
            sx0 = (x0<<3) + smv_x;
            d   = diff;
            
            orig_line = &orig_pic [y0][x0];
            ref_line_p1 = get_line_p1 (ref1_pic, sy0, sx0,  img_height, img_width);
            ref_line_p2 = get_line_p2 (ref2_pic, ry0, rx0,  img_height, img_width);
            weightedpel = clip1a ((((weightSpic * (*(ref_line_p1++ )) + weightRpic * (*(ref_line_p2++ )))
              + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
            *d++ = *orig_line++ - weightedpel;
            weightedpel = clip1a ((((weightSpic * (*(ref_line_p1++ )) + weightRpic * (*(ref_line_p2++ )))
              + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
            *d++ = *orig_line++ - weightedpel;
            weightedpel = clip1a ((((weightSpic * (*(ref_line_p1++ )) + weightRpic * (*(ref_line_p2++ )))
              + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
            *d++ = *orig_line++ - weightedpel;
            weightedpel = clip1a ((((weightSpic * (*(ref_line_p1++ )) + weightRpic * (*(ref_line_p2++ )))
              + 2 * wp_luma_round) >> (luma_log_weight_denom + 1)) + (offsetBi));
            *d++ = *orig_line - weightedpel;
            
            orig_line = &orig_pic [y1][x0];
            ref_line_p1 = get_line_p1 (ref1_pic, sy4, sx0,  img_height, img_width);
            ref_line_p2 = get_line_p2 (ref2_pic, ry4, rx0,  img_height, img_width);
            weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4)))
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
            weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4)))
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
            weightedpel = clip1a ((((weightSpic * (*(ref_line_p1 += 4)) + weightRpic * (*(ref_line_p2 += 4)))
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
              i = (x0&0x7) + (x0>7) * 64 + y_offset;
              for(k=0, j=y0; j<y0 + BLOCK_SIZE; j++, k+=BLOCK_SIZE)
                memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], BLOCK_SIZE*sizeof(int));
            }
          }
        }
        else
        {
          for (x0=0; x0<blocksize_x; x0+=BLOCK_SIZE)
          {
            rx0 = (x0<<3) + cmv_x;
            sx0 = (x0<<3) + smv_x;
            d   = diff;
            
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
            *d++ = *orig_line    -  ((*(ref_line_p1 += 4) + *(ref_line_p2 += 4) + 1) >> 1);
            
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
              i = (x0&0x7) + (x0>7) * 64 + y_offset;
              for(k=0, j=y0; j<y0 + BLOCK_SIZE; j++, k+=BLOCK_SIZE)
                memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], BLOCK_SIZE*sizeof(int));
            }
          }
        }
      }
      if(test8x8transform)
        mcost += find_SATD (c_diff, blocktype);
      
      if (mcost < min_mcost)
      {
        min_mcost = mcost;
        best_pos  = pos;
      }
    }
    if (best_pos)
    {
      *mv_x += spiral_search_x [best_pos];
      *mv_y += spiral_search_y [best_pos];
    }
  }
#endif
  
  //===== return minimum motion cost =====
  return min_mcost;
}

/*!
***********************************************************************
* \brief
*    Motion Cost for Bidirectional modes
***********************************************************************
*/
int BPredPartitionCost (int   blocktype,
                        int   block8x8,
                        short fw_ref,
                        short bw_ref,
                        int   lambda_factor,
                        int   list)
{
  static int  bx0[5][4] = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,2,0,0}, {0,2,0,2}};
  static int  by0[5][4] = {{0,0,0,0}, {0,0,0,0}, {0,2,0,0}, {0,0,0,0}, {0,0,2,2}};
  
  int   diff[64];
  int   curr_blk[MB_BLOCK_SIZE][MB_BLOCK_SIZE]; // ABT pred.error buffer
  int   bsx       = min(input->blc_size[blocktype][0],8);
  int   bsy       = min(input->blc_size[blocktype][1],8);
  
  int   pic_pix_x, pic_pix_y, block_x, block_y;
  int   v, h, mcost, i, j, k;
  int   mvd_bits  = 0;
  int   parttype  = (blocktype<4?blocktype:4);
  int   step_h0   = (input->part_size[ parttype][0]);
  int   step_v0   = (input->part_size[ parttype][1]);
  int   step_h    = (input->part_size[blocktype][0]);
  int   step_v    = (input->part_size[blocktype][1]);
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
      if ((!input->Transform8x8Mode) || (blocktype>4))
      {
        mcost += SATD (diff, input->hadamard);          
      }
    }
  }
  if (input->Transform8x8Mode && (blocktype<=4))  // tchen 4-29-04
  {
    for (byy=0; byy < input->blc_size[parttype][1]; byy+=bsy)
      for (bxx=0; bxx<input->blc_size[parttype][0]; bxx+=bsx)
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
*    Block motion search
***********************************************************************
*/
int                                         //!< minimum motion cost after search
BlockMotionSearch (short     ref,           //!< reference idx
                   int       list,          //!< reference pciture list
                   int       mb_x,          //!< x-coordinate inside macroblock
                   int       mb_y,          //!< y-coordinate inside macroblock
                   int       blocktype,     //!< block type (1-16x16 ... 7-4x4)
                   int       search_range,  //!< 1-d search range for integer-position search
                   int       lambda_factor) //!< lagrangian parameter for determining motion cost
{
  static pel_t   orig_val [256];
  static pel_t  *orig_pic  [16] = {orig_val,     orig_val+ 16, orig_val+ 32, orig_val+ 48,
    orig_val+ 64, orig_val+ 80, orig_val+ 96, orig_val+112,
    orig_val+128, orig_val+144, orig_val+160, orig_val+176,
    orig_val+192, orig_val+208, orig_val+224, orig_val+240};
  
  short     pred_mv_x, pred_mv_y, mv_x, mv_y;
  int       i, j;
#ifdef EIGHTH_PEL
  int       mvshift   = (img->mv_res ? 3 : 2);
#endif
  int       max_value = INT_MAX;
  int       min_mcost = max_value;
  
  int       block_x   = (mb_x>>2);
  int       block_y   = (mb_y>>2);
  
  int       bsx       = input->blc_size[blocktype][0];
  int       bsy       = input->blc_size[blocktype][1];
  
  int       pic_pix_x = img->opix_x + mb_x;
  int       pic_pix_y = img->opix_y + mb_y;
  
  
#ifdef WIN32
  struct _timeb tstruct1;
  struct _timeb tstruct2;
#else
  struct timeb tstruct1;
  struct timeb tstruct2;
#endif
  
  int me_tmp_time;
  short*    pred_mv = img->pred_mv[block_y][block_x][list][ref][blocktype];
  short****** all_mv    = img->all_mv;  
  int list_offset = ((img->MbaffFrameFlag) && (img->mb_data[img->current_mb_nr].mb_field)) ? img->current_mb_nr % 2 ? 4 : 2 : 0;
  int *prevSad = (input->FMEnable == 3)? EPZSDistortion[list + list_offset][blocktype - 1]: NULL;
  
#ifdef SWITCHED_FILTERS
  static int min_mcost_offset;
#endif

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
  
  if(input->FMEnable == 1)
  {
    FME_blocktype = blocktype;
    bipred_flag = 0;
  }
  else if (input->FMEnable == 2)
  {
    simplified_setup_FME(ref, list, block_y, block_x, blocktype, all_mv );
  }
  
  //===========================================
  //=====   GET MOTION VECTOR PREDICTOR   =====
  //===========================================
#ifdef MV_COMPETITION
  if (input->mv_competition > 0)
    SetMotionVectorPredictor_Competition (pred_mv, enc_picture->ref_idx[list], enc_picture->mv[list], ref, list, block_x, block_y, bsx, bsy);
  else
#endif
    SetMotionVectorPredictor (pred_mv, enc_picture->ref_idx[list], enc_picture->mv[list], ref, list, block_x, block_y, bsx, bsy);
  
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
  
  pred_mv_x = pred_mv[0];
  pred_mv_y = pred_mv[1];
  
  //==================================
  //=====   INTEGER-PEL SEARCH   =====
  //==================================
  
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
    
    min_mcost = EPZSPelBlockMotionSearch (orig_pic, ref, list, list_offset, 
      enc_picture->ref_idx, enc_picture->mv,pic_pix_x, pic_pix_y, blocktype,
      pred_mv_x, pred_mv_y, &mv_x, &mv_y, search_range, min_mcost, lambda_factor);
    
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
      {
        if (input->EPZSSubPelME)
#ifdef EIGHTH_PEL
          min_mcost =  EPZSSubPelBlockMotionSearch (orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
          pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9, 9, min_mcost, lambda_factor);            
#else
        min_mcost =  EPZSSubPelBlockMotionSearch (orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
          pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9, min_mcost, lambda_factor);            
#endif
        else
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
  
  //===============================================
  //=====   SET MV'S AND RETURN MOTION COST   =====
  //===============================================
  
#ifdef MV_COMPETITION
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
  for(j = block_y; j < block_y + (bsy>>2); j++)    
  {
    for(i = block_x; i < block_x + (bsx>>2); i++)
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
#ifdef MV_COMPETITION
    if (input->mv_competition > 0)
      SetMotionVectorPredictor_Competition (pred_mv_bi, enc_picture->ref_idx[list ^ 1], enc_picture->mv[(list == LIST_0? LIST_1: LIST_0)], 0, (list == LIST_0? LIST_1: LIST_0), block_x, block_y, bsx, bsy);
    else
#endif
      SetMotionVectorPredictor (pred_mv_bi, enc_picture->ref_idx[list ^ 1], enc_picture->mv[(list == LIST_0? LIST_1: LIST_0)], 0, (list == LIST_0? LIST_1: LIST_0), block_x, block_y, bsx, bsy);
    
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
        min_mcostbi = EPZSBiPredBlockMotionSearch (orig_pic, ref, iterlist, list_offset, 
          enc_picture->ref_idx, enc_picture->mv, 
          pic_pix_x, pic_pix_y, blocktype, 
          pred_mv_x1, pred_mv_y1, pred_mv_x2, pred_mv_y2, 
          &bimv_x, &bimv_y, &tempmv_x, &tempmv_y, 
          input->BiPredMESearchRange, min_mcostbi, lambda_factor);
      }
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
      min_mcostbi =  SubPelBlockSearchBiPred (orig_pic, 0, iterlist, pic_pix_x, pic_pix_y, blocktype,
        pred_mv_x2, pred_mv_y2, &bimv_x, &bimv_y, &mv_x, &mv_y, 9, 9, 9,
        min_mcostbi, lambda_factor);
#else
      min_mcostbi =  SubPelBlockSearchBiPred (orig_pic, 0, iterlist, pic_pix_x, pic_pix_y, blocktype,
        pred_mv_x2, pred_mv_y2, &bimv_x, &bimv_y, &mv_x, &mv_y, 9, 9,
        min_mcostbi, lambda_factor);
#endif
    }

    if (input->BiPredMESubPel==2 && !input->DisableSubpelME)
    {
      if (input->hadamard)
      {
        min_mcostbi = max_value;
      }

#ifdef EIGHTH_PEL
      min_mcostbi =  SubPelBlockSearchBiPred (orig_pic, 0, (iterlist == LIST_0? LIST_1: LIST_0), pic_pix_x, pic_pix_y, blocktype,
        pred_mv_x, pred_mv_y, &mv_x, &mv_y, &bimv_x, &bimv_y, 9, 9, 9,
        min_mcostbi, lambda_factor);
#else
      min_mcostbi =  SubPelBlockSearchBiPred (orig_pic, 0, (iterlist == LIST_0? LIST_1: LIST_0), pic_pix_x, pic_pix_y, blocktype,
        pred_mv_x, pred_mv_y, &mv_x, &mv_y, &bimv_x, &bimv_y, 9, 9,
        min_mcostbi, lambda_factor);
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
              min_mcostbi = EPZSBiPredBlockMotionSearch (orig_pic, ref, iterlist, list_offset, 
                enc_picture->ref_idx, enc_picture->mv, 
                pic_pix_x, pic_pix_y, blocktype, 
                pred_mv_x1, pred_mv_y1, pred_mv_x2, pred_mv_y2, 
                &bimv_x, &bimv_y, &tempmv_x, &tempmv_y, 
                input->BiPredMESearchRange, min_mcostbi, lambda_factor);
            }
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
            min_mcostbi =  SubPelBlockSearchBiPred (orig_pic, 0, iterlist, pic_pix_x, pic_pix_y, blocktype,
              pred_mv_x2, pred_mv_y2, &bimv_x, &bimv_y, &mv_x, &mv_y, 9, 9, 9,
              min_mcostbi, lambda_factor);
#else
            min_mcostbi =  SubPelBlockSearchBiPred (orig_pic, 0, iterlist, pic_pix_x, pic_pix_y, blocktype,
              pred_mv_x2, pred_mv_y2, &bimv_x, &bimv_y, &mv_x, &mv_y, 9, 9,
              min_mcostbi, lambda_factor);
#endif
          }

          if (input->BiPredMESubPel==2 && !input->DisableSubpelME)
          {
            if (input->hadamard)
            {
              min_mcostbi = max_value;
            }

#ifdef EIGHTH_PEL
            min_mcostbi =  SubPelBlockSearchBiPred (orig_pic, 0, (iterlist == LIST_0? LIST_1: LIST_0), pic_pix_x, pic_pix_y, blocktype,
              pred_mv_x, pred_mv_y, &mv_x, &mv_y, &bimv_x, &bimv_y, 9, 9, 9,
              min_mcostbi, lambda_factor);
#else
            min_mcostbi =  SubPelBlockSearchBiPred (orig_pic, 0, (iterlist == LIST_0? LIST_1: LIST_0), pic_pix_x, pic_pix_y, blocktype,
              pred_mv_x, pred_mv_y, &mv_x, &mv_y, &bimv_x, &bimv_y, 9, 9,
              min_mcostbi, lambda_factor);
#endif
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
***********************************************************************
* \brief
*    Motion Cost for Bidirectional modes
***********************************************************************
*/
int BIDPartitionCost (int   blocktype,
                      int   block8x8,
                      short fw_ref,
                      short bw_ref,
                      int   lambda_factor)
{
  static int  bx0[5][4] = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,2,0,0}, {0,2,0,2}};
  static int  by0[5][4] = {{0,0,0,0}, {0,0,0,0}, {0,2,0,0}, {0,0,0,0}, {0,0,2,2}};
  
  int   diff[64];
  int   curr_blk[MB_BLOCK_SIZE][MB_BLOCK_SIZE]; // ABT pred.error buffer
  int   bsx       = min(input->blc_size[blocktype][0],8);
  int   bsy       = min(input->blc_size[blocktype][1],8);
  
  int   pic_pix_x, pic_pix_y, block_x, block_y;
  int   v, h, mcost, i, j, k;
  int   mvd_bits  = 0;
  int   parttype  = (blocktype<4?blocktype:4);
  int   step_h0   = (input->part_size[ parttype][0]);
  int   step_v0   = (input->part_size[ parttype][1]);
  int   step_h    = (input->part_size[blocktype][0]);
  int   step_v    = (input->part_size[blocktype][1]);
  int   bxx, byy;                               // indexing curr_blk
  int   bx = bx0[parttype][block8x8];
  int   by = by0[parttype][block8x8];
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
    for (byy=0; byy < input->blc_size[parttype][1]; byy+=bsy)
      for (bxx=0; bxx<input->blc_size[parttype][0]; bxx+=bsx)
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
*    Get cost for skip mode for an macroblock
************************************************************************
*/
int GetSkipCostMB (int lambda_factor)
{
  int block_y, block_x, pic_pix_y, pic_pix_x, i, j, k;
  int diff[16];
  int cost = 0;
  
  int curr_diff[8][8], diff8x8[64];
  int mb_x, mb_y;
  int block;  
  for(block=0;block<4;block++)
  {
    mb_y    = (block/2)<<3;
    mb_x    = (block%2)<<3;
    for (block_y=mb_y; block_y<mb_y+8; block_y+=4)
    {
      pic_pix_y = img->opix_y + block_y;
      for (block_x=mb_x; block_x<mb_x+8; block_x+=4)
      {
        pic_pix_x = img->opix_x + block_x;
        
        //===== prediction of 4x4 block =====
        LumaPrediction4x4 (block_x, block_y, 0, 0, 0, 0, 0);
        
        //===== get displaced frame difference ======                
        for (k=j=0; j<4; j++)
          for (i=0; i<4; i++, k++)
          {
            diff[k] = curr_diff[block_y-mb_y+j][block_x-mb_x+i] = imgY_org[pic_pix_y+j][pic_pix_x+i] - img->mpr[j+block_y][i+block_x];
          }
          
          if(!((input->rdopt==0)&&(input->Transform8x8Mode)))
            cost += SATD (diff, input->hadamard);
      }
    }
    
    if((input->rdopt==0)&&(input->Transform8x8Mode))
    {      
      for(k=j=0; j<8; j++, k+=8)
        memcpy(&diff8x8[k], &(curr_diff[j]), 8 * sizeof(int));          
      cost += SATD8X8 (diff8x8, input->hadamard);
    }
  }
  
  return cost;
}

/*!
************************************************************************
* \brief
*    Find motion vector for the Skip mode
************************************************************************
*/
void FindSkipModeMotionVector ()
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
  
#ifdef MV_COMPETITION
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
        for (by = 0;by < 4;by++)
          for (bx = 0;bx < 4;bx++)
          {
            memset(all_mv [by][bx][0][0][0], 0, 2* sizeof(short));
          }
          
#ifdef MV_COMPETITION
          if (input->mv_competition > 0)
          {
            mv_comp.mv_pred_skip[prediction_mode_for_skip][0] = 0;
            mv_comp.mv_pred_skip[prediction_mode_for_skip][1] = 0;
          }
#endif
      }
      else
      {
        
        SetMotionVectorPredictor (pmv, enc_picture->ref_idx[LIST_0], enc_picture->mv[LIST_0], 0, LIST_0, 0, 0, 16, 16);
        
        for (by = 0;by < 4;by++)
          for (bx = 0;bx < 4;bx++)
          {
            memcpy(all_mv [by][bx][0][0][0], pmv, 2* sizeof(short));
            //all_mv [by][bx][0][0][0][0] = pmv[0];
            //all_mv [by][bx][0][0][0][1] = pmv[1];
          }
#ifdef MV_COMPETITION
          if (input->mv_competition > 0)
          {
            mv_comp.mv_pred_skip[prediction_mode_for_skip][0] = pmv[0];
            mv_comp.mv_pred_skip[prediction_mode_for_skip][1] = pmv[1];
          }
#endif
      }
#ifdef MV_COMPETITION    
    }
    else   
    {
      PixelPos block_a, block_b, block_c, block_d;
      int mb_nr = img->current_mb_nr;
      int block_x = 0, block_y = 0;
      int blockshape_x = 16;//, blockshape_y = 16;
      int hv, mv_a, mv_b, mv_c;
      int y=(int)img->current_mb_nr/(img->width/16);  // Vertical
      int x=(int)img->current_mb_nr%(img->width/16);  // Horizontal
      //      int ref_frame = 0;
      
      getLuma4x4Neighbour(mb_nr, block_x, block_y,           -1,  0, &block_a);
      getLuma4x4Neighbour(mb_nr, block_x, block_y,            0, -1, &block_b);
      getLuma4x4Neighbour(mb_nr, block_x, block_y, blockshape_x, -1, &block_c);
      getLuma4x4Neighbour(mb_nr, block_x, block_y,           -1, -1, &block_d);
      
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
          for (by = 0;by < 4;by++)
            for (bx = 0;bx < 4;bx++)
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
          for (by = 0;by < 4;by++)
            for (bx = 0;bx < 4;bx++)
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
          for (by = 0;by < 4;by++)
            for (bx = 0;bx < 4;bx++)
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
          for (by = 0;by < 4;by++)
            for (bx = 0;bx < 4;bx++)
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
          for (by = 0;by < 4;by++)
            for (bx = 0;bx < 4;bx++)
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
            for (by = 0;by < 4;by++)
              for (bx = 0;bx < 4;bx++)
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
*    Get cost for direct mode for an 8x8 block
************************************************************************
*/
int Get_Direct_Cost8x8 (int block, int *cost8x8)
{
  int block_y, block_x, pic_pix_y, pic_pix_x, i, j, k;
  int curr_diff[8][8], diff8x8[64];
  int diff[16];
  int cost  = 0;
  int mb_y  = (block/2)<<3;
  int mb_x  = (block%2)<<3;
  
  for (block_y=mb_y; block_y<mb_y+8; block_y+=4)
  {
    pic_pix_y = img->opix_y + block_y;
    
    for (block_x=mb_x; block_x<mb_x+8; block_x+=4)
    {
      pic_pix_x = img->opix_x + block_x;
      
      if (direct_pdir[pic_pix_y>>2][pic_pix_x>>2]<0)
      {
        *cost8x8=INT_MAX;
        return INT_MAX; //mode not allowed 
      }
      
      //===== prediction of 4x4 block =====
      
      LumaPrediction4x4 (block_x, block_y, direct_pdir[pic_pix_y>>2][pic_pix_x>>2], 0, 0, 
        direct_ref_idx[LIST_0][pic_pix_y>>2][pic_pix_x>>2], 
        direct_ref_idx[LIST_1][pic_pix_y>>2][pic_pix_x>>2]);
      
      //===== get displaced frame difference ======                
      for (k=j=0; j<4; j++)
        for (i=0; i<4; i++, k++)
        {
          diff[k] = curr_diff[block_y-mb_y+j][block_x-mb_x+i] =
            imgY_org[pic_pix_y+j][pic_pix_x+i] - img->mpr[j+block_y][i+block_x];
        }
        
        cost += SATD (diff, input->hadamard);
    }
  }
  
  if((input->rdopt==0)&&(input->Transform8x8Mode))
  {
    k=0;
    for(j=0; j<8; j++, k+=8)
      memcpy(&diff8x8[k], &(curr_diff[j]), 8 * sizeof(int));          
    
    *cost8x8 += SATD8X8 (diff8x8, input->hadamard);
  }
  
  return cost;
}



/*!
************************************************************************
* \brief
*    Get cost for direct mode for an macroblock
************************************************************************
*/
int Get_Direct_CostMB (int lambda_factor)
{
  int i;
  int cost = 0;
  int cost8x8 = 0;
  
  for (i=0; i<4; i++)
  {
    cost += Get_Direct_Cost8x8 (i, &cost8x8);
    if (cost8x8 == INT_MAX) return INT_MAX;
  }
  
  switch(input->Transform8x8Mode)
  {
  case 1: // Mixture of 8x8 & 4x4 transform
    if((cost8x8 < cost)||
      !(input->InterSearch8x4 &&
      input->InterSearch4x8 &&
      input->InterSearch4x4)
      )
    {
      cost = cost8x8; //return 8x8 cost
    }
    break;
  case 2: // 8x8 Transform only
    cost = cost8x8;
    break;
  default: // 4x4 Transform only
    break;
  }
  
  return cost;
  // T.Nishi(MEI ACC) 04-28-2004 end
}


/*!
************************************************************************
* \brief
*    Motion search for a partition
************************************************************************
*/
void
PartitionMotionSearch (int    blocktype,
                       int    block8x8,
                       int    lambda_factor
#ifdef RDO_Q
                       ,int   transform8x8
#endif
                       )
{
  static int  bx0[5][4] = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,2,0,0}, {0,2,0,2}};
  static int  by0[5][4] = {{0,0,0,0}, {0,0,0,0}, {0,2,0,0}, {0,0,0,0}, {0,0,2,2}};
  
  char  **ref_array;
  short ***mv_array;
  short *all_mv;
  short ref;
  int   v, h, mcost, search_range, i, j;
  int   pic_block_x, pic_block_y;
  int   bslice    = (img->type==B_SLICE);
  int   parttype  = (blocktype<4?blocktype:4);
  int   step_h0   = (input->part_size[ parttype][0]);
  int   step_v0   = (input->part_size[ parttype][1]);
  int   step_h    = (input->part_size[blocktype][0]);
  int   step_v    = (input->part_size[blocktype][1]);
  int   list;
  int   numlists  = bslice ? 2 : 1;
  int   list_offset = img->mb_data[img->current_mb_nr].list_offset; 
  int   block_y;
  int   *m_cost = NULL;
  int   by = by0[parttype][block8x8];
  int   bx = bx0[parttype][block8x8];
  
#ifdef SWITCHED_FILTERS
  int   mcost_temp;
#endif

  //===== LOOP OVER REFERENCE FRAMES =====
  for (list=0; list<numlists;list++)
  {
    for (ref=0; ref < listXsize[list+list_offset]; ref++)
    {
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
      if(input->UseRDO_Q || img->slice_fractional_quant_flag)
#else
      if(input->UseRDO_Q)
#endif
      {
        if (Motion_Selected == 0)  
          m_cost = &motion_cost[blocktype][list][ref][block8x8];
        else
        {
#ifdef MB32X32
          if(input->UseExtMB > 0)
          {
            if (transform8x8) motion_cost[blocktype][list][ref][block8x8] = motion16_cost8[curr_mbi][blocktype][list][ref][block8x8];
            else motion_cost[blocktype][list][ref][block8x8] = motion16_cost4[curr_mbi][blocktype][list][ref][block8x8];
          }
          else
          {
            if (transform8x8) motion_cost[blocktype][list][ref][block8x8] = motion_cost8[blocktype][list][ref][block8x8];
            else motion_cost[blocktype][list][ref][block8x8] = motion_cost4[blocktype][list][ref][block8x8];
          }
#else
          if (transform8x8) motion_cost[blocktype][list][ref][block8x8] = motion_cost8[blocktype][list][ref][block8x8];
          else motion_cost[blocktype][list][ref][block8x8] = motion_cost4[blocktype][list][ref][block8x8];
#endif
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
      if (!(input->UseRDO_Q  || img->slice_fractional_quant_flag) || Motion_Selected == 0)
#else
      if (!(input->UseRDO_Q) || Motion_Selected == 0)
#endif
#endif
        *m_cost = 0;
      
      //===== LOOP OVER SUB MACRO BLOCK partitions
      for (v=by; v<by + step_v0; v += step_v)
      {
        pic_block_y = img->block_y + v;
        
        for (h=bx; h<bx+step_h0; h+=step_h)
        {
          all_mv = img->all_mv[v][h][list][ref][blocktype];
          pic_block_x = img->block_x + h;
          
          //--- motion search for block ---
          
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
          if(input->UseRDO_Q || img->slice_fractional_quant_flag) 
#else
          if(input->UseRDO_Q)
#endif
          {
            if (Motion_Selected == 0)
            { 
              int       bsx       = input->blc_size[blocktype][0];
              int       bsy       = input->blc_size[blocktype][1];

#ifdef SWITCHED_FILTERS
              img->currMEOffset[list] = 0;
              mcost = BlockMotionSearch     (ref, list, h<<2, v<<2, blocktype, search_range, lambda_factor);
              if(img->filterParam == SIFO_FIRST_PASS_FPO)
              {
                img->min_mcost_offset[list] = mcost;
                if (img->noMEOffsets[list] > 0 && (ref == 0))
                {
                  for (i=0; i<img->noMEOffsets[list]; i++)
                  {
                    img->currMEOffset[list] = img->offsetMETab[list][i];
                    mcost_temp = BlockMotionSearch(ref, list, h<<2, v<<2, blocktype, search_range, lambda_factor); 
                    mcost = (mcost_temp < mcost)? mcost_temp: mcost;
                    img->min_mcost_offset[list] = mcost;
                  }
                }
              }
#else
              mcost = BlockMotionSearch     (ref, list, h<<2, v<<2, blocktype, search_range, lambda_factor);
#endif

              *m_cost += mcost;
              
              if (transform8x8)
              {
                for (j=v; j < v + (bsy>>2); j++)    
                {
                  for (i=h; i < h + (bsx>>2); i++)
                  {
#ifdef MB32X32
                    if(input->UseExtMB > 0)
                    {
                    img->all16_mv8[curr_mbi][j][i][list][ref][blocktype][0] = img->all_mv[j][i][list][ref][blocktype][0];
                    img->all16_mv8[curr_mbi][j][i][list][ref][blocktype][1] = img->all_mv[j][i][list][ref][blocktype][1];

                      img->bipred_mv1_8[curr_mbi][j][i][list][ref][blocktype][0] = img->bipred_mv1[j][i][list][ref][blocktype][0];
                      img->bipred_mv1_8[curr_mbi][j][i][list][ref][blocktype][1] = img->bipred_mv1[j][i][list][ref][blocktype][1];

                      img->bipred_mv2_8[curr_mbi][j][i][list][ref][blocktype][0] = img->bipred_mv2[j][i][list][ref][blocktype][0];
                      img->bipred_mv2_8[curr_mbi][j][i][list][ref][blocktype][1] = img->bipred_mv2[j][i][list][ref][blocktype][1];
                    }
                    else
                    {
                      img->all_mv8[j][i][list][ref][blocktype][0] = img->all_mv[j][i][list][ref][blocktype][0];
                      img->all_mv8[j][i][list][ref][blocktype][1] = img->all_mv[j][i][list][ref][blocktype][1];
                    }
#else
                    img->all_mv8[j][i][list][ref][blocktype][0] = img->all_mv[j][i][list][ref][blocktype][0];
                    img->all_mv8[j][i][list][ref][blocktype][1] = img->all_mv[j][i][list][ref][blocktype][1];
#endif
                  }
                }
              }
              else
              {
                for (j=v; j < v + (bsy>>2); j++)    
                {
                  for (i=h; i < h + (bsx>>2); i++)
                  {
#ifdef MB32X32
                    if(input->UseExtMB>0)
                    {
                    img->all16_mv4[curr_mbi][j][i][list][ref][blocktype][0] = img->all_mv[j][i][list][ref][blocktype][0];
                    img->all16_mv4[curr_mbi][j][i][list][ref][blocktype][1] = img->all_mv[j][i][list][ref][blocktype][1];

                      img->bipred_mv1_4[curr_mbi][j][i][list][ref][blocktype][0] = img->bipred_mv1[j][i][list][ref][blocktype][0];
                      img->bipred_mv1_4[curr_mbi][j][i][list][ref][blocktype][1] = img->bipred_mv1[j][i][list][ref][blocktype][1];

                      img->bipred_mv2_4[curr_mbi][j][i][list][ref][blocktype][0] = img->bipred_mv2[j][i][list][ref][blocktype][0];
                      img->bipred_mv2_4[curr_mbi][j][i][list][ref][blocktype][1] = img->bipred_mv2[j][i][list][ref][blocktype][1];
                    }
                    else
                    {
                      img->all_mv4[j][i][list][ref][blocktype][0] = img->all_mv[j][i][list][ref][blocktype][0];
                      img->all_mv4[j][i][list][ref][blocktype][1] = img->all_mv[j][i][list][ref][blocktype][1];
                    }
#else
                    img->all_mv4[j][i][list][ref][blocktype][0] = img->all_mv[j][i][list][ref][blocktype][0];
                    img->all_mv4[j][i][list][ref][blocktype][1] = img->all_mv[j][i][list][ref][blocktype][1];
#endif
                  }
                }
              }
            }
            else
            {
              int       bsx       = input->blc_size[blocktype][0];
              int       bsy       = input->blc_size[blocktype][1];
              if (transform8x8)
              {
                for (j=v; j < v + (bsy>>2); j++)    
                {
                  for (i=h; i < h + (bsx>>2); i++)
                  {
#ifdef MB32X32
                    if(input->UseExtMB>0)
                    {
                    img->all_mv[j][i][list][ref][blocktype][0] = img->all16_mv8[curr_mbi][j][i][list][ref][blocktype][0];
                    img->all_mv[j][i][list][ref][blocktype][1] = img->all16_mv8[curr_mbi][j][i][list][ref][blocktype][1];

                      img->bipred_mv1[j][i][list][ref][blocktype][0] = img->bipred_mv1_8[curr_mbi][j][i][list][ref][blocktype][0];
                      img->bipred_mv1[j][i][list][ref][blocktype][1] = img->bipred_mv1_8[curr_mbi][j][i][list][ref][blocktype][1];

                      img->bipred_mv2[j][i][list][ref][blocktype][0] = img->bipred_mv2_8[curr_mbi][j][i][list][ref][blocktype][0];
                      img->bipred_mv2[j][i][list][ref][blocktype][1] = img->bipred_mv2_8[curr_mbi][j][i][list][ref][blocktype][1];
                    }
                    else
                    {
                      img->all_mv[j][i][list][ref][blocktype][0] = img->all_mv8[j][i][list][ref][blocktype][0];
                      img->all_mv[j][i][list][ref][blocktype][1] = img->all_mv8[j][i][list][ref][blocktype][1];
                    }
#else
                    img->all_mv[j][i][list][ref][blocktype][0] = img->all_mv8[j][i][list][ref][blocktype][0];
                    img->all_mv[j][i][list][ref][blocktype][1] = img->all_mv8[j][i][list][ref][blocktype][1];
#endif
                  }
                }
              }
              else
              {
                for (j=v; j < v + (bsy>>2); j++)    
                {
                  for (i=h; i < h + (bsx>>2); i++)
                  {
#ifdef MB32X32
                    if(input->UseExtMB>0)
                    {
                    img->all_mv[j][i][list][ref][blocktype][0] = img->all16_mv4[curr_mbi][j][i][list][ref][blocktype][0];
                    img->all_mv[j][i][list][ref][blocktype][1] = img->all16_mv4[curr_mbi][j][i][list][ref][blocktype][1];

                      img->bipred_mv1[j][i][list][ref][blocktype][0] = img->bipred_mv1_4[curr_mbi][j][i][list][ref][blocktype][0];
                      img->bipred_mv1[j][i][list][ref][blocktype][1] = img->bipred_mv1_4[curr_mbi][j][i][list][ref][blocktype][1];

                      img->bipred_mv2[j][i][list][ref][blocktype][0] = img->bipred_mv2_4[curr_mbi][j][i][list][ref][blocktype][0];
                      img->bipred_mv2[j][i][list][ref][blocktype][1] = img->bipred_mv2_4[curr_mbi][j][i][list][ref][blocktype][1];
                    }
                    else
                    {
                      img->all_mv[j][i][list][ref][blocktype][0] = img->all_mv4[j][i][list][ref][blocktype][0];
                      img->all_mv[j][i][list][ref][blocktype][1] = img->all_mv4[j][i][list][ref][blocktype][1];
                    }
#else
                    img->all_mv[j][i][list][ref][blocktype][0] = img->all_mv4[j][i][list][ref][blocktype][0];
                    img->all_mv[j][i][list][ref][blocktype][1] = img->all_mv4[j][i][list][ref][blocktype][1];
#endif
                  }
                }
              }
            }
          }
          else
          {
#ifdef SWITCHED_FILTERS
            img->currMEOffset[list] = 0;
            mcost = BlockMotionSearch(ref, list, h<<2, v<<2, blocktype, search_range, lambda_factor);
            if(img->filterParam == SIFO_FIRST_PASS_FPO)
            {
              img->min_mcost_offset[list] = mcost;
              if(img->noMEOffsets[list] > 0 && ref == 0)
              {
                for(i=0; i<img->noMEOffsets[list]; i++)
                {
                  img->currMEOffset[list] = img->offsetMETab[list][i];
                  mcost_temp = BlockMotionSearch(ref, list, h<<2, v<<2, blocktype, search_range, lambda_factor); 
                  mcost = (mcost_temp < mcost)? mcost_temp: mcost;
                  img->min_mcost_offset[list] = mcost;
                }
              }    
            }
#else
            mcost = BlockMotionSearch     (ref, list, h<<2, v<<2, blocktype, search_range, lambda_factor);
#endif
            *m_cost += mcost;
          }
#else
          mcost = BlockMotionSearch     (ref, list, h<<2, v<<2, blocktype, search_range, lambda_factor);
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
      if(input->UseRDO_Q || img->slice_fractional_quant_flag)
#else
      if(input->UseRDO_Q)
#endif
      {
        if (Motion_Selected == 0)
        {
#ifdef MB32X32

          if (transform8x8)
          {
            if(input->UseExtMB>0)
            {
              motion16_cost8[curr_mbi][blocktype][list][ref][block8x8]= motion_cost[blocktype][list][ref][block8x8];
            }
            else
            {
              motion_cost8[blocktype][list][ref][block8x8]= motion_cost[blocktype][list][ref][block8x8];
            }
          }
          else
          {
            if(input->UseExtMB>0)
            {
              motion16_cost4[curr_mbi][blocktype][list][ref][block8x8]= motion_cost[blocktype][list][ref][block8x8];
            }
            else
            {
              motion_cost4[blocktype][list][ref][block8x8]= motion_cost[blocktype][list][ref][block8x8];
            }
          }

#else
          if (transform8x8)
            motion_cost8[blocktype][list][ref][block8x8]= motion_cost[blocktype][list][ref][block8x8];
          else
            motion_cost4[blocktype][list][ref][block8x8]= motion_cost[blocktype][list][ref][block8x8];
#endif
        }
        else
        {
#ifdef MB32X32
          if(input->UseExtMB>0)
          {
            if (transform8x8)
              motion_cost[blocktype][list][ref][block8x8]= motion16_cost8[curr_mbi][blocktype][list][ref][block8x8];
            else
              motion_cost[blocktype][list][ref][block8x8]= motion16_cost4[curr_mbi][blocktype][list][ref][block8x8];
          }
          else
          {
            if (transform8x8)
              motion_cost[blocktype][list][ref][block8x8]= motion_cost8[blocktype][list][ref][block8x8];
            else
              motion_cost[blocktype][list][ref][block8x8]= motion_cost4[blocktype][list][ref][block8x8];
          }
#else
          if (transform8x8)
            motion_cost[blocktype][list][ref][block8x8]= motion_cost8[blocktype][list][ref][block8x8];
          else
            motion_cost[blocktype][list][ref][block8x8]= motion_cost4[blocktype][list][ref][block8x8];
#endif
        }
      }
#endif
    }
  }
}

/*!
************************************************************************
* \brief
*    Calculate Direct Motion Vectors  *****
************************************************************************
*/
void Get_Direct_Motion_Vectors ()
{
  
  int   block_x, block_y, pic_block_x, pic_block_y, opic_block_x, opic_block_y;
  short ****all_mvs;
  int   mv_scale;
  int refList; 
  int ref_idx;   
  
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
    
    short l0_refA, l0_refB, l0_refD, l0_refC;
    short l1_refA, l1_refB, l1_refD, l1_refC; 
    short l0_refX,l1_refX;
    short pmvfw[2]={0,0},pmvbw[2]={0,0};
    
    PixelPos mb_a, mb_b, mb_d, mb_c;              
    
    getLuma4x4Neighbour(img->current_mb_nr, 0, 0, -1,  0,&mb_a);
    getLuma4x4Neighbour(img->current_mb_nr, 0, 0,  0, -1,&mb_b);
    getLuma4x4Neighbour(img->current_mb_nr, 0, 0, 16, -1,&mb_c);
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
    
    l0_refX = (l0_refA >= 0 && l0_refB >= 0) ? min(l0_refA,l0_refB): max(l0_refA,l0_refB);
    l0_refX = (l0_refX >= 0 && l0_refC >= 0) ? min(l0_refX,l0_refC): max(l0_refX,l0_refC);
    
    l1_refX = (l1_refA >= 0 && l1_refB >= 0) ? min(l1_refA,l1_refB): max(l1_refA,l1_refB);
    l1_refX = (l1_refX >= 0 && l1_refC >= 0) ? min(l1_refX,l1_refC): max(l1_refX,l1_refC);        
    
    if (l0_refX >=0)
      SetMotionVectorPredictor (pmvfw, enc_picture->ref_idx[LIST_0], enc_picture->mv[LIST_0], l0_refX, LIST_0, 0, 0, 16, 16);
    
    if (l1_refX >=0)
      SetMotionVectorPredictor (pmvbw, enc_picture->ref_idx[LIST_1], enc_picture->mv[LIST_1], l1_refX, LIST_1, 0, 0, 16, 16);
    
    for (block_y=0; block_y<4; block_y++)
    {
      pic_block_y  = (img->pix_y  >> 2) + block_y;
      opic_block_y = (img->opix_y >> 2) + block_y;
      
      for (block_x=0; block_x<4; block_x++)
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
}

/*!
************************************************************************
* \brief
*    control the sign of a with b
************************************************************************
*/
int sign(int a,int b)
{
  return (b >= 0) ? absm(a) : -absm(a);
}



