
/*!
 *************************************************************************************
 * \file adaptive_loop_filter.c
 *
 * \brief
 *    adaptive loop filter routines
 *
 * \author
 *    Takeshi Chujoh         <takeshi.chujoh@toshiba.co.jp>             \n
 *    Akiyuki Tanizawa       <akiyuki.tanizawa@toshiba.co.jp>           \n
 *    Goki Yasuda            <goki1.yasuda@toshiba.co.jp>               \n
 *    Naofumi Wada           <naofumi.wada@toshiba.co.jp>               \n
 *    Takashi Watanabe       <takashi39.watanabe@toshiba.co.jp>         \n
 *************************************************************************************
 */

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "adaptive_loop_filter.h"
#include "image.h"
#include "memalloc.h"
#include "vlc.h"
#include "cabac.h"
#include "header.h"

#ifdef ADAPTIVE_LOOP_FILTER

//------- global parameters for adaptive_loop_filter.c
int    num_bit;
int    r_offset;
int    ext_offset;
int    base_type;
double lambda_alf_luma;
double lambda_alf_chroma;
int64  rate, min_rate, orig_rate;
int64  dist, min_dist, orig_dist;
double cost, min_cost, orig_cost;
int64  tmp_luma_rate;
int    temp_pic_bin_count;

alf_parameter_set_rbsp_t *temp_alfps;
double **alf_corr;
double *dbl_alf_coef;
imgpel **imgY_rest;
imgpel ***imgUV_rest;
imgpel **imgY_ext;
imgpel ***imgUV_ext;
imgpel **imgY_temp;
imgpel **imgY_best;
int    **alf_coef_base;
int    ***alf_coef_tap_base;
int    **alf_coef_base_chroma;
Bitstream *tempStream;
Bitstream *alfstream;
unsigned int ****blk_corr;
AlfQtBlock   *qt_blk_data;

int filter_base_luma[42] = {
  0,  0,  0,  0,   0,  0,  0,  0, 0,
  0,  0,  0, -1,   5, -1,  0,  0, 0,
  0,  0,  0,  0,  -7,  0,  0,  0, 0,
  0, -1,  0, -5,  18, -5,  0, -1, 0,
  0,  5, -7, 18, 219,  0
};
int filter_base_chroma[14] = {
   5, -5,  -5, -5,  5,
  -5,  0,  40,  0, -5,
  -5, 40, 160,  0
};

int filter_pos_5x5_in_9x9[14] = {
  20, 21, 22, 23, 24,
  29, 30, 31, 32, 33,
  38, 39, 40, 41
};
int filter_pos_7x7_in_9x9[26] = {
  10, 11, 12, 13, 14, 15, 16,
  19, 20, 21, 22, 23, 24, 25,
  28, 29, 30, 31, 32, 33, 34,
  37, 38, 39, 40, 41
};
int filter_pos_9x9_in_9x9[42] = {
   0,  1,  2,  3,  4,  5,  6,  7,  8, 
   9, 10, 11, 12, 13, 14, 15, 16, 17,
  18, 19, 20, 21, 22, 23, 24, 25, 26,
  27, 28, 29, 30, 31, 32, 33, 34, 35,
  36, 37, 38, 39, 40, 41
};

int symmetric_array9x9[81] = {
   0,  1,  2,  3,  4,  5,  6,  7,  8,
   9, 10, 11, 12, 13, 14, 15, 16, 17,    
  18, 19, 20, 21, 22, 23, 24, 25, 26,
  27, 28, 29, 30, 31, 32, 33, 34, 35,       
  36, 37, 38, 39, 40, 39, 38, 37, 36, 
  35, 34, 33, 32, 31, 30, 29, 28, 27, 
  26, 25, 24, 23, 22, 21, 20, 19, 18,
  17, 16, 15, 14, 13, 12, 11, 10,  9, 
   8,  7,  6,  5,  4,  3,  2,  1,  0
};
int symmetric_array7x7[49] = {
   0,  1,  2,  3,  4,  5,  6,    
   7,  8,  9, 10, 11, 12, 13,
  14, 15, 16, 17, 18, 19, 20,     
  21, 22, 23, 24, 23, 22, 21, 
  20, 19, 18, 17, 16, 15, 14, 
  13, 12, 11, 10,  9,  8,  7,
   6,  5,  4,  3,  2,  1,  0, 
};
int symmetric_array5x5[25] = {    
   0,  1,  2,  3,  4, 
   5,  6,  7,  8,  9,        
  10, 11, 12, 11, 10,  
   9,  8,  7,  6,  5,  
   4,  3,  2,  1,  0,  
};

// for quantization of filter coef.
int symmetric_mag9x9[41] = {
  2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 1
};
int symmetric_mag7x7[25] = {
  2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 1
};
int symmetric_mag5x5[13] = {
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 1
};

int filter_base_flag[3][4];

const int block_scan_order[5][256] = {
  { 0 },  // layer 0
  { 0, 1,
    2, 3 },  // layer 1
  { 0, 1, 4, 5,
    2, 3, 6, 7,
    8, 9,12,13,
   10,11,14,15 },  // layer 2
  { 0, 1, 4, 5,16,17,20,21,
    2, 3, 6, 7,18,19,22,23,
    8, 9,12,13,24,25,28,29,
   10,11,14,15,26,27,30,31,
   32,33,36,37,48,49,52,53,
   34,35,38,39,50,51,54,55,
   40,41,44,45,56,57,60,61,
   42,43,46,47,58,59,62,63 }, // layer 3
  {  0,  1,  4,  5, 16, 17, 20, 21, 64, 65, 68, 69, 80, 81, 84, 85,
     2,  3,  6,  7, 18, 19, 22, 23, 66, 67, 70, 71, 82, 83, 86, 87,
     8,  9, 12, 13, 24, 25, 28, 29, 72, 73, 76, 77, 88, 89, 92, 93,
    10, 11, 14, 15, 26, 27, 30, 31, 74, 75, 78, 79, 90, 91, 94, 95,
    32, 33, 36, 38, 48, 49, 52, 53, 96, 97,100,102,112,113,116,117,
    34, 35, 37, 39, 50, 51, 54, 55, 98, 99,101,103,114,115,118,119,
    40, 41, 44, 45, 56, 57, 60, 61,104,105,108,109,120,121,124,125,
    42, 43, 46, 47, 58, 59, 62, 63,106,107,110,111,122,123,126,127,
   128,129,132,133,144,145,148,149,192,193,196,197,208,209,212,213,
   130,131,134,135,146,147,150,151,194,195,198,199,210,211,214,215,
   136,137,140,141,152,153,156,157,200,201,204,205,216,217,220,221,
   138,139,142,143,154,155,158,159,202,203,206,207,218,219,222,223,
   160,161,164,166,176,177,180,181,224,225,228,230,240,241,244,245,
   162,163,165,167,178,179,182,183,226,227,229,231,242,243,246,247,
   168,169,172,173,184,185,188,189,232,233,236,237,248,249,252,253,
   170,171,174,175,186,187,190,191,234,235,238,239,250,251,254,255 } // layer 4
};

const int pyramid_mem_idx[5][256] = {
  { 0 },  // layer 0
  { 1, 2,
    3, 4 },  // layer 1
  { 5, 6, 7, 8,
    9,10,11,12,
   13,14,15,16,
   17,18,19,20 },  // layer 2
  {21,22,23,24,25,26,27,28,
   29,30,31,32,33,34,35,36,
   37,38,39,40,41,42,43,44,
   45,46,47,48,49,50,51,52,
   53,54,55,56,57,58,59,60,
   61,62,63,64,65,66,67,68,
   69,70,71,72,73,74,75,76,
   77,78,79,80,81,82,83,84 }, // layer 3
  {85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,
  101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,
  117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,
  133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,
  149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,
  165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,
  181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,
  197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,
  213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,
  229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,
  245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,
  261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,
  277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,
  293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,
  309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,
  325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340 } // layer 4
};
//------- end


/*! 
*************************************************************************************
* \brief
*    Allocates memory for an alfps
*
* \return
*    pointer to a alfps
*************************************************************************************
*/
alf_parameter_set_rbsp_t *AllocALFPS (int *memory_size)
{
  alf_parameter_set_rbsp_t *p;

  if ((p=calloc (sizeof (alf_parameter_set_rbsp_t), 1)) == NULL)
    no_mem_exit ("AllocALFPS: ALFPS");
  *memory_size += sizeof(alf_parameter_set_rbsp_t);

  *memory_size += get_mem1Dint(&(p->alf.coeff), MAX_NUM_COEF);
  *memory_size += get_mem1Dint(&(p->alf.coeff_chroma), MAX_NUM_COEF_C);

  if(((p->alf_blk_data) = (AlfBlock *) calloc(img->max_num_alf_block, sizeof(AlfBlock))) == NULL)
    no_mem_exit("AllocALFPS: p->alf_blk_data");
  *memory_size += sizeof(AlfBlock)*img->max_num_alf_block;

  return p;
}

/*! 
*************************************************************************************
* \brief
*    Frees an alfps
*
* \param alfps
*     alfps to be freed
*
* \return
*    none
*************************************************************************************
*/
void FreeALFPS (alf_parameter_set_rbsp_t *p)
{
  assert(p != NULL);

  if(p->alf.coeff != NULL)
  {
    free_mem1Dint (p->alf.coeff);
    p->alf.coeff=NULL;
  }

  if(p->alf.coeff_chroma != NULL)
  {
    free_mem1Dint (p->alf.coeff_chroma);
    p->alf.coeff_chroma=NULL;
  }

  if(p->alf_blk_data != NULL)
  {
    free(p->alf_blk_data);
    p->alf_blk_data=NULL;
  }

  free (p);
  p=NULL;
}

/*!
************************************************************************
* \brief
*    Dynamic memory allocation of frame size related global buffers
*    buffers are defined in global.h, allocated memory must be freed in
*    void FreeALFGlobalBurrers()
*
*  \par Output:
*     Number of allocated bytes
***********************************************************************
*/
int  InitALFGlobalBuffers(void)
{
  int memory_size = 0;
  int idx, num_of_block;
  
  memory_size += get_mem2Ddouble (&(alf_corr),  MAX_NUM_COEF, MAX_NUM_COEF+1);
  memory_size += get_mem1Ddouble (&(dbl_alf_coef), MAX_NUM_COEF);
  memory_size += get_mem2Dpel (&(imgY_rest), input->img_height, input->img_width);
  memory_size += get_mem3Dpel (&(imgUV_rest), 2, input->img_height_cr, input->img_width_cr);
  memory_size += get_mem2Dpel (&(imgY_ext), input->img_height+MAX_NUM_TAP, input->img_width+MAX_NUM_TAP);
  memory_size += get_mem3Dpel (&(imgUV_ext), 2, input->img_height_cr+MAX_NUM_TAP_C, input->img_width_cr+MAX_NUM_TAP_C);
  memory_size += get_mem2Dint (&(alf_coef_base), 4, MAX_NUM_COEF);
  memory_size += get_mem3Dint (&(alf_coef_tap_base), 3, 4, MAX_NUM_COEF);
  memory_size += get_mem2Dint (&(alf_coef_base_chroma), 4, MAX_NUM_COEF_C);
  memory_size += get_mem2Dpel (&(imgY_temp), input->img_height, input->img_width);
  memory_size += get_mem2Dpel (&(imgY_best), input->img_height, input->img_width);

  // for block-based filtering
  for(idx=0; idx<8; idx++)
  {
    num_of_block = CalcNumOfAlfBlocks((1<<idx), input->img_width, input->img_height);
    if( num_of_block<=MAX_NUM_BLOCK ) break;
  }
  img->min_alf_block_size_idx = idx;
  img->min_alf_block_size     = 1<<idx;
  img->max_WidthInAlfBlks     = input->img_width  / img->min_alf_block_size;
  img->max_HeightInAlfBlks    = input->img_height / img->min_alf_block_size;
  if( (input->img_width  % img->min_alf_block_size)!=0 ) img->max_WidthInAlfBlks++;
  if( (input->img_height % img->min_alf_block_size)!=0 ) img->max_HeightInAlfBlks++;
  img->max_num_alf_block      = img->max_WidthInAlfBlks * img->max_HeightInAlfBlks;
  memory_size += get_mem_blk_corr (&(blk_corr), img->max_HeightInAlfBlks, img->max_WidthInAlfBlks, MIN_NUM_COEF, MIN_NUM_COEF+1);

  // allocate bitstream
  if ( (alfstream =  (Bitstream *) calloc(1, sizeof(Bitstream))) == NULL )
    no_mem_exit ("InitALFGlobalBuffers: alfstream");
  memory_size += sizeof(Bitstream);
  if ( (alfstream->streamBuffer = (byte *) calloc(MAX_ALF_PARAM_SIZE, sizeof(byte))) == NULL )
    no_mem_exit ("InitALFGlobalBuffers: alfstream->streamBuffer");
  memory_size += MAX_ALF_PARAM_SIZE*sizeof(byte);
  alfstream->bits_to_go=8;
  alfstream->byte_buf=0;
  alfstream->byte_pos=0;

  // allocate ALF parameter set
  temp_alfps   = AllocALFPS(&memory_size);
  active_alfps = AllocALFPS(&memory_size);

  // for quadtree-based filtering
  if((qt_blk_data = (AlfQtBlock *) calloc(1, sizeof(AlfQtBlock))) == NULL)
    no_mem_exit("InitALFGlobalBuffers: qt_blk_data");
  memory_size += sizeof(AlfQtBlock);

  return memory_size;
}

/*!
************************************************************************
* \brief
*    Free allocated memory of frame size related global buffers
*    buffers are defined in global.h, allocated memory is allocated in
*    int InitALFGlobalBuffers()
************************************************************************
*/
void FreeALFGlobalBurrers(void)
{
  free_mem2Ddouble (alf_corr);
  free_mem1Ddouble (dbl_alf_coef);
  free_mem2Dpel (imgY_rest);
  free_mem3Dpel (imgUV_rest, 2);
  free_mem2Dpel (imgY_ext);
  free_mem3Dpel (imgUV_ext, 2);
  free_mem2Dint (alf_coef_base);
  free_mem3Dint (alf_coef_tap_base, 3);
  free_mem2Dint (alf_coef_base_chroma);
  free_mem2Dpel (imgY_temp);
  free_mem2Dpel (imgY_best);

  free_mem_blk_corr (blk_corr, img->max_HeightInAlfBlks, img->max_WidthInAlfBlks);

  if(alfstream->streamBuffer != NULL)
  {
    free(alfstream->streamBuffer);
    alfstream->streamBuffer=NULL;
  }
  if(alfstream != NULL)
  {
    free(alfstream);
    alfstream=NULL;
  }

  FreeALFPS (temp_alfps);
  FreeALFPS (active_alfps);

  if(qt_blk_data != NULL)
  {
    free(qt_blk_data);
    qt_blk_data=NULL;
  }
}

/*!
 ************************************************************************
 * \brief
 *    set ALF parameters (p -> alfps)
 ************************************************************************
 */
void SetALFParameters( alf_parameter_set_rbsp_t *alfps, alf_parameter_set_rbsp_t *p )
{
  alfps->alf_flag              = p->alf_flag;

  if(alfps->alf_flag)
  {
#ifdef ALF_SPATIAL_PREDICT_COEF
    alfps->pred_coef_mode        = p->pred_coef_mode;
#endif
    alfps->alf.tap               = p->alf.tap;
    alfps->alf.num_coeff         = p->alf.num_coeff;
    memcpy(alfps->alf.coeff, p->alf.coeff, sizeof(int)*(alfps->alf.num_coeff));

    alfps->chroma_idc            = p->chroma_idc;
    if(alfps->chroma_idc)
    {
      alfps->alf.tap_chroma        = p->alf.tap_chroma;
      alfps->alf.num_coeff_chroma  = p->alf.num_coeff_chroma;
      memcpy(alfps->alf.coeff_chroma, p->alf.coeff_chroma, sizeof(int)*(alfps->alf.num_coeff_chroma));
    }

    alfps->block_control_flag    = p->block_control_flag;
    if(alfps->block_control_flag)
    {
      alfps->alf_block_size_idx  = p->alf_block_size_idx;
      alfps->alf_block_size      = p->alf_block_size;
      alfps->WidthInAlfBlks      = p->WidthInAlfBlks;
      alfps->HeightInAlfBlks     = p->HeightInAlfBlks;
      alfps->num_alf_block       = p->num_alf_block;
      memcpy( alfps->alf_blk_data, p->alf_blk_data, sizeof(AlfBlock)*(alfps->num_alf_block));
      alfps->qt_partition_flag   = p->qt_partition_flag;
      if(alfps->qt_partition_flag)
      {
        alfps->max_layer_level     = p->max_layer_level;
        alfps->parent_block_size   = p->parent_block_size;
        alfps->WidthInParentBlks   = p->WidthInParentBlks;
        alfps->HeightInParentBlks  = p->HeightInParentBlks;
        alfps->num_parent_block    = p->num_parent_block;
      }
    }
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    make temporally frame for filtering
 *****************************************************************************************
 */
void MakeExtendFrame(imgpel** ImgDec, imgpel** ImgExt, int width, int height, int offset)
{
  int x, y;

  for(y=0; y<offset; y++)
  {
    for(x=0; x<offset; x++)
    {
      ImgExt[y][x] = ImgDec[0][0];
    }
    for(x=offset; x<width+offset; x++)
    {
      ImgExt[y][x] = ImgDec[0][x-offset];
    }
    for(x=width+offset; x<width+(offset<<1); x++)
    {
      ImgExt[y][x] = ImgDec[0][width-1];
    }
  }
  for(y=offset; y<height+offset; y++)
  {
    for(x=0; x<offset; x++)
    {
      ImgExt[y][x] = ImgDec[y-offset][0];
    }
    for(x=offset; x<width+offset; x++)
    {
      ImgExt[y][x] = ImgDec[y-offset][x-offset];
    }
    for(x=width+offset; x<width+(offset<<1); x++)
    {
      ImgExt[y][x] = ImgDec[y-offset][width-1];
    }
  }
  for(y=height+offset; y<height+(offset<<1); y++)
  {
    for(x=0; x<offset; x++)
    {
      ImgExt[y][x] = ImgDec[height-1][0];
    }
    for(x=offset; x<width+offset; x++)
    {
      ImgExt[y][x] = ImgDec[height-1][x-offset];
    }
    for(x=width+offset; x<width+(offset<<1); x++)
    {
      ImgExt[y][x] = ImgDec[height-1][width-1];
    }
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    set restored frame
 *****************************************************************************************
 */
void SetRestFrame(imgpel** ImgRest, imgpel** ImgDec, int width, int height)
{
  int i;

  for(i=0; i<height; i++)
    memcpy(ImgDec[i], ImgRest[i], sizeof(imgpel)*width);
}

#ifdef ALF_SPATIAL_PREDICT_COEF
 /*! 
 ************************************************************************
 * \brief
 *    PredictALFCoeff: 
 *    Predicts one luma filter coeff based on sum of other luma coeff
 ************************************************************************
 */
void PredictALFCoeff( alf_parameter_set_rbsp_t *alfps )
{
  int i, sum, pred, tap, N;
  int    *pFiltMag;

  tap = alfps->alf.tap;
  switch(tap)
  {
  case 5:
    pFiltMag = symmetric_mag5x5;
    break;
  case 7:
    pFiltMag = symmetric_mag7x7;
    break;
  case 9:
    pFiltMag = symmetric_mag9x9;
    break;
  default:
    printf("not support filter size : %d", tap);
    exit(1);
    break;
  }
  N = (tap*tap+1)>>1;
  sum=0;	
  for(i=0; i<N;i++) 
  {
    sum+=pFiltMag[i]*alfps->alf.coeff[i];
  }
  pred=(1<<num_bit)-(sum-alfps->alf.coeff[N-1]);
  alfps->alf.coeff[N-1]=pred-alfps->alf.coeff[N-1];
}

 /*! 
 ************************************************************************
 * \brief
 *    PredictALFCoeffChroma: 
 *    Predicts one chroma filter coeff based on sum of other chroma coeff
 ************************************************************************
 */
void PredictALFCoeffChroma( alf_parameter_set_rbsp_t *alfps )
{
  int i, sum, pred, tap, N;
  int    *pFiltMag;
  
  tap  = alfps->alf.tap_chroma;
  switch(tap)
  {
  case 5:
    pFiltMag = symmetric_mag5x5;
    break;
  case 7:
    pFiltMag = symmetric_mag7x7;
    break;
  case 9:
    pFiltMag = symmetric_mag9x9;
    break;
  default:
    printf("not support filter size : %d", tap);
    exit(1);
    break;
  }
  N = (tap*tap+1)>>1;
  sum=0;	
  for(i=0; i<N;i++)	
  {	
    sum+=pFiltMag[i]*alfps->alf.coeff_chroma[i];
  }		
  pred=(1<<num_bit)-(sum-alfps->alf.coeff_chroma[N-1]);    
  alfps->alf.coeff_chroma[N-1]=pred-alfps->alf.coeff_chroma[N-1];
}
#endif

/*!
 *****************************************************************************************
 * \brief
 *    get base type for filter subtraction
 *****************************************************************************************
 */
int GetALFBaseType(int slice_type, int ref)
{
  int type=0;

  if ( (slice_type==I_SLICE) || (slice_type==SI_SLICE) )
  {
    type = 0;
  }
  else if ( (slice_type==P_SLICE) || (slice_type==SP_SLICE) )
  {
    type = 1;
  }
  else
  {
    type = (ref>0) ? 2 : 3;
  }
  
  return type;
}

/*!
 *****************************************************************************************
 * \brief
 *    initialize base of filter coefficients
 *****************************************************************************************
 */
void InitFilterBase()
{
  int i, tap_idx, type;

  // luma filter base
  for(tap_idx=0; tap_idx<3; tap_idx++)
  {
    for(type=0; type<4; type++)
    {
      for(i=0; i<MAX_NUM_COEF; i++)
      {
        alf_coef_base[type][i] = filter_base_luma[i];
      }
      filter_base_flag[tap_idx][type] = 0;
    }
  }

  // chroma filter base
  for(type=0; type<4; type++)
  {
    for(i=0; i<MAX_NUM_COEF_C; i++)
    {
      alf_coef_base_chroma[type][i] = filter_base_chroma[i];
    }
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    set difference of filter coefficients 
 *****************************************************************************************
 */
void SetALFDiffCoef(alf_parameter_set_rbsp_t *alfps)
{
  int i, N, tap;
  int *coeff;
  int *filter_pos;
  int *base_flag;
  int **base_coef;

  tap   = alfps->alf.tap;
  N     = (tap*tap+1)>>1;
  coeff = alfps->alf.coeff;

  switch(tap)
  {
  case 5:
    filter_pos = filter_pos_5x5_in_9x9;
    base_flag  = filter_base_flag[0];
    base_coef  = alf_coef_tap_base[0];
    break;
  case 7:
    filter_pos = filter_pos_7x7_in_9x9;
    base_flag  = filter_base_flag[1];
    base_coef  = alf_coef_tap_base[1];
    break;
  case 9:
    filter_pos = filter_pos_9x9_in_9x9;
    base_flag  = filter_base_flag[2];
    base_coef  = alf_coef_tap_base[2];
    break;
  default:
    printf("not support filter size: %d\n", tap);
    exit(1);
    break;
  }

  // set diff coeff
  if(base_flag[base_type] == 0)
  {
    for(i=0; i<N; i++)
    {
      coeff[i] -= alf_coef_base[base_type][filter_pos[i]];
    }
  }
  else
  {
    for(i=0; i<N; i++)
    {
      coeff[i] -= base_coef[base_type][i];
    }
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    set base of filter coefficients 
 *****************************************************************************************
 */
void SetALFBaseCoef(alf_parameter_set_rbsp_t *alfps)
{
  int i, N, tap;
  int *coeff;
  int *filter_pos;
  int *base_flag;
  int **base_coef;

  tap   = alfps->alf.tap;
  N     = (tap*tap+1)>>1;
  coeff = alfps->alf.coeff;

  switch(tap)
  {
  case 5:
    filter_pos = filter_pos_5x5_in_9x9;
    base_flag  = filter_base_flag[0];
    base_coef  = alf_coef_tap_base[0];
    break;
  case 7:
    filter_pos = filter_pos_7x7_in_9x9;
    base_flag  = filter_base_flag[1];
    base_coef  = alf_coef_tap_base[1];
    break;
  case 9:
    filter_pos = filter_pos_9x9_in_9x9;
    base_flag  = filter_base_flag[2];
    base_coef  = alf_coef_tap_base[2];
    break;
  default:
    printf("not support filter size: %d\n", tap);
    exit(1);
    break;
  }

  // set base coeff
  if(base_flag[base_type] == 0)
  {
    for(i=0; i<N; i++)
    {
      alf_coef_base[base_type][filter_pos[i]] += coeff[i];
      base_coef[base_type][i] = alf_coef_base[base_type][filter_pos[i]];
    }
    base_flag[base_type] = 1;
  }
  else
  {
    for(i=0; i<N; i++)
    {
      base_coef[base_type][i] += coeff[i];
    }
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    set difference of filter coefficients for chroma
 *****************************************************************************************
 */
void SetALFDiffCoefChroma(alf_parameter_set_rbsp_t *alfps)
{
  int i, N, tap;

  tap = alfps->alf.tap_chroma;
  N   = (tap*tap+1)>>1;

  for(i=0; i<N; i++)
  {
    alfps->alf.coeff_chroma[i] -= alf_coef_base_chroma[base_type][i];
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    set base of filter coefficients for chroma
 *****************************************************************************************
 */
void SetALFBaseCoefChroma(alf_parameter_set_rbsp_t *alfps)
{
  int i, N, tap;

  tap = alfps->alf.tap_chroma;
  N   = (tap*tap+1)>>1;

  for(i=0; i<N; i++)
  {
    alf_coef_base_chroma[base_type][i] += alfps->alf.coeff_chroma[i];
  }
}


/*!
 *****************************************************************************************
 * \brief
 *    calculation lambda for adaptive loop filter
 *****************************************************************************************
 */
void CalcALFLambda(int qp)
{
  if(img->type == B_SLICE && img->nal_reference_idc == 0)
  {
    lambda_alf_luma   = 0.68 * pow (2, (qp - SHIFT_QP) / 3.0) * 3;
    lambda_alf_chroma = lambda_alf_luma;
  }
  else if(img->type == B_SLICE && img->nal_reference_idc == 1)
  {
    lambda_alf_luma   = 0.68 * pow (2, (qp - SHIFT_QP) / 3.0) * 2;
    lambda_alf_chroma = lambda_alf_luma;
  }
  else if(img->type == P_SLICE)
  {
    lambda_alf_luma   = 0.68 * pow (2, (qp - SHIFT_QP) / 3.0);
    lambda_alf_chroma = 0.68 * pow (2, (qp - SHIFT_QP) / 3.0) * 2;
  }
  else
  {
    lambda_alf_luma   = 0.68 * pow (2, (qp - SHIFT_QP) / 3.0);
    lambda_alf_chroma = lambda_alf_luma;
  }
}

/*!
 ************************************************************************
 * \brief
 *    find distortion for luma
 ************************************************************************
 */
int64 FindDistortionFrame(imgpel** ImgOrg, imgpel** ImgCmp)
{
  int i, j;
  int64  diff_y;

  diff_y = 0;
  for (j = 0; j < input->img_height; ++j)
  {
    for (i = 0; i < input->img_width; ++i)
    {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
      diff_y += SQR_DEPTH(ImgOrg[j][i], ImgCmp[j][i], input->BitDepthLuma, img->BitDepthIncrease);
#else
      diff_y += img->quad[ImgOrg[j][i] - ImgCmp[j][i]];
#endif
    }
  }

  return diff_y;
}

/*!
 ************************************************************************
 * \brief
 *    find distortion for chroma
 ************************************************************************
 */
void FindDistortionFrameChroma(imgpel*** ImgUVOrg, imgpel*** ImgUVCmp, int64 *dist_u, int64 *dist_v)
{
  int i, j;
  int64 diff_u, diff_v;

  diff_u = diff_v = 0;
  if (img->yuv_format != YUV400)
  {
    for (j = 0; j < input->img_height_cr; j++)
    {
      for (i = 0; i < input->img_width_cr; i++)
      {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        diff_u += SQR_DEPTH(ImgUVOrg[0][j][i], ImgUVCmp[0][j][i], input->BitDepthChroma, img->BitDepthIncreaseChroma);
        diff_v += SQR_DEPTH(ImgUVOrg[1][j][i], ImgUVCmp[1][j][i], input->BitDepthChroma, img->BitDepthIncreaseChroma);
#else
        diff_u += img->quad[ImgUVOrg[0][j][i] - ImgUVCmp[0][j][i]];
        diff_v += img->quad[ImgUVOrg[1][j][i] - ImgUVCmp[1][j][i]];
#endif
      }
    }
  }

  *dist_u = diff_u;
  *dist_v = diff_v;
}

/*!
 ************************************************************************
 * \brief
 *    find distortion for block
 ************************************************************************
 */
int64 FindDistortionBlock(int blk_addr, imgpel** ImgOrg, imgpel** ImgCmp, int blk_size, int width_in_alf_blks)
{
  int x, y, blk_x, blk_y, max_blk_x, max_blk_y;
  int64  diff_y;

  GetSamplePosForBALF(blk_addr, &blk_x, &blk_y, blk_size, width_in_alf_blks);

  max_blk_x = blk_x + blk_size;
  max_blk_y = blk_y + blk_size;
  if(max_blk_x > input->img_width)  max_blk_x = input->img_width;
  if(max_blk_y > input->img_height) max_blk_y = input->img_height;

  diff_y = 0;
  for (y = blk_y; y < max_blk_y; y++)
  {
    for (x = blk_x; x < max_blk_x; x++)
    {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
      diff_y += SQR_DEPTH(ImgOrg[y][x], ImgCmp[y][x], input->BitDepthLuma, img->BitDepthIncrease);
#else
      diff_y += img->quad[ImgOrg[y][x] - ImgCmp[y][x]];
#endif
    }
  }

  return diff_y;
}

/*!
 *****************************************************************************************
 * \brief
 *    calculation Rate-Distortion cost
 *****************************************************************************************
 */
void CalcRDCost(imgpel** ImgOrg, imgpel** ImgCmp, alf_parameter_set_rbsp_t *alfps, int64* rate, int64* dist, double* cost)
{
  int *tmp_coef;

  if(alfps != NULL)
  {
    get_mem1Dint    (&(tmp_coef), MAX_NUM_COEF);
    memcpy(tmp_coef, alfps->alf.coeff, sizeof(int)*alfps->alf.num_coeff);
    
#ifdef ALF_SPATIAL_PREDICT_COEF
    if(alfps->pred_coef_mode)
    {
      PredictALFCoeff(alfps);
    }
    else
#endif
    SetALFDiffCoef(alfps);

    if(alfps->block_control_flag && alfps->qt_partition_flag)
    {
      *rate  = send_ALFParameterSet(alfps, tempStream, ImgOrg, enc_picture->imgY, imgY_best);
    }
    else
    {
      *rate  = send_ALFParameterSet(alfps, tempStream, ImgOrg, enc_picture->imgY, ImgCmp);
    }

    memcpy(alfps->alf.coeff, tmp_coef, sizeof(int)*alfps->alf.num_coeff);
    free_mem1Dint(tmp_coef);
  }
  else
  {
    *rate    = 0;
  }
  *dist      = FindDistortionFrame(ImgOrg, ImgCmp);
  *cost      = (double)(*rate) * lambda_alf_luma + (double)(*dist);
}

/*!
 *****************************************************************************************
 * \brief
 *    calculation Rate-Distortion cost for Chroma
 *****************************************************************************************
 */
void CalcRDCostChroma(imgpel*** ImgUVOrg, imgpel*** ImgUVCmp, alf_parameter_set_rbsp_t *alfps, int64* rate, int64* dist, double* cost)
{
  int *tmp_coef;
  int64 dist_u, dist_v;

  if(alfps->chroma_idc)
  {
    get_mem1Dint    (&(tmp_coef), MAX_NUM_COEF_C);
    memcpy(tmp_coef, alfps->alf.coeff_chroma, sizeof(int)*alfps->alf.num_coeff_chroma);
    

#ifdef ALF_SPATIAL_PREDICT_COEF
    if(alfps->pred_coef_mode)
    {
      PredictALFCoeffChroma(alfps);
    }
    else
#endif
    SetALFDiffCoefChroma(alfps);

    if(alfps->block_control_flag && alfps->qt_partition_flag)
    {
      *rate  = send_ALFParameterSet(alfps, tempStream, imgY_org, enc_picture->imgY, imgY_best);
    }
    else
    {
      *rate  = send_ALFParameterSet(alfps, tempStream, imgY_org, enc_picture->imgY, imgY_rest);
    }

    memcpy(alfps->alf.coeff_chroma, tmp_coef, sizeof(int)*alfps->alf.num_coeff_chroma);
    free_mem1Dint(tmp_coef);
  }
  else
  {
    *rate    = tmp_luma_rate;
  }
  FindDistortionFrameChroma(ImgUVOrg, ImgUVCmp, &dist_u, &dist_v);
  *dist      = dist_u + dist_v;
  *cost      = (double)(*rate) * lambda_alf_chroma + (double)(*dist);
}

/*!
 *****************************************************************************************
 * \brief
 *    calculation Rate-Distortion cost for Block-based ALF
 *****************************************************************************************
 */
void CalcRDCostBlock(imgpel** ImgOrg, imgpel** ImgCmp, imgpel** ImgDec, alf_parameter_set_rbsp_t *alfps, int64* rate, int64* dist, double* cost)
{
  int i;
  int *tmp_coef;
  AlfBlock *currAlfBlk;

  if(alfps != NULL)
  {
    get_mem1Dint    (&(tmp_coef), MAX_NUM_COEF);
    memcpy(tmp_coef, alfps->alf.coeff, sizeof(int)*alfps->alf.num_coeff);


#ifdef ALF_SPATIAL_PREDICT_COEF
    if(alfps->pred_coef_mode)
    {
      PredictALFCoeff(alfps);
    }
    else
#endif
    SetALFDiffCoef(alfps);

    *rate   = send_ALFParameterSet(alfps, tempStream, imgY_org, enc_picture->imgY, ImgCmp);

    memcpy(alfps->alf.coeff, tmp_coef, sizeof(int)*alfps->alf.num_coeff);
    free_mem1Dint(tmp_coef);
  }
  else
  {
    *rate   = 0;
  }
  *dist     = 0;
  for (i=0; i<alfps->num_alf_block; i++)
  {
    currAlfBlk = &alfps->alf_blk_data[i];
    if(currAlfBlk->alf_blk_flag)
    {
      *dist += FindDistortionBlock(i, ImgOrg, ImgCmp, alfps->alf_block_size, alfps->WidthInAlfBlks);
    }
    else
    {
      *dist += FindDistortionBlock(i, ImgOrg, ImgDec, alfps->alf_block_size, alfps->WidthInAlfBlks);
    }
  }
  *cost     = (double)(*rate) * lambda_alf_luma + (double)(*dist);
}


/*!
 *****************************************************************************************
 * \brief
 *    calculates correlation
 *****************************************************************************************
 */
void CalcCorrelationFunc(imgpel** ImgOrg, imgpel** ImgDec, double** corr, int tap, int width, int height)
{
  int i, j, x, y, N, offset;
  imgpel *term;
  int *pFiltPos;
  int yy, xx;

  N      = (tap*tap+1)>>1;
  offset = tap>>1;

  switch(tap)
  {
  case 5:
    pFiltPos = symmetric_array5x5;
    break;
  case 7:
    pFiltPos = symmetric_array7x7;
    break;
  case 9:
    pFiltPos = symmetric_array9x9;
    break;
  default:
    printf("not support filter size : %d", tap);
    exit(1);
    break;
  }

  if((term = (imgpel *)calloc(N, sizeof(imgpel))) == NULL)
    no_mem_exit("CalcCorrelationFunc: term");

  for(y=ext_offset; y<height+ext_offset; y++)
  {
    for(x=ext_offset; x<width+ext_offset; x++)
    {
      i=0;
      memset(term, 0, sizeof(imgpel)*N);
      for(yy=y-offset; yy<=y+offset; yy++)
      {
        for(xx=x-offset; xx<=x+offset; xx++)
        {
          term[pFiltPos[i]] += ImgDec[yy][xx];
          i++;
        }
      }

      for(j=0; j<N; j++)
      {
        corr[j][j] += term[j]*term[j];
        for(i=j+1; i<N; i++)
          corr[j][i] += term[j]*term[i];

        // DC offset
        corr[j][N]   += term[j];
        corr[j][N+1] += ImgOrg[y-ext_offset][x-ext_offset]*term[j];
      }
      // DC offset
      for(i=0; i<N; i++)
        corr[N][i] += term[i];
      corr[N][N]   += 1;
      corr[N][N+1] += ImgOrg[y-ext_offset][x-ext_offset];
    }
  }
  for(j=0; j<N-1; j++)
    for(i=j+1; i<N; i++)
      corr[i][j] = corr[j][i];

  free(term);
}

/*!
 *****************************************************************************************
 * \brief
 *    Calculates stored correlation for Block-based ALF
 *****************************************************************************************
 */
void CalcStoredCorrelationFuncBlock(int blk_addr, imgpel** ImgOrg, imgpel** ImgDec, unsigned int** corr, int tap)
{
  int i, j, x, y, N, offset;
  int blk_size, blk_x, blk_y, max_blk_x, max_blk_y;
  unsigned int *term;
  int *pFiltPos;
  int yy, xx;

  N        = (tap*tap+1)>>1;
  offset   = tap>>1;
  blk_size = img->min_alf_block_size;

  switch(tap)
  {
  case 5:
    pFiltPos = symmetric_array5x5;
    break;
  case 7:
    pFiltPos = symmetric_array7x7;
    break;
  case 9:
    pFiltPos = symmetric_array9x9;
    break;
  default:
    printf("not support filter size : %d", tap);
    exit(1);
    break;
  }

  if((term = (unsigned int *)calloc(N, sizeof(unsigned int))) == NULL)
    no_mem_exit("CalcStoredCorrelationFuncBlock: term");

  GetSamplePosForBALF(blk_addr, &blk_x, &blk_y, blk_size, img->max_WidthInAlfBlks);

  max_blk_x = blk_x + blk_size;
  max_blk_y = blk_y + blk_size;
  if(max_blk_x > input->img_width)  max_blk_x = input->img_width;
  if(max_blk_y > input->img_height) max_blk_y = input->img_height;

  for(y=blk_y+ext_offset; y<max_blk_y+ext_offset; y++)
  {
    for(x=blk_x+ext_offset; x<max_blk_x+ext_offset; x++)
    {
      i=0;
      memset(term, 0, sizeof(unsigned int)*N);
      for(yy=y-offset; yy<=y+offset; yy++)
      {
        for(xx=x-offset; xx<=x+offset; xx++)
        {
          term[pFiltPos[i]] += (unsigned int)ImgDec[yy][xx];
          i++;
        }
      }

      for(j=0; j<N; j++)
      {
        corr[j][0] += term[j]*term[j];
        for(i=j+1; i<N; i++)
          corr[j][i-j] += term[j]*term[i];

        // DC offset
        corr[j][N-j]   += term[j];
        corr[j][N-j+1] += (unsigned int)ImgOrg[y-ext_offset][x-ext_offset]*term[j];
      }
      // DC offset
      corr[N][0] += 1;
      corr[N][1] += (unsigned int)ImgOrg[y-ext_offset][x-ext_offset];
    }
  }

  free(term);
}

/*!
 *****************************************************************************************
 * \brief
 *    Calculates correlation for Block-based ALF
 *****************************************************************************************
 */
void CalcCorrelationFuncBlock(int blk_addr, alf_parameter_set_rbsp_t *alfps, imgpel** ImgOrg, imgpel** ImgDec, double** corr)
{
  int i, j, x, y, tap, N, offset;
  int blk_size, blk_x, blk_y, max_blk_x, max_blk_y;
  imgpel *term;
  int *pFiltPos;
  int yy, xx;

  tap      = alfps->alf.tap;
  N        = (tap*tap+1)>>1;
  offset   = tap>>1;
  blk_size = alfps->alf_block_size;

  switch(tap)
  {
  case 5:
    pFiltPos = symmetric_array5x5;
    break;
  case 7:
    pFiltPos = symmetric_array7x7;
    break;
  case 9:
    pFiltPos = symmetric_array9x9;
    break;
  default:
    printf("not support filter size : %d", tap);
    exit(1);
    break;
  }

  if((term = (imgpel *)calloc(N, sizeof(imgpel))) == NULL)
    no_mem_exit("CalcCorrelationFuncBlock: term");

  GetSamplePosForBALF(blk_addr, &blk_x, &blk_y, blk_size, alfps->WidthInAlfBlks);

  max_blk_x = blk_x + blk_size;
  max_blk_y = blk_y + blk_size;
  if(max_blk_x > input->img_width)  max_blk_x = input->img_width;
  if(max_blk_y > input->img_height) max_blk_y = input->img_height;

  for(y=blk_y+ext_offset; y<max_blk_y+ext_offset; y++)
  {
    for(x=blk_x+ext_offset; x<max_blk_x+ext_offset; x++)
    {
      i=0;
      memset(term, 0, sizeof(imgpel)*N);
      for(yy=y-offset; yy<=y+offset; yy++)
      {
        for(xx=x-offset; xx<=x+offset; xx++)
        {
          term[pFiltPos[i]] += ImgDec[yy][xx];
          i++;
        }
      }

      for(j=0; j<N; j++)
      {
        corr[j][j] += (double)(term[j]*term[j]);
        for(i=j+1; i<N; i++)
          corr[j][i] += (double)(term[j]*term[i]);

        // DC offset
        corr[j][N]   += (double)term[j];
        corr[j][N+1] += (double)ImgOrg[y-ext_offset][x-ext_offset]*term[j];
      }
      // DC offset
      for(i=0; i<N; i++)
        corr[N][i] += (double)term[i];
      corr[N][N]   += 1;
      corr[N][N+1] += (double)ImgOrg[y-ext_offset][x-ext_offset];
    }
  }

  free(term);
}

/*!
 *****************************************************************************************
 * \brief
 *    gaussian elimination
 *****************************************************************************************
 */
int Gauss(double **a, int N)
{
  int i, j, k;
  double t;

  for(k=0; k<N-1; k++)
  {
    for(i=k+1;i<N; i++)
    {
      t=a[i][k]/a[k][k];
      for(j=k+1; j<=N; j++)
      {
        a[i][j] -= t * a[k][j];
        if(i==j && fabs(a[i][j])<0.000001) return 1;
      }
    }
  }
  for(i=N-1; i>=0; i--)
  {
    t = a[i][N];
    for(j=i+1; j<N; j++)
      t -= a[i][j] * a[j][N];
    a[i][N] = t / a[i][i];
  }
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *    quick sort filter coefficients
 ************************************************************************
 */
void FilterCoefQuickSort( double *coef_data, int *coef_num, int upper, int lower )
{
  double mid, tmp_data;
  int i, j, tmp_num;

  i = upper;
  j = lower;
  mid = coef_data[(lower+upper)>>1];
  do
  {
    while( coef_data[i] < mid ) i++;
    while( mid < coef_data[j] ) j--;
    if( i <= j )
    {
      tmp_data = coef_data[i];
      tmp_num  = coef_num[i];
      coef_data[i] = coef_data[j];
      coef_num[i]  = coef_num[j];
      coef_data[j] = tmp_data;
      coef_num[j]  = tmp_num;
      i++;
      j--;
    }
  } while( i <= j );
  if ( upper < j ) FilterCoefQuickSort(coef_data, coef_num, upper, j);
  if ( i < lower ) FilterCoefQuickSort(coef_data, coef_num, i, lower);
}

/*!
 *****************************************************************************************
 * \brief
 *    type transition from double to int and quantization 
 *****************************************************************************************
 */
void QuantFilterCoef(double* h, int* qh, int tap, int bit_depth)
{
  int i, N;
  int max_value, min_value;
  double dbl_total_gain;
  int total_gain, q_total_gain;
  int upper, lower;
  double *dh;
  int    *nc;
  int    *pFiltMag;

  switch(tap)
  {
  case 5:
    pFiltMag = symmetric_mag5x5;
    break;
  case 7:
    pFiltMag = symmetric_mag7x7;
    break;
  case 9:
    pFiltMag = symmetric_mag9x9;
    break;
  default:
    printf("not support filter size : %d", tap);
    exit(1);
    break;
  }

  N = (tap*tap+1)>>1;

  get_mem1Ddouble(&(dh), N);
  get_mem1Dint(&(nc), N);

  max_value =   (1<<(1+num_bit))-1;
  min_value = 0-(1<<(1+num_bit));

  dbl_total_gain=0.0;
  q_total_gain=0;
  for(i=0; i<N; i++)
  {
    if(h[i]>=0.0)
      qh[i] =  (int)( h[i]*(1<<num_bit)+0.5);
    else
      qh[i] = -(int)(-h[i]*(1<<num_bit)+0.5);

    dh[i] = (double)qh[i]/(double)(1<<num_bit) - h[i];
    dh[i]*=pFiltMag[i];
    dbl_total_gain += h[i]*pFiltMag[i];
    q_total_gain   += qh[i]*pFiltMag[i];
    nc[i] = i;
  }

  // modification of quantized filter coefficients
  total_gain = (int)(dbl_total_gain*(1<<num_bit)+0.5);

  if( q_total_gain != total_gain )
  {
    FilterCoefQuickSort(dh, nc, 0, N-1);
    if( q_total_gain > total_gain )
    {
      upper = N-1;
      while( q_total_gain > total_gain+1 )
      {
        i = nc[upper%N];
        qh[i]--;
        q_total_gain -= pFiltMag[i];
        upper--;
      }
      if( q_total_gain == total_gain+1 )
      {
        if(dh[N-1]>0)
          qh[N-1]--;
        else
        {
          i=nc[upper%N];
          qh[i]--;
          qh[N-1]++;
        }
      }
    }
    else if( q_total_gain < total_gain )
    {
      lower = 0;
      while( q_total_gain < total_gain-1 )
      {
        i=nc[lower%N];
        qh[i]++;
        q_total_gain += pFiltMag[i];
        lower++;
      }
      if( q_total_gain == total_gain-1 )
      {
        if(dh[N-1]<0)
          qh[N-1]++;
        else
        {
          i=nc[lower%N];
          qh[i]++;
          qh[N-1]--;
        }
      }
    }
  }
  
  // set of filter coefficients
  for(i=0; i<N; i++)
  {
    qh[i] = max(min_value,min(max_value, qh[i]));
  }

  // DC offset
  max_value = min(  (1<<(3+max(img->bitdepth_luma,img->bitdepth_chroma)))-1, (1<<14)-1);
  min_value = max( -(1<<(3+max(img->bitdepth_luma,img->bitdepth_chroma))),  -(1<<14)  );
  qh[N] =  (h[N]>=0.0)? (int)( h[N]*(1<<(num_bit-bit_depth+8)) + 0.5) : -(int)(-h[N]*(1<<(num_bit-bit_depth+8)) + 0.5);
  qh[N] = max(min_value,min(max_value, qh[N]));

  free_mem1Ddouble(dh);
  free_mem1Dint(nc);
}

/*!
 *****************************************************************************************
 * \brief
 *    clear filter coefficients
 *****************************************************************************************
 */
void ClearFilterCoefInt(int* qh, int N)
{
  int i;

  for(i=0; i<N; i++)
  {
    qh[i] = 0;
  }
  // center pos
  qh[N-2]  = 1<<num_bit;
}


//------- filtering for luma
/*!
 *****************************************************************************************
 * \brief
 *    adaptive loop filter main process
 *****************************************************************************************
 */
void AdaptiveLoopFilterProcess(Picture *pic)
{
  int tap, num_coef;
  Slice *currSlice;

  // set global variables
  tap         = MAX_NUM_TAP;
  num_coef    = (tap*tap+1)>>1;
  num_coef    = num_coef + 1; // DC offset
  num_bit     = NUM_BIT_SHIFT;
  r_offset    = 1<<(num_bit-1);
  ext_offset  = tap>>1;
#ifdef ALF_SPATIAL_PREDICT_COEF
  if(input->ALFPredCoefMode==0)
#endif
  base_type   = GetALFBaseType(img->type, img->nal_reference_idc);

  img->currentPicture = pic;
  img->currentSlice   = pic->slices[0];
  img->qp             = pic->slices[0]->qp;

  if(input->symbol_mode == CABAC)
  {
    img->model_number   = pic->slices[0]->model_number;
    temp_pic_bin_count = get_pic_bin_count();
    reset_pic_bin_count();
  }

  tempStream = NULL;

  // initialize base coefficients
#ifdef ALF_SPATIAL_PREDICT_COEF
  if(input->ALFPredCoefMode==0 && img->currentPicture->idr_flag)
#else
  if(img->currentPicture->idr_flag)
#endif
  {
    InitFilterBase();
  }

  // set lambda
  CalcALFLambda(img->qp);

  // extend image for filtering
  MakeExtendFrame( enc_picture->imgY, imgY_ext, input->img_width, input->img_height, ext_offset );

  // set min cost
  min_rate = INT_MAX;
  min_dist = INT_MAX;
  min_cost = DBL_MAX;

  // calc original cost
  CalcRDCost(imgY_org, enc_picture->imgY, NULL, &orig_rate, &orig_dist, &orig_cost);

  // initialize temp_alfps
  temp_alfps->alf_flag           = 1;
  temp_alfps->alf.tap            = tap;
  temp_alfps->alf.num_coeff      = num_coef;
  temp_alfps->chroma_idc         = 0;
  temp_alfps->block_control_flag = 0;
  temp_alfps->qt_partition_flag  = 0;
#ifdef ALF_SPATIAL_PREDICT_COEF
  temp_alfps->pred_coef_mode     = input->ALFPredCoefMode;
#endif

  // adaptive in-loop wiener filtering
  FirstFilteringFrameLuma(imgY_org, imgY_ext, imgY_rest);

  // block-based filter on/off control
  BlockAdaptiveFilterControl(imgY_org, imgY_ext, imgY_rest);

  // adaptive tap-length
  FilterTapDecision(imgY_org, imgY_ext, imgY_rest);

  CalcRDCost(imgY_org, imgY_rest, active_alfps, &rate, &dist, &cost);

  if( cost < orig_cost )
  {
    active_alfps->alf_flag = 1;
  }
  else
  {
    active_alfps->alf_flag = 0;
  }

  if(input->symbol_mode == CABAC)
  {
    reset_pic_bin_count();
    add_pic_bin_count(temp_pic_bin_count);
  }

  if(active_alfps->alf_flag)
  {
    // set alf_flag in slice header
    currSlice = img->currentSlice;
    *(currSlice->alf_flag_stream_pointer) += 1<<(currSlice->alf_flag_bit_to_go-1);

#ifdef ALF_SPATIAL_PREDICT_COEF
    // predict coeff
    if(active_alfps->pred_coef_mode)
    {
      PredictALFCoeff(active_alfps);
    }
    else
    {
#endif
    // set diff coefficients
    SetALFDiffCoef(active_alfps);

    // set base coefficients
    SetALFBaseCoef(active_alfps);
#ifdef ALF_SPATIAL_PREDICT_COEF
    }
#endif

    // adaptive in-loop wiener filtering for chroma
    tmp_luma_rate = rate;
    AdaptiveLoopFilterProcessChroma();

    // make bitstream for ALF
    if(active_alfps->block_control_flag && active_alfps->qt_partition_flag)
    {
      send_ALFParameterSet(active_alfps, alfstream, imgY_org, enc_picture->imgY, imgY_best);
    }
    else
    {
      send_ALFParameterSet(active_alfps, alfstream, imgY_org, enc_picture->imgY, imgY_rest);
    }

    // set filtered frame
    SetRestFrame(imgY_rest, enc_picture->imgY, input->img_width, input->img_height);
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    first filtering frame for luma
 *****************************************************************************************
 */
void FirstFilteringFrameLuma(imgpel** ImgOrg, imgpel** ImgDec, imgpel** ImgRest)
{
  int tap, num_coef;

  tap                 = MIN_NUM_TAP;
  temp_alfps->alf.tap = tap;
  num_coef            = (tap*tap+1)>>1;
  num_coef            = num_coef + 1; // DC offset
  temp_alfps->alf.num_coeff = num_coef;

  FilteringFrameLuma(ImgOrg, ImgDec, ImgRest, 1);

  CalcRDCost(ImgOrg, ImgRest, temp_alfps, &rate, &dist, &cost);

  if( cost < min_cost )
  {
    min_rate = rate;
    min_dist = dist;
    min_cost = cost;
    SetALFParameters(active_alfps, temp_alfps);
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    filtering frame for luma
 *****************************************************************************************
 */
void FilteringFrameLuma(imgpel** ImgOrg, imgpel** ImgDec, imgpel** ImgRest, int store_corr)
{
  int    i, tap, N, err_code;
  double **corr, *h;
  int    *qh;
  int    j, k, bx, by;
  unsigned int  **p_blk_corr;

  tap  = temp_alfps->alf.tap;
  N    = temp_alfps->alf.num_coeff;
  corr = alf_corr;
  h    = dbl_alf_coef;
  qh   = temp_alfps->alf.coeff;

  // initialize correlation
  for(i=0; i<N; i++)
    memset(corr[i], 0, sizeof(double)*(N+1));

  if(store_corr)
  {
    // store correlation per minimum size block
    for (i=0; i<img->max_num_alf_block; i++)
    {
      by = i / img->max_WidthInAlfBlks;
      bx = i % img->max_WidthInAlfBlks;

      p_blk_corr = blk_corr[by][bx];

      for(j=0; j<N; j++)
        memset(p_blk_corr[j], 0, sizeof(unsigned int)*(N+1-j));

      CalcStoredCorrelationFuncBlock(i, ImgOrg, ImgDec, p_blk_corr, tap);

      for(j=0; j<N; j++)
        for(k=j; k<N+1; k++)
          corr[j][k] += (double)p_blk_corr[j][k-j];
    }
    for(j=0; j<N-1; j++)
      for(k=j+1; k<N; k++)
        corr[k][j] = corr[j][k];
  }
  else
  {
    CalcCorrelationFunc(ImgOrg, ImgDec, corr, tap, input->img_width, input->img_height);
  }

  err_code = Gauss(corr, N);

  if(err_code)
  {
    ClearFilterCoefInt(qh, N);
  }
  else
  {
    for(i=0; i<N; i++)
      h[i] = corr[i][N];

    QuantFilterCoef(h, qh, tap, img->bitdepth_luma);
  }

  FilteringProcessLuma(ImgDec, ImgRest, qh, tap);
}

/*!
 *****************************************************************************************
 * \brief
 *    filtering process for luma
 *****************************************************************************************
 */
void FilteringProcessLuma(imgpel** ImgDec, imgpel** ImgRest, int *qh, int tap)
{
  int i, x, y, value, N, offset;
  int *pFiltPos, PixSum[MAX_NUM_COEF];
  int yy, xx;

  N      = (tap*tap+1)>>1;
  offset = tap>>1;

  switch(tap)
  {
  case 5:
    pFiltPos = symmetric_array5x5;
    break;
  case 7:
    pFiltPos = symmetric_array7x7;
    break;
  case 9:
    pFiltPos = symmetric_array9x9;
    break;
  default:
    printf("not support filter size : %d", tap);
    exit(1);
    break;
  }

  for(y=ext_offset; y<input->img_height+ext_offset; y++)
  {
    for(x=ext_offset; x<input->img_width+ext_offset; x++)
    {
      i=0;
      memset(PixSum, 0, sizeof(int)*N);
      for(yy=y-offset; yy<=y+offset; yy++)
      {
        for(xx=x-offset; xx<=x+offset; xx++)
        {
          PixSum[pFiltPos[i]] += (int)ImgDec[yy][xx];
          i++;
        }
      }
      value = 0;
      for(i=0; i<N; i++)
      {
        value += qh[i]*PixSum[i];
      }
      // DC offset
      value += qh[N]<<(img->bitdepth_luma-8);
      value = (value + r_offset)>>num_bit;
      ImgRest[y-ext_offset][x-ext_offset] = (imgpel)Clip1(value);
    }
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    decision of filter-tap by Rate-Distortion Optimization
 *****************************************************************************************
 */
void FilterTapDecision(imgpel** ImgOrg, imgpel** ImgDec, imgpel** ImgRest)
{
  int    tap, N, not_balf_flag;
  int    pblk;

  // restriction for non-stored B-slice
  if(img->type==B_SLICE && img->nal_reference_idc == 0)
  {
    return;
  }

  for(tap=MIN_NUM_TAP+2; tap<=MAX_NUM_TAP; tap+=2)
  {
    SetALFParameters(temp_alfps, active_alfps);

    temp_alfps->alf.tap = tap;
    N = (tap*tap+1)>>1;
    N = N + 1; // DC offset
    temp_alfps->alf.num_coeff = N;

    if(temp_alfps->block_control_flag)
    {
      ReDesignFilterCoeff(ImgOrg, ImgDec, 0);
      FilteringProcessLuma(ImgDec, imgY_temp, temp_alfps->alf.coeff, temp_alfps->alf.tap);

      // set block-based on/off control flag
      if(temp_alfps->qt_partition_flag)
      {
        not_balf_flag = 1;
        for (pblk=0; pblk<temp_alfps->num_parent_block; pblk++)
        {
          SetBlockPartitioningFlag(pblk, temp_alfps, ImgOrg, enc_picture->imgY, imgY_temp);
          if((qt_blk_data->part_flag[0]==TRUE) || (qt_blk_data->alf_blk_flag[0]==TRUE))
          {
            not_balf_flag = 0;
          }
          TransQbpToFbp(pblk, temp_alfps);
        }
      }
      else
      {
        not_balf_flag = SetALFBlockFlag(temp_alfps, ImgOrg, enc_picture->imgY, imgY_temp);
      }
      if(not_balf_flag) continue;

      CalcRDCostBlock(ImgOrg, imgY_temp, enc_picture->imgY, temp_alfps, &rate, &dist, &cost);
    }
    else
    {
      FilteringFrameLuma(ImgOrg, ImgDec, imgY_temp, 0);
      CalcRDCost(ImgOrg, imgY_temp, temp_alfps, &rate, &dist, &cost);
    }

    if( cost < min_cost )
    {
      min_rate = rate;
      min_dist = dist;
      min_cost = cost;
      SetALFParameters(active_alfps, temp_alfps);
      SetRestFrame(imgY_temp, imgY_best, input->img_width, input->img_height);
    }
  }

  if(active_alfps->alf.tap > MIN_NUM_TAP)
  {
    SetRestFrame(imgY_best, ImgRest, input->img_width, input->img_height);
    if(active_alfps->block_control_flag)
    {
      CopyDecToRestBlock(active_alfps, enc_picture->imgY, ImgRest);
      SetBlockSizeParameter(active_alfps, active_alfps->alf_block_size_idx);
    }
  }

  // set temp_alfps
  SetALFParameters(temp_alfps, active_alfps);
}


//------- filtering for chroma
/*!
 *****************************************************************************************
 * \brief
 *    adaptive loop filter main process for chroma
 *****************************************************************************************
 */
void AdaptiveLoopFilterProcessChroma()
{
  int tap, num_coef;
  int64 dist_u, dist_v, dist_u_rest, dist_v_rest;

  // restriction for non-stored B-slice
  if(img->nal_reference_idc == 0)
  {
    active_alfps->chroma_idc = 0;
    return;
  }

  // set global variables for chroma
  tap        = MAX_NUM_TAP_C;
  ext_offset = tap>>1;
  num_coef   = (tap*tap+1)>>1;
  num_coef   = num_coef + 1;  // DC offset
#ifdef ALF_SPATIAL_PREDICT_COEF
  if(input->ALFPredCoefMode==0)
#endif
  base_type  = GetALFBaseType(img->type, img->nal_reference_idc);

  // extend picture for filtering
  MakeExtendFrame( enc_picture->imgUV[0], imgUV_ext[0], input->img_width_cr, input->img_height_cr, ext_offset );
  MakeExtendFrame( enc_picture->imgUV[1], imgUV_ext[1], input->img_width_cr, input->img_height_cr, ext_offset );
  
  // calc. original cost
  SetALFParameters(temp_alfps, active_alfps);
  CalcRDCostChroma(imgUV_org, enc_picture->imgUV, temp_alfps, &min_rate, &min_dist, &min_cost);

  // initialize temp_alfps
  temp_alfps->chroma_idc           = 3;
  temp_alfps->alf.tap_chroma       = tap;
  temp_alfps->alf.num_coeff_chroma = num_coef;

  // Adaptive in-loop wiener filtering for chroma
  FilteringFrameChroma(imgUV_org, imgUV_ext, imgUV_rest);

  // filter on/off decision for chroma
  FindDistortionFrameChroma( imgUV_org, enc_picture->imgUV, &dist_u,      &dist_v      );
  FindDistortionFrameChroma( imgUV_org, imgUV_rest,         &dist_u_rest, &dist_v_rest );

  temp_alfps->chroma_idc = 0;
  if(dist_u > dist_u_rest)
    temp_alfps->chroma_idc += 2;
  if(dist_v  > dist_v_rest )
    temp_alfps->chroma_idc += 1;

  if(temp_alfps->chroma_idc)
  {
    if(temp_alfps->chroma_idc!=3)
    {
      // chroma filter re-design
      FilteringFrameChroma(imgUV_org, imgUV_ext, imgUV_rest);
    }

    CalcRDCostChroma(imgUV_org, imgUV_rest, temp_alfps, &rate, &dist, &cost);

    if( cost < min_cost )
    {
      SetALFParameters(active_alfps, temp_alfps);

#ifdef ALF_SPATIAL_PREDICT_COEF
      // predict coeff
      if(active_alfps->pred_coef_mode)
      {
        PredictALFCoeffChroma(active_alfps);
      }
      else
      {
#endif
      // set diff coefficients
      SetALFDiffCoefChroma(active_alfps);

      // set base coefficients
      SetALFBaseCoefChroma(active_alfps);
#ifdef ALF_SPATIAL_PREDICT_COEF
      }
#endif

      // set filtered frame
      if((active_alfps->chroma_idc>>1)&0x01)
        SetRestFrame(imgUV_rest[0], enc_picture->imgUV[0], input->img_width_cr, input->img_height_cr);
      if(active_alfps->chroma_idc&0x01)
        SetRestFrame(imgUV_rest[1], enc_picture->imgUV[1], input->img_width_cr, input->img_height_cr);
    }
    else
    {
      active_alfps->chroma_idc = 0;
    }
  }
  else
  {
    active_alfps->chroma_idc = 0;
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    filtering frame for chroma
 *****************************************************************************************
 */
void FilteringFrameChroma(imgpel ***ImgUVOrg, imgpel ***ImgUVDec, imgpel ***ImgUVRest)
{
  int i, tap, N, err_code;
  double **corr, *h;
  int    *qh;

  tap  = temp_alfps->alf.tap_chroma;
  N    = temp_alfps->alf.num_coeff_chroma;
  corr = alf_corr;
  h    = dbl_alf_coef;
  qh   = temp_alfps->alf.coeff_chroma;

  // initialize correlation
  for(i=0; i<N; i++)
    memset(corr[i], 0, sizeof(double)*(N+1));

  if((temp_alfps->chroma_idc>>1)&0x01)
  {
    CalcCorrelationFunc(ImgUVOrg[0], ImgUVDec[0], corr, tap, input->img_width_cr, input->img_height_cr);
  }
  if((temp_alfps->chroma_idc)&0x01)
  {
    CalcCorrelationFunc(ImgUVOrg[1], ImgUVDec[1], corr, tap, input->img_width_cr, input->img_height_cr);
  }

  err_code = Gauss(corr, N);

  if(err_code)
  {
    ClearFilterCoefInt(qh, N);
  }
  else
  {
    for(i=0; i<N; i++)
      h[i] = corr[i][N];

    QuantFilterCoef(h, qh, tap, img->bitdepth_chroma);
  }

  if((temp_alfps->chroma_idc>>1)&0x01)
  {
    FilteringProcessChroma(ImgUVDec, ImgUVRest, qh, tap, 0);
  }
  if((temp_alfps->chroma_idc)&0x01)
  {
    FilteringProcessChroma(ImgUVDec, ImgUVRest, qh, tap, 1);
  }

  if(temp_alfps->chroma_idc<3)
  {
    if(temp_alfps->chroma_idc==1)
    {
      SetRestFrame(enc_picture->imgUV[0], ImgUVRest[0], input->img_width_cr, input->img_height_cr);
    }
    if(temp_alfps->chroma_idc==2)
    {
      SetRestFrame(enc_picture->imgUV[1], ImgUVRest[1], input->img_width_cr, input->img_height_cr);
    }
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    filtering process for chroma
 *****************************************************************************************
 */
void FilteringProcessChroma(imgpel ***ImgUVDec, imgpel ***ImgUVRest, int *qh, int tap, int chroma_id)
{
  int i, x, y, value, N, offset;
  int *pFiltPos, PixSum[MAX_NUM_COEF_C];
  int yy, xx;

  N      = (tap*tap+1)>>1;
  offset = tap>>1;

  switch(tap)
  {
  case 5:
    pFiltPos = symmetric_array5x5;
    break;
  case 7:
    pFiltPos = symmetric_array7x7;
    break;
  case 9:
    pFiltPos = symmetric_array9x9;
    break;
  default:
    printf("not support filter size : %d", tap);
    exit(1);
    break;
  }

  for(y=ext_offset; y<input->img_height_cr+ext_offset; y++)
  {
    for(x=ext_offset; x<input->img_width_cr+ext_offset; x++)
    {
      i=0;
      memset(PixSum, 0, sizeof(int)*N);
      for(yy=y-offset; yy<=y+offset; yy++)
      {
        for(xx=x-offset; xx<=x+offset; xx++)
        {
          PixSum[pFiltPos[i]] += (int)ImgUVDec[chroma_id][yy][xx];
          i++;
        }
      }
      value = 0;
      for(i=0; i<N; i++)
      {
        value += qh[i]*PixSum[i];
      }
      // DC offset
      value += qh[N]<<(img->bitdepth_chroma-8);

      value = (value + r_offset)>>num_bit;
      ImgUVRest[chroma_id][y-ext_offset][x-ext_offset] = (imgpel)Clip1_Chr(value);
    }
  }
}


//------- filtering based block
/*!
 *****************************************************************************************
 * \brief
 *    block-based adaptive loop filtering
 *****************************************************************************************
 */
void BlockAdaptiveFilterControl(imgpel** ImgOrg, imgpel** ImgDec, imgpel** ImgRest)
{
  int   idx, min_idx, rd, not_balf_flag;
  int   num_of_redesign = NUM_OF_REDESIGN;
  int   level, max_level, pblk;
  int   i, l;
  int64 temp_min_dist, delta_dist;
  int   temp_block_th, block_th[8][5];

  temp_alfps->block_control_flag = 1;
  min_idx   = img->min_alf_block_size_idx;
  max_level = 4;

  // ET: set threshold
  temp_min_dist = min_dist;
  for(idx=0; idx<8; idx++)
  {
    for(level=0; level<=max_level; level++)
    {
      block_th[idx][level] = MAX_NUM_BLOCK;
    }
  }

  SetRestFrame(ImgRest, imgY_temp, input->img_width, input->img_height);

  // block size decision by Rate-Distortion Optimization
  for(idx=min_idx; idx<8; idx++)
  {
    // set block parameter
    SetBlockSizeParameter(temp_alfps, idx);

    // ET: block size restrictions
    if( temp_alfps->num_alf_block < MIN_NUM_BLOCK ) continue;

    for(level=0; level<=max_level; level++)
    {
      // set quadtree-based block parameter
      temp_alfps->max_layer_level = level;
      SetParentBlockSizeParam(temp_alfps, temp_alfps->alf_block_size);

      // ET: block size restrictions
      if( temp_alfps->num_parent_block > block_th[idx][level] ) continue;

      if(level>0)
      {
        if( temp_alfps->num_parent_block < MIN_NUM_PARENT_BLOCK ) continue;
        if( (temp_alfps->num_parent_block<<1) > block_th[idx][level] ) continue;
      }

      for(rd=0; rd<=num_of_redesign; rd++)
      {
        if(rd>0)
        {
          // re-design filter coefficients
          ReDesignFilterCoeff(ImgOrg, ImgDec, 1);
          FilteringProcessLuma(ImgDec, imgY_temp, temp_alfps->alf.coeff, temp_alfps->alf.tap);
        }

        if(level>0)
        {
          temp_alfps->qt_partition_flag = 1;

          // set quadtree-based on/off control flag
          not_balf_flag = 1;
          for (pblk=0; pblk<temp_alfps->num_parent_block; pblk++)
          {
            SetBlockPartitioningFlag(pblk, temp_alfps, ImgOrg, enc_picture->imgY, imgY_temp);
            if((qt_blk_data->part_flag[0]==TRUE) || (qt_blk_data->alf_blk_flag[0]==TRUE))
            {
              not_balf_flag = 0;
            }
            TransQbpToFbp(pblk, temp_alfps);
          }
        }
        else
        {
          temp_alfps->qt_partition_flag = 0;

          // set block-based on/off control flag
          not_balf_flag = SetALFBlockFlag(temp_alfps, ImgOrg, enc_picture->imgY, imgY_temp);
        }
        if(not_balf_flag) continue;

        CalcRDCostBlock(ImgOrg, imgY_temp, enc_picture->imgY, temp_alfps, &rate, &dist, &cost);

        // ET: set threshold
        delta_dist    = temp_min_dist - dist;
        temp_block_th = (int)((delta_dist/lambda_alf_luma)+0.5);
        if( (rd==0) || (temp_block_th>block_th[idx][level]) )
        {
          block_th[idx][level] = temp_block_th;
        }

        if( cost < min_cost )
        {
          min_rate = rate;
          min_dist = dist;
          min_cost = cost;
          SetALFParameters(active_alfps, temp_alfps);
          SetRestFrame(imgY_temp, imgY_best, input->img_width, input->img_height);
        }
      } // re-design loop

      // ET: set threshold
      for(i=idx; i<8; i++)
      {
        for(l=level; l<=max_level; l++)
        {
          if(block_th[idx][level]<block_th[i][l])
          {
            block_th[i][l] = block_th[idx][level];
          }
        }
      }

    } // layer loop
  } // block size loop

  if(active_alfps->block_control_flag)
  {
    SetRestFrame(imgY_best, ImgRest, input->img_width, input->img_height);
    CopyDecToRestBlock(active_alfps, enc_picture->imgY, ImgRest);
    SetBlockSizeParameter(active_alfps, active_alfps->alf_block_size_idx);
  }

  // set temp_alfps
  SetALFParameters(temp_alfps, active_alfps);
}

/*!
 *****************************************************************************************
 * \brief
 *    re-design filter coefficients
 *****************************************************************************************
 */
void ReDesignFilterCoeff(imgpel** ImgOrg, imgpel** ImgDec, int read_corr)
{  
  int    i, j, k, tap, N, err_code;
  int    x, y, bn, bx, by, max_bx, max_by;
  double **corr, *h;
  int    *qh;
  AlfBlock *currAlfBlk;

  tap  = temp_alfps->alf.tap;
  N    = temp_alfps->alf.num_coeff;
  corr = alf_corr;
  h    = dbl_alf_coef;
  qh   = temp_alfps->alf.coeff;
  bn   = temp_alfps->alf_block_size / img->min_alf_block_size;

  // initialize correlation
  for(i=0; i<N; i++)
    memset(corr[i], 0, sizeof(double)*(N+1));

  for (i=0; i<temp_alfps->num_alf_block; i++)
  {
    currAlfBlk = &temp_alfps->alf_blk_data[i];

    if(currAlfBlk->alf_blk_flag)
    {
      if(read_corr)
      {
        by = (i / temp_alfps->WidthInAlfBlks) * bn;
        bx = (i % temp_alfps->WidthInAlfBlks) * bn;
        max_by = by + bn;
        max_bx = bx + bn;
        if(max_by>img->max_HeightInAlfBlks) max_by = img->max_HeightInAlfBlks;
        if(max_bx>img->max_WidthInAlfBlks ) max_bx = img->max_WidthInAlfBlks;
        for(y=by; y<max_by; y++)
        {
          for(x=bx; x<max_bx; x++)
          {
            for(j=0; j<N; j++)
            {
              for(k=j; k<N+1; k++)
              {
                corr[j][k] += (double)blk_corr[y][x][j][k-j];
              }
            }
          }
        }
      }
      else
      {
        CalcCorrelationFuncBlock(i, temp_alfps, ImgOrg, ImgDec, corr);
      }
    }
  }
  for(j=0; j<N-1; j++)
    for(k=j+1; k<N; k++)
      corr[k][j] = corr[j][k];

  err_code = Gauss(corr,N);

  if(err_code)
  {
    ClearFilterCoefInt(qh, N);
  }
  else
  {
    for(i=0; i<N; i++)
      h[i] = corr[i][N];
    QuantFilterCoef(h, qh, tap, img->bitdepth_luma);
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    set parameters for block-based ALF
 *****************************************************************************************
 */
void SetBlockSizeParameter(alf_parameter_set_rbsp_t *alfps, int size_idx)
{
  alfps->alf_block_size_idx = size_idx;
  alfps->alf_block_size     = 1<<size_idx;
  alfps->WidthInAlfBlks     = input->img_width  / alfps->alf_block_size;
  alfps->HeightInAlfBlks    = input->img_height / alfps->alf_block_size;
  if( (input->img_width  % alfps->alf_block_size)!=0 ) alfps->WidthInAlfBlks++;
  if( (input->img_height % alfps->alf_block_size)!=0 ) alfps->HeightInAlfBlks++;
  alfps->num_alf_block      = alfps->WidthInAlfBlks * alfps->HeightInAlfBlks;

  // set image parameter for CABAC
  img->alf_block_size  = alfps->alf_block_size;
  img->WidthInAlfBlks  = alfps->WidthInAlfBlks;
  img->HeightInAlfBlks = alfps->HeightInAlfBlks;
  img->num_alf_block   = alfps->num_alf_block;
  img->alf_blk_data    = alfps->alf_blk_data;
}

/*!
 *****************************************************************************************
 * \brief
 *    calculate the number of blocks 
 *****************************************************************************************
 */
int CalcNumOfAlfBlocks(int blk_size, int width, int height)
{
  int width_in_blks, height_in_blks, num_of_blks;

  width_in_blks     = width  / blk_size;
  height_in_blks    = height / blk_size;
  if( (width  % blk_size)!=0 ) width_in_blks++;
  if( (height % blk_size)!=0 ) height_in_blks++;
  num_of_blks      = width_in_blks * height_in_blks;

  return num_of_blks;
}

/*!
 *****************************************************************************************
 * \brief
 *    set on/off flag per block
 *****************************************************************************************
 */
int SetALFBlockFlag(alf_parameter_set_rbsp_t *alfps, imgpel** ImgOrg, imgpel** ImgDec, imgpel** ImgRest)
{
  int i;
  int64 diff_dec, diff_rest;
  AlfBlock *currAlfBlk;
  int all_blk_flag_zero=1;

  for (i=0; i<alfps->num_alf_block; i++)
  {
    if(input->symbol_mode == CABAC)
    {
      img->current_alfb_nr = i;
      CheckAvailabilityOfNeighborsAlf();
    }

    currAlfBlk = &alfps->alf_blk_data[i];

    diff_dec  = FindDistortionBlock(i, ImgOrg, ImgDec,  alfps->alf_block_size, alfps->WidthInAlfBlks);
    diff_rest = FindDistortionBlock(i, ImgOrg, ImgRest, alfps->alf_block_size, alfps->WidthInAlfBlks);

    if(diff_rest < diff_dec)
    {
      all_blk_flag_zero = 0;
      currAlfBlk->alf_blk_flag = TRUE;
    }
    else
    {
      currAlfBlk->alf_blk_flag = FALSE;
    }
  }

  return all_blk_flag_zero;
}

/*!
 *****************************************************************************************
 * \brief
 *    copy from dec to rest for filter-off blocks
 *****************************************************************************************
 */
void CopyDecToRestBlock(alf_parameter_set_rbsp_t *alfps, imgpel** ImgDec, imgpel** ImgRest)
{
  int i;
  int x, y, blk_x, blk_y, max_blk_x, max_blk_y;
  AlfBlock *currAlfBlk;

  for(i=0; i<alfps->num_alf_block; i++)
  {
    currAlfBlk = &alfps->alf_blk_data[i];

    if(!currAlfBlk->alf_blk_flag)
    {
      GetSamplePosForBALF(i, &blk_x, &blk_y, alfps->alf_block_size, alfps->WidthInAlfBlks);

      max_blk_x = blk_x + alfps->alf_block_size;
      max_blk_y = blk_y + alfps->alf_block_size;
      if(max_blk_x > input->img_width)  max_blk_x = input->img_width;
      if(max_blk_y > input->img_height) max_blk_y = input->img_height;

      for(y=blk_y; y<max_blk_y; y++)
      {
        for(x=blk_x; x<max_blk_x; x++)
        {
          ImgRest[y][x] = ImgDec[y][x];
        }
      }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    returns the x and y block coordinates for a given AlfBlockAddress
 ************************************************************************
 */
void GetBlockPosForBALF(int blk_addr, int *x, int*y, int width_in_blks)
{
  *x = (blk_addr % width_in_blks);
  *y = (blk_addr / width_in_blks);
}

/*!
 ************************************************************************
 * \brief
 *    returns the x and y sample coordinates for a given AlfBlockAddress
 ************************************************************************
 */
void GetSamplePosForBALF(int blk_addr, int *x, int*y, int blk_size, int width_in_blks)
{
  GetBlockPosForBALF(blk_addr, x, y, width_in_blks);
  
  (*x) *= blk_size;
  (*y) *= blk_size;
}

/*!
 ************************************************************************
 * \brief
 *    Allocate memory array for blk_corr -> unsigned int array4D[height][width][rows][columns]
 *
 * \par Output:
 *    memory size in bytes
 ************************************************************************
 */
int get_mem_blk_corr(unsigned int *****blk_corr, int height, int width, int rows, int columns )
{
  int  i, j, k;
  int  num;       // num of coeff (column + column-1 +  column-2 + ... + 2)
  unsigned int ****ptr_4D, ***ptr_3D;

  num = rows*(columns+2)/2;

  if((*blk_corr = (unsigned int ****) calloc(sizeof(unsigned int***), height))==NULL)
    no_mem_exit("get_mem_blk_corr: blk_corr");

  for(i=0; i<height; i++)
  {
    ptr_4D = *blk_corr+i;
    if((*ptr_4D = (unsigned int ***) calloc(sizeof(unsigned int**), width))==NULL)
      no_mem_exit("get_mem_blk_corr: blk_corr");
    for(j=0; j<width; j++)
    {
      ptr_3D = *ptr_4D+j;
      if((*ptr_3D = (unsigned int **) calloc(sizeof(unsigned int*), rows))==NULL)
        no_mem_exit("get_mem_blk_corr: blk_corr");
      if(((*ptr_3D)[0] = (unsigned int *) calloc(sizeof(unsigned int), num))==NULL)
        no_mem_exit("get_mem_blk_corr: blk_corr");
      for(k=0; k<rows-1; k++)
        (*ptr_3D)[k+1] = (*ptr_3D)[k] + columns - k;
    }
  }

  return height*width*num*sizeof(unsigned int);
}

/*!
 ************************************************************************
 * \brief
 *    free 4D memory array
 *    which was alocated with get_mem_blk_corr()
 ************************************************************************
 */
void free_mem_blk_corr(unsigned int ****blk_corr, int height, int width )
{
  int  i, j;

  if (blk_corr)
  {
    for(j=0;j<height;j++)
    {
      if (blk_corr[j])
      {
        for (i=0;i<width;i++)
        { 
          if (blk_corr[j][i])
          {
            if (blk_corr[j][i][0]) 
              free (blk_corr[j][i][0]);
            else
              error ("free_mem_blk_corr: trying to free unused memory",100);

            free (blk_corr[j][i]);
          }
          else
          {
            error ("free_mem_blk_corr: trying to free unused memory",100);
          }
        }
        free (blk_corr[j]);
      }
      else
      {
        error ("free_mem_blk_corr: trying to free unused memory",100);
      }
    }
    free (blk_corr);
  } else
  {
    error ("free_mem_blk_corr: trying to free unused memory",100);
  }
}

/*!
 ************************************************************************
 * \brief
 *    returns 1 if the ALF block at the given address is available
 ************************************************************************
 */
int ALFBlockIsAvailable(int blk_addr)
{
  if ((blk_addr < 0) || (blk_addr > ((int)img->num_alf_block - 1)))
    return 0;
  
  return 1;
}

/*!
 ************************************************************************
 * \brief
 *    Checks the availability of neighboring ALF blocks of
 *    the current ALF block for prediction and context determination;
 ************************************************************************
 */
void CheckAvailabilityOfNeighborsAlf()
{
  const int alfb_nr  = img->current_alfb_nr;
  AlfBlock *currAlfB = &img->alf_blk_data[alfb_nr];

  // mark all neighbors as unavailable
  currAlfB->alfb_available_up   = NULL;
  currAlfB->alfb_available_left = NULL;

  currAlfB->alfbAddrA = alfb_nr - 1;
  currAlfB->alfbAddrB = alfb_nr - img->WidthInAlfBlks;
  currAlfB->alfbAddrC = alfb_nr - img->WidthInAlfBlks + 1;
  currAlfB->alfbAddrD = alfb_nr - img->WidthInAlfBlks - 1;

  currAlfB->alfbAvailA = ALFBlockIsAvailable(currAlfB->alfbAddrA) && (( alfb_nr    % img->WidthInAlfBlks)!=0);
  currAlfB->alfbAvailB = ALFBlockIsAvailable(currAlfB->alfbAddrB);
  currAlfB->alfbAvailC = ALFBlockIsAvailable(currAlfB->alfbAddrC) && (((alfb_nr+1) % img->WidthInAlfBlks)!=0);
  currAlfB->alfbAvailD = ALFBlockIsAvailable(currAlfB->alfbAddrD) && (( alfb_nr    % img->WidthInAlfBlks)!=0);

  if (currAlfB->alfbAvailA) currAlfB->alfb_available_left = &(img->alf_blk_data[currAlfB->alfbAddrA]);
  if (currAlfB->alfbAvailB) currAlfB->alfb_available_up   = &(img->alf_blk_data[currAlfB->alfbAddrB]);
}

/*!
 ************************************************************************
 * \brief
 *    get neighbouring positions for non-aff coding
 * \param curr_alfb_nr
 *   current ALF block number (decoding order)
 * \param xN
 *    input x position
 * \param yN
 *    input y position
 * \param pix
 *    returns position informations
 ************************************************************************
 */
void GetNonAffNeighbourAlfBlock(unsigned int curr_alfb_nr, int xN, int yN, PixelPos *pix)
{
  AlfBlock *currAlfB = &img->alf_blk_data[curr_alfb_nr];
  int maxW, maxH;

  maxW = img->alf_block_size;
  maxH = img->alf_block_size;

  if ((xN<0)&&(yN<0))
  {
    pix->mb_addr   = currAlfB->alfbAddrD;
    pix->available = currAlfB->alfbAvailD;
  }
  else if ((xN<0)&&((yN>=0)&&(yN<maxH)))
  {
    pix->mb_addr   = currAlfB->alfbAddrA;
    pix->available = currAlfB->alfbAvailA;
  }
  else if (((xN>=0)&&(xN<maxW))&&(yN<0))
  {
    pix->mb_addr   = currAlfB->alfbAddrB;
    pix->available = currAlfB->alfbAvailB;
  }
  else if (((xN>=0)&&(xN<maxW))&&((yN>=0)&&(yN<maxH)))
  {
    pix->mb_addr   = curr_alfb_nr;
    pix->available = 1;
  }
  else if ((xN>=maxW)&&(yN<0))
  {
    pix->mb_addr   = currAlfB->alfbAddrC;
    pix->available = currAlfB->alfbAvailC;
  }
  else 
  {
    pix->available = 0;
  }

  if (pix->available)
  {
    //pix->x = (xN + maxW) % maxW;
    //pix->y = (yN + maxH) % maxH;
    pix->x = (xN + maxW) & (maxW - 1);
    pix->y = (yN + maxH) & (maxH - 1);

    GetSamplePosForBALF(pix->mb_addr, &(pix->pos_x), &(pix->pos_y), img->alf_block_size, img->WidthInAlfBlks);

    pix->pos_x += pix->x;
    pix->pos_y += pix->y;
  }
}

/*!
 ************************************************************************
 * \brief
 *    get neighbouring positions. MB AFF is automatically used from img structure
 * \param curr_alfb_nr
 *   current ALF block number (decoding order)
 * \param xN
 *    input x position
 * \param yN
 *    input y position
 * \param pix
 *    returns position informations
 ************************************************************************
 */
void GetNeighbourAlfBlock(int curr_alfb_nr, int xN, int yN, PixelPos *pix)
{
  if (curr_alfb_nr<0)
    error ("GetNeighbourAlfBlock: invalid block number", 100);

  GetNonAffNeighbourAlfBlock(curr_alfb_nr, xN, yN, pix);
}

/*!
 *****************************************************************************************
 * \brief
 *    SetBlockPartitioningFlag
 *****************************************************************************************
 */
void SetBlockPartitioningFlag(int pblk_addr, alf_parameter_set_rbsp_t *alfps, imgpel** ImgOrg, imgpel** ImgDec, imgpel** ImgRest)
{
  int    x, y, sx, sy, ex, ey;
  int    pl, max_pl, bidx, pidx, pidx1, pidx2, pidx3, pidx4;
  int    pblk_size, pblk_x, pblk_y, pblk_width, num_pblock;
  int    sblk, sblk_size[5], sblk_x, sblk_y, sblk_width[5], num_sblock[5];
  int64  diff_dec, diff_rest;
  int64  d, d0, d1, d2, d3, d_sub;
  int64  r, r0, r1, r2, r3, r_sub;
  double c, c_sub;
  int64  dist_dec[PYRAMIDAL_MEM_SIZE];
  int64  dist_rest[PYRAMIDAL_MEM_SIZE];
  int64  dist_best[PYRAMIDAL_MEM_SIZE];
  int64  rate[PYRAMIDAL_MEM_SIZE];
  
  pblk_size   = alfps->parent_block_size;
  pblk_width  = alfps->WidthInParentBlks;
  num_pblock  = alfps->num_parent_block;
  max_pl      = alfps->max_layer_level;

  GetSamplePosForBALF (pblk_addr, &pblk_x, &pblk_y, pblk_size, pblk_width);

  // initialize partition block data
  memset(qt_blk_data->valid_flag,   FALSE, sizeof(Boolean)*PYRAMIDAL_MEM_SIZE);
  memset(qt_blk_data->part_flag,    FALSE, sizeof(Boolean)*PYRAMIDAL_MEM_SIZE);
  memset(qt_blk_data->alf_blk_flag, FALSE, sizeof(Boolean)*PYRAMIDAL_MEM_SIZE);
  memset(dist_dec,  0, sizeof(int64)*PYRAMIDAL_MEM_SIZE);
  memset(dist_rest, 0, sizeof(int64)*PYRAMIDAL_MEM_SIZE);
  memset(dist_best, 0, sizeof(int64)*PYRAMIDAL_MEM_SIZE);
  memset(rate,      0, sizeof(int64)*PYRAMIDAL_MEM_SIZE);

  // calc. distortion
  for (pl=max_pl; pl>=0; pl--)
  {
    sblk_width[pl] = 1<<pl;
    sblk_size[pl]  = pblk_size>>pl;
    num_sblock[pl] = 1<<(pl<<1);

    for (sblk=0; sblk<num_sblock[pl]; sblk++)
    {
      bidx = block_scan_order[pl][sblk];
      pidx = pyramid_mem_idx[pl][bidx];

      GetSamplePosForBALF (sblk, &sblk_x, &sblk_y, sblk_size[pl], sblk_width[pl]);

      sx = pblk_x + sblk_x;
      sy = pblk_y + sblk_y;
      ex = sx + sblk_size[pl];
      ey = sy + sblk_size[pl];
      if(sx>=input->img_width || sy>=input->img_height) continue;

      if(pl==max_pl)
      {
        if(ex > input->img_width)  ex = input->img_width;
        if(ey > input->img_height) ey = input->img_height;
        diff_dec = diff_rest = 0;
        for (y=sy; y<ey; y++)
        {
          for (x=sx; x<ex; x++)
          {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
            diff_dec  += SQR_DEPTH(ImgOrg[y][x], ImgDec [y][x], input->BitDepthLuma, img->BitDepthIncrease);
            diff_rest += SQR_DEPTH(ImgOrg[y][x], ImgRest[y][x], input->BitDepthLuma, img->BitDepthIncrease);
#else
            diff_dec  += img->quad[ImgOrg[y][x] - ImgDec [y][x]];
            diff_rest += img->quad[ImgOrg[y][x] - ImgRest[y][x]];
#endif
          }
        }
      }
      else
      {
        pidx1 = pyramid_mem_idx[pl+1][(bidx<<2)  ];
        pidx2 = pyramid_mem_idx[pl+1][(bidx<<2)+1];
        pidx3 = pyramid_mem_idx[pl+1][(bidx<<2)+2];
        pidx4 = pyramid_mem_idx[pl+1][(bidx<<2)+3];
        diff_dec  = dist_dec [pidx1] + dist_dec [pidx2] + dist_dec [pidx3] + dist_dec [pidx4];
        diff_rest = dist_rest[pidx1] + dist_rest[pidx2] + dist_rest[pidx3] + dist_rest[pidx4];
      }

      qt_blk_data->valid_flag[pidx] = TRUE;
      dist_dec [pidx] = diff_dec;
      dist_rest[pidx] = diff_rest;
      if(diff_rest<diff_dec)
      {
        qt_blk_data->alf_blk_flag[pidx] = TRUE;
        d = diff_rest;
      }
      else
      {
        qt_blk_data->alf_blk_flag[pidx] = FALSE;
        d = diff_dec;
      }
      dist_best[pidx] = d;
      rate[pidx]      = 1;
    }
  }

  // unification of sub-block data
  for( pl=max_pl-1; pl>=0; pl--)
  {
    for (sblk=0; sblk<num_sblock[pl]; sblk++)
    {
      bidx = block_scan_order[pl][sblk];
      pidx = pyramid_mem_idx[pl][bidx];

      if(qt_blk_data->valid_flag[pidx] == FALSE) continue;

      r  = rate[pidx] + 1;
      d  = dist_best[pidx];
      c  = (double)d + lambda_alf_luma * (double)r;

      pidx1 = pyramid_mem_idx[pl+1][(bidx<<2)  ];
      pidx2 = pyramid_mem_idx[pl+1][(bidx<<2)+1];
      pidx3 = pyramid_mem_idx[pl+1][(bidx<<2)+2];
      pidx4 = pyramid_mem_idx[pl+1][(bidx<<2)+3];
      r0    = rate[pidx1];
      r1    = rate[pidx2];
      r2    = rate[pidx3];
      r3    = rate[pidx4];
      d0    = dist_best[pidx1];
      d1    = dist_best[pidx2];
      d2    = dist_best[pidx3];
      d3    = dist_best[pidx4];
      
      r_sub = r0 + r1 + r2 + r3 + 1;
      d_sub = d0 + d1 + d2 + d3;
      c_sub = (double)d_sub + lambda_alf_luma * (double)r_sub;

      if ( c_sub < c )
      {
        rate[pidx]      = r_sub;
        dist_best[pidx] = d_sub;
        qt_blk_data->part_flag[pidx] = TRUE;
      }
      else
      {
        rate[pidx]      = r;
        dist_best[pidx] = d;
        qt_blk_data->part_flag[pidx] = FALSE;
      }
    }
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    translate quad-tree data into block data
 *****************************************************************************************
 */
void TransQbpToFbp(int pblk_addr, alf_parameter_set_rbsp_t *alfps)
{
  int pl, max_pl, bidx, pidx, pidx1, pidx2, pidx3, pidx4;
  int pblk_size, pblk_x, pblk_y, pblk_width, num_pblock;
  int sblk, sblk_size, sblk_x, sblk_y, sblk_width, num_sblock;
  int blk, blk_x, blk_y;
  AlfBlock  *currAlfBlk;

  pblk_size   = alfps->parent_block_size;
  pblk_width  = alfps->WidthInParentBlks;
  num_pblock  = alfps->num_parent_block;
  max_pl      = alfps->max_layer_level;

  GetBlockPosForBALF (pblk_addr, &pblk_x, &pblk_y, pblk_width);

  for (pl=0; pl<max_pl; pl++)
  {
    num_sblock = 1<<(pl<<1);

    for (bidx=0; bidx<num_sblock; bidx++)
    {
      pidx = pyramid_mem_idx[pl][bidx];
      if(qt_blk_data->part_flag[pidx]==FALSE)
      {
        pidx1 = pyramid_mem_idx[pl+1][(bidx<<2)  ];
        pidx2 = pyramid_mem_idx[pl+1][(bidx<<2)+1];
        pidx3 = pyramid_mem_idx[pl+1][(bidx<<2)+2];
        pidx4 = pyramid_mem_idx[pl+1][(bidx<<2)+3];
        qt_blk_data->part_flag[pidx1]    = FALSE;
        qt_blk_data->part_flag[pidx2]    = FALSE;
        qt_blk_data->part_flag[pidx3]    = FALSE;
        qt_blk_data->part_flag[pidx4]    = FALSE;
        qt_blk_data->alf_blk_flag[pidx1] = qt_blk_data->alf_blk_flag[pidx];
        qt_blk_data->alf_blk_flag[pidx2] = qt_blk_data->alf_blk_flag[pidx];
        qt_blk_data->alf_blk_flag[pidx3] = qt_blk_data->alf_blk_flag[pidx];
        qt_blk_data->alf_blk_flag[pidx4] = qt_blk_data->alf_blk_flag[pidx];
      }
    }
  }

  sblk_width = 1<<max_pl;
  sblk_size  = pblk_size>>max_pl;
  num_sblock = 1<<(max_pl<<1);

  for(sblk=0; sblk<num_sblock; sblk++)
  {
    bidx = block_scan_order[max_pl][sblk];
    pidx = pyramid_mem_idx[max_pl][bidx];

    GetBlockPosForBALF (sblk, &sblk_x, &sblk_y, sblk_width);

    blk_x = pblk_x * sblk_width + sblk_x;
    blk_y = pblk_y * sblk_width + sblk_y;
    if(blk_x >= alfps->WidthInAlfBlks)  continue;
    if(blk_y >= alfps->HeightInAlfBlks) continue;

    blk = blk_y * alfps->WidthInAlfBlks + blk_x;
    currAlfBlk = &alfps->alf_blk_data[blk];
    currAlfBlk->alf_blk_flag = qt_blk_data->alf_blk_flag[pidx];
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    set parent block parameters
 *****************************************************************************************
 */
void SetParentBlockSizeParam(alf_parameter_set_rbsp_t *alfps, int block_size)
{
  alfps->parent_block_size  = block_size<<alfps->max_layer_level;
  alfps->WidthInParentBlks  = input->img_width  / alfps->parent_block_size;
  alfps->HeightInParentBlks = input->img_height / alfps->parent_block_size;
  if( (input->img_width  % alfps->parent_block_size)!=0 ) alfps->WidthInParentBlks++;
  if( (input->img_height % alfps->parent_block_size)!=0 ) alfps->HeightInParentBlks++;
  alfps->num_parent_block   = alfps->WidthInParentBlks * alfps->HeightInParentBlks;
}


//------- send ALF parameters
/*!
 *****************************************************************************************
 * \brief
 *    send ALF parameters
 *****************************************************************************************
 */
int send_ALFParameterSet(alf_parameter_set_rbsp_t *alfps, Bitstream *bitstream, imgpel **ImgOrg, imgpel **ImgDec, imgpel **ImgRest)
{
  int new_bitstream_flag;
  int len = 0;

  new_bitstream_flag = 0;
  if (bitstream==NULL)
  {
    // allocate bitstream
    bitstream =  (Bitstream *) calloc(1, sizeof(Bitstream));
    bitstream->streamBuffer = (byte *) calloc(MAX_ALF_PARAM_SIZE, sizeof(byte));
    bitstream->bits_to_go=8;
    bitstream->byte_buf=0;
    bitstream->byte_pos=0;
    new_bitstream_flag = 1;
  }

  // ALF parameters
#ifdef ALF_SPATIAL_PREDICT_COEF
  len += u_1("ALFPS: pred_coef_mode", alfps->pred_coef_mode, bitstream);
#endif

  // filter parameters for luma
  len += send_ALFLumaParam(alfps, bitstream);

  // filter parameters for chroma
  len += send_ALFChromaParam(alfps, bitstream);

  // block control parameters for luma
  len += send_ALFBlockControlParam(alfps, bitstream, ImgOrg, ImgDec, ImgRest);

  // padding
  if (bitstream->bits_to_go != 8 )
  {
    while(bitstream->bits_to_go != 8)
    {
      len += u_1("ALFPS: padding bits", 1, bitstream); 
    }
  }

  if (new_bitstream_flag)
  {
    // free bitstream
    free(bitstream->streamBuffer);
    bitstream->streamBuffer=NULL;
    free(bitstream);
    bitstream=NULL;
  }

  return len;
}

/*!
 *****************************************************************************************
 * \brief
 *    send ALF parameters for luma
 *****************************************************************************************
 */
int send_ALFLumaParam(alf_parameter_set_rbsp_t *alfps, Bitstream *bitstream)
{
  int len = 0;
  int pos;
  char string[100];

  len += ue_v("ALFPS: alf_tap_size", (alfps->alf.tap-5)/2, bitstream);

  // filter coefficients for luma
  for(pos=0; pos<alfps->alf.num_coeff; pos++)
  {
    sprintf(string, "ALFPS: %d pos", pos);
    len += se_v(string, alfps->alf.coeff[pos], bitstream);
  }

  return len;
}

/*!
 *****************************************************************************************
 * \brief
 *    send ALF parameters for chroma
 *****************************************************************************************
 */
int send_ALFChromaParam(alf_parameter_set_rbsp_t *alfps, Bitstream *bitstream)
{
  int len = 0;
  int pos;
  char string[100];

  len += u_v(2, "ALFPS: alf_chroma_idc", alfps->chroma_idc, bitstream);

  if(alfps->chroma_idc)
  {
    len += ue_v("ALFPS: filter_tap_chroma", (alfps->alf.tap_chroma-5)/2, bitstream);

    // filter coefficients for chroma
    for(pos=0; pos<alfps->alf.num_coeff_chroma; pos++)
    {
      sprintf(string, "ALFPS: [chroma] %d pos", pos);
      len += se_v(string, alfps->alf.coeff_chroma[pos], bitstream);
    }
  }

  return len;
}

/*!
 *****************************************************************************************
 * \brief
 *    send ALF parameters for block-based
 *****************************************************************************************
 */
int send_ALFBlockControlParam(alf_parameter_set_rbsp_t *alfps, Bitstream *bitstream, imgpel **ImgOrg, imgpel **ImgDec, imgpel **ImgRest)
{
  int len = 0;

  len += u_1("ALFPS: block_control_flag", alfps->block_control_flag, bitstream);

  if(alfps->block_control_flag)
  {
    len += u_v(3, "ALFPS: alf_block_size_idx", alfps->alf_block_size_idx, bitstream);

    len += u_1("ALFPS: qt_partition_flag", alfps->qt_partition_flag, bitstream);
    if (alfps->qt_partition_flag)
    {
      len += send_QuadTreeBlockPartitioningParam(alfps, bitstream, ImgOrg, ImgDec, ImgRest);
    }
    else
    {
      // block-based on/off flag
      if (input->symbol_mode == CABAC)
      {
        reset_pic_bin_count();
        len += send_BlockFlagCABAC(alfps, bitstream);
        add_pic_bin_count(temp_pic_bin_count);
      }
      else
      {
        len += send_BlockFlagRunLength(alfps, bitstream);
      }
    }
  }

  return len;
}

/*!
 *****************************************************************************************
 * \brief
 *    send block flag by Run-length coding
 *****************************************************************************************
 */
int send_BlockFlagRunLength(alf_parameter_set_rbsp_t *alfps, Bitstream *bitstream)
{
  int i;
  int len = 0;
  int alf_blk_flag_run=0; 
  int prev_alf_blk_flag=0;
  AlfBlock *currAlfBlk;

  alf_blk_flag_run = 0;

  for (i=0; i<alfps->num_alf_block; i++)
  {
    currAlfBlk = &alfps->alf_blk_data[i];

    if( i==0 )
    {
      len += u_1("ALFPS: alf_blk_flag", currAlfBlk->alf_blk_flag, bitstream);
      alf_blk_flag_run = 0;
      prev_alf_blk_flag = currAlfBlk->alf_blk_flag;
    }
    else
    {
      if( (prev_alf_blk_flag==currAlfBlk->alf_blk_flag) && (i!=alfps->num_alf_block-1) )
      {
        alf_blk_flag_run++;
      }
      else
      {
        if( (prev_alf_blk_flag==currAlfBlk->alf_blk_flag) && (i==alfps->num_alf_block-1) )
        {
          alf_blk_flag_run++;
        }

        len += ue_v("ALFPS: alf_blk_flag_run", alf_blk_flag_run, bitstream);
        alf_blk_flag_run = 0;
        prev_alf_blk_flag = currAlfBlk->alf_blk_flag;
      }
    }
  }

  return len;
}

/*!
 *****************************************************************************************
 * \brief
 *    send block flag by CABAC
 *****************************************************************************************
 */
int send_BlockFlagCABAC(alf_parameter_set_rbsp_t *alfps, Bitstream *bitstream)
{
  int i;
  int len = 0;
  AlfBlock *currAlfBlk;
  SyntaxElement *currSE;
  DataPartition *dataPart;
  EncodingEnvironmentPtr eep=0;
  extern int cabac_encoding;
  Slice *currSlice = img->currentSlice;
#if TRACE
  char string[100];
#endif

  if( ( dataPart = (DataPartition *) calloc(1, sizeof(DataPartition)) ) == NULL )
    no_mem_exit ("calloc DataPartition");
  dataPart->bitstream = bitstream;
  dataPart->ee_cabac  = (currSlice->partArr[0]).ee_cabac;
  eep = &(dataPart->ee_cabac);
  writeVlcByteAlign(bitstream);
  arienco_start_encoding_AlfBlockFlag(eep, bitstream->streamBuffer, &(bitstream->byte_pos));
  init_contexts();
  dataPart->writeSyntaxElement = writeSyntaxElement_CABAC;
  cabac_encoding = 1;

  for (i=0; i<alfps->num_alf_block; i++)
  {
    currAlfBlk = &alfps->alf_blk_data[i];

    img->current_alfb_nr = i;
    currAlfBlk->currSEnr = 0;
    currSE = &img->MB_SyntaxElements[currAlfBlk->currSEnr];

    CheckAvailabilityOfNeighborsAlfCABAC(alfps);
    currSE->value1  = currAlfBlk->alf_blk_flag;
    currSE->type    = SE_HEADER;
    currSE->writing = writeAlfBlockFlagCABAC;
    dataPart->writeSyntaxElement(currSE, dataPart);
#if TRACE
    currSE->bitpattern = eep->Ebuffer;
    sprintf(string, "ALFPS: alf_blk_flag[%4d]", i);
    snprintf(currSE->tracestring, TRACESTRING_SIZE, string);
    trace2out(currSE);
#endif
    len += currSE->len;
    currAlfBlk->currSEnr++;

    // terminate the arithmetic code
    if (i == alfps->num_alf_block-1)
    {
      biari_encode_symbol_final(eep, 1); 
#if TRACE
      fprintf (p_trace, "      CABAC terminating bit = %d\n",1);
#endif
      arienco_done_encoding(eep);
    }
  }
  free(dataPart);
  cabac_encoding = 0;

  return len;
}

/*!
 *****************************************************************************************
 * \brief
 *    send quad-tree block partitioning parameters
 *****************************************************************************************
 */
int send_QuadTreeBlockPartitioningParam(alf_parameter_set_rbsp_t *alfps, Bitstream *bitstream, imgpel **ImgOrg, imgpel **ImgDec, imgpel **ImgRest)
{
  int len = 0;
  int pblk, layer_level, pidx;
  int layer, cblk, num_of_child_block, scblk;
  char string[100];

  layer_level = alfps->max_layer_level;

  len += u_v(2, "ALFPS: max_layer_level", alfps->max_layer_level-1, bitstream);

  for (pblk=0; pblk<alfps->num_parent_block; pblk++)
  {
    SetBlockPartitioningFlag(pblk, alfps, ImgOrg, ImgDec, ImgRest);

    assert(qt_blk_data->valid_flag[0]==TRUE);

    //partition_flag
    for(layer=0; layer<layer_level; layer++)
    {
      num_of_child_block = 1<<(layer<<1);
      for(cblk=0; cblk<num_of_child_block; cblk++)
      {
        pidx = pyramid_mem_idx[layer][cblk];
        if(qt_blk_data->valid_flag[pidx]==TRUE)
        {
          sprintf(string, "ALFPS: block_partitioning_flag[%d]", layer);
          len += u_1(string, qt_blk_data->part_flag[pidx], bitstream);
          if(qt_blk_data->part_flag[pidx]==FALSE)
          {
            for(scblk=0; scblk<4; scblk++)
            {
              pidx = pyramid_mem_idx[layer+1][(cblk<<2)+scblk];
              qt_blk_data->valid_flag[pidx] = FALSE;
            }
          }
        }
        else
        {
          for(scblk=0; scblk<4; scblk++)
          {
            pidx = pyramid_mem_idx[layer+1][(cblk<<2)+scblk];
            qt_blk_data->valid_flag[pidx] = FALSE;
          }
        }
      }
    }

    //onoff_flag
    for(layer=0; layer<=layer_level; layer++)
    {
      num_of_child_block = 1<<(layer<<1);
      for(cblk=0; cblk<num_of_child_block; cblk++)
      {
        pidx  = pyramid_mem_idx[layer][cblk];
        if(qt_blk_data->valid_flag[pidx]==TRUE && qt_blk_data->part_flag[pidx]==FALSE)
        {
          sprintf(string, "ALFPS: alf_blk_flag[%d]", layer);
          len += u_1(string, qt_blk_data->alf_blk_flag[pidx], bitstream);
        }
      }
    }
  }

  return len;
}

/*!
 *****************************************************************************************
 * \brief
 *    terminate slice for ALF
 *****************************************************************************************
 */
int TerminateSliceForALF(Picture* pic)
{
  static int MbWidthC  [4]= { 0, 8, 8,  16};
  static int MbHeightC [4]= { 0, 8, 16, 16};
  
  int bytes_written;
  Bitstream *currStream;
  Slice *currSlice;
  EncodingEnvironmentPtr eep;
  int i, j, slice;
  int byte_pos_before_startcode_emu_prevention;
  int min_num_bytes=0;
  int stuffing_bytes=0;
  int RawMbBits;

  int sh_end_pos, alf_end_pos, sd_end_pos;

  for (slice=0; slice<pic->no_slices; slice++)
  {
    currSlice   = pic->slices[slice];
    sh_end_pos  = currSlice->sh_byte_pos;
    alf_end_pos = sh_end_pos + alfstream->byte_pos;
    for (i=0; i<currSlice->max_part_nr; i++)
    {
      currStream = (currSlice->partArr[i]).bitstream;
      if(active_alfps->alf_flag)
      {
        sd_end_pos = currStream->byte_pos;
        for (j=0; j<sd_end_pos - sh_end_pos; j++)
        {
          currStream->streamBuffer[sd_end_pos+alfstream->byte_pos-j-1] = currStream->streamBuffer[sd_end_pos-j-1];
        }
        for (j=0; j<alfstream->byte_pos; j++)
        {
          currStream->streamBuffer[sh_end_pos+j] = alfstream->streamBuffer[j];
        }
        currStream->byte_pos = sd_end_pos + alfstream->byte_pos;
      }

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
        if (i==((currSlice->max_part_nr-1)))
        {
          RawMbBits = 256 * img->bitdepth_luma + 2 * MbWidthC[active_sps->chroma_format_idc] * MbHeightC[active_sps->chroma_format_idc] * img->bitdepth_chroma;
          min_num_bytes = ((96 * get_pic_bin_count()) - (RawMbBits * (int)img->PicSizeInMbs *3) + 1023) / 1024;
          if (min_num_bytes>img->bytes_in_picture)
          {
            stuffing_bytes = min_num_bytes - img->bytes_in_picture;
            printf ("CABAC stuffing words = %6d\n", stuffing_bytes/3);
          }
        }

        currStream->byte_pos = RBSPtoEBSP(currStream->streamBuffer, 0, currStream->byte_pos, currStream->byte_pos + stuffing_bytes);
        *(stats->em_prev_bits) += (currStream->byte_pos - byte_pos_before_startcode_emu_prevention) * 8;
      }           // CABAC
    }           // partition loop
  }

  alfstream->byte_pos  = 0;

  return 0;
}


#endif // ADAPTIVE_LOOP_FILTER

