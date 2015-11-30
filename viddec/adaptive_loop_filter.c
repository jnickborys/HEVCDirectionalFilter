
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

#include "adaptive_loop_filter.h"
#include "image.h"
#include "memalloc.h"
#include "vlc.h"
#include "cabac.h"
#include "biaridecod.h"
#include "context_ini.h"
#ifdef SWITCHED_FILTERS
#include "header.h"
#include "elements.h"
#include "fmo.h"
#endif

#ifdef ADAPTIVE_LOOP_FILTER

//------- global parameters for adaptive_loop_filter.c
int    num_bit;
int    r_offset;
int    ext_offset;
int    base_type;
int    alf_width;
int    alf_height;
int    alf_width_cr;
int    alf_height_cr;

imgpel **imgY_ext;
imgpel ***imgUV_ext;
int    **alf_coef_base;
int    ***alf_coef_tap_base;
int    **alf_coef_base_chroma;
AlfQtBlock *qt_blk_data;

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

int filter_pos5x5[14] = {
  20, 21, 22, 23, 24,
  29, 30, 31, 32, 33,
  38, 39, 40, 41
};
int filter_pos7x7[26] = {
  10, 11, 12, 13, 14, 15, 16,
  19, 20, 21, 22, 23, 24, 25,
  28, 29, 30, 31, 32, 33, 34,
  37, 38, 39, 40, 41
};
int filter_pos9x9[42] = {
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

#ifdef ALF_SPATIAL_PREDICT_COEF
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
#endif

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
*    Allocates memory for a alfps
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

  if(((p->alf_blk_data) = (AlfBlock *) calloc(img->num_alf_block, sizeof(AlfBlock))) == NULL)
    no_mem_exit("AllocALFPS: p->alf_blk_data");
  *memory_size += sizeof(AlfBlock)*img->num_alf_block;

  return p;
}

/*! 
*************************************************************************************
* \brief
*    Frees a alfps
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
int InitALFGlobalBuffers(void)
{
  int memory_size = 0;
  int idx, num_of_block, min_blk_size;

  memory_size += get_mem2Dpel (&(imgY_ext), alf_height+MAX_NUM_TAP, alf_width+MAX_NUM_TAP);
  memory_size += get_mem3Dpel (&(imgUV_ext), 2, alf_height_cr+MAX_NUM_TAP_C, alf_width_cr+MAX_NUM_TAP_C);
  memory_size += get_mem2Dint (&(alf_coef_base), 4, MAX_NUM_COEF);
  memory_size += get_mem3Dint (&(alf_coef_tap_base), 3, 4, MAX_NUM_COEF);
  memory_size += get_mem2Dint (&(alf_coef_base_chroma), 4, MAX_NUM_COEF_C);

  for(idx=0; idx<8; idx++)
  {
    num_of_block = CalcNumOfAlfBlocks((1<<idx), alf_width, alf_height);
    if( num_of_block<MAX_NUM_BLOCK ) break;
  }
  min_blk_size = 1<<idx;
  img->num_alf_block = CalcNumOfAlfBlocks(min_blk_size, alf_width, alf_height);

  // allocate ALF parameter set
  active_alfps = AllocALFPS(&memory_size);

  if( (qt_blk_data = (AlfQtBlock *) calloc(1, sizeof(AlfQtBlock))) == NULL )
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
  free_mem2Dpel (imgY_ext);
  free_mem3Dpel (imgUV_ext, 2);
  free_mem2Dint (alf_coef_base);
  free_mem3Dint (alf_coef_tap_base, 3);
  free_mem2Dint (alf_coef_base_chroma);

  FreeALFPS (active_alfps);

  if(qt_blk_data != NULL)
  {
    free(qt_blk_data);
    qt_blk_data=NULL;
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    set ALF parameters (alfps -> p)
 *****************************************************************************************
 */
void SetALFParameters(StorablePicture *p, alf_parameter_set_rbsp_t *alfps)
{
  p->alf_flag = alfps->alf_flag;

  if (alfps->alf_flag)
  {
#ifdef ALF_SPATIAL_PREDICT_COEF
    p->pred_coef_mode    = alfps->pred_coef_mode;
#endif
    p->alf.tap           = alfps->alf.tap;
    p->alf.num_coeff     = alfps->alf.num_coeff;
    memcpy(p->alf.coeff, alfps->alf.coeff, sizeof(int)*(p->alf.num_coeff));

    p->alf_chroma_idc = alfps->alf_chroma_idc;
    if (p->alf_chroma_idc)
    {
      p->alf.tap_chroma       = alfps->alf.tap_chroma;
      p->alf.num_coeff_chroma = alfps->alf.num_coeff_chroma;
      memcpy(p->alf.coeff_chroma, alfps->alf.coeff_chroma, sizeof(int)*(p->alf.num_coeff_chroma));
    }

    p->block_control_flag = alfps->block_control_flag;
    if (p->block_control_flag)
    {
      p->alf_block_size_idx = alfps->alf_block_size_idx;
      p->alf_block_size     = alfps->alf_block_size;
      p->WidthInAlfBlks     = alfps->WidthInAlfBlks;
      p->HeightInAlfBlks    = alfps->HeightInAlfBlks;
      p->num_alf_block      = alfps->num_alf_block;
      memcpy( p->alf_blk_data, alfps->alf_blk_data, sizeof(AlfBlock)*(p->num_alf_block));
      p->qt_partition_flag  = alfps->qt_partition_flag;
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

/*!
 *****************************************************************************************
 * \brief
 *    set image size for cropping
 *****************************************************************************************
 */
void SetALFImgSize(ImageParameters *img, seq_parameter_set_rbsp_t *p)
{
  int SubWidthC  [4]= { 1, 2, 2, 1};
  int SubHeightC [4]= { 1, 2, 1, 1};
  int crop_left, crop_right, crop_top, crop_bottom;

  // cropping for luma
  if (p->frame_cropping_flag)
  {
    crop_left   = SubWidthC[p->chroma_format_idc] * p->frame_cropping_rect_left_offset;
    crop_right  = SubWidthC[p->chroma_format_idc] * p->frame_cropping_rect_right_offset;
    crop_top    = SubHeightC[p->chroma_format_idc]*( 2 - p->frame_mbs_only_flag ) *  p->frame_cropping_rect_top_offset;
    crop_bottom = SubHeightC[p->chroma_format_idc]*( 2 - p->frame_mbs_only_flag ) *   p->frame_cropping_rect_bottom_offset;
  }
  else
  {
    crop_left = crop_right = crop_top = crop_bottom = 0;
  }

  alf_width  = img->width  - crop_left - crop_right;
  alf_height = img->height - crop_top - crop_bottom;

  // cropping for chroma
  if (p->frame_cropping_flag)
  {
    crop_left   = p->frame_cropping_rect_left_offset;
    crop_right  = p->frame_cropping_rect_right_offset;
    crop_top    = ( 2 - p->frame_mbs_only_flag ) *  p->frame_cropping_rect_top_offset;
    crop_bottom = ( 2 - p->frame_mbs_only_flag ) *   p->frame_cropping_rect_bottom_offset;
  }
  else
  {
    crop_left = crop_right = crop_top = crop_bottom = 0;
  }

  if ((p->chroma_format_idc==YUV400) && input->write_uv)
  {
    alf_width_cr  = img->width/2;
    alf_height_cr = img->height/2;
  }
  else
  {
    alf_width_cr  = img->width_cr  - crop_left - crop_right;
    alf_height_cr = img->height_cr - crop_top  - crop_bottom;
  }
}

#ifdef ALF_SPATIAL_PREDICT_COEF
 /*! 
 ************************************************************************
 * \brief
 *    PredictALFCoeff: 
 *    Predicts one luma filter coeff based on sum of other luma coeff
 ************************************************************************
 */
void PredictALFCoeff(StorablePicture *p)
{
  int i, sum=0, pred, N;
  int    *pFiltMag;

  switch(p->alf.tap)
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
    printf("not support filter size : %d", p->alf.tap);
    exit(1);
    break;
  }

  N   = (p->alf.tap*p->alf.tap+1)>>1;
  sum=0;
  for(i=0; i<N-1;i++)
  {
    sum+=pFiltMag[i]*p->alf.coeff[i];
  }
  pred=(1<<num_bit)-sum;
  p->alf.coeff[N-1]=pred-p->alf.coeff[N-1];
}

 /*! 
 ************************************************************************
 * \brief
 *    PredictALFCoeff: 
 *    Predicts one chroma filter coeff based on sum of other chroma coeff
 ************************************************************************
 */
void PredictALFCoeffChroma(StorablePicture *p)
{
  int i, pred;
  int N;
  int sum;
  int    *pFiltMag;

  switch(p->alf.tap_chroma)
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
    printf("not support filter size : %d", p->alf.tap_chroma);
    exit(1);
    break;
  }
  N   = (p->alf.tap_chroma*p->alf.tap_chroma+1)>>1;
  sum=0;
  for(i=0; i<N-1;i++)
  {
    sum+=pFiltMag[i]*p->alf.coeff_chroma[i];
  }
  pred=(1<<num_bit)-sum;
  p->alf.coeff_chroma[N-1]=pred-p->alf.coeff_chroma[N-1];
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
 *    set of adaptive loop filter coefficients 
 *****************************************************************************************
 */
void SetALFCoef(StorablePicture *p)
{
  int i, N, tap;
  int *filter_pos;
  int *base_flag;
  int **base_coef;

  tap = p->alf.tap;
  N   = (tap*tap+1)>>1;

  switch(p->alf.tap)
  {
  case 5:
    filter_pos = filter_pos5x5;
    base_flag  = filter_base_flag[0];
    base_coef  = alf_coef_tap_base[0];
    break;
  case 7:
    filter_pos = filter_pos7x7;
    base_flag  = filter_base_flag[1];
    base_coef  = alf_coef_tap_base[1];
    break;
  case 9:
    filter_pos = filter_pos9x9;
    base_flag  = filter_base_flag[2];
    base_coef  = alf_coef_tap_base[2];
    break;
  default:
    printf("not support filter size!!\n");
    exit(1);
    break;
  }

  if(base_flag[base_type] == 0)
  {
    for(i=0; i<N; i++)
    {
      p->alf.coeff[i] += alf_coef_base[base_type][filter_pos[i]];
      alf_coef_base[base_type][filter_pos[i]] = p->alf.coeff[i];
      base_coef[base_type][i]                 = p->alf.coeff[i];
    }
    base_flag[base_type] = 1;
  }
  else
  {
    for(i=0; i<N; i++)
    {
      p->alf.coeff[i] += base_coef[base_type][i];
      base_coef[base_type][i] = p->alf.coeff[i];
    }
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    set of adaptive loop-filter coefficients for chroma
 *****************************************************************************************
 */
void SetALFCoefChroma(StorablePicture *p)
{
  int i, tap, N;

  tap = p->alf.tap_chroma;
  N   = (tap*tap+1)>>1;

  for(i=0; i<N; i++)
  {
    p->alf.coeff_chroma[i] += alf_coef_base_chroma[base_type][i];
    alf_coef_base_chroma[base_type][i] = p->alf.coeff_chroma[i];
  }
}


//------- filtering for luma
/*!
 *****************************************************************************************
 * \brief
 *    adaptive loop filter main process
 *****************************************************************************************
 */
void AdaptiveLoopFilterProcess(StorablePicture *p)
{
  // initialize base coefficients
#ifdef ALF_SPATIAL_PREDICT_COEF
  if(p->pred_coef_mode==0 && p->idr_flag)
#else
  if(p->idr_flag)
#endif
  {
    InitFilterBase();
  }

  if(p->alf_flag)
  {
    // set global variables
    num_bit    = NUM_BIT_SHIFT;
    r_offset   = 1<<(num_bit-1);
    ext_offset = p->alf.tap>>1;
#ifdef ALF_SPATIAL_PREDICT_COEF
    if(p->pred_coef_mode==0)
#endif
    base_type  = GetALFBaseType(p->slice_type, p->used_for_reference);

    // extend image for filtering
    MakeExtendFrame( p->imgY, imgY_ext, alf_width, alf_height, ext_offset );

#ifdef ALF_SPATIAL_PREDICT_COEF
    // predict coeff
    if(p->pred_coef_mode)
    {
      PredictALFCoeff(p);
    }
    else
    {
#endif    
    // set filter coefficients
    SetALFCoef(p);
#ifdef ALF_SPATIAL_PREDICT_COEF
    }
#endif

    // Adaptive in-loop wiener filtering
    FilteringFrame( p, imgY_ext, p->imgY  );
    
    // Adaptive in-loop wiener filtering for chroma
    if(p->alf_chroma_idc)
    {
      AdaptiveLoopFilterProcessChroma(p);
    }
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    filtering frame for luma
 *****************************************************************************************
 */
void FilteringFrame(StorablePicture *p, imgpel** ImgDec, imgpel** ImgRest)
{
  if(p->block_control_flag)
  {
    BlockAdaptiveFilteringProcess( p, ImgDec, ImgRest );
  }
  else
  {
    FilteringProcess( ImgDec, ImgRest, p->alf.coeff, p->alf.tap );
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    filtering process
 *****************************************************************************************
 */
void FilteringProcess(imgpel** ImgDec, imgpel** ImgRest, int *qh, int tap)
{
  int i, x, y, value, N, offset;
  int *pFiltPos, PixSum[MAX_NUM_COEF];
  int yy, xx;

  N      = (tap*tap+1)>>1;
  offset = tap>>1;

  switch(tap)
  {
  case 9:
    pFiltPos = symmetric_array9x9;
    break;
  case 7:
    pFiltPos = symmetric_array7x7;
    break;
  case 5:
    pFiltPos = symmetric_array5x5;
    break;
  default:
    printf("not support filter size : %d", tap);
    exit(1);
    break;
  }

  for(y=ext_offset; y<alf_height+ext_offset; y++)
  {
    for(x=ext_offset; x<alf_width+ext_offset; x++)
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


//------- filtering for chroma
/*!
 *****************************************************************************************
 * \brief
 *    adaptive loop filter main process for chroma
 *****************************************************************************************
 */
void AdaptiveLoopFilterProcessChroma(StorablePicture *p)
{
  // set global variables for chroma
  ext_offset = p->alf.tap_chroma>>1;
#ifdef ALF_SPATIAL_PREDICT_COEF
  if(p->pred_coef_mode==0)
#endif
  base_type  = GetALFBaseType(p->slice_type, p->used_for_reference);

  // extend image for filtering
  MakeExtendFrame( p->imgUV[0], imgUV_ext[0], alf_width_cr, alf_height_cr, ext_offset );
  MakeExtendFrame( p->imgUV[1], imgUV_ext[1], alf_width_cr, alf_height_cr, ext_offset );

#ifdef ALF_SPATIAL_PREDICT_COEF
  // predict coeff
  if(p->pred_coef_mode)
  {
    PredictALFCoeffChroma(p);
  }
  else
  {
#endif
  // set filter coefficients for chroma
  SetALFCoefChroma(p);
#ifdef ALF_SPATIAL_PREDICT_COEF
  }
#endif

  // Adaptive in-loop wiener filtering for chroma
  FilteringFrameChroma( p, imgUV_ext, p->imgUV );
}

/*!
 *****************************************************************************************
 * \brief
 *    filtering frame for chroma
 *****************************************************************************************
 */
void FilteringFrameChroma(StorablePicture *p, imgpel ***ImgUVDec, imgpel ***ImgUVRest)
{
  if((p->alf_chroma_idc>>1)&0x01)
  {
    FilteringProcessChroma(ImgUVDec, ImgUVRest, p->alf.coeff_chroma, p->alf.tap_chroma, 0);
  }
  if(p->alf_chroma_idc&0x01)
  {
    FilteringProcessChroma(ImgUVDec, ImgUVRest, p->alf.coeff_chroma, p->alf.tap_chroma, 1);
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
  case 9:
    pFiltPos = symmetric_array9x9;
    break;
  case 7:
    pFiltPos = symmetric_array7x7;
    break;
  case 5:
    pFiltPos = symmetric_array5x5;
    break;
  default:
    printf("not support filter size : %d", tap);
    exit(1);
    break;
  }

  for(y=ext_offset; y<alf_height_cr+ext_offset; y++)
  {
    for(x=ext_offset; x<alf_width_cr+ext_offset; x++)
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
void BlockAdaptiveFilteringProcess(StorablePicture *p, imgpel** ImgDec, imgpel** ImgRest)
{
  int i;
  AlfBlock *currAlfBlk;

  for (i=0; i<p->num_alf_block; i++)
  { 
    currAlfBlk = &p->alf_blk_data[i];

    if(currAlfBlk->alf_blk_flag)
    {
      FilteringProcessBlock(i, ImgDec, ImgRest, p->alf.coeff, p->alf.tap, p->alf_block_size, p->WidthInAlfBlks);
    }
    else
    {
      CopyDecToRestBlock(i, p->imgY, ImgRest, p->alf_block_size, p->WidthInAlfBlks);
    }
  }
}
/*!
 *****************************************************************************************
 * \brief
 *    filtering process per block
 *****************************************************************************************
 */
void FilteringProcessBlock(int blk_addr, imgpel** ImgDec, imgpel** ImgRest, int *qh, int tap, int blk_size, int width_in_alf_blks)
{
  int i, value, offset, N;
  int x, y, blk_x, blk_y, max_blk_x, max_blk_y;
  int *pFiltPos, PixSum[MAX_NUM_COEF];
  int yy, xx;

  N      = (tap*tap+1)>>1;
  offset = tap>>1;

  switch(tap)
  {
  case 9:
    pFiltPos = symmetric_array9x9;
    break;
  case 7:
    pFiltPos = symmetric_array7x7;
    break;
  case 5:
    pFiltPos = symmetric_array5x5;
    break;
  default:
    printf("not support filter size : %d", tap);
    exit(1);
    break;
  }

  GetSamplePosForBALF(blk_addr, &blk_x, &blk_y, blk_size, width_in_alf_blks);

  max_blk_x = blk_x + blk_size;
  max_blk_y = blk_y + blk_size;
  if(max_blk_x > alf_width)  max_blk_x = alf_width;
  if(max_blk_y > alf_height) max_blk_y = alf_height;

  for(y=blk_y+ext_offset; y<max_blk_y+ext_offset; y++)
  {
    for(x=blk_x+ext_offset; x<max_blk_x+ext_offset; x++)
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
 *    copy from dec to rest
 *****************************************************************************************
 */
void CopyDecToRestBlock(int blk_addr, imgpel** ImgDec, imgpel** ImgRest, int blk_size, int width_in_alf_blks)
{
  int x, y, blk_x, blk_y, max_blk_x, max_blk_y;

  GetSamplePosForBALF(blk_addr, &blk_x, &blk_y, blk_size, width_in_alf_blks);

  max_blk_x = blk_x + blk_size;
  max_blk_y = blk_y + blk_size;
  if(max_blk_x > alf_width)  max_blk_x = alf_width;
  if(max_blk_y > alf_height) max_blk_y = alf_height;

  for(y=blk_y; y<max_blk_y; y++)
  {
    for(x=blk_x; x<max_blk_x; x++)
    {
      ImgRest[y][x] = ImgDec[y][x]; 
    }
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    set parameters for Block-based ALF
 *****************************************************************************************
 */
void SetBlockSizeParameter(alf_parameter_set_rbsp_t *alfps, int size_idx)
{
  alfps->alf_block_size_idx = size_idx;
  alfps->alf_block_size     = 1<<size_idx;
  alfps->WidthInAlfBlks     = alf_width  / alfps->alf_block_size;
  alfps->HeightInAlfBlks    = alf_height / alfps->alf_block_size;
  if( (alf_width  % alfps->alf_block_size)!=0 ) alfps->WidthInAlfBlks++;
  if( (alf_height % alfps->alf_block_size)!=0 ) alfps->HeightInAlfBlks++;
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
void CheckAvailabilityOfNeighborsAlfBlock()
{
  const int alfb_nr = img->current_alfb_nr;
  AlfBlock *currAlfBlk = &img->alf_blk_data[alfb_nr];

  // mark all neighbors as unavailable
  currAlfBlk->alfb_available_up   = NULL;
  currAlfBlk->alfb_available_left = NULL;

  currAlfBlk->alfbAddrA = alfb_nr - 1;
  currAlfBlk->alfbAddrB = alfb_nr - img->WidthInAlfBlks;
  currAlfBlk->alfbAddrC = alfb_nr - img->WidthInAlfBlks + 1;
  currAlfBlk->alfbAddrD = alfb_nr - img->WidthInAlfBlks - 1;

  currAlfBlk->alfbAvailA = ALFBlockIsAvailable(currAlfBlk->alfbAddrA) && (( alfb_nr    % img->WidthInAlfBlks)!=0);
  currAlfBlk->alfbAvailB = ALFBlockIsAvailable(currAlfBlk->alfbAddrB);
  currAlfBlk->alfbAvailC = ALFBlockIsAvailable(currAlfBlk->alfbAddrC) && (((alfb_nr+1) % img->WidthInAlfBlks)!=0);
  currAlfBlk->alfbAvailD = ALFBlockIsAvailable(currAlfBlk->alfbAddrD) && (( alfb_nr    % img->WidthInAlfBlks)!=0);
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
  AlfBlock *currAlfBlk = &img->alf_blk_data[curr_alfb_nr];
  int maxW, maxH;

  maxW = img->alf_block_size;
  maxH = img->alf_block_size;

  if ((xN<0)&&(yN<0))
  {
    pix->mb_addr   = currAlfBlk->alfbAddrD;
    pix->available = currAlfBlk->alfbAvailD;
  }
  else if ((xN<0)&&((yN>=0)&&(yN<maxH)))
  {
    pix->mb_addr   = currAlfBlk->alfbAddrA;
    pix->available = currAlfBlk->alfbAvailA;
  }
  else if (((xN>=0)&&(xN<maxW))&&(yN<0))
  {
    pix->mb_addr   = currAlfBlk->alfbAddrB;
    pix->available = currAlfBlk->alfbAvailB;
  }
  else if (((xN>=0)&&(xN<maxW))&&((yN>=0)&&(yN<maxH)))
  {
    pix->mb_addr   = curr_alfb_nr;
    pix->available = 1;
  }
  else if ((xN>=maxW)&&(yN<0))
  {
    pix->mb_addr   = currAlfBlk->alfbAddrC;
    pix->available = currAlfBlk->alfbAvailC;
  }
  else 
  {
    pix->available = 0;
  }

  if (pix->available)
  {
    pix->x = (xN + maxW) % maxW;
    pix->y = (yN + maxH) % maxH;

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
 *    initialize block partitioning data
 *****************************************************************************************
 */
void InitBlockPartition(int pblk_addr, alf_parameter_set_rbsp_t *alfps)
{
  int    pl, max_pl, bidx, pidx;
  int    pblk_size, pblk_x, pblk_y, pblk_width, num_pblock;
  int    sblk, sblk_size, sblk_x, sblk_y, sblk_width, num_sblock;
  int    sx, sy;

  pblk_size   = alfps->parent_block_size;
  pblk_width  = alfps->WidthInParentBlks;
  num_pblock  = alfps->num_parent_block;
  max_pl      = alfps->max_layer_level;

  GetSamplePosForBALF (pblk_addr, &pblk_x, &pblk_y, pblk_size, pblk_width);

  memset(qt_blk_data->valid_flag,   FALSE, sizeof(Boolean)*PYRAMIDAL_MEM_SIZE);
  memset(qt_blk_data->part_flag,    FALSE, sizeof(Boolean)*PYRAMIDAL_MEM_SIZE);
  memset(qt_blk_data->alf_blk_flag, FALSE, sizeof(Boolean)*PYRAMIDAL_MEM_SIZE);

  for (pl=0; pl<5; pl++)
  {
    sblk_width = 1<<pl;
    sblk_size  = pblk_size>>pl;
    num_sblock = 1<<(pl<<1);

    for (sblk=0; sblk<num_sblock; sblk++)
    {
      GetSamplePosForBALF (sblk, &sblk_x, &sblk_y, sblk_size, sblk_width);

      sx = pblk_x + sblk_x;
      sy = pblk_y + sblk_y;
      if( (sx<alf_width) && (sy<alf_height) )
      {
        bidx = block_scan_order[pl][sblk];
        pidx = pyramid_mem_idx[pl][bidx];
        qt_blk_data->valid_flag[pidx] = TRUE;
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
        pidx1 = pyramid_mem_idx[pl+1][(bidx<<2) ];
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
  alfps->WidthInParentBlks  = alf_width  / alfps->parent_block_size;
  alfps->HeightInParentBlks = alf_height / alfps->parent_block_size;
  if( (alf_width  % alfps->parent_block_size)!=0 ) alfps->WidthInParentBlks++;
  if( (alf_height % alfps->parent_block_size)!=0 ) alfps->HeightInParentBlks++;
  alfps->num_parent_block   = alfps->WidthInParentBlks * alfps->HeightInParentBlks;
}


//------- read ALF parameters
/*!
 *****************************************************************************************
 * \brief
 *    read ALF parameters
 *****************************************************************************************
 */
int read_ALFParameterSet(alf_parameter_set_rbsp_t *alfps, DataPartition *p)
{
  Bitstream *s = p->bitstream;
  
  assert (p != NULL);
  assert (p->bitstream != NULL);
  assert (p->bitstream->streamBuffer != 0);

  // ALF parameters
#ifdef ALF_SPATIAL_PREDICT_COEF
  alfps->pred_coef_mode = u_1("ALFPS: pred_coef_mode", s);
#endif

  // filter parameter for luma
  read_ALFLumaParam(alfps, s);

  // filter parameter for chroma
  read_ALFChromaParam(alfps, s);

  // block control parameter for luma
  read_ALFBlockControlParam(alfps, p);

  // padding
  if (s->frame_bitoffset%8 != 0 )
  {
    while(s->frame_bitoffset%8 != 0)
    {
      u_1("ALFPS: padding bits", s);
    }
  }

  return 1;
}

/*!
 *****************************************************************************************
 * \brief
 *    read ALF parameters for luma
 *****************************************************************************************
 */
int read_ALFLumaParam(alf_parameter_set_rbsp_t *alfps, Bitstream *s)
{
  int pos, N;
  char string[100];

  // ALF parameters for luma
  alfps->alf.tap = ue_v("ALFPS: alf_tap_size", s);
  alfps->alf.tap = (alfps->alf.tap<<1) + 5;

  N = (alfps->alf.tap*alfps->alf.tap+1)>>1;
  N = N + 1;  // DC offset
  alfps->alf.num_coeff = N;

  // filter coefficients for luma
  for(pos=0; pos<N; pos++)
  {
    sprintf(string, "ALFPS: %d pos", pos);
    alfps->alf.coeff[pos] = se_v(string, s);
  }

  return 1;
}

/*!
 *****************************************************************************************
 * \brief
 *    read ALF parameters for chroma
 *****************************************************************************************
 */
int read_ALFChromaParam(alf_parameter_set_rbsp_t *alfps, Bitstream *s)
{
  int pos, N;
  char string[100];

  // ALF parameters for chroma
  alfps->alf_chroma_idc = u_v(2, "ALFPS: alf_chroma_idc", s);

  if(alfps->alf_chroma_idc)
  {
    alfps->alf.tap_chroma = ue_v("ALFPS: filter_tap_chroma", s);
    alfps->alf.tap_chroma = (alfps->alf.tap_chroma<<1) + 5;

    N = (alfps->alf.tap_chroma*alfps->alf.tap_chroma+1)>>1;
    N = N + 1;  // DC offset
    alfps->alf.num_coeff_chroma = N;

    // filter coefficients for chroma
    for(pos=0; pos<N; pos++)
    {
      sprintf(string, "ALFPS: [chroma] %d pos", pos);
      alfps->alf.coeff_chroma[pos] = se_v(string, s);
    }
  }

  return 1;
}

/*!
 *****************************************************************************************
 * \brief
 *    read ALF parameters for block-based
 *****************************************************************************************
 */
int read_ALFBlockControlParam(alf_parameter_set_rbsp_t *alfps, DataPartition *p)
{
  Bitstream *s = p->bitstream;

  // block-based on-off parameter
  alfps->block_control_flag = u_1("ALFPS: block_control_flag", s);

  if(alfps->block_control_flag)
  {
    SetBlockSizeParameter(alfps, u_v(3, "ALFPS: alf_block_size_idx", s));

    alfps->qt_partition_flag  = u_1("ALFPS: qt_partition_flag", s);
    if (alfps->qt_partition_flag)
    {
      read_QuadTreeBlockPartitioningParam(alfps, s);
    }
    else
    {
      if (active_pps->entropy_coding_mode_flag == CABAC)
      {
        read_BlockFlagCABAC(alfps, p);

        if(alfps->block_control_flag)
        {
          s->frame_bitoffset = s->read_len<<3;
        }      
      }
      else
      {
        read_BlockFlagRunLength(alfps, s);
      }
    }
  }

  return 1;
}

/*!
 *****************************************************************************************
 * \brief
 *    read block flag by Run-length coding
 *****************************************************************************************
 */
int read_BlockFlagRunLength(alf_parameter_set_rbsp_t *alfps, Bitstream *s)
{
  int i;
  int alf_blk_flag_run=0;
  Boolean prev_alf_blk_flag=FALSE;
  AlfBlock *currAlfBlk;

  alf_blk_flag_run = 0;

  for (i=0; i<alfps->num_alf_block; i++)
  {
    currAlfBlk = &alfps->alf_blk_data[i];

    if( i==0 )
    {
      currAlfBlk->alf_blk_flag = (Boolean)u_1("ALFPS: alf_blk_flag", s);
      alf_blk_flag_run         = ue_v("ALFPS: alf_blk_flag_run", s);
      prev_alf_blk_flag        = currAlfBlk->alf_blk_flag;
    }
    else
    {
      if(alf_blk_flag_run==0)
      {
        currAlfBlk->alf_blk_flag = (prev_alf_blk_flag) ? FALSE : TRUE;
        if(i!=alfps->num_alf_block-1)
        {
          alf_blk_flag_run = ue_v("ALFPS: alf_blk_flag_run", s);
          prev_alf_blk_flag = currAlfBlk->alf_blk_flag;
        }
      }
      else
      {
        currAlfBlk->alf_blk_flag = prev_alf_blk_flag;
        alf_blk_flag_run--;
      }
    }
  }

  return 1;
}

/*!
 *****************************************************************************************
 * \brief
 *    read block flag by CABAC
 *****************************************************************************************
 */
int read_BlockFlagCABAC(alf_parameter_set_rbsp_t *alfps, DataPartition *p)
{
  int i;
  int len = 0;
  AlfBlock *currAlfBlk;
  SyntaxElement currSE;
  int ByteStartPosition;
#if TRACE
  char string[100];
#endif
  Bitstream *s = p->bitstream;

  ByteStartPosition = s->frame_bitoffset/8;
  if (s->frame_bitoffset%8 != 0)
  {
    ByteStartPosition++;
  }
  arideco_start_decoding_AlfBlkFlag (&p->de_cabac, s->streamBuffer, ByteStartPosition, &s->read_len);
  init_contexts (img);
  p->readSyntaxElement = readSyntaxElement_CABAC;

  for (i=0; i<alfps->num_alf_block; i++)
  {
    currAlfBlk = &alfps->alf_blk_data[i];

    img->current_alfb_nr = i;
    currAlfBlk = &alfps->alf_blk_data[i];
    currSE.type = SE_HEADER;

    CheckAvailabilityOfNeighborsAlfBlock();
    CheckAvailabilityOfNeighborsAlfBlockCABAC();
    currSE.reading = readAlfBlockFlagCABAC;
#if TRACE
    sprintf(string, "ALFPS: alf_blk_flag[%4d]", i);
    strncpy(currSE.tracestring,string,TRACESTRING_SIZE);
#endif
    p->readSyntaxElement(&currSE, img, input, p);
    len+=currSE.len;
    currAlfBlk->alf_blk_flag = (Boolean)currSE.value1;
  }

  return len;
}

/*!
 *****************************************************************************************
 * \brief
 *    read quad-tree block partitioning parameters
 *****************************************************************************************
 */
int read_QuadTreeBlockPartitioningParam(alf_parameter_set_rbsp_t *alfps, Bitstream *s)
{
  int len = 0;
  int pblk, layer_level, pidx;
  int layer, cblk, num_of_child_block, scblk;
  char string[100];

  alfps->max_layer_level = u_v(2, "ALFPS: max_layer_level", s) + 1;
  layer_level            = alfps->max_layer_level;

  SetParentBlockSizeParam(alfps, alfps->alf_block_size);

  for (pblk=0; pblk<alfps->num_parent_block; pblk++)
  {
    InitBlockPartition(pblk, alfps);

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
          qt_blk_data->part_flag[pidx] = u_1(string, s);
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
          qt_blk_data->alf_blk_flag[pidx] = u_1(string, s);
        }
      }
    }
    TransQbpToFbp(pblk, alfps);
  }

  return len;
}

#ifdef SWITCHED_FILTERS
/*!
 *****************************************************************************************
 * \brief
 *    read ALF parameters in SliceHeader
 *****************************************************************************************
 */
void read_ALFHeader()
{
  Slice *currSlice = img->currentSlice;
  int dP_nr = assignSE2partition[currSlice->dp_mode][SE_HEADER];
  DataPartition *partition = &(currSlice->partArr[dP_nr]);
  Bitstream *currStream = partition->bitstream;

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
}
#endif


#endif // ADAPTIVE_LOOP_FILTER
