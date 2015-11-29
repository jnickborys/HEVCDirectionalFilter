
/*!
*************************************************************************************
* \file image.c
*
* \brief
*    Code one image/slice
*
* \author
*    Main contributors (see contributors.h for copyright, address and affiliation details)
*     - Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
*     - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
*     - Jani Lainema                    <jani.lainema@nokia.com>
*     - Sebastian Purreiter             <sebastian.purreiter@mch.siemens.de>
*     - Byeong-Moon Jeon                <jeonbm@lge.com>
*     - Yoon-Seong Soh                  <yunsung@lge.com>
*     - Thomas Stockhammer              <stockhammer@ei.tum.de>
*     - Detlev Marpe                    <marpe@hhi.de>
*     - Guido Heising                   <heising@hhi.de>
*     - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
*     - Ragip Kurceren                  <ragip.kurceren@nokia.com>
*     - Antti Hallapuro                 <antti.hallapuro@nokia.com>
*     - Alexis Michael Tourapis         <alexismt@ieee.org> 
*************************************************************************************
*/
#include "contributors.h"

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <string.h>
#include <memory.h>
#include <assert.h>

#ifdef WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include "global.h"

#include "refbuf.h"
#include "mbuffer.h"
#include "intrarefresh.h"
#include "fmo.h"
#include "sei.h"
#include "memalloc.h"
#include "nalu.h"
#include "ratectl.h"
#include "mb_access.h"
#include "output.h"
#include "cabac.h"
#include "image.h"

#ifdef ADAPTIVE_QUANTIZATION
#include "adaptive_quantization.h"
#endif

#ifdef ADAPTIVE_FILTER
#include "adaptive_filter.h"
#endif

#ifdef ADAPTIVE_LOOP_FILTER
#include "adaptive_loop_filter.h"
#endif

#ifdef MV_COMPETITION
#include "mv_competition.h"
// <FTRD Compatibility with APEC
extern MV_Competition mv_comp;
extern int ******send_index_mv_prediction;
extern int *****predModeMV;
// FTRD>
#endif

#ifdef SWITCHED_FILTERS
#include "switched_filters.h"
#endif  // SWITCHED_FILTERS

extern pic_parameter_set_rbsp_t *PicParSet[MAXPPS];

void code_a_picture(Picture *pic);
void frame_picture (Picture *frame, int method);
void field_picture(Picture *top, Picture *bottom);

static int  writeout_picture(Picture *pic);

static int  picture_structure_decision(Picture *frame, Picture *top, Picture *bot);
static void distortion_fld (float *dis_fld_y, float *dis_fld_u, float *dis_fld_v);
static void find_snr(void);
#ifdef ADAPTIVE_QUANTIZATION
void find_distortion(void);
#else
static void find_distortion(void);
#endif

static void field_mode_buffer(int bit_field, float snr_field_y, float snr_field_u, float snr_field_v);
static void frame_mode_buffer (int bit_frame, float snr_frame_y, float snr_frame_u, float snr_frame_v);

static void init_frame(void);
static void init_field(void);

static void put_buffer_frame(void);
static void put_buffer_top(void);
static void put_buffer_bot(void);

static void copy_motion_vectors_MB(void);

#ifdef ADAPTIVE_FD_SD_CODING
static void store_motion_vectors_FDSD_coding   (void);
static void restore_motion_vectors_FDSD_coding (void);
#endif

static void PaddAutoCropBorders (int org_size_x, int org_size_y, int img_size_x, int img_size_y,
                                 int org_size_x_cr, int org_size_y_cr, int img_size_x_cr, int img_size_y_cr);

static void ReadOneFrame (int FrameNoInFile, int HeaderSize, int xs, int ys, int xs_cr, int ys_cr);

static void writeUnit(Bitstream* currStream ,int partition);
static void rdPictureCoding(void);

#ifdef _ADAPT_LAST_GROUP_
int *last_P_no;
int *last_P_no_frm;
int *last_P_no_fld;
#endif

static void ReportFirstframe(int tmp_time, int me_time);
static void ReportIntra(int tmp_time, int me_time);
static void ReportSP(int tmp_time, int me_time);
static void ReportP(int tmp_time, int me_time);
static void ReportB(int tmp_time, int me_time);
static void ReportNALNonVLCBits(int tmp_time, int me_time);

static int CalculateFrameNumber(void);  // Calculates the next frame number

StorablePicture *enc_picture;
StorablePicture *enc_frame_picture;
StorablePicture *enc_frame_picture2;
StorablePicture *enc_frame_picture3;
StorablePicture *enc_top_picture;
StorablePicture *enc_bottom_picture;
#ifdef ADAPTIVE_FILTER
StorablePicture *enc_frame_picture_aif;
#endif
//Rate control
int    QP;

const int ONE_FOURTH_TAP[3][2] =
{
  {20,20},
  {-5,-4},
  { 1, 0},
};

/*#ifdef USE_HP_FILTER
static double HIGH_PRECISION_H264_FILTER[15][36] = 
{
{ 0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
1.0/64.0  ,   -5.0/64.0  ,   52.0/64.0  ,   20.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0
}, // a_pos
{ 0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
1.0/32.0  ,   -5.0/32.0  ,   20.0/32.0  ,   20.0/32.0  ,   -5.0/32.0  ,   1.0/32.0  ,
0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0
}, // b_pos
{ 0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
1.0/64.0  ,   -5.0/64.0  ,   20.0/64.0  ,   52.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0
}, // c_pos
{ 0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,   52.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0
}, // d_pos
{ 0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
1.0/64.0  ,   -5.0/64.0  ,   40.0/64.0  ,   20.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0
}, // e_pos
{ 1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0,
-5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
52.0/2048.0, -260.0/2048.0, 1040.0/2048.0, 1040.0/2048.0, -260.0/2048.0,  52.0/2048.0,
20.0/2048.0, -100.0/2048.0,  400.0/2048.0,  400.0/2048.0, -100.0/2048.0,  20.0/2048.0,
-5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0
}, // f_pos
{ 0.0       ,    0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,   0.0       ,
1.0/64.0  ,   -5.0/64.0  ,   20.0/64.0  ,   40.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
0.0       ,    0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,   0.0
}, // g_pos
{ 0.0       ,    0.0       ,    1.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,   -5.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,   20.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,   20.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,   -5.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    1.0/32.0  ,    0.0       ,    0.0       ,   0.0
}, // h_pos
{ 1.0/2048.0,   -5.0/2048.0,   52.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0,
-5.0/2048.0,   25.0/2048.0, -260.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
20.0/2048.0, -100.0/2048.0, 1040.0/2048.0,  400.0/2048.0, -100.0/2048.0,  20.0/2048.0,
20.0/2048.0, -100.0/2048.0, 1040.0/2048.0,  400.0/2048.0, -100.0/2048.0,  20.0/2048.0,
-5.0/2048.0,   25.0/2048.0, -260.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
1.0/2048.0,   -5.0/2048.0,   52.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0
}, // i_pos
{ 1.0/1024.0,   -5.0/1024.0,   20.0/1024.0,   20.0/1024.0,   -5.0/1024.0,   1.0/1024.0,
-5.0/1024.0,   25.0/1024.0, -100.0/1024.0, -100.0/1024.0,   25.0/1024.0,  -5.0/1024.0,
20.0/1024.0, -100.0/1024.0,  400.0/1024.0,  400.0/1024.0, -100.0/1024.0,  20.0/1024.0,
20.0/1024.0, -100.0/1024.0,  400.0/1024.0,  400.0/1024.0, -100.0/1024.0,  20.0/1024.0,
-5.0/1024.0,   25.0/1024.0, -100.0/1024.0, -100.0/1024.0,   25.0/1024.0,  -5.0/1024.0,
1.0/1024.0,   -5.0/1024.0,   20.0/1024.0,   20.0/1024.0,   -5.0/1024.0,   1.0/1024.0
}, // j_pos
{ 1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   52.0/2048.0,   -5.0/2048.0,   1.0/2048.0,
-5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -260.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
20.0/2048.0, -100.0/2048.0,  400.0/2048.0, 1040.0/2048.0, -100.0/2048.0,  20.0/2048.0,
20.0/2048.0, -100.0/2048.0,  400.0/2048.0, 1040.0/2048.0, -100.0/2048.0,  20.0/2048.0,
-5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -260.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   52.0/2048.0,   -5.0/2048.0,   1.0/2048.0
}, // k_pos
{ 0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,   52.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0
}, // l_pos
{ 0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
1.0/64.0  ,   -5.0/64.0  ,   40.0/64.0  ,   20.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0
}, // m_pos
{ 1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0,
-5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
20.0/2048.0, -100.0/2048.0,  400.0/2048.0,  400.0/2048.0, -100.0/2048.0,  20.0/2048.0,
52.0/2048.0, -260.0/2048.0, 1040.0/2048.0, 1040.0/2048.0, -260.0/2048.0,  52.0/2048.0,
-5.0/2048.0,   25.0/2048.0, -100.0/2048.0, -100.0/2048.0,   25.0/2048.0,  -5.0/2048.0,
1.0/2048.0,   -5.0/2048.0,   20.0/2048.0,   20.0/2048.0,   -5.0/2048.0,   1.0/2048.0
}, // n_pos
{ 0.0       ,    0.0       ,    0.0     ,      1.0/64.0  ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    0.0     ,     -5.0/64.0  ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    0.0     ,     20.0/64.0  ,    0.0       ,   0.0       ,
1.0/64.0  ,   -5.0/64.0  ,   20.0/64.0,     40.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
0.0       ,    0.0       ,    0.0     ,     -5.0/64.0  ,    0.0       ,   0.0       ,
0.0       ,    0.0       ,    0.0     ,      1.0/64.0  ,    0.0       ,   0.0
} // o_pos
};
#endif*/

void MbAffPostProc(void)
{
  imgpel temp[32][16];

  imgpel ** imgY  = enc_picture->imgY;
  imgpel ***imgUV = enc_picture->imgUV;
  int i, y, x0, y0, uv;

  if (img->yuv_format != YUV400)
  {
    for (i=0; i<(int)img->PicSizeInMbs; i+=2)
    {
      if (enc_picture->mb_field[i])
      {
        get_mb_pos(i, &x0, &y0);
        for (y=0; y<(2*MB_BLOCK_SIZE);y++)
          memcpy(&temp[y],&imgY[y0+y][x0], MB_BLOCK_SIZE * sizeof(imgpel));

        for (y=0; y<MB_BLOCK_SIZE;y++)
        {
          memcpy(&imgY[y0+(2*y)][x0],temp[y], MB_BLOCK_SIZE * sizeof(imgpel));
          memcpy(&imgY[y0+(2*y + 1)][x0],temp[y+ MB_BLOCK_SIZE], MB_BLOCK_SIZE * sizeof(imgpel));
        }

        x0 = x0 / (16/img->mb_cr_size_x);
        y0 = y0 / (16/img->mb_cr_size_y);

        for (uv=0; uv<2; uv++)
        {
          for (y=0; y<(2*img->mb_cr_size_y);y++)    
            memcpy(&temp[y],&imgUV[uv][y0+y][x0], img->mb_cr_size_x * sizeof(imgpel));

          for (y=0; y<img->mb_cr_size_y;y++)
          {
            memcpy(&imgUV[uv][y0+(2*y)][x0],temp[y], img->mb_cr_size_x * sizeof(imgpel));
            memcpy(&imgUV[uv][y0+(2*y + 1)][x0],temp[y+ img->mb_cr_size_y], img->mb_cr_size_x * sizeof(imgpel));

          }
        }
      }
    }
  }
  else
  {
    for (i=0; i<(int)img->PicSizeInMbs; i+=2)
    {
      if (enc_picture->mb_field[i])
      {
        get_mb_pos(i, &x0, &y0);
        for (y=0; y<(2*MB_BLOCK_SIZE);y++)
          memcpy(&temp[y],&imgY[y0+y][x0], MB_BLOCK_SIZE * sizeof(imgpel));

        for (y=0; y<MB_BLOCK_SIZE;y++)
        {
          memcpy(&imgY[y0+(2*y)][x0],temp[y], MB_BLOCK_SIZE * sizeof(imgpel));
          memcpy(&imgY[y0+(2*y + 1)][x0],temp[y+ MB_BLOCK_SIZE], MB_BLOCK_SIZE * sizeof(imgpel));
        }
      }
    }
  }
}

/*!
************************************************************************
* \brief
*    Encodes a picture
*
*    This is the main picture coding loop.. It is called by all this
*    frame and field coding stuff after the img-> elements have been
*    set up.  Not sure whether it is useful for MB-adaptive frame/field
*    coding
************************************************************************
*/
void code_a_picture(Picture *pic)
{
  unsigned int NumberOfCodedMBs = 0;
  int SliceGroup = 0;
  int j;

  img->currentPicture = pic;

  img->currentPicture->idr_flag = ((!IMG_NUMBER) && (!(img->structure==BOTTOM_FIELD))) || (input->idr_enable && (img->type == I_SLICE || img->type==SI_SLICE)&& (!(img->structure==BOTTOM_FIELD)));

  pic->no_slices = 0;
  pic->distortion_u = pic->distortion_v = pic->distortion_y = 0.0;

  RandomIntraNewPicture ();     //! Allocates forced INTRA MBs (even for fields!)

  // The slice_group_change_cycle can be changed here.
  // FmoInit() is called before coding each picture, frame or field
  img->slice_group_change_cycle=1;
  FmoInit(img, active_pps, active_sps);
  FmoStartPicture ();           //! picture level initialization of FMO

  CalculateQuantParam();
  CalculateOffsetParam();

  if(input->Transform8x8Mode)
  {
    CalculateQuant8Param();
    CalculateOffset8Param();
  }

  reset_pic_bin_count();
  img->bytes_in_picture = 0;

  while (NumberOfCodedMBs < img->PicSizeInMbs)       // loop over slices
  {
    // Encode one SLice Group
    while (!FmoSliceGroupCompletelyCoded (SliceGroup))
    {
      // Encode the current slice
      NumberOfCodedMBs += encode_one_slice (SliceGroup, pic, NumberOfCodedMBs);
      FmoSetLastMacroblockInSlice (img->current_mb_nr);
      // Proceed to next slice
      img->current_slice_nr++;
      stats->bit_slice = 0;
    }
    // Proceed to next SliceGroup
    SliceGroup++;
  }
  FmoEndPicture ();

  // Modified for Fast Mode Decision. Inchoon Choi, SungKyunKwan Univ.
  if (input->rdopt == 3 && (img->type != B_SLICE))
    for (j = 0; j < input->NoOfDecoders; j++)
      DeblockFrame (img, decs->decY_best[j], NULL);

  DeblockFrame (img, enc_picture->imgY, enc_picture->imgUV); //comment out to disable loop filter

  if (img->MbaffFrameFlag)
    MbAffPostProc();

#ifdef MV_COMPETITION // bug fix by AT
#ifdef ADAPTIVE_QUANTIZATION
  if(!input->RDPictureDecision && !input->UseAdaptiveQuantMatrix)
#endif
    // <FTRD Compatibility with AIF
#ifdef ADAPTIVE_FILTER
    if(input->UseAdaptiveFilter==0)
#endif
      // FTRD>
      if((img->type == B_SLICE) && (input->successive_Bframe > 1) && (input->mv_competition > 0))
        Copy_MV_B_frame(enc_picture->mv, enc_picture->ref_idx);
#endif

}



/*!
************************************************************************
* \brief
*    Encodes one frame
************************************************************************
*/
int encode_one_frame (void)
{
  static int prev_frame_no = 0; // POC200301
  static int consecutive_non_reference_pictures = 0; // POC200301
  int        FrameNumberInFile;
  int        i, j;

#ifdef _LEAKYBUCKET_
  //extern long Bit_Buffer[20000];
  extern unsigned long total_frame_buffer;
#endif
#ifdef ADAPTIVE_QUANTIZATION
  int rd_qp=0, AQMS_Skip_flag=0, mincost_flg=0;
  sDefaultInputParam_t defInputParam;
#endif

  time_t ltime1;
  time_t ltime2;

#ifdef WIN32
  struct _timeb tstruct1;
  struct _timeb tstruct2;
#else
  struct timeb tstruct1;
  struct timeb tstruct2;
#endif

  int tmp_time;
  int bits_frm = 0, bits_fld = 0;
  float dis_frm = 0.0, dis_frm_y = 0.0, dis_frm_u = 0.0, dis_frm_v = 0.0;
  float dis_fld = 0.0, dis_fld_y = 0.0, dis_fld_u = 0.0, dis_fld_v = 0.0;

  //Rate control
  int pic_type, bits = 0; 

  me_time=0;
  img->rd_pass = 0;
  enc_frame_picture  = NULL;
  enc_frame_picture2 = NULL;
  enc_frame_picture3 = NULL;
#ifdef ADAPTIVE_FILTER
  enc_frame_picture_aif = NULL;
#endif
#ifdef ADAPTIVE_QUANTIZATION
  enc_frame_picture_aqms = NULL;
#endif
#ifdef WIN32
  _ftime (&tstruct1);           // start time ms
#else
  ftime (&tstruct1);
#endif
  time (&ltime1);               // start time s

#ifdef USE_NEW_OFFSET
  if(input->UseNewOffset)
    img->frameOffsetAvail = 0; 
#endif 

  //Rate control 
  img->write_macroblock = 0;
  /*
  //Shankar Regunathan (Oct 2002)
  //Prepare Panscanrect SEI payload
  UpdatePanScanRectInfo ();
  //Prepare Arbitrarydata SEI Payload
  UpdateUser_data_unregistered ();
  //Prepare Registered data SEI Payload
  UpdateUser_data_registered_itu_t_t35 ();
  //Prepare RandomAccess SEI Payload
  UpdateRandomAccess ();
  */

  if (input->ResendPPS && img->number !=0)
  {
    stats->bit_ctr_parametersets_n=write_PPS(0, 0);
    //stats->bit_slice += stats->bit_ctr_parametersets_n;
    stats->bit_ctr_parametersets += stats->bit_ctr_parametersets_n;
  }

  put_buffer_frame ();      // sets the pointers to the frame structures 
  // (and not to one of the field structures)
  init_frame ();
  FrameNumberInFile = CalculateFrameNumber();

  ReadOneFrame (FrameNumberInFile, input->infile_header,
    input->img_width, input->img_height, input->img_width_cr, input->img_height_cr);

  PaddAutoCropBorders (input->img_width, input->img_height, img->width, img->height,
    input->img_width_cr, input->img_height_cr, img->width_cr, img->height_cr);

  // set parameters for direct mode and deblocking filter
  img->direct_spatial_mv_pred_flag     = input->direct_spatial_mv_pred_flag;
  img->LFDisableIdc    = input->LFDisableIdc;
  img->LFAlphaC0Offset = input->LFAlphaC0Offset;
  img->LFBetaOffset    = input->LFBetaOffset;
  img->AdaptiveRounding = input->AdaptiveRounding;

#ifdef ADAPTIVE_QUANTIZATION
  if(input->UseAdaptiveQuantMatrix) // skip consideration
  {
    img->slice_fractional_quant_flag = !(input->UseAdaptiveQuantMatrix==0);
    if(input->ProfileIDC < FREXT_HP)
      img->slice_fractional_quant_flag = (SkipRatio[(img->type)+(img->type==B_SLICE? (img->nal_reference_idc<<2) : 0)] < 90.0);
    AQMS_Skip_flag = (!img->slice_fractional_quant_flag);
  }
#endif
  // Following code should consider optimal coding mode. Currently also does not support
  // multiple slices per frame.
  frame_ctr[img->type]++;
  snr->frame_ctr++;

  if(img->type == SP_SLICE)
  { 
    if(input->sp2_frame_indicator)
    { // switching SP frame encoding
      sp2_frame_indicator=1;
      read_SP_coefficients();
    }
  }
  else
  {
    sp2_frame_indicator=0;
  }

  if (input->PicInterlace == FIELD_CODING)
  {
    //Rate control
    img->FieldControl=1;

    img->field_picture = 1;  // we encode fields
    field_picture (top_pic, bottom_pic);
    img->fld_flag = 1;
  }
  else
  {
    //Rate control
    img->FieldControl=0;

    // For frame coding, turn MB level field/frame coding flag on
    if (input->MbInterlace)
      mb_adaptive = 1;

    img->field_picture = 0; // we encode a frame

    //Rate control
    if(input->RCEnable)
    { 
      /*update the number of MBs in the basic unit for MB adaptive 
      f/f coding*/
      if((input->MbInterlace)&&(input->basicunit<img->Frame_Total_Number_MB)\
        &&(img->type==P_SLICE)&&(img->IFLAG==0))
        img->BasicUnit=input->basicunit*2;
      else
        img->BasicUnit=input->basicunit;

      rc_init_pict(1,0,1); 
      img->qp  = updateQuantizationParameter(0); 


      pic_type = img->type;
      QP =0;
    }

    if( active_sps->frame_mbs_only_flag)
      img->TopFieldFlag=0;

#ifdef ADAPTIVE_QUANTIZATION
    gdwBaseQp = img->qp;
    if (input->UseAdaptiveQuantMatrix )
    {
      SaveInputParam(AQMS_Skip_flag,&defInputParam);
      rd_qp = img->qp;
    }
#endif


    if (input->GenerateMultiplePPS)
      active_pps = PicParSet[0];
#ifdef ADAPTIVE_FILTER
    if(UseAdaptiveFilterForCurrentFrame())
      img->AdaptiveFilterFlag = 0;
    else
      img->AdaptiveFilterFlag = -1;
#endif  
#ifdef DIRECTIONAL_FILTER
    if(input->UseAdaptiveFilter==FILTER_TYPE_EDAIF)
      initEDAIF2();
    else
      initFilterCustom(input->UseAdaptiveFilter);
#endif
    frame_picture (frame_pic_1, 0);

#ifdef ADAPTIVE_FILTER
    if(img->AdaptiveFilterFlag==0)
    {
      if(input->InterpolationDecision)
      {
        if (input->UseAdaptiveFilter == 2) // separable aif
        {
          UnifiedOneForthPixWithNewFilterHor();
          UnifiedOneForthPixWithNewFilterVer();        
        }
#ifdef DIRECTIONAL_FILTER
        if (input->UseAdaptiveFilter == 3)
        {
          if (input->ImpType == IMP_FLOAT32)
            UnifiedOneForthPixWith_1DAIF_float();
          else if (input->ImpType == IMP_INT16)
            UnifiedOneForthPixWith_1DAIF_int16();        
          else 
            printf("DAIF: wrong implementation type!\n");
        }
#endif  // DIRECTIONAL_FILTER
#ifdef E_DAIF
        else if (input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)
        {
          UnifiedOneForthPixWith_EDAIF();
        }
#endif  // E_DAIF
#ifdef EAIF
        else if (input->UseAdaptiveFilter == FILTER_TYPE_EAIF)
        {
          UnifiedOneForthPixWith_EAIF();
        }
#endif
        else 
          UnifiedOneForthPixWithNewFilter();

        if (getDecisionOnAIF_MCP() == 0)
        {
          enc_picture=enc_frame_picture;
          frame_pic = frame_pic_1;
          img->AdaptiveFilterFlag = -1;
          //printf("MCP_De: AVC\n");
        }
        else 
        {
          //printf("MCP_De: AIF\n");
        }
      }
      if(img->AdaptiveFilterFlag==0)
      {
        img->AdaptiveFilterFlag=1;
        frame_picture (frame_pic_aif,3);
        img->AdaptiveFilterFlag=picture_coding_decision(frame_pic_1, frame_pic_aif, img->qp);
        if(img->AdaptiveFilterFlag == 1)
        {
          //printf("RDO_De: AIF\n");
          enc_picture=enc_frame_picture_aif;
          frame_pic = frame_pic_aif;
        }
        else
        {
          //printf("RDO_De: AVC\n");
          enc_picture=enc_frame_picture;
          frame_pic = frame_pic_1;
          img->AdaptiveFilterFlag = -1;
        }
      }      
    }
#endif

#ifdef SWITCHED_FILTERS
    if ((input->RDPictureIntra || img->type != I_SLICE) && input->RDPictureDecision)
    {
      if(input->UseHPFilter == HPF_SIFO)
      {
        if(img->type == P_SLICE)
          rdPictureCodingFilterP();
        else
          rdPictureCodingFilterB();
      }
      else if(input->UseHPFilter == HPF_SIFO_FPO)
      {
        if(img->type == P_SLICE)
          rdPictureCodingFilterP_IPass();

        ComputeFiltersAndOffsets();
        if(img->type == P_SLICE)
          UpdateSequenceFilters_P();
        else
          UpdateSequenceFilters_B();
      }
      else
      {
        rdPictureCoding();
      }
    }
#else
    if ((input->RDPictureIntra || img->type!=I_SLICE) && input->RDPictureDecision)
    {
      rdPictureCoding();
    }
#endif  // SWITCHED_FILTERS

#ifdef ADAPTIVE_QUANTIZATION
    rd_qp = gdwBaseQp = img->qp;
    if(input->UseAdaptiveQuantMatrix)
    {
      rdPictureCodingForIAQMS(AQMS_Skip_flag, &defInputParam, &mincost_flg, rd_qp);
    }
#endif

    // <FTRD Compatibility with AIF
#ifdef MV_COMPETITION // B-slice mv memory override bug fix by AT 2006/11/07
    if( input->RDPictureDecision
#ifdef ADAPTIVE_QUANTIZATION
      || input->UseAdaptiveQuantMatrix 
#endif //ADAPTIV_QUANTIZATION  
#ifdef ADAPTIVE_FILTER
      || (input->UseAdaptiveFilter > 0)
#endif //ADAPTIVE_FILTER 
      )
    {
      if((img->type == B_SLICE) && (input->successive_Bframe > 1))
        Copy_MV_B_frame(enc_picture->mv, enc_picture->ref_idx);
    }
#endif  // MV_COMPETITION
    // FTRD>



    if(img->type==SP_SLICE && si_frame_indicator==0 && input->si_frame_indicator)
    {
      // once the picture has been encoded as a primary SP frame encode as an SI frame
      si_frame_indicator=1;
      frame_picture (frame_pic_2, 0);
    }
    if(img->type==SP_SLICE && input->sp_output_indicator)
    {
      // output the transformed and quantized coefficients (useful for switching SP frames)
      output_SP_coefficients();
    }
    // For field coding, turn MB level field/frame coding flag off
    if (input->MbInterlace)
      mb_adaptive = 0;

    if (input->PicInterlace == ADAPTIVE_CODING)
    {
      //Rate control
      img->FieldControl=1;
      img->write_macroblock = 0;
      img->bot_MB = 0;

      img->field_picture = 1;  // we encode fields
      field_picture (top_pic, bottom_pic);

      //! Note: the distortion for a field coded picture is stored in the top field
      //! the distortion values in the bottom field are dummies
      dis_fld = top_pic->distortion_y + top_pic->distortion_u + top_pic->distortion_v;
      dis_frm = frame_pic_1->distortion_y + frame_pic_1->distortion_u + frame_pic_1->distortion_v;

      if(img->rd_pass==0)
        img->fld_flag = picture_structure_decision (frame_pic_1, top_pic, bottom_pic);
      else if(img->rd_pass==1)
        img->fld_flag = picture_structure_decision (frame_pic_2, top_pic, bottom_pic);
      else
        img->fld_flag = picture_structure_decision (frame_pic_3, top_pic, bottom_pic);

      update_field_frame_contexts (img->fld_flag);

      //Rate control
      if(img->fld_flag==0)
        img->FieldFrame=1;
      /*the current choice is field coding*/
      else
        img->FieldFrame=0;
    }
    else

      img->fld_flag = 0;
  }

  if (img->fld_flag)
    stats->bit_ctr_emulationprevention += stats->em_prev_bits_fld;
  else
    stats->bit_ctr_emulationprevention += stats->em_prev_bits_frm;

  if (img->type != B_SLICE)
  {
    img->pstruct_next_P = img->fld_flag;
  }

  // Here, img->structure may be either FRAME or BOTTOM FIELD depending on whether AFF coding is used
  // The picture structure decision changes really only the fld_flag

  if (img->fld_flag)            // field mode (use field when fld_flag=1 only)
  {
    field_mode_buffer (bits_fld, dis_fld_y, dis_fld_u, dis_fld_v);
    writeout_picture (top_pic);
    writeout_picture (bottom_pic);
  }
  else                          //frame mode
  {
    frame_mode_buffer (bits_frm, dis_frm_y, dis_frm_u, dis_frm_v);

#ifdef ADAPTIVE_QUANTIZATION
    if (input->UseAdaptiveQuantMatrix && mincost_flg == 0)
    {
      writeout_picture (frame_pic_aqms);
    }
    else
    {
      if (input->RDPictureDecision && img->rd_pass == 2)
        writeout_picture (frame_pic_3);
      else if (input->RDPictureDecision && img->rd_pass == 1)
        writeout_picture (frame_pic_2);
      else
        if(img->type==SP_SLICE && si_frame_indicator==1)
        {
          writeout_picture (frame_pic_2);
          si_frame_indicator=0;
        }
        else

#ifdef ADAPTIVE_FILTER
          if(img->AdaptiveFilterFlag == 1)
            writeout_picture(frame_pic_aif);
          else
#endif
            writeout_picture (frame_pic_1);
    }
#else
    if (input->RDPictureDecision && img->rd_pass == 2)
      writeout_picture (frame_pic_3);
    else if (input->RDPictureDecision && img->rd_pass == 1)
      writeout_picture (frame_pic_2);
    else
      if(img->type==SP_SLICE && si_frame_indicator==1)
      {
        writeout_picture (frame_pic_2);
        si_frame_indicator=0;
      }
      else
#ifdef ADAPTIVE_FILTER
        if(img->AdaptiveFilterFlag == 1)
          writeout_picture(frame_pic_aif);
        else
#endif
          writeout_picture (frame_pic_1);
#endif
  }

  if (frame_pic_3)
    free_slice_list(frame_pic_3);  
  if (frame_pic_2)
    free_slice_list(frame_pic_2);  
  if (frame_pic_1)
    free_slice_list(frame_pic_1);
  if (top_pic)
    free_slice_list(top_pic);
  if (bottom_pic)
    free_slice_list(bottom_pic);
#ifdef ADAPTIVE_FILTER
  if(frame_pic_aif)
    free_slice_list(frame_pic_aif);
  if(UseAdaptiveFilterForCurrentFrame() && img->AdaptiveFilterFlag==-1)
    img->AdaptiveFilterFlag=0;      // for the output only
#endif
#ifdef ADAPTIVE_QUANTIZATION
  if(input->UseAdaptiveQuantMatrix)
  {
    if (frame_pic_aqms);
    free_slice_list(frame_pic_aqms);
  }
#endif
  /*
  // Tian Dong (Sept 2002)
  // in frame mode, the newly reconstructed frame has been inserted to the mem buffer
  // and it is time to prepare the spare picture SEI payload.
  if (input->InterlaceCodingOption == FRAME_CODING
  && input->SparePictureOption && img->type != B_SLICE)
  CalculateSparePicture ();
  */

  //Rate control
  if(input->RCEnable)
  {
    bits = (int) (stats->bit_ctr - stats->bit_ctr_n);
    rc_update_pict_frame(bits);
  }

  if (input->PicInterlace == FRAME_CODING)
  {
    if (input->rdopt == 3 && img->type != B_SLICE)
      UpdateDecoders ();      // simulate packet losses and move decoded image to reference buffers

    if (input->RestrictRef)
      UpdatePixelMap ();
  }

  if (input->Verbose != 0)
    find_snr ();
  else
  {
    snr->snr_y = 0.0;
    snr->snr_u = 0.0;
    snr->snr_v = 0.0;
    snr->sse_y = 0.0;
    snr->sse_u = 0.0;
    snr->sse_v = 0.0;
  }

  // redundant pictures: save reconstruction to calculate SNR and replace reference picture
  if(input->redundant_pic_flag && key_frame)
  {
    for(i=0; i<img->width; i++)
    {
      for(j=0; j<img->height; j++)
      {
        imgY_tmp[j][i] = enc_frame_picture->imgY[j][i];
      }
    }
    for(i=0; i<img->width_cr; i++)
    {
      for(j=0; j<img->height_cr; j++)
      {
        imgUV_tmp[0][j][i] = enc_frame_picture->imgUV[0][j][i];
        imgUV_tmp[1][j][i] = enc_frame_picture->imgUV[1][j][i];
      }
    }
  }

  if(input->redundant_pic_flag && redundant_coding)
  {
    for(i=0; i<img->width; i++)
    {
      for(j=0; j<img->height; j++)
      {
        enc_frame_picture->imgY[j][i] = imgY_tmp[j][i];
      }
    }
    for(i=0; i<img->width_cr; i++)
    {
      for(j=0; j<img->height_cr; j++)
      {
        enc_frame_picture->imgUV[0][j][i] = imgUV_tmp[0][j][i];
        enc_frame_picture->imgUV[1][j][i] = imgUV_tmp[1][j][i];
      }
    }
  }

  time (&ltime2);               // end time sec
#ifdef WIN32
  _ftime (&tstruct2);           // end time ms
#else
  ftime (&tstruct2);            // end time ms
#endif

  tmp_time = (int)((ltime2 * 1000 + tstruct2.millitm) - (ltime1 * 1000 + tstruct1.millitm));
  tot_time = tot_time + tmp_time;

  if (input->PicInterlace == ADAPTIVE_CODING)
  {
    if (img->fld_flag)
    {
      // store bottom field
      store_picture_in_dpb(enc_bottom_picture);
      free_storable_picture(enc_frame_picture);
    }
    else
    {
      // replace top with frame
      replace_top_pic_with_frame(enc_frame_picture);
      free_storable_picture(enc_bottom_picture);
    }
  }
  else
  {
    if (img->fld_flag)
    {
      store_picture_in_dpb(enc_bottom_picture);
    }
    else
    {
#ifdef ADAPTIVE_QUANTIZATION
      if ( input->UseAdaptiveQuantMatrix && mincost_flg == 0)
      {
        store_picture_in_dpb(enc_frame_picture_aqms);
        free_storable_picture(enc_frame_picture);
        free_storable_picture(enc_frame_picture2);
        free_storable_picture(enc_frame_picture3);
#ifdef ADAPTIVE_FILTER
        free_storable_picture(enc_frame_picture_aif);
#endif
      }
      else
      {
        if (img->rd_pass==2)
        {
          store_picture_in_dpb(enc_frame_picture3);
          free_storable_picture(enc_frame_picture);
          free_storable_picture(enc_frame_picture2);
#ifdef ADAPTIVE_FILTER
          free_storable_picture(enc_frame_picture_aif);
#endif
        }
        else if (img->rd_pass==1)
        {
          store_picture_in_dpb(enc_frame_picture2);
          free_storable_picture(enc_frame_picture);
          free_storable_picture(enc_frame_picture3);
#ifdef ADAPTIVE_FILTER
          free_storable_picture(enc_frame_picture_aif);
#endif
        }
        else
        {

          if(input->redundant_pic_flag==0)
          {
#ifdef ADAPTIVE_FILTER
            if(img->AdaptiveFilterFlag==1)
            {
              store_picture_in_dpb(enc_frame_picture_aif);
              free_storable_picture(enc_frame_picture);
            }
            else
#endif
            {
              store_picture_in_dpb(enc_frame_picture);
#ifdef ADAPTIVE_FILTER
              free_storable_picture(enc_frame_picture_aif);
#endif
            }
            free_storable_picture(enc_frame_picture2);
            free_storable_picture(enc_frame_picture3);
          }
          else     
          {
            // key picture will be stored in dpb after redundant picture is coded
            if(key_frame==0)
            {
              store_picture_in_dpb(enc_frame_picture);
              free_storable_picture(enc_frame_picture2);
              free_storable_picture(enc_frame_picture3);
            }
          }
        }
        free_storable_picture(enc_frame_picture_aqms);
      }
#else
      if (img->rd_pass==2)
      {
        store_picture_in_dpb(enc_frame_picture3);
        free_storable_picture(enc_frame_picture);
        free_storable_picture(enc_frame_picture2);
#ifdef ADAPTIVE_FILTER
        free_storable_picture(enc_frame_picture_aif);
#endif
      }
      else if (img->rd_pass==1)
      {
        store_picture_in_dpb(enc_frame_picture2);
        free_storable_picture(enc_frame_picture);
        free_storable_picture(enc_frame_picture3);
#ifdef ADAPTIVE_FILTER
        free_storable_picture(enc_frame_picture_aif);
#endif
      }
      else
      {
        if(input->redundant_pic_flag==0)
        {
#ifdef ADAPTIVE_FILTER
          if(img->AdaptiveFilterFlag==1)
          {
            store_picture_in_dpb(enc_frame_picture_aif);
            free_storable_picture(enc_frame_picture);
          }
          else
#endif
          {
            store_picture_in_dpb(enc_frame_picture);
#ifdef ADAPTIVE_FILTER
            free_storable_picture(enc_frame_picture_aif);
#endif
          }
          free_storable_picture(enc_frame_picture2);
          free_storable_picture(enc_frame_picture3);
        }
        else     
        {
          // key picture will be stored in dpb after redundant picture is coded
          if(key_frame==0)
          {
            store_picture_in_dpb(enc_frame_picture);
            free_storable_picture(enc_frame_picture2);
            free_storable_picture(enc_frame_picture3);
          }
        }

      }
#endif
    }
  }

#ifdef _LEAKYBUCKET_
  // Store bits used for this frame and increment counter of no. of coded frames
  Bit_Buffer[total_frame_buffer] = (int) (stats->bit_ctr - stats->bit_ctr_n);
  total_frame_buffer++;
#endif

  // POC200301: Verify that POC coding type 2 is not used if more than one consecutive 
  // non-reference frame is requested or if decoding order is different from output order
  if (img->pic_order_cnt_type == 2)
  {
    if (!img->nal_reference_idc) consecutive_non_reference_pictures++;
    else consecutive_non_reference_pictures = 0;

    if (frame_no < prev_frame_no || consecutive_non_reference_pictures>1)
      error("POC type 2 cannot be applied for the coding pattern where the encoding /decoding order of pictures are different from the output order.\n", -1);
    prev_frame_no = frame_no;
  }

  if (stats->bit_ctr_parametersets_n!=0)
    ReportNALNonVLCBits(tmp_time, me_time);

  if (IMG_NUMBER == 0)
    ReportFirstframe(tmp_time,me_time);
  else
  {
    //Rate control
    if(input->RCEnable)
    {
      if((!input->PicInterlace)&&(!input->MbInterlace))
        bits=(int) (stats->bit_ctr - stats->bit_ctr_n);
      else
      {
        bits = (int) (stats->bit_ctr - Pprev_bits); // used for rate control update */
        Pprev_bits = stats->bit_ctr;
      }
    }

    switch (img->type)
    {
    case I_SLICE:
      stats->bit_ctr_I += stats->bit_ctr - stats->bit_ctr_n;
      ReportIntra(tmp_time,me_time);
      break;
    case SP_SLICE:
      stats->bit_ctr_P += stats->bit_ctr - stats->bit_ctr_n;
      ReportSP(tmp_time,me_time);
      break;
    case B_SLICE:
      stats->bit_ctr_B += stats->bit_ctr - stats->bit_ctr_n;
      ReportB(tmp_time,me_time);
      break;
    default:      // P
      stats->bit_ctr_P += stats->bit_ctr - stats->bit_ctr_n;
      ReportP(tmp_time,me_time);      
    }
  }  


  if (input->Verbose == 0)
  { 
    //for (i = 0; i <= (img->number & 0x0F); i++)
    //printf(".");
    //printf("                              \r");
    printf("Completed Encoding Frame %05d.\r",frame_no);
  }
  // Flush output statistics
  fflush(stdout);

  stats->bit_ctr_n = stats->bit_ctr;

  //Rate control
  if(input->RCEnable) 
  {
    rc_update_pict(bits);

    // update the parameters of quadratic R-D model
    if((img->type==P_SLICE)&&(active_sps->frame_mbs_only_flag))
      updateRCModel();
    else if((img->type==P_SLICE) && (!active_sps->frame_mbs_only_flag) && (img->IFLAG==0))
      updateRCModel();
  }

  stats->bit_ctr_parametersets_n=0;

  if (IMG_NUMBER == 0)
    return 0;
  else
    return 1;
}


/*!
************************************************************************
* \brief
*    This function write out a picture
* \return
*    0 if OK,                                                         \n
*    1 in case of error
*
************************************************************************
*/
static int writeout_picture(Picture *pic)
{
  Bitstream *currStream;
  int partition, slice;
  Slice *currSlice;
  
#ifdef ADAPTIVE_LOOP_FILTER
  if(input->UseAdaptiveLoopFilter)
  {
    AdaptiveLoopFilterProcess(pic);
  }
#endif

  
#ifdef ADAPTIVE_LOOP_FILTER
  if(input->UseAdaptiveLoopFilter)
  {
    TerminateSliceForALF(pic);
  }
#endif

  img->currentPicture=pic;

  for (slice=0; slice<pic->no_slices; slice++)
  {
    currSlice = pic->slices[slice];
    for (partition=0; partition<currSlice->max_part_nr; partition++)
    {
      currStream = (currSlice->partArr[partition]).bitstream;
      assert (currStream->bits_to_go == 8);    //! should always be the case, the 
      //! byte alignment is done in terminate_slice

      // write only if the partition has content
      if (currSlice->partArr[partition].bitstream->write_flag )
        writeUnit (currStream,partition);
    }           // partition loop
  }           // slice loop

  return 0;
}


void copy_params(void)
{
  enc_picture->frame_mbs_only_flag = active_sps->frame_mbs_only_flag;
  enc_picture->frame_cropping_flag = active_sps->frame_cropping_flag;
  enc_picture->chroma_format_idc   = active_sps->chroma_format_idc;

  if (active_sps->frame_cropping_flag)
  {
    enc_picture->frame_cropping_rect_left_offset=active_sps->frame_cropping_rect_left_offset; 
    enc_picture->frame_cropping_rect_right_offset=active_sps->frame_cropping_rect_right_offset; 
    enc_picture->frame_cropping_rect_top_offset=active_sps->frame_cropping_rect_top_offset; 
    enc_picture->frame_cropping_rect_bottom_offset=active_sps->frame_cropping_rect_bottom_offset; 
  }
  else
  {
    enc_picture->frame_cropping_rect_left_offset=0; 
    enc_picture->frame_cropping_rect_right_offset=0; 
    enc_picture->frame_cropping_rect_top_offset=0; 
    enc_picture->frame_cropping_rect_bottom_offset=0; 
  }
}

/*!
************************************************************************
* \brief
*    Encodes a frame picture
************************************************************************
*/
void frame_picture (Picture *frame, int rd_pass)
{
  img->structure = FRAME;
  img->PicSizeInMbs = img->FrameSizeInMbs;

#ifdef SWITCHED_FILTERS
  if((rd_pass == 0) && ((input->UseHPFilter == HPF_SIFO) || (input->UseHPFilter == HPF_SIFO_FPO)))
    img->filterParam = SIFO_FIRST_PASS_FPO;
#endif  // SWITCHED_FILTERS

  if (rd_pass == 2)
  {
#ifdef SWITCHED_FILTERS
    if((enc_frame_picture3) && (input->UseHPFilter == HPF_SIFO))
    {
      free_storable_picture(enc_frame_picture3); 
      free_slice_list(frame_pic_3);
    }
#endif  // SWITCHED_FILTERS

    enc_frame_picture3  = alloc_storable_picture ((PictureStructure)img->structure, img->width, img->height, img->width_cr, img->height_cr);
    img->ThisPOC=enc_frame_picture3->poc=img->framepoc;
    enc_frame_picture3->top_poc    = img->toppoc;
    enc_frame_picture3->bottom_poc = img->bottompoc;

    enc_frame_picture3->frame_poc = img->framepoc;

    enc_frame_picture3->pic_num = img->frame_num;
    enc_frame_picture3->frame_num = img->frame_num;
    enc_frame_picture3->coded_frame = 1;

    enc_frame_picture3->MbaffFrameFlag = img->MbaffFrameFlag = (input->MbInterlace != FRAME_CODING);

    enc_picture=enc_frame_picture3;
    copy_params();
  }
  else if (rd_pass == 1)
  {
#ifdef SWITCHED_FILTERS
    if((enc_frame_picture2) && (input->UseHPFilter == HPF_SIFO))
    {
      free_storable_picture(enc_frame_picture2); 
      free_slice_list(frame_pic_2);
    }
#endif  // SWITCHED_FILTERS

    enc_frame_picture2  = alloc_storable_picture ((PictureStructure)img->structure, img->width, img->height, img->width_cr, img->height_cr);
    img->ThisPOC=enc_frame_picture2->poc=img->framepoc;
    enc_frame_picture2->top_poc    = img->toppoc;
    enc_frame_picture2->bottom_poc = img->bottompoc;

    enc_frame_picture2->frame_poc = img->framepoc;

    enc_frame_picture2->pic_num = img->frame_num;
    enc_frame_picture2->frame_num = img->frame_num;
    enc_frame_picture2->coded_frame = 1;

    enc_frame_picture2->MbaffFrameFlag = img->MbaffFrameFlag = (input->MbInterlace != FRAME_CODING);

    enc_picture=enc_frame_picture2;
    copy_params();
  }
#ifdef ADAPTIVE_FILTER
  else if (rd_pass == 3)    // adaptive filter
  {
    enc_frame_picture_aif  = alloc_storable_picture ((PictureStructure)img->structure, img->width, img->height, img->width_cr, img->height_cr);
    img->ThisPOC=enc_frame_picture_aif->poc=img->framepoc;
    enc_frame_picture_aif->top_poc    = img->toppoc;
    enc_frame_picture_aif->bottom_poc = img->bottompoc;

    enc_frame_picture_aif->frame_poc = img->framepoc;

    enc_frame_picture_aif->pic_num = img->frame_num;
    enc_frame_picture_aif->frame_num = img->frame_num;
    enc_frame_picture_aif->coded_frame = 1;

    enc_frame_picture_aif->MbaffFrameFlag = img->MbaffFrameFlag = (input->MbInterlace != FRAME_CODING);

    enc_picture=enc_frame_picture_aif;
    copy_params();
  }
#endif
#ifdef ADAPTIVE_QUANTIZATION
  else if (rd_pass == 4)    // adaptive quantization matrix
  {
    enc_frame_picture_aqms  = alloc_storable_picture ((PictureStructure)img->structure, img->width, img->height, img->width_cr, img->height_cr);
    img->ThisPOC=enc_frame_picture_aqms->poc=img->framepoc;
    enc_frame_picture_aqms->top_poc    = img->toppoc;
    enc_frame_picture_aqms->bottom_poc = img->bottompoc;

    enc_frame_picture_aqms->frame_poc = img->framepoc;

    enc_frame_picture_aqms->pic_num = img->frame_num;
    enc_frame_picture_aqms->frame_num = img->frame_num;
    enc_frame_picture_aqms->coded_frame = 1;

    enc_frame_picture_aqms->MbaffFrameFlag = img->MbaffFrameFlag = (input->MbInterlace != FRAME_CODING);

    enc_picture=enc_frame_picture_aqms;
    copy_params();
  }
#endif
  else
  {
    enc_frame_picture  = alloc_storable_picture ((PictureStructure)img->structure, img->width, img->height, img->width_cr, img->height_cr);
    img->ThisPOC=enc_frame_picture->poc=img->framepoc;
    enc_frame_picture->top_poc    = img->toppoc;
    enc_frame_picture->bottom_poc = img->bottompoc;

    enc_frame_picture->frame_poc = img->framepoc;

    enc_frame_picture->pic_num = img->frame_num;
    enc_frame_picture->frame_num = img->frame_num;
    enc_frame_picture->coded_frame = 1;

    enc_frame_picture->MbaffFrameFlag = img->MbaffFrameFlag = (input->MbInterlace != FRAME_CODING);

    enc_picture=enc_frame_picture;
    copy_params();
  }


  stats->em_prev_bits_frm = 0;
  stats->em_prev_bits = &stats->em_prev_bits_frm;

  img->fld_flag = 0;

#ifdef ADAPTIVE_FD_SD_CODING
  img->APEC_in_FD_and_SD         = input->APEC_in_FD_and_SD;
#endif
  code_a_picture(frame);

  frame->bits_per_picture = 8 * ((((img->currentSlice)->partArr[0]).bitstream)->byte_pos);

  if (img->structure==FRAME)
  {
    find_distortion ();      
    frame->distortion_y = snr->snr_y;
    frame->distortion_u = snr->snr_u;
    frame->distortion_v = snr->snr_v;
  }
}


/*!
************************************************************************
* \brief
*    Encodes a field picture, consisting of top and bottom field
************************************************************************
*/
void field_picture (Picture *top, Picture *bottom)
{
  //Rate control
  int old_pic_type;              // picture type of top field used for rate control    
  int TopFieldBits;

  //Rate control
  old_pic_type = img->type;

  stats->em_prev_bits_fld = 0;
  stats->em_prev_bits = &stats->em_prev_bits_fld;
  img->number *= 2;
  img->buf_cycle *= 2;
  img->height = (input->img_height+img->auto_crop_bottom) / 2; 
  img->height_cr = img->height_cr_frame / 2;
  img->fld_flag = 1;
  img->PicSizeInMbs = img->FrameSizeInMbs/2;
  // Top field

  enc_top_picture  = alloc_storable_picture ((PictureStructure)img->structure, img->width, img->height, img->width_cr, img->height_cr);
  enc_top_picture->poc=img->toppoc;
  enc_top_picture->frame_poc = img->toppoc;
  enc_top_picture->pic_num = img->frame_num;
  enc_top_picture->frame_num = img->frame_num;
  enc_top_picture->coded_frame = 0;
  enc_top_picture->MbaffFrameFlag = img->MbaffFrameFlag = FALSE;
  img->ThisPOC = img->toppoc;

  img->structure = TOP_FIELD;
  enc_picture = enc_top_picture;
  copy_params();

  put_buffer_top ();
  init_field ();
  if (img->type == B_SLICE)       //all I- and P-frames
    nextP_tr_fld--;


  img->fld_flag = 1;

  //Rate control
  if(input->RCEnable)
  {
    img->BasicUnit=input->basicunit;

    if(input->PicInterlace==FIELD_CODING)
      rc_init_pict(0,1,1); 
    else
      rc_init_pict(0,1,0);

    img->qp  = updateQuantizationParameter(1); 
  }
  img->TopFieldFlag=1;
#ifdef ADAPTIVE_FD_SD_CODING
  img->APEC_in_FD_and_SD         = input->APEC_in_FD_and_SD;
#endif
  code_a_picture(top_pic);
  enc_picture->structure = (PictureStructure)1;

  store_picture_in_dpb(enc_top_picture);

  top->bits_per_picture = 8 * ((((img->currentSlice)->partArr[0]).bitstream)->byte_pos);

  //Rate control
  TopFieldBits=top->bits_per_picture;

  //  Bottom field
  enc_bottom_picture  = alloc_storable_picture ((PictureStructure)img->structure, img->width, img->height, img->width_cr, img->height_cr);
  enc_bottom_picture->poc=img->bottompoc;
  enc_bottom_picture->frame_poc = img->bottompoc;
  enc_bottom_picture->pic_num = img->frame_num;
  enc_bottom_picture->frame_num = img->frame_num;
  enc_bottom_picture->coded_frame = 0;
  enc_bottom_picture->MbaffFrameFlag = img->MbaffFrameFlag = FALSE;
  img->ThisPOC = img->bottompoc;
  img->structure = BOTTOM_FIELD;
  enc_picture = enc_bottom_picture;
  copy_params();
  put_buffer_bot ();
  img->number++;

  init_field ();

  if (img->type == B_SLICE)       //all I- and P-frames
    nextP_tr_fld++;             //check once coding B field

  if (img->type == I_SLICE && input->IntraBottom!=1)
    img->type = (input->BRefPictures == 2) ? B_SLICE : P_SLICE;

  img->fld_flag = 1;

  //Rate control
  if(input->RCEnable)
  {
    setbitscount(TopFieldBits);
    rc_init_pict(0,0,0); 
    img->qp  = updateQuantizationParameter(0); 
  }
  img->TopFieldFlag=0;
#ifdef ADAPTIVE_FD_SD_CODING
  img->APEC_in_FD_and_SD         = input->APEC_in_FD_and_SD;
#endif
  enc_picture->structure = (PictureStructure)2;
  code_a_picture(bottom_pic);

  bottom->bits_per_picture = 8 * ((((img->currentSlice)->partArr[0]).bitstream)->byte_pos);

  // the distortion for a field coded frame (consisting of top and bottom field)
  // lives in the top->distortion variables, the bottom-> are dummies
  distortion_fld (&top->distortion_y, &top->distortion_u, &top->distortion_v);

}


/*!
************************************************************************
* \brief
*    Distortion Field
************************************************************************
*/
static void distortion_fld (float *dis_fld_y, float *dis_fld_u, float *dis_fld_v)
{

  img->number /= 2;
  img->buf_cycle /= 2;
  img->height = (input->img_height+img->auto_crop_bottom);
  img->height_cr = img->height_cr_frame;

  combine_field ();

  imgY_org = imgY_org_frm;
  imgUV_org = imgUV_org_frm;

  find_distortion ();   // find snr from original frame picture

  *dis_fld_y = snr->snr_y;
  *dis_fld_u = snr->snr_u;
  *dis_fld_v = snr->snr_v;
}


/*!
************************************************************************
* \brief
*    Picture Structure Decision
************************************************************************
*/
static int picture_structure_decision (Picture *frame, Picture *top, Picture *bot)
{
  double lambda_picture;
  int bframe = (img->type == B_SLICE);
  float snr_frame, snr_field;
  int bit_frame, bit_field;

  lambda_picture = 0.68 * pow (2, img->bitdepth_lambda_scale + ((img->qp - SHIFT_QP) / 3.0)) * (bframe ? 1 : 1);

  snr_frame = frame->distortion_y + frame->distortion_u + frame->distortion_v;
  //! all distrortions of a field picture are accumulated in the top field
  snr_field = top->distortion_y + top->distortion_u + top->distortion_v;
  bit_field = top->bits_per_picture + bot->bits_per_picture;
  bit_frame = frame->bits_per_picture;

  return decide_fld_frame (snr_frame, snr_field, bit_field, bit_frame, lambda_picture);
}


/*!
************************************************************************
* \brief
*    Field Mode Buffer
************************************************************************
*/
static void field_mode_buffer (int bit_field, float snr_field_y, float snr_field_u, float snr_field_v)
{
  put_buffer_frame ();

  snr->snr_y = snr_field_y;
  snr->snr_u = snr_field_u;
  snr->snr_v = snr_field_v;
}


/*!
************************************************************************
* \brief
*    Frame Mode Buffer
************************************************************************
*/
static void frame_mode_buffer (int bit_frame, float snr_frame_y, float snr_frame_u, float snr_frame_v)
{
  put_buffer_frame ();

  if ((input->PicInterlace != FRAME_CODING)||(input->MbInterlace != FRAME_CODING))
  {
    img->height = img->height / 2;
    img->height_cr = img->height_cr / 2;
    img->number *= 2;

    put_buffer_top ();

    img->number++;
    put_buffer_bot ();

    img->number /= 2;         // reset the img->number to field
    img->height = (input->img_height+img->auto_crop_bottom);
    img->height_cr = img->height_cr_frame;

    snr->snr_y = snr_frame_y;
    snr->snr_u = snr_frame_u;
    snr->snr_v = snr_frame_v;
    put_buffer_frame ();

  }
}


/*!
************************************************************************
* \brief
*    mmco initializations should go here
************************************************************************
*/
static void init_dec_ref_pic_marking_buffer(void)
{
  img->dec_ref_pic_marking_buffer=NULL;
}


/*!
************************************************************************
* \brief
*    Initializes the parameters for a new frame
************************************************************************
*/
static void init_frame (void)
{
  int i;
  int prevP_no, nextP_no;

  last_P_no = last_P_no_frm;

  img->current_mb_nr = 0;
  img->current_slice_nr = 0;
  stats->bit_slice = 0;

  img->mb_y = img->mb_x = 0;
  img->block_y = img->pix_y = img->pix_c_y = 0; 
  img->block_x = img->pix_x = img->block_c_x = img->pix_c_x = 0;

  // The 'slice_nr' of each macroblock is set to -1 here, to guarantee the correct encoding 
  // with FMO (if no FMO, encoding is correct without following assignment), 
  // for which MBs may not be encoded with scan order
  for(i=0;i< ((img->width/MB_BLOCK_SIZE)*(img->height/MB_BLOCK_SIZE));i++)
    img->mb_data[i].slice_nr=-1;

  if (img->b_frame_to_code == 0)
  {
    img->tr = start_tr_in_this_IGOP + IMG_NUMBER * (input->jumpd + 1);

    img->imgtr_last_P_frm = img->imgtr_next_P_frm;
    img->imgtr_next_P_frm = img->tr;

#ifdef _ADAPT_LAST_GROUP_
    if (input->last_frame && img->number + 1 == input->no_frames)
      img->tr = input->last_frame;
#endif

    if (IMG_NUMBER != 0 && input->successive_Bframe != 0)     // B pictures to encode
      nextP_tr_frm = img->tr;

    //Rate control
    if(!input->RCEnable)                  // without using rate control
    {
      if (img->type == I_SLICE)
      {
#ifdef _CHANGE_QP_
        //QP oscillation for secondary SP frames
        if ((input->qp2start > 0 && img->tr >= input->qp2start && input->sp2_frame_indicator==0)||
          ( (input->qp2start > 0) && ( ( (img->tr ) % (2*input->qp2start) ) >=input->qp2start ) && (input->sp2_frame_indicator==1)))
          img->qp = input->qp02;
        else
#endif 
          img->qp = input->qp0;   // set quant. parameter for I-frame
      }
      else
      {
#ifdef _CHANGE_QP_
        //QP oscillation for secondary SP frames
        if ((input->qp2start > 0 && img->tr >= input->qp2start && input->sp2_frame_indicator==0)||
          ( (input->qp2start > 0) && ( ( (img->tr ) % (2*input->qp2start) ) >=input->qp2start ) && (input->sp2_frame_indicator==1)))
          img->qp = input->qpN2 + (img->nal_reference_idc ? 0 : input->DispPQPOffset);
        else
#endif
          img->qp = input->qpN + (img->nal_reference_idc ? 0 : input->DispPQPOffset);

        if (img->type == SP_SLICE)
        {
          if ( (input->qp2start > 0) && ( ( (img->tr ) % (2*input->qp2start) ) >=input->qp2start ))
          { 
            img->qp = input->qpN2-(input->qpN-input->qpsp); 
            img->qpsp = input->qpN2-(input->qpN-input->qpsp_pred);
          }                       
          else
          {
            img->qp = input->qpsp;
            img->qpsp = input->qpsp_pred;
          }
        }   
      }
    }

    img->mb_y_intra = img->mb_y_upd;  //  img->mb_y_intra indicates which GOB to intra code for this frame

    if (input->intra_upd > 0) // if error robustness, find next GOB to update
    {
      img->mb_y_upd = (IMG_NUMBER / input->intra_upd) % (img->height / MB_BLOCK_SIZE);
    }
  }
  else
  {
    img->p_interval = input->jumpd + 1;
    prevP_no = start_tr_in_this_IGOP + (IMG_NUMBER - 1) * img->p_interval;
    nextP_no = start_tr_in_this_IGOP + (IMG_NUMBER) * img->p_interval;

#ifdef _ADAPT_LAST_GROUP_
    last_P_no[0] = prevP_no;
    for (i = 1; i < img->buf_cycle; i++)
      last_P_no[i] = last_P_no[i - 1] - img->p_interval;

    if (input->last_frame && img->number + 1 == input->no_frames)
    {
      nextP_no = input->last_frame;
      img->p_interval = nextP_no - prevP_no;
    }
#endif

    img->b_interval =
      ((double) (input->jumpd + 1) / (input->successive_Bframe + 1.0) );

    if (input->HierarchicalCoding == 3)
      img->b_interval = 1.0;

    if (input->HierarchicalCoding)
      img->tr = prevP_no + (int) (img->b_interval  * (double) (1 + gop_structure[img->b_frame_to_code - 1].display_no));      // from prev_P
    else      
      img->tr = prevP_no + (int) (img->b_interval * (double) img->b_frame_to_code);      // from prev_P


    if (img->tr >= nextP_no)
      img->tr = nextP_no - 1;
    //Rate control
    if(!input->RCEnable && input->HierarchicalCoding == 0)                  // without using rate control   
    {    
#ifdef _CHANGE_QP_
      //QP oscillation for secondary SP frames
      if ((input->qp2start > 0 && img->tr >= input->qp2start && input->sp2_frame_indicator==0)||
        ( (input->qp2start > 0) && ( ( (img->tr ) % (2*input->qp2start) ) >=input->qp2start ) && (input->sp2_frame_indicator==1)))
      {
        img->qp = input->qpB2;
      }
      else
#endif
      {
        img->qp = input->qpB;
      }

      if (img->nal_reference_idc)
      {
#ifdef _CHANGE_QP_
        //QP oscillation for secondary SP frames
        if ((input->qp2start > 0 && img->tr >= input->qp2start && input->sp2_frame_indicator==0)||
          ( (input->qp2start > 0) && ( ( (img->tr ) % (2*input->qp2start) ) >=input->qp2start ) && (input->sp2_frame_indicator==1)))
        {
          img->qp = Clip3(-img->bitdepth_luma_qp_scale,51,input->qpB2 + input->qpBRS2Offset);
        }
        else
#endif
        {
          img->qp = Clip3(-img->bitdepth_luma_qp_scale,51,input->qpB + input->qpBRSOffset);
        }
      }
    }
    else if (input->HierarchicalCoding !=0)  
    {
      // Note that _CHANGE_QP_ does not anymore work for gop_structure. Needs to be fixed
      img->qp =  gop_structure[img->b_frame_to_code - 1].slice_qp;
    }
  }

  UpdateSubseqInfo (img->layer);        // Tian Dong (Sept 2002)
  UpdateSceneInformation ((Boolean)0, 0, 0, -1); // JVT-D099, scene information SEI, nothing included by default

  //! Commented out by StW, needs fixing in SEI.h to keep the trace file clean
  //  PrepareAggregationSEIMessage ();

  img->no_output_of_prior_pics_flag = 0;
  img->long_term_reference_flag = 0;

  init_dec_ref_pic_marking_buffer();
}

/*!
************************************************************************
* \brief
*    Initializes the parameters for a new field
************************************************************************
*/
static void init_field (void)
{
  int i;
  int prevP_no, nextP_no;

  last_P_no = last_P_no_fld;

  img->current_mb_nr = 0;
  img->current_slice_nr = 0;
  stats->bit_slice = 0;

  input->jumpd *= 2;
  input->successive_Bframe *= 2;
  img->number /= 2;
  img->buf_cycle /= 2;

  img->mb_y = img->mb_x = 0;
  img->block_y = img->pix_y = img->pix_c_y = 0; // define vertical positions
  img->block_x = img->pix_x = img->block_c_x = img->pix_c_x = 0;        // define horizontal positions

  if (!img->b_frame_to_code)
  {
    img->tr = img->number * (input->jumpd + 2) + img->fld_type;

    if (!img->fld_type)
    {
      img->imgtr_last_P_fld = img->imgtr_next_P_fld;
      img->imgtr_next_P_fld = img->tr;
    }

#ifdef _ADAPT_LAST_GROUP_
    if (input->last_frame && img->number + 1 == input->no_frames)
      img->tr = input->last_frame;
#endif
    if (img->number != 0 && input->successive_Bframe != 0)    // B pictures to encode
      nextP_tr_fld = img->tr;

    //Rate control
    if(!input->RCEnable)                  // without using rate control
    {
      if (img->type == I_SLICE)
      {
#ifdef _CHANGE_QP_
        //QP oscillation for secondary SP frames
        if ((input->qp2start > 0 && img->tr >= input->qp2start && input->sp2_frame_indicator==0)||
          ( (input->qp2start > 0) && ( ( (img->tr ) % (2*input->qp2start) ) >=input->qp2start ) && (input->sp2_frame_indicator==1)))
          img->qp = input->qp02;
        else
#endif 
          img->qp = input->qp0;   // set quant. parameter for I-frame
      }
      else
      {
#ifdef _CHANGE_QP_
        //QP oscillation for secondary SP frames
        if ((input->qp2start > 0 && img->tr >= input->qp2start && input->sp2_frame_indicator==0)||
          ( (input->qp2start > 0) && ( ( (img->tr ) % (2*input->qp2start) ) >=input->qp2start ) && (input->sp2_frame_indicator==1)))
          img->qp = input->qpN2 + (img->nal_reference_idc ? 0 : input->DispPQPOffset);
        else
#endif
          img->qp = input->qpN + (img->nal_reference_idc ? 0 : input->DispPQPOffset);
        if (img->type == SP_SLICE)
        {
          if ( (input->qp2start > 0) && ( ( (img->tr ) % (2*input->qp2start) ) >=input->qp2start ))
          {
            img->qp = input->qpN2-(input->qpN-input->qpsp); 
            img->qpsp = input->qpN2-(input->qpN-input->qpsp_pred);
          }                       
          else
          {
            img->qp = input->qpsp;
            img->qpsp = input->qpsp_pred;
          }
        }
      }
    }
    img->mb_y_intra = img->mb_y_upd;  //  img->mb_y_intra indicates which GOB to intra code for this frame

    if (input->intra_upd > 0) // if error robustness, find next GOB to update
    {
      img->mb_y_upd =
        (img->number / input->intra_upd) % (img->width / MB_BLOCK_SIZE);
    }
  }
  else
  {
    img->p_interval = input->jumpd + 2;
    prevP_no = (img->number - 1) * img->p_interval + img->fld_type;
    nextP_no = img->number * img->p_interval + img->fld_type;
#ifdef _ADAPT_LAST_GROUP_
    if (!img->fld_type)       // top field
    {
      last_P_no[0] = prevP_no + 1;
      last_P_no[1] = prevP_no;
      for (i = 1; i <= img->buf_cycle; i++)
      {
        last_P_no[2 * i] = last_P_no[2 * i - 2] - img->p_interval;
        last_P_no[2 * i + 1] = last_P_no[2 * i - 1] - img->p_interval;
      }
    }
    else                      // bottom field
    {
      last_P_no[0] = nextP_no - 1;
      last_P_no[1] = prevP_no;
      for (i = 1; i <= img->buf_cycle; i++)
      {
        last_P_no[2 * i] = last_P_no[2 * i - 2] - img->p_interval;
        last_P_no[2 * i + 1] = last_P_no[2 * i - 1] - img->p_interval;
      }
    }

    if (input->last_frame && img->number + 1 == input->no_frames)
    {
      nextP_no = input->last_frame;
      img->p_interval = nextP_no - prevP_no;
    }
#endif
    img->b_interval =
      ((double) (input->jumpd + 1) / (input->successive_Bframe + 1.0) );

    if (input->HierarchicalCoding == 3)
      img->b_interval = 1.0;

    if (input->HierarchicalCoding)
      img->tr = prevP_no + (int) ((img->b_interval + 1.0) * (double) (1 + gop_structure[img->b_frame_to_code - 1].display_no));      // from prev_P
    else      
      img->tr = prevP_no + (int) ((img->b_interval + 1.0) * (double) img->b_frame_to_code);      // from prev_P


    if (img->tr >= nextP_no)
      img->tr = nextP_no - 1; // ?????
    //Rate control
    if(!input->RCEnable && input->HierarchicalCoding == 0)                  // without using rate control
    {
#ifdef _CHANGE_QP_
      //QP oscillation for secondary SP frames
      if ((input->qp2start > 0 && img->tr >= input->qp2start && input->sp2_frame_indicator==0)||
        ( (input->qp2start > 0) && ( ( (img->tr ) % (2*input->qp2start) ) >=input->qp2start ) && (input->sp2_frame_indicator==1)))
      {
        img->qp = input->qpB2;
      }
      else
#endif
        img->qp = input->qpB;
      if (img->nal_reference_idc)
      {
#ifdef _CHANGE_QP_
        //QP oscillation for secondary SP frames
        if ((input->qp2start > 0 && img->tr >= input->qp2start && input->sp2_frame_indicator==0)||
          ( (input->qp2start > 0) && ( ( (img->tr ) % (2*input->qp2start) ) >=input->qp2start ) && (input->sp2_frame_indicator==1)))
        {
          img->qp = Clip3(-img->bitdepth_luma_qp_scale,51,input->qpB2 + input->qpBRS2Offset);
        }
        else
#endif
          img->qp = Clip3(-img->bitdepth_luma_qp_scale,51,input->qpB + input->qpBRSOffset);

      }
    }
    else if (input->HierarchicalCoding != 0)
    {          
      img->qp =  gop_structure[img->b_frame_to_code - 1].slice_qp;
    }
  }
  input->jumpd /= 2;
  input->successive_Bframe /= 2;
  img->buf_cycle *= 2;
  img->number = 2 * img->number + img->fld_type;
}


#define Clip(min,max,val) (((val)<(min))?(min):(((val)>(max))?(max):(val)))


/*!
************************************************************************
* \brief
*    Generate Full Pel Representation
************************************************************************
*/
#ifdef SWITCHED_FILTERS
void GenerateFullPelRepresentation (pel_t ** Fourthpel,
                                    pel_t * Fullpel, int xsize,
                                    int ysize)
#else
static void GenerateFullPelRepresentation (pel_t ** Fourthpel,
                                           pel_t * Fullpel, int xsize,
                                           int ysize)
#endif  // SWITCHED_FILTERS
{
  int x, y, yy , y_pos;

  for (y = 0; y < ysize; y++)
  {
    yy = (y + IMG_PAD_SIZE)<<2;
    y_pos = y * xsize;
    for (x = 0; x < xsize; x++)
      PutPel_11 (Fullpel, y_pos + x, Fourthpel[yy][(x + IMG_PAD_SIZE)<<2]);
  }
}


/*!
************************************************************************
* \brief
*    Upsample 4 times, store them in out4x.  Color is simply copied
*
* \par Input:
*    srcy, srcu, srcv, out4y, out4u, out4v
*
* \par Side Effects_
*    Uses (writes) img4Y_tmp.  This should be moved to a static variable
*    in this module
************************************************************************/
void UnifiedOneForthPix (StorablePicture *s)
{
  int is;
  int i, j, j4;
  int ie2, je2, jj , jpad;
  int maxy = s->size_y + 2 * IMG_PAD_SIZE - 1;
  int ii, i1;
  imgpel **out4Y;
  imgpel  *ref11;
  imgpel  **imgY = s->imgY;
  int size_x_minus1 = s->size_x - 1;
  int size_y_minus1 = s->size_y - 1;
  int ypadded_size = s->size_y + 2 * IMG_PAD_SIZE;

  // don't upsample twice
  if (s->imgY_ups || s->imgY_11)
    return;

  s->imgY_11 = malloc ((s->size_x * s->size_y) * sizeof (imgpel));
  if (NULL == s->imgY_11)
    no_mem_exit("alloc_storable_picture: s->imgY_11");

  get_mem2Dpel (&(s->imgY_ups), (2*IMG_PAD_SIZE + s->size_y)*4, (2*IMG_PAD_SIZE + s->size_x)*4);

#ifdef ADAPTIVE_FILTER       
  //used in adaptive_filter.c
  if(input->UseAdaptiveFilter)
  {
    get_mem2Dpel(&(s->imgY_ups_aif), (2*IMG_PAD_SIZE + s->size_y)*4, (2*IMG_PAD_SIZE + s->size_x)*4);
    if (input->UseAdaptiveFilter == 2) // separable aif
      get_mem2Ddouble(&(s->imgY_ups_aif_hor), 2*IMG_PAD_SIZE + s->size_y, (2*IMG_PAD_SIZE + s->size_x)*4); 
#ifdef E_DAIF
    s->imgY_11_aif = malloc ((s->size_x * s->size_y) * sizeof (imgpel));
    if (NULL == s->imgY_11_aif)
      no_mem_exit("alloc_storable_picture: s->imgY_11_aif");
#endif  // E_DAIF 
  }
#endif  // ADAPTIVE_FILTER  

  if (input->WeightedPrediction || input->WeightedBiprediction || input->GenerateMultiplePPS)
  {
    s->imgY_11_w = malloc ((s->size_x * s->size_y) * sizeof (imgpel));
    if (NULL == s->imgY_11_w)
      no_mem_exit("alloc_storable_picture: s->imgY_11_w");
    get_mem2Dpel (&(s->imgY_ups_w), (2*IMG_PAD_SIZE + s->size_y)*4, (2*IMG_PAD_SIZE + s->size_x)*4);
  }

  out4Y = s->imgY_ups;
  ref11 = s->imgY_11;

  for (j = -IMG_PAD_SIZE; j < s->size_y + IMG_PAD_SIZE; j++)
  {
    jj = max (0, min (size_y_minus1, j));
    jpad = j + IMG_PAD_SIZE;
    for (i = -IMG_PAD_SIZE; i < s->size_x + IMG_PAD_SIZE; i++)
    {      
      is =
        (ONE_FOURTH_TAP[0][0] *
        (imgY[jj][max (0, min (size_x_minus1, i))] +
        imgY[jj][max (0, min (size_x_minus1, i + 1))]) +
        ONE_FOURTH_TAP[1][0] *
        (imgY[jj][max (0, min (size_x_minus1, i - 1))] +
        imgY[jj][max (0, min (size_x_minus1, i + 2))]) +
        ONE_FOURTH_TAP[2][0] *
        (imgY[jj][max (0, min (size_x_minus1, i - 2))] +
        imgY[jj][max (0, min (size_x_minus1, i + 3))]));

      img4Y_tmp[jpad][(i + IMG_PAD_SIZE) * 2] = imgY[jj][max (0, min (size_x_minus1, i))] << 10;    // 1/1 pix pos
      img4Y_tmp[jpad][(i + IMG_PAD_SIZE) * 2 + 1] = is << 5;  // 1/2 pix pos
    }
  }

  for (i = 0; i < (s->size_x + 2 * IMG_PAD_SIZE) * 2; i++)
  {
    ii = i << 1;
#if 1
    // j = 0
    is =
      ( ONE_FOURTH_TAP[0][0] *(img4Y_tmp[0][i] + img4Y_tmp[1][i]) 
      + ONE_FOURTH_TAP[1][0] * (img4Y_tmp[0][i] + img4Y_tmp[2][i]) 
      + ONE_FOURTH_TAP[2][0] * (img4Y_tmp[0][i] + img4Y_tmp[3][i])) >> 5;
    out4Y[0][ii] = (pel_t) Clip3 (0, img->max_imgpel_value, (img4Y_tmp[0][i] + 512) >> 10);  // 1/2 pix
    out4Y[2][ii] = (pel_t) Clip3 (0, img->max_imgpel_value, (is + 512) >> 10);   // 1/2 pix
    // j = 1
    is =
      ( ONE_FOURTH_TAP[0][0] *(img4Y_tmp[1][i] + img4Y_tmp[2][i]) 
      + ONE_FOURTH_TAP[1][0] * (img4Y_tmp[0][i] + img4Y_tmp[3][i]) 
      + ONE_FOURTH_TAP[2][0] * (img4Y_tmp[0][i] + img4Y_tmp[4][i])) >> 5;
    out4Y[4][ii] = (pel_t) Clip3 (0, img->max_imgpel_value, (img4Y_tmp[1][i] + 512) >> 10);  // 1/2 pix
    out4Y[6][ii] = (pel_t) Clip3 (0, img->max_imgpel_value, (is + 512) >> 10);   // 1/2 pix

    for (j = 2; j < maxy - 2 ; j++)
    {
      j4 = j << 2;

      is =
        ( ONE_FOURTH_TAP[0][0] *(img4Y_tmp[j][i] + img4Y_tmp[j + 1][i]) 
        + ONE_FOURTH_TAP[1][0] * (img4Y_tmp[j - 1][i] + img4Y_tmp[j + 2][i]) 
        + ONE_FOURTH_TAP[2][0] * (img4Y_tmp[j - 2][i] + img4Y_tmp[j + 3][i])) >> 5;
      out4Y[j4    ][ii] = (pel_t) Clip3 (0, img->max_imgpel_value, (img4Y_tmp[j][i] + 512) >> 10);  // 1/2 pix
      out4Y[j4 + 2][ii] = (pel_t) Clip3 (0, img->max_imgpel_value, (is + 512) >> 10);   // 1/2 pix
    }

    for (j = maxy - 2; j < ypadded_size ; j++)
    {
      j4 = j << 2;

      is =
        ( ONE_FOURTH_TAP[0][0] *(img4Y_tmp[j][i] + img4Y_tmp[min (maxy, j + 1)][i]) 
        + ONE_FOURTH_TAP[1][0] * (img4Y_tmp[j - 1][i] + img4Y_tmp[min (maxy, j + 2)][i]) 
        + ONE_FOURTH_TAP[2][0] * (img4Y_tmp[j - 2][i] + img4Y_tmp[min (maxy, j + 3)][i])) >> 5;
      out4Y[j4    ][ii] = (pel_t) Clip3 (0, img->max_imgpel_value, (img4Y_tmp[j][i] + 512) >> 10);  // 1/2 pix
      out4Y[j4 + 2][ii] = (pel_t) Clip3 (0, img->max_imgpel_value, (is + 512) >> 10);   // 1/2 pix
    }

#else
    for (j = 0; j < ypadded_size; j++)
      //for (j = 0; j < s->size_y + 2 * IMG_PAD_SIZE; j++)
    {
      j4 = j << 2;

      // change for TML4, use 6 TAP vertical filter
      is =
        ( ONE_FOURTH_TAP[0][0] *(img4Y_tmp[j][i] + img4Y_tmp[min (maxy, j + 1)][i]) 
        + ONE_FOURTH_TAP[1][0] * (img4Y_tmp[max (0, j - 1)][i] + img4Y_tmp[min (maxy, j + 2)][i]) 
        + ONE_FOURTH_TAP[2][0] * (img4Y_tmp[max (0, j - 2)][i] + img4Y_tmp[min (maxy, j + 3)][i])) >> 5;
      //+ ONE_FOURTH_TAP[2][0] * (img4Y_tmp[max (0, j - 2)][i] + img4Y_tmp[min (maxy, j + 3)][i])) / 32;

      //out4Y[j4    ][ii] = (pel_t) Clip3 (0, img->max_imgpel_value, (int) ((img4Y_tmp[j][i] + 512) / 1024));  // 1/2 pix
      //out4Y[j4 + 2][ii] = (pel_t) Clip3 (0, img->max_imgpel_value, (int) ((is + 512) / 1024));   // 1/2 pix
      out4Y[j4    ][ii] = (pel_t) Clip3 (0, img->max_imgpel_value, (int) ((img4Y_tmp[j][i] + 512) >> 10));  // 1/2 pix
      out4Y[j4 + 2][ii] = (pel_t) Clip3 (0, img->max_imgpel_value, (int) ((is + 512) >> 10));   // 1/2 pix
    }
#endif
  }

  /* 1/4 pix */
  /* luma */
  ie2 = (s->size_x + 2 * IMG_PAD_SIZE - 1) * 4 + 2;
  je2 = (s->size_y + 2 * IMG_PAD_SIZE - 1) * 4 + 2;

  for (j = 0; j < je2 + 2; j += 2)
  {
    for (i = 0; i < ie2 - 2; i += 2)
    {
      /*  '-'  */
      //out4Y[j][i+1] = (pel_t) (Clip3 (0, img->max_imgpel_value, (int) (out4Y[j][i] + out4Y[j][i + 2] + 1) >> 1));
      out4Y[j][i+1] = (out4Y[j][i] + out4Y[j][i + 2] + 1) >> 1;
    }

    for (i = ie2 - 2; i < ie2 + 1; i += 2)
    {
      /*  '-'  */
      //out4Y[j][i+1] = (pel_t) (Clip3 (0, img->max_imgpel_value, (int) (out4Y[j][i] + out4Y[j][min (ie2, i + 2)] + 1) >> 1));
      out4Y[j][i+1] = (out4Y[j][i] + out4Y[j][min (ie2, i + 2)] + 1) >> 1;
    }
  }
  for (i = 0; i < ie2 + 2; i++)
  {
    ii = min (ie2, i + 1);
    i1 = i - 1;
    for (j = 0; j < je2 - 2; j += 2)
    {
      if ((i & 0x1) == 0)           /*  '|'  */
      {
        //if((out4Y[j][i] + out4Y[j + 2][i] + 1) >> 1 < 0 || (out4Y[j][i] + out4Y[j + 2][i] + 1) >> 1 > img->max_imgpel_value)
        //printf("Value %d\n",img->max_imgpel_value);
        //out4Y[j + 1][i] = (pel_t) (Clip3 (0, img->max_imgpel_value, (int) (out4Y[j][i] + out4Y[j + 2][i] + 1) >> 1));
        //out4Y[j + 1][i] = (pel_t) ((int) (out4Y[j][i] + out4Y[j + 2][i] + 1) >> 1);
        out4Y[j + 1][i] = (out4Y[j][i] + out4Y[j + 2][i] + 1) >> 1;
      }
      else if (((j & 0x3) == 0 && (i & 0x3) == 1) || ((j & 0x3) == 2 && (i & 0x3) == 3))           /*  '/'  */
      {
        //out4Y[j + 1][i] = (pel_t) (Clip3 (0, img->max_imgpel_value, (int) (out4Y[j][ii] + out4Y[j + 2][i1] + 1) >> 1));
        //out4Y[j + 1][i] = (pel_t) ((int) (out4Y[j][ii] + out4Y[j + 2][i1] + 1) >> 1);
        out4Y[j + 1][i] = (out4Y[j][ii] + out4Y[j + 2][i1] + 1) >> 1;
      }
      else           /*  '\'  */
      {
        //out4Y[j + 1][i] = (pel_t) (Clip3 (0, img->max_imgpel_value, (int) (out4Y[j][i1] + out4Y[j + 2][ii] + 1) >> 1));
        //out4Y[j + 1][i] = (pel_t) ((int) (out4Y[j][i1] + out4Y[j + 2][ii] + 1) >> 1);
        out4Y[j + 1][i] = (out4Y[j][i1] + out4Y[j + 2][ii] + 1) >> 1;
      }
    }

    for (j = je2 - 2; j < je2 + 1; j += 2)
    {
      if ((i & 0x1) == 0)           /*  '|'  */
      {
        //out4Y[j + 1][i] = (pel_t) (Clip3 (0, img->max_imgpel_value, (int) (out4Y[j][i] + out4Y[min (je2, j + 2)][i] + 1) >> 1));
        //out4Y[j + 1][i] = (pel_t) ((int) (out4Y[j][i] + out4Y[min (je2, j + 2)][i] + 1) >> 1);
        out4Y[j + 1][i] = (out4Y[j][i] + out4Y[min (je2, j + 2)][i] + 1) >> 1;
      }
      else if (((j & 0x3) == 0 && (i & 0x3) == 1) || ((j & 0x3) == 2 && (i & 0x3) == 3))           /*  '/'  */
      {
        //out4Y[j + 1][i] = (pel_t) (Clip3 (0, img->max_imgpel_value, (int) (out4Y[j][ii] + out4Y[min (je2, j + 2)][i1] + 1) >> 1));
        //out4Y[j + 1][i] = (pel_t) ((int) (out4Y[j][ii] + out4Y[min (je2, j + 2)][i1] + 1) >> 1);
        out4Y[j + 1][i] = (out4Y[j][ii] + out4Y[min (je2, j + 2)][i1] + 1) >> 1;
      }
      else           /*  '\'  */
      {
        //out4Y[j + 1][i] = (pel_t) (Clip3 (0, img->max_imgpel_value, (int) (out4Y[j][i1] + out4Y[min (je2, j + 2)][ii] + 1) >> 1));
        out4Y[j + 1][i] = (pel_t) ((int) (out4Y[j][i1] + out4Y[min (je2, j + 2)][ii] + 1) >> 1);
      }
    }

  }

  // Generate 1/1th pel representation (used for integer pel MV search)
  GenerateFullPelRepresentation (out4Y, ref11, s->size_x, s->size_y);    
}


#ifdef USE_HP_FILTER
void UnifiedOneForthPixEnhancedFDIF (StorablePicture *s)
{
  imgpel **out4Y;
  imgpel  *ref11;
  imgpel  **imgY = s->imgY;

  int i,j,k,l,i2,j2;
  int sfp;
  static const int f0[4][4] = {{0,5,5,0},{5,22,22,5},{5,22,22,5},{0,5,5,0}};
  static const int f1[6] = {3,-15,111,37,-10,2};
  static const int f2[6] = {3,-17,78,78,-17,3};
  static const int f3[6] = {2,-10,37,111,-15,3};
  int d11,d12,d13,d21,d22,d23;
  int h11,h12,h13;
  int v11,v12,v13;

  // don't upsample twice
  if (s->imgY_ups || s->imgY_11)
    return;

  s->imgY_11 = malloc ((s->size_x * s->size_y) * sizeof (imgpel));
  if (NULL == s->imgY_11)
    no_mem_exit("alloc_storable_picture: s->imgY_11");

  get_mem2Dpel (&(s->imgY_ups), (2*IMG_PAD_SIZE + s->size_y)*4, (2*IMG_PAD_SIZE + s->size_x)*4);
#ifdef ADAPTIVE_FILTER       
  //used in adaptive_filter.c
  if(input->UseAdaptiveFilter)
  {
    get_mem2Dpel(&(s->imgY_ups_aif), (2*IMG_PAD_SIZE + s->size_y)*4, (2*IMG_PAD_SIZE + s->size_x)*4);
    if (input->UseAdaptiveFilter == 2) // separable aif
      get_mem2Ddouble(&(s->imgY_ups_aif_hor), 2*IMG_PAD_SIZE + s->size_y, (2*IMG_PAD_SIZE + s->size_x)*4); 
  }
#endif

  if (input->WeightedPrediction || input->WeightedBiprediction || input->GenerateMultiplePPS)
  {
    s->imgY_11_w = malloc ((s->size_x * s->size_y) * sizeof (imgpel));
    if (NULL == s->imgY_11_w)
      no_mem_exit("alloc_storable_picture: s->imgY_11_w");
    get_mem2Dpel (&(s->imgY_ups_w), (2*IMG_PAD_SIZE + s->size_y)*4, (2*IMG_PAD_SIZE + s->size_x)*4);
  }
  out4Y = s->imgY_ups;
  ref11 = s->imgY_11;

  for (i = 0; i < s->size_y+2*IMG_PAD_SIZE; i++) 
  {            
    for (j = 0; j < s->size_x+2*IMG_PAD_SIZE; j++) 
    {			
      /*! Start Integer position */
      i2 = max(0,min(s->size_y-1,i-IMG_PAD_SIZE));
      j2 = max(0,min(s->size_x-1,j-IMG_PAD_SIZE));
      out4Y[i*4][j*4] = imgY[i2][j2];

      /*! Start Fractional positions */
      h11 = h12 = h13 = v11 = v12 = v13 = d11 = d12 = d13 = d21 = d22 = d23 = 0;

      /*! Strong filter position */
      sfp = 0;
      for (l=0;l<4;l++)
      {
        for (k=0;k<4;k++)
        {
          i2 = max(0,min(s->size_y-1,i-IMG_PAD_SIZE-1+l));
          j2 = max(0,min(s->size_x-1,j-IMG_PAD_SIZE-1+k));
          sfp += f0[l][k] * imgY[i2][j2];
        }
      }

      /*! Temporary variables */
      for (k=0;k<6;k++)
      {
        i2 = max(0,min(s->size_y-1,i-IMG_PAD_SIZE));
        j2 = max(0,min(s->size_x-1,j-IMG_PAD_SIZE-2+k));
        h11 += f1[k] * imgY[i2][j2];
        h12 += f2[k] * imgY[i2][j2];
        h13 += f3[k] * imgY[i2][j2];

        i2 = max(0,min(s->size_y-1,i-IMG_PAD_SIZE-2+k));
        j2 = max(0,min(s->size_x-1,j-IMG_PAD_SIZE));
        v11 += f1[k] * imgY[i2][j2];
        v12 += f2[k] * imgY[i2][j2];
        v13 += f3[k] * imgY[i2][j2];

        i2 = max(0,min(s->size_y-1,i-IMG_PAD_SIZE-2+k));
        j2 = max(0,min(s->size_x-1,j-IMG_PAD_SIZE-2+k));
        d11 += f1[k] * imgY[i2][j2];
        d12 += f2[k] * imgY[i2][j2];
        d13 += f3[k] * imgY[i2][j2];

        i2 = max(0,min(s->size_y-1,i-IMG_PAD_SIZE-2+k));
        j2 = max(0,min(s->size_x-1,j-IMG_PAD_SIZE+3-k));
        d21 += f1[k] * imgY[i2][j2];
        d22 += f2[k] * imgY[i2][j2];
        d23 += f3[k] * imgY[i2][j2];
      }

      /*! Compute fractional position values from temporary variables above */
      out4Y[i*4+0][j*4+1] = (pel_t) Clip3(0,img->max_imgpel_value,(h11    + 64)>>7); //a                           
      out4Y[i*4+0][j*4+2] = (pel_t) Clip3(0,img->max_imgpel_value,(h12    + 64)>>7); //b                          
      out4Y[i*4+0][j*4+3] = (pel_t) Clip3(0,img->max_imgpel_value,(h13    + 64)>>7); //c                          
      out4Y[i*4+1][j*4+0] = (pel_t) Clip3(0,img->max_imgpel_value,(v11    + 64)>>7); //d
      out4Y[i*4+1][j*4+1] = (pel_t) Clip3(0,img->max_imgpel_value,(d11+128+ 64)>>7); //e - Local DC Offset (+)
      out4Y[i*4+1][j*4+2] = (pel_t) Clip3(0,img->max_imgpel_value,(d11+d21+128)>>8); //f
      out4Y[i*4+1][j*4+3] = (pel_t) Clip3(0,img->max_imgpel_value,(sfp    + 64)>>7); //g - Strong Filter Position
      out4Y[i*4+2][j*4+0] = (pel_t) Clip3(0,img->max_imgpel_value,(v12    + 64)>>7); //h
      out4Y[i*4+2][j*4+1] = (pel_t) Clip3(0,img->max_imgpel_value,(d11+d23+128)>>8); //i
      out4Y[i*4+2][j*4+2] = (pel_t) Clip3(0,img->max_imgpel_value,(d12+d22+128)>>8); //j
      out4Y[i*4+2][j*4+3] = (pel_t) Clip3(0,img->max_imgpel_value,(d13+d21+128)>>8); //k
      out4Y[i*4+3][j*4+0] = (pel_t) Clip3(0,img->max_imgpel_value,(v13    + 64)>>7); //l
      out4Y[i*4+3][j*4+1] = (pel_t) Clip3(0,img->max_imgpel_value,(d23-128+ 64)>>7); //m - Local DC Offset (-)   
      out4Y[i*4+3][j*4+2] = (pel_t) Clip3(0,img->max_imgpel_value,(d13+d23+128)>>8); //n    
      out4Y[i*4+3][j*4+3] = (pel_t) Clip3(0,img->max_imgpel_value,(d13    + 64)>>7); //o
    } //for j
  } //for i

  // Generate 1/1th pel representation (used for integer pel MV search)
  GenerateFullPelRepresentation (out4Y, ref11, s->size_x, s->size_y);    
}

short int sub_loc[15][2] =
{ //a     b     c     d     e     f     g     h     i     j     k     l     m     n     o
  {1,0},{2,0},{3,0},{0,1},{1,1},{2,1},{3,1},{0,2},{1,2},{2,2},{3,2},{0,3},{1,3},{2,3},{3,3}
  //0     1     2     3     4     5     6     7     8     9     10    11    12   13     14
};

short int acdl_idx[4]={0, 2, 3, 11};
short int egmo_idx[4]={4, 6, 12, 14};
short int fikn_idx[4]={5, 8, 10, 13};

short int acdl_f[4][2]  = {{0,0},{4,0},{0,0},{0,4}}; //contribution full pel location 
short int acdl_h[4][2]  = {{2,0},{2,0},{0,2},{0,2}}; //contributing half pel location {x,y}
short int egmo[4][2][2] =                            //contributing half pel location
{
  {{2,0},{0,2}},  //"e"
  {{2,0},{4,2}},  //"g"
  {{2,4},{0,2}},  //"m"
  {{2,4},{4,2}}   //"o"
}; 

short int fikn_q[4][2]  = {{2,2},{2,2},{2,2},{2,2}}; //contributing (1/2, 1/2) pel "j"
short int fikn_h[4][2]  = {{2,0},{0,2},{4,2},{2,4}}; //contributing half pel location

int get_mem2Dshortint(short int ***array2D, int rows, int columns)
{
  int i;

  if((*array2D      = (short int**)calloc(rows,        sizeof(short int*))) == NULL)
    no_mem_exit("get_mem2Dshortint: array2D");
  if(((*array2D)[0] = (short int* )calloc(rows*columns,sizeof(short int ))) == NULL)
    no_mem_exit("get_mem2Dshortint: array2D");

  for(i=1 ; i<rows ; i++)
    (*array2D)[i] =  (*array2D)[i-1] + columns  ;

  return rows*columns*sizeof(short int);
}
void free_mem2Dshortint(short int **array2D)
{
  if (array2D)
  {
    if (array2D[0]) 
      free (array2D[0]);
    else error ("free_mem2Dshortint: trying to free unused memory",100);

    free (array2D);

  } else
  {
    error ("free_mem2Dshortint: trying to free unused memory",100);
  }
}

void UnifiedOneForthPixWithHighPrecisionH264Filter_16 (StorablePicture *s)
{
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
  int is, AF[3], full_neighbor, half1_neighbor, half2_neighbor, half_neighbor, temp;
#else
  short int is, AF[3], full_neighbor, half1_neighbor, half2_neighbor, half_neighbor, temp;
#endif
  int i, j, sub_i;
  int ii, jj;
  imgpel **out4Y;
  imgpel  *ref11;
  imgpel  **imgY = s->imgY;
  int pos_y1, pos_y2, pos_x1, pos_x2; 
  int x_sub, y_sub, sub_pos; 
  int x_max, y_max, tmp_x, tmp_y;

#ifdef  INTERNAL_BIT_DEPTH_INCREASE
  int **tmp4Y; 
#else
  short int **tmp4Y; 
#endif

  // don't upsample twice
  if (s->imgY_ups || s->imgY_11)
    return;

  s->imgY_11 = malloc ((s->size_x * s->size_y) * sizeof (imgpel));
  if (NULL == s->imgY_11)
    no_mem_exit("alloc_storable_picture: s->imgY_11");

  get_mem2Dpel (&(s->imgY_ups), (2*IMG_PAD_SIZE + s->size_y)*4, (2*IMG_PAD_SIZE + s->size_x)*4);

#ifdef  INTERNAL_BIT_DEPTH_INCREASE
  get_mem2Dint (&tmp4Y, (2*IMG_PAD_SIZE + s->size_y)*4, (2*IMG_PAD_SIZE + s->size_x)*4);
#else
  get_mem2Dshortint (&tmp4Y, (2*IMG_PAD_SIZE + s->size_y)*4, (2*IMG_PAD_SIZE + s->size_x)*4);
#endif

#ifdef ADAPTIVE_FILTER       
  //used in adaptive_filter.c
  if(input->UseAdaptiveFilter)
  {
#ifdef E_DAIF
    s->imgY_11_aif = malloc ((s->size_x * s->size_y) * sizeof (imgpel));
    if (NULL == s->imgY_11_aif)
      no_mem_exit("alloc_storable_picture: s->imgY_11_aif");
#endif  // E_DAIF 
    get_mem2Dpel(&(s->imgY_ups_aif), (2*IMG_PAD_SIZE + s->size_y)*4, (2*IMG_PAD_SIZE + s->size_x)*4);
    if (input->UseAdaptiveFilter == 2) // separable aif
      get_mem2Ddouble(&(s->imgY_ups_aif_hor), 2*IMG_PAD_SIZE + s->size_y, (2*IMG_PAD_SIZE + s->size_x)*4);
  }
#endif  // ADAPTIVE_FILTER 

  if (input->WeightedPrediction || input->WeightedBiprediction || input->GenerateMultiplePPS)
  {
    s->imgY_11_w = malloc ((s->size_x * s->size_y) * sizeof (imgpel));
    if (NULL == s->imgY_11_w)
      no_mem_exit("alloc_storable_picture: s->imgY_11_w");
    get_mem2Dpel (&(s->imgY_ups_w), (2*IMG_PAD_SIZE + s->size_y)*4, (2*IMG_PAD_SIZE + s->size_x)*4);
  }
  out4Y = s->imgY_ups;
  ref11 = s->imgY_11;

  // interpolate 1D half-pel first (but do not round them down yet)
  for(i=-4*IMG_PAD_SIZE; i < (IMG_PAD_SIZE+s->size_y)*4; i+=2)
    for(j=-4*IMG_PAD_SIZE; j < (IMG_PAD_SIZE+s->size_x)*4; j+=2)
    {
      x_sub = (j+4*IMG_PAD_SIZE)%4;       // x-sub-coordinate in a 4x4block
      y_sub = (i+4*IMG_PAD_SIZE)%4;       // y-sub-coordinate in a 4x4block
      sub_pos = x_sub + 4*y_sub;          // pos 1..15 in a 4x4 block
      is = 0;
      if(sub_pos)
      {
        if(sub_pos/4 == 0)    // horizontal 1/2 position 
        {
          pos_y1 = max(0, min(s->size_y-1, i/4 ) );

          for(jj = 0; jj < 3; jj++)
          {
            pos_x1 = max(0, min(s->size_x -1, j/4-2+jj ) );
            pos_x2 = max(0, min(s->size_x -1, j/4-2+(5-jj) ) );
            if(j<0 && pos_x2 > 0)
              pos_x2 -=1;
            is += ONE_FOURTH_TAP[2-jj][0]*(imgY[pos_y1][pos_x1]+imgY[pos_y1][pos_x2]);
          }
          tmp4Y[i+4*IMG_PAD_SIZE][j+4*IMG_PAD_SIZE] = is;
        }
        else if(sub_pos%4 == 0) // vertical 1/2 position
        {
          pos_x1 = max(0, min(s->size_x -1, j/4 ) );

          for(ii = 0; ii < 3; ii++)
          {
            pos_y1 = max(0, min(s->size_y-1, i/4-2+ii ) );
            pos_y2 = max(0, min(s->size_y-1, i/4-2+(5-ii) ) );
            if(i<0 && pos_y2 > 0)
              pos_y2 -=1;
            is += ONE_FOURTH_TAP[2-ii][0]*(imgY[pos_y1][pos_x1]+imgY[pos_y2][pos_x1]);
          }
          tmp4Y[i+4*IMG_PAD_SIZE][j+4*IMG_PAD_SIZE] = is;
        }
      }
      else
        tmp4Y[i+4*IMG_PAD_SIZE][j+4*IMG_PAD_SIZE] =
        imgY[max(0,min(s->size_y-1,i/4))][max(0,min(s->size_x-1,j/4))];
    }

    // fill the (1/2,1/2) positions
    for(i=-4*IMG_PAD_SIZE+2; i < (IMG_PAD_SIZE+s->size_y)*4; i+=4)
      for(j=-4*IMG_PAD_SIZE+2; j < (IMG_PAD_SIZE+s->size_x)*4; j+=4)
      {
        is = 0;
        //AF[0]=A+F; AF[1]=B+E; AF[2]=C+D;
        pos_x1 = j;
        for(ii = 0; ii < 3; ii++)
        {
          pos_y1 = max(-4*IMG_PAD_SIZE, min((s->size_y+IMG_PAD_SIZE-1)*4, (i-2)-2*4+4*ii ) ); //max(0, min((s->size_y)*4, (i-2)-2*4+4*ii ) );
          pos_y2 = max(-4*IMG_PAD_SIZE, min((s->size_y+IMG_PAD_SIZE-1)*4, (i-2)-2*4+4*(5-ii) ) ); //max(0, min((s->size_y)*4, (i-2)-2*4+4*(5-ii) ) );        

          AF[ii] = tmp4Y[pos_y1+4*IMG_PAD_SIZE][pos_x1+4*IMG_PAD_SIZE]+tmp4Y[pos_y2+4*IMG_PAD_SIZE][pos_x1+4*IMG_PAD_SIZE];
        }

        AF[0] -= AF[1]; //(A+F)-(B+E)
        AF[0] >>= 2;
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        if (img->BitDepthIncrease==0)
        {
          AF[2] = max(-2371, min(18640, AF[2]));
        }
#else
        AF[2] = max(-2371, min(18640, AF[2]));
#endif
        AF[1] = AF[2] - AF[1];   // (C+D)-(B+E)
        AF[0] += AF[1];
        AF[0] >>= 2;
        AF[0] += AF[2];

        tmp4Y[i+4*IMG_PAD_SIZE][j+4*IMG_PAD_SIZE] = AF[0];
      }


      // interpolate quarter-pel and shift down half-pel/full-pel at the same time
      x_max = (s->size_x+2*IMG_PAD_SIZE)*4-1;//(s->size_x+2*IMG_PAD_SIZE)*4-2;
      y_max = (s->size_y+2*IMG_PAD_SIZE)*4-1;//(s->size_y+2*IMG_PAD_SIZE)*4-2;
      for(i=-4*IMG_PAD_SIZE; i < (IMG_PAD_SIZE+s->size_y)*4; i+=4)
        for(j=-4*IMG_PAD_SIZE; j < (IMG_PAD_SIZE+s->size_x)*4; j+=4)
        {
          int full_x, full_y, half1_x, half1_y, half2_x, half2_y, half_x, half_y, sub_x, sub_y;
          // full-pel
          temp = tmp4Y[i+4*IMG_PAD_SIZE][j+4*IMG_PAD_SIZE ];
          out4Y[i+4*IMG_PAD_SIZE][j+4*IMG_PAD_SIZE] = 
            (pel_t)Clip3(0,(1<<img->bitdepth_luma)-1,temp); 

          // a, c, d, l
          for(sub_i=0; sub_i<4; sub_i++)
          {
            sub_x = sub_loc[acdl_idx[sub_i]][0];
            sub_y = sub_loc[acdl_idx[sub_i]][1];

            full_x =acdl_f[sub_i][0];
            full_y =acdl_f[sub_i][1];
            tmp_x = j+4*IMG_PAD_SIZE + full_x;
            if(tmp_x >x_max) 
              tmp_x = x_max -3;
            tmp_y = i+4*IMG_PAD_SIZE + full_y;
            if(tmp_y > y_max)
              tmp_y = y_max -3;

            full_neighbor = tmp4Y[tmp_y][tmp_x];

            half_x = acdl_h[sub_i][0];
            half_y = acdl_h[sub_i][1];
            tmp_x = j+4*IMG_PAD_SIZE + half_x;
            if(tmp_x >x_max)
              tmp_x = x_max -3;
            tmp_y = i+4*IMG_PAD_SIZE + half_y;
            if(tmp_y > y_max)
              tmp_y = y_max -3;

            half_neighbor = tmp4Y[tmp_y][tmp_x];

            temp = ((full_neighbor<<5) + half_neighbor + 32) >> 6;
            out4Y[i+4*IMG_PAD_SIZE +sub_y][j+4*IMG_PAD_SIZE+sub_x] = 
              (pel_t)Clip3(0,(1<<img->bitdepth_luma)-1,temp); 
          }                

          // h,
          temp = (tmp4Y[i+4*IMG_PAD_SIZE + 2][j+4*IMG_PAD_SIZE] + 16)>>5;
          out4Y[i+4*IMG_PAD_SIZE + 2][j+4*IMG_PAD_SIZE] = 
            (pel_t)Clip3(0,(1<<img->bitdepth_luma)-1,temp); 

          // b 
          temp = (tmp4Y[i+4*IMG_PAD_SIZE][j+4*IMG_PAD_SIZE + 2] + 16)>>5;
          out4Y[i+4*IMG_PAD_SIZE][j+4*IMG_PAD_SIZE + 2] = 
            (pel_t)Clip3(0,(1<<img->bitdepth_luma)-1,temp); 

          // j
          temp = (tmp4Y[i+4*IMG_PAD_SIZE + 2][j+4*IMG_PAD_SIZE + 2] + 32)>>6;
          out4Y[i+4*IMG_PAD_SIZE + 2][j+4*IMG_PAD_SIZE + 2] = 
            (pel_t)Clip3(0,(1<<img->bitdepth_luma)-1,temp); 

          // e, g, m, o
          for(sub_i=0; sub_i<4; sub_i++)
          {
            sub_x = sub_loc[egmo_idx[sub_i]][0];
            sub_y = sub_loc[egmo_idx[sub_i]][1];

            half1_x =egmo[sub_i][0][0];
            half1_y =egmo[sub_i][0][1];
            tmp_x = j+4*IMG_PAD_SIZE + half1_x;
            tmp_y = i+4*IMG_PAD_SIZE + half1_y;
            if(tmp_x >x_max)
              tmp_x = x_max -3;
            if(tmp_y > y_max)
              tmp_y = y_max -3;

            half1_neighbor = tmp4Y[tmp_y][tmp_x];

            half2_x = egmo[sub_i][1][0];
            half2_y = egmo[sub_i][1][1];
            tmp_x = j+4*IMG_PAD_SIZE + half2_x;
            tmp_y = i+4*IMG_PAD_SIZE + half2_y;
            if(tmp_x >x_max)
              tmp_x = x_max -3;
            if(tmp_y > y_max)
              tmp_y = y_max -3;
            half2_neighbor = tmp4Y[tmp_y][tmp_x];

            temp = (half1_neighbor + half2_neighbor + 32) >> 6;

            out4Y[i+4*IMG_PAD_SIZE +sub_y][j+4*IMG_PAD_SIZE+sub_x] = 
              (pel_t)Clip3(0,(1<<img->bitdepth_luma)-1,temp); 
          }                

          // i, k,f, n            
          for(sub_i=0; sub_i<4; sub_i++)
          {
            sub_x = sub_loc[fikn_idx[sub_i]][0];
            sub_y = sub_loc[fikn_idx[sub_i]][1];

            half_x = fikn_h[sub_i][0];
            half_y = fikn_h[sub_i][1];
            tmp_x = j+4*IMG_PAD_SIZE + half_x;
            tmp_y = i+4*IMG_PAD_SIZE + half_y;
            if(tmp_x >x_max)
              tmp_x = x_max -3;
            if(tmp_y > y_max)
              tmp_y = y_max -3;

            half_neighbor = tmp4Y[tmp_y][tmp_x];

            temp = (half_neighbor + (tmp4Y[i+4*IMG_PAD_SIZE + 2][j+4*IMG_PAD_SIZE + 2]>>1) + 32) >> 6;

            out4Y[i+4*IMG_PAD_SIZE +sub_y][j+4*IMG_PAD_SIZE+sub_x] = 
              (pel_t)Clip3(0,(1<<img->bitdepth_luma)-1,temp); 
          }                
        }

#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        free_mem2Dint(tmp4Y); 
#else
        free_mem2Dshortint(tmp4Y); 
#endif

        // Generate 1/1th pel representation (used for integer pel MV search)
        GenerateFullPelRepresentation (out4Y, ref11, s->size_x, s->size_y);    
}
#endif

#ifdef USE_NEW_OFFSET
void SetFrameOffset(void)
{
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int x, y, xj, yi, temp, valOrg;
  int mvx=0,  mvy=0;
  int ref_frame=0;
  int subblock=0;
  int x_orig, y_orig;       
  int x_pos, y_pos;
  int out4Y_width  = (img->width+2*IMG_PAD_SIZE)*4-1;  
  int out4Y_height = (img->height+2*IMG_PAD_SIZE)*4-1;

  if(IS_INTRA(currMB)) //intra macroblocks are not used for calculation of the filter coeffs.
    return;

  x_orig = MB_BLOCK_SIZE*(img->block_x/4);
  y_orig = MB_BLOCK_SIZE*(img->block_y/4);
  for(subblock = 0; subblock < 16; subblock++)
  {
    x = x_orig+4*(subblock%4);
    y = y_orig+4*(subblock/4);
    mvx = enc_picture->mv[LIST_0][y/4][x/4][0];
    mvy = enc_picture->mv[LIST_0][y/4][x/4][1];
    ref_frame = enc_picture->ref_idx[LIST_0][y/4][x/4];

    if(ref_frame != -1)
    {
      for(yi = 0; yi < 4; yi++)
      {    //y
        for(xj = 0; xj < 4; xj++)
        {  //x
          valOrg   = imgY_org[y+yi][x+xj];

          y_pos = max(0,min(out4Y_height,4*(y+yi)+4*IMG_PAD_SIZE+mvy));
          x_pos = max(0,min(out4Y_width, 4*(x+xj)+4*IMG_PAD_SIZE+mvx));
          temp=listX[LIST_0][ref_frame]->imgY_ups[y_pos][x_pos];
          img->frameOffsetTotal[LIST_0][ref_frame]+=(valOrg-temp);
          img->frameOffsetCount[LIST_0][ref_frame]++;          

#ifdef USE_LMS     
          img->ref_sum       [LIST_0][ref_frame] += temp;
          img->ref_square_sum[LIST_0][ref_frame] += temp*temp;
          img->ref_cross     [LIST_0][ref_frame] += temp*valOrg;
          img->curr_sum      [LIST_0][ref_frame] += valOrg;
#endif
        }
      }
    }

    mvx = enc_picture->mv[LIST_1][y/4][x/4][0];
    mvy = enc_picture->mv[LIST_1][y/4][x/4][1];
    ref_frame = enc_picture->ref_idx[LIST_1][y/4][x/4];
    if(ref_frame != -1)
    {
      for(yi = 0; yi < 4; yi++)
      {    //y
        for(xj = 0; xj < 4; xj++)
        {  //x
          valOrg   = imgY_org[y+yi][x+xj];

          y_pos = max(0,min(out4Y_height,4*(y+yi)+4*IMG_PAD_SIZE+mvy));
          x_pos = max(0,min(out4Y_width, 4*(x+xj)+4*IMG_PAD_SIZE+mvx));
          temp=listX[LIST_1][ref_frame]->imgY_ups[y_pos][x_pos];
          img->frameOffsetTotal[LIST_1][ref_frame]+=(valOrg-temp);
          img->frameOffsetCount[LIST_1][ref_frame]++;          

#ifdef USE_LMS     
          img->ref_sum        [LIST_1][ref_frame] += temp;
          img->ref_square_sum [LIST_1][ref_frame] += temp*temp;
          img->ref_cross      [LIST_1][ref_frame] += temp*valOrg;
          img->curr_sum       [LIST_1][ref_frame] += valOrg;
#endif
        }
      }
    }
  }
}

void CalcFrameOffset(void)
{
  int frame, list, offset; 
  double temp;
  int numlists  = (img->type == B_SLICE) ? 2 : 1;
#ifdef USE_LMS     
  int default_weight = 1<<6;
  double ref_sum, ref_square_sum, curr_sum, ref_cross;
#endif

  for(list = 0; list < numlists; list++)
    for(frame = 0; frame < listXsize[list]; frame++)
    {
#ifdef USE_LMS     
      img->frameOffset_mse[list][frame] = 0;
      img->frameWeight_mse[list][frame] = default_weight;
#endif

      temp=(double)img->frameOffsetTotal[list][frame];

      if (img->frameOffsetCount[list][frame]>0)
      {
        offset=(int)(fabs(temp)/(double)img->frameOffsetCount[list][frame]+0.5);
        if (img->frameOffsetTotal[list][frame]>=0)
        {
          img->frameOffset[list][frame]=offset;
        }
        else
        {
          img->frameOffset[list][frame]=-offset;
        }

#ifdef USE_LMS     
        ref_sum = (double)img->ref_sum[list][frame] / img->frameOffsetCount[list][frame];
        ref_square_sum = (double)img->ref_square_sum[list][frame]/img->frameOffsetCount[list][frame];
        curr_sum = (double)img->curr_sum[list][frame]/img->frameOffsetCount[list][frame];
        ref_cross = (double)img->ref_cross[list][frame]/img->frameOffsetCount[list][frame];

        temp = ref_sum * ref_sum - ref_square_sum;
        if((int)temp != 0)
        {
          img->frameWeight_mse[list][frame] = (int)(default_weight*(ref_sum * curr_sum - ref_cross)/temp + 0.5);
          img->frameOffset_mse[list][frame] = (int)((-ref_square_sum * curr_sum + ref_sum * ref_cross)/temp + 0.5);
        }
#endif
      }
    }
}

void ResetFrameOffset(void)
{
  int frame, list;

  for(list = 0; list < 2; list++)
    for(frame = 0; frame < MAX_REFERENCE_PICTURES; frame++)
    {
      img->frameOffsetTotal[list][frame] = 0;
      img->frameOffsetCount[list][frame] = 0;

#ifdef USE_LMS     
      img->ref_sum[list][frame] = 0;
      img->ref_square_sum[list][frame] = 0;
      img->ref_cross[list][frame] = 0;
      img->curr_sum[list][frame] = 0;
#endif
    }
}

#endif



#ifdef EIGHTH_PEL
static int h1_25[8] = { -3, 12, -37, 229,  71, -21,  6, -1 };
static int h2_25[8] = { -3, 12, -39, 158, 158, -39, 12, -3 };
static int h3_25[8] = { -1,  6, -21,  71, 229, -37, 12, -3 };
/*!
************************************************************************
* \brief
*    Upsample 4 times for 1/8-pel estimation and store in buffer
*    for multiple reference frames. 1/8-pel resolution is calculated
*    during the motion estimation on the fly with bilinear interpolation.
*
************************************************************************
*/
void oneeighthpix(StorablePicture *s)
{
  int *h1, *h2, *h3;
  int x,y,y1,x4,y4, maxy;
  int IMG_PAD_SIZE_TIMES2 = IMG_PAD_SIZE << 1;

  int i0,i1,i2,i3;

  int rounding = 32768;
  int shifting = 16;

  imgpel  **imgY = s->imgY;
  imgpel  **out4Y;
  imgpel  *ref11;

  int nx_1 = s->size_x - 1;
  int ny_1 = s->size_y - 1;

  h1 = h1_25;
  h2 = h2_25;
  h3 = h3_25;



  // don't upsample twice
  if (s->imgY_ups || s->imgY_11)
    return;

  s->imgY_11 = malloc ((s->size_x * s->size_y) * sizeof (imgpel));
  if (NULL == s->imgY_11)
    no_mem_exit("alloc_storable_picture: s->imgY_11");

  get_mem2Dpel (&(s->imgY_ups), (IMG_PAD_SIZE_TIMES2 + s->size_y) << 2, (IMG_PAD_SIZE_TIMES2 + s->size_x) << 2);

#ifdef ADAPTIVE_FILTER
  //used in adaptive_filter.c
  get_mem2Dpel(&(s->imgY_ups_aif), (IMG_PAD_SIZE_TIMES2 + s->size_y)<< 2, (IMG_PAD_SIZE_TIMES2+ s->size_x)<< 2);
  if (input->UseAdaptiveFilter == 2) // separable aif
    get_mem2Ddouble(&(s->imgY_ups_aif_hor), IMG_PAD_SIZE_TIMES2 + s->size_y, (IMG_PAD_SIZE_TIMES2+ s->size_x)<< 2); 
#endif


  if (input->WeightedPrediction || input->WeightedBiprediction || input->GenerateMultiplePPS)
  {
    s->imgY_11_w = malloc ((s->size_x * s->size_y) * sizeof (imgpel));
    if (NULL == s->imgY_11_w)
      no_mem_exit("alloc_storable_picture: s->imgY_11_w");
    get_mem2Dpel (&(s->imgY_ups_w), (IMG_PAD_SIZE_TIMES2 + s->size_y) << 2, (IMG_PAD_SIZE_TIMES2 + s->size_x) << 2);
  }

  out4Y = s->imgY_ups;
  ref11 = s->imgY_11;

  //horizontal filtering
  for( y = -IMG_PAD_SIZE; y < s->size_y + IMG_PAD_SIZE; y++)
  {
    y1 = Clip3(0, ny_1, y);
    for( x = -IMG_PAD_SIZE; x < s->size_x + IMG_PAD_SIZE; x++)
    {
      i0 =      imgY[y1][Clip3(0, nx_1, x  )];

      i1 = (
        h1[0] * imgY[y1][Clip3(0, nx_1, x-3)] +
        h1[1] * imgY[y1][Clip3(0, nx_1, x-2)] +
        h1[2] * imgY[y1][Clip3(0, nx_1, x-1)] +
        h1[3] * imgY[y1][Clip3(0, nx_1, x  )] +
        h1[4] * imgY[y1][Clip3(0, nx_1, x+1)] +
        h1[5] * imgY[y1][Clip3(0, nx_1, x+2)] +
        h1[6] * imgY[y1][Clip3(0, nx_1, x+3)] +
        h1[7] * imgY[y1][Clip3(0, nx_1, x+4)] );


      i2 = (
        h2[0] * imgY[y1][Clip3(0, nx_1, x-3)] +
        h2[1] * imgY[y1][Clip3(0, nx_1, x-2)] +
        h2[2] * imgY[y1][Clip3(0, nx_1, x-1)] +
        h2[3] * imgY[y1][Clip3(0, nx_1, x  )] +
        h2[4] * imgY[y1][Clip3(0, nx_1, x+1)] +
        h2[5] * imgY[y1][Clip3(0, nx_1, x+2)] +
        h2[6] * imgY[y1][Clip3(0, nx_1, x+3)] +
        h2[7] * imgY[y1][Clip3(0, nx_1, x+4)] );


      i3 = (
        h3[0] * imgY[y1][Clip3(0, nx_1, x-3)] +
        h3[1] * imgY[y1][Clip3(0, nx_1, x-2)] +
        h3[2] * imgY[y1][Clip3(0, nx_1, x-1)] +
        h3[3] * imgY[y1][Clip3(0, nx_1, x  )] +
        h3[4] * imgY[y1][Clip3(0, nx_1, x+1)] +
        h3[5] * imgY[y1][Clip3(0, nx_1, x+2)] +
        h3[6] * imgY[y1][Clip3(0, nx_1, x+3)] +
        h3[7] * imgY[y1][Clip3(0, nx_1, x+4)] );

      x4 = ( x + IMG_PAD_SIZE ) << 2;

      img4Y_tmp[ y + IMG_PAD_SIZE ][ x4++ ] = i0 << 8;
      img4Y_tmp[ y + IMG_PAD_SIZE ][ x4++ ] = i1;
      img4Y_tmp[ y + IMG_PAD_SIZE ][ x4++ ] = i2;
      img4Y_tmp[ y + IMG_PAD_SIZE ][ x4   ] = i3;
    }
  }

  maxy = ny_1 + IMG_PAD_SIZE_TIMES2;

  //vertical filtering
  for( x4=0; x4 < (s->size_x + IMG_PAD_SIZE_TIMES2 ) << 2; x4++ )
  {
    for( y=0; y <= maxy; y++ )
    {
      i0 = ( img4Y_tmp[y][x4] + 128 ) >> 8;

      i1 = (
        h1[0] * img4Y_tmp[max(0   ,y-3)][x4] +
        h1[1] * img4Y_tmp[max(0   ,y-2)][x4] +
        h1[2] * img4Y_tmp[max(0   ,y-1)][x4] +
        h1[3] * img4Y_tmp[y][x4]             +
        h1[4] * img4Y_tmp[min(maxy,y+1)][x4] +
        h1[5] * img4Y_tmp[min(maxy,y+2)][x4] +
        h1[6] * img4Y_tmp[min(maxy,y+3)][x4] +
        h1[7] * img4Y_tmp[min(maxy,y+4)][x4] + rounding ) >> shifting;

      i2 = (
        h2[0] * img4Y_tmp[max(0   ,y-3)][x4] +
        h2[1] * img4Y_tmp[max(0   ,y-2)][x4] +
        h2[2] * img4Y_tmp[max(0   ,y-1)][x4] +
        h2[3] * img4Y_tmp[y][x4]             +
        h2[4] * img4Y_tmp[min(maxy,y+1)][x4] +
        h2[5] * img4Y_tmp[min(maxy,y+2)][x4] +
        h2[6] * img4Y_tmp[min(maxy,y+3)][x4] +
        h2[7] * img4Y_tmp[min(maxy,y+4)][x4] + rounding ) >> shifting;

      i3 = (
        h3[0] * img4Y_tmp[max(0   ,y-3)][x4] +
        h3[1] * img4Y_tmp[max(0   ,y-2)][x4] +
        h3[2] * img4Y_tmp[max(0   ,y-1)][x4] +
        h3[3] * img4Y_tmp[y][x4]             +
        h3[4] * img4Y_tmp[min(maxy,y+1)][x4] +
        h3[5] * img4Y_tmp[min(maxy,y+2)][x4] +
        h3[6] * img4Y_tmp[min(maxy,y+3)][x4] +
        h3[7] * img4Y_tmp[min(maxy,y+4)][x4] + rounding ) >> shifting;

      y4  = y << 2;

      out4Y[ y4++ ][x4] = (pel_t) Clip3( 0, img->max_imgpel_value, i0 );
      out4Y[ y4++ ][x4] = (pel_t) Clip3( 0, img->max_imgpel_value, i1 );
      out4Y[ y4++ ][x4] = (pel_t) Clip3( 0, img->max_imgpel_value, i2 );
      out4Y[ y4   ][x4] = (pel_t) Clip3( 0, img->max_imgpel_value, i3 );
    }
  }

  // Generate 1/1th pel representation (used for integer pel MV search)
  GenerateFullPelRepresentation (out4Y, ref11, s->size_x, s->size_y);

}
#endif

/*!
************************************************************************
* \brief
*    Find SNR for all three components
************************************************************************
*/
static void find_snr (void)
{
  int i, j;
  int64 diff_y=0, diff_u=0, diff_v=0;
  int impix;
  int impix_cr;
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
  unsigned int max_pix_value_sqd = (img->BitDepthIncrease)? ((1<<input->InputBitDepth)-1)*((1<<input->InputBitDepth)-1)
    : img->max_imgpel_value * img->max_imgpel_value;
  unsigned int max_pix_value_sqd_uv = (img->BitDepthIncreaseChroma)? ((1<<input->InputBitDepth)-1)*((1<<input->InputBitDepth)-1)
    : img->max_imgpel_value_uv * img->max_imgpel_value_uv;
#else
  unsigned int max_pix_value_sqd = img->max_imgpel_value * img->max_imgpel_value;
  unsigned int max_pix_value_sqd_uv = img->max_imgpel_value_uv * img->max_imgpel_value_uv;
#endif

  //  Calculate  PSNR for Y, U and V.

  //     Luma.
  impix = input->img_height * input->img_width;
  impix_cr = input->img_height_cr * input->img_width_cr;

  if (img->fld_flag != 0)
  {

    diff_y = 0;
    for (i = 0; i < input->img_width; ++i)
    {
      for (j = 0; j < input->img_height; ++j)
      {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        diff_y += SQR_DEPTH(imgY_org[j][i], imgY_com[j][i], input->BitDepthLuma, img->BitDepthIncrease);
#else
        diff_y += img->quad[imgY_org[j][i] - imgY_com[j][i]];
#endif
      }
    }

    if (img->yuv_format != YUV400)
    {
      //     Chroma.
      diff_u = 0;
      diff_v = 0;

      for (i = 0; i < input->img_width_cr; i++)
      {
        for (j = 0; j < input->img_height_cr; j++)
        {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          diff_u += SQR_DEPTH(imgUV_org[0][j][i], imgUV_com[0][j][i], input->BitDepthChroma, img->BitDepthIncreaseChroma);
          diff_v += SQR_DEPTH(imgUV_org[1][j][i], imgUV_com[1][j][i], input->BitDepthChroma, img->BitDepthIncreaseChroma);
#else
          diff_u += img->quad[imgUV_org[0][j][i] - imgUV_com[0][j][i]];
          diff_v += img->quad[imgUV_org[1][j][i] - imgUV_com[1][j][i]];
#endif
        }
      }
    }
  }
  else
  { 
    imgY_org  = imgY_org_frm;
    imgUV_org = imgUV_org_frm;

    if(input->PicInterlace==ADAPTIVE_CODING)
    {
      enc_picture = enc_frame_picture;
    }  

    diff_y = 0;
    for (i = 0; i < input->img_width; ++i)
    {
      for (j = 0; j < input->img_height; ++j)
      {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        diff_y += SQR_DEPTH(imgY_org[j][i], enc_picture->imgY[j][i], input->BitDepthLuma, img->BitDepthIncrease);
#else
        diff_y += img->quad[imgY_org[j][i] - enc_picture->imgY[j][i]];
#endif
      }
    }

    if (img->yuv_format != YUV400)
    {
      //     Chroma.
      diff_u = 0;
      diff_v = 0;

      for (i = 0; i < input->img_width_cr; i++)
      {
        for (j = 0; j < input->img_height_cr; j++)
        {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          diff_u += SQR_DEPTH(imgUV_org[0][j][i], enc_picture->imgUV[0][j][i], input->BitDepthChroma, img->BitDepthIncreaseChroma);
          diff_v += SQR_DEPTH(imgUV_org[1][j][i], enc_picture->imgUV[1][j][i], input->BitDepthChroma, img->BitDepthIncreaseChroma);
#else
          diff_u += img->quad[imgUV_org[0][j][i] - enc_picture->imgUV[0][j][i]];
          diff_v += img->quad[imgUV_org[1][j][i] - enc_picture->imgUV[1][j][i]];
#endif
        }
      }
    }
  }
  snr->sse_y = (float)diff_y;
  snr->sse_u = (float)diff_u;
  snr->sse_v = (float)diff_v;

#if ZEROSNR
  if (diff_y == 0)
    diff_y = 1;
  if (diff_u == 0)
    diff_u = 1;
  if (diff_v == 0)
    diff_v = 1; 
#endif

  //  Collecting SNR statistics
  if (diff_y != 0)
  {
    snr->snr_y = (float) (10 * log10 (max_pix_value_sqd * (double)((double) impix    / diff_y)));         // luma snr for current frame
    if (img->yuv_format != YUV400)
    {
      snr->snr_u = (float) (10 * log10 (max_pix_value_sqd_uv * (double)((double) impix_cr / diff_u)));   // u croma snr for current frame, 1/4 of luma samples
      snr->snr_v = (float) (10 * log10 (max_pix_value_sqd_uv * (double)((double) impix_cr / diff_v)));   // v croma snr for current frame, 1/4 of luma samples
    }
    else
    {
      snr->snr_u = 0.0;
      snr->snr_v = 0.0;
    }
  }


  if (img->number == 0)
  {
    snr->snr_y1 = snr->snr_y;  // keep luma snr for first frame
    snr->snr_u1 = snr->snr_u;  // keep croma u snr for first frame
    snr->snr_v1 = snr->snr_v;  // keep croma v snr for first frame
    snr->snr_ya = snr->snr_y1;
    snr->snr_ua = snr->snr_u1;
    snr->snr_va = snr->snr_v1;
    // sse stats
    snr->msse_y = snr->sse_y;
    snr->msse_u = snr->sse_u;
    snr->msse_v = snr->sse_v;
    for (i=0; i<5; i++)
    {
      snr->snr_yt[i] = 0.0;
      snr->snr_ut[i] = 0.0;
      snr->snr_vt[i] = 0.0;
    }
  }
  else
  {
    //int total_frames = img->number + frame_ctr[B_SLICE];
    int total_frames = snr->frame_ctr - 1;

    snr->snr_ya = (float) (snr->snr_ya * total_frames + snr->snr_y) / (total_frames + 1); // average snr luma for all frames inc. first
    snr->snr_ua = (float) (snr->snr_ua * total_frames + snr->snr_u) / (total_frames + 1); // average snr u croma for all frames inc. first
    snr->snr_va = (float) (snr->snr_va * total_frames + snr->snr_v) / (total_frames + 1); // average snr v croma for all frames inc. first
    snr->msse_y = (float) (snr->msse_y * total_frames + snr->sse_y) / (total_frames + 1); // average mse luma for all frames inc. first
    snr->msse_u = (float) (snr->msse_u * total_frames + snr->sse_u) / (total_frames + 1); // average mse u croma for all frames inc. first
    snr->msse_v = (float) (snr->msse_v * total_frames + snr->sse_v) / (total_frames + 1); // average mse v croma for all frames inc. first  
  }

  snr->snr_yt[img->type] = (float) (snr->snr_yt[img->type] * (frame_ctr[img->type] - 1) + snr->snr_y) / ( frame_ctr[img->type] );  // average luma snr for img->type coded frames
  snr->snr_ut[img->type] = (float) (snr->snr_ut[img->type] * (frame_ctr[img->type] - 1) + snr->snr_u) / ( frame_ctr[img->type] );  // average chroma u snr for img->type coded frames
  snr->snr_vt[img->type] = (float) (snr->snr_vt[img->type] * (frame_ctr[img->type] - 1) + snr->snr_v) / ( frame_ctr[img->type] );  // average chroma v snr for img->type coded frames
}

/*!
************************************************************************
* \brief
*    Find distortion for all three components
************************************************************************
*/
#ifdef ADAPTIVE_QUANTIZATION
void find_distortion (void)
#else
static void find_distortion (void)
#endif
{
  int i, j;
  int64 diff_y, diff_u, diff_v;
  int impix;

  //  Calculate  PSNR for Y, U and V.

  //     Luma.
  impix = input->img_height * input->img_width;

  if (img->structure!=FRAME)
  {

    diff_y = 0;
    for (i = 0; i < input->img_width; ++i)
    {
      for (j = 0; j < input->img_height; ++j)
      {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        diff_y += SQR_DEPTH(imgY_org[j][i], imgY_com[j][i], input->BitDepthLuma, img->BitDepthIncrease);
#else
        diff_y += img->quad[imgY_org[j][i] - imgY_com[j][i]];
#endif
      }
    }

    diff_u = 0;
    diff_v = 0;

    if (img->yuv_format != YUV400)
    {
      //     Chroma.
      for (i = 0; i < input->img_width_cr; i++)
      {
        for (j = 0; j < input->img_height_cr; j++)
        {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          diff_u += SQR_DEPTH(imgUV_org[0][j][i], imgUV_com[0][j][i], input->BitDepthChroma, img->BitDepthIncreaseChroma);
          diff_v += SQR_DEPTH(imgUV_org[1][j][i], imgUV_com[1][j][i], input->BitDepthChroma, img->BitDepthIncreaseChroma);
#else
          diff_u += img->quad[imgUV_org[0][j][i] - imgUV_com[0][j][i]];
          diff_v += img->quad[imgUV_org[1][j][i] - imgUV_com[1][j][i]];
#endif
        }
      }
    }
  }
  else
  {
    imgY_org   = imgY_org_frm;
    imgUV_org = imgUV_org_frm;

    diff_y = 0;
    for (i = 0; i < input->img_width; ++i)
    {
      for (j = 0; j < input->img_height; ++j)
      {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
        diff_y += SQR_DEPTH(imgY_org[j][i], enc_picture->imgY[j][i], input->BitDepthLuma, img->BitDepthIncrease);
#else
        diff_y += img->quad[imgY_org[j][i] - enc_picture->imgY[j][i]];
#endif
      }
    }

    diff_u = 0;
    diff_v = 0;

    if (img->yuv_format != YUV400)
    {
      //     Chroma.
      for (i = 0; i < input->img_width_cr; i++)
      {
        for (j = 0; j < input->img_height_cr; j++)
        {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
          diff_u += SQR_DEPTH(imgUV_org[0][j][i], enc_picture->imgUV[0][j][i], input->BitDepthChroma, img->BitDepthIncreaseChroma);
          diff_v += SQR_DEPTH(imgUV_org[1][j][i], enc_picture->imgUV[1][j][i], input->BitDepthChroma, img->BitDepthIncreaseChroma);
#else
          diff_u += img->quad[imgUV_org[0][j][i] - enc_picture->imgUV[0][j][i]];
          diff_v += img->quad[imgUV_org[1][j][i] - enc_picture->imgUV[1][j][i]];
#endif
        }
      }
    }
  }
  // Calculate real PSNR at find_snr_avg()
  snr->snr_y = (float) diff_y;
  snr->snr_u = (float) diff_u;
  snr->snr_v = (float) diff_v;
}


/*!
************************************************************************
* \brief
*    Just a placebo
************************************************************************
*/
Boolean dummy_slice_too_big (int bits_slice)
{
  return FALSE;
}


#ifdef ADAPTIVE_FD_SD_CODING
short  enc_picture_mv[2][4][4][2];
char   enc_picture_ref_idx[2][4][4];
int64  enc_picture_ref_pic_id[2][4][4];
int64  enc_picture_ref_pic_num[6][MAX_LIST_SIZE];
double enc_picture_FilterCoefficients[15][36];

static void store_motion_vectors_FDSD_coding (void)
{
  int i,j,k,l;
#ifdef ADAPTIVE_FILTER
  for (i = 0; i < 15; i++)
  {
    for (j = 0; j < 36; j++)
    {
      enc_picture_FilterCoefficients[i][j]=enc_picture->FilterCoefficients[i][j];
    }
  }
#endif
  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < MAX_LIST_SIZE; j++)
    {
      enc_picture_ref_pic_num[i][j]=enc_picture->ref_pic_num[i][j];
    }
  }

  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 4; j++)
    {
      enc_picture_mv[LIST_0][j][i][0]=enc_picture->mv[LIST_0][j+img->block_y][i+img->block_x][0];
      enc_picture_mv[LIST_0][j][i][1]=enc_picture->mv[LIST_0][j+img->block_y][i+img->block_x][1];
      enc_picture_mv[LIST_1][j][i][0]=enc_picture->mv[LIST_1][j+img->block_y][i+img->block_x][0];
      enc_picture_mv[LIST_1][j][i][1]=enc_picture->mv[LIST_1][j+img->block_y][i+img->block_x][1];

      enc_picture_ref_idx[LIST_0][j][i]=enc_picture->ref_idx[LIST_0][j+img->block_y][i+img->block_x];
      enc_picture_ref_idx[LIST_1][j][i]=enc_picture->ref_idx[LIST_1][j+img->block_y][i+img->block_x];

      enc_picture_ref_pic_id[LIST_0][j][i]=enc_picture->ref_pic_id[LIST_0][j+img->block_y][i+img->block_x];
      enc_picture_ref_pic_id[LIST_1][j][i]=enc_picture->ref_pic_id[LIST_1][j+img->block_y][i+img->block_x];

      for (k = 0; k < img->max_num_references; k++)
      {
        for (l = 0; l < 9; l++)
        {
          if (input->BiPredMotionEstimation)
          {
            rdopt_FDSD->bipred_mv1[j][i][LIST_0][k][l][0] = img->bipred_mv1[j][i][LIST_0][k][l][0];
            rdopt_FDSD->bipred_mv1[j][i][LIST_0][k][l][1] = img->bipred_mv1[j][i][LIST_0][k][l][1];
            rdopt_FDSD->bipred_mv2[j][i][LIST_0][k][l][0] = img->bipred_mv2[j][i][LIST_0][k][l][0];
            rdopt_FDSD->bipred_mv2[j][i][LIST_0][k][l][1] = img->bipred_mv2[j][i][LIST_0][k][l][1];

            rdopt_FDSD->bipred_mv1[j][i][LIST_1][k][l][0] = img->bipred_mv1[j][i][LIST_1][k][l][0];
            rdopt_FDSD->bipred_mv1[j][i][LIST_1][k][l][1] = img->bipred_mv1[j][i][LIST_1][k][l][1];
            rdopt_FDSD->bipred_mv2[j][i][LIST_1][k][l][0] = img->bipred_mv2[j][i][LIST_1][k][l][0];
            rdopt_FDSD->bipred_mv2[j][i][LIST_1][k][l][1] = img->bipred_mv2[j][i][LIST_1][k][l][1];
          }
          rdopt_FDSD->all_mv[j][i][LIST_0][k][l][0] = img->all_mv[j][i][LIST_0][k][l][0];
          rdopt_FDSD->all_mv[j][i][LIST_0][k][l][1] = img->all_mv[j][i][LIST_0][k][l][1];

          rdopt_FDSD->all_mv[j][i][LIST_1][k][l][0] = img->all_mv[j][i][LIST_1][k][l][0];
          rdopt_FDSD->all_mv[j][i][LIST_1][k][l][1] = img->all_mv[j][i][LIST_1][k][l][1];

          rdopt_FDSD->pred_mv[j][i][LIST_0][k][l][0] = img->pred_mv[j][i][LIST_0][k][l][0];
          rdopt_FDSD->pred_mv[j][i][LIST_0][k][l][1] = img->pred_mv[j][i][LIST_0][k][l][1];

          rdopt_FDSD->pred_mv[j][i][LIST_1][k][l][0] = img->pred_mv[j][i][LIST_1][k][l][0];
          rdopt_FDSD->pred_mv[j][i][LIST_1][k][l][1] = img->pred_mv[j][i][LIST_1][k][l][1];
        }
      }
    }
  }

  for (i = 0; i < img->height/BLOCK_SIZE; i++)
  {
    for (j = 0; j < img->width/BLOCK_SIZE; j++)
    {
      rdopt_FDSD->ipredmode[i][j] = img->ipredmode[i][j];
    }
  }

  // <FTRD Compatibility with APEC
#ifdef MV_COMPETITION
  if(input->mv_competition>0)
  {
    int maxmode=0;
    int m;
    if (img->type == P_SLICE)
      maxmode = mv_comp.nb_mode_for_mvp;
    else if (img->type == B_SLICE)
      maxmode = mv_comp.nb_mode_for_mvb;

    for(j=0; j<MAXMODE; j++)
      for (i=0; i<img->num_ref_frames; i++)
        for (m=0; m<2; m++)
          for (k=0; k<4; k++)
            for (l=0; l<4; l++)
            {
              rdopt_FDSD->send_index_mv_prediction_RDDATA[j][i][m][k][l] = send_index_mv_prediction[j][i][m][k][l][maxmode];
              rdopt_FDSD->predModeMV_RDDATA[j][i][m][k][l]  =  predModeMV[j][i][m][k][l];
            }
  }
#endif  
  // FTRD>

}


static void restore_motion_vectors_FDSD_coding (void)
{
  int i,j,k,l;
#ifdef ADAPTIVE_FILTER

  for (i = 0; i < 15; i++)
  {
    for (j = 0; j < 36; j++)
    {
      enc_picture->FilterCoefficients[i][j]=enc_picture_FilterCoefficients[i][j];
    }
  }
#endif
  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < MAX_LIST_SIZE; j++)
    {
      enc_picture->ref_pic_num[i][j]=enc_picture_ref_pic_num[i][j];
    }
  }

  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 4; j++)
    {
      enc_picture->mv[LIST_0][j+img->block_y][i+img->block_x][0]=enc_picture_mv[LIST_0][j][i][0];
      enc_picture->mv[LIST_0][j+img->block_y][i+img->block_x][1]=enc_picture_mv[LIST_0][j][i][1];
      enc_picture->mv[LIST_1][j+img->block_y][i+img->block_x][0]=enc_picture_mv[LIST_1][j][i][0];
      enc_picture->mv[LIST_1][j+img->block_y][i+img->block_x][1]=enc_picture_mv[LIST_1][j][i][1];


      enc_picture->ref_idx[LIST_0][j+img->block_y][i+img->block_x]=enc_picture_ref_idx[LIST_0][j][i];
      enc_picture->ref_idx[LIST_1][j+img->block_y][i+img->block_x]=enc_picture_ref_idx[LIST_1][j][i];
      enc_picture->ref_pic_id[LIST_0][j+img->block_y][i+img->block_x]=enc_picture_ref_pic_id[LIST_0][j][i];
      enc_picture->ref_pic_id[LIST_1][j+img->block_y][i+img->block_x]=enc_picture_ref_pic_id[LIST_1][j][i];
      for (k = 0; k < img->max_num_references; k++)
      {
        for (l = 0; l < 9; l++)
        {
          //BUGFIX APEC
          if (input->BiPredMotionEstimation)
          {
            img->bipred_mv1[j][i][LIST_0][k][l][0] = rdopt_FDSD->bipred_mv1[j][i][LIST_0][k][l][0];
            img->bipred_mv1[j][i][LIST_0][k][l][1] = rdopt_FDSD->bipred_mv1[j][i][LIST_0][k][l][1];
            img->bipred_mv1[j][i][LIST_1][k][l][0] = rdopt_FDSD->bipred_mv1[j][i][LIST_1][k][l][0];
            img->bipred_mv1[j][i][LIST_1][k][l][1] = rdopt_FDSD->bipred_mv1[j][i][LIST_1][k][l][1];
            img->bipred_mv2[j][i][LIST_0][k][l][0] = rdopt_FDSD->bipred_mv2[j][i][LIST_0][k][l][0];
            img->bipred_mv2[j][i][LIST_0][k][l][1] = rdopt_FDSD->bipred_mv2[j][i][LIST_0][k][l][1];
            img->bipred_mv2[j][i][LIST_1][k][l][0] = rdopt_FDSD->bipred_mv2[j][i][LIST_1][k][l][0];
            img->bipred_mv2[j][i][LIST_1][k][l][1] = rdopt_FDSD->bipred_mv2[j][i][LIST_1][k][l][1];
          }
          //END BUGFIX APEC

          img->all_mv[j][i][LIST_0][k][l][0] = rdopt_FDSD->all_mv[j][i][LIST_0][k][l][0];
          img->all_mv[j][i][LIST_0][k][l][1] = rdopt_FDSD->all_mv[j][i][LIST_0][k][l][1];

          img->all_mv[j][i][LIST_1][k][l][0] = rdopt_FDSD->all_mv[j][i][LIST_1][k][l][0];
          img->all_mv[j][i][LIST_1][k][l][1] = rdopt_FDSD->all_mv[j][i][LIST_1][k][l][1];

          img->pred_mv[j][i][LIST_0][k][l][0] = rdopt_FDSD->pred_mv[j][i][LIST_0][k][l][0];
          img->pred_mv[j][i][LIST_0][k][l][1] = rdopt_FDSD->pred_mv[j][i][LIST_0][k][l][1];

          img->pred_mv[j][i][LIST_1][k][l][0] = rdopt_FDSD->pred_mv[j][i][LIST_1][k][l][0];
          img->pred_mv[j][i][LIST_1][k][l][1] = rdopt_FDSD->pred_mv[j][i][LIST_1][k][l][1];
        }
      }
    }
  }
  for (i = 0; i < img->height/BLOCK_SIZE; i++)
  {
    for (j = 0; j < img->width/BLOCK_SIZE; j++)
    {
      img->ipredmode[i][j]=rdopt_FDSD->ipredmode[i][j];
    }
  }

  // <FTRD Compatibility with APEC
#ifdef MV_COMPETITION
  if(input->mv_competition>0)
  {
    int maxmode=0;
    int m;
    if (img->type == P_SLICE)
      maxmode = mv_comp.nb_mode_for_mvp;
    else if (img->type == B_SLICE)
      maxmode = mv_comp.nb_mode_for_mvb;

    for(j=0; j<MAXMODE; j++)
      for (i=0; i<img->num_ref_frames; i++)
        for (m=0; m<2; m++)
          for (k=0; k<4; k++)
            for (l=0; l<4; l++)
            {
              send_index_mv_prediction[j][i][m][k][l][maxmode] = rdopt_FDSD->send_index_mv_prediction_RDDATA[j][i][m][k][l];
              predModeMV[j][i][m][k][l] = rdopt_FDSD->predModeMV_RDDATA[j][i][m][k][l] ;
            }
  }
#endif
  // FTRD>

}



void store_rdopt_data_for_FDSD_coding ()
{
  int mb_nr = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];
  int i, j, k;

  int bframe = (img->type == B_SLICE);
  int block_y;

  //===== reconstruction values =====
  for (j = 0; j < MB_BLOCK_SIZE; j++)
    memcpy(rdopt_FDSD->rec_mbY[j],&enc_picture->imgY[img->pix_y + j][img->pix_x], MB_BLOCK_SIZE * sizeof(imgpel));


  if (img->yuv_format != YUV400)
  {
    for (j = 0; j < img->mb_cr_size_y; j++)
    {
      memcpy(rdopt_FDSD->rec_mbU[j],&enc_picture->imgUV[0][img->pix_c_y + j][img->pix_c_x], img->mb_cr_size_x * sizeof(imgpel));
      memcpy(rdopt_FDSD->rec_mbV[j],&enc_picture->imgUV[1][img->pix_c_y + j][img->pix_c_x], img->mb_cr_size_x * sizeof(imgpel));
    }
  }


  for (i = 0; i < 4+img->num_blk8x8_uv; i++)
  {
    for (j = 0; j < 4; j++)
      for (k = 0; k < 2; k++)
        memcpy(rdopt_FDSD->cofAC[i][j][k],img->cofAC[i][j][k], 65 * sizeof(int));
  }

  for (i = 0; i < 3; i++)
  {
    for (k = 0; k < 2; k++)
      memcpy(rdopt_FDSD->cofDC[i][k],img->cofDC[i][k], 18 * sizeof(int));
  }

  rdopt_FDSD->mb_type      = currMB->mb_type;
  rdopt_FDSD->bi_pred_me   = currMB->bi_pred_me;

  memcpy(rdopt_FDSD->b8mode,currMB->b8mode, 4 * sizeof(int));
  memcpy(rdopt_FDSD->b8pdir,currMB->b8pdir, 4 * sizeof(int));

  memcpy(rdopt_FDSD->intra_pred_modes   ,currMB->intra_pred_modes   , MB_BLOCK_PARTITIONS * sizeof(char));
  memcpy(rdopt_FDSD->intra_pred_modes8x8,currMB->intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));

  rdopt_FDSD->cbp          = currMB->cbp;
  rdopt_FDSD->cbp_blk      = currMB->cbp_blk;
  rdopt_FDSD->mode         = currMB->mb_type;
#ifdef ADAPTIVE_QUANTIZATION
  rdopt_FDSD->mb_iaqms_idx         = currMB->mb_iaqms_idx;
#endif

  for (j = 0; j < BLOCK_MULTIPLE; j++)
  {
    block_y = img->block_y + j;
    memcpy( rdopt_FDSD->refar[LIST_0][j],&enc_picture->ref_idx[LIST_0][block_y][img->block_x], BLOCK_MULTIPLE * sizeof(char));
  }

  rdopt_FDSD->i16offset    = img->i16offset;
  rdopt_FDSD->c_ipred_mode = currMB->c_ipred_mode;
  rdopt_FDSD->luma_transform_size_8x8_flag=currMB->luma_transform_size_8x8_flag;
  rdopt_FDSD->qp           = currMB->qp;
  rdopt_FDSD->prev_qp      = currMB->prev_qp;
  rdopt_FDSD->prev_delta_qp= currMB->prev_delta_qp;
  rdopt_FDSD->prev_cbp     = currMB->prev_cbp;
  rdopt_FDSD->delta_qp     = currMB->delta_qp;


  if (bframe)
  {
    for (j = 0; j < BLOCK_MULTIPLE; j++)
    {
      block_y = img->block_y + j;
      memcpy( rdopt_FDSD->refar[LIST_1][j],&enc_picture->ref_idx[LIST_1][block_y][img->block_x], BLOCK_MULTIPLE * sizeof(char));
    }
  }

  rdopt_FDSD->best_SD_Coding_on_off = currMB->SD_Coding_on_off;
  memcpy(rdopt_FDSD->best_quantizer_indices,currMB->quantizer_indices,256*sizeof(int));
  memcpy(rdopt_FDSD->best_SD_or_FD,currMB->SD_or_FD,4*sizeof(int));
  rdopt_FDSD->best_SD_or_FD_t8x8=currMB->SD_or_FD_t8x8;

  store_motion_vectors_FDSD_coding ();
}                           







void restore_rdopt_data_for_FDSD_coding ()
{
  int mb_nr = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];
  int i, j, k;
  int mode;
  int bframe = (img->type == B_SLICE);
  int block_y;
  int list_offset = currMB->list_offset;

  mode                     = rdopt_FDSD->mode;
  currMB->mb_type          = rdopt_FDSD->mb_type;
  currMB->cbp              = rdopt_FDSD->cbp;
  currMB->cbp_blk          = rdopt_FDSD->cbp_blk;
  currMB->bi_pred_me       = rdopt_FDSD->bi_pred_me;
  img->i16offset           = rdopt_FDSD->i16offset;
  currMB->prev_qp          = rdopt_FDSD->prev_qp;
  currMB->prev_delta_qp    = rdopt_FDSD->prev_delta_qp;
  currMB->prev_cbp         = rdopt_FDSD->prev_cbp;
  currMB->delta_qp         = rdopt_FDSD->delta_qp;
  currMB->qp               = rdopt_FDSD->qp;
  currMB->c_ipred_mode     = rdopt_FDSD->c_ipred_mode;
#ifdef ADAPTIVE_QUANTIZATION
  currMB->mb_iaqms_idx     = rdopt_FDSD->mb_iaqms_idx;
#endif
  set_chroma_qp(currMB);

  for (i = 0; i < 4+img->num_blk8x8_uv; i++)
  {
    for (j = 0; j < 4; j++)
      for (k = 0; k < 2; k++)
        memcpy(img->cofAC[i][j][k],rdopt_FDSD->cofAC[i][j][k], 65 * sizeof(int));
  }

  for (i = 0; i < 3; i++)
  {
    for (k = 0; k < 2; k++)
      memcpy(img->cofDC[i][k],rdopt_FDSD->cofDC[i][k], 18 * sizeof(int));
  }

  for (j = 0; j < BLOCK_MULTIPLE; j++)
  {
    block_y = img->block_y + j;
    memcpy(&enc_picture->ref_idx[LIST_0][block_y][img->block_x], rdopt_FDSD->refar[LIST_0][j], BLOCK_MULTIPLE * sizeof(char));
    for (i = 0; i < BLOCK_MULTIPLE; i++)
      enc_picture->ref_pic_id [LIST_0][block_y][img->block_x + i] = 
      enc_picture->ref_pic_num[LIST_0 + list_offset][(short)enc_picture->ref_idx[LIST_0][block_y][img->block_x+i]];
  } 
  if (bframe)
  {
    for (j = 0; j < BLOCK_MULTIPLE; j++)
    {
      block_y = img->block_y + j;
      memcpy(&enc_picture->ref_idx[LIST_1][block_y][img->block_x], rdopt_FDSD->refar[LIST_1][j], BLOCK_MULTIPLE * sizeof(char));
      for (i = 0; i < BLOCK_MULTIPLE; i++)
        enc_picture->ref_pic_id [LIST_1][block_y][img->block_x + i] = 
        enc_picture->ref_pic_num[LIST_1 + list_offset][(short)enc_picture->ref_idx[LIST_1][block_y][img->block_x+i]];
    }
  }

  //===== reconstruction values =====
  for (j = 0; j < MB_BLOCK_SIZE; j++)
    memcpy(&enc_picture->imgY[img->pix_y + j][img->pix_x],rdopt_FDSD->rec_mbY[j], MB_BLOCK_SIZE * sizeof(imgpel));


  currMB->SD_Coding_on_off = rdopt_FDSD->best_SD_Coding_on_off;
  memcpy(currMB->quantizer_indices,rdopt_FDSD->best_quantizer_indices,256*sizeof(int));
  memcpy(currMB->SD_or_FD,rdopt_FDSD->best_SD_or_FD,4*sizeof(int));
  currMB->SD_or_FD_t8x8=rdopt_FDSD->best_SD_or_FD_t8x8;


  if (img->yuv_format != YUV400)
  {
    for (j = 0; j < img->mb_cr_size_y; j++)
    {
      memcpy(&enc_picture->imgUV[0][img->pix_c_y + j][img->pix_c_x],rdopt_FDSD->rec_mbU[j], img->mb_cr_size_x * sizeof(imgpel));
      memcpy(&enc_picture->imgUV[1][img->pix_c_y + j][img->pix_c_x],rdopt_FDSD->rec_mbV[j], img->mb_cr_size_x * sizeof(imgpel));
    }
  }

  memcpy(currMB->b8mode,rdopt_FDSD->b8mode, 4 * sizeof(int));
  memcpy(currMB->b8pdir,rdopt_FDSD->b8pdir, 4 * sizeof(int));

  currMB->luma_transform_size_8x8_flag=rdopt_FDSD->luma_transform_size_8x8_flag;

  memcpy(currMB->intra_pred_modes   ,rdopt_FDSD->intra_pred_modes   , MB_BLOCK_PARTITIONS * sizeof(char));
  memcpy(currMB->intra_pred_modes8x8,rdopt_FDSD->intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));

  restore_motion_vectors_FDSD_coding ();
}                  



void store_rdopt_data_for_FDSD_coding_interlace ()
{
  int i, j, k, l;

  rdopt_FDSD_interlace->min_rdcost   = rdopt->min_rdcost;
  rdopt_FDSD_interlace->mode         = rdopt->mode;
  rdopt_FDSD_interlace->mb_type      = rdopt->mb_type;
  rdopt_FDSD_interlace->cbp          = rdopt->cbp;
  rdopt_FDSD_interlace->cbp_blk      = rdopt->cbp_blk; 
  rdopt_FDSD_interlace->bi_pred_me   = rdopt->bi_pred_me;
  rdopt_FDSD_interlace->i16offset    = rdopt->i16offset;
  rdopt_FDSD_interlace->prev_qp      = rdopt->prev_qp;
  rdopt_FDSD_interlace->prev_delta_qp= rdopt->prev_delta_qp;
  rdopt_FDSD_interlace->prev_cbp     = rdopt->prev_cbp;
  rdopt_FDSD_interlace->delta_qp     = rdopt->delta_qp;
  rdopt_FDSD_interlace->qp           = rdopt->qp;
  rdopt_FDSD_interlace->c_ipred_mode = rdopt->c_ipred_mode;
  rdopt_FDSD_interlace->luma_transform_size_8x8_flag = rdopt->luma_transform_size_8x8_flag;

  for (i = 0; i < 4+img->num_blk8x8_uv; i++)
  {
    for (j = 0; j < 4; j++)
      for (k = 0; k < 2; k++)
        memcpy(rdopt_FDSD_interlace->cofAC[i][j][k],rdopt->cofAC[i][j][k], 65 * sizeof(int));
  }

  for (i = 0; i < 3; i++)
  {
    for (k = 0; k < 2; k++)
      memcpy(rdopt_FDSD_interlace->cofDC[i][k],rdopt->cofDC[i][k], 18 * sizeof(int));
  }

  for (j = 0; j < BLOCK_MULTIPLE; j++)
  {
    memcpy(rdopt_FDSD_interlace->refar[LIST_0][j], rdopt->refar[LIST_0][j], BLOCK_MULTIPLE * sizeof(char));
    memcpy(rdopt_FDSD_interlace->refar[LIST_1][j], rdopt->refar[LIST_1][j], BLOCK_MULTIPLE * sizeof(char));
  } 

  //===== reconstruction values =====
  for (j = 0; j < MB_BLOCK_SIZE; j++)
    memcpy(rdopt_FDSD_interlace->rec_mbY[j],rdopt->rec_mbY[j], MB_BLOCK_SIZE * sizeof(imgpel));

  rdopt_FDSD_interlace->best_SD_Coding_on_off=rdopt->best_SD_Coding_on_off;

  memcpy(rdopt_FDSD_interlace->best_quantizer_indices,rdopt->best_quantizer_indices,256*sizeof(int));
  memcpy(rdopt_FDSD_interlace->best_SD_or_FD,rdopt->best_SD_or_FD,4*sizeof(int));

  rdopt_FDSD_interlace->best_SD_or_FD_t8x8=rdopt->best_SD_or_FD_t8x8;

  if (img->yuv_format != YUV400)
  {
    for (j = 0; j < img->mb_cr_size_y; j++)
    {
      memcpy(rdopt_FDSD_interlace->rec_mbU[j],rdopt->rec_mbU[j], img->mb_cr_size_x * sizeof(imgpel));
      memcpy(rdopt_FDSD_interlace->rec_mbV[j],rdopt->rec_mbV[j], img->mb_cr_size_x * sizeof(imgpel));
    }
  }

  memcpy(rdopt_FDSD_interlace->b8mode,rdopt->b8mode, 4 * sizeof(int));
  memcpy(rdopt_FDSD_interlace->b8pdir,rdopt->b8pdir, 4 * sizeof(int));


  memcpy(rdopt_FDSD_interlace->intra_pred_modes,rdopt->intra_pred_modes, MB_BLOCK_PARTITIONS * sizeof(char));
  memcpy(rdopt_FDSD_interlace->intra_pred_modes8x8,rdopt->intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));

  for (i = 0; i < img->height/BLOCK_SIZE; i++)
  {
    for (j = 0; j < img->width/BLOCK_SIZE; j++)
    {
      rdopt_FDSD_interlace->ipredmode[i][j]=rdopt->ipredmode[i][j];
    }
  }

  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < img->max_num_references; k++)
      {
        for (l = 0; l < 9; l++)
        {
          rdopt_FDSD_interlace->all_mv[j][i][LIST_0][k][l][0] = rdopt->all_mv[j][i][LIST_0][k][l][0];
          rdopt_FDSD_interlace->all_mv[j][i][LIST_0][k][l][1] = rdopt->all_mv[j][i][LIST_0][k][l][1];

          rdopt_FDSD_interlace->all_mv[j][i][LIST_1][k][l][0] = rdopt->all_mv[j][i][LIST_1][k][l][0];
          rdopt_FDSD_interlace->all_mv[j][i][LIST_1][k][l][1] = rdopt->all_mv[j][i][LIST_1][k][l][1];

          rdopt_FDSD_interlace->pred_mv[j][i][LIST_0][k][l][0] = rdopt->pred_mv[j][i][LIST_0][k][l][0];
          rdopt_FDSD_interlace->pred_mv[j][i][LIST_0][k][l][1] = rdopt->pred_mv[j][i][LIST_0][k][l][1];

          rdopt_FDSD_interlace->pred_mv[j][i][LIST_1][k][l][0] = rdopt->pred_mv[j][i][LIST_1][k][l][0];
          rdopt_FDSD_interlace->pred_mv[j][i][LIST_1][k][l][1] = rdopt->pred_mv[j][i][LIST_1][k][l][1];
        }
      }
    }
  }
}


void restore_rdopt_data_for_FDSD_coding_interlace ()
{
  int i, j, k,l;

  rdopt->min_rdcost   = rdopt_FDSD_interlace->min_rdcost;
  rdopt->mode         = rdopt_FDSD_interlace->mode;
  rdopt->mb_type      = rdopt_FDSD_interlace->mb_type;
  rdopt->cbp          = rdopt_FDSD_interlace->cbp;
  rdopt->cbp_blk      = rdopt_FDSD_interlace->cbp_blk; 
  rdopt->bi_pred_me   = rdopt_FDSD_interlace->bi_pred_me;
  rdopt->i16offset    = rdopt_FDSD_interlace->i16offset;
  rdopt->prev_qp      = rdopt_FDSD_interlace->prev_qp;
  rdopt->prev_delta_qp= rdopt_FDSD_interlace->prev_delta_qp;
  rdopt->prev_cbp     = rdopt_FDSD_interlace->prev_cbp;
  rdopt->delta_qp     = rdopt_FDSD_interlace->delta_qp;
  rdopt->qp           = rdopt_FDSD_interlace->qp;

  rdopt->c_ipred_mode = rdopt_FDSD_interlace->c_ipred_mode;

  for (i = 0; i < 4+img->num_blk8x8_uv; i++)
  {
    for (j = 0; j < 4; j++)
      for (k = 0; k < 2; k++)
        memcpy(rdopt->cofAC[i][j][k],rdopt_FDSD_interlace->cofAC[i][j][k], 65 * sizeof(int));
  }

  for (i = 0; i < 3; i++)
  {
    for (k = 0; k < 2; k++)
      memcpy(rdopt->cofDC[i][k],rdopt_FDSD_interlace->cofDC[i][k], 18 * sizeof(int));
  }

  for (j = 0; j < BLOCK_MULTIPLE; j++)
  {
    memcpy(rdopt->refar[LIST_0][j], rdopt_FDSD_interlace->refar[LIST_0][j], BLOCK_MULTIPLE * sizeof(char));
    memcpy(rdopt->refar[LIST_1][j], rdopt_FDSD_interlace->refar[LIST_1][j], BLOCK_MULTIPLE * sizeof(char));
  } 

  //===== reconstruction values =====
  for (j = 0; j < MB_BLOCK_SIZE; j++)
    memcpy(rdopt->rec_mbY[j],rdopt_FDSD_interlace->rec_mbY[j], MB_BLOCK_SIZE * sizeof(imgpel));

  rdopt->best_SD_Coding_on_off=rdopt_FDSD_interlace->best_SD_Coding_on_off;

  memcpy(rdopt->best_quantizer_indices,rdopt_FDSD_interlace->best_quantizer_indices,256*sizeof(int));
  memcpy(rdopt->best_SD_or_FD,rdopt_FDSD_interlace->best_SD_or_FD,4*sizeof(int));

  rdopt->best_SD_or_FD_t8x8=rdopt_FDSD_interlace->best_SD_or_FD_t8x8;

  if (img->yuv_format != YUV400)
  {
    for (j = 0; j < img->mb_cr_size_y; j++)
    {
      memcpy(rdopt->rec_mbU[j],rdopt_FDSD_interlace->rec_mbU[j], img->mb_cr_size_x * sizeof(imgpel));
      memcpy(rdopt->rec_mbV[j],rdopt_FDSD_interlace->rec_mbV[j], img->mb_cr_size_x * sizeof(imgpel));
    }
  }

  memcpy(rdopt->b8mode,rdopt_FDSD_interlace->b8mode, 4 * sizeof(int));
  memcpy(rdopt->b8pdir,rdopt_FDSD_interlace->b8pdir, 4 * sizeof(int));

  rdopt->luma_transform_size_8x8_flag = rdopt_FDSD_interlace->luma_transform_size_8x8_flag;

  memcpy(rdopt->intra_pred_modes,rdopt_FDSD_interlace->intra_pred_modes, MB_BLOCK_PARTITIONS * sizeof(char));
  memcpy(rdopt->intra_pred_modes8x8,rdopt_FDSD_interlace->intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));

  for (i = 0; i < img->height/BLOCK_SIZE; i++)
  {
    for (j = 0; j < img->width/BLOCK_SIZE; j++)
    {
      rdopt->ipredmode[i][j]=rdopt_FDSD_interlace->ipredmode[i][j];
    }
  }

  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < img->max_num_references; k++)
      {
        for (l = 0; l < 9; l++)
        {
          rdopt->all_mv[j][i][LIST_0][k][l][0] = rdopt_FDSD_interlace->all_mv[j][i][LIST_0][k][l][0];
          rdopt->all_mv[j][i][LIST_0][k][l][1] = rdopt_FDSD_interlace->all_mv[j][i][LIST_0][k][l][1];

          rdopt->all_mv[j][i][LIST_1][k][l][0] = rdopt_FDSD_interlace->all_mv[j][i][LIST_1][k][l][0];
          rdopt->all_mv[j][i][LIST_1][k][l][1] = rdopt_FDSD_interlace->all_mv[j][i][LIST_1][k][l][1];

          rdopt->pred_mv[j][i][LIST_0][k][l][0] = rdopt_FDSD_interlace->pred_mv[j][i][LIST_0][k][l][0];
          rdopt->pred_mv[j][i][LIST_0][k][l][1] = rdopt_FDSD_interlace->pred_mv[j][i][LIST_0][k][l][1];

          rdopt->pred_mv[j][i][LIST_1][k][l][0] = rdopt_FDSD_interlace->pred_mv[j][i][LIST_1][k][l][0];
          rdopt->pred_mv[j][i][LIST_1][k][l][1] = rdopt_FDSD_interlace->pred_mv[j][i][LIST_1][k][l][1];
        }
      }
    }
  }
}
#endif






/*! 
***************************************************************************
// For MB level field/frame coding
***************************************************************************
*/
void copy_rdopt_data (int bot_block)
{
  int mb_nr = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];
  int i, j, k;

  int bframe = (img->type == B_SLICE);
  int mode;
  int b8mode, b8pdir;
  int block_y;

  int list_offset = currMB->list_offset;

  mode                = rdopt->mode;
  currMB->mb_type     = rdopt->mb_type;   // copy mb_type 
  currMB->cbp         = rdopt->cbp;   // copy cbp
  currMB->cbp_blk     = rdopt->cbp_blk;   // copy cbp_blk
  currMB->bi_pred_me  = rdopt->bi_pred_me;   // copy biprediction
  img->i16offset      = rdopt->i16offset;

  currMB->prev_qp=rdopt->prev_qp;
  currMB->prev_delta_qp=rdopt->prev_delta_qp;
  currMB->prev_cbp=rdopt->prev_cbp;
  currMB->delta_qp=rdopt->delta_qp;
  currMB->qp = rdopt->qp;
  set_chroma_qp(currMB);
#ifdef ADAPTIVE_QUANTIZATION
  currMB->mb_iaqms_idx = rdopt->mb_iaqms_idx;
#endif
  currMB->c_ipred_mode = rdopt->c_ipred_mode;

  for (i = 0; i < 4+img->num_blk8x8_uv; i++)
  {
    for (j = 0; j < 4; j++)
      for (k = 0; k < 2; k++)
        memcpy(img->cofAC[i][j][k],rdopt->cofAC[i][j][k], 65 * sizeof(int));
  }

  for (i = 0; i < 3; i++)
  {
    for (k = 0; k < 2; k++)
      memcpy(img->cofDC[i][k],rdopt->cofDC[i][k], 18 * sizeof(int));
  }

  for (j = 0; j < BLOCK_MULTIPLE; j++)
  {
    block_y = img->block_y + j;
    memcpy(&enc_picture->ref_idx[LIST_0][block_y][img->block_x], rdopt->refar[LIST_0][j], BLOCK_MULTIPLE * sizeof(char));
    for (i = 0; i < BLOCK_MULTIPLE; i++)
      enc_picture->ref_pic_id [LIST_0][block_y][img->block_x + i] = 
      enc_picture->ref_pic_num[LIST_0 + list_offset][(short)enc_picture->ref_idx[LIST_0][block_y][img->block_x+i]];
  } 
  if (bframe)
  {
    for (j = 0; j < BLOCK_MULTIPLE; j++)
    {
      block_y = img->block_y + j;
      memcpy(&enc_picture->ref_idx[LIST_1][block_y][img->block_x], rdopt->refar[LIST_1][j], BLOCK_MULTIPLE * sizeof(char));
      for (i = 0; i < BLOCK_MULTIPLE; i++)
        enc_picture->ref_pic_id [LIST_1][block_y][img->block_x + i] = 
        enc_picture->ref_pic_num[LIST_1 + list_offset][(short)enc_picture->ref_idx[LIST_1][block_y][img->block_x+i]];
    }
  }

  //===== reconstruction values =====
  for (j = 0; j < MB_BLOCK_SIZE; j++)
    memcpy(&enc_picture->imgY[img->pix_y + j][img->pix_x],rdopt->rec_mbY[j], MB_BLOCK_SIZE * sizeof(imgpel));


#ifdef ADAPTIVE_FD_SD_CODING
  currMB->SD_Coding_on_off=rdopt->best_SD_Coding_on_off;
  memcpy(currMB->quantizer_indices,rdopt->best_quantizer_indices,256*sizeof(int));
  memcpy(currMB->SD_or_FD,rdopt->best_SD_or_FD,4*sizeof(int));
  currMB->SD_or_FD_t8x8=rdopt->best_SD_or_FD_t8x8;
#endif



  if (img->yuv_format != YUV400)
  {
    for (j = 0; j < img->mb_cr_size_y; j++)
    {
      memcpy(&enc_picture->imgUV[0][img->pix_c_y + j][img->pix_c_x],rdopt->rec_mbU[j], img->mb_cr_size_x * sizeof(imgpel));
      memcpy(&enc_picture->imgUV[1][img->pix_c_y + j][img->pix_c_x],rdopt->rec_mbV[j], img->mb_cr_size_x * sizeof(imgpel));
    }
  }

  memcpy(currMB->b8mode,rdopt->b8mode, 4 * sizeof(int));
  memcpy(currMB->b8pdir,rdopt->b8pdir, 4 * sizeof(int));

  currMB->luma_transform_size_8x8_flag = rdopt->luma_transform_size_8x8_flag;

  //==== intra prediction modes ====
  if (mode == P8x8)
  {
    memcpy(currMB->intra_pred_modes,rdopt->intra_pred_modes, MB_BLOCK_PARTITIONS * sizeof(char));
    for (j = img->block_y; j < img->block_y + BLOCK_MULTIPLE; j++)
      memcpy(&img->ipredmode[j][img->block_x],&rdopt->ipredmode[j][img->block_x], BLOCK_MULTIPLE * sizeof(char));
  }
  else if (mode != I4MB && mode != I8MB)
  {
    memset(currMB->intra_pred_modes,DC_PRED, MB_BLOCK_PARTITIONS * sizeof(char));
    for (j = img->block_y; j < img->block_y + BLOCK_MULTIPLE; j++)
      memset(&img->ipredmode[j][img->block_x],DC_PRED, BLOCK_MULTIPLE * sizeof(char));
  }
  else if (mode == I4MB || mode == I8MB)
  {
    memcpy(currMB->intra_pred_modes,rdopt->intra_pred_modes, MB_BLOCK_PARTITIONS * sizeof(char));
    memcpy(currMB->intra_pred_modes8x8,rdopt->intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));
    for (j = img->block_y; j < img->block_y + BLOCK_MULTIPLE; j++)
    {
      memcpy(&img->ipredmode[j][img->block_x],&rdopt->ipredmode[j][img->block_x], BLOCK_MULTIPLE * sizeof(char));
    }
  }

  if (img->MbaffFrameFlag)
  {
    // motion vectors
    copy_motion_vectors_MB ();

    if (!IS_INTRA(currMB))
    {
      for (j = 0; j < 4; j++)
        for (i = 0; i < 4; i++)
        {
          b8mode = currMB->b8mode[i/2+2*(j/2)];
          b8pdir = currMB->b8pdir[i/2+2*(j/2)];

          if (b8pdir!=1)
          {
            enc_picture->mv[LIST_0][j+img->block_y][i+img->block_x][0] = rdopt->all_mv[j][i][LIST_0][(short)rdopt->refar[LIST_0][j][i]][b8mode][0];
            enc_picture->mv[LIST_0][j+img->block_y][i+img->block_x][1] = rdopt->all_mv[j][i][LIST_0][(short)rdopt->refar[LIST_0][j][i]][b8mode][1];
          }
          else
          {
            enc_picture->mv[LIST_0][j+img->block_y][i+img->block_x][0] = 0;
            enc_picture->mv[LIST_0][j+img->block_y][i+img->block_x][1] = 0;
          }
          if (bframe)
          {
            if (b8pdir!=0)
            {
              enc_picture->mv[LIST_1][j+img->block_y][i+img->block_x][0] = rdopt->all_mv[j][i][LIST_1][(short)rdopt->refar[LIST_1][j][i]][b8mode][0];
              enc_picture->mv[LIST_1][j+img->block_y][i+img->block_x][1] = rdopt->all_mv[j][i][LIST_1][(short)rdopt->refar[LIST_1][j][i]][b8mode][1];
            }
            else
            {
              enc_picture->mv[LIST_1][j+img->block_y][i+img->block_x][0] = 0;
              enc_picture->mv[LIST_1][j+img->block_y][i+img->block_x][1] = 0;
            }
          }
        }
    }
    else
    {
      for (j = 0; j < 4; j++)
        memset(enc_picture->mv[LIST_0][j+img->block_y][img->block_x], 0, 2 * BLOCK_MULTIPLE * sizeof(short));
      if (bframe)
      {
        for (j = 0; j < 4; j++)
          memset(enc_picture->mv[LIST_1][j+img->block_y][img->block_x], 0, 2 * BLOCK_MULTIPLE * sizeof(short));
      }
    }
  }  
}                             // end of copy_rdopt_data

static void copy_motion_vectors_MB (void)
{
  int i,j,k,l;

  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < img->max_num_references; k++)
      {
        for (l = 0; l < 9; l++)
        {
          img->all_mv[j][i][LIST_0][k][l][0] = rdopt->all_mv[j][i][LIST_0][k][l][0];
          img->all_mv[j][i][LIST_0][k][l][1] = rdopt->all_mv[j][i][LIST_0][k][l][1];

          img->all_mv[j][i][LIST_1][k][l][0] = rdopt->all_mv[j][i][LIST_1][k][l][0];
          img->all_mv[j][i][LIST_1][k][l][1] = rdopt->all_mv[j][i][LIST_1][k][l][1];

          img->pred_mv[j][i][LIST_0][k][l][0] = rdopt->pred_mv[j][i][LIST_0][k][l][0];
          img->pred_mv[j][i][LIST_0][k][l][1] = rdopt->pred_mv[j][i][LIST_0][k][l][1];

          img->pred_mv[j][i][LIST_1][k][l][0] = rdopt->pred_mv[j][i][LIST_1][k][l][0];
          img->pred_mv[j][i][LIST_1][k][l][1] = rdopt->pred_mv[j][i][LIST_1][k][l][1];
        }
      }
    }
  }
}


static void ReportNALNonVLCBits(int tmp_time, int me_time)
{

  //! Need to add type (i.e. SPS, PPS, SEI etc).
  if (input->Verbose != 0)
    printf ("%04d(NVB)%8d \n", frame_no, stats->bit_ctr_parametersets_n);
}
static void ReportFirstframe(int tmp_time,int me_time)
{
  //Rate control
  int bits;

  if (input->Verbose == 1)
  {
#ifdef ADAPTIVE_FILTER
    printf ("%04d(IDR)%8d   %2d %7.3f %7.3f %7.3f %9d %7d    %3s    %d  %2d\n",
      frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n), 
      img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
      img->fld_flag ? "FLD" : "FRM", img->nal_reference_idc, img->AdaptiveFilterFlag);
#else
    printf ("%04d(IDR)%8d   %2d %7.3f %7.3f %7.3f %9d %7d    %3s    %d\n",
      frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n), 
      img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
      img->fld_flag ? "FLD" : "FRM", img->nal_reference_idc); 
#endif
  }
  else if (input->Verbose == 2)
  {
    printf ("%04d(IDR)%8d %1d %2d %7.3f %7.3f %7.3f %9d %7d    %3s %5d   %2d %2d  %d   %d\n",
      frame_no,(int) (stats->bit_ctr - stats->bit_ctr_n),0,
      img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
      img->fld_flag ? "FLD" : "FRM", intras, img->num_ref_idx_l0_active, img->num_ref_idx_l1_active,img->rd_pass, img->nal_reference_idc);
  }
  //Rate control
  if(input->RCEnable)
  {
    if((!input->PicInterlace)&&(!input->MbInterlace))
      bits = (int) (stats->bit_ctr - stats->bit_ctr_n); // used for rate control update 
    else
    {
      bits = (int)(stats->bit_ctr - Iprev_bits); // used for rate control update 
      Iprev_bits = stats->bit_ctr;
    }
  }

  stats->bit_ctr_I = stats->bit_ctr;
  stats->bit_ctr = 0;

}


static void ReportIntra(int tmp_time, int me_time)
{
  if (input->Verbose == 1)
  {
#ifdef ADAPTIVE_FILTER
    if (img->currentPicture->idr_flag == 1)
      printf ("%04d(IDR)%8d   %2d %7.3f %7.3f %7.3f %9d %7d    %3s    %d  %2d\n",
      frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n), 
      img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
      img->fld_flag ? "FLD" : "FRM", img->nal_reference_idc, img->AdaptiveFilterFlag);
    else
      printf ("%04d(I)  %8d   %2d %7.3f %7.3f %7.3f %9d %7d    %3s    %d   %2d\n",
      frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n), 
      img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
      img->fld_flag ? "FLD" : "FRM", img->nal_reference_idc, img->AdaptiveFilterFlag);
#else
    if (img->currentPicture->idr_flag == 1)
      printf ("%04d(IDR)%8d   %2d %7.3f %7.3f %7.3f %9d %7d    %3s    %d\n",
      frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n), 
      img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
      img->fld_flag ? "FLD" : "FRM", img->nal_reference_idc);
    else
      printf ("%04d(I)  %8d   %2d %7.3f %7.3f %7.3f %9d %7d    %3s    %d\n",
      frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n), 
      img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
      img->fld_flag ? "FLD" : "FRM", img->nal_reference_idc);
#endif
  }
  else if (input->Verbose == 2)
  {
    if (img->currentPicture->idr_flag == 1)
      printf ("%04d(IDR)%8d %1d %2d %7.3f %7.3f %7.3f %9d %7d    %3s %5d   %2d %2d  %d   %d\n",
      frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n), 0,
      img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
      img->fld_flag ? "FLD" : "FRM", intras, img->num_ref_idx_l0_active, img->num_ref_idx_l1_active,img->rd_pass, img->nal_reference_idc); 
    else
      printf ("%04d(I)  %8d %1d %2d %7.3f %7.3f %7.3f %9d %7d    %3s %5d   %2d %2d  %d   %d\n",
      frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n), 0,
      img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
      img->fld_flag ? "FLD" : "FRM", intras, img->num_ref_idx_l0_active, img->num_ref_idx_l1_active,img->rd_pass, img->nal_reference_idc);
  }
}

static void ReportSP(int tmp_time, int me_time)
{
  if (input->Verbose == 1)
  {
    printf ("%04d(SP) %8d   %2d %7.3f %7.3f %7.3f %9d %7d    %3s    %d\n",
      frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n),  
      img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
      img->fld_flag ? "FLD" : "FRM", img->nal_reference_idc);
  }
  else if (input->Verbose == 2)
  {    
    printf ("%04d(SP) %8d %1d %2d %7.3f %7.3f %7.3f %9d %7d    %3s %5d   %2d %2d  %d   %d\n",
      frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n), active_pps->weighted_pred_flag, 
      img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
      img->fld_flag ? "FLD" : "FRM", intras, img->num_ref_idx_l0_active, img->num_ref_idx_l1_active,img->rd_pass, img->nal_reference_idc);
  }
}

static void ReportB(int tmp_time, int me_time)
{
  if (input->Verbose == 1)
  {
#ifdef ADAPTIVE_FILTER
    printf ("%04d(B)  %8d   %2d %7.3f %7.3f %7.3f %9d %7d    %3s    %d  %2d\n",
      frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n),
      img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time,me_time,
      img->fld_flag ? "FLD" : "FRM", img->nal_reference_idc, img->AdaptiveFilterFlag);
#else
    printf ("%04d(B)  %8d   %2d %7.3f %7.3f %7.3f %9d %7d    %3s    %d\n",
      frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n),
      img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time,me_time,
      img->fld_flag ? "FLD" : "FRM", img->nal_reference_idc);
#endif
  }
  else if (input->Verbose == 2)
  {
    printf ("%04d(B)  %8d %1d %2d %7.3f %7.3f %7.3f %9d %7d    %3s %5d %1d %2d %2d  %d   %d\n",
      frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n), active_pps->weighted_bipred_idc,
      img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time,me_time,
      img->fld_flag ? "FLD" : "FRM",intras,img->direct_spatial_mv_pred_flag, img->num_ref_idx_l0_active, img->num_ref_idx_l1_active,img->rd_pass, img->nal_reference_idc);
  }
}


static void ReportP(int tmp_time, int me_time)
{            
  if (input->Verbose == 1)
  {
    if(input->redundant_pic_flag==0)
    {
#ifdef ADAPTIVE_FILTER
      printf ("%04d(P)  %8d   %2d %7.3f %7.3f %7.3f %9d %7d    %3s    %d  %2d\n",
        frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n), 
        img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
        img->fld_flag ? "FLD" : "FRM", img->nal_reference_idc, img->AdaptiveFilterFlag);
#else
      printf ("%04d(P)  %8d   %2d %7.3f %7.3f %7.3f %9d %7d    %3s    %d\n",
        frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n), 
        img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
        img->fld_flag ? "FLD" : "FRM", img->nal_reference_idc);
#endif
    }
    else
    {
      if(!redundant_coding)
      {
        printf ("%04d(P)  %8d   %2d %7.3f %7.3f %7.3f %9d %7d    %3s    %d\n",
          frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n),
          img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
          img->fld_flag ? "FLD" : "FRM", img->nal_reference_idc);
      }
      else
      { // report a redundant picture.
        printf ("    (R)  %8d   %2d %7.3f %7.3f %7.3f %9d %7d    %3s    %d\n",
          (int) (stats->bit_ctr - stats->bit_ctr_n),
          img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
          img->fld_flag ? "FLD" : "FRM", img->nal_reference_idc);
      }
    }
  }
  else if (input->Verbose == 2)
  {
    printf ("%04d(P)  %8d %1d %2d %7.3f %7.3f %7.3f %9d %7d    %3s %5d   %2d %2d  %d   %d\n",
      frame_no, (int) (stats->bit_ctr - stats->bit_ctr_n), active_pps->weighted_pred_flag, 
      img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
      img->fld_flag ? "FLD" : "FRM", intras, img->num_ref_idx_l0_active, img->num_ref_idx_l1_active, img->rd_pass, img->nal_reference_idc);
  }
}

/*!
************************************************************************
* \brief
*    Padding of automatically added border for picture sizes that are not
*     multiples of macroblock/macroblock pair size
*
* \param org_size_x
*    original image horizontal size (luma)
* \param org_size_y
*    original image vertical size (luma)
* \param img_size_x
*    coded image horizontal size (luma)
* \param img_size_y
*    code image vertical size (luma)
* \param org_size_x_cr
*    original image horizontal size (chroma)
* \param org_size_y_cr
*    original image vertical size (chroma)
* \param img_size_x_cr
*    coded image horizontal size (chroma)
* \param img_size_y_cr
*    code image vertical size (chroma)
************************************************************************
*/
static void PaddAutoCropBorders (int org_size_x, int org_size_y, int img_size_x, int img_size_y,
                                 int org_size_x_cr, int org_size_y_cr, int img_size_x_cr, int img_size_y_cr)
{
  int x, y;

  //padding right border
  for (y=0; y<org_size_y; y++)
    for (x=org_size_x; x<img_size_x; x++)
      imgY_org_frm [y][x] = imgY_org_frm [y][x-1];

  //padding bottom border
  for (y=org_size_y; y<img_size_y; y++)
    for (x=0; x<img_size_x; x++)
      imgY_org_frm [y][x] = imgY_org_frm [y-1][x];


  if (img->yuv_format != YUV400)
  {
    //padding right border
    for (y=0; y<org_size_y_cr; y++)
      for (x=org_size_x_cr; x<img_size_x_cr; x++)
      {
        imgUV_org_frm [0][y][x] = imgUV_org_frm [0][y][x-1];
        imgUV_org_frm [1][y][x] = imgUV_org_frm [1][y][x-1];
      }

      //padding bottom border
      for (y=org_size_y_cr; y<img_size_y_cr; y++)
        for (x=0; x<img_size_x_cr; x++)
        {
          imgUV_org_frm [0][y][x] = imgUV_org_frm [0][y-1][x];
          imgUV_org_frm [1][y][x] = imgUV_org_frm [1][y-1][x];
        }
  }
}


/*!
************************************************************************
* \brief
*    Calculates the absolute frame number in the source file out
*    of various variables in img-> and input->
* \return
*    frame number in the file to be read
* \par side effects
*    global variable frame_no updated -- dunno, for what this one is necessary
************************************************************************
*/
static int CalculateFrameNumber(void)
{
  if (img->b_frame_to_code)
  {
    if (input->HierarchicalCoding)
      frame_no = start_tr_in_this_IGOP + (IMG_NUMBER - 1) * (input->jumpd + 1) + (int) (img->b_interval * (double) (1 + gop_structure[img->b_frame_to_code - 1].display_no));
    else      
      frame_no = start_tr_in_this_IGOP + (IMG_NUMBER - 1) * (input->jumpd + 1) + (int) (img->b_interval * (double) img->b_frame_to_code);
  }
  else
  {
    frame_no = start_tr_in_this_IGOP + IMG_NUMBER * (input->jumpd + 1);
#ifdef _ADAPT_LAST_GROUP_
    if (input->last_frame && img->number + 1 == input->no_frames)
      frame_no = input->last_frame;
#endif
  }

  return frame_no;
}


/*!
************************************************************************
* \brief
*    Convert file read buffer to source picture structure
************************************************************************
*/
void buf2img ( imgpel** imgX,           //!< Pointer to image plane
              unsigned char* buf,      //!< Buffer for file output
              int size_x,              //!< horizontal size of picture
              int size_y,              //!< vertical size of picture
              int symbol_size_in_bytes //!< number of bytes in file used for one pixel
              )
{
  int i,j;

  unsigned short tmp16, ui16;
  unsigned long  tmp32, ui32;

  if (symbol_size_in_bytes> sizeof(imgpel))
  {
    error ("Source picture has higher bit depth than imgpel data type. Please recompile with larger data type for imgpel.", 500);
  }

  if (( sizeof(char) == sizeof (imgpel)) && ( sizeof(char) == symbol_size_in_bytes))
  {
    // imgpel == pixel_in_file == 1 byte -> simple copy
    for(j=0;j<size_y;j++)
      memcpy(imgX[j], buf+j*size_x, size_x);
  }
  else
  {
    // sizeof (imgpel) > sizeof(char)
    if (testEndian())
    {
      // big endian
      switch (symbol_size_in_bytes)
      {
      case 1:
        {
          for(j=0;j<size_y;j++)
            for(i=0;i<size_x;i++)
            {
              imgX[j][i]= buf[i+j*size_x];
            }
            break;
        }
      case 2:
        {
          for(j=0;j<size_y;j++)
            for(i=0;i<size_x;i++)
            {
              memcpy(&tmp16, buf+((i+j*size_x)*2), 2);
              ui16  = (tmp16 >> 8) | ((tmp16&0xFF)<<8);
              imgX[j][i] = (imgpel) ui16;
            }
            break;
        }
      case 4:
        {
          for(j=0;j<size_y;j++)
            for(i=0;i<size_x;i++)
            {
              memcpy(&tmp32, buf+((i+j*size_x)*4), 4);
              ui32  = ((tmp32&0xFF00)<<8) | ((tmp32&0xFF)<<24) | ((tmp32&0xFF0000)>>8) | ((tmp32&0xFF000000)>>24);
              imgX[j][i] = (imgpel) ui32;
            }
        }
      default:
        {
          error ("reading only from formats of 8, 16 or 32 bit allowed on big endian architecture", 500);
          break;
        }
      }
    }
    else
    {
      // little endian
      for (j=0; j < size_y; j++)
        for (i=0; i < size_x; i++)
        {
          imgX[j][i]=0;
          memcpy(&(imgX[j][i]), buf +((i+j*size_x)*symbol_size_in_bytes), symbol_size_in_bytes);
        }
    }
  }
}


#ifdef  INTERNAL_BIT_DEPTH_INCREASE
/*!
************************************************************************
* \brief
*    Pixel Bit-Depth Increase
************************************************************************
*/
void ImgBitDepthIncrease ( imgpel** imgX,           //!< Pointer to image plane
                          int size_x,              //!< horizontal size of picture
                          int size_y,              //!< vertical size of picture
                          int input_bit_depth,
                          int bit_depth_increase
                          )
{
  int i,j;

  for(j=0;j<size_y;j++)
    for(i=0;i<size_x;i++)
    {
      imgX[j][i] = imgX[j][i]<<bit_depth_increase;
    }

}
#endif


/*!
************************************************************************
* \brief
*    Reads one new frame from file
*
* \param FrameNoInFile
*    Frame number in the source file
* \param HeaderSize
*    Number of bytes in the source file to be skipped
* \param xs
*    horizontal size of frame in pixels
* \param ys
*    vertical size of frame in pixels
* \param xs_cr
*    horizontal chroma size of frame in pixels
* \param ys_cr
*    vertical chroma size of frame in pixels
************************************************************************
*/
static void ReadOneFrame (int FrameNoInFile, int HeaderSize, int xs, int ys, int xs_cr, int ys_cr)
{
  unsigned int symbol_size_in_bytes = img->pic_unit_size_on_disk/8;

  const int imgsize_y = xs*ys;
  const int imgsize_uv = xs_cr*ys_cr;

  const int bytes_y = imgsize_y * symbol_size_in_bytes;
  const int bytes_uv = imgsize_uv * symbol_size_in_bytes;

  const int64 framesize_in_bytes = bytes_y + 2*bytes_uv;
  unsigned char *buf;

  Boolean rgb_input = (Boolean)(input->rgb_input_flag==1 && input->yuv_format==3);

  assert (p_in != -1);

  // KS: this buffer should actually be allocated only once, but this is still much faster than the previous version
  if (NULL==(buf = malloc (xs*ys * symbol_size_in_bytes))) no_mem_exit("ReadOneFrame: buf");

  // skip Header
  if (lseek (p_in, HeaderSize, SEEK_SET) != HeaderSize)
  {
    error ("ReadOneFrame: cannot fseek to (Header size) in p_in", -1);
  }

  // skip starting frames
  if (lseek (p_in, framesize_in_bytes * input->start_frame, SEEK_CUR) == -1) 
  {
    snprintf(errortext, ET_SIZE, "ReadOneFrame: cannot advance file pointer in p_in beyond frame %d\n", input->start_frame);
    error (errortext,-1);
  } 

  // seek to current frame 
  if (lseek (p_in, framesize_in_bytes * FrameNoInFile, SEEK_CUR) == -1) 
  {
    snprintf(errortext, ET_SIZE, "ReadOneFrame: cannot advance file pointer in p_in beyond frame %d\n", input->start_frame + FrameNoInFile);
    error (errortext,-1);
  }

  // Here we are at the correct position for the source frame in the file.  Now
  // read it.
  if (img->pic_unit_size_on_disk%8 == 0)
  {
    if(rgb_input)
      lseek (p_in, framesize_in_bytes/3, SEEK_CUR);

    if (read(p_in, buf, bytes_y) != bytes_y)
    {
      printf ("ReadOneFrame: cannot read %d bytes from input file, unexpected EOF?, exiting", bytes_y);
      report_stats_on_error();
      exit (-1);
    }

    buf2img(imgY_org_frm, buf, xs, ys, symbol_size_in_bytes);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
    if (img->BitDepthIncrease)
    {
      ImgBitDepthIncrease(imgY_org_frm, xs, ys, input->BitDepthLuma, img->BitDepthIncrease);
    }
#endif

    if (img->yuv_format != YUV400)
    {
      if (read(p_in, buf, bytes_uv) != bytes_uv)
      {
        printf ("ReadOneFrame: cannot read %d bytes from input file, unexpected EOF?, exiting", bytes_y);
        report_stats_on_error();
        exit (-1);
      }
      buf2img(imgUV_org_frm[0], buf, xs_cr, ys_cr, symbol_size_in_bytes);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
      if (img->BitDepthIncreaseChroma)
      {
        ImgBitDepthIncrease(imgUV_org_frm[0], xs_cr, ys_cr, input->BitDepthChroma, img->BitDepthIncreaseChroma);
      }
#endif

      if(rgb_input)
        lseek (p_in, -framesize_in_bytes, SEEK_CUR);

      if (read(p_in, buf, bytes_uv) != bytes_uv)
      {
        printf ("ReadOneFrame: cannot read %d bytes from input file, unexpected EOF?, exiting", bytes_y);
        report_stats_on_error();
        exit (-1);
      }
      buf2img(imgUV_org_frm[1], buf, xs_cr, ys_cr, symbol_size_in_bytes);
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
      if (img->BitDepthIncreaseChroma)
      {
        ImgBitDepthIncrease(imgUV_org_frm[1], xs_cr, ys_cr, input->BitDepthChroma, img->BitDepthIncreaseChroma);
      }
#endif

      if(rgb_input)
        lseek (p_in, framesize_in_bytes*2/3, SEEK_CUR);
    }
  }
  else
  {
    printf ("ReadOneFrame (NOT IMPLEMENTED): pic unit size on disk must be divided by 8");
    exit (-1);
  }
  free (buf);

}


/*!
************************************************************************
* \brief
*    point to frame coding variables 
************************************************************************
*/
static void put_buffer_frame(void)
{
  imgY_org  = imgY_org_frm;
  imgUV_org = imgUV_org_frm;  
}

/*!
************************************************************************
* \brief
*    point to top field coding variables 
************************************************************************
*/
static void put_buffer_top(void)
{
  img->fld_type = 0;

  imgY_org = imgY_org_top;
  imgUV_org = imgUV_org_top;
}

/*!
************************************************************************
* \brief
*    point to bottom field coding variables 
************************************************************************
*/
static void put_buffer_bot(void)
{
  img->fld_type = 1;

  imgY_org = imgY_org_bot;
  imgUV_org = imgUV_org_bot;
}

/*!
************************************************************************
* \brief
*    Writes a NAL unit of a partition or slice
************************************************************************
*/

static void writeUnit(Bitstream* currStream, int partition)
{
  const int buffer_size = 500 + img->FrameSizeInMbs * (128 + 256 * img->bitdepth_luma + 512 * img->bitdepth_chroma);
  // KS: this is approx. max. allowed code picture size
  NALU_t *nalu;
  assert (currStream->bits_to_go == 8);
  nalu = AllocNALU(buffer_size);
  nalu->startcodeprefix_len = 2+(img->current_mb_nr == 0?ZEROBYTES_SHORTSTARTCODE+1:ZEROBYTES_SHORTSTARTCODE);
  //printf ("nalu->startcodeprefix_len %d\n", nalu->startcodeprefix_len);
  nalu->len = currStream->byte_pos +1;            // add one for the first byte of the NALU
  //printf ("nalu->len %d\n", nalu->len);
  memcpy (&nalu->buf[1], currStream->streamBuffer, nalu->len-1);
  if (img->currentPicture->idr_flag)
  {
    nalu->nal_unit_type = NALU_TYPE_IDR;
    nalu->nal_reference_idc = NALU_PRIORITY_HIGHEST;
  }
  else if (img->type == B_SLICE)
  {
    //different nal header for different partitions
    if(input->partition_mode == 0)
    {
      nalu->nal_unit_type = NALU_TYPE_SLICE;
    }
    else
    {
      nalu->nal_unit_type = NALU_TYPE_DPA +  partition;
    }

    if (img->nal_reference_idc !=0)
    {
      nalu->nal_reference_idc = NALU_PRIORITY_HIGH;
    }
    else
    {
      nalu->nal_reference_idc = NALU_PRIORITY_DISPOSABLE;
    }
  }
  else   // non-b frame, non IDR slice
  {
    //different nal header for different partitions
    if(input->partition_mode == 0)
    {
      nalu->nal_unit_type = NALU_TYPE_SLICE;
    }
    else
    {
      nalu->nal_unit_type = NALU_TYPE_DPA +  partition;
    }
    if (img->nal_reference_idc !=0)
    {
      nalu->nal_reference_idc = NALU_PRIORITY_HIGH;
    }
    else
    {
      nalu->nal_reference_idc = NALU_PRIORITY_DISPOSABLE;
    }
  }      
  nalu->forbidden_bit = 0;
  stats->bit_ctr += WriteNALU (nalu);

  FreeNALU(nalu);
}

/*!
************************************************************************
* \brief
*    performs multi-pass encoding of same picture using different 
*    coding conditions
************************************************************************
*/
static void rdPictureCoding(void)
{
  int second_qp = img->qp, rd_qp = img->qp;
  int previntras = intras;
  int prevtype = img->type;
  int skip_encode = 0;
  pic_parameter_set_rbsp_t *sec_pps;
#ifdef AIF_BUG_FIX
  int change_aif_flag = 0;
#endif
#ifdef SIMPLIFIED_RDPIC_DECISION
  int pps_id;
#endif

  if (img->type!=I_SLICE && input->GenerateMultiplePPS)
  {
    if (img->type==P_SLICE)
    {
      if (test_wp_P_slice(0) == 1)
      {
        active_pps = PicParSet[1];
      }
      else
      {
        skip_encode = input->RDPSliceWeightOnly;
        active_pps = PicParSet[0];
        if (!img->AdaptiveRounding)
          img->qp-=1;
      }
    }
    else
    {
#ifdef SIMPLIFIED_RDPIC_DECISION
      if(input->SimplifiedRDPicDecision)
      {
        if (test_wp_B_slice(1) == 1)//Implicit wp based on POC...skip this pass if weights are (0.5,0.5)
        {
          active_pps = PicParSet[2];		
        }
        else
        {
          skip_encode = 1;
        }
      }
      else
        active_pps = PicParSet[2];
#else
      active_pps = PicParSet[2];
#endif
    }
  }
  else        
  {
    if (!img->AdaptiveRounding)
      img->qp-=1;
  }

  sec_pps = active_pps;
  second_qp = img->qp;

  img->write_macroblock = 0;

  if (skip_encode)
  {
    img->rd_pass = 0;
    enc_frame_picture2 = NULL;
  }
  else
  {
    frame_picture (frame_pic_2,1);
#ifdef ADAPTIVE_FILTER
    if(img->AdaptiveFilterFlag==1)
      img->rd_pass=picture_coding_decision(frame_pic_aif, frame_pic_2, rd_qp);
    else
#endif
      img->rd_pass=picture_coding_decision(frame_pic_1, frame_pic_2, rd_qp);
  }
  //      update_rd_picture_contexts (img->rd_pass); 
  if (img->rd_pass==0)
  {
    enc_picture=enc_frame_picture;
    if (img->type!=I_SLICE && input->GenerateMultiplePPS)
    { 
      img->qp=rd_qp;
      active_pps = PicParSet[0];
    }
    else       
    {
      img->qp=rd_qp;
    }
    intras = previntras;
#ifdef AIF_BUG_FIX
    if(img->AdaptiveFilterFlag == 1)
      frame_pic = frame_pic_aif;
    else
      frame_pic = frame_pic_1;
#else
    frame_pic = frame_pic_1;
#endif
  }
  else
  {
    previntras = intras;
    frame_pic = frame_pic_2;
  }

#ifdef SIMPLIFIED_RDPIC_DECISION
  pps_id = active_pps->pic_parameter_set_id;
#endif

  // Final Encoding pass - note that we should 
  // make this more flexible in a later version.

  if (img->type!=I_SLICE && input->GenerateMultiplePPS)
  {
    skip_encode = 0;
    img->qp    = rd_qp;

#ifndef USE_NEW_OFFSET
    if (img->type == P_SLICE && input->GenerateMultiplePPS && (intras * 100 )/img->FrameSizeInMbs >=75)
    {
      img->type=I_SLICE;
      active_pps = PicParSet[0];
#ifdef AIF_BUG_FIX
      if(img->AdaptiveFilterFlag == 1)
      {
        img->AdaptiveFilterFlag = -1;
        change_aif_flag = 1;
      }
#endif
    }
    else if (img->type==P_SLICE)
#else
    if (!input->UseNewOffset && (img->type == P_SLICE && input->GenerateMultiplePPS && (intras * 100 )/img->FrameSizeInMbs >=75))
    {
      img->type=I_SLICE;
      active_pps = PicParSet[0];
#ifdef AIF_BUG_FIX
      if(img->AdaptiveFilterFlag == 1)
      {
        img->AdaptiveFilterFlag = -1;
        change_aif_flag = 1;
      }
#endif
    }
    else if (img->type==P_SLICE)
#endif 
    {
      if (test_wp_P_slice(1) == 1)
      {
        active_pps = PicParSet[1];
      }
      else if (input->RDPSliceBTest && active_sps->profile_idc != 66)
      {
        img->type = B_SLICE;
        active_pps = PicParSet[0];        
      }
      else
      {
        skip_encode = input->RDPSliceWeightOnly;
        active_pps = PicParSet[0];
        if (!img->AdaptiveRounding)
          img->qp+=1;
      }
    }
    else
    {
      if (test_wp_B_slice(0) == 1)
      {
        active_pps = PicParSet[1];
      }
      else
      {
        skip_encode = input->RDBSliceWeightOnly;
#ifdef SIMPLIFIED_RDPIC_DECISION
        if(input->SimplifiedRDPicDecision)
        {
          if(!input->UseRDO_Q)
            img->qp = rd_qp + (img->nal_reference_idc ? - 1 : 1);
        }
        else
          img->qp = rd_qp + (img->nal_reference_idc ? - 1 : 1);
#else
        img->qp = rd_qp + (img->nal_reference_idc ? - 1 : 1);
#endif
      }      
    }
  }
  else 
  {
    active_pps = PicParSet[0];
    if (!img->AdaptiveRounding)
      img->qp    = (rd_qp + 1);
  }

  img->write_macroblock = 0;

#ifdef SIMPLIFIED_RDPIC_DECISION
	if (skip_encode || ((input->SimplifiedRDPicDecision==1) && (img->type==B_SLICE) && (pps_id==active_pps->pic_parameter_set_id) && (input->UseRDO_Q)))
#else
  if (skip_encode)
#endif
  {
    enc_frame_picture3 = NULL;
    img->qp = rd_qp;
  }
  else
  {
    frame_picture (frame_pic_3,2);

#ifdef ADAPTIVE_FILTER
#ifdef AIF_BUG_FIX
    if(change_aif_flag == 1)
    {
      img->AdaptiveFilterFlag = 1;
    }
#endif

    if (img->AdaptiveFilterFlag==1)
    {
      if (img->rd_pass==0)
        img->rd_pass  = 2*picture_coding_decision(frame_pic_aif, frame_pic_3, rd_qp);
      else
        img->rd_pass +=   picture_coding_decision(frame_pic_2, frame_pic_3, rd_qp);
    }
    else
    {
      if (img->rd_pass==0)
        img->rd_pass  = 2*picture_coding_decision(frame_pic_1, frame_pic_3, rd_qp);
      else
        img->rd_pass +=   picture_coding_decision(frame_pic_2, frame_pic_3, rd_qp);
    }
#else    
    if (img->rd_pass==0)
      img->rd_pass  = 2*picture_coding_decision(frame_pic_1, frame_pic_3, rd_qp);
    else
      img->rd_pass +=   picture_coding_decision(frame_pic_2, frame_pic_3, rd_qp);
#endif
  }

  //update_rd_picture_contexts (img->rd_pass); 
  if (img->rd_pass==0)
  {
#ifdef AIF_BUG_FIX
    if(img->AdaptiveFilterFlag == 1)
      enc_picture = enc_frame_picture_aif;
    else
      enc_picture = enc_frame_picture;
#else
    enc_picture = enc_frame_picture;
#endif
    img->type   = prevtype;
    active_pps  = PicParSet[0];
    img->qp     = rd_qp;
    intras      = previntras;
  }
  else if (img->rd_pass==1)
  {
    enc_picture = enc_frame_picture2;
    img->type   = prevtype;
    active_pps  = sec_pps;
    img->qp     = second_qp;
    intras      = previntras;
  }       
}



/*!
*************************************************************************************
* Brief
*     Output SP frames coefficients
*************************************************************************************
*/
void output_SP_coefficients()
{
  int i,k;
  FILE *SP_coeff_file;
  if(number_sp2_frames==0) 
  {
    if ((SP_coeff_file = fopen(input->sp_output_filename,"w")) == NULL)
    {
      printf ("Fatal: cannot open SP output file '%s', exit (-1)\n", input->sp_output_filename);
      exit (-1);
    }
    number_sp2_frames++;
  }
  else
  {
    if ((SP_coeff_file = fopen(input->sp_output_filename,"a")) == NULL)
    {
      printf ("Fatal: cannot open SP output file '%s', exit (-1)\n", input->sp_output_filename);
      exit (-1);
    }
  }

  for(i=0;i<img->height;i++)
  {
    fwrite(lrec[i],sizeof(int),img->width,SP_coeff_file);
  }
  for(k=0;k<2;k++)
  {
    for(i=0;i<img->height_cr;i++)
    {
      fwrite(lrec_uv[k][i],sizeof(int),img->width_cr,SP_coeff_file);
    }
  }
  fclose(SP_coeff_file);
}

/*!
*************************************************************************************
* Brief
*     Read SP frames coefficients
*************************************************************************************
*/
void read_SP_coefficients()
{
  int i,k;
  FILE *SP_coeff_file;

  if ( (input->qp2start > 0) && ( ( (img->tr ) % (2*input->qp2start) ) >=input->qp2start ))
  {
    if ((SP_coeff_file = fopen(input->sp2_input_filename1,"r")) == NULL)
    {
      printf ("Fatal: cannot open SP input file '%s', exit (-1)\n", input->sp2_input_filename2);
      exit (-1);
    }   
  }
  else
  {
    if ((SP_coeff_file = fopen(input->sp2_input_filename2,"r")) == NULL)
    {
      printf ("Fatal: cannot open SP input file '%s', exit (-1)\n", input->sp2_input_filename1);
      exit (-1);
    }   
  }

  if (0 != fseek (SP_coeff_file, img->width*img->height*3/2*number_sp2_frames*sizeof(int), SEEK_SET))
  {
    printf ("Fatal: cannot seek in SP input file, exit (-1)\n");
    exit (-1);
  }     
  number_sp2_frames++;

  for(i=0;i<img->height;i++) 
  {
    if(img->width!=(int)fread(lrec[i],sizeof(int),img->width,SP_coeff_file))
    { 
      printf ("Fatal: cannot read in SP input file, exit (-1)\n");
      exit (-1);
    } 
  }

  for(k=0;k<2;k++)
  {
    for(i=0;i<img->height_cr;i++) 
    {
      if(img->width_cr!=(int)fread(lrec_uv[k][i],sizeof(int),img->width_cr,SP_coeff_file))
      { 
        printf ("Fatal: cannot read in SP input file, exit (-1)\n");
        exit (-1);
      } 
    }
  }
  fclose(SP_coeff_file);
}



