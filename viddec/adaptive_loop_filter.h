
/*!
 *************************************************************************************
 * \file adaptive_loop_filter.h
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

#ifndef _ADAPTIVE_LOOP_FILTER_H_
#define _ADAPTIVE_LOOP_FILTER_H_

#include "global.h"
#include "mbuffer.h"

#ifdef ADAPTIVE_LOOP_FILTER

//------- defines
#define MAX_NUM_TAP    9
#define MAX_NUM_TAP_C  5
#define NUM_BIT_SHIFT  8
#define MAX_NUM_COEF   42
#define MAX_NUM_COEF_C 14
#define MAX_ALF_PARAM_SIZE 50000
#define MAX_NUM_BLOCK  60000
#define PYRAMIDAL_MEM_SIZE 341
typedef struct
{
  Boolean      valid_flag[PYRAMIDAL_MEM_SIZE];
  Boolean      part_flag[PYRAMIDAL_MEM_SIZE];
  Boolean      alf_blk_flag[PYRAMIDAL_MEM_SIZE];
} AlfQtBlock;
//------- end


//------- prototype
alf_parameter_set_rbsp_t *AllocALFPS (int *memory_size);
void   FreeALFPS (alf_parameter_set_rbsp_t *alfps);
int    InitALFGlobalBuffers(void);
void   FreeALFGlobalBurrers(void);
void   SetALFParameters(StorablePicture *p, alf_parameter_set_rbsp_t *alfps);
void   MakeExtendFrame( imgpel** ImgDec, imgpel** ImgExt, int width, int height, int offset );
void   SetRestFrame(imgpel** ImgRest, imgpel** ImgDec, int width, int height);
void   SetALFImgSize(ImageParameters *img, seq_parameter_set_rbsp_t *p);

//------- prediction of coeff
#ifdef ALF_SPATIAL_PREDICT_COEF
void   PredictALFCoeff(StorablePicture *p);
void   PredictALFCoeffChroma(StorablePicture *p);
#endif
int    GetALFBaseType(int slice_type, int ref);
void   InitFilterBase();
void   SetALFCoef(StorablePicture *p);
void   SetALFCoefChroma(StorablePicture *p);


//------- filtering for luma
void   AdaptiveLoopFilterProcess(StorablePicture *p);
void   FilteringFrame(StorablePicture *p, imgpel** ImgDec, imgpel** ImgRest);
void   FilteringProcess(imgpel** ImgDec, imgpel** ImgRest, int *qh, int tap);

//------- filtering for chroma
void   AdaptiveLoopFilterProcessChroma(StorablePicture *p);
void   FilteringFrameChroma(StorablePicture *p, imgpel ***ImgUVDec, imgpel ***ImgUVRest);
void   FilteringProcessChroma(imgpel ***ImgUVDec, imgpel ***ImgUVRest, int *qh, int tap, int chroma_id);

//------- filtering based block
void   BlockAdaptiveFilteringProcess(StorablePicture *p, imgpel** ImgDec, imgpel** ImgRest);
void   FilteringProcessBlock(int blk_addr, imgpel** ImgDec, imgpel** ImgRest, int *qh, int tap, int blk_size, int width_in_alf_blks);
void   CopyDecToRestBlock(int blk_addr, imgpel** ImgDec, imgpel** ImgRest, int blk_size, int width_in_alf_blks);
void   SetBlockSizeParameter(alf_parameter_set_rbsp_t *alfps, int size_idx);
int    CalcNumOfAlfBlocks(int blk_size, int width, int height);
void   GetBlockPosForBALF(int blk_addr, int *x, int *y, int width_in_blks);
void   GetSamplePosForBALF(int blk_addr, int *x, int *y, int blk_size, int width_in_blks);
int    ALFBlockIsAvailable(int blk_addr);
void   CheckAvailabilityOfNeighborsAlfBlock();
void   GetNeighbourAlfBlock(int curr_alfb_nr, int xN, int yN, PixelPos *pix);
void   InitBlockPartition(int pblk_addr, alf_parameter_set_rbsp_t *alfps);
void   TransQbpToFbp(int pblk_addr, alf_parameter_set_rbsp_t *alfps);
void   SetParentBlockSizeParam(alf_parameter_set_rbsp_t *alfps, int block_size);

//------- read ALF parameters
int    read_ALFParameterSet(alf_parameter_set_rbsp_t *alfps, DataPartition *p);
int    read_ALFLumaParam(alf_parameter_set_rbsp_t *alfps, Bitstream *s);
int    read_ALFChromaParam(alf_parameter_set_rbsp_t *alfps, Bitstream *s);
int    read_ALFBlockControlParam(alf_parameter_set_rbsp_t *alfps,  DataPartition *p);
int    read_BlockFlagRunLength(alf_parameter_set_rbsp_t *alfps, Bitstream *s);
int    read_BlockFlagCABAC(alf_parameter_set_rbsp_t *alfps, DataPartition *p);
int    read_QuadTreeBlockPartitioningParam(alf_parameter_set_rbsp_t *alfps, Bitstream *s);
#ifdef SWITCHED_FILTERS
void   read_ALFHeader();
#endif


#endif // ADAPTIVE_LOOP_FILTER

#endif // _ADAPTIVE_LOOP_FILTER_H_
