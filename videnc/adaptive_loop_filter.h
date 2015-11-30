
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

#ifdef ADAPTIVE_LOOP_FILTER

//------- defines
#define MIN_NUM_TAP    5
#define MAX_NUM_TAP    9
#define MAX_NUM_TAP_C  5
#define NUM_BIT_SHIFT  8
#define MAX_NUM_COEF   42
#define MIN_NUM_COEF   14
#define MAX_NUM_COEF_C 14
#define MAX_ALF_PARAM_SIZE 50000
#define NUM_OF_REDESIGN 3
#define MIN_NUM_BLOCK  10
#define MAX_NUM_BLOCK  60000
#define PYRAMIDAL_MEM_SIZE 341
#define MIN_NUM_PARENT_BLOCK 7
typedef struct
{
  Boolean      valid_flag[PYRAMIDAL_MEM_SIZE];
  Boolean      part_flag[PYRAMIDAL_MEM_SIZE];
  Boolean      alf_blk_flag[PYRAMIDAL_MEM_SIZE];
} AlfQtBlock;
//------- end

//------- prototype
alf_parameter_set_rbsp_t *AllocALFPS(int *memory_size);
void   FreeALFPS (alf_parameter_set_rbsp_t *alfps);
int    InitALFGlobalBuffers(void);
void   FreeALFGlobalBurrers(void);
void   SetALFParameters(alf_parameter_set_rbsp_t *alfps, alf_parameter_set_rbsp_t *p);
void   MakeExtendFrame( imgpel** ImgDec, imgpel** ImgExt, int width, int height, int offset );
void   SetRestFrame(imgpel** ImgRest, imgpel** ImgDec, int width, int height);

//------- prediction of coeff
#ifdef ALF_SPATIAL_PREDICT_COEF
void   PredictALFCoeff( alf_parameter_set_rbsp_t *alfps );
void   PredictALFCoeffChroma( alf_parameter_set_rbsp_t *alfps );
#endif
int    GetALFBaseType(int slice_type, int ref);
void   InitFilterBase();
void   SetALFDiffCoef(alf_parameter_set_rbsp_t *alfps);
void   SetALFBaseCoef(alf_parameter_set_rbsp_t *alfps);
void   SetALFDiffCoefChroma(alf_parameter_set_rbsp_t *alfps);
void   SetALFBaseCoefChroma(alf_parameter_set_rbsp_t *alfps);

void   CalcALFLambda(int qp);
int64  FindDistortionFrame(imgpel** ImgOrg, imgpel** ImgCmp);
void   FindDistortionFrameChroma(imgpel*** ImgUVOrg, imgpel*** ImgUVCmp, int64 *dist_u, int64 *dist_v);
int64  FindDistortionBlock(int blk_addr, imgpel** ImgOrg, imgpel** ImgCmp, int blk_size, int width_in_alf_blks);
void   CalcRDCost(imgpel** ImgOrg, imgpel** ImgCmp, alf_parameter_set_rbsp_t *alfps, int64* rate, int64* dist, double* cost);
void   CalcRDCostChroma(imgpel*** ImgUVOrg, imgpel*** ImgUVCmp, alf_parameter_set_rbsp_t *alfps, int64* rate, int64* dist, double* cost);
void   CalcRDCostBlock(imgpel** ImgOrg, imgpel** ImgCmp, imgpel** ImgDec, alf_parameter_set_rbsp_t *alfps, int64* rate, int64* dist, double* cost);
void   CalcCorrelationFunc(imgpel** ImgOrg, imgpel** ImgDec, double** corr, int tap, int width, int height);
void   CalcStoredCorrelationFuncBlock(int blk_addr, imgpel** ImgOrg, imgpel** ImgDec, unsigned int** corr, int tap);
void   CalcCorrelationFuncBlock(int blk_addr, alf_parameter_set_rbsp_t *alfps, imgpel** ImgOrg, imgpel** ImgDec, double** corr);
int    Gauss(double **a, int N);
void   FilterCoefQuickSort( double *coef_data, int *coef_num, int upper, int lower );
void   QuantFilterCoef(double* h, int* qh, int tap, int bitdepth);
void   ClearFilterCoefInt(int* qh, int N);

//------- filtering for luma
void   AdaptiveLoopFilterProcess(Picture *pic);
void   FirstFilteringFrameLuma(imgpel** ImgOrg, imgpel** ImgDec, imgpel** ImgRest);
void   FilteringFrameLuma(imgpel** ImgOrg, imgpel** ImgDec, imgpel** ImgRest, int store_corr);
void   FilteringProcessLuma(imgpel** ImgDec, imgpel** ImgRest, int *qh, int tap);
void   FilterTapDecision(imgpel** ImgOrg, imgpel** ImgDec, imgpel** ImgRest);

//------- filtering for chroma
void   AdaptiveLoopFilterProcessChroma();
void   FilteringFrameChroma(imgpel ***ImgUVOrg, imgpel ***ImgUVDec, imgpel ***ImgUVRest);
void   FilteringProcessChroma(imgpel ***ImgUVDec, imgpel ***ImgUVRest, int *qh, int tap, int chroma_id);

//------- filtering based block
void   BlockAdaptiveFilterControl(imgpel** ImgOrg, imgpel** ImgDec, imgpel** ImgRest);
void   ReDesignFilterCoeff(imgpel** ImgOrg, imgpel** ImgDec, int read_corr);
void   SetBlockSizeParameter(alf_parameter_set_rbsp_t *alfps, int size_idx);
int    CalcNumOfAlfBlocks(int blk_size, int width, int height);
int    SetALFBlockFlag(alf_parameter_set_rbsp_t *alfps, imgpel** ImgOrg, imgpel** ImgDec, imgpel** ImgRest);
void   CopyDecToRestBlock(alf_parameter_set_rbsp_t *alfps, imgpel** ImgDec, imgpel** ImgRest);
void   GetBlockPosForBALF(int blk_addr, int *x, int *y, int width_in_blks);
void   GetSamplePosForBALF(int blk_addr, int *x, int *y, int blk_size, int width_in_blks);
int    get_mem_blk_corr(unsigned int *****blk_corr, int height, int width, int rows, int columns);
void   free_mem_blk_corr(unsigned int ****blk_corr, int height, int width );
int    ALFBlockIsAvailable(int blk_addr);
void   CheckAvailabilityOfNeighborsAlf();
void   GetNeighbourAlfBlock(int curr_alfb_nr, int xN, int yN, PixelPos *pix);
void   SetBlockPartitioningFlag(int pblk_addr, alf_parameter_set_rbsp_t *alfps, imgpel **ImgOrg, imgpel **ImgDec, imgpel **ImgRest);
void   TransQbpToFbp(int pblk_addr, alf_parameter_set_rbsp_t *alfps);
void   SetParentBlockSizeParam(alf_parameter_set_rbsp_t *alfps, int block_size);

//------- send ALF parameters
int    send_ALFParameterSet(alf_parameter_set_rbsp_t *alfps, Bitstream *bitstream, imgpel **ImgOrg, imgpel **ImgDec, imgpel **ImgRest);
int    send_ALFLumaParam(alf_parameter_set_rbsp_t *alfps, Bitstream *bitstream);
int    send_ALFChromaParam(alf_parameter_set_rbsp_t *alfps, Bitstream *bitstream);
int    send_ALFBlockControlParam(alf_parameter_set_rbsp_t *alfps, Bitstream *bitstream, imgpel **ImgOrg, imgpel **ImgDec, imgpel **ImgRest);
int    send_BlockFlagRunLength(alf_parameter_set_rbsp_t *alfps, Bitstream *bitstream);
int    send_BlockFlagCABAC(alf_parameter_set_rbsp_t *alfps, Bitstream *bitstream);
int    send_QuadTreeBlockPartitioningParam(alf_parameter_set_rbsp_t *alfps, Bitstream *bitstream, imgpel **ImgOrg, imgpel **ImgDec, imgpel **ImgRest);
int    TerminateSliceForALF(Picture* pic);


#endif // ADAPTIVE_LOOP_FILTER

#endif // _ADAPTIVE_LOOP_FILTER_H_
