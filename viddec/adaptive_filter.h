/*!
************************************************************************
* \file adaptive_filter.h
*
* \brief
*    adaptive filter routines
*
* \author
*  Yuri Vatis               <vatis@tnt.uni-hannover.de>    \n
*
* This software may be used for research purposes only.
* Universitaet Hannover would also appreciate that all those using the 
* software or any extension of it in any way comply with the 
* following condition:
*    Including in any technical report, conference or journal paper
*    produced which uses the software or any extension of it in
*    any way, at least one of bibliographic references:
*    Yuri Vatis, Bernd Edler, Dieu Thanh Nguyen, Jörn Ostermann,
*    "Motion-and Aliasing-compensated Prediction using a two-dimensional non-separable Adaptive Wiener Interpolation Filter",
*    Proc. ICIP 2005, IEEE International Conference on Image Processing, Genova, Italy, September 2005.
*    Yuri Vatis, Bernd Edler, Ingolf Wassermann, Dieu Thanh Nguyen, Jörn Ostermann,
*    "Coding of Coefficients of two-dimensional non-separable Adaptive Wiener Interpolation Filter", 
*    Proc. VCIP 2005, SPIE Visual Communication & Image Processing, Beijing, China, July 2005.
************************************************************************
*/


#ifndef _ADAPTIVE_FILTER_H
#define _ADAPTIVE_FILTER_H
#include "mbuffer.h"
//#ifdef EIGHTH_PEL
#include "global.h"
#ifdef EDAIF2
#include "aif_common.h" // we collect in a single file defines related to encoder and decoder
#endif
//#endif

#define SQR(a)  ((a)*(a))

#ifdef E_DAIF
#define FILTER_SIZE_INT   5
#define SQR_FILTER_INT    (FILTER_SIZE_INT * FILTER_SIZE_INT + 1)

#ifndef EDAIF2
#define SQR_FILTER        37
#endif
#else
#define SQR_FILTER        36
#endif

#ifndef EDAIF2
#define DEFAULT_QUANT 9
#endif

/*
X       X       X       X       X       X



X       X       X       X       X       X



X       X       X a b c X       X       X
                d e f g
                h i j k
                l m n o
X       X       X       X       X       X



X       X       X       X       X       X



X       X       X       X       X       X
*/
#ifndef EDAIF2
#define a_pos  0         //  1
#define b_pos  1         //  2
#define c_pos  2         //  3
#define d_pos  3         //  4
#define e_pos  4         //  5
#define f_pos  5         //  6
#define g_pos  6         //  7
#define h_pos  7         //  8
#define i_pos  8         //  9
#define j_pos  9         // 10
#define k_pos 10         // 11
#define l_pos 11         // 12
#define m_pos 12         // 13
#define n_pos 13         // 14
#define o_pos 14         // 15
#endif

#define NUMBER_OF_DIF_COEFFS 54
// 54 =    3 for b_pos +
//         6 for a_pos +
//        18 for f_pos +
//         21 for e_pos +
//         6  for j_pos
/*
Ematrix*FilterCoef=CVector;
*/
#ifdef DIRECTIONAL_FILTER
#define FILTER_TYPE_2D_NS  1
#define FILTER_TYPE_2D_S   2
#define FILTER_TYPE_1D     3
#define FILTER_TYPE_EDAIF  4
#define FILTER_TYPE_EAIF   5

#define IMP_FLOAT32        0
#define IMP_INT16          1
#endif

void InitAdaptiveFilter(void);
void FreeAdaptiveFilter(void);
void ResetAdaptiveFilter(void);
int  UseAdaptiveFilterForCurrentFrame(void);
void QuantizeFilterCoef(double **QCoef, double **Coef, int N, int quant);
void PrintFilterCoef(void);
void ReadFilterCoeffitients(void);
void GetBlockWith2DAIF(int ref_idx, StorablePicture **list, int x_pos, int y_pos, int mvx, int mvy,
                       int img_width, int img_height, int block[BLOCK_SIZE][BLOCK_SIZE]);
void GetBlockWithSepAIF(int ref_idx, StorablePicture **list, int x_pos, int y_pos, int xoff4x4, int yoff4x4, int mvx, int mvy,
                        int img_width, int img_height, int block_x, int block_y, int block[MB_BLOCK_SIZE][MB_BLOCK_SIZE]);
#ifdef EIGHTH_PEL
void GetBlockWith2DAIF_quarter_pel(int ref_idx, StorablePicture **list, int x_pos, int y_pos, int mvx, int mvy,
                                   int img_width, int img_height, int block[BLOCK_SIZE][BLOCK_SIZE]);
void GetBlockWith2DAIF_eighth_pel(int ref_idx, StorablePicture **list, int x_pos, int y_pos, int mvx, int mvy,
                                  int img_width, int img_height, int block[BLOCK_SIZE][BLOCK_SIZE]);
void GetBlockWithSepAIF_quarter_pel(int ref_idx, StorablePicture **list, int x_pos, int y_pos, int xoff4x4, int yoff4x4, int mvx, int mvy,
                                    int img_width, int img_height, int block_size_x, int block_size_y, int block[MB_BLOCK_SIZE][MB_BLOCK_SIZE]); 
void GetBlockWithSepAIF_eighth_pel(int ref_idx, StorablePicture **list, int x_pos, int y_pos, int xoff4x4, int yoff4x4, int mvx, int mvy,
                                   int img_width, int img_height, int block_size_x, int block_size_y, int block[MB_BLOCK_SIZE][MB_BLOCK_SIZE]); 
#endif


int  FindPosition(int Dim, int cur_pos, int Offset, int mv);
void WriteSubBlockData(int mvx, int mvy, int cur_mb_nr, int block8x8, int block4x4);
void WriteSubBlockValues(int aif_std, int tmp_block[BLOCK_SIZE][BLOCK_SIZE]);
int  CountWithNewFilter(int mvx, int mvy);
void WriteDifferencesIntoTestFile(int cur_number, int img_height, int img_width,
                                  byte **imgY, byte **imgY_ref);
void CalculateFilterCoefficients(void);
void CalculateFilterCoefficientsSep(void); // separable aif
void ExtendFilterCoefficientsFloat(int sub_pos, double FCoef[15][SQR_FILTER]);
void ExtendFilterCoefficientsInt  (int sub_pos, int   FCoef[15][SQR_FILTER]);
void DequantizeFilterCoefficients(int QFCoef[SQR_FILTER], double FCoef[SQR_FILTER]);
void StoreFilterCoefficients(StorablePicture *p);
void PrintFilterCoefFloat(double FCoef[15][SQR_FILTER]);
void PrintFilterCoefInt(int FCoef[15][SQR_FILTER]);
#ifdef DIRECTIONAL_FILTER
void ExtendFilterCoefficientsShort(int sub_pos, short FCoef[15][SQR_FILTER]);

void initFilterCustom(int filterType);
void CalculateFilterCoefficients1DAIF(void);

void GetBlockWith1DAIF( int             ref_idx,
                       StorablePicture **list,
                       int             x_pos,
                       int             y_pos,
                       int             mvx,
                       int             mvy,
                       int             img_width,
                       int             img_height,
                       int             block[BLOCK_SIZE][BLOCK_SIZE]);

void GetBlockWith1DAIF_quarter_pel(int             ref_idx,
                                   StorablePicture **list,
                                   int             x_pos,
                                   int             y_pos,
                                   int             mvx,
                                   int             mvy,
                                   int             img_width,
                                   int             img_height,
                                   int             block[BLOCK_SIZE][BLOCK_SIZE]); 

void GetBlockWith1DAIF_16bits__quarter_pel(int             ref_idx,
                                           StorablePicture **list,
                                           int             x_pos,
                                           int             y_pos,
                                           int             mvx,
                                           int             mvy,
                                           int             img_width,
                                           int             img_height,
                                           int             block[BLOCK_SIZE][BLOCK_SIZE]); 
#endif


#ifdef E_DAIF
void CalculateFilterCoefficients_EDAIF(void);
void DequantizeFilterCoefficientsExt(int QFCoef[SQR_FILTER], int *numQBits, int QFOffsetI, int QFOffsetF, int QFOffsetSign, double FCoef[SQR_FILTER]);
int  FindPositionInt(int Dim, int cur_pos, int Offset, int mv); 
#endif
#ifdef EAIF
void GetBlockWithEAIF(int ref_idx, StorablePicture **list, int x_pos, int y_pos, int mvx, int mvy,
											 int img_width, int img_height, int block[BLOCK_SIZE][BLOCK_SIZE]);
void CalculateFilterCoefficients_EAIF(void);
#endif

#ifdef EDAIF2
void SetCalc2FilterIdx();
void CalculateFilterCoefficients_EDAIF2(void);
void ArrangeSymmetryEDAIF2(void);
void PrintFilterCoefFloatFile(double FCoef[15][SQR_FILTER]);
void DequantizeFilterOffsetEDAIF2(int QFOffsetI, int QFOffsetF, int QFOffsetSign, double FCoef[SQR_FILTER]);
void DequantizeFilterCoefficientsEDAIF2(int QFCoef[SQR_FILTER], double FCoef[SQR_FILTER], int sub_pos);
#endif
#endif
