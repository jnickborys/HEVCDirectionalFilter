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

#include "global.h"
#include "mbuffer.h"

#ifdef EDAIF2
#include "aif_common.h"
#endif

#define SQR(a) ((a)*(a))
#ifndef EDAIF2
#define DEFAULT_QUANT 9
#endif

#ifdef E_DAIF
#define FILTER_SIZE_INT   5
#define SQR_FILTER_INT    (FILTER_SIZE_INT * FILTER_SIZE_INT + 1)
#define SQR_FILTER        37
#else
#define SQR_FILTER        36
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
/*#define a_pos  0         //  1
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
#define o_pos 14         // 15*/

#define NUMBER_OF_DIF_COEFFS 54
// 54 =      3 for b_pos +
//           6 for a_pos +
//          18 for f_pos +
//          21 for e_pos +
//          6 for j_pos

/*
Ematrix*FilterCoef=CVector;
*/

#ifdef DIRECTIONAL_FILTER

#define FILTER_TYPE_2D_NS 1
#define FILTER_TYPE_2D_S  2
#define FILTER_TYPE_1D    3
#define FILTER_TYPE_EDAIF 4
#define FILTER_TYPE_EAIF  5

#define IMP_FLOAT32    0
#define IMP_INT16      1

#define MAX_NUM_SUBPELS_AIF 15
#define MAX_AIF_SUPPORT 36
void UnifiedOneForthPixWith_1DAIF_float(void);// Upsampling with 1D-AIF, floating point arithm. implementation
void UnifiedOneForthPixWith_1DAIF_int16(void);// Upsampling with 1D-AIF, floating point arithm. implementation
void initFilterCustom(int filterID);
void PredictFilterCoefficients1DAIF(void);

void RepresentCoefsIn16bits(int ImpType);
void ApplySumRestriction(int ImpType, int sub_pos);
#endif
#ifdef E_DAIF
void UnifiedOneForthPixWith_EDAIF(void);
#endif
#ifdef EAIF
void PredictFilterCoefficientsEAIF(void);
void UnifiedOneForthPixWith_EAIF(void);
#endif

#ifdef EDAIF2
typedef struct
{
  int g_aif_ready;
  int g_ImpType;
  int	UseAllSubpelPositions;									// 1 if FilterFlag for all independent positions is 1
  int SubpelPositionsPattern;
  int RealTimeDecomposition;
  int realSymCommands[4][MAX_NUM_AIF];
  int FilterFlag[MAX_NUM_AIF];													// Flags defining if Filter at the particular position calculated

  int Calc2Filt_Indexes[MAX_NUM_AIF][SQR_FILTER];
  int Filt_Indexes_offset[MAX_NUM_AIF];
  int FILTER_NEW_SUB_POS[MAX_NUM_AIF];
  double STANDARD_2D_FILTER[MAX_NUM_AIF][SQR_FILTER];
  int TwoDEquationPattern[MAX_NUM_AIF][SQR_FILTER];
  double CalculatedFilter[MAX_NUM_AIF][SQR_FILTER];							// temporal storage
  double FilterCoef[MAX_NUM_AIF][SQR_FILTER];									// estimated filter coefficients
  short int FilterCoef16bits[MAX_NUM_AIF][SQR_FILTER];	

  double DiffFilterCoef[MAX_NUM_AIF][SQR_FILTER];									
  int DiffQFilterCoef[MAX_NUM_AIF][SQR_FILTER];				// differences to be transmitted
  int POS_EQUATION_NUMBER[MAX_NUM_AIF];								// number of different coefficietns for each sub-pel position
  int SymmetryPosition[MAX_NUM_AIF];	// if 0, the position is copied from another one
  int Is1DPosition[MAX_NUM_AIF];// if 0, the position is a 2D one
  int IsDiagonal1D[MAX_NUM_AIF]; // if 1, the filter alighned NW-SE, 2 - NE-SW
  int nBitsIntRepresent[MAX_NUM_AIF];

  double NxN_Matrix[MAX_NUM_AIF][SQR_FILTER][SQR_FILTER];			// SQR(PEL)-1 X SQR(FILTER_SIZE) X SQR(FILTER_SIZE) - matrix
  double N_Vector[MAX_NUM_AIF][SQR_FILTER];										// SQR(PEL)-1 X SQR(FILTER_SIZE) - vector
}FilterEntity;


extern FilterEntity SetAIF[3];
extern double costMCP[4][MAX_NUM_AIF+1];

void PredictFilterCoefficientsDAIF2(int mode);
void ExtendDiagonalFilter_DAIF2(int sub_pos, double FCoef[MAX_NUM_AIF][SQR_FILTER]);
void PreProcessFilterCoef(void);
void initEDAIF2(void);
void SwapAIF(int filterID, int dir);
#endif
void InitAdaptiveFilter(InputParameters *input);
void FreeAdaptiveFilter(void);
int  FindPosition(int Dim, int cur_pos, int Offset, int mv);
void ReorderPosition(int* filter_pos, int* new_eq,  int* new_sub_pos,
                     int cur_xy, int equation, int sub_pos);
void FindFilterCoef(void);
void FindFilterCoefHor(void); // separable aif
void FindFilterCoefVer(void); // separable aif
void ExtendFilterCoefficientsDouble(int sub_pos, double FCoef[15][SQR_FILTER]);
void ExtendFilterCoefficientsInt(int sub_pos, int ICoef[15][SQR_FILTER]);
void QuantizeFilterCoefficients(double FCoef[SQR_FILTER], int QFCoef[SQR_FILTER]);
void PredictFilterCoefficients(void);
void PredictFilterCoefficientsHor(void);
void PredictFilterCoefficientsVer(void);
void ResetAdaptiveFilter(void);
int  UseAdaptiveFilterForCurrentFrame(void);
void Position1DFilter(int sub_pos, double FCoef[SQR_FILTER]);
void SetAdaptiveFilter(void);
void SetAdaptiveFilterHor(void); // separable aif
void SetAdaptiveFilterVer(void); // separable aif
void SetAdaptiveFilterEighthpel(void);
void SetAdaptiveFilterEighthpelHor(void); // separable aif
void SetAdaptiveFilterEighthpelVer(void); // separable aif
void SetMotionVectorsAIF(void);
void PrintFilterCoefDouble(double FCoef[15][SQR_FILTER]);
void PrintFilterCoefInt(int FCoef[15][SQR_FILTER]);
void PrintFilterCoefDoubleSep (double FCoef[15][SQR_FILTER]); 
void UnifiedOneForthPixWithNewFilter(void);
void UnifiedOneForthPixWithNewFilterHor (void); // separable aif
void UnifiedOneForthPixWithNewFilterVer (void); // separable aif
void SwapUpsampledFrames(void);
void FindSAD(void);
int getDecisionOnAIF_MCP(void);

#ifdef E_DAIF
void AIFDecision(int *FilterFlagInt, int *FilterFlag, int *SymmetryPosition);
void PredictFilterCoefficients_EDAIF(void);
int  FindPositionInt(int Dim, int cur_pos, int Offset, int mv);
void QuantizeIntegerFilter();
void QuantizeFilterOffsets();
#endif  // E_DAIF

#ifdef EDAIF2
void FindFilterCoef2(void);
void GetErrorMCP(int *FilterFlagInt,int filterID);
void AIFDecision2(int *FilterFlagInt, int *FilterFlag, int *SymmetryPosition, int *FILTER_NEW_SUB_POS);
void PrintFilterCoefFloatFile(double FCoef[15][SQR_FILTER]);
void QuantizeFilterCoefficientsEDAIF2(double FCoef[SQR_FILTER], int QFCoef[SQR_FILTER], int sub_pos);
void SetCalc2FilterIdx();
#endif
#endif

