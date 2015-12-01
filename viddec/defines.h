
/*!
**************************************************************************
* \file defines.h
*
* \brief
*    Headerfile containing some useful global definitions
*
* \author
*    Detlev Marpe  
*    Copyright (C) 2000 HEINRICH HERTZ INSTITUTE All Rights Reserved.
*
* \date
*    21. March 2001
**************************************************************************
*/

#ifndef _DEFINES_H_
#define _DEFINES_H_

#if defined _DEBUG
#define TRACE           0                   //!< 0:Trace off 1:Trace on 2:detailed CABAC context information
#else
#define TRACE           0                   //!< 0:Trace off 1:Trace on 2:detailed CABAC context information
#endif

//*** KTA Adoption techniques ********************************************************

// Adaptive Filter
#define ADAPTIVE_FILTER

#ifdef ADAPTIVE_FILTER
#define DIRECTIONAL_FILTER
#define E_DAIF      // Enhanced DAIF
#define EAIF     
#define EDAIF2  // DAIF scheme with variable symmetry and filter structure
#endif  // ADAPTIVE_FILTER

// 1/8 pel motion compensation
#define EIGHTH_PEL
// Adaptive Prediction Error Coding in Frequency and Spatial Domain
#define ADAPTIVE_FD_SD_CODING
// mv competition predictor
#define MV_COMPETITION
// Adaptive Quanziation Matirx Selection (IAQMS)
#define ADAPTIVE_QUANTIZATION
// Directional transform for intra coding
#define USE_INTRA_MDDT      
// use high precision interpolation and prediction (including adaptive bidirectional prediction rounding)
#define USE_HP_FILTER

#ifdef USE_HP_FILTER
// Use switched filters
#define SWITCHED_FILTERS 
#define HPF_COMBO
#endif
#define MB32X32
#ifdef MB32X32
#ifdef MV_COMPETITION
#define MB32X32_MVC
#endif
#define MAX_MB_EXT_LEVEL 2 
#endif




// Quadtree-based Adaptive Loop Filter (QALF)
#define ADAPTIVE_LOOP_FILTER

#ifdef ADAPTIVE_LOOP_FILTER
#define ALF_SPATIAL_PREDICT_COEF  // C129
#endif // ADAPTIVE_LOOP_FILTER

//*** KTA Adoption techniques as reference tools *************************************

// Internal bit-depth increase
#define INTERNAL_BIT_DEPTH_INCREASE

// Bug fix in JM11.0
#define BUG_FIX_FOR_FREXT

//************************************************************************************

//#define PAIR_FIELDS_IN_OUTPUT

//#define MAX_NUM_SLICES 150
#define MAX_NUM_SLICES 50

//FREXT Profile IDC definitions
#define FREXT_HP        100      //!< YUV 4:2:0/8 "High"
#define FREXT_Hi10P     110      //!< YUV 4:2:0/10 "High 10"
#define FREXT_Hi422     122      //!< YUV 4:2:2/10 "High 4:2:2"
#define FREXT_Hi444     144      //!< YUV 4:4:4/12 "High 4:4:4"

#define YUV400 0
#define YUV420 1
#define YUV422 2
#define YUV444 3


#define ZEROSNR 0

// CAVLC
#define LUMA              0
#define LUMA_INTRA16x16DC 1
#define LUMA_INTRA16x16AC 2

#define TOTRUN_NUM    15
#define RUNBEFORE_NUM  7


//--- block types for CABAC ----
#define LUMA_16DC       0
#define LUMA_16AC       1
#define LUMA_8x8        2
#define LUMA_8x4        3
#define LUMA_4x8        4
#define LUMA_4x4        5
#define CHROMA_DC       6
#define CHROMA_AC       7
#define CHROMA_DC_2x4   8
#define CHROMA_DC_4x4   9
#ifdef USE_INTRA_MDDT
#define LUMA_16x16      10
#ifdef MB32X32
#define LUMA_16x16P     11
#define LUMA_16x8P      12
#define LUMA_8x16P      13
#define NUM_BLOCK_TYPES 14
#else
#define NUM_BLOCK_TYPES 11
#endif
#else
#define NUM_BLOCK_TYPES 10
#endif 


#define MAX_CODED_FRAME_SIZE 8000000         //!< bytes for one frame

//#define _LEAKYBUCKET_

#ifdef ADAPTIVE_FILTER 
enum {
  PSKIP        =  0,
  BSKIP_DIRECT =  0,
  P16x16       =  1,
  P16x8        =  2,
  P8x16        =  3,
  SMB8x8       =  4,
  SMB8x4       =  5,
  SMB4x8       =  6,
  SMB4x4       =  7,
  P8x8         =  8,
  I4MB         =  9,
  I16MB        = 10,
  IBLOCK       = 11,
  SI4MB        = 12,
  I8MB         = 13,
  IPCM         = 14,
  MAXMODE      = 15
} MBModeTypes;
#endif

#define absm(A) ((A)<(0) ? (-(A)):(A))      //!< abs macro, faster than procedure

#define Clip1(a)            ((a)>img->max_imgpel_value?img->max_imgpel_value:((a)<0?0:(a)))
#define Clip1_Chr(a)        ((a)>img->max_imgpel_value_uv?img->max_imgpel_value_uv:((a)<0?0:(a)))
#define Clip3(min,max,val) (((val)<(min))?(min):(((val)>(max))?(max):(val)))

#define P8x8    8
#define I4MB    9
#define I16MB   10
#define IBLOCK  11
#define SI4MB   12
#define I8MB    13
#define IPCM    14
#define MAXMODE 15

#define IS_INTRA(MB)    ((MB)->mb_type==I4MB  || (MB)->mb_type==I16MB ||(MB)->mb_type==IPCM || (MB)->mb_type==I8MB || (MB)->mb_type==SI4MB)
#define IS_NEWINTRA(MB) ((MB)->mb_type==I16MB  || (MB)->mb_type==IPCM)
#define IS_OLDINTRA(MB) ((MB)->mb_type==I4MB)

#define IS_INTER(MB)    ((MB)->mb_type!=I4MB  && (MB)->mb_type!=I16MB && (MB)->mb_type!=I8MB  && (MB)->mb_type!=IPCM)
#define IS_INTERMV(MB)  ((MB)->mb_type!=I4MB  && (MB)->mb_type!=I16MB && (MB)->mb_type!=I8MB  && (MB)->mb_type!=0 && (MB)->mb_type!=IPCM)
#define IS_DIRECT(MB)   ((MB)->mb_type==0     && (img->type==B_SLICE ))
#define IS_COPY(MB)     ((MB)->mb_type==0     && (img->type==P_SLICE || img->type==SP_SLICE))
#define IS_P8x8(MB)     ((MB)->mb_type==P8x8)

#ifdef  INTERNAL_BIT_DEPTH_INCREASE
#define SQR_DEPTH(rec,ref,depth,increase) ((increase)? img->quad[abs((ref)-Clip3(0,((1<<(depth))-1),((rec)+(1<<(increase-1)))>>increase))] : img->quad[abs((ref)-(rec))])
#endif

// Quantization parameter range

#define MIN_QP          0
#define MAX_QP          51

#define BLOCK_SIZE      4
#define MB_BLOCK_SIZE   16


#define NO_INTRA_PMODE  9        //!< #intra prediction modes
/* 4x4 intra prediction modes */
#define VERT_PRED             0
#define HOR_PRED              1
#define DC_PRED               2
#define DIAG_DOWN_LEFT_PRED   3
#define DIAG_DOWN_RIGHT_PRED  4
#define VERT_RIGHT_PRED       5
#define HOR_DOWN_PRED         6
#define VERT_LEFT_PRED        7
#define HOR_UP_PRED           8

// 16x16 intra prediction modes
#define VERT_PRED_16    0
#define HOR_PRED_16     1
#define DC_PRED_16      2
#define PLANE_16        3

// 8x8 chroma intra prediction modes
#define DC_PRED_8       0
#define HOR_PRED_8      1
#define VERT_PRED_8     2
#define PLANE_8         3

#define EOS             1                       //!< End Of Sequence
#define SOP             2                       //!< Start Of Picture
#define SOS             3                       //!< Start Of Slice

#define DECODING_OK     0
#define SEARCH_SYNC     1
#define PICTURE_DECODED 2

#define MAX_REFERENCE_PICTURES 32               //!< H264 allows 32 fields

#define INVALIDINDEX  (-135792468)

#if !defined(WIN32) || defined(__GNUC__)
#define max(a, b)      ((a) > (b) ? (a) : (b))  //!< Macro returning max value
#define min(a, b)      ((a) < (b) ? (a) : (b))  //!< Macro returning min value
#endif


#define MVPRED_MEDIAN   0
#define MVPRED_L        1
#define MVPRED_U        2
#define MVPRED_UR       3

#define DECODE_COPY_MB  0
#define DECODE_MB       1
//#define DECODE_MB_BFRAME 2

#define BLOCK_MULTIPLE      (MB_BLOCK_SIZE/BLOCK_SIZE)

//Start code and Emulation Prevention need this to be defined in identical manner at encoder and decoder
#define ZEROBYTES_SHORTSTARTCODE 2 //indicates the number of zero bytes in the short start-code prefix

#ifdef ADAPTIVE_QUANTIZATION
#define NUM_OF_IAQMS_MAT 4
#endif

#endif

