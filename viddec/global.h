
/*!
************************************************************************
*  \file
*     global.h
*  \brief
*     global definitions for for H.264 decoder.
*  \author
*     Copyright (C) 1999  Telenor Satellite Services,Norway
*                         Ericsson Radio Systems, Sweden
*
*     Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
*
*     Telenor Satellite Services
*     Keysers gt.13                       tel.:   +47 23 13 86 98
*     N-0130 Oslo,Norway                  fax.:   +47 22 77 79 80
*
*     Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
*
*     Ericsson Radio Systems
*     KI/ERA/T/VV
*     164 80 Stockholm, Sweden
*
************************************************************************
*/
#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <stdio.h>                              //!< for FILE
#include <time.h>
#include <sys/timeb.h>
#include "defines.h"
#include "parsetcommon.h"

// Tool defines (AIF, 1/8pel, MVC) is moved to define.h by AT 2006/11/01

#ifdef MV_COMPETITION
#define MAX_MV_PREDICTOR  8
#endif

#ifdef WIN32
#define  snprintf _snprintf
#define  open     _open
#define  close    _close
#define  read     _read
#define  write    _write
#define  lseek    _lseeki64
#define  fsync    _commit
#define  OPENFLAGS_WRITE _O_WRONLY|_O_CREAT|_O_BINARY|_O_TRUNC
#define  OPEN_PERMISSIONS _S_IREAD | _S_IWRITE
#define  OPENFLAGS_READ  _O_RDONLY|_O_BINARY
#else
#define  OPENFLAGS_WRITE O_WRONLY|O_CREAT|O_TRUNC
#define  OPENFLAGS_READ  O_RDONLY
#define  OPEN_PERMISSIONS S_IRUSR | S_IWUSR
#endif


typedef unsigned char   byte;                   //!<  8 bit unsigned
typedef int             int32;
typedef unsigned int    u_int32;

#define imgpel unsigned short

#if defined(WIN32) && !defined(__GNUC__)
typedef __int64   int64;
# define FORMAT_OFF_T "I64d"
#ifndef INT64_MIN
# define INT64_MIN        (-9223372036854775807i64 - 1i64)
#endif
#else
typedef long long int64;
# define FORMAT_OFF_T "lld"
#ifndef INT64_MIN
# define INT64_MIN        (-9223372036854775807LL - 1LL)
#endif
#endif





pic_parameter_set_rbsp_t *active_pps;
seq_parameter_set_rbsp_t *active_sps;


// global picture format dependend buffers, mem allocation in decod.c ******************
int  **refFrArr;                                //!< Array for reference frames of each block

imgpel **imgY_ref;                                //!< reference frame find snr
imgpel ***imgUV_ref;

int  ReMapRef[20];
// B pictures
int  Bframe_ctr;
int  frame_no;

int  g_nFrame;

// For MB level frame/field coding
int  TopFieldForSkip_Y[16][16];
int  TopFieldForSkip_UV[2][16][16];

int  InvLevelScale4x4Luma_Intra[6][4][4];
int  InvLevelScale4x4Chroma_Intra[2][6][4][4];

int  InvLevelScale4x4Luma_Inter[6][4][4];
int  InvLevelScale4x4Chroma_Inter[2][6][4][4];

int  InvLevelScale8x8Luma_Intra[6][8][8];

int  InvLevelScale8x8Luma_Inter[6][8][8];

int  *qmatrix[8];

#ifdef ADAPTIVE_QUANTIZATION
int  *ScalingList4x4[6], *ScalingList8x8[2];
int  *ScalingList4x4_IAQMS[NUM_OF_IAQMS_MAT][6], *ScalingList8x8_IAQMS[NUM_OF_IAQMS_MAT][6];
int  InvLevelScale4x4Luma_Intra_IAQMS[NUM_OF_IAQMS_MAT][6][4][4];
int  InvLevelScale4x4Chroma_Intra_IAQMS[NUM_OF_IAQMS_MAT][2][6][4][4];
int  InvLevelScale4x4Luma_Inter_IAQMS[NUM_OF_IAQMS_MAT][6][4][4];
int  InvLevelScale4x4Chroma_Inter_IAQMS[NUM_OF_IAQMS_MAT][2][6][4][4];
int  InvLevelScale8x8Luma_Intra_IAQMS[NUM_OF_IAQMS_MAT][6][8][8];
int  InvLevelScale8x8Luma_Inter_IAQMS[NUM_OF_IAQMS_MAT][6][8][8];
int  *qmatrixIAQMS[NUM_OF_IAQMS_MAT][12];
#endif

#ifdef ADAPTIVE_FILTER
double **tmp_coef; 
#endif

#define ET_SIZE 300      //!< size of error text buffer
char errortext[ET_SIZE]; //!< buffer for error message for exit with error()

                         /***********************************************************************
                         * T y p e    d e f i n i t i o n s    f o r    T M L
                         ***********************************************************************
*/

//! Data Partitioning Modes
typedef enum
{
  PAR_DP_1,    //!< no data partitioning is supported
    PAR_DP_3,    //!< data partitioning with 3 partitions
} PAR_DP_TYPE;


//! Output File Types
typedef enum
{
  PAR_OF_ANNEXB,   //!< Current TML description
    PAR_OF_RTP,   //!< RTP Packet Output format
    //  PAR_OF_IFF    //!< Interim File Format
} PAR_OF_TYPE;

//! Boolean Type
/*typedef enum {
FALSE,
TRUE
} Boolean;
*/
//! definition of H.264 syntax elements
typedef enum {
  SE_HEADER,
    SE_PTYPE,
    SE_MBTYPE,
    SE_REFFRAME,
    SE_INTRAPREDMODE,
    SE_MVD,
#ifdef MV_COMPETITION
    SE_MV_PREDICTOR,
#endif
    SE_CBP_INTRA,
    SE_LUM_DC_INTRA,
    SE_CHR_DC_INTRA,
    SE_LUM_AC_INTRA,
    SE_CHR_AC_INTRA,
    SE_CBP_INTER,
    SE_LUM_DC_INTER,
    SE_CHR_DC_INTER,
    SE_LUM_AC_INTER,
    SE_CHR_AC_INTER,
    SE_DELTA_QUANT_INTER,
    SE_DELTA_QUANT_INTRA,
    SE_BFRAME,
    SE_EOS,
    SE_MAX_ELEMENTS //!< number of maximum syntax elements, this MUST be the last one!
} SE_type;        // substituting the definitions in element.h


typedef enum {
  INTER_MB,
    INTRA_MB_4x4,
    INTRA_MB_16x16
} IntraInterDecision;

typedef enum {
  BITS_TOTAL_MB,
    BITS_HEADER_MB,
    BITS_INTER_MB,
    BITS_CBP_MB,
    BITS_COEFF_Y_MB,
    BITS_COEFF_UV_MB,
    MAX_BITCOUNTER_MB
} BitCountType;

typedef enum {
  NO_SLICES,
    FIXED_MB,
    FIXED_RATE,
    CALLBACK,
    FMO
} SliceMode;


typedef enum {
  UVLC,
    CABAC
} SymbolMode;

typedef enum {
  LIST_0=0,
    LIST_1=1
} Lists;


typedef enum {
  FRAME,
    TOP_FIELD,
    BOTTOM_FIELD
} PictureStructure;           //!< New enum for field processing


typedef enum {
  P_SLICE = 0,
    B_SLICE,
    I_SLICE,
    SP_SLICE,
    SI_SLICE
} SliceType;

/***********************************************************************
* D a t a    t y p e s   f o r  C A B A C
***********************************************************************
*/

//! struct to characterize the state of the arithmetic coding engine
typedef struct
{
  unsigned int    Dlow, Drange;
  unsigned int    Dvalue;
  unsigned int    Dbuffer;
  int             Dbits_to_go;
  byte            *Dcodestrm;
  int             *Dcodestrm_len;
} DecodingEnvironment;

typedef DecodingEnvironment *DecodingEnvironmentPtr;

//! struct for context management
typedef struct
{
  unsigned short state;         // index into state-table CP  
  unsigned char  MPS;           // Least Probable Symbol 0/1 CP
} BiContextType;

typedef BiContextType *BiContextTypePtr;


/**********************************************************************
* C O N T E X T S   F O R   T M L   S Y N T A X   E L E M E N T S
**********************************************************************
*/

#define NUM_MB_TYPE_CTX  11
#define NUM_B8_TYPE_CTX  9
#define NUM_MV_RES_CTX   10
#define NUM_REF_NO_CTX   6
#define NUM_DELTA_QP_CTX 4
#define NUM_MB_AFF_CTX 4
#define NUM_TRANSFORM_SIZE_CTX 3
#ifdef ADAPTIVE_QUANTIZATION
#define NUM_MODULATED_QUANTIZATION_CTX 2
#endif
#ifdef ADAPTIVE_LOOP_FILTER
#define NUM_ALF_BLK_CTX  3
#endif

typedef struct
{
  BiContextType mb_type_contexts [4][NUM_MB_TYPE_CTX];
  BiContextType b8_type_contexts [2][NUM_B8_TYPE_CTX];
  BiContextType mv_res_contexts  [2][NUM_MV_RES_CTX];
#ifdef MV_COMPETITION
  BiContextType mv_predictor_skip_contexts  [MAX_MV_PREDICTOR];
  BiContextType mv_predictor_mvp_contexts   [MAX_MV_PREDICTOR];
  BiContextType mv_predictor_mvb_contexts   [MAX_MV_PREDICTOR];
#endif
  BiContextType ref_no_contexts  [2][NUM_REF_NO_CTX];
  BiContextType delta_qp_contexts[NUM_DELTA_QP_CTX];
  BiContextType mb_aff_contexts  [NUM_MB_AFF_CTX];
  BiContextType transform_size_contexts [NUM_TRANSFORM_SIZE_CTX];
#ifdef ADAPTIVE_QUANTIZATION
  BiContextType modulated_quantization_contexts [NUM_MODULATED_QUANTIZATION_CTX];
#endif
#ifdef ADAPTIVE_LOOP_FILTER
  BiContextType alf_blk_flag_contexts   [NUM_ALF_BLK_CTX];
#endif
} MotionInfoContexts;

#define NUM_IPR_CTX    2
#define NUM_CIPR_CTX   4
#define NUM_CBP_CTX    4
#define NUM_BCBP_CTX   4
#define NUM_MAP_CTX   15
#define NUM_LAST_CTX  15
#define NUM_ONE_CTX    5
#define NUM_ABS_CTX    5

#ifdef ADAPTIVE_FD_SD_CODING
#define NUM_SCAN_CTX        4
#define NUM_MB_adap_CTX     4
#define NUM_BCBP_FD_SD_CTX  4
#define NUM_MAP8_CTX       63
#define NUM_LAST8_CTX      63
#endif

typedef struct
{
#ifdef ADAPTIVE_FD_SD_CODING
  BiContextType  MB_adaptive_SD_context[NUM_MB_adap_CTX];
  BiContextType  bcbp_contexts_FD_SD   [NUM_BCBP_FD_SD_CTX];
  BiContextType  bcbp8_contexts_FD_SD  [NUM_BCBP_FD_SD_CTX];
  BiContextType  one_contexts_SD       [NUM_ONE_CTX];
  BiContextType  one8_contexts_SD      [NUM_ONE_CTX];
  BiContextType  abs_contexts_SD       [NUM_ABS_CTX];
  BiContextType  abs8_contexts_SD      [NUM_ABS_CTX];
  BiContextType  map_contexts_SD       [NUM_MAP_CTX];
  BiContextType  map8_contexts_SD      [NUM_MAP8_CTX];
  BiContextType  last_contexts_SD      [NUM_LAST_CTX];
  BiContextType  last8_contexts_SD     [NUM_LAST8_CTX];
#endif
  BiContextType  ipr_contexts [NUM_IPR_CTX];
  BiContextType  cipr_contexts[NUM_CIPR_CTX]; 
  BiContextType  cbp_contexts [3][NUM_CBP_CTX];
  BiContextType  bcbp_contexts[NUM_BLOCK_TYPES][NUM_BCBP_CTX];
  BiContextType  map_contexts [NUM_BLOCK_TYPES][NUM_MAP_CTX];
  BiContextType  last_contexts[NUM_BLOCK_TYPES][NUM_LAST_CTX];
  BiContextType  one_contexts [NUM_BLOCK_TYPES][NUM_ONE_CTX];
  BiContextType  abs_contexts [NUM_BLOCK_TYPES][NUM_ABS_CTX];
  BiContextType  fld_map_contexts [NUM_BLOCK_TYPES][NUM_MAP_CTX];
  BiContextType  fld_last_contexts[NUM_BLOCK_TYPES][NUM_LAST_CTX];
} TextureInfoContexts;


//*********************** end of data type definition for CABAC *******************

/***********************************************************************
* N e w   D a t a    t y p e s   f o r    T M L
***********************************************************************
*/

struct img_par;
struct inp_par;
struct stat_par;

/*! Buffer structure for decoded referenc picture marking commands */
typedef struct DecRefPicMarking_s
{
  int memory_management_control_operation;
  int difference_of_pic_nums_minus1;
  int long_term_pic_num;
  int long_term_frame_idx;
  int max_long_term_frame_idx_plus1;
  struct DecRefPicMarking_s *Next;
} DecRefPicMarking_t;

//! Syntaxelement
typedef struct syntaxelement
{
  int           type;                  //!< type of syntax element for data part.
  int           value1;                //!< numerical value of syntax element
  int           value2;                //!< for blocked symbols, e.g. run/level
#ifdef MB32X32
  int           mb_ext_level;
#endif
  int           len;                   //!< length of code
  int           inf;                   //!< info part of UVLC code
  unsigned int  bitpattern;            //!< UVLC bitpattern
  int           context;               //!< CABAC context
  int           k;                     //!< CABAC context for coeff_count,uv
#ifdef ADAPTIVE_FD_SD_CODING
  int           writing_b4;
  int           writing_b8;
#endif
#if TRACE
#define       TRACESTRING_SIZE 100           //!< size of trace string
  char          tracestring[TRACESTRING_SIZE]; //!< trace string
#endif
  
  //! for mapping of UVLC to syntaxElement
  void    (*mapping)(int len, int info, int *value1, int *value2);
  //! used for CABAC: refers to actual coding method of each individual syntax element type
  void  (*reading)(struct syntaxelement *, struct inp_par *, struct img_par *, DecodingEnvironmentPtr);
  
} SyntaxElement;

//! Macroblock
typedef struct macroblock
{
  int           qp;                  //!< QP luma
  int           qpc[2];              //!< QP chroma
  
  int           slice_nr;
  int           delta_quant;          //!< for rate control
  
  struct macroblock   *mb_available_up;   //!< pointer to neighboring MB (CABAC)
  struct macroblock   *mb_available_left; //!< pointer to neighboring MB (CABAC)
  
  // some storage of macroblock syntax elements for global access
  int           mb_type;
  int           mvd[2][BLOCK_MULTIPLE][BLOCK_MULTIPLE][2];      //!< indices correspond to [forw,backw][block_y][block_x][x,y]
  int           cbp;
  int64         cbp_blk ;
  int64         cbp_bits;
#ifdef ADAPTIVE_FD_SD_CODING
  int64         FD_or_SD_bits;
  int           written_SD_Coding_on_off;
  int           SD_Coding_on_off;
  int           SD_or_FD           [2][2];    //y,x
  int           SD_or_FD_t8x8;                //y,x
  int           quantizer_indices  [16][16];  //y,x
#endif  
  int           is_skip;
#ifdef ADAPTIVE_QUANTIZATION
  int           mb_iaqms_idx;
#endif
  int           i16mode;
  int           b8mode[4];
  int           b8pdir[4];
  int           ei_flag;
  
  int           LFDisableIdc;
  int           LFAlphaC0Offset;
  int           LFBetaOffset;
  
  int           c_ipred_mode;       //!< chroma intra prediction mode
  int           mb_field;
  
  int           skip_flag;
  
  int mbAddrA, mbAddrB, mbAddrC, mbAddrD;
  int mbAvailA, mbAvailB, mbAvailC, mbAvailD;
  
  int           luma_transform_size_8x8_flag;
  int           NoMbPartLessThan8x8Flag;
} Macroblock;

//! Bitstream
typedef struct
{
  // CABAC Decoding
  int           read_len;           //!< actual position in the codebuffer, CABAC only
  int           code_len;           //!< overall codebuffer length, CABAC only
  // UVLC Decoding
  int           frame_bitoffset;    //!< actual position in the codebuffer, bit-oriented, UVLC only
  int           bitstream_length;   //!< over codebuffer lnegth, byte oriented, UVLC only
  // ErrorConcealment
  byte          *streamBuffer;      //!< actual codebuffer for read bytes
  int           ei_flag;            //!< error indication, 0: no error, else unspecified error
} Bitstream;

//! DataPartition
typedef struct datapartition
{
  
  Bitstream           *bitstream;
  DecodingEnvironment de_cabac;
  
  int     (*readSyntaxElement)(SyntaxElement *, struct img_par *, struct inp_par *, struct datapartition *);
  /*!< virtual function;
  actual method depends on chosen data partition and
  entropy coding method  */
} DataPartition;

//! Slice
typedef struct
{
  int                 ei_flag;       //!< 0 if the partArr[0] contains valid information
  int                 qp;
  int                 slice_qp_delta;
  int                 picture_type;  //!< picture type
  PictureStructure    structure;     //!< Identify picture structure type
  int                 start_mb_nr;   //!< MUST be set by NAL even in case of ei_flag == 1
  int                 max_part_nr;
  int                 dp_mode;       //!< data partioning mode
  int                 next_header;
  //  int                 last_mb_nr;    //!< only valid when entropy coding == CABAC
  DataPartition       *partArr;      //!< array of partitions
  MotionInfoContexts  *mot_ctx;      //!< pointer to struct of context models for use in CABAC
  TextureInfoContexts *tex_ctx;      //!< pointer to struct of context models for use in CABAC
  
  int                 ref_pic_list_reordering_flag_l0;
  int                 *reordering_of_pic_nums_idc_l0;
  int                 *abs_diff_pic_num_minus1_l0;
  int                 *long_term_pic_idx_l0;
  int                 ref_pic_list_reordering_flag_l1;
  int                 *reordering_of_pic_nums_idc_l1;
  int                 *abs_diff_pic_num_minus1_l1;
  int                 *long_term_pic_idx_l1;
  
  int     (*readSlice)(struct img_par *, struct inp_par *);
  
  int                 LFDisableIdc;     //!< Disable loop filter on slice
  int                 LFAlphaC0Offset;  //!< Alpha and C0 offset for filtering slice
  int                 LFBetaOffset;     //!< Beta offset for filtering slice
  
  int                 pic_parameter_set_id;   //!<the ID of the picture parameter set the slice is reffering to
} Slice;

#ifdef ADAPTIVE_LOOP_FILTER
typedef struct
{
  int          tap;
  int          num_coeff;
  int          *coeff;
  int          tap_chroma;
  int          num_coeff_chroma;
  int          *coeff_chroma;
} ALFParam;

typedef struct alfblock
{
  Boolean          alf_blk_flag;
  struct alfblock  *alfb_available_up;   //!< pointer to neighboring Block (CABAC)
  struct alfblock  *alfb_available_left; //!< pointer to neighboring Block (CABAC)
  int              alfbAddrA, alfbAddrB, alfbAddrC, alfbAddrD;
  int              alfbAvailA, alfbAvailB, alfbAvailC, alfbAvailD;
} AlfBlock;

typedef struct
{
  int          alf_flag;
#ifdef ALF_SPATIAL_PREDICT_COEF
  int          pred_coef_mode;
#endif
  int          alf_chroma_idc;
  ALFParam     alf;
  int          block_control_flag;
  int          alf_block_size_idx;
  int          alf_block_size;
  int          WidthInAlfBlks;
  int          HeightInAlfBlks;
  int          num_alf_block;
  AlfBlock     *alf_blk_data;
  int          qt_partition_flag;
  int          max_layer_level;
  int          parent_block_size;
  int          WidthInParentBlks;
  int          HeightInParentBlks;
  int          num_parent_block;
} alf_parameter_set_rbsp_t;

alf_parameter_set_rbsp_t *active_alfps;
#endif // ADAPTIVE_LOOP_FILTER

//****************************** ~DM ***********************************

// image parameters
typedef struct img_par
{
  int number;                                 //!< frame number
  unsigned int current_mb_nr; // bitstream order
  unsigned int num_dec_mb;
  int current_slice_nr;
  int *intra_block;
  int tr;                                     //!< temporal reference, 8 bit, wrapps at 255
  int qp;                                     //!< quant for the current frame
  int qpsp;                                   //!< quant for SP-picture predicted frame
  int sp_switch;                              //!< 1 for switching sp, 0 for normal sp
  int direct_spatial_mv_pred_flag;            //!< 1 for Spatial Direct, 0 for Temporal
  int type;                                   //!< image type INTER/INTRA
  int width;
  int height;
  int width_cr;                               //!< width chroma
  int height_cr;                              //!< height chroma
  int mb_y;
  int mb_x;
  int block_y;
  int pix_y;
  int pix_x;
  int pix_c_y;
  int block_x;
  int pix_c_x;
  
  int allrefzero;
  imgpel mpr[16][16];                         //!< predicted block
  int mvscale[6][MAX_REFERENCE_PICTURES];
  
  // <FTRD : Compatibility with hierarchical B frames
#ifdef MV_COMPETITION
  int mvscale_hb[2][2][MAX_REFERENCE_PICTURES][MAX_REFERENCE_PICTURES]; 
#endif 
  // FTRD>
  
  int m7[16][16];                             //!< final 4x4 block. Extended to 16x16 for ABT
  int cof[4][12][4][4];                       //!< correction coefficients from predicted   
  int cofu[16];                                                                             
  int **ipredmode;                            //!< prediction type [90][74]
  int *quad;
  int ***nz_coeff;
  int **siblock;
  int cod_counter;                            //!< Current count of number of skipped macroblocks in a row
  
  int newframe;
  
  int structure;                               //!< Identify picture structure type
  int pstruct_next_P;
  
  // B pictures
  Slice       *currentSlice;                   //!< pointer to current Slice data struct
  Macroblock          *mb_data;                //!< array containing all MBs of a whole frame
  int subblock_x;
  int subblock_y;
  int is_intra_block;
  int is_v_block;
  
  // For MB level frame/field coding
  int MbaffFrameFlag;
  
  // for signalling to the neighbour logic that this is a deblocker call
  int DeblockCall;
  
  int **field_anchor;
#ifdef EIGHTH_PEL
  int mv_res;     //!< motion vector resolution 0 = 1/4-pel;  1 = 1/8-pel
#endif
#ifdef ADAPTIVE_FD_SD_CODING
  int Allow_SD_Coding;
  int SD_Quantizer;
#endif
  DecRefPicMarking_t *dec_ref_pic_marking_buffer;                    //!< stores the memory management control operations
  
  int num_ref_idx_l0_active;             //!< number of forward reference
  int num_ref_idx_l1_active;             //!< number of backward reference
  
  int slice_group_change_cycle;
  
  int redundant_pic_cnt; 
  
  int explicit_B_prediction;
  
  unsigned int pre_frame_num;           //!< store the frame_num in the last decoded slice. For detecting gap in frame_num.
  
  // End JVT-D101
  // POC200301: from unsigned int to int
  int toppoc;      //poc for this top field // POC200301
  int bottompoc;   //poc of bottom field of frame
  int framepoc;    //poc of this frame // POC200301
  unsigned int frame_num;   //frame_num for this frame
  unsigned int field_pic_flag;
  unsigned int bottom_field_flag;
  
  //the following is for slice header syntax elements of poc
  // for poc mode 0.
  unsigned int pic_order_cnt_lsb;
  int delta_pic_order_cnt_bottom;
  // for poc mode 1.
  int delta_pic_order_cnt[3];
  
  // ////////////////////////
  // for POC mode 0:
  signed int PrevPicOrderCntMsb;
  unsigned int PrevPicOrderCntLsb;
  signed int PicOrderCntMsb;
  
  // for POC mode 1:
  unsigned int AbsFrameNum;
  signed int ExpectedPicOrderCnt, PicOrderCntCycleCnt, FrameNumInPicOrderCntCycle;
  unsigned int PreviousFrameNum, FrameNumOffset;
  int ExpectedDeltaPerPicOrderCntCycle;
  int PreviousPOC, ThisPOC;
  int PreviousFrameNumOffset;
  // /////////////////////////
  
  //weighted prediction
  unsigned int luma_log2_weight_denom;
  unsigned int chroma_log2_weight_denom;
  int ***wp_weight;  // weight in [list][index][component] order
  int ***wp_offset;  // offset in [list][index][component] order
  int ****wbp_weight; //weight in [list][fw_index][bw_index][component] order
  int wp_round_luma;
  int wp_round_chroma;
  unsigned int apply_weights;
  
  int idr_flag;
  int nal_reference_idc;                       //!< nal_reference_idc from NAL unit
#ifdef  USE_HP_FILTER
  int bipred_rounding_control;
#endif

  int idr_pic_id;
  
  int MaxFrameNum;
  
  unsigned int PicWidthInMbs;
  unsigned int PicHeightInMapUnits;
  unsigned int FrameHeightInMbs;
  unsigned int PicHeightInMbs;
  unsigned int PicSizeInMbs;
  unsigned int FrameSizeInMbs;
  unsigned int oldFrameSizeInMbs;
  
  int no_output_of_prior_pics_flag;
  int long_term_reference_flag;
  int adaptive_ref_pic_buffering_flag;
  
  int last_has_mmco_5;
  int last_pic_bottom_field;
  
  int model_number;
  
  // Fidelity Range Extensions Stuff
  int pic_unit_bitsize_on_disk;
  int bitdepth_luma;
  int bitdepth_chroma;
  int bitdepth_luma_qp_scale;
  int bitdepth_chroma_qp_scale;
  unsigned int dc_pred_value_luma;            //!< luma value for DC prediction (depends on luma pel bit depth)
  unsigned int dc_pred_value_chroma;          //!< chroma value for DC prediction (depends on chroma pel bit depth)
  int max_imgpel_value;                       //!< max value that one luma picture element (pixel) can take (depends on pic_unit_bitdepth)
  int max_imgpel_value_uv;                    //!< max value that one chroma picture element (pixel) can take (depends on pic_unit_bitdepth)
  int Transform8x8Mode;        
  int profile_idc;              
  int yuv_format;
  int lossless_qpprime_flag;
  int num_blk8x8_uv;
  int num_cdc_coeff;
  int mb_cr_size_x;
  int mb_cr_size_y;
  
  int idr_psnr_number;
  int psnr_number;
  
  time_t ltime_start;               // for time measurement
  time_t ltime_end;                 // for time measurement
  
  // Residue Color Transform
  int residue_transform_flag;
#ifdef ADAPTIVE_QUANTIZATION
  int slice_fractional_quant_flag;
  int slice_mqm_signaling_flag;
  int slice_modeling_qm_param0;
  int slice_modeling_qm_param1;
#endif
#ifdef WIN32
  struct _timeb tstruct_start;
  struct _timeb tstruct_end;
#else
  struct timeb tstruct_start;
  struct timeb tstruct_end;
#endif
  
  // picture error concealment
#ifdef ADAPTIVE_FILTER
  int UseAdaptiveFilter;
  int AdaptiveFilterFlag;
#endif
#ifdef DIRECTIONAL_FILTER
  int ImpType;
#endif
#ifdef USE_INTRA_MDDT
  int UseIntraMDDT; 
#endif
#ifdef USE_HP_FILTER    
  int use_high_precision_flag;
#endif
  int last_ref_pic_poc;
  int ref_poc_gap;
  int poc_gap;
  int conceal_mode;
  int earlier_missing_poc;
  unsigned int frame_to_conceal;
  int IDR_concealment_flag;
  int conceal_slice_type;
  
  // random access point decoding
  int recovery_point;           
  int recovery_point_found;
  int recovery_frame_cnt;
  int recovery_frame_num;
  int recovery_poc;
  
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
  int BitDepthIncrease;
  int BitDepthIncreaseChroma;
#endif
  
#ifdef SWITCHED_FILTERS
  int filterFrame[16];                          // Filter selection computed on the sequence
  int filterSequence[16];                       // Filter selection computed on current frame
  int filter[16];                               // Filter selection
  int subpelOffset[2][16];                      // Sub-pel offsets
  int imgOffset[2][16];                         // Image offsets
  int filterParam;
#endif

#ifdef USE_INTRA_MDDT
  byte scanOrder4x4[9][16][2];
  byte scanOrder8x8[9][64][2];
  int  scanStats4x4[9][16];
  int  scanStats8x8[9][64];
  int  update4x4[9];
  int  update8x8[9];
  int  update4x4Count[9];
  int  update8x8Count[9];
  int  update4x4Thres[9];
  int  update8x8Thres[9];
  
  int  cof16x16[16][16];
  byte scanOrder16x16[4][256][2];
#endif 
#ifdef MB32X32
  int mb_ext_level;
  int    dct_coeff16I[16][16];
  int    dct_coeff8I [16][16];
  int    cof16x16P[16][16];
  int    cof16x8[2][16][16];
  int    cof8x16[2][16][16];
#endif

#ifdef ADAPTIVE_LOOP_FILTER
  int  current_alfb_nr;
  int  alf_block_size;
  int  WidthInAlfBlks;
  int  HeightInAlfBlks;
  int  num_alf_block;
  AlfBlock *alf_blk_data;
#endif
} ImageParameters;

extern ImageParameters *img;
extern struct snr_par  *snr;

// signal to noise ratio parameters
struct snr_par
{
  int   frame_ctr;
  float snr_y;                                 //!< current Y SNR
  float snr_u;                                 //!< current U SNR
  float snr_v;                                 //!< current V SNR
  float snr_y1;                                //!< SNR Y(dB) first frame
  float snr_u1;                                //!< SNR U(dB) first frame
  float snr_v1;                                //!< SNR V(dB) first frame
  float snr_ya;                                //!< Average SNR Y(dB) remaining frames
  float snr_ua;                                //!< Average SNR U(dB) remaining frames
  float snr_va;                                //!< Average SNR V(dB) remaining frames
  float sse_y;                                 //!< SSE Y
  float sse_u;                                 //!< SSE U
  float sse_v;                                 //!< SSE V
  float msse_y;                                //!< Average SSE Y
  float msse_u;                                //!< Average SSE U
  float msse_v;                                //!< Average SSE V
};

int tot_time;

#define FILE_NAME_SIZE 200
// input parameters from configuration file
struct inp_par
{
  char infile[FILE_NAME_SIZE];                       //!< H.264 inputfile
  char outfile[FILE_NAME_SIZE];                      //!< Decoded YUV 4:2:0 output
  char reffile[FILE_NAME_SIZE];                      //!< Optional YUV 4:2:0 reference file for SNR measurement
  int FileFormat;                         //!< File format of the Input file, PAR_OF_ANNEXB or PAR_OF_RTP
  int ref_offset;
  int poc_scale;
  int write_uv;
  
#ifdef _LEAKYBUCKET_
  unsigned long R_decoder;                //!< Decoder Rate in HRD Model
  unsigned long B_decoder;                //!< Decoder Buffer size in HRD model
  unsigned long F_decoder;                //!< Decoder Initial buffer fullness in HRD model
  char LeakyBucketParamFile[FILE_NAME_SIZE];         //!< LeakyBucketParamFile
#endif
  
  // picture error concealment
  int conceal_mode;
  int ref_poc_gap;
  int poc_gap;
  
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
  int output_bitdepth;
#endif
};

extern struct inp_par *input;

typedef struct pix_pos
{
  int available;
  int mb_addr;
  int x;
  int y;
  int pos_x;
  int pos_y;
} PixelPos;

typedef struct old_slice_par
{
  unsigned field_pic_flag;
  unsigned bottom_field_flag;
  unsigned frame_num;
  int nal_ref_idc;
  unsigned pic_oder_cnt_lsb;
  int delta_pic_oder_cnt_bottom;
  int delta_pic_order_cnt[2];
  int idr_flag;
  int idr_pic_id;
  int pps_id;
} OldSliceParams;

extern OldSliceParams old_slice;

// files
int p_out;                    //!< file descriptor to output YUV file
//FILE *p_out2;                    //!< pointer to debug output YUV file
int p_ref;                    //!< pointer to input original reference YUV file file

FILE *p_log;                    //!< SNR file

#if TRACE
FILE *p_trace;
#endif

// Residue Color Transform
int mprRGB[3][16][16];
int rec_res[3][16][16];

// Redundant slices
int previous_frame_num;          //!< frame number of previous slice
int ref_flag[17];                //!< 0: i-th previous frame is incorrect 
//!< non-zero: i-th previous frame is correct
int Is_primary_correct;          //!< if primary frame is correct, 0: incorrect
int Is_redundant_correct;        //!< if redundant frame is correct, 0:incorrect
int redundant_slice_ref_idx;     //!< reference index of redundant slice
void Error_tracking(void);       

// prototypes
void init_conf(struct inp_par *inp, char *config_filename);
void report(struct inp_par *inp, struct img_par *img, struct snr_par *snr);
void init(struct img_par *img);

void malloc_slice(struct inp_par *inp, struct img_par *img);
void free_slice(struct inp_par *inp, struct img_par *img);

int  decode_one_frame(struct img_par *img,struct inp_par *inp, struct snr_par *snr);
void init_picture(struct img_par *img, struct inp_par *inp);
void exit_picture();

int  read_new_slice();
void decode_one_slice(struct img_par *img,struct inp_par *inp);

void start_macroblock(struct img_par *img,struct inp_par *inp, int CurrentMBInScanOrder);
int  read_one_macroblock(struct img_par *img,struct inp_par *inp);
void read_ipred_modes(struct img_par *img,struct inp_par *inp);
#ifdef MB32X32
int  decode_one_macroblock(struct img_par *img,struct inp_par *inp, int mb_ext_level, int mb_count);
#else
int  decode_one_macroblock(struct img_par *img,struct inp_par *inp);
#endif
#ifdef ADAPTIVE_FILTER
int  decode_one_macroblock_sepAIF(struct img_par *img,struct inp_par *inp); 
#endif
int  exit_macroblock(struct img_par *img,struct inp_par *inp, int eos_bit);
void decode_ipcm_mb(struct img_par *img);


void readMotionInfoFromNAL (struct img_par *img,struct inp_par *inp);
#ifdef MB32X32
void calcCoeff16x16I(int coeff[16][16]);
void calcCoeff8x8I  (int coeff[16][16]);
void itrans16x16(struct img_par *img);
void itrans16x8 (struct img_par *img, int blk);
void itrans8x16 (struct img_par *img, int blk);
void readLumaCoeffBigBlock_CABAC (struct img_par *img,struct inp_par *inp, int mb32);
void readCBPandCoeffsFromNAL(struct img_par *img,struct inp_par *inp, int mb32, int mbcount);
#else
void readCBPandCoeffsFromNAL(struct img_par *img,struct inp_par *inp);
#endif
void readIPCMcoeffsFromNAL(struct img_par *img, struct inp_par *inp, struct datapartition *dP);

void readLumaCoeff8x8_CABAC (struct img_par *img,struct inp_par *inp, int b8);
void itrans8x8(struct img_par *img, int ioff, int joff);

void copyblock_sp(struct img_par *img,int block_x,int block_y);
void itrans_sp_chroma(struct img_par *img,int ll);
void itrans(struct img_par *img,int ioff,int joff,int i0,int j0, int chroma);
void itrans_sp(struct img_par *img,int ioff,int joff,int i0,int j0);
int  intrapred(struct img_par *img,int ioff,int joff,int i4,int j4);
void itrans_2(struct img_par *img);
int  intrapred_luma_16x16(struct img_par *img,int predmode);
void intrapred_chroma(struct img_par *img, int uv);
int  sign(int a , int b);

// SLICE function pointers
int  (*nal_startcode_follows) ();

// NAL functions TML/CABAC bitstream
int  uvlc_startcode_follows();
int  cabac_startcode_follows();
void free_Partition(Bitstream *currStream);

// ErrorConcealment
void reset_ec_flags();

void error(char *text, int code);
int  is_new_picture();
void init_old_slice();

// dynamic mem allocation
int  init_global_buffers();
void free_global_buffers();

void frame_postprocessing(struct img_par *img, struct inp_par *inp);
void field_postprocessing(struct img_par *img, struct inp_par *inp);
int  bottom_field_picture(struct img_par *img,struct inp_par *inp);
void decode_slice(struct img_par *img,struct inp_par *inp, int current_header);

int RBSPtoSODB(byte *streamBuffer, int last_byte_pos);
int EBSPtoRBSP(byte *streamBuffer, int end_bytepos, int begin_bytepos);

// For MB level frame/field coding
void init_super_macroblock(struct img_par *img,struct inp_par *inp);
void exit_super_macroblock(struct img_par *img,struct inp_par *inp);
int  decode_super_macroblock(struct img_par *img,struct inp_par *inp);
void decode_one_Copy_topMB(struct img_par *img,struct inp_par *inp);

void SetOneRefMV(struct img_par* img);
int peekSyntaxElement_UVLC(SyntaxElement *sym, struct img_par *img, struct inp_par *inp, struct datapartition *dP);

void fill_wp_params(struct img_par *img);

void reset_wp_params(struct img_par *img);

void FreePartition (DataPartition *dp, int n);
DataPartition *AllocPartition();

void tracebits2(const char *trace_str, int len, int info);

void init_decoding_engine_IPCM(struct img_par *img);
void readIPCMBytes_CABAC(SyntaxElement *sym, Bitstream *currStream);

unsigned CeilLog2( unsigned uiVal);

// For Q-matrix
void AssignQuantParam(pic_parameter_set_rbsp_t* pps, seq_parameter_set_rbsp_t* sps);
void CalculateQuantParam(void);
void CalculateQuant8Param(void);

void SetMotionVectorPredictor (struct img_par  *img,
                               short           *pmv_x,
                               short           *pmv_y,
                               char            ref_frame,
                               byte            list,
                               char            ***refPic,
                               short           ****tmp_mv,
                               int             block_x,
                               int             block_y,
                               int             blockshape_x,
                               int             blockshape_y);

#ifdef DIRECTIONAL_FILTER
int FILTCOEF_BITS;
#endif

#ifdef USE_INTRA_MDDT
void precompute_all_inner_product16x16();
void precompute_all_inner_product8x8();
void itransklt16x16_sep_fast(struct img_par *img,int ipredmode);
void itrans8x8klt_sep_fast(struct img_par *img, int ioff, int joff, int ipmode);
void itransklt_sep(struct img_par *img,int ioff,int joff,int i0,int j0, int ipredmode);

void updateScanOrder(int first);
void normalizeScanStats();
#endif 
#endif


