
//#include "424_functions.h"


#define FILE_NAME_SIZE 200
typedef unsigned char byte;    //!< byte type definition
typedef struct
{
  int ProfileIDC;               //!< profile idc
  int LevelIDC;                 //!< level idc
  
  int no_frames;                //!< number of frames to be encoded
  int qp0;                      //!< QP of first frame
  int qpN;                      //!< QP of remaining frames
  int jumpd;                    //!< number of frames to skip in input sequence (e.g 2 takes frame 0,3,6,9...)
#ifdef EIGHTH_PEL
  int mv_res;                   //!< motion vector resolution: 0 = 1/4-pel, 1 = 1/8-pel
#endif
#ifdef MV_COMPETITION
  int mv_competition;
  char predictors_skip[1024];
  char predictors_mvp[1024];
  char predictors_mvb[1024];
#endif
  
  int hadamard;                 /*!< 0: 'normal' SAD in sub pixel search.  1: use 4x4 Hadamard transform and '
  Sum of absolute transform difference' in sub pixel search                   */
  int DisableSubpelME;          //!< Disable Subpixel Motion Estimation
  int search_range;             /*!< search range - integer pel search and 16x16 blocks.  The search window is
                                generally around the predicted vector. Max vector is 2xmcrange.  For 8x8
  and 4x4 block sizes the search range is 1/2 of that for 16x16 blocks.       */
  int num_ref_frames;           //!< number of reference frames to be used
  int P_List0_refs;
  int B_List1_refs;
  int B_List0_refs;
  int Log2MaxFNumMinus4;
  int Log2MaxPOCLsbMinus4;
  int ResendPPS;
  int GenerateMultiplePPS;
  
  int img_width;                //!< image width  (must be a multiple of 16 pels)
  int img_height;               //!< image height (must be a multiple of 16 pels)
  int yuv_format;               //!< YUV format (0=4:0:0, 1=4:2:0, 2=4:2:2, 3=4:4:4)
  int intra_upd;                /*!< For error robustness. 0: no special action. 1: One GOB/frame is intra coded
                                as regular 'update'. 2: One GOB every 2 frames is intra coded etc.
                                In connection with this intra update, restrictions is put on motion vectors
  to prevent errors to propagate from the past                                */
  int blc_size[8][2];           //!< array for different block sizes
  int part_size[8][2];          //!< array for different partition sizes
  int slice_mode;               //!< Indicate what algorithm to use for setting slices
  int slice_argument;           //!< Argument to the specified slice algorithm
  int UseConstrainedIntraPred;  //!< 0: Inter MB pixels are allowed for intra prediction 1: Not allowed
  int  infile_header;           //!< If input file has a header set this to the length of the header
  char infile[FILE_NAME_SIZE];             //!< YUV 4:2:0 input format
  char outfile[FILE_NAME_SIZE];            //!< H.264 compressed output bitstream
  char ReconFile[FILE_NAME_SIZE];          //!< Reconstructed Pictures
  char TraceFile[FILE_NAME_SIZE];          //!< Trace Outputs
  char QmatrixFile[FILE_NAME_SIZE];        //!< Q matrix cfg file
  int intra_period;             //!< Random Access period though intra
  int EnableOpenGOP;            //!< support for open gops.
  
  int idr_enable;         //!< Encode intra slices as IDR
  int start_frame;        //!< Encode sequence starting from Frame start_frame
  
  // B pictures
  int successive_Bframe;        //!< number of B frames that will be used
  int qpB;                      //!< QP for non-reference B slice coded pictures
  int qpBRSOffset;                     //!< QP for reference B slice coded pictures
  int direct_spatial_mv_pred_flag;              //!< Direct Mode type to be used (0: Temporal, 1: Spatial)
  int directInferenceFlag;      //!< Direct Inference Flag
  
  int BiPredMotionEstimation;
  int BiPredMERefinements;
  int BiPredMESearchRange;
  int BiPredMESubPel;
  
  
  // SP Pictures
  int sp_periodicity;           //!< The periodicity of SP-pictures
  int qpsp;                     //!< SP Picture QP for prediction error
  int qpsp_pred;                //!< SP Picture QP for predicted block
  
  int si_frame_indicator;       //!< Flag indicating whether SI frames should be encoded rather than SP frames (0: not used, 1: used)
  int sp2_frame_indicator;      //!< Flag indicating whether switching SP frames should be encoded rather than SP frames (0: not used, 1: used)
  int sp_output_indicator;      //!< Flag indicating whether coefficients are output to allow future encoding of switchin SP frames (0: not used, 1: used)
  char sp_output_filename[FILE_NAME_SIZE]; //!<Filename where SP coefficients are output
  char sp2_input_filename1[FILE_NAME_SIZE]; //!<Filename of coefficients of the first bitstream when encoding SP frames to switch bitstreams 
  char sp2_input_filename2[FILE_NAME_SIZE]; //!<Filenames of coefficients of the second bitstream when encoding SP frames to switch bitstreams
  
  int WeightedPrediction;        //!< Weighted prediciton for P frames (0: not used, 1: explicit)
  int WeightedBiprediction;      //!< Weighted prediciton for B frames (0: not used, 1: explicit, 2: implicit)
  int UseWeightedReferenceME;    //!< Use Weighted Reference for ME.
  int RDPictureDecision;         //!< Perform RD optimal decision between various coded versions of same picture
  int RDPictureIntra;            //!< Enabled RD pic decision for intra as well.
  int RDPSliceWeightOnly;        //!< If enabled, does not check QP variations for P slices.
  int RDPSliceBTest;             //!< Tests B slice replacement for P.      
  int RDBSliceWeightOnly;        //!< If enabled, does not check QP variations for B slices.
  int SkipIntraInInterSlices;    //!< Skip intra type checking in inter slices if best_mode is skip/direct
  int BRefPictures;              //!< B coded reference pictures replace P pictures (0: not used, 1: used)
  int HierarchicalCoding;
  int HierarchyLevelQPEnable;
  char ExplicitHierarchyFormat[1024];  //!< Explicit GOP format (HierarchicalCoding==3). 
  int ReferenceReorder;          //!< Reordering based on Poc distances
  int PocMemoryManagement;       //!< Memory management based on Poc distances for hierarchical coding
  
  int symbol_mode;              //!< Specifies the mode the symbols are mapped on bits
  int of_mode;                  //!< Specifies the mode of the output file
  int partition_mode;           //!< Specifies the mode of data partitioning
  
#ifdef MB32X32
  int InterSearchSkip;
#endif
  int InterSearch16x16;
  int InterSearch16x8;
  int InterSearch8x16;
  int InterSearch8x8;
  int InterSearch8x4;
  int InterSearch4x8;
  int InterSearch4x4;
  
  int IntraDisableInterOnly;
  int Intra4x4ParDisable;
  int Intra4x4DiagDisable;
  int Intra4x4DirDisable;
  int Intra16x16ParDisable;
  int Intra16x16PlaneDisable;
  int ChromaIntraDisable;
  
  int EnableIPCM;
  
  double FrameRate;
  
  int EPZSPattern;
  int EPZSDual;
  int EPZSFixed;
  int EPZSTemporal;
  int EPZSSpatialMem;
  int EPZSMinThresScale;
  int EPZSMaxThresScale;
  int EPZSMedThresScale;
  int EPZSSubPelME;
  int EPZSSubPelThresScale;
  
  int chroma_qp_index_offset;
#ifdef _FULL_SEARCH_RANGE_
  int full_search;
#endif
#ifdef _ADAPT_LAST_GROUP_
  int last_frame;
#endif
#ifdef _CHANGE_QP_
  int qpN2, qpB2, qp2start;
  int qp02, qpBRS2Offset;
#endif
  int rdopt;
  int disthres;
  int nobskip;
  
#ifdef _LEAKYBUCKET_
  int NumberLeakyBuckets;
  char LeakyBucketRateFile[FILE_NAME_SIZE];
  char LeakyBucketParamFile[FILE_NAME_SIZE];
#endif
  
  int PicInterlace;           //!< picture adaptive frame/field
  int MbInterlace;            //!< macroblock adaptive frame/field
  
  int IntraBottom;            //!< Force Intra Bottom at GOP periods.
  
  int LossRateA;              //!< assumed loss probablility of partition A (or full slice), in per cent, used for loss-aware R/D optimization
  int LossRateB;              //!< assumed loss probablility of partition B, in per cent, used for loss-aware R/D 
  int LossRateC;              //!< assumed loss probablility of partition C, in per cent, used for loss-aware R/D 
  int NoOfDecoders;
  int RestrictRef;
  int NumFramesInELSubSeq;
  int NumFrameIn2ndIGOP;
  
  int RandomIntraMBRefresh;     //!< Number of pseudo-random intra-MBs per picture
  
  int LFSendParameters;
  int LFDisableIdc;
  int LFAlphaC0Offset;
  int LFBetaOffset;
  
  int SparePictureOption;
  int SPDetectionThreshold;
  int SPPercentageThreshold;
  
  // FMO
  char SliceGroupConfigFileName[FILE_NAME_SIZE];    //!< Filename for config info fot type 0, 2, 6
  int num_slice_groups_minus1;           //!< "FmoNumSliceGroups" in encoder.cfg, same as FmoNumSliceGroups, which should be erased later
  int slice_group_map_type; 
  
  int *top_left;                         //!< top_left and bottom_right store values indicating foregrounds
  int *bottom_right; 
  byte *slice_group_id;                   //!< slice_group_id is for slice group type being 6  
  int *run_length_minus1;                //!< run_length_minus1 is for slice group type being 0
  
  int slice_group_change_direction_flag;
  int slice_group_change_rate_minus1;
  int slice_group_change_cycle;
  
  int redundant_pic_flag;   //! encoding of redundant pictures
  int pic_order_cnt_type;   //! POC type
  
  int context_init_method;
  int model_number;
  int Transform8x8Mode;
  int ReportFrameStats;
  int DisplayEncParams;
  int Verbose;
  
  //! Rate Control on JVT standard 
  int RCEnable;    
  int bit_rate;
  int SeinitialQP;
  int basicunit;
  int channel_type;
  
  int ScalingMatrixPresentFlag;
  int ScalingListPresentFlag[8];
  
  // FastME enable
  int FMEnable;
  
  // Adaptive filter enable
#ifdef ADAPTIVE_FILTER
  int UseAdaptiveFilter;
  int InterpolationDecision;
#ifdef DIRECTIONAL_FILTER
  int ImpType;//Type of implementation: 0 - floating point 2D, 1- floating point 1D-Diag, 2 - integer 32 bits, 3 - Russanovskyy 16bits, 2 - Kemal's 16bits
#endif
  
#endif
#ifdef ADAPTIVE_FD_SD_CODING
  int  APEC_in_FD_and_SD;
  int  SD_Quantizer;
#endif
#ifdef ADAPTIVE_QUANTIZATION
  int UseAdaptiveQuantMatrix;
#endif
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
  int InputBitDepth;
#endif
#ifdef USE_INTRA_MDDT
  int UseIntraMDDT; 
#endif
#ifdef RDO_Q
  int UseRDO_Q;
#endif
#ifdef USE_HP_FILTER
  int UseHPFilter;
#endif
#ifdef USE_NEW_OFFSET
  int UseNewOffset;
#endif
#ifdef ADAPTIVE_LOOP_FILTER
  int UseAdaptiveLoopFilter;
#ifdef ALF_SPATIAL_PREDICT_COEF
  int ALFPredCoefMode;
#endif
#endif
#ifdef SIMPLIFIED_RDPIC_DECISION
  int SimplifiedRDPicDecision;
#endif
#ifdef MB32X32
  int UseExtMB; //0:16x16MB, 2:64x64MB
#endif
  int DynamicSearchRange;//!< Dynamic Search Range
  int FMEScale;
  //////////////////////////////////////////////////////////////////////////
  // Fidelity Range Extensions
  int BitDepthLuma;
  int BitDepthChroma;
  int img_height_cr;
  int img_width_cr;
  int rgb_input_flag;
  int cb_qp_index_offset;
  int cr_qp_index_offset;
  
  // Lossless Coding
  int lossless_qpprime_y_zero_flag;
  
  //Residue Color Transform
  int residue_transform_flag;
  
  // Lambda Params
  int UseExplicitLambdaParams;
  double LambdaWeight[6];
  
  char QOffsetMatrixFile[FILE_NAME_SIZE];        //!< Quantization Offset matrix cfg file
  int  OffsetMatrixPresentFlag;                  //!< Enable Explicit Quantization Offset Matrices
  
  int AdaptiveRounding;                          //!< Adaptive Rounding parameter based on JVT-N011
  int AdaptRndPeriod;                            //!< Set period for adaptive rounding of JVT-N011 in MBs
  int AdaptRndChroma;
  int AdaptRndWFactor[2][5];                     //!< Weighting factors based on reference indicator and slice type 
  // Fast Mode Decision
  int EarlySkipEnable;
  int SelectiveIntraEnable;
  int DisposableP;
  int DispPQPOffset;
  
  //Redundant picture
  int NumRedundantHierarchy;   //!< number of entries to allocate redundant pictures in a GOP
  int PrimaryGOPLength;        //!< GOP length of primary pictures  
  int NumRefPrimary;           //!< number of reference frames for primary picture
} InputParameters;