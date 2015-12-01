
/*!
 ***************************************************************************
 * \file
 *    cabac.h
 *
 * \brief
 *    Headerfile for entropy coding routines
 *
 * \author
 *    Detlev Marpe                                                         \n
 *    Copyright (C) 2000 HEINRICH HERTZ INSTITUTE All Rights Reserved.
 *
 * \date
 *    21. Oct 2000 (Changes by Tobias Oelbaum 28.08.2001)
 ***************************************************************************
 */

#ifndef _CABAC_H_
#define _CABAC_H_

#include "global.h"

MotionInfoContexts* create_contexts_MotionInfo(void);
TextureInfoContexts* create_contexts_TextureInfo(void);
void init_contexts_MotionInfo(struct img_par *img, MotionInfoContexts *enco_ctx);
void init_contexts_TextureInfo(struct img_par *img, TextureInfoContexts *enco_ctx);
void delete_contexts_MotionInfo(MotionInfoContexts *enco_ctx);
void delete_contexts_TextureInfo(TextureInfoContexts *enco_ctx);

void cabac_new_slice();

void readMB_typeInfo_CABAC(SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readB8_typeInfo_CABAC(SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readIntraPredMode_CABAC(SyntaxElement *se, struct inp_par *inp,struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readRefFrame_CABAC(SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readMVD_CABAC(SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readCBP_CABAC(SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readRunLevel_CABAC(SyntaxElement *se, struct inp_par *inp, struct img_par *img,  DecodingEnvironmentPtr dep_dp);
#ifdef USE_INTRA_MDDT
void readRunLevel16x16_CABAC(SyntaxElement *se, struct inp_par *inp, struct img_par *img,  DecodingEnvironmentPtr dep_dp);
#endif 
void readDquant_CABAC(SyntaxElement *se,struct inp_par *inp,struct img_par *img,DecodingEnvironmentPtr dep_dp);
void readCIPredMode_CABAC(SyntaxElement *se,struct inp_par *inp,struct img_par *img,DecodingEnvironmentPtr dep_dp);
void readMB_skip_flagInfo_CABAC( SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readFieldModeInfo_CABAC(SyntaxElement *se,struct inp_par *inp,struct img_par *img,DecodingEnvironmentPtr dep_dp); 

void readMB_transform_size_flag_CABAC( SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);

int  readSyntaxElement_CABAC(SyntaxElement *se, struct img_par *img, struct inp_par *inp, DataPartition *this_dataPart);

int  check_next_mb_and_get_field_mode_CABAC(SyntaxElement *se,struct img_par *img,struct inp_par *inp,DataPartition  *act_dp);
void CheckAvailabilityOfNeighborsCABAC();

#ifdef ADAPTIVE_FD_SD_CODING
void read_adaptive_prederror_coding_flag_for_MB( SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void set_bit(int *dest, int bit, int value);
int  read_spatial_domain_coding_flag_for_4x4block (Macroblock *currMB, DecodingEnvironmentPtr  dep_dp, struct img_par *img, int type);
int  read_significance_map_SD (Macroblock *currMB, DecodingEnvironmentPtr  dep_dp, struct img_par *img, int type, int coeff[]);
void read_significant_coefficients_SD (Macroblock *currMB, DecodingEnvironmentPtr dep_dp, struct img_par *img, int type, int coeff[]);
int  read_spatial_domain_coding_flag_for_8x8block (Macroblock *currMB, DecodingEnvironmentPtr dep_dp, struct img_par *img, int type, int writing_b8);
int  read_significance_map_SD8(Macroblock *currMB, DecodingEnvironmentPtr dep_dp, struct img_par *img, int type, int coeff[]);
void read_significant_coefficients_SD8(Macroblock *currMB, DecodingEnvironmentPtr dep_dp,struct img_par *img, int type, int coeff[]);
#endif

#ifdef ADAPTIVE_QUANTIZATION
void readIAQMS_idx_CABAC( SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
#endif

#ifdef ADAPTIVE_LOOP_FILTER
void readAlfBlockFlagCABAC(SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void CheckAvailabilityOfNeighborsAlfBlockCABAC(void);
#endif

#endif  // _CABAC_H_

