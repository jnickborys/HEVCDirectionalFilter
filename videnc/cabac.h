
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

#ifdef RDO_Q
#define MAX_PREC_COEFF    25
#endif
// CABAC
int get_pic_bin_count(void);
void reset_pic_bin_count(void);

void arienco_start_encoding(EncodingEnvironmentPtr eep, unsigned char *code_buffer, int *code_len);
int  arienco_bits_written(EncodingEnvironmentPtr eep);
void arienco_done_encoding(EncodingEnvironmentPtr eep);
void biari_init_context (BiContextTypePtr ctx, const int* ini);
void rescale_cum_freq(BiContextTypePtr bi_ct);
#ifdef RDO_Q
int biari_no_bits(EncodingEnvironmentPtr eep, signed short symbol, BiContextTypePtr bi_ct );
int biari_state(EncodingEnvironmentPtr eep, signed short symbol, BiContextTypePtr bi_ct );
int est_unary_exp_golomb_level_encode(unsigned int symbol, int ctx, int type);
void precalculate_unary_exp_golomb_level();
#endif
void biari_encode_symbol(EncodingEnvironmentPtr eep, signed short symbol, BiContextTypePtr bi_ct );
void biari_encode_symbol_eq_prob(EncodingEnvironmentPtr eep, signed short symbol);
void biari_encode_symbol_final(EncodingEnvironmentPtr eep, signed short symbol);
MotionInfoContexts* create_contexts_MotionInfo(void);
TextureInfoContexts* create_contexts_TextureInfo(void);
void init_contexts_MotionInfo (MotionInfoContexts  *enco_ctx);
void init_contexts_TextureInfo(TextureInfoContexts *enco_ctx);
void delete_contexts_MotionInfo(MotionInfoContexts *enco_ctx);
void delete_contexts_TextureInfo(TextureInfoContexts *enco_ctx);
void writeHeaderToBuffer(void);
int  writeSyntaxElement_CABAC(SyntaxElement *se, DataPartition *this_dataPart);
void writeMB_typeInfo_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeIntraPredMode_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeB8_typeInfo_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeRefFrame2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeRefFrame_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeMVD_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeCBP_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeDquant_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeRunLevel_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
#ifdef USE_INTRA_MDDT
void writeRunLevel16x16_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
#endif 
#ifdef MB32X32
void writeRunLevel16x16P_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeRunLevel16x8_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeRunLevel8x16_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
#endif
void writeBiDirBlkSize_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeCIPredMode_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void print_ctx_TextureInfo(TextureInfoContexts *enco_ctx);
void writeMB_skip_flagInfo_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeFieldModeInfo_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp); //GB
void writeCBP_BIT_CABAC (int b8, int bit, int cbp, Macroblock* currMB, int inter, EncodingEnvironmentPtr eep_dp);
void cabac_new_slice(void);
void CheckAvailabilityOfNeighborsCABAC(void);

void writeMB_transform_size_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);

#ifdef ADAPTIVE_QUANTIZATION
void writeIAQMS_idx_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
#endif
#ifdef ADAPTIVE_FD_SD_CODING
void write_adaptive_prederror_coding_flag_for_MB(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void set_bit(int *dest, int bit, int value);
#endif
#ifdef ADAPTIVE_LOOP_FILTER
void add_pic_bin_count(int val);
void arienco_start_encoding_AlfBlockFlag(EncodingEnvironmentPtr eep, unsigned char *code_buffer, int *code_len);
void writeAlfBlockFlagCABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void CheckAvailabilityOfNeighborsAlfCABAC();
#endif
#endif  // CABAC_H

