
/*!
 *************************************************************************************
 * \file cabac.c
 *
 * \brief
 *    CABAC entropy coding routines
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Detlev Marpe                    <marpe@hhi.de>
 **************************************************************************************
 */

#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include "global.h"

#include "cabac.h"
#include "image.h"
#include "mb_access.h"
#ifdef ADAPTIVE_LOOP_FILTER
#include "adaptive_loop_filter.h"
#endif
#ifdef ADAPTIVE_FD_SD_CODING
#include "spatial_domain_coding.h"

extern int cabac_encoding;
#endif

#ifdef RDO_Q
estBitsCabacStruct estBitsCabac[NUM_BLOCK_TYPES];
int precalcUnaryLevelTab[128][MAX_PREC_COEFF];
extern int entropyBits[128];
#endif

static const int type2ctx_bcbp[] = { 0,  1,  2,  2,  3,  4,  5,  6,  5,  5}; // 7
#ifdef USE_INTRA_MDDT
#ifdef MB32X32
static const int type2ctx_map [] = { 0,  1,  2,  3,  4,  5,  6,  7,  6,  6, 8, 9, 10, 11}; 
static const int type2ctx_last[] = { 0,  1,  2,  3,  4,  5,  6,  7,  6,  6, 8, 9, 10, 11}; 
static const int type2ctx_one [] = { 0,  1,  2,  3,  3,  4,  5,  6,  5,  5, 7, 8,  9, 10}; 
static const int type2ctx_abs [] = { 0,  1,  2,  3,  3,  4,  5,  6,  5,  5, 7, 8,  9, 10}; 
static const int max_c2       [] = { 4,  4,  4,  4,  4,  4,  3,  4,  3,  3, 4, 4,  4,  4}; 
#else
static const int type2ctx_map [] = { 0,  1,  2,  3,  4,  5,  6,  7,  6,  6, 8}; // 9
static const int type2ctx_last[] = { 0,  1,  2,  3,  4,  5,  6,  7,  6,  6, 8}; // 9
static const int type2ctx_one [] = { 0,  1,  2,  3,  3,  4,  5,  6,  5,  5, 7}; // 8
static const int type2ctx_abs [] = { 0,  1,  2,  3,  3,  4,  5,  6,  5,  5, 7}; // 8
static const int max_c2       [] = { 4,  4,  4,  4,  4,  4,  3,  4,  3,  3, 4}; // 9
#endif
#else
static const int type2ctx_map [] = { 0,  1,  2,  3,  4,  5,  6,  7,  6,  6}; // 8
static const int type2ctx_last[] = { 0,  1,  2,  3,  4,  5,  6,  7,  6,  6}; // 8
static const int type2ctx_one [] = { 0,  1,  2,  3,  3,  4,  5,  6,  5,  5}; // 7
static const int type2ctx_abs [] = { 0,  1,  2,  3,  3,  4,  5,  6,  5,  5}; // 7
static const int max_c2       [] = { 4,  4,  4,  4,  4,  4,  3,  4,  3,  3}; // 9
#endif

#define BIT_SET(x,n)  ((int)(((x)&((int64)1<<(n)))>>(n)))

int last_dquant = 0;

/***********************************************************************
 * L O C A L L Y   D E F I N E D   F U N C T I O N   P R O T O T Y P E S
 ***********************************************************************
 */


void unary_bin_encode(EncodingEnvironmentPtr eep_frame,
                      unsigned int symbol,
                      BiContextTypePtr ctx,
                      int ctx_offset);

void unary_bin_max_encode(EncodingEnvironmentPtr eep_frame,
                          unsigned int symbol,
                          BiContextTypePtr ctx,
                          int ctx_offset,
                          unsigned int max_symbol);

void unary_exp_golomb_level_encode( EncodingEnvironmentPtr eep_dp,
                                   unsigned int symbol,
                                   BiContextTypePtr ctx);

void unary_exp_golomb_mv_encode(EncodingEnvironmentPtr eep_dp,
                                unsigned int symbol,
                                BiContextTypePtr ctx,
                                unsigned int max_bin);


void cabac_new_slice(void)
{
  last_dquant=0;
}


/*!
 ************************************************************************
 * \brief
 *    Check for available neighbouring blocks
 *    and set pointers in current macroblock
 ************************************************************************
 */
void CheckAvailabilityOfNeighborsCABAC(void)
{
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  PixelPos up, left;
  
  getNeighbour(img->current_mb_nr, -1,  0, 1, &left);
  getNeighbour(img->current_mb_nr,  0, -1, 1, &up);
  
  if (up.available)
    currMB->mb_available_up = &img->mb_data[up.mb_addr];
  else
    currMB->mb_available_up = NULL;
  
  if (left.available)
    currMB->mb_available_left = &img->mb_data[left.mb_addr];
  else
    currMB->mb_available_left = NULL;
}

/*!
 ************************************************************************
 * \brief
 *    Allocation of contexts models for the motion info
 *    used for arithmetic encoding
 ************************************************************************
 */
MotionInfoContexts* create_contexts_MotionInfo(void)
{
  MotionInfoContexts* enco_ctx;

  enco_ctx = (MotionInfoContexts*) calloc(1, sizeof(MotionInfoContexts) );
  if( enco_ctx == NULL )
    no_mem_exit("create_contexts_MotionInfo: enco_ctx");

  return enco_ctx;
}


/*!
 ************************************************************************
 * \brief
 *    Allocates of contexts models for the texture info
 *    used for arithmetic encoding
 ************************************************************************
 */
TextureInfoContexts* create_contexts_TextureInfo(void)
{
  TextureInfoContexts*  enco_ctx;

  enco_ctx = (TextureInfoContexts*) calloc(1, sizeof(TextureInfoContexts) );
  if( enco_ctx == NULL )
    no_mem_exit("create_contexts_TextureInfo: enco_ctx");

  return enco_ctx;
}




/*!
 ************************************************************************
 * \brief
 *    Frees the memory of the contexts models
 *    used for arithmetic encoding of the motion info.
 ************************************************************************
 */
void delete_contexts_MotionInfo(MotionInfoContexts *enco_ctx)
{
  if( enco_ctx == NULL )
    return;

  free( enco_ctx );

  return;
}

/*!
 ************************************************************************
 * \brief
 *    Frees the memory of the contexts models
 *    used for arithmetic encoding of the texture info.
 ************************************************************************
 */
void delete_contexts_TextureInfo(TextureInfoContexts *enco_ctx)
{
  if( enco_ctx == NULL )
    return;

  free( enco_ctx );

  return;
}


/*!
 **************************************************************************
 * \brief
 *    generates arithmetic code and passes the code to the buffer
 **************************************************************************
 */
int writeSyntaxElement_CABAC(SyntaxElement *se, DataPartition *this_dataPart)
{
  EncodingEnvironmentPtr eep_dp = &(this_dataPart->ee_cabac);
  int curr_len = arienco_bits_written(eep_dp);

  // perform the actual coding by calling the appropriate method
  se->writing(se, eep_dp);

  if(se->type != SE_HEADER)
    this_dataPart->bitstream->write_flag = 1;

  return (se->len = (arienco_bits_written(eep_dp) - curr_len));
}

/*!
 ***************************************************************************
 * \brief
 *    This function is used to arithmetically encode the field
 *    mode info of a given MB  in the case of mb-based frame/field decision
 ***************************************************************************
 */
void writeFieldModeInfo_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  int a,b,act_ctx;
  MotionInfoContexts *ctx         = (img->currentSlice)->mot_ctx;
  Macroblock         *currMB      = &img->mb_data[img->current_mb_nr];
  int                mb_field     = se->value1;
  
  a = currMB->mbAvailA ? img->mb_data[currMB->mbAddrA].mb_field : 0;
  b = currMB->mbAvailB ? img->mb_data[currMB->mbAddrB].mb_field : 0;
  
  act_ctx = a + b;
  
  biari_encode_symbol(eep_dp, (signed short) (mb_field != 0),&ctx->mb_aff_contexts[act_ctx]);
  
  se->context = act_ctx;
  
  return;
}

#ifdef ADAPTIVE_FD_SD_CODING
void write_adaptive_prederror_coding_flag_for_MB(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  int a, b;
  int ctx = 0;
  int act_sym;

  Macroblock         *currMB      = &img->mb_data[img->current_mb_nr];

  b = (currMB->mb_available_up   == NULL) ? 0 : currMB->mb_available_up  ->SD_Coding_on_off;
  a = (currMB->mb_available_left == NULL) ? 0 : currMB->mb_available_left->SD_Coding_on_off;

  ctx         = a + 2 * b;
  act_sym     = se->value1;

  biari_encode_symbol (eep_dp, (short)act_sym, img->currentSlice->tex_ctx->MB_adaptive_SD_context+ctx);
}

void write_spatial_domain_coding_flag_for_4x4block (Macroblock* currMB, EncodingEnvironmentPtr eep_dp, int type, int FD_SD_bit)
{
  int y_ac        = (type==LUMA_16AC || type==LUMA_8x8 || type==LUMA_8x4 || type==LUMA_4x8 || type==LUMA_4x4);
  int y_dc        = (type==LUMA_16DC);
  int u_ac        = (type==CHROMA_AC && !img->is_v_block);
  int v_ac        = (type==CHROMA_AC &&  img->is_v_block);
  int chroma_dc   = (type==CHROMA_DC || type==CHROMA_DC_2x4 || type==CHROMA_DC_4x4);
  int u_dc        = (chroma_dc && !img->is_v_block);
  int v_dc        = (chroma_dc &&  img->is_v_block);
  int j           = (y_ac || u_ac || v_ac ? img->subblock_y : 0);
  int i           = (y_ac || u_ac || v_ac ? img->subblock_x : 0);
  int bit         = (y_dc ? 0 : y_ac ? 1 : u_dc ? 17 : v_dc ? 18 : u_ac ? 19 : 23);
  int default_bit = (img->is_intra_block ? 1 : 0);
  int upper_bit   = default_bit;
  int left_bit    = default_bit;
  int ctx;

  int bit_pos_a   = 0;
  int bit_pos_b   = 0;

  PixelPos block_a, block_b;

  if (y_ac || y_dc)
  {
    getLuma4x4Neighbour(img->current_mb_nr, i, j, -1,  0, &block_a);
    getLuma4x4Neighbour(img->current_mb_nr, i, j,  0, -1, &block_b);
    if (y_ac)
    {
      if (block_a.available)
        bit_pos_a = 4*block_a.y + block_a.x;
      if (block_b.available)
        bit_pos_b = 4*block_b.y + block_b.x;
    }
  }
  else
  {
    getChroma4x4Neighbour(img->current_mb_nr, i, j, -1,  0, &block_a);
    getChroma4x4Neighbour(img->current_mb_nr, i, j,  0, -1, &block_b);
    if (u_ac||v_ac)
    {
      if (block_a.available)
        bit_pos_a = 4*block_a.y + block_a.x;
      if (block_b.available)
        bit_pos_b = 4*block_b.y + block_b.x;
    }
  }

  bit = (y_dc ? 0 : y_ac ? 1+4*j+i : u_dc ? 17 : v_dc ? 18 : u_ac ? 19+4*j+i : 35+4*j+i);
  //--- set bits for current block ---
  if (FD_SD_bit)
  {
    if (type==LUMA_8x8)
    {
      currMB->FD_or_SD_bits   |= ((int64)1<< bit   );
      currMB->FD_or_SD_bits   |= ((int64)1<<(bit+1));
      currMB->FD_or_SD_bits   |= ((int64)1<<(bit+4));
      currMB->FD_or_SD_bits   |= ((int64)1<<(bit+5));
    }
    else if (type==LUMA_8x4)
    {
      currMB->FD_or_SD_bits   |= ((int64)1<< bit   );
      currMB->FD_or_SD_bits   |= ((int64)1<<(bit+1));
    }
    else if (type==LUMA_4x8)
    {
      currMB->FD_or_SD_bits   |= ((int64)1<< bit   );
      currMB->FD_or_SD_bits   |= ((int64)1<<(bit+4));
    }
    else
    {
      currMB->FD_or_SD_bits   |= ((int64)1<<bit);
    }
  }

  bit = (y_dc ? 0 : y_ac ? 1 : u_dc ? 17 : v_dc ? 18 : u_ac ? 19 : 35);
  if (type!=LUMA_8x8)
  {
    if (block_b.available)
    {
      if(img->mb_data[block_b.mb_addr].mb_type==IPCM)
        upper_bit=1;
      else
        upper_bit = BIT_SET(img->mb_data[block_b.mb_addr].FD_or_SD_bits,bit+bit_pos_b);
    }


    if (block_a.available)
    {
      if(img->mb_data[block_a.mb_addr].mb_type==IPCM)
        left_bit=1;
      else
        left_bit = BIT_SET(img->mb_data[block_a.mb_addr].FD_or_SD_bits,bit+bit_pos_a);
    }

    ctx = 2*upper_bit+left_bit;

    //===== encode symbol =====
    biari_encode_symbol (eep_dp, (short)FD_SD_bit, img->currentSlice->tex_ctx->bcbp_contexts_FD_SD+ctx);
  }
}


void write_spatial_domain_coding_flag_for_8x8block  (Macroblock* currMB, EncodingEnvironmentPtr eep_dp, int type, int FD_SD_bit, int writing_b8)
{
  int ctx=0;

  int upper_bit   = 0;
  int left_bit    = 0;

  if (writing_b8==0)
  {
    upper_bit = (currMB->mb_available_up   == NULL) ? 0 : (((currMB->mb_available_up  ->SD_or_FD_t8x8) & 4)==0?0:1);
    left_bit  = (currMB->mb_available_left == NULL) ? 0 : (((currMB->mb_available_left->SD_or_FD_t8x8) & 2)==0?0:1);
  }
  else if (writing_b8==1)
  {
    upper_bit = (currMB->mb_available_up   == NULL) ? 0 : (((currMB->mb_available_up  ->SD_or_FD_t8x8) & 8)==0?0:1);
    left_bit  =                                                              (((currMB->SD_or_FD_t8x8) & 1)==0?0:1);
  }
  else if (writing_b8==2)
  {
    upper_bit =                                                              (((currMB->SD_or_FD_t8x8) & 1)==0?0:1);
    left_bit  = (currMB->mb_available_left == NULL) ? 0 : (((currMB->mb_available_left->SD_or_FD_t8x8) & 8)==0?0:1);
  }
  else
  {
    upper_bit =                                                              (((currMB->SD_or_FD_t8x8) & 2)==0?0:1);
    left_bit  =                                                              (((currMB->SD_or_FD_t8x8) & 4)==0?0:1);
  }

  ctx=2*upper_bit+left_bit;

  //===== encode symbol =====
  biari_encode_symbol (eep_dp, (short)FD_SD_bit, img->currentSlice->tex_ctx->bcbp8_contexts_FD_SD+ctx);
}
#endif




/*!
***************************************************************************
* \brief
*    This function is used to arithmetically encode the mb_skip_flag.
***************************************************************************
*/
void writeMB_skip_flagInfo_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  int a,b,act_ctx;
  int bframe   = (img->type==B_SLICE);
  MotionInfoContexts *ctx         = (img->currentSlice)->mot_ctx;
  Macroblock         *currMB      = &img->mb_data[img->current_mb_nr];
  int                curr_mb_type = se->value1;
  
  if (bframe)
  {
    if (currMB->mb_available_up == NULL)
      b = 0;
    else
      b = (currMB->mb_available_up->skip_flag==0 ? 1 : 0);
    if (currMB->mb_available_left == NULL)
      a = 0;
    else
      a = (currMB->mb_available_left->skip_flag==0 ? 1 : 0);
    
    act_ctx = 7 + a + b;

    if (se->value1==0 && se->value2==0) // DIRECT mode, no coefficients
      biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][act_ctx]);
    else
      biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[2][act_ctx]);   

    currMB->skip_flag = (se->value1==0 && se->value2==0)?1:0;
  }
  else
  {
    if (currMB->mb_available_up == NULL)
      b = 0;
    else
      b = (( (currMB->mb_available_up)->skip_flag == 0) ? 1 : 0 );
    if (currMB->mb_available_left == NULL)
      a = 0;
    else
      a = (( (currMB->mb_available_left)->skip_flag == 0) ? 1 : 0 );

    act_ctx = a + b;

    if (curr_mb_type==0) // SKIP
      biari_encode_symbol(eep_dp, 1,&ctx->mb_type_contexts[1][act_ctx]);
    else
      biari_encode_symbol(eep_dp, 0,&ctx->mb_type_contexts[1][act_ctx]);

    currMB->skip_flag = (curr_mb_type==0)?1:0;
  }
  se->context = act_ctx;

  return;
}

/*!
***************************************************************************
* \brief
*    This function is used to arithmetically encode the macroblock
*    intra_pred_size flag info of a given MB.
***************************************************************************
*/

void writeMB_transform_size_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  int a, b;
  int act_ctx = 0;
  int act_sym;
  
  MotionInfoContexts *ctx         = (img->currentSlice)->mot_ctx;
  Macroblock         *currMB      = &img->mb_data[img->current_mb_nr];
  
#ifdef MB32X32
  if(input->UseExtMB  > 0)
  {
    if(currMB->mb_type >= 0 && currMB->mb_type <= 3)
    {
      b = (currMB->mb_available_up == NULL) ? 0 : currMB->mb_available_up->luma_transform_size_8x8_flag ? 1:0;
      a = (currMB->mb_available_left == NULL) ? 0 :currMB->mb_available_left->luma_transform_size_8x8_flag ? 1:0;

      act_ctx     = a + b;
      assert(act_ctx <=2);
      act_sym     = currMB->luma_transform_size_8x8_flag ? 1:0;
      se->context = act_ctx; // store context
      biari_encode_symbol(eep_dp, (signed short) (act_sym != 0), ctx->transform_size_contexts + act_ctx );

      // send another bit to indicate either 8x8 or 16x16
      if(currMB->luma_transform_size_8x8_flag)
      {
        b = (currMB->mb_available_up == NULL) ? 0 : currMB->mb_available_up->luma_transform_size_8x8_flag ? 1:0;
        a = (currMB->mb_available_left == NULL) ? 0 :currMB->mb_available_left->luma_transform_size_8x8_flag ? 1:0;

        act_ctx     = a + b;
        assert(act_ctx <=2);
        act_sym     = (currMB->luma_transform_size_8x8_flag == 2) ? 1:0;
        se->context = act_ctx; // store context
        biari_encode_symbol(eep_dp, (signed short) (act_sym != 0), ctx->transform_size_contexts + act_ctx );
      }
    }
    else
    {
      b = (currMB->mb_available_up == NULL) ? 0 : currMB->mb_available_up->luma_transform_size_8x8_flag ? 1:0;
      a = (currMB->mb_available_left == NULL) ? 0 :currMB->mb_available_left->luma_transform_size_8x8_flag ? 1:0;

      act_ctx     = a + b;
      act_sym     = currMB->luma_transform_size_8x8_flag;
      se->context = act_ctx; // store context
      biari_encode_symbol(eep_dp, (signed short) (act_sym != 0), ctx->transform_size_contexts + act_ctx );
    }
  }
  else
  {
#endif
    b = (currMB->mb_available_up == NULL) ? 0 : currMB->mb_available_up->luma_transform_size_8x8_flag;
    a = (currMB->mb_available_left == NULL) ? 0 :currMB->mb_available_left->luma_transform_size_8x8_flag;
      
    act_ctx     = a + b;
    act_sym     = currMB->luma_transform_size_8x8_flag;
    se->context = act_ctx; // store context
    biari_encode_symbol(eep_dp, (signed short) (act_sym != 0), ctx->transform_size_contexts + act_ctx );  
#ifdef MB32X32
  }
#endif
}

/*!
 ***************************************************************************
 * \brief
 *    This function is used to arithmetically encode the macroblock
 *    type info of a given MB.
 ***************************************************************************
 */

void writeMB_typeInfo_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  int a, b;
  int act_ctx = 0;
  int act_sym;
  signed short csym;
  int bframe   = (img->type==B_SLICE);
  int mode_sym = 0;
  int mode16x16;


  MotionInfoContexts *ctx         = (img->currentSlice)->mot_ctx;
  Macroblock         *currMB      = &img->mb_data[img->current_mb_nr];
  int                curr_mb_type = se->value1;

  if(img->type == I_SLICE)  // INTRA-frame
  {
    if (currMB->mb_available_up == NULL)
      b = 0;
    else 
      b = ((currMB->mb_available_up->mb_type != I4MB &&  currMB->mb_available_up->mb_type != I8MB) ? 1 : 0 );

    if (currMB->mb_available_left == NULL)
      a = 0;
    else 
      a = ((currMB->mb_available_left->mb_type != I4MB &&  currMB->mb_available_left->mb_type != I8MB) ? 1 : 0 );
    
    act_ctx     = a + b;
    act_sym     = curr_mb_type;
    se->context = act_ctx; // store context

    if (act_sym==0) // 4x4 Intra
    {
      biari_encode_symbol(eep_dp, 0, ctx->mb_type_contexts[0] + act_ctx );
    }
    else if( act_sym == 25 ) // PCM-MODE
    {
      biari_encode_symbol(eep_dp, 1, ctx->mb_type_contexts[0] + act_ctx );
      biari_encode_symbol_final(eep_dp, 1);
    }
    else // 16x16 Intra
    {
      biari_encode_symbol(eep_dp, 1, ctx->mb_type_contexts[0] + act_ctx );

      biari_encode_symbol_final(eep_dp, 0);

      mode_sym = act_sym-1; // Values in the range of 0...23
      act_ctx  = 4;
      act_sym  = mode_sym/12;
      biari_encode_symbol(eep_dp, (signed short) act_sym, ctx->mb_type_contexts[0] + act_ctx ); // coding of AC/no AC
      mode_sym = mode_sym % 12;
      act_sym  = mode_sym / 4; // coding of cbp: 0,1,2
      act_ctx  = 5;
      if (act_sym==0)
      {
        biari_encode_symbol(eep_dp, 0, ctx->mb_type_contexts[0] + act_ctx );
      }
      else
      {
        biari_encode_symbol(eep_dp, 1, ctx->mb_type_contexts[0] + act_ctx );
        act_ctx=6;
        biari_encode_symbol(eep_dp, (signed short) (act_sym!=1), ctx->mb_type_contexts[0] + act_ctx );
      }
      mode_sym = mode_sym & 0x03; // coding of I pred-mode: 0,1,2,3
      act_sym  = mode_sym >> 1;
      act_ctx  = 7;
      biari_encode_symbol(eep_dp, (signed short) act_sym, ctx->mb_type_contexts[0] + act_ctx );
      act_ctx  = 8;
      act_sym  = mode_sym & 0x01;
      biari_encode_symbol(eep_dp, (signed short) act_sym, ctx->mb_type_contexts[0] + act_ctx );
    }
  }
  else // INTER
  {
    
    if (bframe)
    {
      if (currMB->mb_available_up == NULL)
        b = 0;
      else
        b = ((currMB->mb_available_up->mb_type != 0) ? 1 : 0 );

      if (currMB->mb_available_left == NULL)
        a = 0;
      else
        a = ((currMB->mb_available_left->mb_type != 0) ? 1 : 0 );
      act_ctx = a + b;
      se->context = act_ctx; // store context
    }
    act_sym = curr_mb_type;

    if (act_sym>=(mode16x16=(bframe?24:7)))
    {
      mode_sym = act_sym-mode16x16;
      act_sym  = mode16x16; // 16x16 mode info
    }

    if (!bframe)
    {
      switch (act_sym)
      {
      case 0:
        break;
      case 1:
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][4]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][5]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][6]);
        break;
      case 2:
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][4]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][5]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][7]);
        break;
      case 3:
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][4]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][5]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][7]);
        break;
      case 4:
      case 5:
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][4]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][5]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][6]);
        break;
      case 6:
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][4]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][7]);
        break;
      case 7:
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][4]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][7]);
        break;
      default:
        printf ("Unsupported MB-MODE in writeMB_typeInfo_CABAC!\n");
        exit (1);
      }
    }
    else //===== B-FRAMES =====
    {
      if (act_sym==0)
      {
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[2][act_ctx]);
      }
      else if (act_sym<=2)
      {
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][act_ctx]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[2][4]);
        csym = (act_sym-1 != 0);
        biari_encode_symbol (eep_dp, csym, &ctx->mb_type_contexts[2][6]);
      }
      else if (act_sym<=10)
      {
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][act_ctx]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][4]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[2][5]);
        csym=(((act_sym-3)>>2)&0x01) != 0;
        biari_encode_symbol (eep_dp, csym, &ctx->mb_type_contexts[2][6]);
        csym=(((act_sym-3)>>1)&0x01) != 0;
        biari_encode_symbol (eep_dp, csym, &ctx->mb_type_contexts[2][6]);
        csym=((act_sym-3)&0x01) != 0;
        biari_encode_symbol (eep_dp, csym, &ctx->mb_type_contexts[2][6]);
      }
      else if (act_sym==11 || act_sym==22)
      {
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][act_ctx]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][4]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][5]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][6]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][6]);
        csym = (act_sym != 11);
        biari_encode_symbol (eep_dp, csym, &ctx->mb_type_contexts[2][6]);
      }
      else
      {
        if (act_sym > 22) act_sym--;
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][act_ctx]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][4]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][5]);
        csym=(((act_sym-12)>>3)&0x01) != 0;
        biari_encode_symbol (eep_dp, csym, &ctx->mb_type_contexts[2][6]);
        csym=(((act_sym-12)>>2)&0x01) != 0;
        biari_encode_symbol (eep_dp, csym, &ctx->mb_type_contexts[2][6]);
        csym=(((act_sym-12)>>1)&0x01) != 0;
        biari_encode_symbol (eep_dp, csym, &ctx->mb_type_contexts[2][6]);
        csym=((act_sym-12)&0x01) != 0;
        biari_encode_symbol (eep_dp, csym, &ctx->mb_type_contexts[2][6]);         
        if (act_sym >=22) act_sym++;
      }
    }

    if(act_sym==mode16x16) // additional info for 16x16 Intra-mode
    {
      if( mode_sym==25 )
      {
        biari_encode_symbol_final(eep_dp, 1 );
        return;
      }
      biari_encode_symbol_final(eep_dp, 0 );

      act_ctx = 8;
      act_sym = mode_sym/12;
      biari_encode_symbol(eep_dp, (signed short) act_sym, ctx->mb_type_contexts[1] + act_ctx ); // coding of AC/no AC
      mode_sym = mode_sym % 12;

      act_sym = mode_sym / 4; // coding of cbp: 0,1,2
      act_ctx = 9;
      if (act_sym==0)
      {
        biari_encode_symbol(eep_dp, 0, ctx->mb_type_contexts[1] + act_ctx );
      }
      else
      {
        biari_encode_symbol(eep_dp, 1, ctx->mb_type_contexts[1] + act_ctx );
        biari_encode_symbol(eep_dp, (signed short) (act_sym!=1), ctx->mb_type_contexts[1] + act_ctx );
      }

      mode_sym = mode_sym % 4; // coding of I pred-mode: 0,1,2,3
      act_ctx  = 10;
      act_sym  = mode_sym/2;
      biari_encode_symbol(eep_dp, (signed short) act_sym, ctx->mb_type_contexts[1] + act_ctx );
      act_sym  = mode_sym%2;
      biari_encode_symbol(eep_dp, (signed short) act_sym, ctx->mb_type_contexts[1] + act_ctx );
    }
  }
}


/*!
 ***************************************************************************
 * \brief
 *    This function is used to arithmetically encode the 8x8 block
 *    type info
 ***************************************************************************
 */
void writeB8_typeInfo_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  int act_ctx;
  int act_sym;
  signed short csym;
  int bframe=(img->type==B_SLICE);

  MotionInfoContexts *ctx = (img->currentSlice)->mot_ctx;

  act_sym = se->value1;
  act_ctx = 0;

  if (!bframe)  
  {
    switch (act_sym)
    {
    case 0:
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[0][1]);
      break;
    case 1:
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[0][1]);
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[0][3]);
      break;
    case 2:
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[0][1]);
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[0][3]);
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[0][4]);
      break;
    case 3:
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[0][1]);
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[0][3]);
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[0][4]);
      break;
    }
  }
  else //===== B-FRAME =====
  {
    if (act_sym==0)
    {
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[1][0]);
      return;
    }
    else
    {
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][0]);
      act_sym--;
    }
    if (act_sym<2)
    {
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[1][1]);
      biari_encode_symbol (eep_dp, (signed short) (act_sym!=0), &ctx->b8_type_contexts[1][3]);
    }
    else if (act_sym<6)
    {
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][1]);
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[1][2]);
      csym=(((act_sym-2)>>1)&0x01) != 0;
      biari_encode_symbol (eep_dp, csym, &ctx->b8_type_contexts[1][3]);
      csym=((act_sym-2)&0x01) != 0;
      biari_encode_symbol (eep_dp, csym, &ctx->b8_type_contexts[1][3]);
    }
    else
    {
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][1]);
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][2]);
      csym=(((act_sym-6)>>2)&0x01);
      if (csym)
      {
        biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][3]);
        csym=((act_sym-6)&0x01) != 0;
        biari_encode_symbol (eep_dp, csym, &ctx->b8_type_contexts[1][3]);
      }
      else
      {
        biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[1][3]);
        csym=(((act_sym-6)>>1)&0x01) != 0;
        biari_encode_symbol (eep_dp, csym, &ctx->b8_type_contexts[1][3]);
        csym=((act_sym-6)&0x01) != 0;
        biari_encode_symbol (eep_dp, csym, &ctx->b8_type_contexts[1][3]);
      }
    }
  }
}


/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode a pair of
 *    intra prediction modes of a given MB.
 ****************************************************************************
 */
void writeIntraPredMode_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  TextureInfoContexts *ctx = img->currentSlice->tex_ctx;

  // use_most_probable_mode
  if (se->value1 == -1)
    biari_encode_symbol(eep_dp, 1, ctx->ipr_contexts);
  else
  {
    biari_encode_symbol(eep_dp, 0, ctx->ipr_contexts);
        
    // remaining_mode_selector
    biari_encode_symbol(eep_dp,(signed short)( se->value1 & 0x1    ), ctx->ipr_contexts+1);
    biari_encode_symbol(eep_dp,(signed short)((se->value1 & 0x2)>>1), ctx->ipr_contexts+1);
    biari_encode_symbol(eep_dp,(signed short)((se->value1 & 0x4)>>2), ctx->ipr_contexts+1);
  }
}
/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode the reference
 *    parameter of a given MB.
 ****************************************************************************
 */
void writeRefFrame_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  MotionInfoContexts  *ctx    = img->currentSlice->mot_ctx;
  Macroblock          *currMB = &img->mb_data[img->current_mb_nr];
  int                 addctx  = 0;

  int   a, b;
  int   act_ctx;
  int   act_sym;
  char** refframe_array = enc_picture->ref_idx[se->value2];

  int bslice = (img->type==B_SLICE);

  int   b8a, b8b;

  PixelPos block_a, block_b;
  
  getLuma4x4Neighbour(img->current_mb_nr, img->subblock_x, img->subblock_y, -1,  0, &block_a);
  getLuma4x4Neighbour(img->current_mb_nr, img->subblock_x, img->subblock_y,  0, -1, &block_b);

  b8a=((block_a.x >> 1) & 0x01)+2*((block_a.y >> 1) & 0x01);
  b8b=((block_b.x >> 1) & 0x01)+2*((block_b.y >> 1) & 0x01);

  
  if (!block_b.available)
    b=0;
  //else if (IS_DIRECT(&img->mb_data[block_b.mb_addr]) || (img->mb_data[block_b.mb_addr].b8mode[b8b]==0 && bslice))
  else if ((IS_DIRECT(&img->mb_data[block_b.mb_addr]) && !giRDOpt_B8OnlyFlag) || (img->mb_data[block_b.mb_addr].b8mode[b8b]==0 && bslice))
    b=0;
  else
  {
    if (img->MbaffFrameFlag && (currMB->mb_field == 0) && (img->mb_data[block_b.mb_addr].mb_field == 1))
      b = (refframe_array[block_b.pos_y][block_b.pos_x] > 1 ? 1 : 0);
    else
      b = (refframe_array[block_b.pos_y][block_b.pos_x] > 0 ? 1 : 0);
  }

  if (!block_a.available)
    a=0;
  // else if (IS_DIRECT(&img->mb_data[block_a.mb_addr]) || (img->mb_data[block_a.mb_addr].b8mode[b8a]==0 && bslice))
  else if ((IS_DIRECT(&img->mb_data[block_a.mb_addr]) && !giRDOpt_B8OnlyFlag) || (img->mb_data[block_a.mb_addr].b8mode[b8a]==0 && bslice))
    a=0;
  else 
  {
    if (img->MbaffFrameFlag && (currMB->mb_field == 0) && (img->mb_data[block_a.mb_addr].mb_field == 1))
      a = (refframe_array[block_a.pos_y][block_a.pos_x] > 1 ? 1 : 0);
    else
      a = (refframe_array[block_a.pos_y][block_a.pos_x] > 0 ? 1 : 0);
  }

  act_ctx     = a + 2*b; 
  se->context = act_ctx; // store context
  act_sym     = se->value1;

  if (act_sym==0)
  {
    biari_encode_symbol(eep_dp, 0, ctx->ref_no_contexts[addctx] + act_ctx );
  }
  else
  {
    biari_encode_symbol(eep_dp, 1, ctx->ref_no_contexts[addctx] + act_ctx);
    act_sym--;
    act_ctx=4;
    unary_bin_encode(eep_dp, act_sym,ctx->ref_no_contexts[addctx]+act_ctx,1);
  }
}

/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode the coded
 *    block pattern of a given delta quant.
 ****************************************************************************
 */
void writeDquant_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  MotionInfoContexts *ctx = img->currentSlice->mot_ctx;

  int act_ctx;
  int act_sym;
  int dquant = se->value1;
  int sign=0;

  Macroblock     *currMB       = &img->mb_data[img->current_mb_nr];

  last_dquant=currMB->prev_delta_qp;

  if (dquant <= 0)
    sign = 1;
  act_sym = absm(dquant) << 1;

  act_sym += sign;
  act_sym --;

  act_ctx = ( (last_dquant != 0) ? 1 : 0);

  if (act_sym==0)
  {
    biari_encode_symbol(eep_dp, 0, ctx->delta_qp_contexts + act_ctx );
  }
  else
  {
    biari_encode_symbol(eep_dp, 1, ctx->delta_qp_contexts + act_ctx);
    act_ctx=2;
    act_sym--;
    unary_bin_encode(eep_dp, act_sym,ctx->delta_qp_contexts+act_ctx,1);
  }
}

/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode the motion
 *    vector data of a B-frame MB.
 ****************************************************************************
 */
void writeMVD_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  int i = img->subblock_x;
  int j = img->subblock_y;
  int a, b;
  int act_ctx;
  int act_sym;
  int mv_pred_res;
  int mv_local_err;
  int mv_sign;
  int list_idx = se->value2 & 0x01;
  int k = (se->value2>>1); // MVD component

  PixelPos block_a, block_b;

  MotionInfoContexts  *ctx    = img->currentSlice->mot_ctx;
  Macroblock          *currMB = &img->mb_data[img->current_mb_nr];

  getLuma4x4Neighbour(img->current_mb_nr, i, j, -1,  0, &block_a);
  getLuma4x4Neighbour(img->current_mb_nr, i, j,  0, -1, &block_b);

  if (block_b.available)
  {
    b = absm(img->mb_data[block_b.mb_addr].mvd[list_idx][block_b.y][block_b.x][k]);
    if (img->MbaffFrameFlag && (k==1)) 
    {
      if ((currMB->mb_field==0) && (img->mb_data[block_b.mb_addr].mb_field==1))
        b *= 2;
      else if ((currMB->mb_field==1) && (img->mb_data[block_b.mb_addr].mb_field==0))
        b /= 2;
    }
  }
  else
    b=0;
          
  if (block_a.available)
  {
    a = absm(img->mb_data[block_a.mb_addr].mvd[list_idx][block_a.y][block_a.x][k]);
    if (img->MbaffFrameFlag && (k==1)) 
    {
      if ((currMB->mb_field==0) && (img->mb_data[block_a.mb_addr].mb_field==1))
        a *= 2;
      else if ((currMB->mb_field==1) && (img->mb_data[block_a.mb_addr].mb_field==0))
        a /= 2;
    }
  }
  else
    a = 0;

  if ((mv_local_err=a+b)<3)
    act_ctx = 5*k;
  else
  {
    if (mv_local_err>32)
      act_ctx=5*k+3;
    else
      act_ctx=5*k+2;
  }

  mv_pred_res = se->value1;
  se->context = act_ctx;

  act_sym = absm(mv_pred_res);

  if (act_sym == 0)
    biari_encode_symbol(eep_dp, 0, &ctx->mv_res_contexts[0][act_ctx] );
  else
  {
    biari_encode_symbol(eep_dp, 1, &ctx->mv_res_contexts[0][act_ctx] );
    act_sym--;
    act_ctx=5*k;
    unary_exp_golomb_mv_encode(eep_dp,act_sym,ctx->mv_res_contexts[1]+act_ctx,3);
    mv_sign = (mv_pred_res<0) ? 1: 0;
    biari_encode_symbol_eq_prob(eep_dp, (signed short) mv_sign);
  }
}


/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode the chroma
 *    intra prediction mode of an 8x8 block
 ****************************************************************************
 */
void writeCIPredMode_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  TextureInfoContexts *ctx     = img->currentSlice->tex_ctx;
  Macroblock          *currMB  = &img->mb_data[img->current_mb_nr];
  int                 act_ctx,a,b;
  int                 act_sym  = se->value1;

  if (currMB->mb_available_up == NULL) b = 0;
  else  b = ( ((currMB->mb_available_up)->c_ipred_mode != 0) ? 1 : 0);

  if (currMB->mb_available_left == NULL) a = 0;
  else  a = ( ((currMB->mb_available_left)->c_ipred_mode != 0) ? 1 : 0);

  act_ctx = a+b;

  if (act_sym==0) 
    biari_encode_symbol(eep_dp, 0, ctx->cipr_contexts + act_ctx );
  else
  {
    biari_encode_symbol(eep_dp, 1, ctx->cipr_contexts + act_ctx );
    unary_bin_max_encode(eep_dp,(unsigned int) (act_sym-1),ctx->cipr_contexts+3,0,2);
  }
}


/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode the coded
 *    block pattern of an 8x8 block
 ****************************************************************************
 */
void writeCBP_BIT_CABAC (int b8, int bit, int cbp, Macroblock* currMB, int inter, EncodingEnvironmentPtr eep_dp)
{
  PixelPos block_a;
  int a, b;
  
  int mb_x=(b8 & 0x01)<<1;
  int mb_y=(b8 >> 1)<<1;

  if (mb_y == 0)
  {
    if (currMB->mb_available_up == NULL)
      b = 0;
    else
    {
      if((currMB->mb_available_up)->mb_type==IPCM)
        b=0;
      else
        b = (( ((currMB->mb_available_up)->cbp & (1<<(2+(mb_x>>1)))) == 0) ? 1 : 0);   //VG-ADD
    }    
    
  }
  else
    b = ( ((cbp & (1<<(mb_x/2))) == 0) ? 1: 0);
  
  if (mb_x == 0)
  {
    getLuma4x4Neighbour(img->current_mb_nr, mb_x, mb_y, -1, 0, &block_a);
    if (block_a.available)
    {
      {
        if(img->mb_data[block_a.mb_addr].mb_type==IPCM)
          a=0;
        else
          a = (( (img->mb_data[block_a.mb_addr].cbp & (1<<(2*(block_a.y>>1)+1))) == 0) ? 1 : 0); //VG-ADD
      }
      
    }
    else
      a=0;
  }
  else
    a = ( ((cbp & (1<<mb_y)) == 0) ? 1: 0);
  
  //===== WRITE BIT =====
  biari_encode_symbol (eep_dp, (signed short) bit,
    img->currentSlice->tex_ctx->cbp_contexts[0] + a+2*b);
}

/*!
****************************************************************************
* \brief
*    This function is used to arithmetically encode the coded
*    block pattern of a macroblock
****************************************************************************
*/
void writeCBP_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  TextureInfoContexts *ctx = img->currentSlice->tex_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  
  int a, b;
  int curr_cbp_ctx, curr_cbp_idx;
  int cbp = se->value1; // symbol to encode
  int cbp_bit;
  int b8;
  
#ifdef MB32X32
  int b8step = 1;
  int b8num = 4;
  
  if(input->UseExtMB > 0)
  {
    if ((currMB->mb_type == 1 || se->value2) && currMB->luma_transform_size_8x8_flag == 2) // se->value2 = mb32;
      b8step = 4;
    else if (currMB->mb_type == 2 && currMB->luma_transform_size_8x8_flag == 2)
      b8step = 2;
    else if (currMB->mb_type == 3 && currMB->luma_transform_size_8x8_flag == 2)
    {
      b8step = 1;   
      b8num = 2;
    }
  }
#endif



#ifdef MB32X32
  for (b8=0; b8<b8num; b8+=b8step)
#else
  for (b8=0; b8<4; b8++)
#endif
  {
    curr_cbp_idx = (currMB->b8mode[b8] == IBLOCK ? 0 : 1);
    writeCBP_BIT_CABAC (b8, cbp&(1<<b8), cbp, currMB, curr_cbp_idx, eep_dp);
  }

  if (img->yuv_format != YUV400)
  {
    // coding of chroma part
    b = 0;
    if (currMB->mb_available_up != NULL)
    {
      if((currMB->mb_available_up)->mb_type==IPCM)
        b=1;
      else
        b = ((currMB->mb_available_up)->cbp > 15) ? 1 : 0;
    }
    
    
    a = 0;
    if (currMB->mb_available_left != NULL)
    {
      if((currMB->mb_available_left)->mb_type==IPCM)
        a=1;
      else
        a = ((currMB->mb_available_left)->cbp > 15) ? 1 : 0;
    }
    
    curr_cbp_ctx = a+2*b;
    cbp_bit = (cbp > 15 ) ? 1 : 0;
    biari_encode_symbol(eep_dp, (signed short) cbp_bit, ctx->cbp_contexts[1] + curr_cbp_ctx );
    
    if (cbp > 15)
    {
      b = 0;
      if (currMB->mb_available_up != NULL)
      {
        if((currMB->mb_available_up)->mb_type==IPCM)
          b=1;
        else
          if ((currMB->mb_available_up)->cbp > 15)
            b = (( ((currMB->mb_available_up)->cbp >> 4) == 2) ? 1 : 0);
      }
      
      
      a = 0;
      if (currMB->mb_available_left != NULL)
      {
        if((currMB->mb_available_left)->mb_type==IPCM)
          a=1;
        else
          if ((currMB->mb_available_left)->cbp > 15)
            a = (( ((currMB->mb_available_left)->cbp >> 4) == 2) ? 1 : 0);
      }
      
      curr_cbp_ctx = a+2*b;
      cbp_bit = ((cbp>>4) == 2) ? 1 : 0;
      biari_encode_symbol(eep_dp, (signed short) cbp_bit, ctx->cbp_contexts[2] + curr_cbp_ctx );
    }
  }
}

#ifdef USE_INTRA_MDDT
#ifdef MB32X32
static const int maxpos       [] = {16, 15, 64, 32, 32, 16,  4, 15,  8, 16, 256, 256, 128, 128};
static const int c1isdc       [] = { 1,  0,  1,  1,  1,  1,  1,  0,  1,  1,   1,   1,   1,   1};
#else
static const int maxpos       [] = {16, 15, 64, 32, 32, 16,  4, 15,  8, 16, 256};
static const int c1isdc       [] = { 1,  0,  1,  1,  1,  1,  1,  0,  1,  1,   1};
#endif
#else
static const int maxpos       [] = {16, 15, 64, 32, 32, 16,  4, 15,  8, 16};
static const int c1isdc       [] = { 1,  0,  1,  1,  1,  1,  1,  0,  1,  1};
#endif

#ifdef USE_INTRA_MDDT
/*! 
*************************************************************************************
* \brief
*   Set coded-block-pattern bits for INTRA16x16 MB 
*
* \para set_CBP_block_bit_16x16()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
void set_CBP_block_bit_16x16 (Macroblock* currMB)
{
  int bit;

  for(bit = 1; bit < 17; bit++)
  {
    currMB->cbp_bits   |= ((int64)1<< bit   );
  }
}
#endif 

/*!
 ****************************************************************************
 * \brief
 *    Write CBP4-BIT
 ****************************************************************************
 */
void write_and_store_CBP_block_bit (Macroblock* currMB, EncodingEnvironmentPtr eep_dp, int type, int cbp_bit)
{
  int y_ac        = (type==LUMA_16AC || type==LUMA_8x8 || type==LUMA_8x4 || type==LUMA_4x8 || type==LUMA_4x4);
  int y_dc        = (type==LUMA_16DC);
  int u_ac        = (type==CHROMA_AC && !img->is_v_block);
  int v_ac        = (type==CHROMA_AC &&  img->is_v_block);
  int chroma_dc   = (type==CHROMA_DC || type==CHROMA_DC_2x4 || type==CHROMA_DC_4x4);
  int u_dc        = (chroma_dc && !img->is_v_block);
  int v_dc        = (chroma_dc &&  img->is_v_block);
  int j           = (y_ac || u_ac || v_ac ? img->subblock_y : 0);
  int i           = (y_ac || u_ac || v_ac ? img->subblock_x : 0);
  int bit         = (y_dc ? 0 : y_ac ? 1 : u_dc ? 17 : v_dc ? 18 : u_ac ? 19 : 23);
  int default_bit = (img->is_intra_block ? 1 : 0);
  int upper_bit   = default_bit;
  int left_bit    = default_bit;
  int ctx;

  int bit_pos_a   = 0;
  int bit_pos_b   = 0;

  PixelPos block_a, block_b;

  if (y_ac || y_dc)
  {
    getLuma4x4Neighbour(img->current_mb_nr, i, j, -1,  0, &block_a);
    getLuma4x4Neighbour(img->current_mb_nr, i, j,  0, -1, &block_b);
    if (y_ac)
    {
      if (block_a.available)
        bit_pos_a = 4*block_a.y + block_a.x;
      if (block_b.available)
        bit_pos_b = 4*block_b.y + block_b.x;
    }
  }
  else
  {
    getChroma4x4Neighbour(img->current_mb_nr, i, j, -1,  0, &block_a);
    getChroma4x4Neighbour(img->current_mb_nr, i, j,  0, -1, &block_b);
    if (u_ac||v_ac)
    {
      if (block_a.available)
        bit_pos_a = 4*block_a.y + block_a.x;
      if (block_b.available)
        bit_pos_b = 4*block_b.y + block_b.x;
    }
  }

  bit = (y_dc ? 0 : y_ac ? 1+4*j+i : u_dc ? 17 : v_dc ? 18 : u_ac ? 19+4*j+i : 35+4*j+i);
  //--- set bits for current block ---
  if (cbp_bit)
  {
    if (type==LUMA_8x8)
    {
      currMB->cbp_bits   |= ((int64)1<< bit   );
      currMB->cbp_bits   |= ((int64)1<<(bit+1));
      currMB->cbp_bits   |= ((int64)1<<(bit+4));
      currMB->cbp_bits   |= ((int64)1<<(bit+5));
    }
    else if (type==LUMA_8x4)
    {
      currMB->cbp_bits   |= ((int64)1<< bit   );
      currMB->cbp_bits   |= ((int64)1<<(bit+1));
    }
    else if (type==LUMA_4x8)
    {
      currMB->cbp_bits   |= ((int64)1<< bit   );
      currMB->cbp_bits   |= ((int64)1<<(bit+4));
    }
    else
    {
      currMB->cbp_bits   |= ((int64)1<<bit);
    }
  }

  bit = (y_dc ? 0 : y_ac ? 1 : u_dc ? 17 : v_dc ? 18 : u_ac ? 19 : 35);
  if (type!=LUMA_8x8)
  {
    if (block_b.available)
    {
      if(img->mb_data[block_b.mb_addr].mb_type==IPCM)
        upper_bit=1;
      else
        upper_bit = BIT_SET(img->mb_data[block_b.mb_addr].cbp_bits,bit+bit_pos_b);
    }

    
    if (block_a.available)
    {
      if(img->mb_data[block_a.mb_addr].mb_type==IPCM)
        left_bit=1;
      else
        left_bit = BIT_SET(img->mb_data[block_a.mb_addr].cbp_bits,bit+bit_pos_a);
    }

    ctx = 2*upper_bit+left_bit;

    //===== encode symbol =====
    biari_encode_symbol (eep_dp, (short)cbp_bit, img->currentSlice->tex_ctx->bcbp_contexts[type2ctx_bcbp[type]]+ctx);
  }
}


#ifdef RDO_Q
/*!
 ****************************************************************************
 * \brief
 *    estimate CABAC CBP bits
 ****************************************************************************
 */
int est_write_and_store_CBP_block_bit(Macroblock* currMB, int type) // marta - CBP
{
  int j           = img->subblock_y;
  int i           = img->subblock_x;
  int bit, default_bit = (IS_INTRA(currMB) ? 1 : 0);
  int upper_bit   = default_bit;
  int left_bit    = default_bit;
  int ctx, estBits;

  int bit_pos_a   = 0;
  int bit_pos_b   = 0;

  if (type!=LUMA_8x8)
  {
    PixelPos block_a, block_b;

    getLuma4x4Neighbour(img->current_mb_nr, i, j, -1,  0, &block_a);
    getLuma4x4Neighbour(img->current_mb_nr, i, j,  0, -1, &block_b);

    if (block_a.available)
      bit_pos_a = 4*block_a.y + block_a.x;
    if (block_b.available)
      bit_pos_b = 4*block_b.y + block_b.x;

    bit = 1; // 4x4: bit=1

    if (block_b.available)
    {
      if(img->mb_data[block_b.mb_addr].mb_type==IPCM)
        upper_bit=1;
      else
        upper_bit = BIT_SET(img->mb_data[block_b.mb_addr].cbp_bits,bit+bit_pos_b);
    }


    if (block_a.available)
    {
      if(img->mb_data[block_a.mb_addr].mb_type==IPCM)
        left_bit=1;
      else
        left_bit = BIT_SET(img->mb_data[block_a.mb_addr].cbp_bits,bit+bit_pos_a);
    }

    ctx = 2*upper_bit+left_bit;
    //===== encode symbol =====
    estBits=estBitsCabac[type].blockCbpBits[ctx][0]-estBitsCabac[type].blockCbpBits[ctx][1];
  }
  else
  {
    estBits=0;
  }
  return(estBits);
}
/*!
 ****************************************************************************
 * \brief
 *    estimate bit cost for each CBP bit
 ****************************************************************************
 */
void est_CBP_block_bit (Macroblock* currMB, EncodingEnvironmentPtr eep_dp, int type)
{
  int ctx;
  short cbp_bit;

  for (ctx=0; ctx<=3; ctx++)
  {
    cbp_bit=0;
    estBitsCabac[type].blockCbpBits[ctx][cbp_bit]=biari_no_bits(eep_dp, cbp_bit, img->currentSlice->tex_ctx->bcbp_contexts[type2ctx_bcbp[type]]+ctx);

    cbp_bit=1;
    estBitsCabac[type].blockCbpBits[ctx][cbp_bit]=biari_no_bits(eep_dp, cbp_bit, img->currentSlice->tex_ctx->bcbp_contexts[type2ctx_bcbp[type]]+ctx);
  }
}
#endif


//===== position -> ctx for MAP =====
//--- zig-zag scan ----
static const int  pos2ctx_map8x8 [] = { 0,  1,  2,  3,  4,  5,  5,  4,  4,  3,  3,  4,  4,  4,  5,  5,
                                        4,  4,  4,  4,  3,  3,  6,  7,  7,  7,  8,  9, 10,  9,  8,  7,
                                        7,  6, 11, 12, 13, 11,  6,  7,  8,  9, 14, 10,  9,  8,  6, 11,
#ifdef RDO_Q
                                       12, 13, 11,  6,  9, 14, 10,  9, 11, 12, 13, 11 ,14, 10, 12, 15}; // 15 CTX
#else
                                       12, 13, 11,  6,  9, 14, 10,  9, 11, 12, 13, 11 ,14, 10, 12, 14}; // 15 CTX
#endif
static const int  pos2ctx_map8x4 [] = { 0,  1,  2,  3,  4,  5,  7,  8,  9, 10, 11,  9,  8,  6,  7,  8,
#ifdef RDO_Q
                                        9, 10, 11,  9,  8,  6, 12,  8,  9, 10, 11,  9, 13, 13, 14, 15}; // 15 CTX
#else
                                        9, 10, 11,  9,  8,  6, 12,  8,  9, 10, 11,  9, 13, 13, 14, 14}; // 15 CTX
#endif
#ifdef RDO_Q
static const int  pos2ctx_map4x4 [] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15}; // 15 CTX
#else
static const int  pos2ctx_map4x4 [] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 14}; // 15 CTX
#endif
static const int  pos2ctx_map2x4c[] = { 0,  0,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2}; // 15 CTX
static const int  pos2ctx_map4x4c[] = { 0,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2}; // 15 CTX
#ifndef USE_INTRA_MDDT
static const int* pos2ctx_map    [] = {pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map8x8, pos2ctx_map8x4,
                                       pos2ctx_map8x4, pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map4x4,
                                       pos2ctx_map2x4c, pos2ctx_map4x4c};
#else
static const int  pos2ctx_map16x16[]= { 
 0,  1,  2,  2,  1,  1,  1,  1,  2,  2,  3,  2,  4,  1,  5,  5, 
 5,  4,  4,  3,  3,  3,  3,  4,  4,  4,  5,  5,  5,  5,  4,  4, 
 4,  4,  3,  3,  3,  3,  4,  4,  4,  4,  4,  5,  5,  5,  5,  4, 
 4,  4,  4,  4,  4,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4, 
 5,  5,  5,  5,  4,  4,  4,  4,  4,  4,  4,  4,  3,  3,  3,  3, 
 6,  4,  7,  4,  7,  4,  7,  4,  8,  5,  9,  9,  9,  8,  8,  7, 
 7,  7,  7,  7,  7,  6,  6,  3,  3, 12,  3, 11,  6,  6,  7,  7, 
 7,  7,  7,  8,  8,  9,  9, 10, 10, 10,  9,  9,  8,  8,  7,  7, 
 7,  7,  6,  6, 11, 11, 12, 12, 12, 13, 11, 11,  6,  6,  7,  7, 
 7,  8,  8,  9,  9, 14, 10, 14, 14,  9,  9,  8,  8,  7,  7,  6, 
 6, 11, 11, 13, 13, 13, 12, 11, 11,  6,  6,  7,  8,  8,  9,  9, 
10, 14, 10, 10,  9,  9,  8,  8,  6,  6, 11, 11, 12, 12, 12, 13, 
11, 11,  6,  6,  8,  9,  9, 14, 10, 14, 14,  9,  9,  6,  6, 11, 
11, 13, 13, 13, 12, 11, 11,  6,  9,  9, 10, 14, 10, 10,  9,  9, 
11, 11, 12, 12, 12, 13, 11, 11,  9, 14, 10, 14, 14, 11, 11, 13, 
13, 13, 12, 11, 10, 14, 10, 10, 12, 12, 12, 14, 10, 14, 14, 15, 
}; // 15   CTX
#ifdef MB32X32
static const int  pos2ctx_map16x8[]= { 
   0, 2, 1, 1, 2, 3, 3, 3, 4, 1, 5, 4, 4, 3, 3, 3, 
   3, 4, 4, 4, 5, 5, 4, 4, 4, 4, 3, 3,12, 3, 6, 4, 
   4, 4, 4, 5, 5, 4, 4, 4, 7, 6,11,12,13,11, 6, 7, 
   7, 4, 4, 5, 5, 4, 7, 7, 7, 6,11,13,12,11, 6, 7, 
   7, 7, 8, 5, 9, 8, 8, 7, 7, 6,11,12,13,11, 6, 7, 
   8, 8, 9, 9,10, 9, 9, 8, 8, 6,11,13,12,11, 6, 8, 
   9, 9,14,10,14,10, 9, 9, 6,11,12,13,11, 9, 9,14, 
  10,14,10, 9,11,13,12,11,14,10,14,10,12,10,10,15
}; // 15 CTX

static const int* pos2ctx_map    [] = {pos2ctx_map4x4,  pos2ctx_map4x4, pos2ctx_map8x8, pos2ctx_map8x4,
                                       pos2ctx_map8x4,  pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map4x4,
                                       pos2ctx_map2x4c, pos2ctx_map4x4c,pos2ctx_map16x16, 
                                       pos2ctx_map16x16, pos2ctx_map16x8, pos2ctx_map16x8};

#else
static const int* pos2ctx_map    [] = {pos2ctx_map4x4,  pos2ctx_map4x4, pos2ctx_map8x8, pos2ctx_map8x4,
                                       pos2ctx_map8x4,  pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map4x4,
                                       pos2ctx_map2x4c, pos2ctx_map4x4c,pos2ctx_map16x16};
#endif
#endif

//--- interlace scan ----
//Taken from ABT
static const int  pos2ctx_map8x8i[] = { 0,  1,  1,  2,  2,  3,  3,  4,  5,  6,  7,  7,  7,  8,  4,  5,
                                        6,  9, 10, 10,  8, 11, 12, 11,  9,  9, 10, 10,  8, 11, 12, 11,
                                        9,  9, 10, 10,  8, 11, 12, 11,  9,  9, 10, 10,  8, 13, 13,  9,
                                        9, 10, 10,  8, 13, 13,  9,  9, 10, 10, 14, 14, 14, 14, 14, 14}; // 15 CTX

static const int  pos2ctx_map8x4i[] = { 0,  1,  2,  3,  4,  5,  6,  3,  4,  5,  6,  3,  4,  7,  6,  8,
                                        9,  7,  6,  8,  9, 10, 11, 12, 12, 10, 11, 13, 13, 14, 14, 14}; // 15 CTX
static const int  pos2ctx_map4x8i[] = { 0,  1,  1,  1,  2,  3,  3,  4,  4,  4,  5,  6,  2,  7,  7,  8,
                                        8,  8,  5,  6,  9, 10, 10, 11, 11, 11, 12, 13, 13, 14, 14, 14}; // 15 CTX
static const int* pos2ctx_map_int[] = {pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map8x8i,pos2ctx_map8x4i,
                                       pos2ctx_map4x8i,pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map4x4,
                                       pos2ctx_map2x4c, pos2ctx_map4x4c};


//===== position -> ctx for LAST =====
static const int  pos2ctx_last8x8 [] = { 0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                                         2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
                                         3,  3,  3,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,
#ifdef RDO_Q
                                         5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  15}; //  9 CTX
#else
                                         5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  8}; //  9 CTX
#endif
static const int  pos2ctx_last8x4 [] = { 0,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,
#ifdef RDO_Q
                                         3,  3,  3,  3,  4,  4,  4,  4,  5,  5,  6,  6,  7,  7,  8,  15}; //  9 CTX
#else
                                         3,  3,  3,  3,  4,  4,  4,  4,  5,  5,  6,  6,  7,  7,  8,  8}; //  9 CTX
#endif
static const int  pos2ctx_last4x4 [] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15}; // 15 CTX
static const int  pos2ctx_last2x4c[] = { 0,  0,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2}; // 15 CTX
static const int  pos2ctx_last4x4c[] = { 0,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2}; // 15 CTX
#ifndef USE_INTRA_MDDT
static const int* pos2ctx_last    [] = {pos2ctx_last4x4, pos2ctx_last4x4, pos2ctx_last8x8, pos2ctx_last8x4,
                                        pos2ctx_last8x4, pos2ctx_last4x4, pos2ctx_last4x4, pos2ctx_last4x4,
                                        pos2ctx_last2x4c, pos2ctx_last4x4c};
#else
static const int  pos2ctx_last16x16[]= {  
    0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, 
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, 
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, 
    1,  1,  1,  1,  1,  1,  1,  2,  1,  2,  1,  2,  1,  2,  1,  2, 
    1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 
    2,  2,  2,  2,  2,  2,  2,  2,  2,  1,  2,  2,  2,  2,  2,  2, 
    2,  2,  2,  2,  2,  2,  2,  2,  2,  3,  2,  3,  2,  3,  2,  3, 
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 
    3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3, 
    2,  4,  2,  4,  2,  4,  2,  4,  4,  4,  4,  4,  4,  3,  3,  3, 
    3,  3,  3,  3,  3,  3,  5,  3,  4,  3,  4,  3,  4,  4,  4,  4, 
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  5,  5,  5,  5, 
    4,  5,  4,  5,  4,  6,  4,  6,  4,  6,  6,  6,  6,  5,  5,  5, 
    5,  5,  5,  5,  7,  5,  7,  5,  6,  6,  6,  6,  6,  6,  6,  6, 
    7,  7,  7,  7,  7,  7,  7,  7,  6,  8,  6,  8,  8,  7,  7,  7, 
    7,  7,  8,  7,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8, 15, 
}; //  9 CTX

#ifdef MB32X32

static const int  pos2ctx_last16x8[]= {  
    0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2, 
    1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  3,  2,  2,  2, 
    2,  1,  1,  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3,  3,  2, 
    2,  2,  2,  1,  1,  2,  2,  2,  3,  3,  3,  3,  5,  3,  3,  3, 
    2,  2,  2,  1,  2,  2,  2,  2,  3,  3,  4,  5,  5,  4,  4,  3, 
    4,  2,  2,  2,  2,  2,  4,  4,  4,  4,  5,  5,  7,  5,  5,  4, 
    4,  4,  4,  2,  4,  4,  4,  6,  5,  7,  7,  7,  7,  6,  6,  6, 
    4,  6,  6,  6,  7,  7,  8,  7,  8,  6,  8,  8,  8,  8,  8, 15, 
}; 

static const int* pos2ctx_last    [] = {pos2ctx_last4x4,  pos2ctx_last4x4,  pos2ctx_last8x8, pos2ctx_last8x4,
                                        pos2ctx_last8x4,  pos2ctx_last4x4,  pos2ctx_last4x4, pos2ctx_last4x4,
                                        pos2ctx_last2x4c, pos2ctx_last4x4c, pos2ctx_last16x16, 
                                        pos2ctx_last16x16, pos2ctx_last16x8, pos2ctx_last16x8};

#else
static const int* pos2ctx_last    [] = {pos2ctx_last4x4,  pos2ctx_last4x4,  pos2ctx_last8x8, pos2ctx_last8x4,
                                        pos2ctx_last8x4,  pos2ctx_last4x4,  pos2ctx_last4x4, pos2ctx_last4x4,
                                        pos2ctx_last2x4c, pos2ctx_last4x4c, pos2ctx_last16x16};
#endif
#endif




/*!
****************************************************************************
* \brief
*    Write Significance MAP
****************************************************************************
*/
void write_significance_map (Macroblock* currMB, EncodingEnvironmentPtr eep_dp, int type, int coeff[], int coeff_ctr)
{
  int   k;
  unsigned short sig, last;
  int   k0      = 0;
  int   k1      = maxpos[type]-1;
  
  int               fld       = ( img->structure!=FRAME || currMB->mb_field );
  BiContextTypePtr  map_ctx   = ( fld ? img->currentSlice->tex_ctx->fld_map_contexts[type2ctx_map [type]]
    : img->currentSlice->tex_ctx->map_contexts[type2ctx_map [type]] );
  BiContextTypePtr  last_ctx  = ( fld ? img->currentSlice->tex_ctx->fld_last_contexts[type2ctx_last[type]]
    : img->currentSlice->tex_ctx->last_contexts[type2ctx_last[type]] );
  
  if (!c1isdc[type])
  {
    k0++; k1++; coeff--;
  }
  
  if (!fld)
  {
    for (k=k0; k<k1; k++) // if last coeff is reached, it has to be significant
    {
      sig   = (coeff[k] != 0);      
      biari_encode_symbol  (eep_dp, sig,  map_ctx+pos2ctx_map     [type][k]);
      if (sig)
      {
        last = (--coeff_ctr == 0);
        
        biari_encode_symbol(eep_dp, last, last_ctx+pos2ctx_last[type][k]);
        if (last) return;
      }
    }
    return;
  }
  else
  {
    for (k=k0; k<k1; k++) // if last coeff is reached, it has to be significant
    {
      sig   = (coeff[k] != 0);
      
      biari_encode_symbol  (eep_dp, sig,  map_ctx+pos2ctx_map_int [type][k]);
      if (sig)
      {
        last = (--coeff_ctr == 0);
        
        biari_encode_symbol(eep_dp, last, last_ctx+pos2ctx_last[type][k]);
        if (last) return;
      }
    }
  }
}

#ifdef RDO_Q
/*!
 ****************************************************************************
 * \brief
 *    estimate CABAC bit cost for significant coefficient map
 ****************************************************************************
 */
void est_significance_map(Macroblock* currMB, EncodingEnvironmentPtr eep_dp, int type)
{
  int   k;
  unsigned short sig, last;
  int   k1      = maxpos[type]-1;
  int               fld       = ( img->structure!=FRAME || currMB->mb_field );
  BiContextTypePtr  map_ctx   = ( fld ? img->currentSlice->tex_ctx->fld_map_contexts[type2ctx_map [type]]
    : img->currentSlice->tex_ctx->map_contexts[type2ctx_map [type]] );
  BiContextTypePtr  last_ctx  = ( fld ? img->currentSlice->tex_ctx->fld_last_contexts[type2ctx_last[type]]
    : img->currentSlice->tex_ctx->last_contexts[type2ctx_last[type]] );
  
  for (k=0; k<k1; k++) // if last coeff is reached, it has to be significant
  {
    sig   = 0;     
    estBitsCabac[type].significantBits[pos2ctx_map[type][k]][sig]=biari_no_bits  (eep_dp, sig,  map_ctx+pos2ctx_map     [type][k]);

    sig   = 1;     
    estBitsCabac[type].significantBits[pos2ctx_map[type][k]][sig]=biari_no_bits  (eep_dp, sig,  map_ctx+pos2ctx_map     [type][k]);

    last=0;
    estBitsCabac[type].lastBits[pos2ctx_last[type][k]][last]=biari_no_bits(eep_dp, last, last_ctx+pos2ctx_last[type][k]);

    last=1;
    estBitsCabac[type].lastBits[pos2ctx_last[type][k]][last]=biari_no_bits(eep_dp, last, last_ctx+pos2ctx_last[type][k]);
  }
  // if last coeff is reached, it has to be significant
  estBitsCabac[type].significantBits[pos2ctx_map[type][k1]][0]=0;
  estBitsCabac[type].significantBits[pos2ctx_map[type][k1]][1]=0;
  estBitsCabac[type].lastBits[pos2ctx_last[type][k1]][0]=0;
  estBitsCabac[type].lastBits[pos2ctx_last[type][k1]][1]=0;
}
#endif

/*!
****************************************************************************
* \brief
*    Write Significance MAP
****************************************************************************
*/
#ifdef ADAPTIVE_FD_SD_CODING
void write_significance_map_SD (Macroblock* currMB, EncodingEnvironmentPtr eep_dp, int type, int coeff[], int coeff_ctr)
{
  int   k;
  unsigned short sig, last;
  int   k0      = 0;
  int   k1      = maxpos[type]-1;

  BiContextTypePtr  map_ctx       = img->currentSlice->tex_ctx->map_contexts_SD;
  BiContextTypePtr  last_ctx      = img->currentSlice->tex_ctx->last_contexts_SD;

  if (!c1isdc[type])
  {
    k0++; k1++; coeff--;
  }

  for (k=k0; k<k1; k++) // if last coeff is reached, it has to be significant
  {
    sig   = (coeff[k] != 0);
    biari_encode_symbol  (eep_dp, sig,  map_ctx+pos2ctx_map_int [type][k]);

    if (sig)
    {
      last = (--coeff_ctr == 0);

      biari_encode_symbol(eep_dp, last, last_ctx+pos2ctx_last[type][k]);

      if (last) return;
    }
  }
}




void write_significance_map_SD8 (Macroblock* currMB, EncodingEnvironmentPtr eep_dp, int type, int coeff[], int coeff_ctr)
{
  int   k;
  unsigned short sig, last;
  int   k0      = 0;
  int   k1      = maxpos[type]-1;

  BiContextTypePtr  map_ctx       = img->currentSlice->tex_ctx->map8_contexts_SD;
  BiContextTypePtr  last_ctx      = img->currentSlice->tex_ctx->last8_contexts_SD;

  if (!c1isdc[type])
  {
    k0++; k1++; coeff--;
  }

  for (k=k0; k<k1; k++) // if last coeff is reached, it has to be significant
  {
    sig   = (coeff[k] != 0);
    biari_encode_symbol    (eep_dp, sig, map_ctx+k);

    if (sig)
    {
      last = (--coeff_ctr == 0);

      biari_encode_symbol(eep_dp, last, last_ctx+k);

      if (last) return;
    }
  }
  return;
}
#endif





/*!
 ****************************************************************************
 * \brief
 *    Write Levels
 ****************************************************************************
 */
void write_significant_coefficients (Macroblock* currMB, EncodingEnvironmentPtr eep_dp, int type, int coeff[])
{
  int   i;
  int   absLevel;
  int   ctx;
  short sign;
  short greater_one;
  int   c1 = 1;
  int   c2 = 0;
  
  for (i=maxpos[type]-1; i>=0; i--)
  {
    if (coeff[i]!=0)
    {
      if (coeff[i]>0) 
      {
        absLevel =  coeff[i];
        sign = 0;
      }
      else            
      {
        absLevel = -coeff[i];
        sign = 1;
      }

      greater_one = (absLevel>1);

      //--- if coefficient is one ---
      ctx = min(c1,4);    
      biari_encode_symbol (eep_dp, greater_one, img->currentSlice->tex_ctx->one_contexts[type2ctx_one[type]] + ctx);

      if (greater_one)
      {
        ctx = min(c2, max_c2[type]);
        unary_exp_golomb_level_encode(eep_dp, absLevel-2, img->currentSlice->tex_ctx->abs_contexts[type2ctx_abs[type]] + ctx);
        c1 = 0;
        c2++;
      }
      else if (c1)
      {
        c1++;
      }
      biari_encode_symbol_eq_prob (eep_dp, sign);
    }
  }
}

#ifdef RDO_Q
/*!
 ****************************************************************************
 * \brief
 *    estimate bit cost of significant coefficient
 ****************************************************************************
 */
void est_significant_coefficients (Macroblock* currMB, EncodingEnvironmentPtr eep_dp,  int type)
{
  int   ctx;
  short greater_one;
  int maxCtx=min(4, max_c2[type]);

  for (ctx=0; ctx<=4; ctx++)
  {    
    greater_one=0;
    estBitsCabac[type].greaterOneBits[0][ctx][greater_one]=
      biari_no_bits (eep_dp, greater_one, img->currentSlice->tex_ctx->one_contexts[type2ctx_one[type]] + ctx);

    greater_one=1;
    estBitsCabac[type].greaterOneBits[0][ctx][greater_one]=
      biari_no_bits (eep_dp, greater_one, img->currentSlice->tex_ctx->one_contexts[type2ctx_one[type]] + ctx);
  }

  for (ctx=0; ctx<=maxCtx; ctx++)
  {
    estBitsCabac[type].greaterOneBits[1][ctx][0]=
      biari_no_bits(eep_dp, 0, img->currentSlice->tex_ctx->abs_contexts[type2ctx_abs[type]] + ctx);
  
    estBitsCabac[type].greaterOneState[ctx]=biari_state(eep_dp, 0, img->currentSlice->tex_ctx->abs_contexts[type2ctx_abs[type]] + ctx);

    estBitsCabac[type].greaterOneBits[1][ctx][1]=
      biari_no_bits(eep_dp, 1, img->currentSlice->tex_ctx->abs_contexts[type2ctx_abs[type]] + ctx);
  }
}
#endif

#ifdef ADAPTIVE_FD_SD_CODING
/*!
 ****************************************************************************
 * \brief
 *    Write Levels
 ****************************************************************************
 */
void write_significant_coefficients_SD (Macroblock* currMB, EncodingEnvironmentPtr eep_dp, int type, int coeff[])
{
  int   i;
  int   absLevel;
  int   ctx;
  short sign;
  short greater_one;
  int   c1 = 1;
  int   c2 = 0;


  for (i=maxpos[type]-1; i>=0; i--)
  {
    if (coeff[i]!=0)
    {

      if (coeff[i]>0) 
      {
        absLevel =  coeff[i];
        sign = 0;
      }
      else            
      {
        absLevel = -coeff[i];
        sign = 1;
      }

      greater_one = (absLevel>1);

      //--- if coefficient is one ---
      ctx = min(c1,4);

      biari_encode_symbol (eep_dp, greater_one, img->currentSlice->tex_ctx->one_contexts_SD + ctx);

      if (greater_one)
      {
        ctx = min(c2, max_c2[type]);

        unary_exp_golomb_level_encode(eep_dp, absLevel-2, img->currentSlice->tex_ctx->abs_contexts_SD + ctx);

        c1 = 0;
        c2++;
      }
      else if (c1)
      {
        c1++;
      }

      biari_encode_symbol_eq_prob (eep_dp, sign);
    }
  }
}

void write_significant_coefficients_SD8 (Macroblock* currMB, EncodingEnvironmentPtr eep_dp, int type, int coeff[])
{
  int   i;
  int   absLevel;
  int   ctx;
  short sign;
  short greater_one;
  int   c1 = 1;
  int   c2 = 0;

  for (i=maxpos[type]-1; i>=0; i--)
  {
    if (coeff[i]!=0)
    {
      if (coeff[i]>0) 
      {
        absLevel =  coeff[i];
        sign = 0;
      }
      else            
      {
        absLevel = -coeff[i];
        sign = 1;
      }

      greater_one = (absLevel>1);

      //--- if coefficient is one ---
      ctx = min(c1,4);

      biari_encode_symbol (eep_dp, greater_one, img->currentSlice->tex_ctx->one8_contexts_SD+ ctx);

      if (greater_one)
      {
        ctx = min(c2, max_c2[type]);

        unary_exp_golomb_level_encode(eep_dp, absLevel-2, img->currentSlice->tex_ctx->abs8_contexts_SD + ctx);

        c1 = 0;
        c2++;
      }
      else if (c1)
      {
        c1++;
      }
      biari_encode_symbol_eq_prob (eep_dp, sign);
    }
  }
}
#endif




#ifdef ADAPTIVE_FD_SD_CODING
/*!
 ****************************************************************************
 * \brief
 *    Write Block-Transform Coefficients if coded in frequency domain or
 *    samples if coded in spatial domain
 ****************************************************************************
 */
void writeRunLevel_CABAC (SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  int         writing_b4=se->writing_b4;
  int         writing_b8=se->writing_b8;
#ifdef USE_INTRA_MDDT
  static int  coeff    [256];
#else
  static int  coeff    [64];
#endif
  static int  coeff_ctr = 0;
  static int  pos       = 0;
  int         i;
  int         DF_or_SD=0;
  int         b8_y=((writing_b8/ 2)<<3)+((writing_b4/ 2)<<2);
  int         b8_x=((writing_b8% 2)<<3)+((writing_b4% 2)<<2);
  Macroblock* currMB    = &img->mb_data[img->current_mb_nr];
  int         need_SD_FD_bit     =(currMB->b8mode[writing_b8] != IBLOCK && currMB->luma_transform_size_8x8_flag==0 && se->context==LUMA_4x4 && !IS_INTRA (currMB));
  int         need_SD_FD_bit_t8x8=(currMB->b8mode[writing_b8] != IBLOCK && currMB->luma_transform_size_8x8_flag==1 && se->context==LUMA_8x8 && !IS_INTRA (currMB));
  int SD;

  if (img->APEC_in_FD_and_SD==0)
  {
    need_SD_FD_bit=0;
    need_SD_FD_bit_t8x8=0;
  }
  if (img->APEC_in_FD_and_SD==1 && currMB->SD_Coding_on_off==0)
  {
    need_SD_FD_bit=0;
    need_SD_FD_bit_t8x8=0;
  }

  if (cabac_encoding && (currMB->b8mode[writing_b8] == IBLOCK || currMB->luma_transform_size_8x8_flag==0 || se->context==LUMA_4x4 || IS_INTRA (currMB)))
    set_bit ( &(currMB->SD_or_FD_t8x8), writing_b8, 0);

  SD= (need_SD_FD_bit      && (currMB->SD_or_FD[writing_b8/ 2][writing_b8% 2] & (1<<writing_b4))!=0) ||
          (need_SD_FD_bit_t8x8 && ((currMB->SD_or_FD_t8x8) & (1<<writing_b8))!=0);

  //--- accumulate run-level information ---
  if (se->value1 != 0)
  {
    if (SD==0)
    {
      pos += se->value2;
      coeff[pos++] = se->value1;
      coeff_ctr++;
      //return;
    }
  }
  else
  {
    //===== encode CBP-BIT =====
    if (coeff_ctr>0 || need_SD_FD_bit || need_SD_FD_bit_t8x8)
    {
      if (need_SD_FD_bit)
      {
        set_bit ( &(currMB->SD_or_FD_t8x8), writing_b8, 0);

        DF_or_SD=0;

        if      (currMB->SD_or_FD[writing_b8/ 2][writing_b8% 2] & (1<<writing_b4) )
        {
          DF_or_SD=1;
        }

        if (DF_or_SD)
        {
          coeff_ctr   = 0;
          for (i=0; i<16; i++)
          {
            coeff    [i] = currMB->quantizer_indices[b8_y+(i/4)][b8_x+(i%4)];
            if( coeff[i] != 0 ) coeff_ctr++;
          }
          if (coeff_ctr>0)
          {
            write_and_store_CBP_block_bit                (currMB, eep_dp, se->context, 1);
            write_spatial_domain_coding_flag_for_4x4block(currMB, eep_dp, se->context, DF_or_SD);
            //===== encode significance map =====
            write_significance_map_SD         (currMB, eep_dp, se->context, coeff,coeff_ctr);
            //===== encode significant coefficients =====
            write_significant_coefficients_SD (currMB, eep_dp, se->context, coeff);
          }
          else
          {
            write_and_store_CBP_block_bit  (currMB, eep_dp, se->context, 0);
          }
        }
        else
        {//!DF_or_SD
          if (coeff_ctr>0)
          {
            write_and_store_CBP_block_bit                     (currMB, eep_dp, se->context, 1);
            write_spatial_domain_coding_flag_for_4x4block     (currMB, eep_dp, se->context, DF_or_SD);
            //===== encode significance map =====
            write_significance_map         (currMB, eep_dp, se->context, coeff, coeff_ctr);
            //===== encode significant coefficients =====
            write_significant_coefficients (currMB, eep_dp, se->context, coeff);
          }
          else
          {
            write_and_store_CBP_block_bit          (currMB, eep_dp, se->context, 0);
          }
        }
      }//need_SD_FD_bit
      else if (need_SD_FD_bit_t8x8)
      {
        DF_or_SD=((currMB->SD_or_FD_t8x8) & (1<<writing_b8))!=0?1:0;

        if (DF_or_SD)
        {
          coeff_ctr   = 0;
          for (i=0; i<64; i++)
          {
            coeff   [i] = currMB->quantizer_indices[((writing_b8/ 2)<<3)+(i/8)][((writing_b8% 2)<<3)+(i%8)];
            if( coeff[i] != 0 )  coeff_ctr++;
          }

          if (coeff_ctr>0)
          {
            write_and_store_CBP_block_bit                 (currMB, eep_dp, se->context, 1);
            write_spatial_domain_coding_flag_for_8x8block (currMB, eep_dp, se->context, DF_or_SD, writing_b8);

            //===== encode significance map =====
            write_significance_map_SD8         (currMB, eep_dp, se->context, coeff, coeff_ctr);
            //===== encode significant coefficients =====
            write_significant_coefficients_SD8 (currMB, eep_dp, se->context, coeff);
          }
          else
          {
            write_and_store_CBP_block_bit  (currMB, eep_dp, se->context, 0);
          }
        }
        else
        {//!DF_or_SD
          if (coeff_ctr>0)
          {
            write_and_store_CBP_block_bit                    (currMB, eep_dp, se->context, 1);
            write_spatial_domain_coding_flag_for_8x8block    (currMB, eep_dp, se->context, DF_or_SD, writing_b8);
            //===== encode significance map =====
            write_significance_map         (currMB, eep_dp, se->context, coeff, coeff_ctr);
            //===== encode significant coefficients =====
            write_significant_coefficients (currMB, eep_dp, se->context, coeff);
          }
          else
          {
            write_and_store_CBP_block_bit          (currMB, eep_dp, se->context, 0);
          }
        }
      }
      else
      {
        write_and_store_CBP_block_bit  (currMB, eep_dp, se->context, 1);
        //===== encode significance map =====
        write_significance_map         (currMB, eep_dp, se->context, coeff, coeff_ctr);
        //===== encode significant coefficients =====
        write_significant_coefficients (currMB, eep_dp, se->context, coeff);
      }
    }
    else
      write_and_store_CBP_block_bit  (currMB, eep_dp, se->context, 0);

    //--- reset counters ---
    pos = coeff_ctr = 0;
#ifdef USE_INTRA_MDDT
    memset(coeff, 0 , 256 * sizeof(int));
#else
    memset(coeff, 0 , 64 * sizeof(int));
#endif
  }
}
#else
/*!
 ****************************************************************************
 * \brief
 *    Write Block-Transform Coefficients
 ****************************************************************************
 */
void writeRunLevel_CABAC (SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
#ifdef USE_INTRA_MDDT
  static int  coeff[256];
#else
  static int  coeff[64];
#endif
  static int  coeff_ctr = 0;
  static int  pos       = 0;
    
  //--- accumulate run-level information ---
  if (se->value1 != 0)
  {
    pos += se->value2;
    coeff[pos++] = se->value1; 
    coeff_ctr++;
    //return;
  }
  else
  {
    Macroblock* currMB    = &img->mb_data[img->current_mb_nr];
    //===== encode CBP-BIT =====
    if (coeff_ctr>0)
    {      
      write_and_store_CBP_block_bit  (currMB, eep_dp, se->context, 1);      
      //===== encode significance map =====
      write_significance_map         (currMB, eep_dp, se->context, coeff, coeff_ctr);      
      //===== encode significant coefficients =====
      write_significant_coefficients (currMB, eep_dp, se->context, coeff);
    }
    else
      write_and_store_CBP_block_bit  (currMB, eep_dp, se->context, 0);
    
    //--- reset counters ---
    pos = coeff_ctr = 0;  
#ifdef USE_INTRA_MDDT
    memset(coeff, 0 , 256 * sizeof(int));
#else
    memset(coeff, 0 , 64 * sizeof(int));
#endif
  }
}
#endif

#ifdef RDO_Q
/*!
 ****************************************************************************
 * \brief
 *   estimate bit cost for CBP, significant map and significant coefficients
 ****************************************************************************
 */
void estRunLevel_CABAC (int context) // marta - writes CABAC run/level 
{
  Macroblock* currMB    = &img->mb_data[img->current_mb_nr];
  Slice*          currSlice = img->currentSlice;
  DataPartition*  dataPart = &(currSlice->partArr[0]); // assumed that no DP is used (table assignSE2partition_NoDP)
  EncodingEnvironmentPtr eep_dp = &(dataPart->ee_cabac); 

  est_CBP_block_bit  (currMB, eep_dp, context);      
  //===== encode significance map =====
  est_significance_map         (currMB, eep_dp, context);      
  //===== encode significant coefficients =====
  est_significant_coefficients (currMB, eep_dp, context);
}

#define SIGN_BITS    1
/*!
 ****************************************************************************
 * \brief
 *    Rate distortion optimized trellis quantization
 ****************************************************************************
 */
void est_writeRunLevel_CABAC(levelDataStruct levelData[], int levelTabMin[], int type, double lambda, int kInit, int kStop, 
                              int noCoeff, int estCBP)
{
  int   k, i;
  int estBits;
  double lagr, lagrMin=0, lagrTabMin, lagrTab;
  int   c1 = 1, c2 = 0, c1Tab[3], c2Tab[3];
  int   iBest, levelTab[256];
  int   ctx, greater_one, last, maxK;
  double   lagrAcc, lagrLastMin=0, lagrLast;
  int      kBest=0, kStart, first;

  maxK=maxpos[type];
  for (k=0; k<maxK; k++)
  {
    levelTabMin[k]=0;
  }

  if (noCoeff>0)
  {
    if (noCoeff>1)
    {
      kStart=kInit; kBest=0; first=1; 

      lagrAcc=0; 
      for (k=kStart; k<=kStop; k++)
      {
        lagrAcc+=levelData[k].errLevel[0];
      }

      if (levelData[kStart].noLevels>2)
      { 
        lagrAcc-=levelData[kStart].errLevel[0];
        lagrLastMin=lambda*(estBitsCabac[type].lastBits[pos2ctx_last[type][kStart]][1]-estBitsCabac[type].lastBits[pos2ctx_last[type][kStart]][0])+lagrAcc;

        kBest=kStart;
        kStart=kStart+1;
        first=0;
      }

      for (k=kStart; k<=kStop; k++)
      {

        lagrMin=levelData[k].errLevel[0]+lambda*estBitsCabac[type].significantBits[pos2ctx_map[type][k]][0];

        lagrAcc-=levelData[k].errLevel[0];
        if (levelData[k].noLevels>1)
        { 

          estBits=SIGN_BITS+estBitsCabac[type].significantBits[pos2ctx_map[type][k]][1]+
            estBitsCabac[type].greaterOneBits[0][4][0];

          lagrLast=levelData[k].errLevel[1]+lambda*(estBits+estBitsCabac[type].lastBits[pos2ctx_last[type][k]][1])+lagrAcc;
          lagr=levelData[k].errLevel[1]+lambda*(estBits+estBitsCabac[type].lastBits[pos2ctx_last[type][k]][0]);

          lagrMin=(lagr<lagrMin)?lagr:lagrMin;

          if (lagrLast<lagrLastMin || first==1)
          {
            kBest=k;
            first=0;
            lagrLastMin=lagrLast;
          }

        }
        lagrAcc+=lagrMin;
      }

      kStart=kBest;
    }
    else
    {
      kStart=kStop;
    }

    lagrTabMin=0;
    for (k=0; k<=kStart; k++)
    {
      lagrTabMin+=levelData[k].errLevel[0];
    }
    // Initial Lagrangian calculation
    lagrTab=0;

    //////////////////////////

    lagrTabMin+=(lambda*estCBP);
    iBest=0; first=1;
    for (k=kStart; k>=0; k--)
    {

      last=(k==kStart);
      if (!last)
      {
        lagrMin=levelData[k].errLevel[0]+lambda*estBitsCabac[type].significantBits[pos2ctx_map[type][k]][0];
        iBest=0;
        first=0;
      }

      for (i=1; i<levelData[k].noLevels; i++)
      {

        estBits=SIGN_BITS+estBitsCabac[type].significantBits[pos2ctx_map[type][k]][1];
        estBits+=estBitsCabac[type].lastBits[pos2ctx_last[type][k]][last];

        // greater than 1
        greater_one = (levelData[k].level[i]>1);

        c1Tab[i]=c1;   c2Tab[i]=c2;

        ctx = min(c1Tab[i],4);  
        estBits+=estBitsCabac[type].greaterOneBits[0][ctx][greater_one];

        // magnitude if greater than 1
        if (greater_one)
        {
          ctx = min(c2Tab[i], max_c2[type]);
          if ((levelData[k].level[i]-2)<MAX_PREC_COEFF)
          {
            estBits+=precalcUnaryLevelTab[estBitsCabac[type].greaterOneState[ctx]][levelData[k].level[i]-2];
          }
          else
          {
            estBits+=est_unary_exp_golomb_level_encode((unsigned int)levelData[k].level[i]-2, ctx, type);
          }

          c1Tab[i] = 0;
          c2Tab[i]++;
        }
        else if (c1Tab[i])
        {
          c1Tab[i]++;
        }

        lagr=levelData[k].errLevel[i]+lambda*estBits;
        if (lagr<lagrMin || first==1)
        {
          iBest=i;
          lagrMin=lagr;
          first=0;
        }
      }

      if (iBest>0)
      {
        c1=c1Tab[iBest]; c2=c2Tab[iBest];
      }

      levelTab[k] = (int)levelData[k].level[iBest];
      lagrTab+=lagrMin;
    }
    ///////////////////////////////////

    if (lagrTab<lagrTabMin)
    {
      for (k=0; k<=kStart; k++)
      {
        levelTabMin[k]=levelTab[k];
      }
    }
  }
}
#endif

#ifdef USE_INTRA_MDDT
/*! 
*************************************************************************************
* \brief
*   write all the runs and levels of an INTRA16x16 macroblock
*
* \para writeRunLevel16x16_CABAC()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
void writeRunLevel16x16_CABAC (SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  static int  coeff[256];
  static int  coeff_ctr = 0;
  static int  pos       = 0;
    
  //--- accumulate run-level information ---
  if (se->value1 != 0)
  {
    pos += se->value2;
    coeff[pos++] = se->value1; 
    coeff_ctr++;
    //return;
  }
  else
  {
    Macroblock* currMB    = &img->mb_data[img->current_mb_nr];
    //===== encode CBP-BIT =====
    if (coeff_ctr>0)
    {      
      // for 16x16, just set all block CBP to 1
      set_CBP_block_bit_16x16  (currMB);      
      //===== encode significance map =====
      write_significance_map         (currMB, eep_dp, se->context, coeff, coeff_ctr);      
      //===== encode significant coefficients =====
      write_significant_coefficients (currMB, eep_dp, se->context, coeff);
    }
//    else
//      write_and_store_CBP_block_bit  (currMB, eep_dp, se->context, 0);
    
    //--- reset counters ---
    pos = coeff_ctr = 0;  
    memset(coeff, 0 , 256 * sizeof(int));
  }
}
#endif

#ifdef MB32X32
void set_CBP_block_bit_16x8 (Macroblock* currMB)
{
  int bit;

  for(bit = img->subblock_y*4+1; bit < img->subblock_y*4+9; bit++)
  {
    currMB->cbp_bits   |= ((int64)1<< bit   );
  }
}
void set_CBP_block_bit_8x16 (Macroblock* currMB)
{
  int bit;

  for(bit = img->subblock_x+1; bit < img->subblock_x+3; bit++)
  {
    currMB->cbp_bits   |= ((int64)1<< bit   );
  }
  for(bit = img->subblock_x+5; bit < img->subblock_x+7; bit++)
  {
    currMB->cbp_bits   |= ((int64)1<< bit   );
  }
  for(bit = img->subblock_x+9; bit < img->subblock_x+11; bit++)
  {
    currMB->cbp_bits   |= ((int64)1<< bit   );
  }
  for(bit =img->subblock_x+13; bit < img->subblock_x+15; bit++)
  {
    currMB->cbp_bits   |= ((int64)1<< bit   );
  }
}
/*! 
*************************************************************************************
* \brief
*   write all the runs and levels of an INTRA16x16 macroblock
*
* \para writeRunLevel16x16P_CABAC()
*    <paragraph>
*
* \author
*    - Yan Ye                      <yye@qualcomm.com>
*************************************************************************************
*/
void writeRunLevel16x16P_CABAC (SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  static int  coeff[256];
  static int  coeff_ctr = 0;
  static int  pos       = 0;
    
  //--- accumulate run-level information ---
  if (se->value1 != 0)
  {
    pos += se->value2;
    coeff[pos++] = se->value1; 
    coeff_ctr++;
    //return;
  }
  else
  {
    Macroblock* currMB    = &img->mb_data[img->current_mb_nr];
    //===== encode CBP-BIT =====
    if (coeff_ctr>0)
    {      
      // for 16x16, just set all block CBP to 1
      set_CBP_block_bit_16x16  (currMB);      
      //===== encode significance map =====
      write_significance_map         (currMB, eep_dp, se->context, coeff, coeff_ctr);      
      //===== encode significant coefficients =====
      write_significant_coefficients (currMB, eep_dp, se->context, coeff);
    }
//    else
//      write_and_store_CBP_block_bit  (currMB, eep_dp, se->context, 0);
    
    //--- reset counters ---
    pos = coeff_ctr = 0;  
    memset(coeff, 0 , 256 * sizeof(int));
  }
}

void writeRunLevel16x8_CABAC (SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  static int  coeff[256];
  static int  coeff_ctr = 0;
  static int  pos       = 0;
    
  //--- accumulate run-level information ---
  if (se->value1 != 0)
  {
    pos += se->value2;
    coeff[pos++] = se->value1; 
    coeff_ctr++;
    //return;
  }
  else
  {
    Macroblock* currMB    = &img->mb_data[img->current_mb_nr];
    //===== encode CBP-BIT =====
    if (coeff_ctr>0)
    {      
      // for 16x16, just set all block CBP to 1
      set_CBP_block_bit_16x8  (currMB);      
      //===== encode significance map =====
      write_significance_map         (currMB, eep_dp, se->context, coeff, coeff_ctr);      
      //===== encode significant coefficients =====
      write_significant_coefficients (currMB, eep_dp, se->context, coeff);
    }
//    else
//      write_and_store_CBP_block_bit  (currMB, eep_dp, se->context, 0);
    
    //--- reset counters ---
    pos = coeff_ctr = 0;  
    memset(coeff, 0 , 256 * sizeof(int));
  }
}

void writeRunLevel8x16_CABAC (SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  static int  coeff[256];
  static int  coeff_ctr = 0;
  static int  pos       = 0;
    
  //--- accumulate run-level information ---
  if (se->value1 != 0)
  {
    pos += se->value2;
    coeff[pos++] = se->value1; 
    coeff_ctr++;
    //return;
  }
  else
  {
    Macroblock* currMB    = &img->mb_data[img->current_mb_nr];
    //===== encode CBP-BIT =====
    if (coeff_ctr>0)
    {      
      // for 16x16, just set all block CBP to 1
      set_CBP_block_bit_8x16  (currMB);      
      //===== encode significance map =====
      write_significance_map         (currMB, eep_dp, se->context, coeff, coeff_ctr);      
      //===== encode significant coefficients =====
      write_significant_coefficients (currMB, eep_dp, se->context, coeff);
    }
//    else
//      write_and_store_CBP_block_bit  (currMB, eep_dp, se->context, 0);
    
    //--- reset counters ---
    pos = coeff_ctr = 0;  
    memset(coeff, 0 , 256 * sizeof(int));
  }
}
#endif

/*!
 ************************************************************************
 * \brief
 *    Unary binarization and encoding of a symbol by using
 *    one or two distinct models for the first two and all
 *    remaining bins
*
************************************************************************/
void unary_bin_encode(EncodingEnvironmentPtr eep_dp,
                      unsigned int symbol,
                      BiContextTypePtr ctx,
                      int ctx_offset)
{
  unsigned int l;
  BiContextTypePtr ictx;

  if (symbol==0)
  {
    biari_encode_symbol(eep_dp, 0, ctx );
    return;
  }
  else
  {
    biari_encode_symbol(eep_dp, 1, ctx );
    l = symbol;
    ictx = ctx+ctx_offset;
    while ((--l)>0)
      biari_encode_symbol(eep_dp, 1, ictx);
    biari_encode_symbol(eep_dp, 0, ictx);
  }
  return;
}

/*!
 ************************************************************************
 * \brief
 *    Unary binarization and encoding of a symbol by using
 *    one or two distinct models for the first two and all
 *    remaining bins; no terminating "0" for max_symbol
 *    (finite symbol alphabet)
 ************************************************************************
 */
void unary_bin_max_encode(EncodingEnvironmentPtr eep_dp,
                          unsigned int symbol,
                          BiContextTypePtr ctx,
                          int ctx_offset,
                          unsigned int max_symbol)
{
  unsigned int l;
  BiContextTypePtr ictx;

  if (symbol==0)
  {
    biari_encode_symbol(eep_dp, 0, ctx );
    return;
  }
  else
  {
    biari_encode_symbol(eep_dp, 1, ctx );
    l=symbol;
    ictx=ctx+ctx_offset;
    while ((--l)>0)
      biari_encode_symbol(eep_dp, 1, ictx);
    if (symbol<max_symbol)
      biari_encode_symbol(eep_dp, 0, ictx);
  }
  return;
}



/*!
 ************************************************************************
 * \brief
 *    Exp Golomb binarization and encoding
 ************************************************************************
 */
void exp_golomb_encode_eq_prob( EncodingEnvironmentPtr eep_dp,
                                unsigned int symbol,
                                int k) 
{
  while(1)
  {
    if (symbol >= (unsigned int)(1<<k))   
    {
      biari_encode_symbol_eq_prob(eep_dp, 1);   //first unary part
      symbol = symbol - (1<<k);
      k++;
    }
    else                  
    {
      biari_encode_symbol_eq_prob(eep_dp, 0);   //now terminated zero of unary part
      while (k--)                               //next binary part
        biari_encode_symbol_eq_prob(eep_dp, (signed short)((symbol>>k)&1)); 
      break;
    }
  }

  return;
}

#ifdef RDO_Q
/*!
 ****************************************************************************
 * \brief
 *    estimate exp golomb bit cost 
 ****************************************************************************
 */
int est_exp_golomb_encode_eq_prob(unsigned int symbol)
{
  int k=0, estBits=0;
  
  while(1)
  {
    if (symbol >= (unsigned int)(1<<k))   
    {
      estBits++;
      symbol = symbol - (1<<k);
      k++;
    }
    else                  
    {
      estBits++;  
      while (k--)  
      {
        estBits++;
      }
      break;
    }
  }
  return(estBits << 15);
}
#endif

/*!
 ************************************************************************
 * \brief
 *    Exp-Golomb for Level Encoding
*
************************************************************************/
void unary_exp_golomb_level_encode( EncodingEnvironmentPtr eep_dp,
                                    unsigned int symbol,
                                    BiContextTypePtr ctx)
{
  unsigned int l,k;
  unsigned int exp_start = 13; // 15-2 : 0,1 level decision always sent

  if (symbol==0)
  {
    biari_encode_symbol(eep_dp, 0, ctx );
    return;
  }
  else
  {
    biari_encode_symbol(eep_dp, 1, ctx );
    l=symbol;
    k=1;
    while (((--l)>0) && (++k <= exp_start))
      biari_encode_symbol(eep_dp, 1, ctx);
    if (symbol < exp_start) biari_encode_symbol(eep_dp, 0, ctx);
    else exp_golomb_encode_eq_prob(eep_dp,symbol-exp_start,0);
  }
  return;
}

#ifdef RDO_Q
/*!
 ****************************************************************************
 * \brief
 *    estimate unary exp golomb bit cost
 ****************************************************************************
 */
int est_unary_exp_golomb_level_encode(unsigned int symbol, int ctx, int type)
{
  unsigned int l,k;
  unsigned int exp_start = 13; // 15-2 : 0,1 level decision always sent
  int estBits;

  if (symbol==0)
  {
    estBits=estBitsCabac[type].greaterOneBits[1][ctx][0];
    return (estBits);
  }
  else
  {
    estBits=estBitsCabac[type].greaterOneBits[1][ctx][1];
    l=symbol;
    k=1;
    while (((--l)>0) && (++k <= exp_start))
    {
      estBits+=estBitsCabac[type].greaterOneBits[1][ctx][1];
    }
    if (symbol < exp_start)
    {
      estBits+=estBitsCabac[type].greaterOneBits[1][ctx][0];
    }
    else 
    {
      estBits+=est_exp_golomb_encode_eq_prob(symbol-exp_start);
    }
  }
  return(estBits);
}

int est_unary_exp_golomb_level_bits(unsigned int symbol, int bits0, int bits1)
{
  unsigned int l,k;
  unsigned int exp_start = 13; // 15-2 : 0,1 level decision always sent
  int estBits;

  if (symbol==0)
  {
    return (bits0);
  }
  else
  {
    estBits=bits1;
    l=symbol;
    k=1;
    while (((--l)>0) && (++k <= exp_start))
    {
      estBits+=bits1;
    }
    if (symbol < exp_start)
    {
      estBits+=bits0;
    }
    else 
    {
      estBits+=est_exp_golomb_encode_eq_prob(symbol-exp_start);
    }
  }
  return(estBits);
}

void precalculate_unary_exp_golomb_level()
{
  int state, ctx_state0, ctx_state1, estBits0, estBits1, symbol;

  for (state=0; state<=63; state++)
  {
    // symbol 0 is MPS
    ctx_state0=64+state;
    estBits0=entropyBits[127-ctx_state0];
    ctx_state1=63-state;
    estBits1=entropyBits[127-ctx_state1];

    for (symbol=0; symbol<MAX_PREC_COEFF; symbol++)
    {
      precalcUnaryLevelTab[ctx_state0][symbol]=est_unary_exp_golomb_level_bits(symbol, estBits0, estBits1);

      // symbol 0 is LPS
      precalcUnaryLevelTab[ctx_state1][symbol]=est_unary_exp_golomb_level_bits(symbol, estBits1, estBits0);
    }
  }
}

#endif


/*!
 ************************************************************************
 * \brief
 *    Exp-Golomb for MV Encoding
*
************************************************************************/
void unary_exp_golomb_mv_encode(EncodingEnvironmentPtr eep_dp,
                                unsigned int symbol,
                                BiContextTypePtr ctx,
                                unsigned int max_bin)
{
  unsigned int l,k;
  unsigned int bin=1;
  BiContextTypePtr ictx=ctx;
  unsigned int exp_start = 8; // 9-1 : 0 mvd decision always sent

  if (symbol==0)
  {
    biari_encode_symbol(eep_dp, 0, ictx );
    return;
  }
  else
  {
    biari_encode_symbol(eep_dp, 1, ictx );
    l=symbol;
    k=1;
    ictx++;
    while (((--l)>0) && (++k <= exp_start))
    {
      biari_encode_symbol(eep_dp, 1, ictx  );
      if ((++bin)==2) ictx++;
      if (bin==max_bin) ictx++;
    }
    if (symbol < exp_start) biari_encode_symbol(eep_dp, 0, ictx);
    else exp_golomb_encode_eq_prob(eep_dp,symbol-exp_start,3);
  }
  return;
}


#ifdef ADAPTIVE_FD_SD_CODING
void set_bit(int *dest, int bit, int value)
{
  (*dest)=(value==0?(*dest)&(~(1<<bit)) : (*dest)|(1<<bit));
}
#endif

#ifdef ADAPTIVE_QUANTIZATION
/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode the iaqms_idx
 ****************************************************************************
 */
void writeIAQMS_idx_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  MotionInfoContexts *ctx = img->currentSlice->mot_ctx;

  // remaining_mode_selector
  biari_encode_symbol(eep_dp,(signed short)( se->value1 & 0x1    ), ctx->modulated_quantization_contexts+1);
  if( !(img->type==B_SLICE || (active_sps->profile_idc < FREXT_HP && img->type==P_SLICE)) )
    biari_encode_symbol(eep_dp,(signed short)((se->value1 & 0x2)>>1), ctx->modulated_quantization_contexts+1);
}
#endif

#ifdef ADAPTIVE_LOOP_FILTER
/*!
 ***************************************************************************
 * \brief
 *    This function is used to arithmetically encode the alf_blk_flag
 *    info of a given Block.
 ***************************************************************************
 */
void writeAlfBlockFlagCABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  int a, b;
  int act_ctx = 0;
  int act_sym;
  
  MotionInfoContexts *ctx        = (img->currentSlice)->mot_ctx;
  AlfBlock           *currAlfBlk = &img->alf_blk_data[img->current_alfb_nr];
  
  b = (currAlfBlk->alfb_available_up   == NULL) ? 0 : currAlfBlk->alfb_available_up->alf_blk_flag;
  a = (currAlfBlk->alfb_available_left == NULL) ? 0 : currAlfBlk->alfb_available_left->alf_blk_flag;

  act_ctx     = a + b;
  act_sym     = se->value1;
  se->context = act_ctx; // store context

  biari_encode_symbol(eep_dp, (short)act_sym, ctx->alf_blk_flag_contexts + act_ctx );  
}
/*!
 ************************************************************************
 * \brief
 *    Check for available neighbouring blocks
 *    and set pointers in current block
 ************************************************************************
 */
void CheckAvailabilityOfNeighborsAlfCABAC()
{
  AlfBlock *currAlfBlk = &img->alf_blk_data[img->current_alfb_nr];
  PixelPos up, left;
  
  GetNeighbourAlfBlock(img->current_alfb_nr, -1,  0, &left);
  GetNeighbourAlfBlock(img->current_alfb_nr,  0, -1, &up);
  
  if (up.available)
    currAlfBlk->alfb_available_up = &img->alf_blk_data[up.mb_addr];
  else
    currAlfBlk->alfb_available_up = NULL;
  
  if (left.available)
    currAlfBlk->alfb_available_left = &img->alf_blk_data[left.mb_addr];
  else
    currAlfBlk->alfb_available_left = NULL;
}
#endif
