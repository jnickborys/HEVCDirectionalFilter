
/*!
 *************************************************************************************
 * \file biariencode.c
 *
 * \brief
 *    Routines for binary arithmetic encoding
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Detlev Marpe                    <marpe@hhi.de>
 *    - Gabi Blaettermann               <blaetter@hhi.de>
 *************************************************************************************
 */

#include <stdlib.h>
#include <stdio.h>

#include "global.h"
#include "biariencode.h"
#ifdef RDO_Q
#include <math.h>

int entropyBits[128]=
{
     895,    943,    994,   1048,   1105,   1165,   1228,   1294, 
    1364,   1439,   1517,   1599,   1686,   1778,   1875,   1978, 
    2086,   2200,   2321,   2448,   2583,   2725,   2876,   3034, 
    3202,   3380,   3568,   3767,   3977,   4199,   4435,   4684, 
    4948,   5228,   5525,   5840,   6173,   6527,   6903,   7303, 
    7727,   8178,   8658,   9169,   9714,  10294,  10914,  11575, 
   12282,  13038,  13849,  14717,  15650,  16653,  17734,  18899, 
   20159,  21523,  23005,  24617,  26378,  28306,  30426,  32768, 
   32768,  35232,  37696,  40159,  42623,  45087,  47551,  50015, 
   52479,  54942,  57406,  59870,  62334,  64798,  67262,  69725, 
   72189,  74653,  77117,  79581,  82044,  84508,  86972,  89436, 
   91900,  94363,  96827,  99291, 101755, 104219, 106683, 109146, 
  111610, 114074, 116538, 119002, 121465, 123929, 126393, 128857, 
  131321, 133785, 136248, 138712, 141176, 143640, 146104, 148568, 
  151031, 153495, 155959, 158423, 160887, 163351, 165814, 168278, 
  170742, 173207, 175669, 178134, 180598, 183061, 185525, 187989
};
#endif

int binCount = 0;

/*!
 ************************************************************************
 * Macro for writing bytes of code
 ***********************************************************************
 */

#define put_byte() { \
                     Ecodestrm[(*Ecodestrm_len)++] = Ebuffer; \
                     Ebits_to_go = 8; \
                     while (eep->C > 7) { \
                       eep->C-=8; \
                       eep->E++; \
                     } \
                    } 

#define put_one_bit(b) { \
                         Ebuffer <<= 1; Ebuffer |= (b); \
                         if (--Ebits_to_go == 0) \
                           put_byte(); \
                       }

#define put_one_bit_plus_outstanding(b) { \
                                          put_one_bit(b); \
                                          while (Ebits_to_follow > 0) \
                                          { \
                                            Ebits_to_follow--; \
                                            put_one_bit(!(b)); \
                                          } \
                                         }

int pic_bin_count;

void reset_pic_bin_count(void)
{
  pic_bin_count = 0;
}

int get_pic_bin_count(void)
{
  return pic_bin_count;
}



/*!
 ************************************************************************
 * \brief
 *    Allocates memory for the EncodingEnvironment struct
 ************************************************************************
 */
EncodingEnvironmentPtr arienco_create_encoding_environment(void)
{
  EncodingEnvironmentPtr eep;

  if ( (eep = (EncodingEnvironmentPtr) calloc(1,sizeof(EncodingEnvironment))) == NULL)
    no_mem_exit("arienco_create_encoding_environment: eep");

  return eep;
}



/*!
 ************************************************************************
 * \brief
 *    Frees memory of the EncodingEnvironment struct
 ************************************************************************
 */
void arienco_delete_encoding_environment(EncodingEnvironmentPtr eep)
{
  if (eep == NULL)
  {
    snprintf(errortext, ET_SIZE, "Error freeing eep (NULL pointer)");
    error (errortext, 200);
  }
  else
    free(eep);
}



/*!
 ************************************************************************
 * \brief
 *    Initializes the EncodingEnvironment for the arithmetic coder
 ************************************************************************
 */
void arienco_start_encoding(EncodingEnvironmentPtr eep,
                            unsigned char *code_buffer,
                            int *code_len )
{
  Elow = 0;
  Ebits_to_follow = 0;
  Ebuffer = 0;
  Ebits_to_go = 9; // to swallow first redundant bit

  Ecodestrm = code_buffer;
  Ecodestrm_len = code_len;

  Erange = HALF-2;

  eep->C = 0;
  eep->E = 0;

}

/*!
 ************************************************************************
 * \brief
 *    Returns the number of currently written bits
 ************************************************************************
 */
int arienco_bits_written(EncodingEnvironmentPtr eep)
{
   return (8 * (*Ecodestrm_len) + Ebits_to_follow + 8  - Ebits_to_go);
}


/*!
 ************************************************************************
 * \brief
 *    Terminates the arithmetic codeword, writes stop bit and stuffing bytes (if any)
 ************************************************************************
 */
void arienco_done_encoding(EncodingEnvironmentPtr eep)
{
  put_one_bit_plus_outstanding((unsigned char) ((Elow >> (B_BITS-1)) & 1));
  put_one_bit((unsigned char) (Elow >> (B_BITS-2))&1);
  put_one_bit((unsigned char) 1);

  stats->bit_use_stuffingBits[img->type]+=(8-Ebits_to_go);

  while (Ebits_to_go != 8)
    put_one_bit(0);

  pic_bin_count += eep->E*8 + eep->C; // no of processed bins
}

extern int cabac_encoding;

/*!
 ************************************************************************
 * \brief
 *    Actually arithmetic encoding of one binary symbol by using
 *    the probability estimate of its associated context model
 ************************************************************************
 */
void biari_encode_symbol(EncodingEnvironmentPtr eep, signed short symbol, BiContextTypePtr bi_ct )
{
  register unsigned int range = Erange;
  register unsigned int low = Elow;
  unsigned int rLPS = rLPS_table_64x4[bi_ct->state][(range>>6) & 3];
  
#if (2==TRACE)
  if (cabac_encoding)
    fprintf(p_trace, "%d  0x%04x  %d  %d\n", binCount++, Erange , bi_ct->state, bi_ct->MPS );
#endif
  
  range -= rLPS;  
  bi_ct->count += cabac_encoding;

  /* covers all cases where code does not bother to shift down symbol to be 
   * either 0 or 1, e.g. in some cases for cbp, mb_Type etc the code simply 
   * masks off the bit position and passes in the resulting value */
  symbol = (short) (symbol != 0);

  if (symbol != bi_ct->MPS) 
  {
    low += range;
    range = rLPS;
    
    if (!bi_ct->state)
      bi_ct->MPS = (unsigned char) (bi_ct->MPS ^ 0x01);               // switch LPS if necessary
    bi_ct->state = AC_next_state_LPS_64[bi_ct->state]; // next state
  } 
  else 
    bi_ct->state = AC_next_state_MPS_64[bi_ct->state]; // next state

  /* renormalisation */    
  while (range < QUARTER)
  {
    if (low >= HALF)
    {
      put_one_bit_plus_outstanding(1);
      low -= HALF;
    }
    else if (low < QUARTER)
    {
      put_one_bit_plus_outstanding(0);
    }
    else
    {
      Ebits_to_follow++;
      low -= QUARTER;
    }
    low <<= 1;
    range <<= 1;
  }
  Erange = range;
  Elow = low;
  eep->C++;
}

#ifdef RDO_Q
int biari_no_bits(EncodingEnvironmentPtr eep, signed short symbol, BiContextTypePtr bi_ct )
{
  int ctx_state, estBits;

  symbol = (short) (symbol != 0);

  ctx_state=(symbol==bi_ct->MPS)?64+bi_ct->state:63-bi_ct->state;
  estBits=entropyBits[127-ctx_state];

  return(estBits);
}

int biari_state(EncodingEnvironmentPtr eep, signed short symbol, BiContextTypePtr bi_ct )
{ 
  int ctx_state;

  symbol = (short) (symbol != 0);
  ctx_state=(symbol==bi_ct->MPS)?64+bi_ct->state:63-bi_ct->state;

  return(ctx_state);
}

#endif

/*!
 ************************************************************************
 * \brief
 *    Arithmetic encoding of one binary symbol assuming 
 *    a fixed prob. distribution with p(symbol) = 0.5
 ************************************************************************
 */
void biari_encode_symbol_eq_prob(EncodingEnvironmentPtr eep, signed short symbol)
{
  register unsigned int low = (Elow<<1);
  
#if (2==TRACE)
  extern int cabac_encoding;
  if (cabac_encoding)
    fprintf(p_trace, "%d  0x%04x\n", binCount++, Erange );
#endif
  
  if (symbol != 0)
    low += Erange;

  /* renormalisation as for biari_encode_symbol; 
     note that low has already been doubled */ 
  if (low >= ONE)
  {
    put_one_bit_plus_outstanding(1);
    low -= ONE;
  }
  else 
    if (low < HALF)
    {
      put_one_bit_plus_outstanding(0);
    }
    else
    {
      Ebits_to_follow++;
      low -= HALF;
    }
    Elow = low;
    eep->C++;    
}

/*!
 ************************************************************************
 * \brief
 *    Arithmetic encoding for last symbol before termination
 ************************************************************************
 */
void biari_encode_symbol_final(EncodingEnvironmentPtr eep, signed short symbol)
{
  register unsigned int range = Erange-2;
  register unsigned int low = Elow;
  
#if (2==TRACE)
  extern int cabac_encoding;
  if (cabac_encoding)
    fprintf(p_trace, "%d  0x%04x\n", binCount++, Erange);
#endif
  
  if (symbol)
  {
    low += range;
    range = 2;
  }
  
  while (range < QUARTER)
  {
    if (low >= HALF)
    {
      put_one_bit_plus_outstanding(1);
      low -= HALF;
    }
    else 
      if (low < QUARTER)
      {
        put_one_bit_plus_outstanding(0);
      }
      else
      {
        Ebits_to_follow++;
        low -= QUARTER;
      }
      low <<= 1;
      range <<= 1;
  }
  Erange = range;
  Elow = low;
  eep->C++;
}


/*!
 ************************************************************************
 * \brief
 *    Initializes a given context with some pre-defined probability state
 ************************************************************************
 */
void biari_init_context (BiContextTypePtr ctx, const int* ini)
{
  int pstate;

  pstate = ((ini[0]* max(0, img->currentSlice->qp)) >> 4) + ini[1];
  pstate = min (max ( 1, pstate), 126);

  if ( pstate >= 64 )
  {
    ctx->state  = (unsigned short) (pstate - 64);
    ctx->MPS    = 1;
  }
  else
  {
    ctx->state  = (unsigned short) (63 - pstate);
    ctx->MPS    = 0;
  }
  
  ctx->count = 0;
}

#ifdef ADAPTIVE_LOOP_FILTER
void add_pic_bin_count(int val)
{
  pic_bin_count += val;
}
/*!
 ************************************************************************
 * \brief
 *    Initializes the EncodingEnvironment for the arithmetic coder
 ************************************************************************
 */
void arienco_start_encoding_AlfBlockFlag(EncodingEnvironmentPtr eep,
                            unsigned char *code_buffer,
                            int *code_len )
{
  Elow = 0;
  Ebits_to_follow = 0;
  Ebuffer = 0;
  Ebits_to_go = 9; // to swallow first redundant bit

  Ecodestrm = code_buffer;
  Ecodestrm_len = code_len;

  Erange = HALF-2;

  eep->C = 0;
  eep->E = 0;

}
#endif
