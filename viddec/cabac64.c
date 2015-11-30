#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "global.h"
#include "image.h"
#include "elements.h"
#include "errorconcealment.h"
//#include "macroblock.h"
#include "fmo.h"
#include "cabac.h"
#include "vlc.h"
#include "image.h"
#include "mb_access.h"
#include "biaridecod.h"
#include "transform8x8.h"


#ifdef MB32X32
#if TRACE
#define TRACE_STRING(s) strncpy(currSE.tracestring, s, TRACESTRING_SIZE)
#else
#define TRACE_STRING(s) // do nothing
#endif

extern int last_dquant;


void readCBP_CABAC_luma1bit(SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp); 
void readCBP_CABAC_chroma(SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
unsigned int unary_bin_decode(DecodingEnvironmentPtr dep_dp, BiContextTypePtr ctx, int ctx_offset);
unsigned int unary_exp_golomb_mv_decode(DecodingEnvironmentPtr dep_dp, BiContextTypePtr ctx, unsigned int max_bin);
int set_CBP_block_bit_16x16 (Macroblock *currMB);
int read_and_store_CBP_block_bit (Macroblock *currMB, DecodingEnvironmentPtr dep_dp, struct img_par *img, int type);
void read_significant_coefficients (Macroblock *currMB, DecodingEnvironmentPtr  dep_dp, struct img_par *img, int type, int coeff[]);
int read_significance_map (Macroblock *currMB, DecodingEnvironmentPtr  dep_dp, struct img_par *img, int type, int coeff[]);

void readLumaTransform_size(struct inp_par *inp)
{
  //int i,j,k;
  //int level;
  int mb_nr = img->current_mb_nr;
  //int ii,jj;
  //int m2,jg2;// i1,j1;
  Macroblock *currMB = &img->mb_data[mb_nr];
  //int cbp;
  SyntaxElement currSE;
  Slice *currSlice = img->currentSlice;
  DataPartition *dP;
  int *partMap = assignSE2partition[currSlice->dp_mode];
  //int iii,jjj;
  //int coef_ctr, i0, j0, b8;
  //int ll;
  //int block_x,block_y;
  //int start_scan;
  //int run, len;
  //int levarr[16], runarr[16], numcoeff;
  
  int smb       = ((img->type==SP_SLICE) && IS_INTER (currMB)) || (img->type == SI_SLICE && currMB->mb_type == SI4MB);
  
  //int uv;
  
  //int intra     = IS_INTRA (currMB);
  //int temp[4];
  
  //int b4;
  //int yuv = dec_picture->chroma_format_idc-1;
  
  //int need_transform_size_flag;
  
  
  if(img->type==SP_SLICE  && currMB->mb_type!=I16MB )
    smb=1;
  
    
    //if (need_transform_size_flag)
    {
      currSE.type   =  SE_HEADER;
      dP = &(currSlice->partArr[partMap[SE_HEADER]]);
      currSE.reading = readMB_transform_size_flag_CABAC;
      TRACE_STRING("transform size 8x8 flag");
      
      // read UVLC transform_size_8x8_flag
      if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
      {
        currSE.len = 1;
        readSyntaxElement_FLC(&currSE, dP->bitstream);
      } else 
      {
        dP->readSyntaxElement(&currSE,img,inp,dP);
      }
      currMB->luma_transform_size_8x8_flag = currSE.value1;

    }                 
}

void readlumaCBPTransform_size(struct inp_par *inp)
{
  //int i,j,k;
  //int level;
  int mb_nr = img->current_mb_nr;
  //int ii,jj;
  //int m2,jg2;// i1,j1;
  Macroblock *currMB = &img->mb_data[mb_nr];
  int cbp;
  SyntaxElement currSE;
  Slice *currSlice = img->currentSlice;
  DataPartition *dP;
  int *partMap = assignSE2partition[currSlice->dp_mode];
  //int iii,jjj;
  //int coef_ctr, i0, j0, b8;
  //int ll;
  //int block_x,block_y;
  //int start_scan;
  //int run, len;
  //int levarr[16], runarr[16], numcoeff;
  
  int smb       = ((img->type==SP_SLICE) && IS_INTER (currMB)) || (img->type == SI_SLICE && currMB->mb_type == SI4MB);
  
  //int uv;
  
  //int intra     = IS_INTRA (currMB);
  //int temp[4];
  
  //int b4;
  //int yuv = dec_picture->chroma_format_idc-1;
  
  int need_transform_size_flag;
  
  
  if(img->type==SP_SLICE  && currMB->mb_type!=I16MB )
    smb=1;
  
  
#ifdef ADAPTIVE_FD_SD_CODING
  currMB->SD_Coding_on_off=0;
#endif
  
    //=====   C B P   =====
    //---------------------
    if (IS_OLDINTRA (currMB) || currMB->mb_type == SI4MB || currMB->mb_type == I8MB)   currSE.type = SE_CBP_INTRA;
    else                        currSE.type = SE_CBP_INTER;
    
    dP = &(currSlice->partArr[partMap[currSE.type]]);
    
    if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
    {
      if (IS_OLDINTRA (currMB) || currMB->mb_type == SI4MB || currMB->mb_type == I8MB)  currSE.mapping = linfo_cbp_intra;
      else                        currSE.mapping = linfo_cbp_inter;
    }
    else
    {
      currSE.reading = readCBP_CABAC_luma1bit;
    }
    
    TRACE_STRING("coded_block_pattern");
    
    dP->readSyntaxElement(&currSE,img,inp,dP);
    currMB->cbp = cbp = currSE.value1;
    
#ifdef ADAPTIVE_FD_SD_CODING
    currMB->SD_Coding_on_off=1;
    if (cbp!=0 && !IS_INTRA (currMB) && img->Allow_SD_Coding)
    {
      currSE.type = SE_HEADER;
      dP = &(currSlice->partArr[partMap[currSE.type]]);
      currSE.reading = read_adaptive_prederror_coding_flag_for_MB;
      TRACE_STRING("adaptive_prederror_coding_flag");
      dP->readSyntaxElement(&currSE,img,inp,dP);
      currMB->SD_Coding_on_off = currSE.value1;
      currMB->written_SD_Coding_on_off=1;
    }
#endif
    

    //============= Transform size flag for INTER MBs =============
    //-------------------------------------------------------------
    need_transform_size_flag = (((currMB->mb_type >= 1 && currMB->mb_type <= 3)||
      (IS_DIRECT(currMB) && active_sps->direct_8x8_inference_flag) ||
      (currMB->NoMbPartLessThan8x8Flag))
      && currMB->mb_type != I8MB && currMB->mb_type != I4MB
      && (currMB->cbp&15)
      && img->Transform8x8Mode);
    
    if (need_transform_size_flag)
    {
      currSE.type   =  SE_HEADER;
      dP = &(currSlice->partArr[partMap[SE_HEADER]]);
      currSE.reading = readMB_transform_size_flag_CABAC;
      TRACE_STRING("transform size 8x8 flag");
      
      // read UVLC transform_size_8x8_flag
      if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
      {
        currSE.len = 1;
        readSyntaxElement_FLC(&currSE, dP->bitstream);
      } else 
      {
        dP->readSyntaxElement(&currSE,img,inp,dP);
      }
      currMB->luma_transform_size_8x8_flag = currSE.value1;

    }                 
}

int read_chroma_cbp(struct inp_par *inp)
{
  //int i,j,k;
  //int level;
  int mb_nr = img->current_mb_nr;
  //int ii,jj;
  //int m2,jg2;// i1,j1;
  Macroblock *currMB = &img->mb_data[mb_nr];
  int cbp;
  SyntaxElement currSE;
  Slice *currSlice = img->currentSlice;
  DataPartition *dP;
  int *partMap = assignSE2partition[currSlice->dp_mode];
  //int iii,jjj;
  //int coef_ctr, i0, j0, b8;
  //int ll;
  //int block_x,block_y;
  //int start_scan;
  //int run, len;
  //int levarr[16], runarr[16], numcoeff;
  
  int smb       = ((img->type==SP_SLICE) && IS_INTER (currMB)) || (img->type == SI_SLICE && currMB->mb_type == SI4MB);
  
  //int uv;
  
  //int intra     = IS_INTRA (currMB);
  //int temp[4];
  
  //int b4;
  //int yuv = dec_picture->chroma_format_idc-1;
  
  //int need_transform_size_flag;
  
  
  if(img->type==SP_SLICE  && currMB->mb_type!=I16MB )
    smb=1;
  
  
#ifdef ADAPTIVE_FD_SD_CODING
  currMB->SD_Coding_on_off=0;
#endif
  
    //=====   C B P   =====
    //---------------------
    if (IS_OLDINTRA (currMB) || currMB->mb_type == SI4MB || currMB->mb_type == I8MB)   currSE.type = SE_CBP_INTRA;
    else                        currSE.type = SE_CBP_INTER;
    
    dP = &(currSlice->partArr[partMap[currSE.type]]);
    
    if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
    {
      if (IS_OLDINTRA (currMB) || currMB->mb_type == SI4MB || currMB->mb_type == I8MB)  currSE.mapping = linfo_cbp_intra;
      else                        currSE.mapping = linfo_cbp_inter;
    }
    else
    {
      currSE.reading = readCBP_CABAC_chroma;
    }
    
    TRACE_STRING("coded_block_pattern");
    
    dP->readSyntaxElement(&currSE,img,inp,dP);
    cbp = currSE.value1;
    


    return cbp;
}

void readCBP_CABAC_luma1bit(SyntaxElement *se,
                   struct inp_par *inp,
                   struct img_par *img,
                   DecodingEnvironmentPtr dep_dp)
{
  TextureInfoContexts *ctx = img->currentSlice->tex_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int mb_x, mb_y;
  int a, b;
  int curr_cbp_ctx, curr_cbp_idx;
  int cbp = 0;
  int cbp_bit;
  int mask;
  PixelPos block_a;

  //  coding of luma part (bit by bit)
#ifdef MB32X32
  for (mb_y=0; mb_y < 2; mb_y += 2)
#else
  for (mb_y=0; mb_y < 4; mb_y += 2)
#endif
  {
#ifdef MB32X32
    for (mb_x=0; mb_x < 2; mb_x += 2)
#else
    for (mb_x=0; mb_x < 4; mb_x += 2)
#endif
    {
      if (currMB->b8mode[mb_y+(mb_x/2)]==IBLOCK)
        curr_cbp_idx = 0;
      else
        curr_cbp_idx = 1;

      if (mb_y == 0)
      {
        if (currMB->mb_available_up == NULL)
          b = 0;
        else
        {
          if((currMB->mb_available_up)->mb_type==IPCM)
            b=0;
          else
            b = (( ((currMB->mb_available_up)->cbp & (1<<(2+mb_x/2))) == 0) ? 1 : 0);
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
              a = (( (img->mb_data[block_a.mb_addr].cbp & (1<<(2*(block_a.y/2)+1))) == 0) ? 1 : 0);
          }
          
        }
        else
          a=0;
      }
      else
        a = ( ((cbp & (1<<mb_y)) == 0) ? 1: 0);

      curr_cbp_ctx = a+2*b;
      mask = (1<<(mb_y+mb_x/2));
      cbp_bit = biari_decode_symbol(dep_dp, ctx->cbp_contexts[0] + curr_cbp_ctx );
      if (cbp_bit) cbp += mask;
    }
  }


 /* if (dec_picture->chroma_format_idc != YUV400)
  {
    // coding of chroma part
    // CABAC decoding for BinIdx 0
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
    cbp_bit = biari_decode_symbol(dep_dp, ctx->cbp_contexts[1] + curr_cbp_ctx );
    
    // CABAC decoding for BinIdx 1 
    if (cbp_bit) // set the chroma bits
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
      cbp_bit = biari_decode_symbol(dep_dp, ctx->cbp_contexts[2] + curr_cbp_ctx );
      cbp += (cbp_bit == 1) ? 32 : 16;
    }
  }*/

  se->value1 = cbp;


#if TRACE
  fprintf(p_trace, "@%d %s\t\t\t%d\n",symbolCount++, se->tracestring, se->value1);
  fflush(p_trace);
#endif
}

void readCBP_CABAC_chroma(SyntaxElement *se,
                   struct inp_par *inp,
                   struct img_par *img,
                   DecodingEnvironmentPtr dep_dp)
{
  TextureInfoContexts *ctx = img->currentSlice->tex_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  //int mb_x, mb_y;
  int a, b;
  int curr_cbp_ctx;
  int cbp = 0;
  int cbp_bit;
  //int mask;
  //PixelPos block_a;

  if (dec_picture->chroma_format_idc != YUV400)
  {
    // coding of chroma part
    // CABAC decoding for BinIdx 0
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
    cbp_bit = biari_decode_symbol(dep_dp, ctx->cbp_contexts[1] + curr_cbp_ctx );
    
    // CABAC decoding for BinIdx 1 
    if (cbp_bit) // set the chroma bits
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
      cbp_bit = biari_decode_symbol(dep_dp, ctx->cbp_contexts[2] + curr_cbp_ctx );
      cbp += (cbp_bit == 1) ? 32 : 16;
    }
  }

  se->value1 = cbp;



#if TRACE
  fprintf(p_trace, "@%d %s\t\t\t%d\n",symbolCount++, se->tracestring, se->value1);
  fflush(p_trace);
#endif
}

/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the coded
 *    block pattern of a given MB.
 ************************************************************************
 */
void readCBP_CABAC32(SyntaxElement *se,
                   struct inp_par *inp,
                   struct img_par *img,
                   DecodingEnvironmentPtr dep_dp)
{
  TextureInfoContexts *ctx = img->currentSlice->tex_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int mb_x, mb_y;
  int a, b;
  int curr_cbp_ctx, curr_cbp_idx;
  int cbp = 0;
  int cbp_bit;
  int mask;
  PixelPos block_a;

  //  coding of luma part (bit by bit)
#ifdef MB32X32
  for (mb_y=0; mb_y < 2; mb_y += 2)
#else
  for (mb_y=0; mb_y < 4; mb_y += 2)
#endif
  {
#ifdef MB32X32
    for (mb_x=0; mb_x < 2; mb_x += 2)
#else
    for (mb_x=0; mb_x < 4; mb_x += 2)
#endif
    {
      if (currMB->b8mode[mb_y+(mb_x/2)]==IBLOCK)
        curr_cbp_idx = 0;
      else
        curr_cbp_idx = 1;

      if (mb_y == 0)
      {
        if (currMB->mb_available_up == NULL)
          b = 0;
        else
        {
          if((currMB->mb_available_up)->mb_type==IPCM)
            b=0;
          else
            b = (( ((currMB->mb_available_up)->cbp & (1<<(2+mb_x/2))) == 0) ? 1 : 0);
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
              a = (( (img->mb_data[block_a.mb_addr].cbp & (1<<(2*(block_a.y/2)+1))) == 0) ? 1 : 0);
          }
          
        }
        else
          a=0;
      }
      else
        a = ( ((cbp & (1<<mb_y)) == 0) ? 1: 0);

      curr_cbp_ctx = a+2*b;
      mask = (1<<(mb_y+mb_x/2));
      cbp_bit = biari_decode_symbol(dep_dp, ctx->cbp_contexts[0] + curr_cbp_ctx );
      if (cbp_bit) cbp += mask;
    }
  }


 /* if (dec_picture->chroma_format_idc != YUV400)
  {
    // coding of chroma part
    // CABAC decoding for BinIdx 0
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
    cbp_bit = biari_decode_symbol(dep_dp, ctx->cbp_contexts[1] + curr_cbp_ctx );
    
    // CABAC decoding for BinIdx 1 
    if (cbp_bit) // set the chroma bits
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
      cbp_bit = biari_decode_symbol(dep_dp, ctx->cbp_contexts[2] + curr_cbp_ctx );
      cbp += (cbp_bit == 1) ? 32 : 16;
    }
  }*/

  se->value1 = cbp;

  if (!cbp)
  {
    last_dquant=0;
  }

#if TRACE
  fprintf(p_trace, "@%d %s\t\t\t%d\n",symbolCount++, se->tracestring, se->value1);
  fflush(p_trace);
#endif
}

/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the reference
 *    parameter of a given MB.
 ************************************************************************
 */
void readRefFrame_CABAC32( SyntaxElement *se,
                         struct inp_par *inp,
                         struct img_par *img,
                         DecodingEnvironmentPtr dep_dp)
{
  MotionInfoContexts *ctx = img->currentSlice->mot_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int   addctx  = 0;
  int   a, b;
  int   act_ctx;
  int   act_sym;
  char** refframe_array = dec_picture->ref_idx[se->value2];
  int   b8a, b8b;

  PixelPos block_a, block_b;
#ifdef MB32X32
  int block_a_is_direct=0, block_b_is_direct=0, block_a_b8mode_is_direct=0, block_b_b8mode_is_direct=0;
  int mb_nr=0, subblock_x=0, subblock_y=0;
  int mb_ext_level = se->mb_ext_level;
  if(img->subblock_x==0 && img->subblock_y==0 )
  {
    subblock_x = 0;
    subblock_y = 0;
    mb_nr      = img->current_mb_nr;
  }
  else if(img->subblock_x==0 && img->subblock_y==(2<<mb_ext_level) )
  {
    subblock_x = 0;
    subblock_y = 0;
    mb_nr      = img->current_mb_nr + (img->PicWidthInMbs<<(mb_ext_level-1));
  }
  else if(img->subblock_x==(2<<mb_ext_level) && img->subblock_y==0 )
  {
    subblock_x = 0;
    subblock_y = 0;
    mb_nr      = img->current_mb_nr + (1<<(mb_ext_level-1));
  }
  else
  {
    error("subblock_x and subblock_y are out of range\n", -1);
  }

  getLuma4x4Neighbour(mb_nr, subblock_x, subblock_y, -1,  0, &block_a);
  getLuma4x4Neighbour(mb_nr, subblock_x, subblock_y,  0, -1, &block_b);
#else
  
  getLuma4x4Neighbour(img->current_mb_nr, img->subblock_x, img->subblock_y, -1,  0, &block_a);
  getLuma4x4Neighbour(img->current_mb_nr, img->subblock_x, img->subblock_y,  0, -1, &block_b);
#endif
  b8a=((block_a.x/2)%2)+2*((block_a.y/2)%2);
  b8b=((block_b.x/2)%2)+2*((block_b.y/2)%2);


#ifdef MB32X32
  if(block_b.available && (mb_nr == img->current_mb_nr || mb_nr == img->current_mb_nr+(1<<(mb_ext_level-1))) )//block_b should be out of the current 32x32 block, otherwise it is impossible to be direct mode
  {
    //int bslice = (img->type==B_SLICE);
    block_b_is_direct = IS_DIRECT(&img->mb_data[block_b.mb_addr]);
    block_b_b8mode_is_direct = img->mb_data[block_b.mb_addr].b8mode[b8b]==0 && img->mb_data[block_b.mb_addr].b8pdir[b8b]==2;
  }
  if(block_a.available && (mb_nr == img->current_mb_nr || mb_nr == img->current_mb_nr + (img->PicWidthInMbs << (mb_ext_level-1))) )//block_a should be out of the current 32x32 block, otherwise it is impossible to be direct mode
  {     
    //int bslice = (img->type==B_SLICE);
    block_a_is_direct = IS_DIRECT(&img->mb_data[block_a.mb_addr]);
    block_a_b8mode_is_direct = img->mb_data[block_a.mb_addr].b8mode[b8a]==0 && img->mb_data[block_a.mb_addr].b8pdir[b8a]==2;
  }
#endif

  if (!block_b.available)
    b=0;
#ifdef MB32X32
  else if ( block_b_is_direct || block_b_b8mode_is_direct )
#else
  else if ( (img->mb_data[block_b.mb_addr].mb_type==IPCM) || IS_DIRECT(&img->mb_data[block_b.mb_addr]) || (img->mb_data[block_b.mb_addr].b8mode[b8b]==0 && img->mb_data[block_b.mb_addr].b8pdir[b8b]==2))
#endif
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
#ifdef MB32X32
  else if ( block_a_is_direct || block_a_b8mode_is_direct)
#else
  else if ((img->mb_data[block_a.mb_addr].mb_type==IPCM) || IS_DIRECT(&img->mb_data[block_a.mb_addr]) || (img->mb_data[block_a.mb_addr].b8mode[b8a]==0 && img->mb_data[block_a.mb_addr].b8pdir[b8a]==2))
#endif
    a=0;
  else 
  {
    if (img->MbaffFrameFlag && (currMB->mb_field == 0) && (img->mb_data[block_a.mb_addr].mb_field == 1))
      a = (refframe_array[block_a.pos_y][block_a.pos_x] > 1 ? 1 : 0);
    else
      a = (refframe_array[block_a.pos_y][block_a.pos_x] > 0 ? 1 : 0);
  }

  act_ctx = a + 2*b;
  se->context = act_ctx; // store context

  act_sym = biari_decode_symbol(dep_dp,ctx->ref_no_contexts[addctx] + act_ctx );

  if (act_sym != 0)
  {
    act_ctx = 4;
    act_sym = unary_bin_decode(dep_dp,ctx->ref_no_contexts[addctx]+act_ctx,1);
    act_sym++;
  }
  se->value1 = act_sym;


#if TRACE
  fprintf(p_trace, "@%d %s\t\t\t%d \n",symbolCount++, se->tracestring, se->value1);
//  fprintf(p_trace," c: %d :%d \n",ctx->ref_no_contexts[addctx][act_ctx].cum_freq[0],ctx->ref_no_contexts[addctx][act_ctx].cum_freq[1]);
  fflush(p_trace);
#endif
}

/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the motion
 *    vector data of a B-frame MB.
 ************************************************************************
 */
void readMVD_CABAC32( SyntaxElement *se,
                    struct inp_par *inp,
                    struct img_par *img,
                    DecodingEnvironmentPtr dep_dp)
{
  //int i = img->subblock_x;
  //int j = img->subblock_y;
  int a, b;
  int act_ctx;
  int act_sym;
  int mv_local_err;
  int mv_sign;
  int list_idx = se->value2 & 0x01;
  int k = (se->value2>>1); // MVD component

  PixelPos block_a, block_b;

  MotionInfoContexts *ctx = img->currentSlice->mot_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

#ifdef MB32X32
  int mb_nr=0, subblock_x=0, subblock_y=0; 
  int mb_ext_level = se->mb_ext_level;
  if(img->subblock_x==0 && img->subblock_y==0 )
  {
    subblock_x = 0;
    subblock_y = 0;
    mb_nr      = img->current_mb_nr;
  }
  else if(img->subblock_x==0 && img->subblock_y==(2<<mb_ext_level) )
  {
    subblock_x = 0;
    subblock_y = 0;
    mb_nr      = img->current_mb_nr + (img->PicWidthInMbs<<(mb_ext_level-1));
  }
  else if(img->subblock_x==(2<<mb_ext_level) && img->subblock_y==0 )
  {
    subblock_x = 0;
    subblock_y = 0;
    mb_nr      = img->current_mb_nr + (1<<(mb_ext_level-1));
  }
  else
  {
    error("subblock_x and subblock_y are out of the range\n", -1);
  }

  getLuma4x4Neighbour(mb_nr, subblock_x, subblock_y, -1,  0, &block_a);
  getLuma4x4Neighbour(mb_nr, subblock_x, subblock_y,  0, -1, &block_b);
#else
  getLuma4x4Neighbour(img->current_mb_nr, i, j, -1,  0, &block_a);
  getLuma4x4Neighbour(img->current_mb_nr, i, j,  0, -1, &block_b);
#endif

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
  se->context = act_ctx;

  act_sym = biari_decode_symbol(dep_dp,&ctx->mv_res_contexts[0][act_ctx] );

  if (act_sym != 0)
  {
    act_ctx=5*k;
    act_sym = unary_exp_golomb_mv_decode(dep_dp,ctx->mv_res_contexts[1]+act_ctx,3);
    act_sym++;
    mv_sign = biari_decode_symbol_eq_prob(dep_dp);

    if(mv_sign)
      act_sym = -act_sym;
  }
  se->value1 = act_sym;


#if TRACE
  fprintf(p_trace, "@%d %s\t\t\t%d \n",symbolCount++, se->tracestring, se->value1);
  fflush(p_trace);
#endif
}

void readRunLevel_CABAC16x16 (SyntaxElement  *se,
                         struct inp_par *inp,
                         struct img_par *img,
                         DecodingEnvironmentPtr dep_dp)
{
#ifdef USE_INTRA_MDDT
  static int  coeff[256]; // one more for EOB
#else
  static int  coeff[64]; // one more for EOB
#endif
  static int  coeff_ctr = -1;
  static int  pos       =  0;
  int         DF_or_SD=0;
  int         i;
  int         writing_b8=min(3,max(0,se->writing_b8));
  int         writing_b4=min(3,max(0,se->writing_b4));

  Macroblock  *currMB = &img->mb_data[img->current_mb_nr];

  int         isBigBlock; 

  isBigBlock = currMB->luma_transform_size_8x8_flag == 2 && (se->context != CHROMA_AC && se->context != CHROMA_DC) ? 1:0;

  if(isBigBlock)
  {
#ifdef MB32X32
    set_CBP_block_bit_16x16(currMB);
#else
    if(currMB->mb_type == 0 || currMB->mb_type == P16x16)
      set_CBP_block_bit_16x16(currMB);
    else if(currMB->mb_type == P16x8)
      set_CBP_block_bit_16x8 (currMB);
    else
      set_CBP_block_bit_8x16 (currMB);
#endif

  }

  //--- read coefficients for whole block ---
  if (coeff_ctr < 0)
  {
    //===== decode CBP-BIT =====
    if (isBigBlock || (coeff_ctr = read_and_store_CBP_block_bit (currMB, dep_dp, img, se->context)))
    {
      if (currMB->SD_Coding_on_off==1 && img->Allow_SD_Coding && currMB->b8mode[writing_b8] != IBLOCK && currMB->luma_transform_size_8x8_flag==0 && se->context==LUMA_4x4 && !IS_INTRA (currMB) )
      {
        DF_or_SD = read_spatial_domain_coding_flag_for_4x4block (currMB, dep_dp, img, se->context);

        if (DF_or_SD)
        {
          //===== decode significance map =====
          coeff_ctr = read_significance_map_SD (currMB, dep_dp, img, se->context, coeff);

          //===== decode significant coefficients =====
          read_significant_coefficients_SD     (currMB, dep_dp, img, se->context, coeff);

          for (i=0; i<16; i++)
          {
            currMB->quantizer_indices[((writing_b8/ 2)<<3)+((writing_b4/ 2)<<2)+(i/4)][((writing_b8% 2)<<3)+((writing_b4% 2)<<2)+(i%4)]=coeff[i];
          }
          currMB->SD_or_FD[writing_b8/ 2][writing_b8% 2]+=(1<<writing_b4);
        }
        else
        {
          //===== decode significance map =====
          coeff_ctr = read_significance_map (currMB, dep_dp, img, se->context, coeff);

          //===== decode significant coefficients =====
          read_significant_coefficients     (currMB, dep_dp, img, se->context, coeff);
        }
      }

      else if (currMB->SD_Coding_on_off==1 && img->Allow_SD_Coding && currMB->b8mode[writing_b8] != IBLOCK && currMB->luma_transform_size_8x8_flag==1 && se->context==LUMA_8x8 && !IS_INTRA (currMB) )
      {
        if ((DF_or_SD = read_spatial_domain_coding_flag_for_8x8block (currMB, dep_dp, img, se->context,writing_b8)))
        {
          //===== decode significance map =====
          coeff_ctr = read_significance_map_SD8 (currMB, dep_dp, img, se->context, coeff);

          //===== decode significant coefficients =====
          read_significant_coefficients_SD8     (currMB, dep_dp, img, se->context, coeff);
          for (i=0; i<64; i++)
          {
            currMB->quantizer_indices[((writing_b8/ 2)<<3)+(i/8)][((writing_b8% 2)<<3)+(i%8)]=coeff[i];
          }
          currMB->SD_or_FD_t8x8+=(1<<writing_b8);
        }
        else
        {
          //===== decode significance map =====
          coeff_ctr = read_significance_map (currMB, dep_dp, img, se->context, coeff);

          //===== decode significant coefficients =====
          read_significant_coefficients     (currMB, dep_dp, img, se->context, coeff);
        }
      }

      else
      {
        //===== decode significance map =====
        coeff_ctr = read_significance_map (currMB, dep_dp, img, se->context, coeff);

        //===== decode significant coefficients =====
        read_significant_coefficients     (currMB, dep_dp, img, se->context, coeff);
      }
    }
  }

  //--- set run and level ---
  if (coeff_ctr)
  {
    //--- set run and level (coefficient) ---
    for (se->value2=0; coeff[pos]==0; pos++, se->value2++);
    se->value1=coeff[pos++];
  }
  else
  {
    //--- set run and level (EOB) ---
    se->value1 = se->value2 = 0;
  }
  //--- decrement coefficient counter and re-set position ---
  if (coeff_ctr-- == 0) pos=0;

#if TRACE
  fprintf(p_trace, "@%d %s\t\t\t%d\t%d\n",symbolCount++, se->tracestring, se->value1,se->value2);
  fflush(p_trace);
#endif
}

#endif


