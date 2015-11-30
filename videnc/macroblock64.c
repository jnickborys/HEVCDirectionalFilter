
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include "global.h"

#include "elements.h"
#include "refbuf.h"
#include "fmo.h"
#include "vlc.h"
#include "image.h"
#include "mb_access.h"
#include "ratectl.h"              // head file for rate control
#include "cabac.h"
#include "transform8x8.h"
#include "header.h"
#ifdef ADAPTIVE_FILTER
#include "adaptive_filter.h"
#endif
#include <math.h>
#ifdef MV_COMPETITION
#include "mv_competition.h"
extern MV_Competition mv_comp;
#ifdef MB32X32_MVC
extern int ******send_index_mv_prediction;        // int predModeMV2[mode][10][list][4][4][MAX_MODEPREDMV+1]; 
#endif
#endif

#ifdef MB32X32
#define DCT16PREC   7

#define SIMPLIFY_CODE

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#ifdef ADAPTIVE_FD_SD_CODING
extern int cabac_encoding;
#endif

extern Macroblock  MB32;


extern Macroblock  MB64;//MB64X64
extern int   cbp_32[1<<(2*(MAX_MB_EXT_LEVEL-1))], cbp_32_RDCost[1<<(2*(MAX_MB_EXT_LEVEL-1))];//MB64X64
extern int write_mb64_mbtype;//MB64X64

extern int   cbp16[1<<(2*MAX_MB_EXT_LEVEL)];
extern int64 cbp_blk16[1<<(2*MAX_MB_EXT_LEVEL)];
extern int   ****cofAC32[1<<(2*MAX_MB_EXT_LEVEL)];     // [8x8block][4x4block][level/run][scan_pos]
extern int   ***cofDC32[1<<(2*MAX_MB_EXT_LEVEL)];                       // [yuv][level/run][scan_pos]
extern int   luma_transform_size_8x8_flag32[1<<(2*MAX_MB_EXT_LEVEL)];
extern int    **cof32AC16x16[1<<(2*MAX_MB_EXT_LEVEL)];

extern int bestmode32;
extern int bestmode64;

extern Macroblock  MB64;//MB64X64
int mb64_stat[2][5];//MB64X64
extern int mb32type[1<<(2*(MAX_MB_EXT_LEVEL-1))];//MB64X64


extern int delta_qp_sent;
extern int eos_bit_32;
extern int mb32_mvd[2][BLOCK_MULTIPLE<<MAX_MB_EXT_LEVEL][BLOCK_MULTIPLE<<MAX_MB_EXT_LEVEL][2];          //!< indices correspond to [forw,backw][block_y][block_x][x,y]
void start_macroblock32(int mb_addr, int mb_field, int mb32);
void encode_one_macroblock32 (int level);
int write_one_macroblock32 (int eos_bit, int mb32, int mbi);
int write_one_macroblock64 (int eos_bit, int mb32i, int mb16i);
void unary_bin_encode(EncodingEnvironmentPtr eep_frame, unsigned int symbol, BiContextTypePtr ctx, int ctx_offset);
void unary_exp_golomb_mv_encode(EncodingEnvironmentPtr eep_dp, unsigned int symbol, BiContextTypePtr ctx, unsigned int max_bin);
int MBType2Value (Macroblock* currMB);
int writeIntraModes(void);
int ZeroRef (Macroblock* currMB);
int BType2CtxRef (int btype);

/*!
 ************************************************************************
 * \brief
 *    level           : 1-32x32; 2-64x64
 *    CurrentMbAddr   : 16x16MB index in the raster scan order starting from 0
 *    CurrentMbAddr_x : 16x16MB x-axis index
 *    CurrentMbAddr_y : 16x16MB y_axis index
 *    extend_x        : the horizonal next 16x16MB x-axis index in the current cluster
 *    extend_y        : the vertical  next 16x16MB y-axis index in the current cluster
 *    num_MB          : number of 16x16 MB in the cluster
 *
 ************************************************************************
 */
void find_Mb_cluster(int level, int CurrentMbAddr, int *CurrentMbAddr_x, int *CurrentMbAddr_y, int *extend_x, int *extend_y, int *num_MB, int *nextMbAddr)
{
  int extend = 1<<level;

  *CurrentMbAddr_x = CurrentMbAddr % img->PicWidthInMbs;
  *CurrentMbAddr_y = CurrentMbAddr / img->PicWidthInMbs;

  *extend_x = min((unsigned)(*CurrentMbAddr_x+extend-1), img->PicWidthInMbs-1);
  *extend_y = min((unsigned)(*CurrentMbAddr_y+extend-1), img->FrameHeightInMbs-1);

  *num_MB = (*extend_x - *CurrentMbAddr_x + 1) * (*extend_y - *CurrentMbAddr_y + 1);

  if((unsigned)(*extend_x+1) <= img->PicWidthInMbs-1)
    *nextMbAddr = CurrentMbAddr + extend;
  else if((unsigned)(*extend_y+1) <= img->FrameHeightInMbs-1)
    *nextMbAddr = (*extend_y+1) * img->PicWidthInMbs;
  else
    *nextMbAddr = -1;

}

void SetFrameOffset32(int CurrentMbAddr)
{
    int i, j;
    int encodeMbAddr;
    int CurrentMbAddr_x, CurrentMbAddr_y, extend_x, extend_y, num_MB, nextMbAddr;
 
    find_Mb_cluster(input->UseExtMB, CurrentMbAddr, &CurrentMbAddr_x, &CurrentMbAddr_y, &extend_x, &extend_y, &num_MB, &nextMbAddr);//MB64X64

    for(i=CurrentMbAddr_y; i<=extend_y; i++)
    {
      for(j=CurrentMbAddr_x; j<=extend_x; j++)
      {
        encodeMbAddr = i*img->PicWidthInMbs + j;
        set_MB_parameters(encodeMbAddr);
        SetFrameOffset();
      }
    }
}
#ifdef ADAPTIVE_FILTER 
void SetAdaptiveFilter(void);
void SetAdaptiveFilterHor(void);
void SetAdaptiveFilter32(int CurrentMbAddr)
{
    int i, j;
    int encodeMbAddr;
    int CurrentMbAddr_x, CurrentMbAddr_y, extend_x, extend_y, num_MB, nextMbAddr;
    find_Mb_cluster(input->UseExtMB, CurrentMbAddr, &CurrentMbAddr_x, &CurrentMbAddr_y, &extend_x, &extend_y, &num_MB, &nextMbAddr);

    for(i=CurrentMbAddr_y; i<=extend_y; i++)
    {
      for(j=CurrentMbAddr_x; j<=extend_x; j++)
      {
        encodeMbAddr = i*img->PicWidthInMbs + j;
        set_MB_parameters (encodeMbAddr);

#ifdef EAIF
          if ((input->UseAdaptiveFilter == 1) || (input->UseAdaptiveFilter == 3) || (input->UseAdaptiveFilter == FILTER_TYPE_EDAIF) || (input->UseAdaptiveFilter == FILTER_TYPE_EAIF))
#else
          if ((input->UseAdaptiveFilter == 1) || ((input->UseAdaptiveFilter == 3) || (input->UseAdaptiveFilter == FILTER_TYPE_EDAIF)))
#endif
            SetAdaptiveFilter();
          else
            SetAdaptiveFilterHor(); // separable aif
      }
    }
}
#endif

void prepare_data_for_write(int encodeMbAddr, Macroblock *MB, int mbcount, int mode)
{
  int ****i4p, ***i3p;
  Macroblock  *currMB  = &img->mb_data[encodeMbAddr];
  int **i2p;

  currMB->cbp = cbp16[mbcount];
  currMB->cbp_blk = cbp_blk16[mbcount];
  currMB->mb_type = mode;
  currMB->luma_transform_size_8x8_flag = luma_transform_size_8x8_flag32[mbcount];

  i4p=cofAC32[mbcount]; cofAC32[mbcount]=img->cofAC; img->cofAC=i4p;
  i3p=cofDC32[mbcount]; cofDC32[mbcount]=img->cofDC; img->cofDC=i3p;
  i2p=cof32AC16x16[mbcount]; cof32AC16x16[mbcount]=img->cofAC16x16; img->cofAC16x16=i2p;

  memcpy(currMB->b8mode,MB->b8mode, BLOCK_MULTIPLE * sizeof(int));
  memcpy(currMB->b8pdir,MB->b8pdir, BLOCK_MULTIPLE * sizeof(int));    
  currMB->bi_pred_me = MB->bi_pred_me;
}

void assign_qp_info(int mbcount, int encodeMbAddr, int saved_cluster_qp, int saved_prev_qp, int saved_prev_delta_qp)
{
  //in RDOQ, currMB->qp and currMB->delta_qp may change if cbp=0;
//#ifdef MB32_DELTA_QP
    img->mb_data[encodeMbAddr].qp = img->qp = saved_cluster_qp;
    set_chroma_qp(&img->mb_data[encodeMbAddr]);
//#endif

    img->mb_data[encodeMbAddr].prev_qp = saved_prev_qp;
    img->mb_data[encodeMbAddr].prev_delta_qp = saved_prev_delta_qp; //in writeDquant_CABAC, last_dquant=currMB->prev_delta_qp
    img->mb_data[encodeMbAddr].delta_qp = img->mb_data[encodeMbAddr].qp - img->mb_data[encodeMbAddr].prev_qp;
}

void propagate_qp_info(int delta_qp_sent_flag, int CurrentMbAddr, int CurrentMbAddr_x, int CurrentMbAddr_y, int extend_x, int extend_y,
                       int saved_cluster_qp, int saved_prev_qp, int saved_prev_delta_qp)      
{
  int i, j, encodeMbAddr;

  for(i=CurrentMbAddr_y; i<=extend_y; i++)
  {
    for(j=CurrentMbAddr_x; j<=extend_x; j++)
    {
      encodeMbAddr = i*img->PicWidthInMbs + j;

      img->mb_data[encodeMbAddr].prev_qp = saved_prev_qp;
      img->mb_data[encodeMbAddr].prev_delta_qp = saved_prev_delta_qp;

      if(delta_qp_sent_flag == 1)
      {
        img->mb_data[encodeMbAddr].qp = img->qp = saved_cluster_qp;
      }
      else
      {
        img->mb_data[encodeMbAddr].qp = saved_prev_qp;
      }
      set_chroma_qp(&img->mb_data[encodeMbAddr]);

      img->mb_data[encodeMbAddr].delta_qp = img->mb_data[encodeMbAddr].qp - img->mb_data[encodeMbAddr].prev_qp;
    }
  }

}
      
void propagate_mb_info(Macroblock* currMB32, int CurrentMbAddr, int CurrentMbAddr_x, int CurrentMbAddr_y, int extend_x, int extend_y,
                       int saved_cluster_qp, int saved_prev_qp, int saved_prev_delta_qp)      
{
  int i, j, encodeMbAddr;
  //set 16x16 block type 
  for(i=CurrentMbAddr_y; i<=extend_y; i++)
  {
    for(j=CurrentMbAddr_x; j<=extend_x; j++)
    {
      encodeMbAddr = i*img->PicWidthInMbs + j;

      if(img->type==P_SLICE && currMB32->mb_type == 0)
        assert(currMB32->cbp==0);

      if(currMB32->mb_type == 0 && currMB32->cbp==0)
        img->mb_data[encodeMbAddr].skip_flag = 1;
      else
        img->mb_data[encodeMbAddr].skip_flag = 0;

      if(img->mb_data[CurrentMbAddr].mb_type == 0)
        img->mb_data[encodeMbAddr].mb_type = 0;
      else
      {
        if( CurrentMbAddr != encodeMbAddr)
        {
          img->mb_data[encodeMbAddr].mb_type = 1;
        }
      }
//#ifdef MB32_DELTA_QP
  //reset mb qp
      img->mb_data[encodeMbAddr].prev_qp = saved_prev_qp;
      img->mb_data[encodeMbAddr].prev_delta_qp = saved_prev_delta_qp;

      if(currMB32->cbp != 0)
      {
        img->mb_data[encodeMbAddr].qp = img->qp = saved_cluster_qp;            
      }
      else
      {
        img->mb_data[encodeMbAddr].qp = saved_prev_qp;
      }         
      set_chroma_qp(&img->mb_data[encodeMbAddr]);

      img->mb_data[encodeMbAddr].delta_qp = img->mb_data[encodeMbAddr].qp - img->mb_data[encodeMbAddr].prev_qp;
//#endif
    }//j
  }//i
}

int write_macroblock_cluster32(int *NumberOfCodedMBs, int CurrentMbAddr, int mb32count, 
                               int saved_prev_qp, int saved_prev_delta_qp, int saved_cluster_qp, int *delta_qp_sent)
{    
    int no_bits=0;
    int i, j;
    int encodeMbAddr;
    int CurrentMbAddr_x, CurrentMbAddr_y, extend_x, extend_y, num_MB, nextMbAddr;

    CalculateQuantParam();
    CalculateOffsetParam();

    if(input->Transform8x8Mode)
    {
      CalculateQuant8Param();
      CalculateOffset8Param();
    }
 
    find_Mb_cluster(1, CurrentMbAddr, &CurrentMbAddr_x, &CurrentMbAddr_y, &extend_x, &extend_y, &num_MB, &nextMbAddr);

    //prepare 32x32 block type and 32x32 block cbp
    //bestmode32 = 1;//P8x8;

#ifdef RDO_Q
    if(input->UseRDO_Q)
      bestmode32 = saveBestMode.MBsize32[mb32count];
#endif          
    MB32.mb_type = bestmode32;
    MB32.cbp=0;


    //32x32, 32x16, 16x32
   //if(0)
   if(num_MB == 4 && bestmode32 != P8x8)
    {      
      int mbcount;
      short save_valid_mode[4];

      save_valid_mode[0]     = input->InterSearchSkip;      
      save_valid_mode[1]     = input->InterSearch16x16;//use for 32x32
      save_valid_mode[2]     = input->InterSearch16x8; //use for 32x16
      save_valid_mode[3]     = input->InterSearch8x16; //use for 16x32

      // if mb32count != 0 && bestmode32 == 0 we have to redo mode decision,
      // because the neighboring mv could change so the mv of skipped MB may be messed up.
      if(mb32count == 0 || bestmode32 != 0) 
      {
        input->InterSearchSkip  = 0;
        input->InterSearch16x16 = 0;
        input->InterSearch16x8  = 0;
        input->InterSearch8x16  = 0;
      }
      switch(bestmode32)
      {
      case 0:
        input->InterSearchSkip = 1;
        break;
      case 1:
        input->InterSearch16x16 = 1;
        break;
      case 2:
        input->InterSearch16x8 = 1;
        break;
      case 3:
        input->InterSearch8x16 = 1;
        break;
      default: 
        error("this mode is not correct!\n", -1);
        break;
      }

      start_macroblock32 (CurrentMbAddr, FALSE, 1);

      //qp handling
      assign_qp_info(mb32count,  CurrentMbAddr, saved_cluster_qp, saved_prev_qp, saved_prev_delta_qp);//MB32_DELTA_QP

      encode_one_macroblock32 (1);

      bestmode32=img->mb_data[img->current_mb_nr].mb_type; 


      mbcount=0;

      memset(mb32_mvd, 0, 2*(BLOCK_MULTIPLE<<MAX_MB_EXT_LEVEL)*(BLOCK_MULTIPLE<<MAX_MB_EXT_LEVEL)*2 * sizeof(int));//mb32_mvd was assigned only on predicted direction

      for(i=CurrentMbAddr_y; i<=extend_y; i++)
      {
        for(j=CurrentMbAddr_x; j<=extend_x; j++)
        {
          int b4_x,b4_y, list, k;
          //Macroblock *currMB;
          encodeMbAddr = i*img->PicWidthInMbs + j;

          start_macroblock32 (encodeMbAddr, FALSE, 0);//reset img->current_mb_nr and other global variables

          /******************************************
           *for the upper left 16x16 block, set b8mode
           *b8pdir for the whole 32x32 block in order
           *to write motion info in 1 shot
           *****************************************/
          prepare_data_for_write(encodeMbAddr, &MB32, mbcount, bestmode32);

          no_bits += write_one_macroblock32(1, 1, mbcount);

          if(MB32.mb_type != 0)
          {
            for(list=0; list<2; list++)            
              for(b4_y=0; b4_y<4; b4_y++)
                for(b4_x=0; b4_x<4; b4_x++)
                  for(k=0; k<2; k++)
                  {                    
                    img->mb_data[encodeMbAddr].mvd[list][b4_y][b4_x][k] = mb32_mvd[list][(mbcount>>1)*4+b4_y][(mbcount&1)*4+b4_x][k];
                  }
          }

          mbcount ++;

        }
     }
      if(MB32.cbp != 0)//MB32_DELTA_QP
        *delta_qp_sent = 1;
      /*
       * propagate information to each 16x16 of 32x32 block
       */      
      propagate_mb_info(&MB32, CurrentMbAddr, CurrentMbAddr_x,CurrentMbAddr_y,extend_x, extend_y, 
        saved_cluster_qp, saved_prev_qp, saved_prev_delta_qp);
      
      /*
       * recover input->InterSearch
       */
      input->InterSearchSkip  = save_valid_mode[0];      
      input->InterSearch16x16 = save_valid_mode[1];
      input->InterSearch16x8  = save_valid_mode[2];
      input->InterSearch8x16  = save_valid_mode[3];

      img->currentSlice->num_mb += 4;
      (*NumberOfCodedMBs) += 4;

    }//mode 32x32, 32x16, 16x32
    else
    {
      int mbcount=0;
      //int mb32_16x16_cbp_2nd[4]; //each 16x16 block cbp of 32x32 in the 2nd pass coding
      
      //16x16 mode
      for(i=CurrentMbAddr_y; i<=extend_y; i++)
      {
        for(j=CurrentMbAddr_x; j<=extend_x; j++)
        {
          /****************
           * initialization
           ****************/
          encodeMbAddr = i*img->PicWidthInMbs + j;
          start_macroblock (encodeMbAddr, FALSE);          
          assign_qp_info(mbcount,  encodeMbAddr, saved_cluster_qp, saved_prev_qp, saved_prev_delta_qp);

          /*********
           * coding
           *********/
          encode_one_macroblock ();

          //mb32_16x16_cbp_2nd[mbcount] = img->mb_data[encodeMbAddr].cbp;
          //if(img->mb_data[encodeMbAddr].cbp != 0)
          //{
          //  MB32.cbp = 1;
          //}
          if(num_MB == 4)
          {            
            no_bits += write_one_macroblock32(1, 0, mbcount); //32x32 block at 16x16 mode
          }
          else
            write_one_macroblock(1);


          proceed2nextMacroblock();
          img->currentSlice->num_mb++;
          (*NumberOfCodedMBs)++;

          if(img->mb_data[encodeMbAddr].cbp > 0 || img->mb_data[encodeMbAddr].mb_type == I16MB)//MB32_DELTA_QP
          {
            *delta_qp_sent = 1;
          }
          mbcount++;
          //total_cost 
        }//j
      }//i
         

      /*if(num_MB ==4 && MB32.cbp != MB32_16x16.cbp)
      {
        int mbi;
        printf("============32x32 block 16x16 mode mb32_cbp(2nd pass) != MB32_16x16.cbp (1st pass)\n");
        for(mbi=0; mbi<4; mbi++)
        {
          printf("mbi %d 1stpass %d 2ndpass %d\n", mbi, mb32_16x16_cbp[mbi], mb32_16x16_cbp_2nd[mbi]);
        }
      }*/

      /*************
       * reset mb qp
       *************/
      propagate_qp_info(*delta_qp_sent, CurrentMbAddr, CurrentMbAddr_x,CurrentMbAddr_y,extend_x, extend_y, 
        saved_cluster_qp, saved_prev_qp, saved_prev_delta_qp);//MB32_DELTA_QP
    }


    img->currentSlice->num_mb--;//currSlice->num_mb will be updated for the last mb of this cluster in terminate_macroblock  

    return nextMbAddr;
}

 
int write_macroblock_cluster64(int *NumberOfCodedMBs, int CurrentMbAddr, 
                               int saved_prev_qp, int saved_prev_delta_qp, int saved_cluster_qp, int *delta_qp_sent)//MB64X64
{
    static int first_time=1;
    int no_bits=0;
    int i, j, ii, jj;
    int encodeMbAddr;
    int CurrentMbAddr_x, CurrentMbAddr_y, extend_x, extend_y, num_MB, nextMbAddr;

#ifndef SIMPLIFY_CODE
    CalculateQuantParam();
    CalculateOffsetParam();

    if(input->Transform8x8Mode)
    {
      CalculateQuant8Param();
      CalculateOffset8Param();
    }
#endif

    find_Mb_cluster(2, CurrentMbAddr, &CurrentMbAddr_x, &CurrentMbAddr_y, &extend_x, &extend_y, &num_MB, &nextMbAddr);

    //bestmode64 = P8x8;
//SPEEDUP
#ifdef RDO_Q
    if(input->UseRDO_Q)
      bestmode64 = saveBestMode.MBsize64;
#endif

    MB64.mb_type = bestmode64;

    //collect statistics
    if(first_time==1)
    {
      memset(mb64_stat[0], 0, sizeof(int)*5);
      memset(mb64_stat[1], 0, sizeof(int)*5);
      first_time++;
    }
    if(num_MB==16)
    {
      mb64_stat[img->type][bestmode64<P8x8?bestmode64:4]++;
    }

    //64x64, 64x32, 32x64
   //if(0)
   if(num_MB == 16 && bestmode64 != P8x8)
    {      
      int mb16count, mb32count;
      short save_valid_mode[4];

      save_valid_mode[0]     = input->InterSearchSkip;      
      save_valid_mode[1]     = input->InterSearch16x16;//use for 32x32
      save_valid_mode[2]     = input->InterSearch16x8; //use for 32x16
      save_valid_mode[3]     = input->InterSearch8x16; //use for 16x32

      input->InterSearchSkip  = 0;
      input->InterSearch16x16 = 0;
      input->InterSearch16x8  = 0;
      input->InterSearch8x16  = 0;
      
      switch(bestmode64)
      {
      case 0:
        input->InterSearchSkip = 1;
        break;
      case 1:
        input->InterSearch16x16 = 1;
        break;
      case 2:
        input->InterSearch16x8 = 1;
        break;
      case 3:
        input->InterSearch8x16 = 1;
        break;
      default: 
        error("this mode is not correct!\n", -1);
        break;
      }

      start_macroblock32 (CurrentMbAddr, FALSE, 2);

      assign_qp_info(0,  CurrentMbAddr, saved_cluster_qp, saved_prev_qp, saved_prev_delta_qp);//overwrite what have been done in start_macroblock for qp info//MB32_DELTA_QP

      encode_one_macroblock32 (2);



      memset(mb32_mvd, 0, 2*(BLOCK_MULTIPLE<<MAX_MB_EXT_LEVEL)*(BLOCK_MULTIPLE<<MAX_MB_EXT_LEVEL)*2 * sizeof(int));//mb32_mvd was assigned only on predicted direction

      mb32count=0;
      for(i=CurrentMbAddr_y; i<=extend_y; i+=2)
      {
        for(j=CurrentMbAddr_x; j<=extend_x; j+=2)
        {
          MB32.cbp = cbp_32[mb32count];
          mb16count=0;
          for(ii=0; ii<2; ii++)
          {
            for(jj=0; jj<2; jj++)
            {
              int b4_x,b4_y, list, k;
              //Macroblock *currMB;
              encodeMbAddr = (i+ii)*img->PicWidthInMbs + (j+jj);

              start_macroblock32 (encodeMbAddr, FALSE, 0);//reset img->current_mb_nr and other global variables

              /******************************************
               *for the upper left 16x16 block, set b8mode
               *b8pdir for the whole 32x32 block in order
               *to write motion info in 1 shot
               *****************************************/               
              prepare_data_for_write(encodeMbAddr, &MB64, mb32count*4+mb16count, bestmode64);

              no_bits += write_one_macroblock64(1, mb32count, mb16count);

              if(MB64.mb_type != 0)
              {
                int mb_y = (mb32count>>1)*2 + (mb16count>>1);
                int mb_x = (mb32count&1) *2 + (mb16count&1);
                for(list=0; list<2; list++)            
                  for(b4_y=0; b4_y<4; b4_y++)
                    for(b4_x=0; b4_x<4; b4_x++)
                      for(k=0; k<2; k++)
                      {                    
                        img->mb_data[encodeMbAddr].mvd[list][b4_y][b4_x][k] = mb32_mvd[list][mb_y*4+b4_y][mb_x*4+b4_x][k];
                      }
              }

              mb16count ++;
            }
          }
          mb32count++;
        }
      }

      /*
       * propagate information to each 16x16 of 32x32 block
       */      
      propagate_mb_info(&MB64, CurrentMbAddr, CurrentMbAddr_x,CurrentMbAddr_y,extend_x, extend_y, 
        saved_cluster_qp, saved_prev_qp, saved_prev_delta_qp);
      /*
       * recover input->InterSearch
       */
      input->InterSearchSkip  = save_valid_mode[0];      
      input->InterSearch16x16 = save_valid_mode[1];
      input->InterSearch16x8  = save_valid_mode[2];
      input->InterSearch8x16  = save_valid_mode[3];

      img->currentSlice->num_mb += 16;
      (*NumberOfCodedMBs) += 16;

    }//mode 64x64, 64x32, 32x64
    else
    {
      int mb32count=0;
      
      //32x32 mode
      for(i=CurrentMbAddr_y; i<=extend_y; i+=2)
      {
        for(j=CurrentMbAddr_x; j<=extend_x; j+=2)
        {
          /****************
           * initialization
           ****************/
          encodeMbAddr = i*img->PicWidthInMbs + j;
          start_macroblock (encodeMbAddr, FALSE);

          assign_qp_info(mb32count,  encodeMbAddr, saved_cluster_qp, saved_prev_qp, saved_prev_delta_qp);//overwrite what have been done in start_macroblock for qp info

          if(i+1<=extend_y && j+1<=extend_x)
          {
// SPEEDUP
#ifdef RDO_Q
          if(input->UseRDO_Q)
            bestmode32 = saveBestMode.MBsize32[mb32count];
          else
#endif          
            bestmode32 = mb32type[mb32count];
          }
          /*************
           * coding
           *************/
          if(num_MB==16 && mb32count == 0)
            write_mb64_mbtype = 1;

          if(num_MB==16 && mb32count != 0)
            eos_bit_32 = 0;

          write_macroblock_cluster32(NumberOfCodedMBs, encodeMbAddr, mb32count, saved_prev_qp, saved_prev_delta_qp, saved_cluster_qp, delta_qp_sent);

          if(num_MB==16 && mb32count != 0)
            eos_bit_32 = 1;

          img->currentSlice->num_mb ++; //inside write_macroblock_cluster32, currSlice->num_mb--;

          write_mb64_mbtype = 0;

          mb32count++;
          //total_cost 
        }//j
      }//i
         

      /*************
       * reset mb qp
       *************/
      propagate_qp_info(*delta_qp_sent, CurrentMbAddr, CurrentMbAddr_x,CurrentMbAddr_y,extend_x, extend_y, 
        saved_cluster_qp, saved_prev_qp, saved_prev_delta_qp);//MB32_DELTA_QP
    }


    img->currentSlice->num_mb--;//currSlice->num_mb will be updated for the last mb of this cluster in terminate_macroblock  

    return nextMbAddr;
}

int write_macroblock_cluster(int *NumberOfCodedMBs, int CurrentMbAddr)
{
  int saved_prev_qp, saved_prev_delta_qp, saved_cluster_qp;
  if(img->type == I_SLICE)
  {
    write_one_macroblock(1);

    (*NumberOfCodedMBs)++;
    proceed2nextMacroblock();

    if(CurrentMbAddr+1 == img->PicSizeInMbs)
      return -1;
    else
      return (CurrentMbAddr+1);
  }
  else
  {
    int nextMbAddr;
//#ifdef MB32_DELTA_QP
    //qp handling
    start_macroblock(CurrentMbAddr, FALSE); // the purpose of calling this function is to set img->mb_data[CurrentMbAddr].(prev_qp, prev_delta_qp, qp)
    saved_prev_qp = img->mb_data[CurrentMbAddr].prev_qp;
    saved_prev_delta_qp = img->mb_data[CurrentMbAddr].prev_delta_qp;
    saved_cluster_qp = img->mb_data[CurrentMbAddr].qp;
    delta_qp_sent = 0;
//#endif

    if(input->UseExtMB == 2)
      nextMbAddr = write_macroblock_cluster64(NumberOfCodedMBs, CurrentMbAddr, saved_prev_qp, saved_prev_delta_qp, saved_cluster_qp, &delta_qp_sent);//MB64X64
    else
      nextMbAddr = write_macroblock_cluster32(NumberOfCodedMBs, CurrentMbAddr, 0, saved_prev_qp, saved_prev_delta_qp, saved_cluster_qp, &delta_qp_sent);

    return nextMbAddr;
  } 
}

/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode the reference
 *    parameter of a given MB.
 ****************************************************************************
 */
void writeRefFrame_CABAC32(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  MotionInfoContexts  *ctx    = img->currentSlice->mot_ctx;
  //Macroblock          *currMB = &img->mb_data[img->current_mb_nr];
  int                 addctx  = 0;

  int   a, b;
  int   act_ctx;
  int   act_sym;
  char** refframe_array = enc_picture->ref_idx[se->value2];

  int bslice = (img->type==B_SLICE);

  int   b8a, b8b;

  PixelPos block_a, block_b;

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
    mb_nr      = img->current_mb_nr + (img->PicWidthInMbs << (mb_ext_level-1));
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

  b8a=((block_a.x >> 1) & 0x01)+2*((block_a.y >> 1) & 0x01);
  b8b=((block_b.x >> 1) & 0x01)+2*((block_b.y >> 1) & 0x01);
  

  if(block_b.available && (mb_nr == img->current_mb_nr || mb_nr == img->current_mb_nr+(1<<(mb_ext_level-1))) )//block_b should be out of the current 32x32 block
  {
    block_b_is_direct = IS_DIRECT(&img->mb_data[block_b.mb_addr]);
    block_b_b8mode_is_direct = img->mb_data[block_b.mb_addr].b8mode[b8b]==0 && bslice;
  }
  if(block_a.available && (mb_nr == img->current_mb_nr || mb_nr == img->current_mb_nr + (img->PicWidthInMbs << (mb_ext_level-1))) )//block_a should be out of the current 32x32 block
  {
    block_a_is_direct = IS_DIRECT(&img->mb_data[block_a.mb_addr]);
    block_a_b8mode_is_direct = img->mb_data[block_a.mb_addr].b8mode[b8a]==0 && bslice;
  }

  if (!block_b.available)
    b=0;
  else if (( block_b_is_direct && !giRDOpt_B8OnlyFlag) || block_b_b8mode_is_direct )
    b=0;
  else
  {
#ifndef SIMPLIFY_CODE
    if (img->MbaffFrameFlag && (currMB->mb_field == 0) && (img->mb_data[block_b.mb_addr].mb_field == 1))
      b = (refframe_array[block_b.pos_y][block_b.pos_x] > 1 ? 1 : 0);
    else
#endif
      b = (refframe_array[block_b.pos_y][block_b.pos_x] > 0 ? 1 : 0);
  }

  if (!block_a.available)
    a=0;
  else if (( block_a_is_direct && !giRDOpt_B8OnlyFlag) || block_a_b8mode_is_direct)
    a=0;
  else 
  {
#ifndef SIMPLIFY_CODE
    if (img->MbaffFrameFlag && (currMB->mb_field == 0) && (img->mb_data[block_a.mb_addr].mb_field == 1))
      a = (refframe_array[block_a.pos_y][block_a.pos_x] > 1 ? 1 : 0);
    else
#endif
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
 *    This function is used to arithmetically encode the motion
 *    vector data of a B-frame MB.
 ****************************************************************************
 */
void writeMVD_CABAC32(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  //int i = img->subblock_x;
  //int j = img->subblock_y;
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
  //Macroblock          *currMB = &img->mb_data[img->current_mb_nr];

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
    mb_nr      = img->current_mb_nr + (img->PicWidthInMbs << (mb_ext_level-1));
  }
  else if(img->subblock_x==(2<<mb_ext_level) && img->subblock_y==0 )
  {
    subblock_x = 0;
    subblock_y = 0;
    mb_nr      = img->current_mb_nr + (1<<(mb_ext_level- 1));
  }
  else
  {
    error("subblock_x and subblock_y are out of the range\n", -1);
  }

  getLuma4x4Neighbour(mb_nr, subblock_x, subblock_y, -1,  0, &block_a);
  getLuma4x4Neighbour(mb_nr, subblock_x, subblock_y,  0, -1, &block_b);
  if (block_b.available)
  {
    b = absm(img->mb_data[block_b.mb_addr].mvd[list_idx][block_b.y][block_b.x][k]);
#ifndef SIMPLIFY_CODE
    if (img->MbaffFrameFlag && (k==1)) 
    {
      if ((currMB->mb_field==0) && (img->mb_data[block_b.mb_addr].mb_field==1))
        b *= 2;
      else if ((currMB->mb_field==1) && (img->mb_data[block_b.mb_addr].mb_field==0))
        b /= 2;
    }
#endif
  }
  else
    b=0;
          
  if (block_a.available)
  {
    a = absm(img->mb_data[block_a.mb_addr].mvd[list_idx][block_a.y][block_a.x][k]);
#ifndef SIMPLIFY_CODE
    if (img->MbaffFrameFlag && (k==1)) 
    {
      if ((currMB->mb_field==0) && (img->mb_data[block_a.mb_addr].mb_field==1))
        a *= 2;
      else if ((currMB->mb_field==1) && (img->mb_data[block_a.mb_addr].mb_field==0))
        a /= 2;
    }
#endif
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

void writeCBP_CABAC32(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  //TextureInfoContexts *ctx = img->currentSlice->tex_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  
  //int a, b;
  int curr_cbp_idx;
  int cbp = se->value1; // symbol to encode
  //int cbp_bit;
  //int b8;
  
  curr_cbp_idx = 1;
  writeCBP_BIT_CABAC (0, cbp, cbp, currMB, curr_cbp_idx, eep_dp);

#ifndef SIMPLIFY_CODE
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
#endif
}

void writeCBP_CABAC_luma1bit(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  //TextureInfoContexts *ctx = img->currentSlice->tex_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  
  //int a, b;
  int curr_cbp_idx;
  int cbp = se->value1; // symbol to encode
  //int cbp_bit;
  //int b8;
  
  curr_cbp_idx = 1;

  if((cbp&15)!=0)
    writeCBP_BIT_CABAC (0, 1, 1, currMB, curr_cbp_idx, eep_dp);
  else
    writeCBP_BIT_CABAC (0, 0, 0, currMB, curr_cbp_idx, eep_dp);


#ifndef SIMPLIFY_CODE
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
#endif
}

void SetMotionVectorPredictor32 (short  pmv[2],
                               char   **refPic,
                               short  ***tmp_mv,
                               short  ref_frame,
                               int    list,
                               int    block_x,
                               int    block_y,
                               int    blockshape_x,
                               int    blockshape_y,
                               int    mb_ext_level);

void writeMVD_CABAC32(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
/*!
************************************************************************
* \brief
*    Writes motion vectors of an 16x16 block
************************************************************************
*/
int write32MotionVector8x8 (int  i0,
                          int  j0,
                          int  i1,
                          int  j1,
                          int  refframe,
                          int  list_idx,
                          int  mv_mode,
                          int  mb_ext_level)
{
  int            i, j, k, l, m;
  int            curr_mvd;
  DataPartition* dataPart;
  
  int            rate       = 0;
  int            step_h     = input->part_size[mv_mode][0] << mb_ext_level;
  int            step_v     = input->part_size[mv_mode][1] << mb_ext_level;
  Macroblock*    currMB     = &img->mb_data[img->current_mb_nr];
  SyntaxElement* currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*         currSlice  = img->currentSlice;
  int*           bitCount   = currMB->bitcounter;
  const int*     partMap    = assignSE2partition[input->partition_mode];
  int            refindex   = refframe;
  
  short******    all_mv     = img->all_mv;
  short******    pred_mv    = img->pred_mv;
  
  if (currMB->bi_pred_me && currMB->b8pdir[0]==2 && mv_mode == 1 && refindex == 0)
    all_mv = currMB->bi_pred_me == 1? img->bipred_mv1 : img->bipred_mv2 ;  

  for (j=j0; j<j1; j+=step_v)
  {
    for (i=i0; i<i1; i+=step_h)
    {
#ifdef MB32X32_MVC
#ifdef RDO_Q
      if((input->Transform8x8Mode)&&(pass_with_writing == 1) && (input->mv_competition > 0) &&(!input->UseRDO_Q))
#endif
      {
        short best_prediction[2];
        SetMotionVectorPredictor_Competition32 (best_prediction, enc_picture->ref_idx[list_idx], enc_picture->mv[list_idx], (short)refindex, list_idx, i, j,(input->blc_size[mv_mode][0]<<mb_ext_level),(input->blc_size[mv_mode][1]<<mb_ext_level),mb_ext_level);
        GetBestPredVector(best_prediction, all_mv[j][i][list_idx][refindex][mv_mode][0],all_mv[j][i][list_idx][refindex][mv_mode][1],img->lambda_mf[img->type][img->qp],refindex,list_idx,mv_mode,i*4,j*4); 
        pred_mv[j][i][list_idx][refindex][mv_mode][0]=best_prediction[0];
        pred_mv[j][i][list_idx][refindex][mv_mode][1]=best_prediction[1];
      }
#endif
      
#ifdef RDO_Q
#ifdef ADAPTIVE_QUANTIZATION
      if(input->UseRDO_Q || img->slice_fractional_quant_flag)
#else
      if(input->UseRDO_Q)
#endif
      {
#ifdef MB32X32_MVC
        if( (input->mv_competition > 0))    
        {
          short best_prediction[2];
          SetMotionVectorPredictor_Competition32 (best_prediction, enc_picture->ref_idx[list_idx], enc_picture->mv[list_idx], (short)refindex, list_idx, i, j,(input->blc_size[mv_mode][0]<<mb_ext_level),(input->blc_size[mv_mode][1]<<mb_ext_level),mb_ext_level);
          GetBestPredVector(best_prediction, all_mv[j][i][list_idx][refindex][mv_mode][0],all_mv[j][i][list_idx][refindex][mv_mode][1],img->lambda_mf[img->type][img->qp],refindex,list_idx,mv_mode,i*4,j*4); 
          pred_mv[j][i][list_idx][refindex][mv_mode][0]=best_prediction[0];
          pred_mv[j][i][list_idx][refindex][mv_mode][1]=best_prediction[1];
        }
        else
        {
#endif
          short pred_mv_stored[2];
          SetMotionVectorPredictor32 (pred_mv_stored, enc_picture->ref_idx[list_idx], enc_picture->mv[list_idx], (short)refindex, list_idx, i, j, (input->blc_size[mv_mode][0]<<mb_ext_level),(input->blc_size[mv_mode][1]<<mb_ext_level), mb_ext_level);
          pred_mv[j][i][list_idx][refindex][mv_mode][0]=pred_mv_stored[0];
          pred_mv[j][i][list_idx][refindex][mv_mode][1]=pred_mv_stored[1];
#ifdef MB32X32_MVC
        }  
#endif
      }
#endif   //RDO_Q
      for (k=0; k<2; k++) 
      {        
        curr_mvd = all_mv[j][i][list_idx][refindex][mv_mode][k] - pred_mv[j][i][list_idx][refindex][mv_mode][k];
        
        //--- store (oversampled) mvd ---
        //only save mvd for the upper left 16x16 block for the context of the 2nd partition of 32x16 and 16x32        
        if(j==0 && i==0)
        {
          for (l=0; l < 4; l++)
          for (m=0; m < 4; m++)
          {
            currMB->mvd[list_idx][j+l][i+m][k] = curr_mvd; //mvd[2][BLOCK_MULTIPLE][BLOCK_MULTIPLE][2]
          }
        }
        //the mvd for the whole 32x32 is saved in mb32_mvd and will be copied to each 16x16 block's mvd in write_macroblock_cluster
        for (l=0; l < step_v; l++)
          for (m=0; m < step_h; m++)
          {
            mb32_mvd[list_idx][j+l][i+m][k] = curr_mvd;
          }
          currSE->value1 = curr_mvd;
          currSE->value2 = 0;
          currSE->type   = SE_MVD;
          if (input->symbol_mode == UVLC)
          {
            currSE->mapping = se_linfo;
          }
          else
          {
            img->subblock_x = i; // position used for context determination
            img->subblock_y = j; // position used for context determination
            currSE->value2  = 2*k+list_idx; // identifies the component and the direction; only used for context determination

            currSE->mb_ext_level  = mb_ext_level;
            currSE->writing = writeMVD_CABAC32;
          }  
          dataPart = &(currSlice->partArr[partMap[SE_MVD]]);
          dataPart->writeSyntaxElement (currSE, dataPart);

#if TRACE
          if (!list_idx)
          {
            snprintf(currSE->tracestring, TRACESTRING_SIZE, "mvd_l0 (%d) = %3d  (org_mv %3d pred_mv %3d)",k, curr_mvd, all_mv[j][i][list_idx][refindex][mv_mode][k], pred_mv[j][i][list_idx][refindex][mv_mode][k]);
          }
          else
          {
            snprintf(currSE->tracestring, TRACESTRING_SIZE, "mvd_l1 (%d) = %3d  (org_mv %3d pred_mv %3d)",k, curr_mvd, all_mv[j][i][list_idx][refindex][mv_mode][k], pred_mv[j][i][list_idx][refindex][mv_mode][k]);
          }
          
#endif
          bitCount[BITS_INTER_MB] += currSE->len;
          rate                    += currSE->len;
          currSE++;  
          currMB->currSEnr++;
      }
#ifdef MB32X32_MVC
      if (input->mv_competition > 0)
      { 
        int maxmode=0;
        if (img->type == P_SLICE)
          maxmode = mv_comp.nb_mode_for_mvp;
        else if (img->type == B_SLICE)
          maxmode = mv_comp.nb_mode_for_mvb;
        //Rahul ----- reverse logic....0 means "send it", 1 means "dont send"
        if (send_index_mv_prediction[mv_mode][refindex][list_idx][j][i][maxmode] == 0)
          rate+=writeIndexForMotionVectorPrediction (refframe, list_idx, mv_mode, i, j);
      }
#endif
    }
  }
  return rate;
}


/*!
************************************************************************
* \brief
*    Codes the reference frame
************************************************************************
*/
int writeReferenceFrame32 (int mode, int i, int j, int fwd_flag, int  ref, int  mb_ext_level)
{
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  int*            bitCount  = currMB->bitcounter;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int             rate      = 0;
  DataPartition*  dataPart;
  int             num_ref   = ( fwd_flag ? listXsize[LIST_0+currMB->list_offset]: listXsize[LIST_1+currMB->list_offset]);
  int             flag_mode = 0;
  
  if( num_ref == 1 )
  {
    return 0;
  }
  
  if ( num_ref == 2 )
  {
    flag_mode = 1;
  }
  
  currSE->value1 = ref;
  currSE->value2  = 0;
  currSE->type   = SE_REFFRAME;
  
  dataPart = &(currSlice->partArr[partMap[currSE->type]]);
  if (input->symbol_mode == UVLC)
  {
    if( flag_mode )
    {
      currSE->bitpattern = 1 - currSE->value1;
      currSE->len = 1;
      writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
    }
    else
    {
      currSE->mapping = ue_linfo;
      dataPart->writeSyntaxElement (currSE, dataPart);
    }
  }
  else
  {
    currSE->context = BType2CtxRef (mode);
    img->subblock_x = i; // position used for context determination
    img->subblock_y = j; // position used for context determination
    currSE->mb_ext_level = mb_ext_level;
    currSE->writing = writeRefFrame_CABAC32;
    currSE->value2 = (fwd_flag)? LIST_0:LIST_1;
    dataPart->writeSyntaxElement (currSE, dataPart);
  }
  bitCount[BITS_INTER_MB] += currSE->len;
  rate                    += currSE->len;
#if TRACE
  if (fwd_flag)
  {
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "ref_idx_l0 = %d", currSE->value1);
  }
  else
  {
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "ref_idx_l1 = %d", currSE->value1);
  }
#endif
  currSE++;
  currMB->currSEnr++;
  
  return rate;
}


/*!
************************************************************************
* \brief
*    Writes motion info
************************************************************************
*/
int writeMotionInfo2NAL32 (int mb_ext_level)
{
  int k, j0, i0, refframe;
  int jj;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  int             no_bits   = 0;
  int   bframe          = (img->type==B_SLICE);
  int   step_h0         = ((input->blc_size[IS_P8x8(currMB) ? 4 : currMB->mb_type][0]<<mb_ext_level) >> 2);
  int   step_v0         = ((input->blc_size[IS_P8x8(currMB) ? 4 : currMB->mb_type][1]<<mb_ext_level) >> 2);

  //=== If multiple ref. frames, write reference frame for the MB ===
  if (IS_INTERMV (currMB))
  {
    // if UVLC is turned on, a 8x8 macroblock with all ref=0 in a P-frame is signalled in macroblock mode
    if (!IS_P8x8 (currMB) || !ZeroRef (currMB) || input->symbol_mode==CABAC || bframe)
    {
      for (j0=0; j0<(4<<mb_ext_level); j0+=step_v0)
      {
        jj = img->block_y+j0;
        for (i0=0; i0<(4<<mb_ext_level); i0+=step_h0)
        {
          k=(j0>>mb_ext_level)+(i0 >> (1+mb_ext_level));
          
          if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)//has forward vector
          {
            no_bits += writeReferenceFrame32 (currMB->b8mode[k], i0, j0, 1, enc_picture->ref_idx[LIST_0][jj][img->block_x+i0], mb_ext_level);
          }
        }
      }

      for (j0=0; j0<(4<<mb_ext_level); j0+=step_v0)
      {
        jj = img->block_y+j0;
        for (i0=0; i0<(4<<mb_ext_level); i0+=step_h0)
        {
          k=(j0>>mb_ext_level)+(i0 >> (1+mb_ext_level));
          if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)//has backward vector
          {
            no_bits += writeReferenceFrame32 (currMB->b8mode[k], i0, j0, 0, enc_picture->ref_idx[LIST_1][jj][img->block_x+i0], mb_ext_level);
          }
        }
      }
    }
  }
  
  //===== write forward motion vectors =====
  if (IS_INTERMV (currMB))
  {
    for (j0=0; j0<(4<<mb_ext_level); j0+=step_v0)
    {
      jj = img->block_y+j0;
      for (i0=0; i0<(4<<mb_ext_level); i0+=step_h0)
      {
        k=(j0>>mb_ext_level)+(i0 >> (1+mb_ext_level));
        if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)//has forward vector
        {
          refframe  = enc_picture->ref_idx[LIST_0][jj][img->block_x+i0];
          no_bits  += write32MotionVector8x8 (i0, j0, i0+step_h0, j0+step_v0, refframe, LIST_0, currMB->b8mode[k], mb_ext_level);
        }
      }
    }
  }
  
  
  //===== write backward motion vectors =====
  if (IS_INTERMV (currMB) && bframe)
  {
    for (j0=0; j0<(4<<mb_ext_level); j0+=step_v0)
    {
      jj = img->block_y+j0;
      for (i0=0; i0<(4<<mb_ext_level); i0+=step_h0)
      {
        k=(j0>>mb_ext_level)+(i0 >> (1+mb_ext_level));
        if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)//has backward vector
        {
          refframe  = enc_picture->ref_idx[LIST_1][jj][img->block_x+i0];
          no_bits  += write32MotionVector8x8 (i0, j0, i0+step_h0, j0+step_v0, refframe, LIST_1, currMB->b8mode[k], mb_ext_level);
        }
      }
    }
  }
  return no_bits;
}
// writeMB_typeInfo_CABAC32 is based on writeMB_typeInfo_CABAC
// the first bit of the P slice mb type is omitted.
void writeMB_typeInfo_CABAC32(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
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
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][5]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][6]);
        break;
      case 2:
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][5]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][7]);
        break;
      case 3:
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][5]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][7]);
        break;
      case 4:
      case 5:
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

int writeCBP32(int rdopt, int cbp, int dquant)
{
  //int             mb_x, mb_y, i, j, k;
  //int             level, run;
  int             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  int*            bitCount  = currMB->bitcounter;
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  //int             cbp       = currMB->cbp;
  DataPartition*  dataPart;
  //int             need_transform_size_flag;   //ADD-VG-24062004
  
//  int   b8, b4;
//  int*  DCLevel = img->cofDC[0][0];
//  int*  DCRun   = img->cofDC[0][1];
//  int*  ACLevel;
//  int*  ACRun;
  
  if (1)//(!IS_NEWINTRA (&MB32))
  {
    //=====   C B P   =====
    //---------------------
    currSE->value1 = cbp;

#ifndef SIMPLIFY_CODE
    if (IS_OLDINTRA (currMB) || currMB->mb_type == SI4MB ||  currMB->mb_type == I8MB)
    {
      if (input->symbol_mode == UVLC)  currSE->mapping = cbp_linfo_intra;
      currSE->type = SE_CBP_INTRA;
    }
    else
#endif 
    {
      if (input->symbol_mode == UVLC)  currSE->mapping = cbp_linfo_inter;
      currSE->type = SE_CBP_INTER;
    }

    if (input->symbol_mode == CABAC)   currSE->writing = writeCBP_CABAC32;
    
    // choose the appropriate data partition
    dataPart = &(currSlice->partArr[partMap[currSE->type]]);
    
    dataPart->writeSyntaxElement(currSE, dataPart);
    bitCount[BITS_CBP_MB] += currSE->len;
    rate                  += currSE->len;
#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "CBP (%2d,%2d) = %3d",img->mb_x, img->mb_y, cbp);
#endif

    // proceed to next SE
    currSE++;  
    currMB->currSEnr++;        
  }
  
//#ifdef MB32_DELTA_QP
  if(cbp != 0 && dquant)
  {
        currSE->value1 = currMB->delta_qp;
      
      if (input->symbol_mode==UVLC)   currSE->mapping = se_linfo;
      else                            currSE->writing = writeDquant_CABAC;
      
      if (IS_INTER (currMB))  currSE->type = SE_DELTA_QUANT_INTER;
      else                    currSE->type = SE_DELTA_QUANT_INTRA;
      
      
      // choose the appropriate data partition
      dataPart = &(img->currentSlice->partArr[partMap[currSE->type]]);
      dataPart->writeSyntaxElement(  currSE, dataPart);
      bitCount[BITS_DELTA_QUANT_MB] += currSE->len;
      rate                          += currSE->len;
  #if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "Delta QP (%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->delta_qp);
  #endif



      // proceed to next SE
      currSE++;
      currMB->currSEnr++;
  }
//#endif

  return rate;
}


/* mb32: (0~3) indext of 32x32 block
 * mb16: (0~3) index of 16x16 block inside a 32x32 block
 */
int writeMBLayer64 (int rdopt, int *coeff_rate, int mb32i, int mb16i)//MB64X64
{
  int             i,j;
  int             mb_nr      = img->current_mb_nr;
  Macroblock*     currMB     = &img->mb_data[mb_nr];
  SyntaxElement  *currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
  int*            bitCount   = currMB->bitcounter;
  Slice*          currSlice  = img->currentSlice;
  DataPartition*  dataPart;
  const int*      partMap    = assignSE2partition[input->partition_mode];
  int             no_bits    = 0;
  int             mb_type;  
  //int             mb_field_tmp;
  //Macroblock      *topMB = NULL;
  
  int             WriteFrameFieldMBInHeader = 0;
#ifndef SIMPLIFY_CODE
  if (img->MbaffFrameFlag)
  {
    if (0==(mb_nr & 0x01))
    {
      WriteFrameFieldMBInHeader = 1; // top field
      
      prevMbSkipped = 0;
    }
    else
    {
      if (prevMB->mb_type ? 0:((img->type == B_SLICE) ? !prevMB->cbp:1))
      {
        WriteFrameFieldMBInHeader = 1; // bottom, if top was skipped
      }
      
      topMB= &img->mb_data[prev_mb_nr];
      prevMbSkipped = topMB->skip_flag;
    }
  }
#endif
  currMB->IntraChromaPredModeFlag = IS_INTRA(currMB);
  
  // choose the appropriate data partition
  dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);
  
  if(img->type == I_SLICE)
  {
    //========= write mb_aff (I_SLICE) =========
    if(WriteFrameFieldMBInHeader)
    {
      currSE->value1 = currMB->mb_field;
      currSE->value2 = 0;
      currSE->type   = SE_MBTYPE;      
      
#if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "mb_field_decoding_flag");
#endif
      if( input->symbol_mode==UVLC)
      {
        currSE->mapping = ue_linfo;
        currSE->bitpattern = (currMB->mb_field ? 1 : 0);
        currSE->len = 1;
        writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
      }
      else
      {
        currSE->writing = writeFieldModeInfo_CABAC;
        dataPart->writeSyntaxElement(currSE, dataPart);
      }
      
      bitCount[BITS_MB_MODE] += currSE->len;
      no_bits                += currSE->len;
      currSE++;
      currMB->currSEnr++;
    }
    
    //========= write mb_type (I_SLICE) =========
    currSE->value1  = MBType2Value (currMB);
    currSE->value2  = 0;
    currSE->type    = SE_MBTYPE;
    
    if (input->symbol_mode == UVLC)  
      currSE->mapping = ue_linfo;
    else
      currSE->writing = writeMB_typeInfo_CABAC;
    
    dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE,   "mb_type (I_SLICE) (%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
#endif
    bitCount[BITS_MB_MODE] += currSE->len;
    no_bits                += currSE->len;
    currSE++;
    currMB->currSEnr++;
    
  }
  // not I_SLICE, CABAC
  else if (input->symbol_mode == CABAC)
  {
#ifndef SIMPLIFY_CODE
    if (img->MbaffFrameFlag && ((img->current_mb_nr & 0x01) == 0||prevMbSkipped))
    {
      mb_field_tmp = currMB->mb_field;
      currMB->mb_field = field_flag_inference();
      CheckAvailabilityOfNeighborsCABAC();
      currMB->mb_field = mb_field_tmp;
    }
#endif

    /********************************************
     *MB type for mode 32x32, 32x16, 16x32, 16x16
     *******************************************/
    if(mb32i==0 && mb16i==0)
    {

      //========= write mb_skip_flag (CABAC) =========
      mb_type         = MBType2Value (&MB64); //MB32 set in write_macroblock_cluster

      currSE->value1  = mb_type;

      currSE->value2  = MB64.cbp; //MB32.cbp set in write_macroblock_cluster

      currSE->type    = SE_MBTYPE;
      currSE->writing = writeMB_skip_flagInfo_CABAC;
      dataPart->writeSyntaxElement( currSE, dataPart);
  #if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "mb_skip_flag");
  #endif
      bitCount[BITS_MB_MODE] += currSE->len;
      no_bits                += currSE->len;
      currSE++;
      currMB->currSEnr++;


      CheckAvailabilityOfNeighborsCABAC();

      //========= write mb_type (CABAC) =========
      if (MB64.mb_type != 0 || ((img->type == B_SLICE) && MB64.cbp != 0))
      {
        currSE->value1  = mb_type;
        currSE->value2  = 0;
        currSE->type    = SE_MBTYPE;
        currSE->writing = writeMB_typeInfo_CABAC32;
        dataPart->writeSyntaxElement( currSE, dataPart);
  #if TRACE
        if (img->type == B_SLICE) 
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "mb_type (B_SLICE) (%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
        else                      
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "mb_type (P_SLICE) (%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
  #endif
        bitCount[BITS_MB_MODE] += currSE->len;
        no_bits                += currSE->len;
        currSE++;
        currMB->currSEnr++;
        
      }
    }


    CheckAvailabilityOfNeighborsCABAC();
  }  
  
#ifdef MB32X32_MVC
  if (input->mv_competition > 0)
  {
    if((rdopt==0) && (input->symbol_mode == CABAC) && (img->type==P_SLICE) && (mv_comp.nb_mode_for_skip > 1))  
    {  
      //make sure that this function is called only once...whether it is big blocks or 16x16
      //64x64 P-Skip
      if(mb32i==0 && mb16i==0 && MB64.mb_type==0)
          no_bits  += write_predictor_index_for_skip_mode64();
    }
  }
#endif

  
  //init NoMbPartLessThan8x8Flag
  currMB->NoMbPartLessThan8x8Flag = (IS_DIRECT(currMB) && !(active_sps->direct_8x8_inference_flag))? 0: 1;
  
  if (currMB->mb_type == IPCM)
  {
    int jj, uv;
    if (input->symbol_mode == CABAC) 
    {
      int len;
      EncodingEnvironmentPtr eep = &dataPart->ee_cabac;
      len = arienco_bits_written(eep);
      arienco_done_encoding(eep); // This pads to byte
      len = arienco_bits_written(eep) - len;
      currSE--; // Need to append more bits to mb_type
      currSE->len += len;
      no_bits += len;
      currSE++;
      // Now restart the encoder
      arienco_start_encoding(eep, dataPart->bitstream->streamBuffer, &(dataPart->bitstream->byte_pos));
      reset_pic_bin_count();
    }
    if (dataPart->bitstream->bits_to_go < 8)
    {
      // This will only happen in the CAVLC case, CABAC is already padded
      currSE->mapping = ue_linfo;
      currSE->len = dataPart->bitstream->bits_to_go;
      no_bits += currSE->len;
      bitCount[BITS_COEFF_Y_MB]+= currSE->len;
      currSE->bitpattern = 0;
#if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "IPCM alignment bits = %d", currSE->len);
#endif
      writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
      currSE++;
      currMB->currSEnr++;
    }
    for (j=0;j<MB_BLOCK_SIZE;j++)
    {
      jj = img->pix_y+j;
      for (i=0;i<MB_BLOCK_SIZE;i++)
      {
        currSE->mapping = ue_linfo;
        currSE->len = img->bitdepth_luma;  
        currSE->type    = SE_MBTYPE;
        no_bits += currSE->len;
        currSE->bitpattern = enc_picture->imgY[jj][img->pix_x+i];       
        currSE->value1 = currSE->bitpattern;
        currSE->type = SE_MBTYPE;
        bitCount[BITS_COEFF_Y_MB]+=currSE->len;
#if TRACE
        snprintf(currSE->tracestring, TRACESTRING_SIZE, "IPCM Luma (%d %d) = %d", j,i,currSE->bitpattern);
#endif
        writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);        
        currSE++;
        currMB->currSEnr++;
      }
    }
    if (img->yuv_format != YUV400)
    {
      for (uv = 0; uv < 2; uv ++)
      {
        for (j=0;j<img->mb_cr_size_y;j++)
        {
          jj = img->pix_c_y+j;
          for (i=0;i<img->mb_cr_size_x;i++)
          {
            currSE->mapping = ue_linfo;
            currSE->len = img->bitdepth_chroma;
            currSE->type    = SE_MBTYPE;
            no_bits += currSE->len;
            currSE->bitpattern = enc_picture->imgUV[uv][jj][img->pix_c_x+i];
            currSE->value1 = currSE->bitpattern;
            currSE->type = SE_MBTYPE;
            bitCount[BITS_COEFF_UV_MB]+=currSE->len;
#if TRACE
            snprintf(currSE->tracestring, TRACESTRING_SIZE, "IPCM chroma(%d) (%d %d) = %d", uv, j,i,currSE->bitpattern);
#endif
            writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
            currSE++;
            currMB->currSEnr++;
          }
        }
      }
    }
    return no_bits;
  }
  
  //===== BITS FOR 8x8 SUB-PARTITION MODES =====
  if (IS_P8x8 (currMB))
  {
    dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);
    
    for (i=0; i<4; i++)
    {
      if (input->symbol_mode==UVLC)   
        currSE->mapping = ue_linfo;
      else
        currSE->writing = writeB8_typeInfo_CABAC;
      
      currSE->value1  = B8Mode2Value (currMB->b8mode[i], currMB->b8pdir[i]);
      currSE->value2  = 0;
      currSE->type    = SE_MBTYPE;
      dataPart->writeSyntaxElement (currSE, dataPart);
#if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "8x8 mode/pdir(%2d) = %3d/%d", i, currMB->b8mode[i], currMB->b8pdir[i]);
#endif
      bitCount[BITS_MB_MODE]+= currSE->len;
      no_bits               += currSE->len;
      currSE++;
      currMB->currSEnr++;
      
      //set NoMbPartLessThan8x8Flag for P8x8 mode
      currMB->NoMbPartLessThan8x8Flag &= (currMB->b8mode[i]==0 && active_sps->direct_8x8_inference_flag) || 
        (currMB->b8mode[i]==4);
    }
#ifdef MB32X32_MVC
    if (input->mv_competition > 0)
      pass_with_writing = !rdopt;
#endif
    no_bits += writeMotionInfo2NAL  ();
    currSE   = &img->MB_SyntaxElements[currMB->currSEnr];
  }
  
  //============= Transform size flag for INTRA MBs =============
  //-------------------------------------------------------------
  //transform size flag for INTRA_4x4 and INTRA_8x8 modes
  if ((currMB->mb_type == I8MB || currMB->mb_type == I4MB) && input->Transform8x8Mode)
  {
    currSE->value1 = currMB->luma_transform_size_8x8_flag;
    currSE->type   = SE_HEADER;
    
    if( input->symbol_mode==UVLC)
    {
      currSE->mapping = ue_linfo;
      currSE->type    = SE_MBTYPE;
      currSE->bitpattern = currMB->luma_transform_size_8x8_flag;
      currSE->len = 1;
      writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
    }
    else
    {
      currSE->writing = writeMB_transform_size_CABAC;
      dataPart->writeSyntaxElement(currSE, dataPart);
    }
#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "transform_size_8x8_flag = %3d", currMB->luma_transform_size_8x8_flag);
#endif
    
    bitCount[BITS_MB_MODE] += currSE->len;
    no_bits                += currSE->len;
    currSE++;
    currMB->currSEnr++;
  }
  
  
  //===== BITS FOR INTRA PREDICTION MODES ====
  no_bits += writeIntraModes();
  //===== BITS FOR CHROMA INTRA PREDICTION MODE ====
  if (currMB->IntraChromaPredModeFlag && img->yuv_format != YUV400)
    no_bits += writeChromaIntraPredMode();
  else if(!rdopt) //GB CHROMA !!!!!
    currMB->c_ipred_mode = DC_PRED_8; //setting c_ipred_mode to default is not the right place here
  //resetting in rdopt.c (but where ??)
  //with cabac and bframes maybe it could crash without this default
  //since cabac needs the right neighborhood for the later MBs
  
  //----- motion information -----
  if(mb32i==0 && mb16i==0)
  {

    if(MB64.mb_type != 0)
    {
#ifdef MB32X32_MVC
      if (input->mv_competition > 0)
        pass_with_writing = !rdopt;
#endif
      no_bits  += writeMotionInfo2NAL32  (2);
    }
  }


  if(MB64.mb_type!=0 || (img->type==B_SLICE && MB64.cbp!=0))
  {
    if(mb32i==0 && mb16i==0)
    {

      no_bits += writeCBP32(rdopt, MB64.cbp, 1);//with dquant
    }

    if(MB64.cbp!=0)
    {
      if(mb16i == 0)
      {

        no_bits += writeCBP32(rdopt, MB32.cbp, 0);//with dquant
      }

      if(MB32.cbp != 0)
      {
        *coeff_rate = writeCBPandLumaCoeff (rdopt, 1, mb16i);

        if (img->yuv_format != YUV400)
          *coeff_rate  += writeChromaCoeff ();

        no_bits  += *coeff_rate;           
      }   
    }
  }
#ifdef ADAPTIVE_QUANTIZATION
  if(!rdopt)
  {
    dwSkipMbCount+=(currMB->cbp==0);
    dwSkipEnableMbCount++;
    
    dwSkipEnableMbCountForSlice[(img->type)+(img->type==B_SLICE? (img->nal_reference_idc<<2) : 0)]++;
    dwSkipMbCountForSlice[(img->type)+(img->type==B_SLICE? (img->nal_reference_idc<<2) : 0)]+=(currMB->cbp==0);
  }
#endif
#ifdef ADAPTIVE_FD_SD_CODING
  else
  {
    if (cabac_encoding)
    {
      currMB->SD_or_FD_t8x8=0;
    }
  }  
#endif
  return no_bits;
}

/*!
************************************************************************
* \brief
*    Passes the chosen syntax elements to the NAL for 16x16 block inside of 64x64 block
* \param coeff_rate
*    bitrate of Luma and Chroma coeff
* \param (mb32i,mb16i)
*    mb32i: index of 32x32 block inside 64x64
*    mb16i: index of 16x16 block inside 32x32
************************************************************************
*/
int write_one_macroblock64 (int eos_bit, int mb32i, int mb16i)
{
  Macroblock* currMB   = &img->mb_data[img->current_mb_nr];
  int*        bitCount = currMB->bitcounter;
  int i;
  int no_bits;
  
  extern int cabac_encoding;



#ifndef SIMPLIFY_CODE
  //--- constrain intra prediction ---
  if(input->UseConstrainedIntraPred && (img->type==P_SLICE || img->type==B_SLICE))
  {
    if( !IS_INTRA(currMB) )
    {
      img->intra_block[img->current_mb_nr] = 0;
    }
    else
    {
      img->intra_block[img->current_mb_nr] = 1;
    }
  }
#endif
  //===== init and update number of intra macroblocks =====
  if (img->current_mb_nr==0)
    intras=0;
  
  if (IS_INTRA(currMB))
    intras++;
  
  //--- write non-slice termination symbol if the macroblock is not the first one in its slice ---
  if (input->symbol_mode==CABAC && img->current_mb_nr!=img->currentSlice->start_mb_nr && eos_bit && mb32i==0 && mb16i==0)
  {
    write_terminating_bit (0);
  }
  
  cabac_encoding = 1;
  
#ifdef ADAPTIVE_FD_SD_CODING
  if (cabac_encoding) currMB->written_SD_Coding_on_off=0;
#endif
#ifdef ADAPTIVE_QUANTIZATION
  // Condieration for no residual
  if ( img->slice_fractional_quant_flag && currMB->mb_type==IPCM)
    currMB->mb_iaqms_idx = 0;
  if(img->slice_fractional_quant_flag && (currMB->cbp==0 && !IS_NEWINTRA (currMB)) )
    currMB->mb_iaqms_idx = 0;
#endif
  
  //--- write macroblock ---
  no_bits = writeMBLayer64 (0, &i, mb32i, mb16i);  // i is temporary

  if (!((currMB->mb_type !=0 ) || ((img->type==B_SLICE) && currMB->cbp != 0) ))
  { 
    for (i=0; i < 4; i++)
      memset(img->nz_coeff [img->current_mb_nr][i], 0, (4 + img->num_blk8x8_uv) * sizeof(int));  // CAVLC
  }
    
  //--- set total bit-counter ---
  
#ifdef MB32X32_MVC
  if (input->mv_competition > 0)
  {
    bitCount[BITS_TOTAL_MB] = bitCount[BITS_MB_MODE]  + bitCount[BITS_COEFF_Y_MB]     
      + bitCount[BITS_INTER_MB] + bitCount[BITS_CBP_MB]  
      + bitCount[BITS_DELTA_QUANT_MB] + bitCount[BITS_COEFF_UV_MB]
      + bitCount[BITS_SKIP_PRED]
      + bitCount[BITS_MV_PRED];
  }
  else
#endif
  {
    bitCount[BITS_TOTAL_MB] = bitCount[BITS_MB_MODE]  + bitCount[BITS_COEFF_Y_MB]     
      + bitCount[BITS_INTER_MB] + bitCount[BITS_CBP_MB]  
      + bitCount[BITS_DELTA_QUANT_MB] + bitCount[BITS_COEFF_UV_MB];
    
    //Rate control
    img->NumberofMBHeaderBits=bitCount[BITS_MB_MODE]   + bitCount[BITS_INTER_MB]
      + bitCount[BITS_CBP_MB]  + bitCount[BITS_DELTA_QUANT_MB];
  }  
  ;
  img->NumberofMBTextureBits= bitCount[BITS_COEFF_Y_MB]+ bitCount[BITS_COEFF_UV_MB];
  img->NumberofTextureBits +=img->NumberofMBTextureBits;
  img->NumberofHeaderBits +=img->NumberofMBHeaderBits;
  /*basic unit layer rate control*/
  if(img->BasicUnit<img->Frame_Total_Number_MB)
  {
    img->NumberofBasicUnitHeaderBits +=img->NumberofMBHeaderBits;
    img->NumberofBasicUnitTextureBits +=img->NumberofMBTextureBits;
  }
  /*record the total number of MBs*/
  img->NumberofCodedMacroBlocks++;
  
  stats->bit_slice += bitCount[BITS_TOTAL_MB];
  
  cabac_encoding = 0;



  return no_bits;
}


/*!
************************************************************************
* \brief
*    Codes macroblock header
* \param rdopt
*    true for calls during RD-optimization
* \param coeff_rate
*    bitrate of Luma and Chroma coeff
* \param (mb32,mbi)
*        (1,0) - 32x32, 32x16, 16x32 mode info, mv info, coeff of upper left 16x16, dequant
*        (1,1~)-                                         coeff of other      16x16
*        (0,0) - 16x16 mode, upper left 16x16 partition mode, mv, coeff, dequant
*        (0,1~)-                        16x16 partition mode, mv, coeff
************************************************************************
*/
int writeMBLayer32 (int rdopt, int *coeff_rate, int mb32, int mbi)
{
  int             i,j;
  int             mb_nr      = img->current_mb_nr;
  Macroblock*     currMB     = &img->mb_data[mb_nr];
  SyntaxElement  *currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
  int*            bitCount   = currMB->bitcounter;
  Slice*          currSlice  = img->currentSlice;
  DataPartition*  dataPart;
  const int*      partMap    = assignSE2partition[input->partition_mode];
  int             no_bits    = 0;
  int             mb_type=0;    
  int             WriteFrameFieldMBInHeader = 0;
#ifndef SIMPLIFY_CODE
  if (img->MbaffFrameFlag)
  {
    if (0==(mb_nr & 0x01))
    {
      WriteFrameFieldMBInHeader = 1; // top field
      
      prevMbSkipped = 0;
    }
    else
    {
      if (prevMB->mb_type ? 0:((img->type == B_SLICE) ? !prevMB->cbp:1))
      {
        WriteFrameFieldMBInHeader = 1; // bottom, if top was skipped
      }
      
      topMB= &img->mb_data[prev_mb_nr];
      prevMbSkipped = topMB->skip_flag;
    }
  }
#endif
  currMB->IntraChromaPredModeFlag = IS_INTRA(currMB);
  
  // choose the appropriate data partition
  dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);
  
  if(img->type == I_SLICE)
  {
    //========= write mb_aff (I_SLICE) =========
    if(WriteFrameFieldMBInHeader)
    {
      currSE->value1 = currMB->mb_field;
      currSE->value2 = 0;
      currSE->type   = SE_MBTYPE;      
      
#if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "mb_field_decoding_flag");
#endif
      if( input->symbol_mode==UVLC)
      {
        currSE->mapping = ue_linfo;
        currSE->bitpattern = (currMB->mb_field ? 1 : 0);
        currSE->len = 1;
        writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
      }
      else
      {
        currSE->writing = writeFieldModeInfo_CABAC;
        dataPart->writeSyntaxElement(currSE, dataPart);
      }
      
      bitCount[BITS_MB_MODE] += currSE->len;
      no_bits                += currSE->len;
      currSE++;
      currMB->currSEnr++;
    }
    
    //========= write mb_type (I_SLICE) =========
    currSE->value1  = MBType2Value (currMB);
    currSE->value2  = 0;
    currSE->type    = SE_MBTYPE;
    
    if (input->symbol_mode == UVLC)  
      currSE->mapping = ue_linfo;
    else
      currSE->writing = writeMB_typeInfo_CABAC;
    
    dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE,   "mb_type (I_SLICE) (%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
#endif
    bitCount[BITS_MB_MODE] += currSE->len;
    no_bits                += currSE->len;
    currSE++;
    currMB->currSEnr++;
    
  }
  // not I_SLICE, CABAC
  else if (input->symbol_mode == CABAC)
  {
#ifndef SIMPLIFY_CODE
    if (img->MbaffFrameFlag && ((img->current_mb_nr & 0x01) == 0||prevMbSkipped))
    {
      mb_field_tmp = currMB->mb_field;
      currMB->mb_field = field_flag_inference();
      CheckAvailabilityOfNeighborsCABAC();
      currMB->mb_field = mb_field_tmp;
    }
#endif


    if(input->UseExtMB == 2)//MB64X64
    {
      if(write_mb64_mbtype && mbi==0)
      {
        Macroblock mb32x32;
        mb32x32.mb_type = P8x8;

        //========= write mb_skip_flag (CABAC) =========
        mb_type         = MBType2Value (&mb32x32);
        currSE->value1  = mb_type;
        currSE->value2  = currMB->cbp;
        currSE->type    = SE_MBTYPE;
        currSE->writing = writeMB_skip_flagInfo_CABAC;
        dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
        snprintf(currSE->tracestring, TRACESTRING_SIZE, "mb_skip_flag");
#endif
        bitCount[BITS_MB_MODE] += currSE->len;
        no_bits                += currSE->len;
        currSE++;
        currMB->currSEnr++;
        
        CheckAvailabilityOfNeighborsCABAC();
            
        //========= write mb_type (CABAC) =========
        if(1)
        {
          currSE->value1  = mb_type;
          currSE->value2  = 0;
          currSE->type    = SE_MBTYPE;
          currSE->writing = writeMB_typeInfo_CABAC32;
          dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
          if (img->type == B_SLICE) 
            snprintf(currSE->tracestring, TRACESTRING_SIZE, "mb_type (B_SLICE) (%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
          else                      
            snprintf(currSE->tracestring, TRACESTRING_SIZE, "mb_type (P_SLICE) (%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
#endif
          bitCount[BITS_MB_MODE] += currSE->len;
          no_bits                += currSE->len;
          currSE++;
          currMB->currSEnr++;
          
        }
      }
    }

    /********************************************
     *MB type for mode 32x32, 32x16, 16x32, 16x16
     *******************************************/
    if(mbi==0)
    {
      //========= write mb_skip_flag (CABAC) =========
      mb_type         = MBType2Value (&MB32); //MB32 set in write_macroblock_cluster

      currSE->value1  = mb_type;

      currSE->value2  = MB32.cbp; //MB32.cbp set in write_macroblock_cluster

      currSE->type    = SE_MBTYPE;
      currSE->writing = writeMB_skip_flagInfo_CABAC;
      dataPart->writeSyntaxElement( currSE, dataPart);
  #if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "mb_skip_flag");
  #endif
      bitCount[BITS_MB_MODE] += currSE->len;
      no_bits                += currSE->len;
      currSE++;
      currMB->currSEnr++;


      CheckAvailabilityOfNeighborsCABAC();

      //========= write mb_type (CABAC) =========
      if (MB32.mb_type != 0 || ((img->type == B_SLICE) && MB32.cbp != 0))
      {
        currSE->value1  = mb_type;
        currSE->value2  = 0;
        currSE->type    = SE_MBTYPE;
        currSE->writing = writeMB_typeInfo_CABAC32;
        dataPart->writeSyntaxElement( currSE, dataPart);
  #if TRACE
        if (img->type == B_SLICE) 
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "mb_type (B_SLICE) (%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
        else                      
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "mb_type (P_SLICE) (%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
  #endif
        bitCount[BITS_MB_MODE] += currSE->len;
        no_bits                += currSE->len;
        currSE++;
        currMB->currSEnr++;
        
      }
    }


    if(mb32==0)
    {
      //========= write mb_skip_flag (CABAC) =========
      mb_type         = MBType2Value (currMB);
      currSE->value1  = mb_type;
      currSE->value2  = currMB->cbp;
      currSE->type    = SE_MBTYPE;
      currSE->writing = writeMB_skip_flagInfo_CABAC;
      dataPart->writeSyntaxElement( currSE, dataPart);
  #if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "mb_skip_flag");
  #endif
      bitCount[BITS_MB_MODE] += currSE->len;
      no_bits                += currSE->len;
      currSE++;
      currMB->currSEnr++;
    }

    CheckAvailabilityOfNeighborsCABAC();
    
#ifndef SIMPLIFY_CODE
    //========= write mb_aff (CABAC) =========
    if(img->MbaffFrameFlag && !skip) // check for copy mode
    {
      if(WriteFrameFieldMBInHeader)
      {
        currSE->value1 = currMB->mb_field;
        currSE->value2 = 0;
        currSE->type   =  SE_MBTYPE;
        
        currSE->writing = writeFieldModeInfo_CABAC;
        dataPart->writeSyntaxElement(currSE, dataPart);
#if TRACE
        snprintf(currSE->tracestring, TRACESTRING_SIZE, "mb_field_decoding_flag");
#endif
        bitCount[BITS_MB_MODE] += currSE->len;
        no_bits                += currSE->len;
        currSE++;
        currMB->currSEnr++;
      }
    }
#endif

    //each 16x16 block mode
    if(mb32==0)
    {    
      //========= write mb_type (CABAC) =========
      if (currMB->mb_type != 0 || ((img->type == B_SLICE) && currMB->cbp != 0))
      {
        currSE->value1  = mb_type;
        currSE->value2  = 0;
        currSE->type    = SE_MBTYPE;
        currSE->writing = writeMB_typeInfo_CABAC;
        dataPart->writeSyntaxElement( currSE, dataPart);
  #if TRACE
        if (img->type == B_SLICE) 
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "mb_type (B_SLICE) (%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
        else                      
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "mb_type (P_SLICE) (%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
  #endif
        bitCount[BITS_MB_MODE] += currSE->len;
        no_bits                += currSE->len;
        currSE++;
        currMB->currSEnr++;
        
      }
    }
  }

  
#ifdef MV_COMPETITION
#ifdef MB32X32_MVC
  if (input->mv_competition > 0)
  {
    if((rdopt==0) && (input->symbol_mode == CABAC) && (img->type==P_SLICE) && (mv_comp.nb_mode_for_skip > 1))  
    {  
      //make sure that this function is called only once...whether it is big blocks or 16x16
      //32x32 P-Skip
      if(mb32==1 && mbi==0 && MB32.mb_type==0)
          no_bits  += write_predictor_index_for_skip_mode32();
      //16x16 P-Skip
      if(mb32==0 && currMB->mb_type==0)
          no_bits  += write_predictor_index_for_skip_mode  ();
    }
  }
#else
  if (input->mv_competition > 0 && !mb32)
  {
    if((rdopt==0)&&(input->symbol_mode == CABAC))  
    {  
      if ((currMB->mb_type == 0) && (img->type==P_SLICE) && (mv_comp.nb_mode_for_skip > 1))
        no_bits  += write_predictor_index_for_skip_mode();
    }
  }
#endif
#endif
  
  
  //init NoMbPartLessThan8x8Flag
  currMB->NoMbPartLessThan8x8Flag = (IS_DIRECT(currMB) && !(active_sps->direct_8x8_inference_flag))? 0: 1;
  
  if (currMB->mb_type == IPCM)
  {
    int jj, uv;
    if (input->symbol_mode == CABAC) 
    {
      int len;
      EncodingEnvironmentPtr eep = &dataPart->ee_cabac;
      len = arienco_bits_written(eep);
      arienco_done_encoding(eep); // This pads to byte
      len = arienco_bits_written(eep) - len;
      currSE--; // Need to append more bits to mb_type
      currSE->len += len;
      no_bits += len;
      currSE++;
      // Now restart the encoder
      arienco_start_encoding(eep, dataPart->bitstream->streamBuffer, &(dataPart->bitstream->byte_pos));
      reset_pic_bin_count();
    }
    if (dataPart->bitstream->bits_to_go < 8)
    {
      // This will only happen in the CAVLC case, CABAC is already padded
      currSE->mapping = ue_linfo;
      currSE->len = dataPart->bitstream->bits_to_go;
      no_bits += currSE->len;
      bitCount[BITS_COEFF_Y_MB]+= currSE->len;
      currSE->bitpattern = 0;
#if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "IPCM alignment bits = %d", currSE->len);
#endif
      writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
      currSE++;
      currMB->currSEnr++;
    }
    for (j=0;j<MB_BLOCK_SIZE;j++)
    {
      jj = img->pix_y+j;
      for (i=0;i<MB_BLOCK_SIZE;i++)
      {
        currSE->mapping = ue_linfo;
        currSE->len = img->bitdepth_luma;  
        currSE->type    = SE_MBTYPE;
        no_bits += currSE->len;
        currSE->bitpattern = enc_picture->imgY[jj][img->pix_x+i];       
        currSE->value1 = currSE->bitpattern;
        currSE->type = SE_MBTYPE;
        bitCount[BITS_COEFF_Y_MB]+=currSE->len;
#if TRACE
        snprintf(currSE->tracestring, TRACESTRING_SIZE, "IPCM Luma (%d %d) = %d", j,i,currSE->bitpattern);
#endif
        writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);        
        currSE++;
        currMB->currSEnr++;
      }
    }
    if (img->yuv_format != YUV400)
    {
      for (uv = 0; uv < 2; uv ++)
      {
        for (j=0;j<img->mb_cr_size_y;j++)
        {
          jj = img->pix_c_y+j;
          for (i=0;i<img->mb_cr_size_x;i++)
          {
            currSE->mapping = ue_linfo;
            currSE->len = img->bitdepth_chroma;
            currSE->type    = SE_MBTYPE;
            no_bits += currSE->len;
            currSE->bitpattern = enc_picture->imgUV[uv][jj][img->pix_c_x+i];
            currSE->value1 = currSE->bitpattern;
            currSE->type = SE_MBTYPE;
            bitCount[BITS_COEFF_UV_MB]+=currSE->len;
#if TRACE
            snprintf(currSE->tracestring, TRACESTRING_SIZE, "IPCM chroma(%d) (%d %d) = %d", uv, j,i,currSE->bitpattern);
#endif
            writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
            currSE++;
            currMB->currSEnr++;
          }
        }
      }
    }
    return no_bits;
  }
  
  //===== BITS FOR 8x8 SUB-PARTITION MODES =====
  if (IS_P8x8 (currMB))
  {
    dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);
    
    for (i=0; i<4; i++)
    {
      if (input->symbol_mode==UVLC)   
        currSE->mapping = ue_linfo;
      else
        currSE->writing = writeB8_typeInfo_CABAC;
      
      currSE->value1  = B8Mode2Value (currMB->b8mode[i], currMB->b8pdir[i]);
      currSE->value2  = 0;
      currSE->type    = SE_MBTYPE;
      dataPart->writeSyntaxElement (currSE, dataPart);
#if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "8x8 mode/pdir(%2d) = %3d/%d", i, currMB->b8mode[i], currMB->b8pdir[i]);
#endif
      bitCount[BITS_MB_MODE]+= currSE->len;
      no_bits               += currSE->len;
      currSE++;
      currMB->currSEnr++;
      
      //set NoMbPartLessThan8x8Flag for P8x8 mode
      currMB->NoMbPartLessThan8x8Flag &= (currMB->b8mode[i]==0 && active_sps->direct_8x8_inference_flag) || 
        (currMB->b8mode[i]==4);
    }
#ifdef MV_COMPETITION
#ifdef MB32X32_MVC
    if (input->mv_competition > 0)
#else
    if (input->mv_competition > 0 && !mb32)
#endif
      pass_with_writing = !rdopt;
#endif
    no_bits += writeMotionInfo2NAL  ();
    currSE   = &img->MB_SyntaxElements[currMB->currSEnr];
  }
  
  //============= Transform size flag for INTRA MBs =============
  //-------------------------------------------------------------
  //transform size flag for INTRA_4x4 and INTRA_8x8 modes
  if ((currMB->mb_type == I8MB || currMB->mb_type == I4MB) && input->Transform8x8Mode)
  {
    currSE->value1 = currMB->luma_transform_size_8x8_flag;
    currSE->type   = SE_HEADER;
    
    if( input->symbol_mode==UVLC)
    {
      currSE->mapping = ue_linfo;
      currSE->type    = SE_MBTYPE;
      currSE->bitpattern = currMB->luma_transform_size_8x8_flag;
      currSE->len = 1;
      writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
    }
    else
    {
      currSE->writing = writeMB_transform_size_CABAC;
      dataPart->writeSyntaxElement(currSE, dataPart);
    }
#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "transform_size_8x8_flag = %3d", currMB->luma_transform_size_8x8_flag);
#endif
    
    bitCount[BITS_MB_MODE] += currSE->len;
    no_bits                += currSE->len;
    currSE++;
    currMB->currSEnr++;
  }
  
  
  //===== BITS FOR INTRA PREDICTION MODES ====
  no_bits += writeIntraModes();
  //===== BITS FOR CHROMA INTRA PREDICTION MODE ====
  if (currMB->IntraChromaPredModeFlag && img->yuv_format != YUV400)
    no_bits += writeChromaIntraPredMode();
  else if(!rdopt) //GB CHROMA !!!!!
    currMB->c_ipred_mode = DC_PRED_8; //setting c_ipred_mode to default is not the right place here
  //resetting in rdopt.c (but where ??)
  //with cabac and bframes maybe it could crash without this default
  //since cabac needs the right neighborhood for the later MBs
  
  //----- motion information -----
  if(mb32==1 && mbi==0)
  {
    if(MB32.mb_type != 0)
    {
#ifdef MB32X32_MVC
      if (input->mv_competition > 0)
        pass_with_writing = !rdopt;
#endif
      no_bits  += writeMotionInfo2NAL32  (1);
    }
  }
  else if(mb32==0 && currMB->mb_type !=0 && currMB->mb_type !=P8x8)
  {
#ifdef MV_COMPETITION
#ifdef MB32X32_MVC
    if (input->mv_competition > 0)
#else
    if (input->mv_competition > 0 && !mb32)
#endif
      pass_with_writing = !rdopt;
#endif

    no_bits  += writeMotionInfo2NAL  ();
  }

  if(  (mb32 && (MB32.mb_type!=0 || (img->type==B_SLICE && MB32.cbp!=0)))
    || (!mb32 && /*MB32.cbp!=0 && */(currMB->mb_type!=0 || (img->type==B_SLICE && currMB->cbp!=0)))  )
  {
    if(mb32 && mbi==0)
    {
      no_bits += writeCBP32(rdopt, MB32.cbp, !delta_qp_sent);//with dquant
    }

    if( !mb32 || (mb32 && MB32.cbp!=0) )
    {
#ifdef USE_INTRA_MDDT
      *coeff_rate = writeCBPandLumaCoeff (rdopt, mb32, mbi);
#else
      *coeff_rate = writeCBPandLumaCoeff ();
#endif
      if (img->yuv_format != YUV400)
        *coeff_rate  += writeChromaCoeff ();
    
      no_bits  += *coeff_rate;  
    }
  }
#ifdef ADAPTIVE_QUANTIZATION
  if(!rdopt)
  {
    dwSkipMbCount+=(currMB->cbp==0);
    dwSkipEnableMbCount++;
    
    dwSkipEnableMbCountForSlice[(img->type)+(img->type==B_SLICE? (img->nal_reference_idc<<2) : 0)]++;
    dwSkipMbCountForSlice[(img->type)+(img->type==B_SLICE? (img->nal_reference_idc<<2) : 0)]+=(currMB->cbp==0);
  }
#endif
#ifdef ADAPTIVE_FD_SD_CODING
  else
  {
    if (cabac_encoding)
    {
      currMB->SD_or_FD_t8x8=0;
    }
  }  
#endif
  return no_bits;
}

/*!
************************************************************************
* \brief
*    Passes the chosen syntax elements to the NAL for 16x16 block inside of 32x32 block
* \param coeff_rate
*    bitrate of Luma and Chroma coeff
* \param (mb32,mbi)
*        (1,0) - 32x32, 32x16, 16x32 mode info, mv info, coeff of upper left 16x16, dequant
*        (1,1~)-                                         coeff of other      16x16
*        (0,0) - 16x16 mode, upper left 16x16 partition mode, mv, coeff, dequant
*        (0,1~)-                        16x16 partition mode, mv, coeff
************************************************************************
*/
int write_one_macroblock32 (int eos_bit, int mb32, int mbi)
{
  Macroblock* currMB   = &img->mb_data[img->current_mb_nr];
  int*        bitCount = currMB->bitcounter;
  int i;
  int no_bits;
  
  extern int cabac_encoding;

  //===== init and update number of intra macroblocks =====
  if (img->current_mb_nr==0)
    intras=0;
  
  if (IS_INTRA(currMB))
    intras++;
  
  //--- write non-slice termination symbol if the macroblock is not the first one in its slice ---
  if (input->symbol_mode==CABAC && img->current_mb_nr!=img->currentSlice->start_mb_nr && eos_bit && mbi==0 && eos_bit_32)
  {
    write_terminating_bit (0);
  }
  
  cabac_encoding = 1;
  
#ifdef ADAPTIVE_FD_SD_CODING
  if (cabac_encoding) currMB->written_SD_Coding_on_off=0;
#endif
#ifdef ADAPTIVE_QUANTIZATION
  // Condieration for no residual
  if ( img->slice_fractional_quant_flag && currMB->mb_type==IPCM)
    currMB->mb_iaqms_idx = 0;
  if(img->slice_fractional_quant_flag && (currMB->cbp==0 && !IS_NEWINTRA (currMB)) )
    currMB->mb_iaqms_idx = 0;
#endif
  
  //--- write macroblock ---
  no_bits = writeMBLayer32 (0, &i, mb32, mbi);  // i is temporary
  if (!((currMB->mb_type !=0 ) || ((img->type==B_SLICE) && currMB->cbp != 0) ))
  { 
    for (i=0; i < 4; i++)
      memset(img->nz_coeff [img->current_mb_nr][i], 0, (4 + img->num_blk8x8_uv) * sizeof(int));  // CAVLC
  }
    
  //--- set total bit-counter ---
  
  
  
#ifdef MV_COMPETITION
#ifdef MB32X32_MVC
  if (input->mv_competition > 0)
#else
  if (input->mv_competition > 0 && !mb32)
#endif
  {
    bitCount[BITS_TOTAL_MB] = bitCount[BITS_MB_MODE]  + bitCount[BITS_COEFF_Y_MB]     
      + bitCount[BITS_INTER_MB] + bitCount[BITS_CBP_MB]  
      + bitCount[BITS_DELTA_QUANT_MB] + bitCount[BITS_COEFF_UV_MB]
      + bitCount[BITS_SKIP_PRED]
      + bitCount[BITS_MV_PRED];
    
    //Rate control
    img->NumberofMBHeaderBits=bitCount[BITS_MB_MODE]   + bitCount[BITS_INTER_MB]
      + bitCount[BITS_CBP_MB]  + bitCount[BITS_DELTA_QUANT_MB]
      + bitCount[BITS_SKIP_PRED]
      + bitCount[BITS_MV_PRED];
  }
  else
#endif
    
  {
    bitCount[BITS_TOTAL_MB] = bitCount[BITS_MB_MODE]  + bitCount[BITS_COEFF_Y_MB]     
      + bitCount[BITS_INTER_MB] + bitCount[BITS_CBP_MB]  
      + bitCount[BITS_DELTA_QUANT_MB] + bitCount[BITS_COEFF_UV_MB];
    
    //Rate control
    img->NumberofMBHeaderBits=bitCount[BITS_MB_MODE]   + bitCount[BITS_INTER_MB]
      + bitCount[BITS_CBP_MB]  + bitCount[BITS_DELTA_QUANT_MB];
  }  
  ;
  img->NumberofMBTextureBits= bitCount[BITS_COEFF_Y_MB]+ bitCount[BITS_COEFF_UV_MB];
  img->NumberofTextureBits +=img->NumberofMBTextureBits;
  img->NumberofHeaderBits +=img->NumberofMBHeaderBits;
  /*basic unit layer rate control*/
  if(img->BasicUnit<img->Frame_Total_Number_MB)
  {
    img->NumberofBasicUnitHeaderBits +=img->NumberofMBHeaderBits;
    img->NumberofBasicUnitTextureBits +=img->NumberofMBTextureBits;
  }
  /*record the total number of MBs*/
  img->NumberofCodedMacroBlocks++;
  
  stats->bit_slice += bitCount[BITS_TOTAL_MB];
  
  cabac_encoding = 0;



  return no_bits;
}



int writeLumaCoeff16x16P_CABAC ()
{
  int             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int*            bitCount  = currMB->bitcounter;
  DataPartition*  dataPart;
  
  int   level, run;
  int   k;
  int*  ACLevel = img->cofAC16x16[0];
  int*  ACRun   = img->cofAC16x16[1];
  
  img->subblock_x = 0;//img->pix_x;  // horiz. position for coeff_count context
  img->subblock_y = 0;//img->pix_y;  // vert.  position for coeff_count context



  level=1; // get inside loop
  for(k=0; k<=256 && level !=0; k++)
  {
    level = currSE->value1 = ACLevel[k]; // level
    run   = currSE->value2 = ACRun  [k]; // run
    
    currSE->writing = writeRunLevel16x16P_CABAC;
    
    currSE->context     = LUMA_16x16P; 
    currSE->type        = SE_LUM_AC_INTER;
    img->is_intra_block = 0;    // only applicable to inter blocks
    
    // choose the appropriate data partition
    dataPart = &(currSlice->partArr[partMap[img->type != B_SLICE ? currSE->type : SE_BFRAME]]);
    dataPart->writeSyntaxElement (currSE, dataPart);
    bitCount[BITS_COEFF_Y_MB] += currSE->len;
    rate                      += currSE->len;



#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "Luma8x8 sng(%2d) level =%3d run =%2d", k, level,run);
#endif
    /* proceed to next SE */
    currSE++;  
    currMB->currSEnr++;
    
  }

  
  return rate;
}

int writeLumaCoeff16x16P ()
{
  int  rate = 0;
  
  if( input->symbol_mode == UVLC) // allow here if 4x4 or UVLC
  {
    assert(0);
  }
  else 
    rate += writeLumaCoeff16x16P_CABAC ();
  
  return rate;
}

int writeLumaCoeff16x8_CABAC (int blk)
{
  int             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int*            bitCount  = currMB->bitcounter;
  DataPartition*  dataPart;
  
  int   level, run;
  int   k;
  int*  ACLevel = img->cofAC16x8[blk][0];
  int*  ACRun   = img->cofAC16x8[blk][1];
  
  img->subblock_x = 0;               // horiz. position for coeff_count context
  img->subblock_y = (blk == 0)?0:2;  // vert.  position for coeff_count context
  
  level=1; // get inside loop
  for(k=0; k<=128 && level !=0; k++)
  {
    level = currSE->value1 = ACLevel[k]; // level
    run   = currSE->value2 = ACRun  [k]; // run
    
    currSE->writing = writeRunLevel16x8_CABAC;
    
    currSE->context     = LUMA_16x8P;

    currSE->type        = SE_LUM_AC_INTER;
    img->is_intra_block = 0;    // only applicable to inter blocks
    
    // choose the appropriate data partition
    dataPart = &(currSlice->partArr[partMap[img->type != B_SLICE ? currSE->type : SE_BFRAME]]);
    dataPart->writeSyntaxElement (currSE, dataPart);
    bitCount[BITS_COEFF_Y_MB] += currSE->len;
    rate                      += currSE->len;
#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "Luma8x8 sng(%2d) level =%3d run =%2d", k, level,run);
#endif
    /* proceed to next SE */
    currSE++;  
    currMB->currSEnr++;
    
  }
  
  return rate;
}

int writeLumaCoeff16x8 ( int blk16x8 )
{
  int  rate = 0;
  
  if( input->symbol_mode == UVLC) // allow here if 4x4 or UVLC
  {
    assert(0);
  }
  else 
    rate += writeLumaCoeff16x8_CABAC (blk16x8);
  
  return rate;
}

int writeLumaCoeff8x16_CABAC (int blk)
{
  int             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int*            bitCount  = currMB->bitcounter;
  DataPartition*  dataPart;
  
  int   level, run;
  int   k;
  int*  ACLevel = img->cofAC8x16[blk][0];
  int*  ACRun   = img->cofAC8x16[blk][1];
  
  img->subblock_x = (blk == 0)?0:2;     // horiz. position for coeff_count context
  img->subblock_y = 0;                  // vert.  position for coeff_count context
  
  level=1; // get inside loop
  for(k=0; k<=128 && level !=0; k++)
  {
    level = currSE->value1 = ACLevel[k]; // level
    run   = currSE->value2 = ACRun  [k]; // run
    
    currSE->writing = writeRunLevel8x16_CABAC;
    
    currSE->context     = LUMA_8x16P;

    currSE->type        = SE_LUM_AC_INTER;
    img->is_intra_block = 0;    // only applicable to inter blocks
    
    // choose the appropriate data partition
    dataPart = &(currSlice->partArr[partMap[img->type != B_SLICE ? currSE->type : SE_BFRAME]]);
    dataPart->writeSyntaxElement (currSE, dataPart);
    bitCount[BITS_COEFF_Y_MB] += currSE->len;
    rate                      += currSE->len;
#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "Luma8x8 sng(%2d) level =%3d run =%2d", k, level,run);
#endif
    /* proceed to next SE */
    currSE++;  
    currMB->currSEnr++;
    
  }
  
  return rate;
}

int writeLumaCoeff8x16 ( int blk8x16 )
{
  int  rate = 0;
  
  if( input->symbol_mode == UVLC) 
  {
    assert(0);
  }
  else 
    rate += writeLumaCoeff8x16_CABAC (blk8x16);
  
  return rate;
}
int writeLumaCoefBigBlock(Macroblock *currMB, int mb32, int rdopt)
{
  int rate = 0;

  assert(currMB->mb_type >= 0 && currMB->mb_type <= 3);
  assert(currMB->luma_transform_size_8x8_flag);
  if((currMB->mb_type < 2 || mb32)&& currMB->luma_transform_size_8x8_flag == 2)  // 16x16 partition
  {
    if(currMB->cbp & 0xf)
      rate += writeLumaCoeff16x16P(); 
  }
  else if(currMB->mb_type == 2 && currMB->luma_transform_size_8x8_flag == 2)  // 16x8 partition
  {
    if(currMB->cbp & 0x3)
      rate += writeLumaCoeff16x8(0);
    if(currMB->cbp & 0xc)
      rate += writeLumaCoeff16x8(1);
  }
  else if(currMB->mb_type == 3 && currMB->luma_transform_size_8x8_flag == 2)  // 16x8 partition
  {
    if(currMB->cbp & 0x5)
      rate += writeLumaCoeff8x16(0);
    if(currMB->cbp & 0xa)
      rate += writeLumaCoeff8x16(1);
  }

  return rate;
}


int writelumaTransform_size()
{
  //int             mb_x, mb_y, i, j, k;
  //int             level, run;
  int             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  int*            bitCount  = currMB->bitcounter;
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  DataPartition*  dataPart;
  //int             need_transform_size_flag;   //ADD-VG-24062004
  

    if (IS_OLDINTRA (currMB) || currMB->mb_type == SI4MB ||  currMB->mb_type == I8MB)
    {
      if (input->symbol_mode == UVLC)  currSE->mapping = cbp_linfo_intra;
      currSE->type = SE_CBP_INTRA;
    }
    else
    {
      if (input->symbol_mode == UVLC)  currSE->mapping = cbp_linfo_inter;
      currSE->type = SE_CBP_INTER;
    }    

    // choose the appropriate data partition
    dataPart = &(currSlice->partArr[partMap[currSE->type]]);
    
    //if (need_transform_size_flag)
    {
      currSE->value1 = currMB->luma_transform_size_8x8_flag;
      currSE->type   = SE_HEADER;
      
      if( input->symbol_mode==UVLC)
      {
        currSE->mapping = ue_linfo;
        currSE->bitpattern = currMB->luma_transform_size_8x8_flag;
        currSE->len = 1;
        writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
      }
      else
      {
        currSE->writing = writeMB_transform_size_CABAC;
        dataPart->writeSyntaxElement(currSE, dataPart);
      }
#if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "transform size 8x8 flag = %3d", currMB->luma_transform_size_8x8_flag);
#endif
      
      bitCount[BITS_MB_MODE] += currSE->len;
      rate                   += currSE->len;
      currSE++;
      currMB->currSEnr++;
    }

    return rate;
}

// write 1 bit indicating luma cbp
int writelumaCBPTransform_size()
{
  //int             mb_x, mb_y, i, j, k;
  //int             level, run;
  int             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  int*            bitCount  = currMB->bitcounter;
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int             cbp       = currMB->cbp;
  DataPartition*  dataPart;
  int             need_transform_size_flag;   //ADD-VG-24062004
  

    //=====   C B P   =====
    //---------------------
    currSE->value1 = cbp;

    if (IS_OLDINTRA (currMB) || currMB->mb_type == SI4MB ||  currMB->mb_type == I8MB)
    {
      if (input->symbol_mode == UVLC)  currSE->mapping = cbp_linfo_intra;
      currSE->type = SE_CBP_INTRA;
    }
    else
    {
      if (input->symbol_mode == UVLC)  currSE->mapping = cbp_linfo_inter;
      currSE->type = SE_CBP_INTER;
    }
    if (input->symbol_mode == CABAC)   currSE->writing = writeCBP_CABAC_luma1bit;
    
    // choose the appropriate data partition
    dataPart = &(currSlice->partArr[partMap[currSE->type]]);
    
    dataPart->writeSyntaxElement(currSE, dataPart);
    bitCount[BITS_CBP_MB] += currSE->len;
    rate                  += currSE->len;
#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "CBP (%2d,%2d) = %3d",img->mb_x, img->mb_y, cbp);
#endif
    // proceed to next SE
    currSE++;  
    currMB->currSEnr++;

    //============= Transform Size Flag for INTER MBs =============
    //-------------------------------------------------------------
    need_transform_size_flag = (((currMB->mb_type >= 1 && currMB->mb_type <= 3)||
      (IS_DIRECT(currMB) && active_sps->direct_8x8_inference_flag) ||
      (currMB->NoMbPartLessThan8x8Flag))
      && currMB->mb_type != I8MB && currMB->mb_type != I4MB
      && (currMB->cbp&15)
      && input->Transform8x8Mode);
    
    if (need_transform_size_flag)
    {
      currSE->value1 = currMB->luma_transform_size_8x8_flag;
      currSE->type   = SE_HEADER;
      
      if( input->symbol_mode==UVLC)
      {
        currSE->mapping = ue_linfo;
        currSE->bitpattern = currMB->luma_transform_size_8x8_flag;
        currSE->len = 1;
        writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
      }
      else
      {
        currSE->writing = writeMB_transform_size_CABAC;
        dataPart->writeSyntaxElement(currSE, dataPart);
      }
#if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "transform size 8x8 flag = %3d", currMB->luma_transform_size_8x8_flag);
#endif
      
      bitCount[BITS_MB_MODE] += currSE->len;
      rate                   += currSE->len;
      currSE++;
      currMB->currSEnr++;
    }

    return rate;
}
void writeCBP_CABAC_chroma(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  TextureInfoContexts *ctx = img->currentSlice->tex_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  
  int a, b;
  int curr_cbp_ctx;
  int cbp = se->value1; // symbol to encode
  int cbp_bit;
  //int b8;
  
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

int write_chroma_cbp()
{
  //int             mb_x, mb_y, i, j, k;
  //int             level, run;
  int             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  int*            bitCount  = currMB->bitcounter;
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int             cbp       = currMB->cbp;
  DataPartition*  dataPart;
  //int             need_transform_size_flag;   //ADD-VG-24062004
  

    //=====   C B P   =====
    //---------------------
    currSE->value1 = cbp;
 
    if (IS_OLDINTRA (currMB) || currMB->mb_type == SI4MB ||  currMB->mb_type == I8MB)
    {
      if (input->symbol_mode == UVLC)  currSE->mapping = cbp_linfo_intra;
      currSE->type = SE_CBP_INTRA;
    }
    else
    {
      if (input->symbol_mode == UVLC)  currSE->mapping = cbp_linfo_inter;
      currSE->type = SE_CBP_INTER;
    }
    if (input->symbol_mode == CABAC)   currSE->writing = writeCBP_CABAC_chroma;
    
    // choose the appropriate data partition
    dataPart = &(currSlice->partArr[partMap[currSE->type]]);
    
    dataPart->writeSyntaxElement(currSE, dataPart);
    bitCount[BITS_CBP_MB] += currSE->len;
    rate                  += currSE->len;
#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "CBP (%2d,%2d) = %3d",img->mb_x, img->mb_y, cbp);
#endif
    // proceed to next SE
    currSE++;  
    currMB->currSEnr++;

    return rate;
}


void calcCoeff16x16I (int coeffI[16][16])
{
  int i, j;
  double a = M_PI/32.;
  double b = 1./(sqrt(2.));
  double c = 1./sqrt(8.);
  double coeff[16][16];
  int factor = 1 << DCT16PREC;

  for(j = 0; j < 16; j++)
    for(i = 0; i < 16; i++)
      coeff[j][i] = cos((2*i+1)*j*a)*c; 

  for(i = 0; i < 16; i++)
    coeff[0][i] *= b; 

  for(j = 0; j < 16; j++)
    for(i = 0; i < 16; i++)
    {
      coeffI[j][i] = (int)(fabs(coeff[j][i])*factor+0.5); 
      if(coeff[j][i] < 0.)
        coeffI[j][i] = - coeffI[j][i];
    }
}

void calcCoeff8x8I(int coeffI[16][16])
{
  int i, j;
  double a = M_PI/16.;
  double b = 1./(sqrt(2.));
  double c = 1./sqrt(4.);
  double coeff[16][16];
  int factor = 1 << DCT16PREC;

  for(j = 0; j < 8; j++)
    for(i = 0; i < 8; i++)
      coeff[j][i] = cos((2*i+1)*j*a)*c; 

  for(i = 0; i < 8; i++)
    coeff[0][i] *= b; 

  for(j = 0; j < 8; j++)
    for(i = 0; i < 8; i++)
    {
      coeffI[j][i] = (int)(fabs(coeff[j][i])*factor+0.5); 
      if(coeff[j][i] < 0.)
        coeffI[j][i] = - coeffI[j][i];
    }
}


static void forward16x16I (int (*block) [16], int (*tblock)[16])
{
  int x, y;
  int k; 
  int temp[16][16];
  int tmp; 

  // horizontal 
  for(y = 0; y < 16; y++)
  {
    for(x = 0; x < 16; x++)
    {
      temp[y][x] = 0;
      for(k = 0; k < 16; k++)
        temp[y][x] += img->dct_coeff16I[y][k]*block[x][k];
    }
  }
  
  // vertical 
  for(y = 0; y < 16; y++)
  {
    for(x = 0; x < 16; x++)
    {
      tmp = 0;
      for(k = 0; k < 16; k++)
        tmp += img->dct_coeff16I[y][k]*temp[x][k];
      tblock[y][x] = tmp;
    }
  }
}

static void inverse16x16I (int (*tblock) [16], int (*block)[16])
{
  int x, y, k;
  int temp[16][16];
  int tmp;

  // horizontal 
  for(y = 0; y < 16; y++)
  {
    for(x = 0; x < 16; x++)
    {
      temp[y][x] = 0;
      for(k = 0; k < 16; k++)
        temp[y][x] += img->dct_coeff16I[k][y]*tblock[x][k];
    }
  }
  
  // vertical 
  for(y = 0; y < 16; y++)
  {
    for(x = 0; x < 16; x++)
    {
      tmp = 0;
      for(k = 0; k < 16; k++)
        tmp += img->dct_coeff16I[k][y]*temp[x][k];
      block[y][x] = tmp;
    }
  }
}


const byte SNGL_SCAN16x16[256][2] = 
{
  { 0, 0}, { 1, 0}, { 0, 1}, { 0, 2}, { 1, 1}, { 2, 0}, { 3, 0}, { 2, 1}, { 1, 2}, { 0, 3}, { 0, 4}, { 1, 3}, { 2, 2}, { 3, 1}, { 4, 0}, { 5, 0},
  { 4, 1}, { 3, 2}, { 2, 3}, { 1, 4}, { 0, 5}, { 0, 6}, { 1, 5}, { 2, 4}, { 3, 3}, { 4, 2}, { 5, 1}, { 6, 0}, { 7, 0}, { 6, 1}, { 5, 2}, { 4, 3},
  { 3, 4}, { 2, 5}, { 1, 6}, { 0, 7}, { 0, 8}, { 1, 7}, { 2, 6}, { 3, 5}, { 4, 4}, { 5, 3}, { 6, 2}, { 7, 1}, { 8, 0}, { 9, 0}, { 8, 1}, { 7, 2}, 
  { 6, 3}, { 5, 4}, { 4, 5}, { 3, 6}, { 2, 7}, { 1, 8}, { 0, 9}, { 0,10}, { 1, 9}, { 2, 8}, { 3, 7}, { 4, 6}, { 5, 5}, { 6, 4}, { 7, 3}, { 8, 2},
  { 9, 1}, {10, 0}, {11, 0}, {10, 1}, { 9, 2}, { 8, 3}, { 7, 4}, { 6, 5}, { 5, 6}, { 4, 7}, { 3, 8}, { 2, 9}, { 1,10}, { 0,11}, { 0,12}, { 1,11},
  { 2,10}, { 3, 9}, { 4, 8}, { 5, 7}, { 6, 6}, { 7, 5}, { 8, 4}, { 9, 3}, {10, 2}, {11, 1}, {12, 0}, {13, 0}, {12, 1}, {11, 2}, {10, 3}, { 9, 4},
  { 8, 5}, { 7, 6}, { 6, 7}, { 5, 8}, { 4, 9}, { 3,10}, { 2,11}, { 1,12}, { 0,13}, { 0,14}, { 1,13}, { 2,12}, { 3,11}, { 4,10}, { 5, 9}, { 6, 8},
  { 7, 7}, { 8, 6}, { 9, 5}, {10, 4}, {11, 3}, {12, 2}, {13, 1}, {14, 0}, {15, 0}, {14, 1}, {13, 2}, {12, 3}, {11, 4}, {10, 5}, { 9, 6}, { 8, 7},
  { 7, 8}, { 6, 9}, { 5,10}, { 4,11}, { 3,12}, { 2,13}, { 1,14}, { 0,15}, { 1,15}, { 2,14}, { 3,13}, { 4,12}, { 5,11}, { 6,10}, { 7, 9}, { 8, 8},
  { 9, 7}, {10, 6}, {11, 5}, {12, 4}, {13, 3}, {14, 2}, {15, 1}, {15, 2}, {14, 3}, {13, 4}, {12, 5}, {11, 6}, {10, 7}, { 9, 8}, { 8, 9}, { 7,10},
  { 6,11}, { 5,12}, { 4,13}, { 3,14}, { 2,15}, { 3,15}, { 4,14}, { 5,13}, { 6,12}, { 7,11}, { 8,10}, { 9, 9}, {10, 8}, {11, 7}, {12, 6}, {13, 5},
  {14, 4}, {15, 3}, {15, 4}, {14, 5}, {13, 6}, {12, 7}, {11, 8}, {10, 9}, { 9,10}, { 8,11}, { 7,12}, { 6,13}, { 5,14}, { 4,15}, { 5,15}, { 6,14},
  { 7,13}, { 8,12}, { 9,11}, {10,10}, {11, 9}, {12, 8}, {13, 7}, {14, 6}, {15, 5}, {15, 6}, {14, 7}, {13, 8}, {12, 9}, {11,10}, {10,11}, { 9,12},
  { 8,13}, { 7,14}, { 6,15}, { 7,15}, { 8,14}, { 9,13}, {10,12}, {11,11}, {12,10}, {13, 9}, {14, 8}, {15, 7}, {15, 8}, {14, 9}, {13,10}, {12,11},
  {11,12}, {10,13}, { 9,14}, { 8,15}, { 9,15}, {10,14}, {11,13}, {12,12}, {13,11}, {14,10}, {15, 9}, {15,10}, {14,11}, {13,12}, {12,13}, {11,14},
  {10,15}, {11,15}, {12,14}, {13,13}, {14,12}, {15,11}, {15,12}, {14,13}, {13,14}, {12,15}, {13,15}, {14,14}, {15,13}, {15,14}, {14,15}, {15,15}
};


int dct_luma_16x16_new (void)
{
  int nonzero = FALSE;  
  int i,j,coeff_ctr;
  int scan_pos,run;

  int*  ACLevel = img->cofAC16x16[0];
  int*  ACRun   = img->cofAC16x16[1];

  imgpel **img_enc      = enc_picture->imgY;
  int    (*curr_res)[MB_BLOCK_SIZE] = img->m7; 

  //Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  //int   max_imgpel_value = img->max_imgpel_value;
  const byte (*pos_scan)[2] = SNGL_SCAN16x16;
  int qp_per    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)/6; //CHA-VG01
  int qp_rem    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)%6; //CHA-VG01
  //int q_bits = Q_BITS + qp_per;
  int DCTprec = 2*DCT16PREC; 

  int R[6] = {40, 45, 50, 57, 63, 71};
  int Q[6] = {410,   364,   328,   287,   260,   231};
  int qpShift, qpRound;
  int level;
  int M1[MB_BLOCK_SIZE][MB_BLOCK_SIZE];

  int k;
  levelDataStruct levelData[256];
  double lambda_md = 0.0;
  double err;
  int lowerInt;
  int levelRDOQ[256];
  int kStart = 0;
  int kStop = 0;
  int noCoeff;
  double normFact = pow(2, 2 * (DCTprec + 14)-15); 

  //  Forward 16x16 transform
  forward16x16I(curr_res, M1);

  // Quantization process
  qpShift=qp_per+DCTprec+8;
  qpRound=(1<<(qpShift))/6;
  

  // quantization and zig-zag scanning 
  if(input->UseRDO_Q)
  {
    if ((img->type == B_SLICE) && img->nal_reference_idc)
    {
      lambda_md = img->lambda_md[5][img->masterQP];  
    }
    else
    {
      lambda_md = img->lambda_md[img->type][img->masterQP];
    }

    noCoeff = 0;
    for (coeff_ctr = 0; coeff_ctr < 256; coeff_ctr++)
    {
      i = pos_scan[coeff_ctr][0];
      j = pos_scan[coeff_ctr][1];

      levelData[coeff_ctr].levelDouble = (int64)absm(M1[j][i]) * Q[qp_rem];
      level = (int)(levelData[coeff_ctr].levelDouble >> qpShift);
      lowerInt = ((levelData[coeff_ctr].levelDouble - (level << qpShift)) < ((int64)1 << (qpShift - 1)))? 1: 0;

      levelData[coeff_ctr].level[0] = 0;
      if (level == 0 && lowerInt == 1)
      {
        levelData[coeff_ctr].noLevels = 1;
      }
      else if (level == 0 && lowerInt == 0)
      {
        levelData[coeff_ctr].level[1] = level + 1;
        levelData[coeff_ctr].noLevels = 2;
        kStop = coeff_ctr;
        noCoeff++;
      }
      else if (level > 0 && lowerInt == 1)
      {
        levelData[coeff_ctr].level[1] = level;
        levelData[coeff_ctr].noLevels = 2;
        kStop = coeff_ctr;
        noCoeff++;
      }
      else
      {
        levelData[coeff_ctr].level[1] = level;
        levelData[coeff_ctr].level[2] = level + 1;
        levelData[coeff_ctr].noLevels = 3;
        kStop = coeff_ctr;
        kStart = coeff_ctr;
        noCoeff++;
      }

      for (k = 0; k < levelData[coeff_ctr].noLevels; k++)
      {      
        err = (double)(levelData[coeff_ctr].level[k] << qpShift) - (double) levelData[coeff_ctr].levelDouble;
        err = err * R[qp_rem];
        levelData[coeff_ctr].errLevel[k] = (err * err) / normFact; 
      }
    }

    est_writeRunLevel_CABAC(levelData, levelRDOQ, LUMA_16x16P, lambda_md, kStart, kStop, noCoeff, 0);
  }

  nonzero=FALSE;
  run=-1;
  scan_pos=0;

  for (coeff_ctr=0;coeff_ctr < 256;coeff_ctr++)
  {
    i=pos_scan[coeff_ctr][0];
    j=pos_scan[coeff_ctr][1];

    run++;

    if(input->UseRDO_Q)
    {
      level = levelRDOQ[coeff_ctr];
    }
    else
    {
      level=(abs(M1[j][i])*Q[qp_rem]+qpRound);
      level = level >>(qpShift);
    }

    if (level != 0)
    {
      nonzero=TRUE;

      ACLevel[scan_pos] = sign(level, M1[j][i]);
      ACRun  [scan_pos] = run;
      ++scan_pos;
      run=-1;                     // reset zero level counter

      M1[j][i] = sign(level, M1[j][i]) * (R[qp_rem] << qp_per); 
    }
    else 
      M1[j][i]=0;
  }
  ACLevel[scan_pos] = 0;

  if (nonzero)
  {
    int shift, shiftBy2;
    // Inverse 16x16 transform
    inverse16x16I(M1, curr_res);

    // store the encoded image
    shift=DCTprec+6; 
    shiftBy2=1<<(shift-1);
    for (j=0; j < MB_BLOCK_SIZE; j++)
    {
      for (i=0; i < MB_BLOCK_SIZE; i++)
      {
        short temp;
        temp=(short)((curr_res[j][i]+((int)img->mpr[j][i]<<shift)+shiftBy2)>>shift);
        temp=min(img->max_imgpel_value,max(0,temp));
        img_enc[img->pix_y+j][img->pix_x+i]=temp;

        //printf("%d(%d) ", img->mpr[j][i], curr_res[j][i]);

      }
      //printf("\n");
    }
  }
  else // if (nonzero) => No transformed residual. Just use prediction.
  {
    for (j=0; j < MB_BLOCK_SIZE; j++)
    {

      memcpy(&(img_enc[img->pix_y + j][img->pix_x]),&(img->mpr[j][0]), MB_BLOCK_SIZE * sizeof(imgpel));
    }
  }

  return nonzero;
}

static void forward16x8I(int (*block) [16], int (*tblock)[16])
{
  int x, y;
  int k; 
  int temp[16][16];
  int tmp;
  int bw, bh;

  bw = 16; bh = 8;

  // horizontal 
  for(y = 0; y < bw; y++)
  {
    for(x = 0; x < bh; x++)
    {
      temp[x][y] = 0;
      for(k = 0; k < bw; k++)
        temp[x][y] += img->dct_coeff16I[y][k]*block[x][k];
    }
  }
  
  // vertical 
  for(y = 0; y < bh; y++)
  {
    for(x = 0; x < bw; x++)
    {
      tmp = 0;
      for(k = 0; k < bh; k++)
        tmp += img->dct_coeff8I[y][k]*temp[k][x];
      tblock[y][x] = tmp;
    }
  }
}


static void inverse16x8I(int (*tblock) [16], int (*block)[16])
{
  int x, y, k;
  int temp[16][16];
  int tmp;
  int bw, bh;

  bw = 16; bh = 8;
  // horizontal
  for(y = 0; y < bw; y++)
  {
    for(x = 0; x < bh; x++)
    {
      temp[x][y] = 0;
      for(k = 0; k < bw; k++)
        temp[x][y] += img->dct_coeff16I[k][y]*tblock[x][k];
    }
  }
  
  // vertical 
  for(y = 0; y < bh; y++)
  {
    for(x = 0; x < bw; x++)
    {
      tmp = 0;
      for(k = 0; k < bh; k++)
        tmp += img->dct_coeff8I[k][y]*temp[k][x];
      block[y][x] = tmp;
    }
  }
}

const byte SNGL_SCAN16x8[128][2] = 
{
  { 0, 0}, { 0, 1}, { 1, 0}, { 2, 0}, { 1, 1}, { 0, 2}, { 0, 3}, { 1, 2}, { 2, 1}, { 3, 0}, { 4, 0}, { 3, 1}, { 2, 2}, { 1, 3}, { 0, 4}, { 0, 5}, 
  { 1, 4}, { 2, 3}, { 3, 2}, { 4, 1}, { 5, 0}, { 6, 0}, { 5, 1}, { 4, 2}, { 3, 3}, { 2, 4}, { 1, 5}, { 0, 6}, { 0, 7}, { 1, 6}, { 2, 5}, { 3, 4}, 
  { 4, 3}, { 5, 2}, { 6, 1}, { 7, 0}, { 8, 0}, { 7, 1}, { 6, 2}, { 5, 3}, { 4, 4}, { 3, 5}, { 2, 6}, { 1, 7}, { 2, 7}, { 3, 6}, { 4, 5}, { 5, 4}, 
  { 6, 3}, { 7, 2}, { 8, 1}, { 9, 0}, {10, 0}, { 9, 1}, { 8, 2}, { 7, 3}, { 6, 4}, { 5, 5}, { 4, 6}, { 3, 7}, { 4, 7}, { 5, 6}, { 6, 5}, { 7, 4}, 
  { 8, 3}, { 9, 2}, {10, 1}, {11, 0}, {12, 0}, {11, 1}, {10, 2}, { 9, 3}, { 8, 4}, { 7, 5}, { 6, 6}, { 5, 7}, { 6, 7}, { 7, 6}, { 8, 5}, { 9, 4}, 
  {10, 3}, {11, 2}, {12, 1}, {13, 0}, {14, 0}, {13, 1}, {12, 2}, {11, 3}, {10, 4}, { 9, 5}, { 8, 6}, { 7, 7}, { 8, 7}, { 9, 6}, {10, 5}, {11, 4}, 
  {12, 3}, {13, 2}, {14, 1}, {15, 0}, {15, 1}, {14, 2}, {13, 3}, {12, 4}, {11, 5}, {10, 6}, { 9, 7}, {10, 7}, {11, 6}, {12, 5}, {13, 4}, {14, 3},
  {15, 2}, {15, 3}, {14, 4}, {13, 5}, {12, 6}, {11, 7}, {12, 7}, {13, 6}, {14, 5}, {15, 4}, {15, 5}, {14, 6}, {13, 7}, {14, 7}, {15, 6}, {15, 7}
};

int dct_luma_16x8 (int blk16x8)
{
  int nonzero = FALSE;  
  int i,j,coeff_ctr;
  int scan_pos,run;

  int*  ACLevel = img->cofAC16x8[blk16x8][0];
  int*  ACRun   = img->cofAC16x8[blk16x8][1];

  imgpel **img_enc      = enc_picture->imgY;
  int    (*curr_res)[MB_BLOCK_SIZE] = &img->m7 [blk16x8*MB_BLOCK_SIZE/2]; 

  //Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  //int   max_imgpel_value = img->max_imgpel_value;
  const byte (*pos_scan)[2] = SNGL_SCAN16x8;
  int qp_per    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)/6; //CHA-VG01
  int qp_rem    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)%6; //CHA-VG01
  //int q_bits = Q_BITS + qp_per;
  int DCTprec = 2*DCT16PREC; 

  int R[6] = {40, 45, 50, 57, 63, 71};
  int Q[6] = {410,   364,   328,   287,   260,   231};
  int qpShift, qpRound;
  int level;
  int M1[MB_BLOCK_SIZE][MB_BLOCK_SIZE];

  int k;
  levelDataStruct levelData[256];
  double lambda_md = 0.0;
  double err;
  int lowerInt;
  int levelRDOQ[256];
  int kStart = 0;
  int kStop = 0;
  int noCoeff;
  double normFact = pow(2, 2 * (DCTprec + 14)-15); 

  //  Forward 16x8 transform
  forward16x8I(curr_res, M1);

  // Quantization process
  qpShift=qp_per+DCTprec+8;
  qpRound=(1<<(qpShift))/6;
  
  // quantization and zig-zag scanning 
  if(input->UseRDO_Q)
  {
    if ((img->type == B_SLICE) && img->nal_reference_idc)
    {
      lambda_md = img->lambda_md[5][img->masterQP];  
    }
    else
    {
      lambda_md = img->lambda_md[img->type][img->masterQP];
    }

    noCoeff = 0;
    for (coeff_ctr = 0; coeff_ctr < 128; coeff_ctr++)
    {
      i = pos_scan[coeff_ctr][0];
      j = pos_scan[coeff_ctr][1];

      levelData[coeff_ctr].levelDouble = (int64)absm(M1[j][i]) * Q[qp_rem];
      level = (int)(levelData[coeff_ctr].levelDouble >> qpShift);
      lowerInt = ((levelData[coeff_ctr].levelDouble - (level << qpShift)) < ((int64)1 << (qpShift - 1)))? 1: 0;

      levelData[coeff_ctr].level[0] = 0;
      if (level == 0 && lowerInt == 1)
      {
        levelData[coeff_ctr].noLevels = 1;
      }
      else if (level == 0 && lowerInt == 0)
      {
        levelData[coeff_ctr].level[1] = level + 1;
        levelData[coeff_ctr].noLevels = 2;
        kStop = coeff_ctr;
        noCoeff++;
      }
      else if (level > 0 && lowerInt == 1)
      {
        levelData[coeff_ctr].level[1] = level;
        levelData[coeff_ctr].noLevels = 2;
        kStop = coeff_ctr;
        noCoeff++;
      }
      else
      {
        levelData[coeff_ctr].level[1] = level;
        levelData[coeff_ctr].level[2] = level + 1;
        levelData[coeff_ctr].noLevels = 3;
        kStop = coeff_ctr;
        kStart = coeff_ctr;
        noCoeff++;
      }

      for (k = 0; k < levelData[coeff_ctr].noLevels; k++)
      {      
        err = (double)(levelData[coeff_ctr].level[k] << qpShift) - (double) levelData[coeff_ctr].levelDouble;
        err = err * R[qp_rem];
        levelData[coeff_ctr].errLevel[k] = (err * err) / normFact; 
      }
    }

    est_writeRunLevel_CABAC(levelData, levelRDOQ, LUMA_16x8P, lambda_md, kStart, kStop, noCoeff, 0);
  }

  run=-1;
  scan_pos=0;

  for (coeff_ctr=0;coeff_ctr < 128;coeff_ctr++)
  {
    i=pos_scan[coeff_ctr][0];
    j=pos_scan[coeff_ctr][1];
    
    run++;

    if(input->UseRDO_Q)
    {
      level = levelRDOQ[coeff_ctr];
    }
    else
    {
      level=(abs(M1[j][i])*Q[qp_rem]+qpRound);
      level = level >>(qpShift);
    }

    if (level != 0)
    {
      nonzero=TRUE;

      ACLevel[scan_pos] = sign(level, M1[j][i]);
      ACRun  [scan_pos] = run;
      ++scan_pos;
      run=-1;                     // reset zero level counter

      M1[j][i] = sign(level, M1[j][i]) * (R[qp_rem] << qp_per); 
    }
    else 
      M1[j][i]=0;
  }
  ACLevel[scan_pos] = 0;

  if (nonzero)
  {
    int shift, shiftBy2;
    // Inverse 16x8 transform
    inverse16x8I(M1, curr_res);

    // store the encoded image
    shift=DCTprec+6; 
    shiftBy2=1<<(shift-1);
    for (j=0; j < MB_BLOCK_SIZE/2; j++)
      for (i=0; i < MB_BLOCK_SIZE; i++)
      {
        short temp;

        temp=(short)((curr_res[j][i]+((int)img->mpr[blk16x8*MB_BLOCK_SIZE/2+j][i]<<shift)+shiftBy2)>>shift);
        temp=min(img->max_imgpel_value,max(0,temp));
        img_enc[img->pix_y+blk16x8*MB_BLOCK_SIZE/2+j][img->pix_x+i]=temp;
      }
  }
  else // if (nonzero) => No transformed residual. Just use prediction.
  {
    for (j=0; j < MB_BLOCK_SIZE/2; j++)
    {
      memcpy(&(img_enc[img->pix_y + blk16x8*MB_BLOCK_SIZE/2 + j][img->pix_x]),&(img->mpr[blk16x8*MB_BLOCK_SIZE/2+j][0]), MB_BLOCK_SIZE * sizeof(imgpel));
    }
  }

  return nonzero;
}


static void forward8x16I(int (*block) [16], int (*tblock)[16])
{
  int x, y;
  int k; 
  int temp[16][16];
  int tmp;
  int bw, bh;

  bh = 16; bw = 8;
  // horizontal 
  for(y = 0; y < bw; y++)
  {
    for(x = 0; x < bh; x++)
    {
      temp[x][y] = 0;
      for(k = 0; k < bw; k++)
        temp[x][y] += img->dct_coeff8I[y][k]*block[x][k];
    }
  }
  
  // vertical 
  for(y = 0; y < bh; y++)
  {
    for(x = 0; x < bw; x++)
    {
      tmp = 0;
      for(k = 0; k < bh; k++)
        tmp += img->dct_coeff16I[y][k]*temp[k][x];
      tblock[y][x] = tmp;
    }
  }
}


static void inverse8x16I(int (*tblock) [16], int (*block)[16])
{
  int x, y, k;
  int temp[16][16];
  int tmp;
  int bw, bh;

  bh = 16; bw = 8;
  // horizontal
  for(y = 0; y < bw; y++)
  {
    for(x = 0; x < bh; x++)
    {
      temp[x][y] = 0;
      for(k = 0; k < bw; k++)
        temp[x][y] += img->dct_coeff8I[k][y]*tblock[x][k];
    }
  }
  
  // vertical 
  for(y = 0; y < bh; y++)
  {
    for(x = 0; x < bw; x++)
    {
      tmp = 0;
      for(k = 0; k < bh; k++)
        tmp += img->dct_coeff16I[k][y]*temp[k][x];
      block[y][x] = tmp;
    }
  }
}

int dct_luma_8x16 (int blk8x16)
{
  int nonzero = FALSE;  
  int i,j,coeff_ctr;
  int scan_pos,run;

  int*  ACLevel = img->cofAC8x16[blk8x16][0];
  int*  ACRun   = img->cofAC8x16[blk8x16][1];

  imgpel **img_enc      = enc_picture->imgY;
  int    res[MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  int    (*curr_res)[MB_BLOCK_SIZE] = res;

  //Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  //int   max_imgpel_value = img->max_imgpel_value;
  const byte (*pos_scan)[2] = SNGL_SCAN16x8;
  int qp_per    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)/6; //CHA-VG01
  int qp_rem    = (img->qp + img->bitdepth_luma_qp_scale - MIN_QP)%6; //CHA-VG01
  //int q_bits = Q_BITS + qp_per;
  int DCTprec = 2*DCT16PREC; 

  int R[6] = {40, 45, 50, 57, 63, 71};
  int Q[6] = {410,   364,   328,   287,   260,   231};
  int qpShift, qpRound;
  int shift, shiftBy2;
  int level;
  int M1[MB_BLOCK_SIZE][MB_BLOCK_SIZE];

  int k;
  levelDataStruct levelData[256];
  double lambda_md = 0.0;
  double err;
  int lowerInt;
  int levelRDOQ[256];
  int kStart = 0;
  int kStop = 0;
  int noCoeff;
  double normFact = pow(2, 2 * (DCTprec + 14)-15); 

  //  Forward 8x16 transform
  for(j = 0; j < MB_BLOCK_SIZE; j++)
    for(i = 0; i < MB_BLOCK_SIZE/2; i++)
      curr_res[j][i] = img->m7[j][i+blk8x16*MB_BLOCK_SIZE/2];

  forward8x16I(curr_res, M1);

  // Quantization process
  qpShift=qp_per+DCTprec+8;
  qpRound=(1<<(qpShift))/6;
  
  // quantization and zig-zag scanning 
  if(input->UseRDO_Q)
  {
    if ((img->type == B_SLICE) && img->nal_reference_idc)
    {
      lambda_md = img->lambda_md[5][img->masterQP];  
    }
    else
    {
      lambda_md = img->lambda_md[img->type][img->masterQP];
    }

    noCoeff = 0;
    for (coeff_ctr = 0; coeff_ctr < 128; coeff_ctr++)
    {
      j = pos_scan[coeff_ctr][0];
      i = pos_scan[coeff_ctr][1];

      levelData[coeff_ctr].levelDouble = (int64)absm(M1[j][i]) * Q[qp_rem];
      level = (int)(levelData[coeff_ctr].levelDouble >> qpShift);
      lowerInt = ((levelData[coeff_ctr].levelDouble - (level << qpShift)) < ((int64)1 << (qpShift - 1)))? 1: 0;

      levelData[coeff_ctr].level[0] = 0;
      if (level == 0 && lowerInt == 1)
      {
        levelData[coeff_ctr].noLevels = 1;
      }
      else if (level == 0 && lowerInt == 0)
      {
        levelData[coeff_ctr].level[1] = level + 1;
        levelData[coeff_ctr].noLevels = 2;
        kStop = coeff_ctr;
        noCoeff++;
      }
      else if (level > 0 && lowerInt == 1)
      {
        levelData[coeff_ctr].level[1] = level;
        levelData[coeff_ctr].noLevels = 2;
        kStop = coeff_ctr;
        noCoeff++;
      }
      else
      {
        levelData[coeff_ctr].level[1] = level;
        levelData[coeff_ctr].level[2] = level + 1;
        levelData[coeff_ctr].noLevels = 3;
        kStop = coeff_ctr;
        kStart = coeff_ctr;
        noCoeff++;
      }

      for (k = 0; k < levelData[coeff_ctr].noLevels; k++)
      {      
        err = (double)(levelData[coeff_ctr].level[k] << qpShift) - (double) levelData[coeff_ctr].levelDouble;
        err = err * R[qp_rem];
        levelData[coeff_ctr].errLevel[k] = (err * err) / normFact; 
      }
    }

    est_writeRunLevel_CABAC(levelData, levelRDOQ, LUMA_8x16P, lambda_md, kStart, kStop, noCoeff, 0);
  }

  nonzero=FALSE;

  run=-1;
  scan_pos=0;

  for (coeff_ctr=0;coeff_ctr < 128;coeff_ctr++)
  {
    j=pos_scan[coeff_ctr][0];
    i=pos_scan[coeff_ctr][1];

    run++;

    if(input->UseRDO_Q)
    {
      level = levelRDOQ[coeff_ctr];
    }
    else
    {
      level=(abs(M1[j][i])*Q[qp_rem]+qpRound);
      level = level >>(qpShift);
    }

    if (level != 0)
    {
      nonzero=TRUE;

      ACLevel[scan_pos] = sign(level, M1[j][i]);
      ACRun  [scan_pos] = run;
      ++scan_pos;
      run=-1;                     // reset zero level counter

      M1[j][i] = sign(level, M1[j][i]) * (R[qp_rem] << qp_per); 
    }
    else 
      M1[j][i]=0;
  }
  ACLevel[scan_pos] = 0;

  if (nonzero)
  {
    // Inverse 8x16 transform
    inverse8x16I(M1, curr_res);

    // store the encoded image
    shift=DCTprec+6; 
    shiftBy2=1<<(shift-1);
    for (j=0; j < MB_BLOCK_SIZE; j++)
      for (i=0; i < MB_BLOCK_SIZE/2; i++)
      {
        short temp;

        temp=(short)((curr_res[j][i]+((int)img->mpr[j][i+blk8x16*MB_BLOCK_SIZE/2]<<shift)+shiftBy2)>>shift);
        temp=min(img->max_imgpel_value,max(0,temp));
        img_enc[img->pix_y+j][img->pix_x+blk8x16*MB_BLOCK_SIZE/2+i]=temp;
      }
  }
  else // if (nonzero) => No transformed residual. Just use prediction.
  {
    for (j=0; j < MB_BLOCK_SIZE; j++)
    {
      memcpy(&(img_enc[img->pix_y + j][img->pix_x+ blk8x16*MB_BLOCK_SIZE/2]),&(img->mpr[j][blk8x16*MB_BLOCK_SIZE/2]), MB_BLOCK_SIZE/2 * sizeof(imgpel));
    }
  }

  return nonzero;
}


#endif 
